// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_POROUS_MEDIA_EROSION_MODEL_HH
#define DUMUX_POROUS_MEDIA_EROSION_MODEL_HH

#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/evalgradients.hh>

namespace Dumux::Properties::TTag {
struct PorousMediaErosionModel {};
} // end namespace Dumux::Properties::TTag

namespace Dumux {

// variables
template <class Traits>
class PorousMediaErosionModelVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());
public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;

    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        priVars_ = elemSol[scv.indexInElement()];
        const auto solidity = priVars_[Indices::solidityIdx];
        const auto porosity = 1.0 - solidity;
        permeability_ = problem.permeability(scv.center(), porosity);
    }

    Scalar pressure() const
    { return priVars_[Indices::pressureIdx]; }

    Scalar solidity() const
    { return priVars_[Indices::solidityIdx]; }

    Scalar g() const
    { return priVars_[Indices::gIdx]; }

    Scalar permeability() const
    { return permeability_; }

    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    const PrimaryVariables& priVars() const
    { return priVars_; }

    Scalar extrusionFactor() const
    { return 1.0; }

private:
    PrimaryVariables priVars_;
    Scalar permeability_;
};
} // end namespace Dumux


// equations
namespace Dumux {
template<class TypeTag>
class PorousMediaErosionModelLocalResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    static constexpr int dimWorld = GridView::dimensionworld;
public:
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage;
        const auto kappa = problem.fluidCompressibility();
        const auto gRelaxationRate = problem.fieldRelaxationRate();
        storage[Indices::massBalanceEqIdx] = kappa*volVars.pressure();
        storage[Indices::solidityEqIdx] = volVars.solidity();
        storage[Indices::blurEqIdx] = volVars.g()*gRelaxationRate;
        return storage;
    }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        static_assert(DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>,
            "This local residual is hard-coded to control-volume finite element schemes");

        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradPressure(0.0);
        Dune::FieldVector<Scalar, dimWorld> gradG(0.0);
        Scalar permeability = 0.0;
        const auto& shapeValues = fluxVarCache.shapeValues();
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            gradPressure.axpy(volVars.pressure(), fluxVarCache.gradN(scv.indexInElement()));
            gradG.axpy(volVars.g(), fluxVarCache.gradN(scv.indexInElement()));
            permeability += volVars.permeability()*shapeValues[scv.indexInElement()][0];
        }

        NumEqVector flux;
        flux[Indices::massBalanceEqIdx] = -1.0*vtmv(scvf.unitOuterNormal(), permeability, gradPressure)*scvf.area();
        const auto B = problem.fieldDiffusivity();
        flux[Indices::blurEqIdx] = -1.0*vtmv(scvf.unitOuterNormal(), B, gradG)*scvf.area();

        const auto D = problem.solidDiffusivity();
        if (D > 0.0)
        {
            Dune::FieldVector<Scalar, dimWorld> gradSolidity(0.0);
            for (const auto& scv : scvs(fvGeometry))
                gradSolidity.axpy(elemVolVars[scv].solidity(), fluxVarCache.gradN(scv.indexInElement()));
            flux[Indices::solidityEqIdx] = -1.0*vtmv(scvf.unitOuterNormal(), D, gradSolidity)*scvf.area();
        }

        return flux;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);
        source[Indices::massBalanceEqIdx] = 0.0;
        const auto& volVars = elemVolVars[scv];
        static const Scalar omega = getParam<Scalar>("ModelParameters.Omega");
        static const Scalar phiStar = getParam<Scalar>("ModelParameters.PhiStar");
        static const auto H = [&](const Scalar phi) { return 0.5*(1.0 + std::tanh(omega*(phi - phiStar))); };
        static const Scalar H0 = H(0.0);
        static const Scalar H1 = H(1.0);
        const Scalar psi = (H(volVars.g())-H0)/(H1 - H0);

        const auto geometry = element.geometry();
        const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
        const auto gradP = evalGradients(element, geometry, fvGeometry.gridGeometry(), elemSol, scv.dofPosition())[Indices::pressureIdx];
        const auto pSquared = gradP*gradP;

        const auto erosionRateFactor = problem.erosionRateFactor(scv.center());
        source[Indices::solidityEqIdx] = -volVars.solidity()*std::max(0.0, pSquared - psi)*erosionRateFactor;

        source[Indices::blurEqIdx] = volVars.solidity() - volVars.g();
        source += problem.source(element, fvGeometry, elemVolVars, scv);
        return source;
    }
};
} // end namespace Dumux


// configure the PorousMediaErosionModel
namespace Dumux::Properties {

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PorousMediaErosionModel>
{ using type = PorousMediaErosionModelLocalResidual<TypeTag>; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::PorousMediaErosionModel>
{ using type = double; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PorousMediaErosionModel>
{
    struct type
    {
        struct Indices
        {
            static constexpr int pressureIdx = 0;
            static constexpr int solidityIdx = 1;
            static constexpr int gIdx = 2;

            static constexpr int massBalanceEqIdx = 0;
            static constexpr int solidityEqIdx = 1;
            static constexpr int blurEqIdx = 2;
        };

        static constexpr int numEq() { return 3; }
    };
};

template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::PorousMediaErosionModel>
{
    using type = Dune::FieldVector<
        GetPropType<TypeTag, Properties::Scalar>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PorousMediaErosionModel>
{
    struct Traits
    {
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    };

    using type = PorousMediaErosionModelVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
