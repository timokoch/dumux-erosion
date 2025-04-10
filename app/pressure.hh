// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_POROUS_MEDIA_EROSION_PRESSURE_SOLVER_HH
#define DUMUX_POROUS_MEDIA_EROSION_PRESSURE_SOLVER_HH

#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/defaultlocaloperator.hh>

namespace Dumux::Properties::TTag {
struct PorousMediaErosionDarcyModel {};
} // end namespace Dumux::Properties::TTag

namespace Dumux {

// variables
template <class Traits>
class PorousMediaErosionDarcyModelVolumeVariables
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
        const auto solidity = problem.solidity(element, scv);
        const auto porosity = 1.0 - solidity;
        permeability_ = problem.permeability(scv.center(), porosity);
    }

    Scalar pressure() const
    { return priVars_[Indices::pressureIdx]; }

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
class PorousMediaErosionDarcyModelLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
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
        storage[Indices::massBalanceEqIdx] = kappa*volVars.pressure();
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
        Scalar permeability = 0.0;
        const auto& shapeValues = fluxVarCache.shapeValues();
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            gradPressure.axpy(volVars.pressure(), fluxVarCache.gradN(scv.indexInElement()));
            permeability += volVars.permeability()*shapeValues[scv.indexInElement()][0];
        }

        NumEqVector flux;
        flux[Indices::massBalanceEqIdx] = -1.0*vtmv(scvf.unitOuterNormal(), permeability, gradPressure)*scvf.area();
        return flux;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);
        source += problem.source(element, fvGeometry, elemVolVars, scv);
        return source;
    }
};
} // end namespace Dumux


// configure the PorousMediaErosionDarcyModel
namespace Dumux::Properties {

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PorousMediaErosionDarcyModel>
{ using type = PorousMediaErosionDarcyModelLocalResidual<TypeTag>; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::PorousMediaErosionDarcyModel>
{ using type = double; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::PorousMediaErosionDarcyModel>
{
    struct type
    {
        struct Indices
        {
            static constexpr int pressureIdx = 0;
            static constexpr int massBalanceEqIdx = 0;
        };

        static constexpr int numEq() { return 1; }
    };
};

template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::PorousMediaErosionDarcyModel>
{
    using type = Dune::FieldVector<
        GetPropType<TypeTag, Properties::Scalar>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PorousMediaErosionDarcyModel>
{
    struct Traits
    {
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    };

    using type = PorousMediaErosionDarcyModelVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties

namespace Dumux {

// solve for pressure with given solidity field to make
// solution compatible with the Darcy model. This is required
// to solve the pure Neumann problem with regularization
template<class TypeTag>
class PorousMediaErosionPressureSolverProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
public:
    PorousMediaErosionPressureSolverProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
       obstaclePermeability_ = 1e-9;
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    // only called on Dirichlet boundaries
    PrimaryVariables dirichlet(const Element& element, const SubControlVolume& scv) const
    { return boundaryPressure_[scv.dofIndex()]; }

    Scalar permeability(const GlobalPosition& globalPos, const Scalar porosity) const
    {
        if (insideObstacle_(globalPos))
            return obstaclePermeability_;

        const auto solidFraction = 1.0 - porosity;
        return porosity*porosity*porosity/(solidFraction*solidFraction); // Kozeny-Carman
    }

    Scalar solidity(const Element& element, const SubControlVolume& scv) const
    { return solidity_[scv.dofIndex()]; }

    template<class FullSolutionVector>
    void setBoundaryAndInitialConditions(const FullSolutionVector& sol)
    {
        boundaryPressure_.resize(sol.size());
        solidity_.resize(sol.size());
        for (int n = 0; n < sol.size(); ++n)
        {
            boundaryPressure_[n] = sol[n][0];
            solidity_[n] = sol[n][1];
        }
    }

    template<class Func>
    requires std::predicate<Func, const Dune::FieldVector<Scalar, dimWorld>&>
    void setObstables(Func&& obs)
    { insideObstacle_ = std::forward<Func>(obs); }

    Scalar fluidCompressibility() const { return 0.0; }

private:
    std::vector<Scalar> boundaryPressure_;
    std::vector<Scalar> solidity_;
    std::function<bool(const Dune::FieldVector<Scalar, dimWorld>&)> insideObstacle_;
    Scalar obstaclePermeability_;
};

} // end namespace Dumux

namespace Dumux::Properties::TTag {

#ifndef GRID_TYPE
#define GRID_TYPE Dune::YaspGrid<2,Dune::EquidistantOffsetCoordinates<double,2>>
#endif

struct PorousMediaErosionPressureSolver
{
    using InheritsFrom = std::tuple<PorousMediaErosionDarcyModel, BoxModel>;
    using Grid = GRID_TYPE;

    template<class TypeTag>
    using Problem = PorousMediaErosionPressureSolverProblem<TypeTag>;

    using EnableGridGeometryCache = std::true_type;
    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
};

} // end namespace Dumux::Properties

namespace Dumux {

template<class FullSolutionVector, class GridGeometry, class ObstacleFunc>
void reinitializePressure(FullSolutionVector& sol, std::shared_ptr<GridGeometry> gg, ObstacleFunc&& obs)
{
    using TypeTag = Dumux::Properties::TTag::PorousMediaErosionPressureSolver;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    auto problem = std::make_shared<Problem>(gg);
    problem->setBoundaryAndInitialConditions(sol);
    problem->setObstables(std::forward<ObstacleFunc>(obs));

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector solPressure(gg->numDofs());

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gg);
    gridVariables->init(solPressure);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Solver = Dumux::NewtonSolver<Assembler, LinearSolver>;

    auto assembler = std::make_shared<Assembler>(problem, gg, gridVariables);
    auto linearSolver = std::make_shared<LinearSolver>(gg->gridView(), gg->dofMapper());
    Solver solver(assembler, linearSolver);

    solver.solve(solPressure);

    for (int n = 0; n < sol.size(); ++n)
        sol[n][0] = solPressure[n][0];
}

} // end namespace Dumux

#endif
