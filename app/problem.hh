// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_POROUS_MEDIA_EROSION_MODEL_PROBLEM_HH
#define DUMUX_POROUS_MEDIA_EROSION_MODEL_PROBLEM_HH

#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

template<class TypeTag>
class PorousMediaErosionTestProblem : public FVProblem<TypeTag>
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

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    enum class Scenario
    {
        hyperbolic_symmetric, hyperbolic_evaporation,
        half_sphere_boundary, retina,
        sphere_symmetric,
        heart,
        rectangular_inflow_outflow, rectangular_evaporation,
        cuboid_inflow_outflow, cuboid_evaporation
    };

public:
    PorousMediaErosionTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // lens
        lensLowerLeft_ = getParam<GlobalPosition>("BoundaryConditions.LensLowerLeft", GlobalPosition(-101.0));
        lensUpperRight_ = getParam<GlobalPosition>("BoundaryConditions.LensUpperRight", GlobalPosition(-100.0));

        // model parameters
        fluidCompressibility_ = getParam<Scalar>("ModelParameters.FluidCompressibility");
        fieldRelaxationRate_ = getParam<Scalar>("ModelParameters.FieldRelaxationRate");
        solidDiffusivity_ = getParam<Scalar>("ModelParameters.SolidDiffusivity", 0.0);
        fieldDiffusivity_ = getParam<Scalar>("ModelParameters.FieldDiffusivity");

        // boundary conditions
        rampUpTime_ = Chrono::toSeconds(getParam("BoundaryConditions.RampUpTime")).count();
        const auto flowRate = getParam<Scalar>("BoundaryConditions.FlowRate");
        const bool evaporation = getParam<bool>("BoundaryConditions.Evaporation");

        scenario_ = [&]() -> Scenario {
            if constexpr (dim == 2 && dimWorld == 3)
            {
                const auto msh = getParam<std::string>("Grid.File");
                if (msh.find("halfsphere") != std::string::npos)
                    return evaporation ? Scenario::retina : Scenario::half_sphere_boundary;
                else if (msh.find("sphere") != std::string::npos)
                    return Scenario::sphere_symmetric;
                else if (msh.find("hyperbolic") != std::string::npos)
                    return evaporation ? Scenario::hyperbolic_evaporation : Scenario::hyperbolic_symmetric;
                else if (msh.find("heart") != std::string::npos)
                    return Scenario::heart;
            }
            else if (dim == 2)
            {
                return evaporation ? Scenario::rectangular_evaporation
                                   : Scenario::rectangular_inflow_outflow;
            }
            else if (dim == 3)
            {
                return evaporation ? Scenario::cuboid_evaporation
                                   : Scenario::cuboid_inflow_outflow;
            }

            DUNE_THROW(Dune::IOError, "Unknown scenario");
        }();

        for (int i = 0; i < dim; ++i)
            domainSize_[i] = this->gridGeometry().bBoxMax()[i] - this->gridGeometry().bBoxMin()[i];

        // inflow_outflow and evaporation scenarios for 2D and 3D cubes
        if constexpr (dim == dimWorld)
        {
            const auto inflowLength = getParam<Scalar>("BoundaryConditions.InflowLength");
            const auto outflowLength = getParam<Scalar>("BoundaryConditions.OutflowLength");
            const auto inflowArea = dim == 2 ? inflowLength : inflowLength*inflowLength;
            const auto outflowArea = dim == 2 ? outflowLength : outflowLength*outflowLength;
            const auto bottomTopArea = dim == 2 ? domainSize_[0] : domainSize_[0]*domainSize_[1];

            inflowRatePerArea_ = flowRate/inflowArea;
            if (inflowArea < bottomTopArea)
            {
                inflowBottom_ = 0.5*inflowLength;
                inflowSide_ = 0.0;
            }
            else
            {
                inflowBottom_ = 0.5*domainSize_[0];
                inflowSide_ = dim == 2 ? (inflowArea - bottomTopArea)/2.0
                                       : ((inflowArea - bottomTopArea)/4.0)/domainSize_[0];
            }

            outflowRatePerArea_ = evaporation ? 0.0 : flowRate/outflowArea;
            if (outflowArea < bottomTopArea)
            {
                outflowTop_ = 0.5*outflowLength;
                outflowSide_ = 0.0;
            }
            else
            {
                outflowTop_ = 0.5*domainSize_[0];
                outflowSide_ = dim == 2 ? (outflowArea - bottomTopArea)/2.0
                                        : ((outflowArea - bottomTopArea)/4.0)/domainSize_[0];
            }

            const auto domainVolume = dim == 2 ? domainSize_[0]*domainSize_[1] : domainSize_[0]*domainSize_[1]*domainSize_[2];
            evapRatePerVolume_ = evaporation ? flowRate/domainVolume : 0.0;

            if (this->gridGeometry().gridView().comm().rank() == 0)
            {
                std::cout << "\nBoundary conditions:\n";
                std::cout << "-- inflow: " << inflowRatePerArea_*inflowArea << "\n";
                std::cout << "-- outflow: " << outflowRatePerArea_*outflowArea << "\n";
                std::cout << "-- evaporation: " << evapRatePerVolume_*domainVolume << "\n";
            }

            if (inflowRatePerArea_*inflowArea - outflowRatePerArea_*outflowArea - evapRatePerVolume_*domainVolume > 1e-6)
                DUNE_THROW(Dune::Exception, "Mass balance in the boundary conditions not satisfied.");
        }

        // surface geometries for 2D (in 3D) scenarios
        else if constexpr (dim == 2 && dimWorld == 3)
        {
            if (scenario_ == Scenario::half_sphere_boundary)
            {
                const auto r = 5.0;
                const auto a = r*std::sin(M_PI/6.0);
                const auto outflowRingArea = 2*M_PI*a;
                inflowRatePerArea_ = 0.0;
                outflowRatePerArea_ = flowRate/outflowRingArea;

                heightSphericalCapIn_ = 0.1;
                const auto areaSpehricalCap = 2.0*M_PI*heightSphericalCapIn_*r;
                evapRatePerAreaIn_ = flowRate/areaSpehricalCap;
                evapRatePerVolume_ = 0.0;
            }
            else if (scenario_ == Scenario::retina)
            {
                const auto r = 5.0;
                const auto h = r*(1.0 - std::cos(M_PI/6.0));
                const auto outflowMissingSphericalCap = 2.0*M_PI*h*r;

                inflowRatePerArea_ = 0.0;
                outflowRatePerArea_ = 0.0;
                heightSphericalCapIn_ = 0.1;
                const auto areaSpehricalCap = 2.0*M_PI*heightSphericalCapIn_*r;
                const auto sphereArea = 4.0*M_PI*r*r;
                evapRatePerAreaIn_ = flowRate/areaSpehricalCap;
                evapRatePerVolume_ = flowRate/(sphereArea - outflowMissingSphericalCap);
            }
            else if (scenario_ == Scenario::sphere_symmetric)
            {
                const auto r = 5.0;
                heightSphericalCapIn_ = 0.1;
                const auto areaSpehricalCap = 2.0*M_PI*heightSphericalCapIn_*r;

                inflowRatePerArea_ = flowRate/areaSpehricalCap;
                outflowRatePerArea_ = flowRate/areaSpehricalCap;
                evapRatePerAreaIn_ = 0.0;
                evapRatePerVolume_ = 0.0;
            }
            else if (scenario_ == Scenario::hyperbolic_symmetric)
            {
                const auto r = 5.0;
                const auto ringArea = 2.0*M_PI*r;
                inflowRatePerArea_ = flowRate/ringArea;
                outflowRatePerArea_ = flowRate/ringArea;
            }
            else if (scenario_ == Scenario::hyperbolic_evaporation)
                DUNE_THROW(Dune::NotImplemented, "Scenario::hyperbolic_evaporation");
        }
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        // TODO Dirichlet scenario
        // if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6)
        // {
        //     values.setDirichlet(0);
        //     values.setNeumann(1);
        //     values.setNeumann(2);
        // }
        // else
        values.setAllNeumann();

        return values;
    }

    // only called on Dirichlet boundaries
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        return { 0.0, 0.0, 0.0 }; // second and third should be unused
    }

    // only called on Dirichlet boundaries
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);

        if (isIntializationPhase_)
            return values;

        const auto rampFactor = std::min(1.0, 1.0/rampUpTime_*time_);

        // inflow_outflow and evaporation scenarios for 2D and 3D cubes
        if constexpr (dim == dimWorld)
        {
            if (globalPos[1] < this->gridGeometry().bBoxMin()[1] + 1e-6 && globalPos[0] < inflowBottom_ && globalPos[0] > -inflowBottom_)
                values[0] = -inflowRatePerArea_*rampFactor;
            else if (onSide_(globalPos) && globalPos[1] < this->gridGeometry().bBoxMin()[1] + inflowSide_)
                values[0] = -inflowRatePerArea_*rampFactor;
            else if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6 && globalPos[0] < outflowTop_ && globalPos[0] > -outflowTop_)
                values[0] = outflowRatePerArea_*rampFactor;
            else if (onSide_(globalPos) && globalPos[1] > this->gridGeometry().bBoxMax()[1] - outflowSide_)
                values[0] = outflowRatePerArea_*rampFactor;
        }

        // surface geometries for 2D (in 3D) scenarios
        else if constexpr (dim == 2 && dimWorld == 3)
        {
            if (globalPos[1] < this->gridGeometry().bBoxMin()[1] + 1e-6)
                values[0] = -inflowRatePerArea_*rampFactor;
            else if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6)
                values[0] = outflowRatePerArea_*rampFactor;
        }

        return values;
    }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        if (isIntializationPhase_)
            return { 0.0, 0.0, 0.0 };

        const auto rampFactor = std::min(1.0, 1.0/rampUpTime_*time_);

        // inflow_outflow and evaporation scenarios for 2D and 3D cubes
        if constexpr (dim == dimWorld)
            return { -evapRatePerVolume_*rampFactor, 0.0, 0.0 };

        // surface geometries for 2D (in 3D) scenarios
        else if constexpr (dim == 2 && dimWorld == 3)
        {
            assert(heightSphericalCapIn_ > 0.0);
            if (globalPos[1] < this->gridGeometry().bBoxMin()[1] + heightSphericalCapIn_)
                return { evapRatePerAreaIn_*rampFactor, 0.0, 0.0 };
            else
                return { -evapRatePerVolume_*rampFactor, 0.0, 0.0 };
        }
    }

    Scalar permeability(const GlobalPosition& globalPos, const Scalar porosity) const
    {
        const auto permeabilityFactor = inLens(globalPos) ? eps_ : 1.0;
        const auto solidFraction = 1.0 - porosity;
        return permeabilityFactor*porosity*porosity*porosity/(solidFraction*solidFraction); // Kozeny-Carman
    }

    Scalar erosionRateFactor(const GlobalPosition& globalPos) const
    { return inLens(globalPos) ? 0.0 : 1.0; }

    Scalar fluidCompressibility() const { return fluidCompressibility_; }
    Scalar fieldRelaxationRate() const { return fieldRelaxationRate_; }
    Scalar solidDiffusivity() const { return solidDiffusivity_; }
    Scalar fieldDiffusivity() const { return fieldDiffusivity_; }

    void stopInitializationPhase() { isIntializationPhase_ = false; }
    void startInitializationPhase() { isIntializationPhase_ = true; }

    void setSolidDiffusivity(const Scalar D) { solidDiffusivity_ = D; }
    void setTime(Scalar time) { time_ = time; }

    // lens occlusion
    bool inLens(const GlobalPosition& globalPos) const
    {
        bool inLens = true;
        for (int i = 0; i < dim; ++i)
            inLens = inLens && globalPos[i] > lensLowerLeft_[i] - eps_
                            && globalPos[i] < lensUpperRight_[i] + eps_;
        return inLens;
    }

    // boundary conditions
    bool onSide_(const GlobalPosition& globalPos) const
    {
        bool onSide = false;
        for (int i = 0; i < dim-1; ++i)
            onSide = onSide || globalPos[i] > this->gridGeometry().bBoxMin()[i] + eps_
                            || globalPos[i] < this->gridGeometry().bBoxMax()[i] - eps_;
        return onSide;
    }

private:

    Scalar time_ = 0.0;
    Scalar rampUpTime_;
    Scalar inflowRatePerArea_, inflowBottom_, inflowSide_;
    Scalar outflowRatePerArea_, outflowTop_, outflowSide_;
    Scalar evapRatePerVolume_;
    Scalar evapRatePerAreaIn_ = 0.0;
    Scalar heightSphericalCapIn_ = 0.0;

    GlobalPosition domainSize_;

    GlobalPosition lensUpperRight_;
    GlobalPosition lensLowerLeft_;

    // model parameters
    Scalar fluidCompressibility_;
    Scalar fieldRelaxationRate_;
    Scalar solidDiffusivity_;
    Scalar fieldDiffusivity_;

    bool isIntializationPhase_ = true;

    Scenario scenario_;

    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
