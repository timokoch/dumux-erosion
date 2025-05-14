// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_POROUS_MEDIA_EROSION_MODEL_PROBLEM_HH
#define DUMUX_POROUS_MEDIA_EROSION_MODEL_PROBLEM_HH

#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <algorithm>
#include <functional>
#include <array>
#include <cmath>
#include <iostream>
#include <concepts>

#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/geometry/distance.hh>

#include <dumux/io/json.hh>

namespace Dumux {

template<class Scalar, int dimWorld>
class Obstacles {
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
public:
    Obstacles() = default;

    template<class T>
    requires std::predicate<T, GlobalPosition>
    void addObstacle(T&& t)
    { obstacles_.emplace_back(std::forward<T>(t)); }

    bool inside(const GlobalPosition& point) const {
        return std::any_of(obstacles_.begin(), obstacles_.end(),
            [&](const auto& query) { return query(point); }
        );
    }

    std::vector<std::function<bool(const GlobalPosition&)>> obstacles_;
};

} // end namespace Dumux

namespace Dumux::Obstacle {

template<class Scalar, int dimWorld>
void addBall(const Dune::FieldVector<Scalar, dimWorld>& center,
             const Scalar radius,
             Obstacles<Scalar, dimWorld>& obstacles)
{
    obstacles.addObstacle([=](const auto& point){
        const auto dist2 = (center-point).two_norm2();
        return dist2 <= radius*radius;
    });
}

template<class Scalar, int dimWorld>
void addBox(const Dune::FieldVector<Scalar, dimWorld>& lowerLeft,
            const Dune::FieldVector<Scalar, dimWorld>& upperRight,
            Obstacles<Scalar, dimWorld>& obstacles)
{
    obstacles.addObstacle([=](const auto& point){
        for (int i = 0; i < dimWorld; ++i)
            if (point[i] < lowerLeft[i] + 1e-7 || point[i] > upperRight[i] - 1e-7)
                return false;
        return true;
    });
}

template<class Scalar, int dimWorld>
void addCylinder(const Dune::FieldVector<Scalar, dimWorld>& lowerCenter,
                 const Dune::FieldVector<Scalar, dimWorld>& upperCenter,
                 const Scalar radius,
                 Obstacles<Scalar, dimWorld>& obstacles)
{
    obstacles.addObstacle([=](const auto& point){
        const auto sdf = [&]{
            const auto ba = upperCenter - lowerCenter;
            const auto pa = point - lowerCenter;
            const Scalar baba = ba*ba;
            const Scalar paba = pa*ba;
            const Scalar x = (pa*baba - ba*paba).two_norm() - radius*baba;
            const Scalar y = std::abs(paba-baba*0.5) - baba*0.5;
            const Scalar x2 = x*x;
            const Scalar y2 = y*y*baba;
            const Scalar d = (std::max(x,y) < 0.0) ? -std::min(x2, y2) : (((x > 0.0) ? x2 : 0.0) + ((y > 0.0) ? y2 : 0.0));
            return std::copysign(std::sqrt(std::abs(d))/baba, d);
        }();
        return sdf <= 0.0;
    });
}

} // end namespace Dumux::Obstacle

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
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    enum class Scenario
    {
        hyperbolic_symmetric, hyperbolic_evaporation,
        half_sphere_boundary, retina,
        sphere_symmetric,
        torus_symmetric, torus_evaporation,
        heart,
        rectangular_inflow_outflow, rectangular_evaporation,
        cuboid_inflow_outflow, cuboid_evaporation,
        wavy_surface_inflow_outflow, wavy_surface_evaporation
    };

public:
    PorousMediaErosionTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // obstacles
        const auto obstacleType = getParam<std::string>("BoundaryConditions.ObstacleType", "None");
        if (obstacleType == "Ball")
        {
            const auto center = getParam<GlobalPosition>("BoundaryConditions.BallObstacleCenter");
            const auto radius = getParam<Scalar>("BoundaryConditions.BallObstacleRadius");
            Obstacle::addBall(center, radius, obstacles_);
        }
        else if (obstacleType == "Box")
        {
            const auto lowerLeft = getParam<GlobalPosition>("BoundaryConditions.BoxObstacleLowerLeft");
            const auto upperRight = getParam<GlobalPosition>("BoundaryConditions.BoxUpperRight");
            Obstacle::addBox(lowerLeft, upperRight, obstacles_);
        }
        else if (obstacleType == "Json")
        {
            const auto filename = getParam<std::string>("BoundaryConditions.JsonObstacleFile");
            auto data = Json::JsonTree::parse(std::ifstream{filename});
            for (const auto& obstacle : data)
            {
                const auto type = obstacle["type"].template get<std::string>();
                if (type == "Ball")
                {
                    const auto center = obstacle["center"].template get<GlobalPosition>();
                    const auto radius = obstacle["radius"].template get<Scalar>();
                    Obstacle::addBall(center, radius, obstacles_);
                }
                else if (type == "Box")
                {
                    const auto lowerLeft = obstacle["lower_left"].template get<GlobalPosition>();
                    const auto upperRight = obstacle["upper_right"].template get<GlobalPosition>();
                    Obstacle::addBox(lowerLeft, upperRight, obstacles_);
                }
                else if (type == "Cylinder")
                {
                    const auto lowerCenter = obstacle["lower_center"].template get<GlobalPosition>();
                    const auto upperCenter = obstacle["upper_center"].template get<GlobalPosition>();
                    const auto radius = obstacle["radius"].template get<Scalar>();
                    Obstacle::addCylinder(lowerCenter, upperCenter, radius, obstacles_);
                }
                else
                {
                    DUNE_THROW(Dune::NotImplemented, "Unknown obstacle type: " << type);
                }
            }
        }
        else if (obstacleType != "None")
        {
            DUNE_THROW(Dune::NotImplemented, "Unknown obstacle type: " << obstacleType);
        }

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
                else if (msh.find("torus") != std::string::npos)
                    return evaporation ? Scenario::torus_evaporation : Scenario::torus_symmetric;
                else if (msh.find("hyperbolic") != std::string::npos)
                    return evaporation ? Scenario::hyperbolic_evaporation : Scenario::hyperbolic_symmetric;
                else if (msh.find("heart") != std::string::npos)
                    return Scenario::heart;
                else if (msh.find("surface") != std::string::npos)
                    return evaporation ? Scenario::wavy_surface_evaporation : Scenario::wavy_surface_inflow_outflow;
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
                heightSphericalCapIn_ = 0.025;
                auto fvGeometry = localView(this->gridGeometry());
                Scalar areaSphericalCapInflow = 0.0, areaRingOutflow = 0.0;
                for (const auto& element : elements(this->gridGeometry().gridView(), Dune::Partitions::interior))
                {
                    fvGeometry.bindElement(element);
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        if (scv.center()[1] < this->gridGeometry().bBoxMin()[1] + heightSphericalCapIn_)
                            areaSphericalCapInflow += scv.volume();
                    }

                    for (const auto& scvf : scvfs(fvGeometry))
                    {
                        if (scvf.boundary() && scvf.center()[1] > this->gridGeometry().bBoxMax()[1] - eps_)
                            areaRingOutflow += scvf.area();
                    }
                }

                areaSphericalCapInflow = this->gridGeometry().gridView().comm().sum(areaSphericalCapInflow);
                areaRingOutflow = this->gridGeometry().gridView().comm().sum(areaRingOutflow);

                inflowRatePerArea_ = 0.0;
                outflowRatePerArea_ = flowRate/areaRingOutflow;
                evapRatePerAreaIn_ = flowRate/areaSphericalCapInflow;
                evapRatePerAreaOut_ = 0.0;
                evapRatePerVolume_ = 0.0;
            }
            else if (scenario_ == Scenario::wavy_surface_inflow_outflow)
            {
                Scalar inflowCurveArea = 0.0, outflowCurveArea = 0.0;
                auto fvGeometry = localView(this->gridGeometry());
                for (const auto& element : elements(this->gridGeometry().gridView(), Dune::Partitions::interior))
                {
                    fvGeometry.bindElement(element);
                    for (const auto& scvf : scvfs(fvGeometry))
                    {
                        if (scvf.boundary())
                        {
                            if (scvf.ipGlobal()[1] < this->gridGeometry().bBoxMin()[1] + 1e-6)
                                inflowCurveArea += scvf.area();
                            if (scvf.ipGlobal()[1] > this->gridGeometry().bBoxMax()[1] - 1e-6)
                                outflowCurveArea += scvf.area();
                        }
                    }
                }

                inflowCurveArea = this->gridGeometry().gridView().comm().sum(inflowCurveArea);
                outflowCurveArea = this->gridGeometry().gridView().comm().sum(outflowCurveArea);

                inflowRatePerArea_ = flowRate/inflowCurveArea;
                outflowRatePerArea_ = flowRate/outflowCurveArea;
                evapRatePerAreaIn_ = 0.0;
                evapRatePerAreaOut_ = 0.0;
                evapRatePerVolume_ = 0.0;
            }
            else if (scenario_ == Scenario::wavy_surface_evaporation)
            {
                Scalar inflowCurveArea = 0.0, evaporationArea = 0.0;
                auto fvGeometry = localView(this->gridGeometry());
                for (const auto& element : elements(this->gridGeometry().gridView(), Dune::Partitions::interior))
                {
                    fvGeometry.bindElement(element);
                    for (const auto& scvf : scvfs(fvGeometry))
                    {
                        if (scvf.boundary())
                        {
                            if (scvf.ipGlobal()[1] < this->gridGeometry().bBoxMin()[1] + 1e-6)
                                inflowCurveArea += scvf.area();
                        }
                    }

                    evaporationArea += element.geometry().volume();
                }

                inflowCurveArea = this->gridGeometry().gridView().comm().sum(inflowCurveArea);
                evaporationArea = this->gridGeometry().gridView().comm().sum(evaporationArea);

                inflowRatePerArea_ = flowRate/inflowCurveArea;
                outflowRatePerArea_ = 0.0;
                evapRatePerAreaIn_ = 0.0;
                evapRatePerAreaOut_ = 0.0;
                evapRatePerVolume_ = flowRate/evaporationArea;
            }
            else if (scenario_ == Scenario::sphere_symmetric || scenario_ == Scenario::torus_symmetric)
            {
                heightSphericalCapIn_ = 0.025;
                auto fvGeometry = localView(this->gridGeometry());
                Scalar areaSphericalCapInflow = 0.0, areaSphericalCapOutflow = 0.0;
                for (const auto& element : elements(this->gridGeometry().gridView(), Dune::Partitions::interior))
                {
                    fvGeometry.bindElement(element);
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        if (scv.center()[1] < this->gridGeometry().bBoxMin()[1] + heightSphericalCapIn_)
                            areaSphericalCapInflow += scv.volume();
                        if (scv.center()[1] > this->gridGeometry().bBoxMax()[1] - heightSphericalCapIn_)
                            areaSphericalCapOutflow += scv.volume();
                    }
                }

                areaSphericalCapInflow = this->gridGeometry().gridView().comm().sum(areaSphericalCapInflow);
                areaSphericalCapOutflow = this->gridGeometry().gridView().comm().sum(areaSphericalCapOutflow);

                inflowRatePerArea_ = 0.0;
                outflowRatePerArea_ = 0.0;
                evapRatePerAreaIn_ = flowRate/areaSphericalCapInflow;
                evapRatePerAreaOut_ = flowRate/areaSphericalCapOutflow;
                evapRatePerVolume_ = 0.0;
            }
            else if (scenario_ == Scenario::hyperbolic_symmetric)
            {
                Scalar inflowRingArea = 0.0, outflowRingArea = 0.0;
                auto fvGeometry = localView(this->gridGeometry());
                for (const auto& element : elements(this->gridGeometry().gridView(), Dune::Partitions::interior))
                {
                    fvGeometry.bindElement(element);
                    for (const auto& scvf : scvfs(fvGeometry))
                    {
                        if (scvf.boundary())
                        {
                            if (scvf.ipGlobal()[1] < this->gridGeometry().bBoxMin()[1] + 1e-6)
                                inflowRingArea += scvf.area();
                            if (scvf.ipGlobal()[1] > this->gridGeometry().bBoxMax()[1] - 1e-6)
                                outflowRingArea += scvf.area();
                        }
                    }
                }

                inflowRingArea = this->gridGeometry().gridView().comm().sum(inflowRingArea);
                outflowRingArea = this->gridGeometry().gridView().comm().sum(outflowRingArea);

                inflowRatePerArea_ = flowRate/inflowRingArea;
                outflowRatePerArea_ = flowRate/outflowRingArea;
                evapRatePerAreaIn_ = 0.0;
                evapRatePerAreaOut_ = 0.0;
                evapRatePerVolume_ = 0.0;
            }
            else if (scenario_ == Scenario::hyperbolic_evaporation)
                DUNE_THROW(Dune::NotImplemented, "Scenario::hyperbolic_evaporation");

            else if (scenario_ == Scenario::heart || scenario_ == Scenario::torus_evaporation || scenario_ == Scenario::retina)
            {
                const auto strengths = getParam<std::vector<Scalar>>("BoundaryConditions.PointSourcesStrength");
                const auto positions = getParam<std::vector<Scalar>>("BoundaryConditions.PointSourcePositions");
                std::vector<std::pair<Scalar, GlobalPosition>> pointSourceCandidates;
                for (int i = 0; i < strengths.size(); ++i)
                    pointSourceCandidates.push_back(std::make_pair(strengths[i], GlobalPosition{positions[3*i], positions[3*i + 1], positions[3*i + 2]} ));

                for (const auto& p : pointSourceCandidates)
                {
                    const auto [dist, eIdx] = closestEntity(p.second, this->gridGeometry().boundingBoxTree(), 0.05*0.05);
                    const int found = static_cast<int>(eIdx != 0);
                    const int totalFound = this->gridGeometry().gridView().comm().sum(found);
                    if (found)
                    {
                        if (this->gridGeometry().element(eIdx).partitionType() == Dune::PartitionType::InteriorEntity)
                        {
                            const auto rank = this->gridGeometry().gridView().comm().rank();
                            std::cout << "Added point source on rank " << rank << " at position " << p.second << " on element " << eIdx << "\n";
                            pointSources_.push_back(std::make_pair(eIdx, p.first/static_cast<double>(totalFound)));
                        }
                    }
                }

                const auto totalStrength = std::accumulate(strengths.begin(), strengths.end(), 0.0);
                const auto totalArea = [&]{
                    Scalar area = 0.0;
                    for (const auto& element : elements(this->gridGeometry().gridView(), Dune::Partitions::interior))
                        if (const auto geo = element.geometry(); !insideObstacle(geo.center()))
                            area += geo.volume();
                    area = this->gridGeometry().gridView().comm().sum(area);
                    return area;
                }();

                inflowRatePerArea_ = 0.0;
                outflowRatePerArea_ = 0.0;
                evapRatePerAreaIn_ = 0.0;
                evapRatePerVolume_ = totalStrength/totalArea;
            }

            else
                DUNE_THROW(Dune::NotImplemented, "Unknown scenario");
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
            if constexpr (dim == 2)
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
            if constexpr (dim == 3)
            {
                if (globalPos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + 1e-6 && globalPos[0] < inflowBottom_ && globalPos[0] > -inflowBottom_ && globalPos[dimWorld-2] < inflowBottom_ && globalPos[dimWorld-2] > -inflowBottom_)
                    values[0] = -inflowRatePerArea_*rampFactor;
                else if (onSide_(globalPos) && globalPos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + inflowSide_)
                    values[0] = -inflowRatePerArea_*rampFactor;
                else if (globalPos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - 1e-6 && globalPos[0] < outflowTop_ && globalPos[0] > -outflowTop_ && globalPos[dimWorld-2] < outflowTop_ && globalPos[dimWorld-2] > -outflowTop_)
                    values[0] = outflowRatePerArea_*rampFactor;
                else if (onSide_(globalPos) && globalPos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - outflowSide_)
                    values[0] = outflowRatePerArea_*rampFactor;
            }
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

    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        if (isIntializationPhase_)
            return { 0.0, 0.0, 0.0 };

        const auto& globalPos = scv.center();
        const auto rampFactor = std::min(1.0, 1.0/rampUpTime_*time_);

        // inflow_outflow and evaporation scenarios for 2D and 3D cubes
        if constexpr (dim == dimWorld)
            return { -evapRatePerVolume_*rampFactor, 0.0, 0.0 };

        // surface geometries for 2D (in 3D) scenarios
        else if constexpr (dim == 2 && dimWorld == 3)
        {
            if (scenario_ == Scenario::heart || scenario_ == Scenario::torus_evaporation || scenario_ == Scenario::retina)
            {
                for (const auto& [eIdx, strength] : pointSources_)
                    if (fvGeometry.elementIndex() == eIdx)
                        return { rampFactor*(strength/element.geometry().volume() - evapRatePerVolume_), 0.0, 0.0 };

                if (!insideObstacle(globalPos))
                    return { -evapRatePerVolume_*rampFactor, 0.0, 0.0 };

                return { 0.0, 0.0, 0.0 };
            }

            assert(heightSphericalCapIn_ > 0.0);
            if (scenario_ == Scenario::sphere_symmetric || scenario_ == Scenario::torus_symmetric || scenario_ == Scenario::half_sphere_boundary)
                if (globalPos[1] < this->gridGeometry().bBoxMin()[1] + heightSphericalCapIn_)
                    return { (evapRatePerAreaIn_-evapRatePerVolume_)*rampFactor, 0.0, 0.0 };

            if (scenario_ == Scenario::sphere_symmetric || scenario_ == Scenario::torus_symmetric)
                if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - heightSphericalCapIn_)
                    return { (-evapRatePerAreaOut_-evapRatePerVolume_)*rampFactor, 0.0, 0.0 };

            return { -evapRatePerVolume_*rampFactor, 0.0, 0.0 };
        }
    }

    Scalar permeability(const GlobalPosition& globalPos, const Scalar porosity) const
    {
        const auto poro = std::max(porosity, 1e-2);
        const auto solidFraction = 1.0 - poro;
        return poro*poro*poro/(solidFraction*solidFraction); // Kozeny-Carman
    }

    Scalar erosionRateFactor(const GlobalPosition& globalPos) const
    { return insideObstacle(globalPos) ? 0.0 : 1.0; }

    Scalar fluidCompressibility() const { return fluidCompressibility_; }
    Scalar fieldRelaxationRate() const { return fieldRelaxationRate_; }
    Scalar solidDiffusivity() const { return solidDiffusivity_; }
    Scalar fieldDiffusivity() const { return fieldDiffusivity_; }

    void stopInitializationPhase() { isIntializationPhase_ = false; }
    void startInitializationPhase() { isIntializationPhase_ = true; }

    void setSolidDiffusivity(const Scalar D) { solidDiffusivity_ = D; }
    void setTime(Scalar time) { time_ = time; }

    // obstacle occlusion
    bool insideObstacle(const GlobalPosition& globalPos) const
    { return obstacles_.inside(globalPos); }

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
    Scalar evapRatePerAreaOut_ = 0.0;
    Scalar heightSphericalCapIn_ = 0.0;

    GlobalPosition domainSize_;

    Obstacles<Scalar, dimWorld> obstacles_;

    // model parameters
    Scalar fluidCompressibility_;
    Scalar fieldRelaxationRate_;
    Scalar solidDiffusivity_;
    Scalar fieldDiffusivity_;

    std::vector<std::pair<std::size_t, Scalar>> pointSources_;

    bool isIntializationPhase_ = true;

    Scenario scenario_;

    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
