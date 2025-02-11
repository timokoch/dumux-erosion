// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>
#include <type_traits>
#include <random>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>

#include <dumux/discretization/box.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dune/grid/yaspgrid.hh>
#include <dumux/io/chrono.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#if HAVE_DUNE_ALUGRID
#include <dumux/io/grid/gridmanager_alu.hh>
#endif

#include "model.hh"

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

public:
    PorousMediaErosionTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // boundary conditions
        rampUpTime_ = Chrono::toSeconds(getParam("BoundaryConditions.RampUpTime")).count();
        const auto flowRate = getParam<Scalar>("BoundaryConditions.FlowRate");
        const bool evaporation = getParam<bool>("BoundaryConditions.Evaporation");
        const auto xLength = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        const auto yLength = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];

        if constexpr (dim == dimWorld)
        {
            const auto inflowLength = getParam<Scalar>("BoundaryConditions.InflowLength");
            const auto outflowLength = getParam<Scalar>("BoundaryConditions.OutflowLength");

            inflowRatePerLength_ = flowRate/inflowLength;
            if (inflowLength < xLength)
            {
                inflowBottom_ = 0.5*inflowLength;
                inflowSide_ = 0.0;
            }
            else
            {
                inflowBottom_ = 0.5*xLength;
                inflowSide_ = 0.5*(inflowLength - xLength);
            }

            outflowRatePerLength_ = evaporation ? 0.0 : flowRate/outflowLength;
            if (outflowLength < xLength)
            {
                outflowBottom_ = 0.5*outflowLength;
                outflowSide_ = 0.0;
            }
            else
            {
                outflowBottom_ = 0.5*xLength;
                outflowSide_ = 0.5*(outflowLength - xLength);
            }

            evapRatePerArea_ = evaporation ? flowRate/(yLength*xLength) : 0.0;

            if (this->gridGeometry().gridView().comm().rank() == 0)
            {
                std::cout << "\nBoundary conditions:\n";
                std::cout << "-- inflow: " << inflowRatePerLength_*inflowLength << "\n";
                std::cout << "-- outflow: " << outflowRatePerLength_*outflowLength << "\n";
                std::cout << "-- evaporation: " << evapRatePerArea_*yLength*xLength << "\n";
            }

            if (inflowRatePerLength_*inflowLength - outflowRatePerLength_*outflowLength - evapRatePerArea_*yLength*xLength > 1e-6)
                DUNE_THROW(Dune::Exception, "Mass balance in the boundary conditions not satisfied.");
        }
        else
        {
            inflowRatePerLength_ = 0.0;
            outflowRatePerLength_ = 0.0;
            evapRatePerArea_ = flowRate/(yLength*0.1*xLength*0.1);
        }

        // model parameters
        fluidCompressibility_ = getParam<Scalar>("ModelParameters.FluidCompressibility");
        fieldRelaxationRate_ = getParam<Scalar>("ModelParameters.FieldRelaxationRate");
        solidDiffusivity_ = getParam<Scalar>("ModelParameters.SolidDiffusivity");
        fieldDiffusivity_ = getParam<Scalar>("ModelParameters.FieldDiffusivity");
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);

        if (isIntializationPhase_)
            return values;

        const auto rampFactor = std::min(1.0, 1.0/rampUpTime_*time_);
        if (globalPos[1] < this->gridGeometry().bBoxMin()[1] + 1e-6 && globalPos[0] < inflowBottom_ && globalPos[0] > -inflowBottom_)
            values[0] = -inflowRatePerLength_*rampFactor;
        else if (onSide_(globalPos) && globalPos[1] < this->gridGeometry().bBoxMin()[1] + inflowSide_)
            values[0] = -inflowRatePerLength_*rampFactor;
        else if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6)
            values[0] = outflowRatePerLength_*rampFactor;
        else if (onSide_(globalPos) && globalPos[1] > this->gridGeometry().bBoxMax()[1] - outflowSide_)
            values[0] = outflowRatePerLength_*rampFactor;

        return values;
    }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        if (isIntializationPhase_)
            return { 0.0, 0.0, 0.0 };

        const auto rampFactor = std::min(1.0, 1.0/rampUpTime_*time_);
        if constexpr (dim == dimWorld)
            return { -evapRatePerArea_*rampFactor, 0.0, 0.0 };
        else
        {
            if (globalPos[1] < this->gridGeometry().bBoxMin()[1] + 0.5)
                return { evapRatePerArea_*rampFactor, 0.0, 0.0 };
            else if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - 0.5)
                return { -evapRatePerArea_*rampFactor, 0.0, 0.0 };
            else
                return { 0.0, 0.0, 0.0 };
        }
    }

    Scalar fluidCompressibility() const { return fluidCompressibility_; }
    Scalar fieldRelaxationRate() const { return fieldRelaxationRate_; }
    Scalar solidDiffusivity() const { return solidDiffusivity_; }
    Scalar fieldDiffusivity() const { return fieldDiffusivity_; }

    void stopInitializationPhase() { isIntializationPhase_ = false; }
    void startInitializationPhase() { isIntializationPhase_ = true; }

    void setSolidDiffusivity(const Scalar D) { solidDiffusivity_ = D; }
    void setTime(Scalar time) { time_ = time; }

private:
    // boundary conditions
    bool onSide_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < this->gridGeometry().bBoxMin()[0] + 1e-6 ||
               globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-6;
    }
    Scalar time_ = 0.0;
    Scalar rampUpTime_;
    Scalar inflowRatePerLength_, inflowBottom_, inflowSide_;
    Scalar outflowRatePerLength_, outflowBottom_, outflowSide_;
    Scalar evapRatePerArea_;

    // model parameters
    Scalar fluidCompressibility_;
    Scalar fieldRelaxationRate_;
    Scalar solidDiffusivity_;
    Scalar fieldDiffusivity_;

    bool isIntializationPhase_ = true;
};

} // end namespace Dumux

namespace Dumux::Properties::TTag {

#ifndef GRID_TYPE
#define GRID_TYPE Dune::YaspGrid<2,Dune::EquidistantOffsetCoordinates<double,2>>
#endif

struct PorousMediaErosionTest
{
    using InheritsFrom = std::tuple<PorousMediaErosionModel, BoxModel>;

    using Grid = GRID_TYPE;

    template<class TypeTag>
    using Problem = PorousMediaErosionTestProblem<TypeTag>;

    using EnableGridGeometryCache = std::true_type;
    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
};
} // end namespace Dumux::Properties

// functor for data communication with MPI
struct MinScatter
{
    template<class A, class B>
    static void apply(A& a, const B& b)
    { a[1] = std::min(a[1], b[1]); }
};

// create the random initial solution
template<class SolutionVector, class GridGeometry>
SolutionVector createInitialSolution(const GridGeometry& gg)
{
    SolutionVector sol(gg.numDofs());

    // Generate random number and add processor offset
    // For sequential runs `rank` always returns `0`.
    std::mt19937 gen(0); // seed is 0 for deterministic results
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int n = 0; n < sol.size(); ++n)
    {
        sol[n][0] = 0.0;
        sol[n][1] = 0.8 + 0.1*(0.5-dis(gen)) + gg.gridView().comm().rank();
        sol[n][2] = 0.8;
    }

    // We take the value of the processor with the minimum rank
    // and subtract the rank offset
    if (gg.gridView().comm().size() > 1)
    {
        Dumux::VectorCommDataHandle<
            typename GridGeometry::VertexMapper,
            SolutionVector,
            GridGeometry::GridView::dimension,
            MinScatter
        > minHandle(gg.vertexMapper(), sol);
        gg.gridView().communicate(minHandle, Dune::All_All_Interface, Dune::ForwardCommunication);

        // Remove processor offset
        for (int n = 0; n < sol.size(); ++n)
            sol[n][1] -= std::floor(sol[n][1]);
    }
    return sol;
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::PorousMediaErosionTest;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    GridManager<Grid> gridManager;
    gridManager.init();
    const auto rank = gridManager.grid().leafGridView().comm().rank();

    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);
    auto sol = createInitialSolution<SolutionVector>(*gridGeometry);
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.pressure(); }, "p");
    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.solidity(); }, "phi");
    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.g(); }, "g");
    vtkWriter.write(-1.0);

    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(
        Chrono::toSeconds("-1.0s"), Chrono::toSeconds("0.1s"), Chrono::toSeconds(getParam("TimeLoop.TEnd"))
    );
    timeLoop->setMaxTimeStepSize(0.1);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Solver = Dumux::NewtonSolver<Assembler, LinearSolver>;

    auto oldSol = sol;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, oldSol);
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
    Solver solver(assembler, linearSolver);

    // initialization phase idea: we start with a random field and then diffusive such that we get
    // a Gaussian random field for the solidity to start with
    const auto correlationLength = getParam<Scalar>("ModelParameters.CorrelationLength");
    const auto D = correlationLength*correlationLength/4.0; // Running diffusion for 1 time unit
    problem->startInitializationPhase();
    problem->setSolidDiffusivity(D);

    if (rank == 0)
        std::cout << "\nStarting initialization phase with solid diffusivity " << problem->solidDiffusivity() << "\n" << std::endl;

    timeLoop->start(); do {
        problem->setTime(timeLoop->time() + timeLoop->timeStepSize());
        solver.solve(sol, *timeLoop);
        oldSol = sol;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();
        if (timeLoop->timeStepIndex() != 10)
            vtkWriter.write(timeLoop->time());
        timeLoop->reportTimeStep();
    } while (timeLoop->timeStepIndex() < 10);

    // re-initialize timeloop and problem
    timeLoop->setTime(0.0);
    timeLoop->setTimeStepSize(Chrono::toSeconds(getParam("TimeLoop.InitialTimeStepSize")));
    timeLoop->setMaxTimeStepSize(Chrono::toSeconds(getParam("TimeLoop.MaxTimeStepSize")));
    problem->stopInitializationPhase();
    problem->setSolidDiffusivity(getParam<Scalar>("ModelParameters.SolidDiffusivity", 0.0));

    // re-initialize solution
    for (int n = 0; n < sol.size(); ++n)
    {
        sol[n][0] = 0.0;
        sol[n][2] = 0.8;
    }

    oldSol = sol;
    gridVariables->init(sol);
    vtkWriter.write(0.0);

    if (rank == 0)
        std::cout << "\nStarting main phase with solid diffusivity " << problem->solidDiffusivity() << "\n" << std::endl;

    do {
        problem->setTime(timeLoop->time() + timeLoop->timeStepSize());
        solver.solve(sol, *timeLoop);
        oldSol = sol;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();
        vtkWriter.write(timeLoop->time());
        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(solver.suggestTimeStepSize(timeLoop->timeStepSize()));
    } while (!timeLoop->finished());

    // print the final report
    timeLoop->finalize(gridGeometry->gridView().comm());
    return 0;
}
