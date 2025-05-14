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

#include <dumux/discretization/box.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dune/grid/yaspgrid.hh>
#include <dumux/io/chrono.hh>
#include <dumux/io/gridwriter.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/loadsolution.hh>

#if HAVE_DUNE_ALUGRID
#include <dumux/io/grid/gridmanager_alu.hh>
#endif

#include "model.hh"
#include "problem.hh"

namespace Dumux {

template<class Grid>
struct OutputModuleTraits
{
    static constexpr auto format = IO::Format::pvd_with(IO::Format::vti.with({
        .encoder = IO::Encoding::raw,
        .compressor = IO::Compression::none,
    }));
};

#if HAVE_DUNE_ALUGRID
template<int dim, int dimWorld>
struct OutputModuleTraits<Dune::ALUGrid<dim, dimWorld, Dune::simplex, Dune::conforming>>
{
    static constexpr auto format = IO::Format::pvd_with(IO::Format::vtu.with({
        .encoder = IO::Encoding::raw,
        .compressor = IO::Compression::none,
    }));
};
#endif

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
template<class SolutionVector, class GridGeometry, class Scalar>
SolutionVector createInitialSolution(const GridGeometry& gg, const Scalar initStddev = 0.1)
{
    SolutionVector sol(gg.numDofs());

    // Generate random number and add processor offset
    // For sequential runs `rank` always returns `0`.
    const int rank = gg.gridView().comm().rank();
    std::mt19937 gen(rank);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int n = 0; n < sol.size(); ++n)
    {
        sol[n][0] = 0.0;
        sol[n][1] = 0.8 + initStddev*(0.5-dis(gen)) + 100*rank;
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
            sol[n][1] -= std::floor(sol[n][1]/100)*100;
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
    const auto initStddev = getParam<Scalar>("ModelParameters.InitialStddev", 0.1);
    const auto initialSolutionFile = getParam<std::string>("Problem.ReadInitialSolutionFrom", "");
    SolutionVector sol = [&]{
        if (initialSolutionFile.empty())
            return createInitialSolution<SolutionVector>(*gridGeometry, initStddev);
        else
        {
            SolutionVector sol(gridGeometry->numDofs());
            loadSolution(sol, initialSolutionFile,
                [](int pvIdx, int){
                    constexpr std::array<std::string_view, 3> n{{"p", "phi", "g"}};
                    return std::string{n[pvIdx]};
                }, *gridGeometry
            );
            return sol;
        }
    }();

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    std::vector<int> insideObstacle(gridGeometry->gridView().size(0), 0);
    for (const auto& element : elements(gridGeometry->gridView(), Dune::Partitions::interior))
        if (problem->insideObstacle(element.geometry().center()))
            insideObstacle[gridGeometry->elementMapper().index(element)] = 1;

    IO::OutputModule vtkWriter(OutputModuleTraits<Grid>::format, *gridVariables, sol, problem->name());

    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.pressure(); }, "p");
    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.solidity(); }, "phi");
    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.g(); }, "g");
    vtkWriter.addField(insideObstacle, "obstacle");

    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(
        Chrono::toSeconds("-1.0s"), Chrono::toSeconds("0.1s"), Chrono::toSeconds(getParam("TimeLoop.TEnd"))
    );
    timeLoop->setMaxTimeStepSize(0.1);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Solver = Dumux::NewtonSolver<Assembler, LinearSolver>;

    auto oldSol = sol;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, oldSol);
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
    Solver solver(assembler, linearSolver);

    problem->startInitializationPhase();
    if (initialSolutionFile.empty())
    {
        // initialization phase idea: we start with a random field and then diffusive such that we get
        // a Gaussian random field for the solidity to start with
        const auto correlationLength = getParam<Scalar>("ModelParameters.CorrelationLength");
        const auto D = correlationLength*correlationLength/(2.0*Grid::dimension); // Running diffusion for 1 time unit
        problem->setSolidDiffusivity(D);

        if (rank == 0)
            std::cout << "\nStarting initialization phase with solid diffusivity " << problem->solidDiffusivity() << "\n" << std::endl;

        timeLoop->start(); do {
            problem->setTime(timeLoop->time() + timeLoop->timeStepSize());
            solver.solve(sol, *timeLoop);
            oldSol = sol;
            gridVariables->advanceTimeStep();
            timeLoop->advanceTimeStep();
            timeLoop->reportTimeStep();
        } while (timeLoop->timeStepIndex() < 10);

        // re-initialize solution
        const Scalar obtacleVolumeFraction = getParam<Scalar>("ModelParameters.ObstacleVolumeFraction", 0.8);
        auto fvGeometry = localView(*gridGeometry);
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);
            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
            {
                const auto dofIndex = scv.dofIndex();
                if (insideObstacle[eIdx] == 1)
                    sol[dofIndex][1] = obtacleVolumeFraction;

                sol[dofIndex][0] = 0.0;
                sol[dofIndex][2] = 0.8;
            }
        }
    }
    else
    {
        // re-initialize pressure
        for (int n = 0; n < sol.size(); ++n)
            sol[n][0] = 0.0;
    }

    // re-initialize timeloop and problem
    timeLoop->setTime(0.0);
    timeLoop->setTimeStepSize(Chrono::toSeconds(getParam("TimeLoop.InitialTimeStepSize")));
    timeLoop->setMaxTimeStepSize(Chrono::toSeconds(getParam("TimeLoop.MaxTimeStepSize")));
    problem->stopInitializationPhase();
    problem->setSolidDiffusivity(getParam<Scalar>("ModelParameters.SolidDiffusivity", 0.0));
    timeLoop->setPeriodicCheckPoint(Chrono::toSeconds(getParam("TimeLoop.CheckPointPeriod", "2.0s")));

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
        if (timeLoop->isCheckPoint())
            vtkWriter.write(timeLoop->time());
        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(solver.suggestTimeStepSize(timeLoop->timeStepSize()));
    } while (!timeLoop->finished());

    // print the final report
    timeLoop->finalize(gridGeometry->gridView().comm());
    return 0;
}
