
#!/usr/bin/env python
# SPDX-FileCopyrightText: Copyright © Timo Koch
# SPDX-License-Identifier: GPL-3.0-or-later

import subprocess
import click
import os
from copy import deepcopy

SCENARIOS = {
    "retina": {
        "Executable": "run_erosion_2d_in_3d",
        "ParameterFile": "run_erosion_2d_in_3d.input",
        "Parameters": {
            "Problem.Name": "retina",
            "Grid.File": "halfsphere_fine.msh",
            "BoundaryConditions.Evaporation": "true",
            "BoundaryConditions.PointSourcesStrength": "0.7",
            "BoundaryConditions.PointSourcePositions": "2.5 -4.33012701892 0.0",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "Json",
            "BoundaryConditions.JsonObstacleFile": "obstacles_retina.json",
            "ModelParameters.CorrelationLength": "0.0025",
            "ModelParameters.FieldDiffusivity": "0.1",
            "TimeLoop.CheckPointPeriod": "0.25s",
            "TimeLoop.TEnd": "20.0s",
        },
    },
    "heart": {
        "Executable": "run_erosion_2d_in_3d",
        "ParameterFile": "run_erosion_2d_in_3d.input",
        "Parameters": {
            "Problem.Name": "heart",
            "Grid.File": "../../data/heart.msh",
            "BoundaryConditions.Evaporation": "true",
            "BoundaryConditions.PointSourcesStrength": "4 4",
            "BoundaryConditions.PointSourcePositions": "-14.6428 21.1399 37.5829 -23.9962 9.63808 23.1206",
            "BoundaryConditions.RampUpTime": "40.0",
            "BoundaryConditions.ObstacleType": "None",
            "ModelParameters.CorrelationLength": "0.6",
            "ModelParameters.FieldDiffusivity": "20.0",
            "TimeLoop.CheckPointPeriod": "0.25s",
            "TimeLoop.TEnd": "40.0s",
        },
    },
    "hyperbolic": {
        "Executable": "run_erosion_2d_in_3d",
        "ParameterFile": "run_erosion_2d_in_3d.input",
        "Parameters": {
            "Problem.Name": "hyperbolic",
            "Grid.File": "hyperbolic.msh",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.FlowRate": "0.4",
            "ModelParameters.CorrelationLength": "0.05",
            "ModelParameters.FieldDiffusivity": "0.2",
            "TimeLoop.CheckPointPeriod": "0.25s",
            "TimeLoop.TEnd": "20.0s",
        },
    },
    "sphere": {
        "Executable": "run_erosion_2d_in_3d",
        "ParameterFile": "run_erosion_2d_in_3d.input",
        "Parameters": {
            "Problem.Name": "sphere",
            "Grid.File": "sphere.msh",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.FlowRate": "0.4",
            "ModelParameters.CorrelationLength": "0.05",
            "ModelParameters.FieldDiffusivity": "0.2",
            "TimeLoop.CheckPointPeriod": "0.25s",
            "TimeLoop.TEnd": "20.0s",
        },
    },
    "torus": {
        "Executable": "run_erosion_2d_in_3d",
        "ParameterFile": "run_erosion_2d_in_3d.input",
        "Parameters": {
            "Problem.Name": "torus",
            "Grid.File": "torus.msh",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.FlowRate": "0.4",
            "ModelParameters.CorrelationLength": "0.05",
            "ModelParameters.FieldDiffusivity": "0.2",
            "TimeLoop.CheckPointPeriod": "0.25s",
            "TimeLoop.TEnd": "20.0s",
        },
    },
    "surface1": {
        "Executable": "run_erosion_2d_in_3d",
        "ParameterFile": "run_erosion_2d_in_3d.input",
        "Parameters": {
            "Problem.Name": "surface1",
            "Grid.File": "surface.msh",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.FlowRate": "0.5",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.FieldDiffusivity": "0.2",
            "ModelParameters.InitialStddev": "0.2",
            "LinearSolver.Iterations": "2000",
            "TimeLoop.TEnd": "20.0s",
            "TimeLoop.CheckPointPeriod": "0.25s",
        },
    },
    "surface2": {
        "Executable": "run_erosion_2d_in_3d",
        "ParameterFile": "run_erosion_2d_in_3d.input",
        "Parameters": {
            "Problem.Name": "surface2",
            "Grid.File": "surface_2.msh",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.FlowRate": "0.5",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.FieldDiffusivity": "0.2",
            "ModelParameters.InitialStddev": "0.2",
            "LinearSolver.Iterations": "2000",
            "TimeLoop.TEnd": "20.0s",
            "TimeLoop.CheckPointPeriod": "0.25s",
        },
    },
    "surface3": {
        "Executable": "run_erosion_2d_in_3d",
        "ParameterFile": "run_erosion_2d_in_3d.input",
        "Parameters": {
            "Problem.Name": "surface3",
            "Grid.File": "surface_3.msh",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.FlowRate": "0.5",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.FieldDiffusivity": "0.2",
            "ModelParameters.InitialStddev": "0.2",
            "LinearSolver.Iterations": "2000",
            "TimeLoop.TEnd": "20.0s",
            "TimeLoop.CheckPointPeriod": "0.25s",
        },
    },
    "surface4": {
        "Executable": "run_erosion_2d_in_3d",
        "ParameterFile": "run_erosion_2d_in_3d.input",
        "Parameters": {
            "Problem.Name": "surface4",
            "Grid.File": "surface_4.msh",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.FlowRate": "0.5",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.FieldDiffusivity": "0.2",
            "ModelParameters.InitialStddev": "0.2",
            "LinearSolver.Iterations": "2000",
            "TimeLoop.TEnd": "20.0s",
            "TimeLoop.CheckPointPeriod": "0.25s",
        },
    },
    "2d_line_line": {
        "Executable": "run_erosion_2d",
        "ParameterFile": "run_erosion_2d.input",
        "Parameters": {
            "Problem.Name": "2d_line_line",
            "Grid.LowerLeft": "-5 -5",
            "Grid.UpperRight": "5 5",
            "Grid.Cells": "200 200",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.InflowLength": "10.0",
            "BoundaryConditions.OutflowLength": "10.0",
            "BoundaryConditions.FlowRate": "0.8",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.Omega": "6.5",
            "ModelParameters.PhiStar": "0.7",
            "ModelParameters.FieldDiffusivity": "0.3",
            "ModelParameters.InitialStddev": "0.1",
            "LinearSolver.Iterations": "2000",
            "TimeLoop.TEnd": "20.0s",
            "TimeLoop.CheckPointPeriod": "1.0s",
        },
    },
    "2d_point_point": {
        "Executable": "run_erosion_2d",
        "ParameterFile": "run_erosion_2d.input",
        "Parameters": {
            "Problem.Name": "2d_point_point",
            "Grid.LowerLeft": "-5 -5",
            "Grid.UpperRight": "5 5",
            "Grid.Cells": "200 200",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.InflowLength": "0.1",
            "BoundaryConditions.OutflowLength": "0.1",
            "BoundaryConditions.FlowRate": "0.8",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.Omega": "6.5",
            "ModelParameters.PhiStar": "0.7",
            "ModelParameters.FieldDiffusivity": "0.3",
            "ModelParameters.InitialStddev": "0.1",
            "LinearSolver.Iterations": "2000",
            "TimeLoop.TEnd": "20.0s",
            "TimeLoop.CheckPointPeriod": "1.0s",
        },
    },
    "2d_point_line": {
        "Executable": "run_erosion_2d",
        "ParameterFile": "run_erosion_2d.input",
        "Parameters": {
            "Problem.Name": "2d_point_line",
            "Grid.LowerLeft": "-5 -5",
            "Grid.UpperRight": "5 5",
            "Grid.Cells": "200 200",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.InflowLength": "0.1",
            "BoundaryConditions.OutflowLength": "10.0",
            "BoundaryConditions.FlowRate": "0.8",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.Omega": "6.5",
            "ModelParameters.PhiStar": "0.7",
            "ModelParameters.FieldDiffusivity": "0.3",
            "ModelParameters.InitialStddev": "0.1",
            "LinearSolver.Iterations": "2000",
            "TimeLoop.TEnd": "40.0s",
            "TimeLoop.CheckPointPeriod": "1.0s",
        },
    },
    "2d_point_long_line": {
        "Executable": "run_erosion_2d",
        "ParameterFile": "run_erosion_2d.input",
        "Parameters": {
            "Problem.Name": "2d_point_long_line",
            "Grid.LowerLeft": "-5 -5",
            "Grid.UpperRight": "5 5",
            "Grid.Cells": "200 200",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.InflowLength": "0.1",
            "BoundaryConditions.OutflowLength": "20.0",
            "BoundaryConditions.FlowRate": "0.8",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.Omega": "6.5",
            "ModelParameters.PhiStar": "0.7",
            "ModelParameters.FieldDiffusivity": "0.3",
            "ModelParameters.InitialStddev": "0.1",
            "LinearSolver.Iterations": "2000",
            "TimeLoop.TEnd": "20.0s",
            "TimeLoop.CheckPointPeriod": "1.0s",
        },
    },
    "2d_point_evaporation": {
        "Executable": "run_erosion_2d",
        "ParameterFile": "run_erosion_2d.input",
        "Parameters": {
            "Problem.Name": "2d_point_evaporation",
            "Grid.LowerLeft": "-5 -5",
            "Grid.UpperRight": "5 5",
            "Grid.Cells": "200 200",
            "BoundaryConditions.Evaporation": "true",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.InflowLength": "0.1",
            "BoundaryConditions.FlowRate": "0.8",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.Omega": "6.5",
            "ModelParameters.PhiStar": "0.7",
            "ModelParameters.FieldDiffusivity": "0.3",
            "ModelParameters.InitialStddev": "0.1",
            "LinearSolver.Iterations": "2000",
            "TimeLoop.TEnd": "20.0s",
            "TimeLoop.CheckPointPeriod": "1.0s",
        },
    },
    "2d_point_point_long": {
        "Executable": "run_erosion_2d",
        "ParameterFile": "run_erosion_2d.input",
        "Parameters": {
            "Problem.Name": "2d_point_point_long",
            "Grid.LowerLeft": "-5 -20",
            "Grid.UpperRight": "5 20",
            "Grid.Cells": "200 1000",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "BoundaryConditions.ObstacleType": "None",
            "BoundaryConditions.InflowLength": "0.1",
            "BoundaryConditions.OutflowLength": "0.1",
            "BoundaryConditions.FlowRate": "0.8",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.Omega": "6.5",
            "ModelParameters.PhiStar": "0.7",
            "ModelParameters.FieldDiffusivity": "0.3",
            "ModelParameters.InitialStddev": "0.1",
            "LinearSolver.Iterations": "2000",
            "TimeLoop.TEnd": "20.0s",
            "TimeLoop.CheckPointPeriod": "1.0s",
        },
    },
}

_more_scenarios = {}
for k, v in SCENARIOS.items():
    if "surface" in k:
        _more_scenarios[k + "_evaporation"] = deepcopy(v)
        _more_scenarios[k + "_evaporation"]["Parameters"]["BoundaryConditions.Evaporation"] = "true"
        _more_scenarios[k + "_evaporation"]["Parameters"]["Problem.Name"] = v["Parameters"]["Problem.Name"] + "_evaporation"

for q in ["0.1", "0.3", "0.6", "0.8", "1.0", "1.2", "1.5", "2.0", "10.0"]:
    _more_scenarios[f"2d_point_long_line_q{q}"] = deepcopy(SCENARIOS["2d_point_long_line"])
    _more_scenarios[f"2d_point_long_line_q{q}"]["Parameters"]["BoundaryConditions.FlowRate"] = q

SCENARIOS.update(_more_scenarios)

@click.command()
@click.argument("scenario", type=click.Choice(SCENARIOS.keys()), required=True)
@click.option("-p", "--num_processes", type=int, default=8)
def run_scenario(scenario, num_processes):
    s = SCENARIOS[scenario]
    command = ["mpiexec", "-np", f"{num_processes}", s["Executable"], s["ParameterFile"]]
    params = s.get("Parameters", {})
    for k, v in params.items():
        command.extend([f"-{k}", v])
    print(f"Running scenario: {scenario}")
    print(f"-- command: {' '.join(command)}")
    cmd_env = os.environ.copy()
    cmd_env["DUMUX_NUM_THREADS"] = str(8//num_processes)
    cmd_env["OMP_NUM_THREADS"] = str(8//num_processes)
    subprocess.run(command, env=cmd_env)

if __name__ == "__main__":
    run_scenario()
