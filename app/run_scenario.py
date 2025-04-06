#!/usr/bin/env python

import subprocess
import click
import os

SCENARIOS = {
    "retina": {
        "Executable": "test_erosion_md_2d_in_3d",
        "ParameterFile": "params_2d3d.input",
        "Parameters": {
            "Problem.Name": "retina",
            "Grid.File": "halfsphere.msh",
            "BoundaryConditions.Evaporation": "true",
            "BoundaryConditions.PointSourcesStrength": "0.4",
            "BoundaryConditions.PointSourcePosition": "2.5 -4.33012701892 0.0",
            "BoundaryConditions.RampUpTime": "10.0",
            "ModelParameters.CorrelationLength": "0.0025",
            "ModelParameters.FieldDiffusivity": "0.001",
        },
    },
    "heart": {
        "Executable": "test_erosion_md_2d_in_3d",
        "ParameterFile": "params_2d3d.input",
        "Parameters": {
            "Problem.Name": "heart",
            "Grid.File": "heart.msh",
            "BoundaryConditions.Evaporation": "true",
            "BoundaryConditions.PointSourcesStrength": "4 4",
            "BoundaryConditions.PointSourcePosition": "-14.6428 21.1399 37.5829 -23.9962 9.63808 23.1206",
            "BoundaryConditions.RampUpTime": "20.0",
            "ModelParameters.CorrelationLength": "0.8",
            "ModelParameters.FieldDiffusivity": "0.4",
        },
    },
    "surface1": {
        "Executable": "test_erosion_md_2d_in_3d",
        "ParameterFile": "params_2d3d.input",
        "Parameters": {
            "Problem.Name": "surface",
            "Grid.File": "surface.msh",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.FieldDiffusivity": "0.001",
        },
    },
    "surface2": {
        "Executable": "test_erosion_md_2d_in_3d",
        "ParameterFile": "params_2d3d.input",
        "Parameters": {
            "Problem.Name": "surface2",
            "Grid.File": "surface_2.msh",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.FieldDiffusivity": "0.001",
        },
    },
    "surface3": {
        "Executable": "test_erosion_md_2d_in_3d",
        "ParameterFile": "params_2d3d.input",
        "Parameters": {
            "Problem.Name": "surface3",
            "Grid.File": "surface_3.msh",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.FieldDiffusivity": "0.001",
        },
    },
    "surface4": {
        "Executable": "test_erosion_md_2d_in_3d",
        "ParameterFile": "params_2d3d.input",
        "Parameters": {
            "Problem.Name": "surface4",
            "Grid.File": "surface_4.msh",
            "BoundaryConditions.Evaporation": "false",
            "BoundaryConditions.RampUpTime": "10.0",
            "ModelParameters.CorrelationLength": "0.1",
            "ModelParameters.FieldDiffusivity": "0.001",
        },
    },
}

_more_scenarios = {}
for k, v in SCENARIOS.items():
    if "surface" in k:
        _more_scenarios[k + "_evaporation"] = v.copy()
        _more_scenarios[k + "_evaporation"]["Parameters"]["BoundaryConditions.Evaporation"] = "true"
        _more_scenarios[k + "_evaporation"]["Parameters"]["Problem.Name"] = v["Parameters"]["Problem.Name"] + "_evaporation"

SCENARIOS.update(_more_scenarios)

@click.command()
@click.argument("scenario", type=click.Choice(SCENARIOS.keys()), required=True)
@click.option("-p", "--num_processes", type=int, default=8)
def run_scenario(scenario, num_processes):
    s = SCENARIOS[scenario]
    command = ["mpirun", "-np", f"{num_processes}", s["Executable"], s["ParameterFile"]]
    params = s.get("Parameters", {})
    for k, v in params.items():
        command.extend([f"-{k}", v])
    subprocess.run(command, env=dict(os.environment, OMP_NUM_THREADS=1))

if __name__ == "__main__":
    run_scenario()
