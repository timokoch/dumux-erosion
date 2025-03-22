# SPDX-FileCopyrightInfo: Copyright © Timo Koch
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES
  params.input params_3d.input params_2d3d.input
  sphere.msh parabol.msh halfsphere.msh
  surface.msh surface_2.msh surface_3.msh surface_4.msh torus.msh
  run_scenario.py
)

dumux_add_test(
  NAME test_erosion_2d
  COMPILE_DEFINITIONS
    GRID_TYPE=Dune::YaspGrid<2,Dune::EquidistantOffsetCoordinates<double,2>>
  SOURCES main.cc
)

dumux_add_test(
  NAME test_erosion_3d
  COMPILE_DEFINITIONS
    GRID_TYPE=Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>
  SOURCES main.cc
)

dumux_add_test(
  NAME test_erosion_md_2d_in_3d
  COMPILE_DEFINITIONS
    GRID_TYPE=Dune::ALUGrid<2,3,Dune::simplex,Dune::conforming>
  SOURCES main.cc
)
