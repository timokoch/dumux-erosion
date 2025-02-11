# SPDX-FileCopyrightInfo: Copyright © Timo Koch
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input params_3d.input)

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
