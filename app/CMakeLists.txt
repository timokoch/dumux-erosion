# SPDX-FileCopyrightInfo: Copyright © Timo Koch
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB MSH_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.msh)
file(GLOB INPUT_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.input)
file(GLOB JSON_INPUT_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.json)
message(STATUS "Found MSH_FILES: ${MSH_FILES}")
message(STATUS "Found INPUT_FILES: ${INPUT_FILES}")
message(STATUS "Found JSON_INPUT_FILES: ${JSON_INPUT_FILES}")
dune_symlink_to_source_files(FILES
  ${MSH_FILES}
  ${INPUT_FILES}
  ${JSON_INPUT_FILES}
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
