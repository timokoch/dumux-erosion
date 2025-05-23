#!/usr/bin/env bash

# reconfigure & build libs
./dune-common/bin/dunecontrol --opts=dumux-erosion/cmake.opts bexec rm -r CMakeCache.txt CMakeFiles
./dune-common/bin/dunecontrol --opts=dumux-erosion/cmake.opts cmake
./dune-common/bin/dunecontrol --opts=dumux-erosion/cmake.opts make -j
