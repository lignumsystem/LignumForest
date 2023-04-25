# cmake file to compile LignumForest
# For Unix Makefile build system:
# mkdir build
# cd  build
# cmake .. -DCMAKE_BUILD_TYPE=Debug/Release
# Explicitely set c++ compiler
# cmake ..  -DCMAKE_BUILD_TYPE=Debug  -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS=-stdlib=libc++
# make
# For Xcode IDE:
# mkdir xcode
# cd xcode
# cmake .. -G Xcode
# open lignum-forest.xcodeproj   #Opens Xcode IDE
# Build the 'lignum-forest' Product in Xcode. See also Product (menu bar) -> Scheme is set to 'lignum-forest' to allow Run to debug.
# Copy necessary *.fun files and *.txt parameter files to xcode/Debug where 'lignum-forest´ is located.
# Otherwise hard coded files are not found in the program.
# Set command line parameters for 'lignum-forest' in Product (menu bar) -> Scheme -> Edit Schmeme -> Arguments.
# Divide the command line into practical parts for debugging from Arguments -> '+'.
cmake_minimum_required(VERSION 3.12)
# cmake 3.23 for Xcode 13 recommends the policy CMP0114 (warning output)
cmake_policy(SET CMP0114 NEW)
project(lignum-forest DESCRIPTION "Single tree in homogeneous forest")

# Lignum core projects as an ExternalProject. 'lignum-core' target is compiled first
# Note: lignum-core  has to be configured with cmake first. See lignum-core/CMakeLists.txt
include(ExternalProject)
ExternalProject_add(lignum-core
SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..
DOWNLOAD_COMMAND "")

# The '*.cmake' files are produced with 'cmake' configuration for 'lignum-core' (see lignum-core/CMakeLists.txt).  
include(../c++adt/lib/cxxadt.cmake)
include(../Firmament/lib/sky.cmake)
include(../stl-lignum/lib/LGM.cmake)
include(../stl-voxelspace/lib/voxel.cmake)
include(../LEngine/lib/L.cmake)
include(../LEngine/lib/l2c.cmake)
include(../Graphics/lib/Visual.cmake)

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
find_package(Qt5 COMPONENTS Xml Core REQUIRED)
#find_package(OpenGL REQUIRED)
# cmake 3.23 output erroneusly prints 'glut not found'. It is found. 
#find_package(GLUT REQUIRED)

set(SOURCE_FILES lignum-forest.cc
branchfunctor.cc
forest-trees.cc
generate-tree-locations.cc
harvest-stand.cc
pine-em98.cpp
branchfunctor.cc
src/borderforest.cc
src/space.cc
../CrownDensity/src/ScotsPine.cc
)

# Rebuild L system compiler with custom command for 'l2c'
add_custom_command(
OUTPUT ${CMAKE_CURRENT_LIST_DIR}/pine-em98.cpp
DEPENDS ${CMAKE_CURRENT_LIST_DIR}/pine-em98.L l2c
COMMAND ${CMAKE_CURRENT_LIST_DIR}/../LEngine/bin/l2c ${CMAKE_CURRENT_LIST_DIR}/pine-em98.L  ${CMAKE_CURRENT_LIST_DIR}/pine-em98.cpp
)

# Main target, the executable 'lignum-forest'.
add_executable(lignum-forest ${SOURCE_FILES})
target_include_directories(lignum-forest PUBLIC include
${CMAKE_CURRENT_LIST_DIR}/
${CMAKE_CURRENT_LIST_DIR}/include
${CMAKE_CURRENT_LIST_DIR}/../CrownDensity/include
${CMAKE_CURRENT_LIST_DIR}/../stl-voxelspace/include
${CMAKE_CURRENT_LIST_DIR}/../Pine
${CMAKE_CURRENT_LIST_DIR}/../XMLTree
${CMAKE_CURRENT_LIST_DIR}/../LEngine/include
)

# This add_dependencies for lignum-forest seem to generate time stamp based necessary automatic rebuilds
# of 'lignum-forest' and 'lignum-core' projects for Unix Makefile system, i.e. if a file is changed then
# that file and files depending on it will be recompiled.
# Note: 'make clean' is required first before 'make'.
#
# Xcode project created for 'lignum-forest' can compile and rebuild 'lignum-forest' (Product in Xcode) when needed.
# However 'lignum-core' dependency (Targets in Xcode) is compiled once ignoring after
# that time stamp based changes in Lignum core project files when rebuilding lignum-forest.
add_dependencies(lignum-forest lignum-core)

target_link_libraries(lignum-forest
HDF5::HDF5 Qt5::Xml Qt5::Core
L voxel sky LGM cxxadt
)

install(TARGETS lignum-forest DESTINATION "${CMAKE_CURRENT_LIST_DIR}")
