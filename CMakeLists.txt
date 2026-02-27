# -----------------------------------------
# CMAKE file to compile 'lignum-forest' and
# install it in LignumForest directory
# -----------------------------------------
# For Unix Makefile build system for Debug:
# mkdir debug
# cd  debug
# cmake .. -DCMAKE_BUILD_TYPE=Debug
# Optionally explicitely set c++ compiler
# cmake ..  -DCMAKE_BUILD_TYPE=Debug  -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS=-stdlib=libc++
# make install
# ------------------------------------------
# For Unix Makefile build system for Release:
# mkdir release
# cd release
# cmake .. -DCMAKE_BUILD_TYPE=Release
# make install
# ---------------------------------------------
# For Xcode IDE first create Xcode project file:
# mkdir xcode
# cd xcode
# cmake .. -G Xcode
# -------------------------------------
# Open Xcode IDE
# open lignum-forest.xcodeproj
# -------------------------------------
# Build the 'lignum-forest' Product in Xcode for debugging.
# See also Xcode -> Product (in menu bar) -> Scheme is set to 'lignum-forest' to allow Run to debug.
# Copy necessary *.fun files and *.txt parameter files to xcode/Debug where 'lignum-forest´ is located.
# Otherwise hard coded files are not found in the program.
# Set command line parameters for 'lignum-forest' in Product (menu bar) -> Scheme -> Edit Schmeme -> Arguments.
# Divide the command line into practical parts for debugging from Arguments -> '+'.
# -------------------------------------
cmake_minimum_required(VERSION 3.30)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# cmake 3.23 for Xcode 13 recommends the policy CMP0114 (warning output)
#cmake_policy(SET CMP0114 NEW)
#cmake_policy(SET CMP0175 NEW)
#cmake_policy(SET CMP0177 NEW)
project(LignumForest DESCRIPTION "Individual trees in single species Scots pine forest")
message("Project  ${CMAKE_PROJECT_NAME}: ${LignumForest_DESCRIPTION}")
message("---------------------------------------------------------------------------")
# Lignum core projects as an ExternalProject. 'lignum-core' target is compiled first
# Note: lignum-core  has to be configured with cmake first. See lignum-core/CMakeLists.txt
include(ExternalProject)
# Comment out other L-systems
#Note other L-systems than  pine-em98-branch-C.L may generate "dupclicate symbols" errors
#due to restructuring CrownDensity and LignumForest and rigorous use of using namespaces.
#set(LSYSTEMFILE pine-em98.L)
#set(LSYSTEMSRC  pine-em98.cpp)
#set(LSYSTEMFILE pine-em98-branch-A.L)
#set(LSYSTEMSRC  pine-em98-branch-A.cpp)
#set(LSYSTEMFILE pine-em98-branch-B.L)
#set(LSYSTEMSRC  pine-em98-branch-B.cpp)
set(LSYSTEMFILE pine-em98-branch-C.L)
set(LSYSTEMSRC  pine-em98-branch-C.cpp)
ExternalProject_Add(lignum-core
SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..
DOWNLOAD_COMMAND "")
#The '*.cmake' files for dependencies are preproduced with 'cmake' configuration
#for 'lignum-core' (see lignum-core/CMakeLists.txt).  
include(../c++adt/lib/cxxadt.cmake)
include(../Firmament/lib/sky.cmake)
include(../stl-lignum/lib/LGM.cmake)
include(../stl-voxelspace/lib/voxel.cmake)
include(../LEngine/lib/L.cmake)
include(../LEngine/lib/l2c.cmake)

if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
  message("LignumForest")
  message("------------")
  message("Apple macOS (Darwin) Compiling Visual")
  include(../Graphics/lib/Visual.cmake)
endif()

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
  message("LignumForest")
  message("------------")
  message("Apple macOS (Darwin) platform")
  #The following lines find Qt4 headers and libraries, qmake must be in PATH
  #find_package(Qt4 COMPONENTS QtOpenGL QtXml QtGui QtCore REQUIRED)
  # The following lines find Qt5 headers and libraries, qmake must be in PATH
  find_package(Qt5 COMPONENTS Xml Core  REQUIRED)
  #find_package(Qt6 COMPONENTS Xml Core  REQUIRED)
  #HSF5 for  collecting simulation results
  find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
  #Simple visualization tool in macOS
  find_package(OpenGL REQUIRED) 
  # cmake 3.23 output erroneusly prints 'glut not found'. It is found. 
  find_package(GLUT REQUIRED)
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
  message("LignumForest")
  message("------------")
  message("Linux platform")
  find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
  find_package(Qt5 COMPONENTS Xml Core  REQUIRED)
else()
  message(${CMAKE_HOST_SYSTEM_NAME} " not supported")
  message(FATAL_ERROR "Exiting cmake")
endif()

set(SOURCE_FILES lignum-forest.cc
  src/globalvariables.cc
  ${LSYSTEMSRC}
  branchfunctor.cc
  forest-trees.cc
  generate-tree-locations.cc
  harvest-stand.cc
  branchfunctor.cc
  src/metabolism.cc
  src/borderforest.cc
  src/space.cc
  src/CreateHDF5Files.cc
  src/CreateVoxelSpaceData.cc
)

# Rebuild L system with custom command for 'l2c'
add_custom_command(
OUTPUT ${CMAKE_CURRENT_LIST_DIR}/${LSYSTEMSRC}
DEPENDS ${CMAKE_CURRENT_LIST_DIR}/${LSYSTEMFILE}
COMMAND ${CMAKE_CURRENT_LIST_DIR}/../LEngine/bin/l2c ${CMAKE_CURRENT_LIST_DIR}/${LSYSTEMFILE}  ${CMAKE_CURRENT_LIST_DIR}/${LSYSTEMSRC}
)

# Main target, the executable 'lignum-forest'
add_executable(lignum-forest ${SOURCE_FILES})
target_include_directories(lignum-forest PUBLIC . include
  ${CMAKE_CURRENT_LIST_DIR}/
  ${CMAKE_CURRENT_LIST_DIR}/include
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

if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
  message("LignumForest")
  message("------------")
  message("Apple macOS (Darwin) libraries")
  target_link_libraries(lignum-forest
    #HDF5::HDF5 Qt4::QtOpenGL Qt4::QtXml Qt4::QtGui Qt4::QtCore
    HDF5::HDF5 Qt5::Xml Qt5::Core
    #HDF5::HDF5 Qt6::Xml Qt6::Core
    GLUT::GLUT OpenGL::GLU OpenGL::GL
    Visual L voxel sky LGM cxxadt
  )
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
  message("LignumForest")
  message("------------")
  message("Linux libraries")
  target_link_libraries(lignum-forest
    HDF5::HDF5 Qt5::Xml Qt5::Core
    L voxel sky LGM cxxadt
  )
else()
  message(${CMAKE_HOST_SYSTEM_NAME} " not supported")
  message(FATAL_ERROR "Exiting cmake")
endif()
add_custom_command(TARGET lignum-forest
  COMMAND echo "Compilation for LignumForest done"
  POST_BUILD
  COMMENT "REMEMBER to type \"make install\" to install \"lignum-forest\" into the LignumForest working directory above."
  )
install(TARGETS lignum-forest DESTINATION "${CMAKE_CURRENT_LIST_DIR}")

