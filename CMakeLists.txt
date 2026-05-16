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
project(LignumForest DESCRIPTION "Individual trees in single species Scots pine forest" LANGUAGES CXX)
message("Project  ${CMAKE_PROJECT_NAME}: ${LignumForest_DESCRIPTION}")
message("---------------------------------------------------------------------------")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
include(CMakePrintHelpers)

#Using "add_subdirectory" aligns with modern CMake as it supports namespaced targets
#and there is no need to use explicit directory paths.
#Also, the maintenance is easier - to refer to "lignum-core" relative to LignumForest
#change the first argument "../" to point to "lignum-core" CMakeLists.txt location.
#"lignum-core_build" is the local build directory for lignum-core,
#located inside the LignumForest build directory.
add_subdirectory(../ lignum-core_build)

# Choose one L-system and comment out other L-systems
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

cmake_print_variables(CMAKE_CXX_COMPILER)
cmake_print_variables(CMAKE_CXX_COMPILER_ID) 
cmake_print_variables(CMAKE_CXX_STANDARD)

#HDF5 is mandatory
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)

#Select one of Qt4, Qt5 or Qt6 for XML parsing.
#----------------------------------------------
#Set also libray linking accordingly (target_link_libray).
#Qt4 is obsolete and only for macOS with Intel processors
#find_package(Qt4 COMPONENTS QtOpenGL QtXml QtGui QtCore REQUIRED)
#find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
#Qt5 package
find_package(Qt5 COMPONENTS Xml Core  REQUIRED)
#Qt6 package
#find_package(Qt6 COMPONENTS Xml Core  REQUIRED)
if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
  message("LignumForest")
  message("Apple macOS (Darwin) platform")
  message("-----------------------------")
  #Simple visualization tool in macOS
  find_package(OpenGL REQUIRED) 
  # cmake 3.x output erroneusly prints 'glut not found'. It is found. 
  find_package(GLUT REQUIRED)
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
  message("LignumForest")
  message("Linux platform")
  message("--------------")
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

# Rebuild L system with custom command 'l2c'
add_custom_command(
OUTPUT ${CMAKE_CURRENT_LIST_DIR}/${LSYSTEMSRC}
DEPENDS ${CMAKE_CURRENT_LIST_DIR}/${LSYSTEMFILE}
COMMAND l2c::l2c  ${CMAKE_CURRENT_LIST_DIR}/${LSYSTEMFILE}  ${CMAKE_CURRENT_LIST_DIR}/${LSYSTEMSRC}
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

if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
  message("LignumForest")
  message("------------")
  message("Apple macOS (Darwin) libraries")
  target_link_libraries(lignum-forest
    HDF5::HDF5
    #Obsolete Qt4 
    #Qt4::QtOpenGL Qt4::QtXml Qt4::QtGui Qt4::QtCore
    #Qt5 linking
    Qt5::Xml Qt5::Core
    #Qt6 linking
    #Qt6::Xml Qt6::Core
    GLUT::GLUT OpenGL::GLU OpenGL::GL
    Visual::Visual L::L voxel::voxel sky::sky LGM::LGM cxxadt::cxxadt
  )
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
  message("LignumForest")
  message("------------")
  message("Linux libraries")
  target_link_libraries(lignum-forest
    HDF5::HDF5
    #Obsolete Qt4 
    #Qt4::QtOpenGL Qt4::QtXml Qt4::QtGui Qt4::QtCore
    #Qt5 linking
    Qt5::Xml Qt5::Core
    #Qt6 linking
    #Qt6::Xml Qt6::Core
    L::L voxel::voxel sky::sky LGM::LGM cxxadt::cxxadt
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

