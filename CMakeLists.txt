# $Id: CMakeLists.txt 68058 2013-03-13 14:47:43Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4Medipix)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "-g -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11")


if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

# Build multithreaded
set(G4MULTITHREADED ON)



set(useHDF5 true)        # use true or false


#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()


       
# include boost
find_package(Boost COMPONENTS system filesystem REQUIRED)

# CADMesh
# find_package(cadmesh)

# Include the HDF5 library

find_package(HDF5 COMPONENTS CXX)
if(HDF5_FOUND)
  option(WITH_HDF5 "HDF5 output enabled" ${HDF5_FOUND})
  add_definitions(-DWITH_HDF5)
  include_directories(${HDF5_CXX_INCLUDE_DIR})
  #add_definitions(${HDF5_hdf5_LIBRARY})
  set(hdf5_libs hdf5_serial hdf5_cpp hdf5_hl_cpp)
#  set(LIBS ${LIBS} ${HDF5_hdf5_cpp_LIBRARY})
endif()

#find_package(HDF5 COMPONENTS CXX)
#include_directories(${HDF5_CXX_INCLUDE_DIR})



#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(G4Medipix G4Medipix.cc ${sources} ${headers} include/OutputMessenger.hh)
if(HDF5_FOUND)
  target_link_libraries(G4Medipix ${HDF5_hdf5_cpp_LIBRARY})
  target_link_libraries(G4Medipix ${HDF5_hdf5_LIBRARY})
endif()
target_link_libraries(G4Medipix ${Geant4_LIBRARIES})
target_link_libraries(G4Medipix ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries(G4Medipix ${Boost_SYSTEM_LIBRARY})
#target_link_libraries(G4Medipix ${cadmesh_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4c. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(G4MEDIPIX_SCRIPTS
  icons.mac
  gui.mac
  init.mac
  init_vis.mac
  vis.mac
#  mpx.mac
#  mpx_vis.mac
#  dna.mac
  am241.mac
  alpha.mac
  dacscan.mac
#  test_erik.mac
  gamma.mac
  proton.mac
  electron.mac
  plotHisto.C
  gamma_calibration.mac
  fluorescence.mac
  )

foreach(_script ${G4MEDIPIX_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/macros/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

configure_file(
    ${PROJECT_SOURCE_DIR}/run.png
    ${PROJECT_BINARY_DIR}/run.png
    COPYONLY
)

include(${PROJECT_SOURCE_DIR}/tubes_list.txt)

#----------------------------------------------------------------------------
# copy the potential maps as well, if they don't exist already
set(POTENTIAL_MAPS
  3d_regular_grid_weighting_field110_300_fine.bin
  3d_regular_grid_weighting_field110_500_fine.bin
  3d_regular_grid_weighting_field110_1000_fine.bin
  3d_regular_grid_weighting_field55_300_fine.bin
  3d_regular_grid_weighting_field55_500_fine.bin
  3d_regular_grid_weighting_field55_1000_fine.bin
  )

foreach(_script ${POTENTIAL_MAPS})
  get_filename_component(name ${_script} NAME)
  if (NOT EXISTS ${PROJECT_BINARY_DIR}/potential/${_script} )
    configure_file(
      ${PROJECT_SOURCE_DIR}/potential/${_script}
      ${PROJECT_BINARY_DIR}/potential/${_script}
      COPYONLY
      )
  endif (NOT EXISTS ${PROJECT_BINARY_DIR}/potential/${_script} )
endforeach()

#----------------------------------------------------------------------------
# copy other files as .ini
set(FilesToCopy
  DetectorConfig.ini
  )

foreach(_script ${FilesToCopy})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
#install(TARGETS G4Medipix DESTINATION bin)
