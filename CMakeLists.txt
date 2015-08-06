cmake_minimum_required(VERSION 2.8)
project(magudi Fortran C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
string(TOLOWER ${CMAKE_BUILD_TYPE} cmake_build_type_tolower)
if(NOT cmake_build_type_tolower STREQUAL "debug" AND NOT cmake_build_type_tolower STREQUAL "release")
  message(FATAL_ERROR "Unknown build type ${CMAKE_BUILD_TYPE}. Allowed values are: Debug, Release (case-insensitive).")
endif()

set(magudi_SCALAR_TYPES real32_iso real64_iso real128_iso binary128_IEEE754)
if(NOT SCALAR_TYPE)
  set(SCALAR_TYPE real64_iso)
endif()
list(FIND magudi_SCALAR_TYPES ${SCALAR_TYPE} SCALAR_TYPE_IS_VALID)
if(SCALAR_TYPE_IS_VALID EQUAL -1)
  message(FATAL_ERROR "Unknown scalar type ${SCALAR_TYPE}. Allowed values are: ${magudi_SCALAR_TYPES} (case-sensitive).")
endif()
add_definitions(-DSCALAR_TYPE_IS_${SCALAR_TYPE})

option(USE_EXTENDED_PLOT3D "Use an extended PLOT3D format with 64-bit integer byte counts." OFF)
if(USE_EXTENDED_PLOT3D)
  add_definitions(-DUSE_EXTENDED_PLOT3D)
endif()

find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})

set(CMAKE_Fortran_FORMAT "FREE")
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

configure_file(${PROJECT_SOURCE_DIR}/include/config.h.in ${PROJECT_BINARY_DIR}/include/config.h)
include_directories(${PROJECT_BINARY_DIR}/include)

set(magudi_SINGLETONS
  
  ${PROJECT_SOURCE_DIR}/src/RandomNumber.f90
  ${PROJECT_SOURCE_DIR}/src/ErrorHandler.f90
  ${PROJECT_SOURCE_DIR}/src/MPIHelper.f90
  ${PROJECT_SOURCE_DIR}/src/MPITimingsHelper.f90
  ${PROJECT_SOURCE_DIR}/src/PLOT3DHelper.f90
  ${PROJECT_SOURCE_DIR}/src/CNSHelper.f90
  ${PROJECT_SOURCE_DIR}/src/InputHelper.f90
  ${PROJECT_SOURCE_DIR}/src/Solver.f90
  ${PROJECT_SOURCE_DIR}/src/RhsHelper.f90
  
  )

set(magudi_FACTORIES
  
  ${PROJECT_SOURCE_DIR}/src/PatchFactory.f90
  ${PROJECT_SOURCE_DIR}/src/ControllerFactory.f90
  ${PROJECT_SOURCE_DIR}/src/FunctionalFactory.f90
  ${PROJECT_SOURCE_DIR}/src/TimeIntegratorFactory.f90
  ${PROJECT_SOURCE_DIR}/src/ReverseMigratorFactory.f90
  
  )

set(magudi_OBJECTS

  ${PROJECT_SOURCE_DIR}/src/SimulationFlags.f90
  ${PROJECT_SOURCE_DIR}/src/SolverOptions.f90
  ${PROJECT_SOURCE_DIR}/src/StencilOperator.f90
  ${PROJECT_SOURCE_DIR}/src/Grid.f90
  ${PROJECT_SOURCE_DIR}/src/State.f90
  ${PROJECT_SOURCE_DIR}/src/Region.f90
  ${PROJECT_SOURCE_DIR}/src/AcousticSource.f90
  ${PROJECT_SOURCE_DIR}/src/ResidualManager.f90

  )

set(magudi_ABSTRACT_TYPES

  ${PROJECT_SOURCE_DIR}/src/Patch.f90
  ${PROJECT_SOURCE_DIR}/src/Controller.f90
  ${PROJECT_SOURCE_DIR}/src/Functional.f90
  ${PROJECT_SOURCE_DIR}/src/TimeIntegrator.f90
  ${PROJECT_SOURCE_DIR}/src/ReverseMigrator.f90

  )

set(magudi_EXTENDED_TYPES
  
  ${PROJECT_SOURCE_DIR}/src/RK4Integrator.f90
  ${PROJECT_SOURCE_DIR}/src/SpongePatch.f90
  ${PROJECT_SOURCE_DIR}/src/InflowOutflowPatch.f90
  ${PROJECT_SOURCE_DIR}/src/ImpenetrableWall.f90
  ${PROJECT_SOURCE_DIR}/src/IsothermalWall.f90
  ${PROJECT_SOURCE_DIR}/src/CostTargetPatch.f90
  ${PROJECT_SOURCE_DIR}/src/ActuatorPatch.f90
  ${PROJECT_SOURCE_DIR}/src/BlockInterfacePatch.f90
  ${PROJECT_SOURCE_DIR}/src/UniformCheckpointing.f90
  ${PROJECT_SOURCE_DIR}/src/ThermalActuator.f90
  ${PROJECT_SOURCE_DIR}/src/AcousticNoise.f90
  ${PROJECT_SOURCE_DIR}/src/ProbePatch.f90
  ${PROJECT_SOURCE_DIR}/src/SolenoidalExcitationPatch.f90

  )

set(magudiObj_FORTRAN_SOURCES ${magudi_SINGLETONS} ${magudi_FACTORIES} ${magudi_OBJECTS} ${magudi_ABSTRACT_TYPES} ${magudi_EXTENDED_TYPES})
set(magudiObj_C_SOURCES ${PROJECT_SOURCE_DIR}/src/PLOT3DFormat.c)
set(magudiObj_SOURCES ${magudiObj_FORTRAN_SOURCES} ${magudiObj_C_SOURCES})

add_library(magudiObj OBJECT ${magudiObj_SOURCES})
add_executable(magudi $<TARGET_OBJECTS:magudiObj> ${PROJECT_SOURCE_DIR}/src/main.f90)

target_link_libraries(magudi ${MPI_Fortran_LIBRARIES})

if(${CMAKE_COMPILER_IS_GNUG77})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fmodule-private -fimplicit-none -pedantic-errors -std=f2008")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fbacktrace -Wall -Wextra -Waliasing -ffree-line-length-none -fall-intrinsics -fcheck=all")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -ftree-vectorize -ffast-math -funroll-loops")
  if(cmake_build_type_tolower STREQUAL "debug")
    set_source_files_properties(${magudi_EXTENDED_TYPES} PROPERTIES COMPILE_FLAGS "-Wno-unused-dummy-argument -Wno-unused-parameter")
  endif()
else()
  get_filename_component(CMAKE_Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
  message(FATAL_ERROR "${CMAKE_Fortran_COMPILER_NAME}: Fortran compiler is not supported.")
endif()

set(magudi_BUILD_FILES

  # CMake-generated files
  ${PROJECT_BINARY_DIR}/CMakeCache.txt
  ${PROJECT_BINARY_DIR}/CMakeFiles
  ${PROJECT_BINARY_DIR}/cmake_install.cmake
  ${PROJECT_BINARY_DIR}/Makefile

  # Build directories
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_BINARY_DIR}/modules

  )

set(magudi_BUILD_FILES_STRING "")
foreach(file ${magudi_BUILD_FILES})
  file(TO_NATIVE_PATH ${file} file)
  set(magudi_BUILD_FILES_STRING "${magudi_BUILD_FILES_STRING} ${file}")
endforeach()
add_custom_target(purge make clean COMMAND echo ${magudi_BUILD_FILES_STRING} | xargs rm -rf VERBATIM)
