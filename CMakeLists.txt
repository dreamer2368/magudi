cmake_minimum_required(VERSION 2.8)
project(magudi Fortran C)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if (NOT cmake_build_type_tolower STREQUAL "debug" AND NOT cmake_build_type_tolower STREQUAL "release" AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()

if (cmake_build_type_tolower STREQUAL "debug")
  add_definitions(-DDEBUG)
endif()

if (NOT SCALAR_TYPE)
  set(SCALAR_TYPE "real64_iso")
endif()

if (NOT SCALAR_TYPE STREQUAL "real32_iso" AND 
    NOT SCALAR_TYPE STREQUAL "real64_iso" AND 
    NOT SCALAR_TYPE STREQUAL "real128_iso" AND
    NOT SCALAR_TYPE STREQUAL "binary128_IEEE754" AND
    NOT SCALAR_TYPE STREQUAL "complex" AND 
    NOT SCALAR_TYPE STREQUAL "double_complex")
  message(FATAL_ERROR "Unknown scalar type \"${SCALAR_TYPE}\". Allowed values are: real32_iso, real64_iso, real128_iso, binary128_IEEE754, complex, double_complex (case-sensitive).")
endif()

set(SCALAR_TYPE_MACRO "SCALAR_TYPE_IS_${SCALAR_TYPE}")

if (NOT PLOT3D_FORMAT)
  set(PLOT3D_FORMAT "regular")
endif()

if (NOT PLOT3D_FORMAT STREQUAL "regular" AND 
    NOT PLOT3D_FORMAT STREQUAL "extended")
  message(FATAL_ERROR "Unknown PLOT3D format \"${PLOT3D_FORMAT}\". Allowed values are: regular, extended (case-sensitive).")
endif()

if (PLOT3D_FORMAT STREQUAL "extended")
  add_definitions(-DUSE_EXTENDED_PLOT3D)
endif()

find_program(MPIF90_FOUND mpif90)
if(NOT MPIF90_FOUND)
  message(FATAL_ERROR "mpif90 not found.")
endif()
set(CMAKE_Fortran_COMPILER "mpif90")
execute_process(COMMAND mpif90 -show COMMAND sed "s/^\\([^[:space:]]*\\).*/\\1/g" OUTPUT_VARIABLE CMAKE_Fortran_COMPILER_PATH)
execute_process(COMMAND basename ${CMAKE_Fortran_COMPILER_PATH} OUTPUT_VARIABLE CMAKE_Fortran_COMPILER_NAME OUTPUT_STRIP_TRAILING_WHITESPACE)

if(CMAKE_Fortran_COMPILER_NAME STREQUAL "gfortran")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -fmodule-private -fimplicit-none -cpp -pedantic-errors -std=f2008")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fbacktrace -Wall -Wextra -Waliasing -ffree-line-length-none -fall-intrinsics -fcheck=all -g")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ftree-vectorize -ffast-math -funroll-loops")
elseif(CMAKE_Fortran_COMPILER_NAME STREQUAL "ifort")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -FR -u -fpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-traceback -gen-interfaces -check all -debug all -debug-parameters all -ftrapuv -fpe0 -warn all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -parallel -no-prec-div -ipo -gen-interfaces")
else()
  message(FATAL_ERROR "${CMAKE_Fortran_COMPILER_NAME}: Fortran compiler could not be determined correctly or is unsupported.")
endif()
set(CMAKE_C_FLAGS_DEBUG "-std=c89 -pedantic-errors")

set(CMAKE_Fortran_FORMAT "FREE")
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
configure_file(${PROJECT_SOURCE_DIR}/include/config.h.in ${PROJECT_BINARY_DIR}/include/config.h)
include_directories(${PROJECT_BINARY_DIR}/include)

set(magudiObj_SOURCES

  # Singleton classes.
  ${PROJECT_SOURCE_DIR}/include/RandomNumber.f90
  ${PROJECT_SOURCE_DIR}/src/RandomNumberImpl.f90
  ${PROJECT_SOURCE_DIR}/include/ErrorHandler.f90
  ${PROJECT_SOURCE_DIR}/src/ErrorHandlerImpl.f90
  ${PROJECT_SOURCE_DIR}/include/MPIHelper.f90
  ${PROJECT_SOURCE_DIR}/src/MPIHelperImpl.f90
  ${PROJECT_SOURCE_DIR}/include/MPITimingsHelper.f90
  ${PROJECT_SOURCE_DIR}/src/MPITimingsHelperImpl.f90
  ${PROJECT_SOURCE_DIR}/include/PLOT3DHelper.f90
  ${PROJECT_SOURCE_DIR}/src/PLOT3DHelperImpl.f90
  ${PROJECT_SOURCE_DIR}/include/CNSHelper.f90
  ${PROJECT_SOURCE_DIR}/src/CNSHelperImpl.f90
  ${PROJECT_SOURCE_DIR}/include/InputHelper.f90
  ${PROJECT_SOURCE_DIR}/src/InputHelperImpl.f90
  ${PROJECT_SOURCE_DIR}/include/Solver.f90
  ${PROJECT_SOURCE_DIR}/src/SolverImpl.f90

  # Object-oriented classes.
  ${PROJECT_SOURCE_DIR}/include/SimulationFlags.f90
  ${PROJECT_SOURCE_DIR}/src/SimulationFlagsImpl.f90
  ${PROJECT_SOURCE_DIR}/include/SolverOptions.f90
  ${PROJECT_SOURCE_DIR}/src/SolverOptionsImpl.f90
  ${PROJECT_SOURCE_DIR}/include/StencilOperator.f90
  ${PROJECT_SOURCE_DIR}/src/StencilOperatorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/Grid.f90
  ${PROJECT_SOURCE_DIR}/src/GridImpl.f90
  ${PROJECT_SOURCE_DIR}/include/State.f90
  ${PROJECT_SOURCE_DIR}/src/StateImpl.f90
  ${PROJECT_SOURCE_DIR}/include/Patch.f90
  ${PROJECT_SOURCE_DIR}/src/PatchImpl.f90
  ${PROJECT_SOURCE_DIR}/include/PatchDescriptor.f90
  ${PROJECT_SOURCE_DIR}/src/PatchDescriptorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/Region.f90
  ${PROJECT_SOURCE_DIR}/src/RegionImpl.f90
  ${PROJECT_SOURCE_DIR}/include/AcousticSource.f90
  ${PROJECT_SOURCE_DIR}/src/AcousticSourceImpl.f90
  ${PROJECT_SOURCE_DIR}/include/SolenoidalExcitation.f90
  ${PROJECT_SOURCE_DIR}/src/SolenoidalExcitationImpl.f90
  ${PROJECT_SOURCE_DIR}/include/TimeIntegrator.f90
  ${PROJECT_SOURCE_DIR}/src/TimeIntegratorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/RK4Integrator.f90
  ${PROJECT_SOURCE_DIR}/src/RK4IntegratorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/ReverseMigrator.f90
  ${PROJECT_SOURCE_DIR}/src/ReverseMigratorImpl.f90

  # C source files.
  ${PROJECT_SOURCE_DIR}/src/PLOT3DFormat.c

  )

add_library(magudiObj OBJECT ${magudiObj_SOURCES})
add_executable(magudi $<TARGET_OBJECTS:magudiObj> ${PROJECT_SOURCE_DIR}/src/main.f90)
enable_testing()
add_subdirectory(test)
add_subdirectory(utils)
