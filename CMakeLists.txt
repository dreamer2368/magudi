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
execute_process(COMMAND mpif90 -show COMMAND awk "{print $1}" OUTPUT_VARIABLE CMAKE_Fortran_COMPILER_PATH)
execute_process(COMMAND basename ${CMAKE_Fortran_COMPILER_PATH} OUTPUT_VARIABLE CMAKE_Fortran_COMPILER_NAME OUTPUT_STRIP_TRAILING_WHITESPACE)

if(CMAKE_Fortran_COMPILER_NAME STREQUAL "gfortran")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -fmodule-private -fimplicit-none -cpp -pedantic-errors -std=f2008")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fbacktrace -Wall -Wextra -Waliasing -Wno-unused-dummy-argument -Wno-unused-parameter -ffree-line-length-none -fall-intrinsics -fcheck=all -g")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ftree-vectorize -ffast-math -funroll-loops")
elseif(CMAKE_Fortran_COMPILER_NAME STREQUAL "ifort")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -FR -u -fpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-traceback -gen-interfaces -check all -debug all -debug-parameters all -ftrapuv -fpe0 -warn all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -parallel -no-prec-div -ipo -gen-interfaces")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -xHost -parallel -no-prec-div -ipo")
else()
  message(FATAL_ERROR "${CMAKE_Fortran_COMPILER_NAME}: Fortran compiler could not be determined correctly or is unsupported.")
endif()
set(CMAKE_Fortran_FORMAT "FREE")
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

set(CMAKE_C_FLAGS_DEBUG "-std=c89 -pedantic-errors")

configure_file(${PROJECT_SOURCE_DIR}/include/config.h.in ${PROJECT_BINARY_DIR}/include/config.h)
include_directories(${PROJECT_BINARY_DIR}/include)

set(magudiObj_SOURCES

  # Singletons.
  ${PROJECT_SOURCE_DIR}/include/RandomNumber.f90
  ${PROJECT_SOURCE_DIR}/src/RandomNumberImpl.f90
  ${PROJECT_SOURCE_DIR}/include/ErrorHandler.f90
  ${PROJECT_SOURCE_DIR}/src/ErrorHandlerImpl.f90
  ${PROJECT_SOURCE_DIR}/include/ThermoChemistry.f90
  ${PROJECT_SOURCE_DIR}/src/ThermoChemistryImpl.f90
  ${PROJECT_SOURCE_DIR}/include/MPIHelper.f90
  ${PROJECT_SOURCE_DIR}/src/MPIHelperImpl.f90
  ${PROJECT_SOURCE_DIR}/include/MPITimingsHelper.f90
  ${PROJECT_SOURCE_DIR}/src/MPITimingsHelperImpl.f90
  ${PROJECT_SOURCE_DIR}/include/PLOT3DHelper.f90
  ${PROJECT_SOURCE_DIR}/src/PLOT3DHelperImpl.f90
  ${PROJECT_SOURCE_DIR}/include/InputHelper.f90
  ${PROJECT_SOURCE_DIR}/src/InputHelperImpl.f90
  ${PROJECT_SOURCE_DIR}/include/Solver.f90
  ${PROJECT_SOURCE_DIR}/src/SolverImpl.f90
  ${PROJECT_SOURCE_DIR}/include/CNSHelper.f90
  ${PROJECT_SOURCE_DIR}/src/CNSHelperImpl.f90
  ${PROJECT_SOURCE_DIR}/include/RhsHelper.f90
  ${PROJECT_SOURCE_DIR}/src/RhsHelperImpl.f90
  ${PROJECT_SOURCE_DIR}/include/InterfaceHelper.f90
  ${PROJECT_SOURCE_DIR}/src/InterfaceHelperImpl.f90
  ${PROJECT_SOURCE_DIR}/include/EnsightHelper.f90
  ${PROJECT_SOURCE_DIR}/src/EnsightHelperImpl.f90

  # Factories.
  ${PROJECT_SOURCE_DIR}/include/PatchFactory.f90
  ${PROJECT_SOURCE_DIR}/src/PatchFactoryImpl.f90
  ${PROJECT_SOURCE_DIR}/include/ControllerFactory.f90
  ${PROJECT_SOURCE_DIR}/src/ControllerFactoryImpl.f90
  ${PROJECT_SOURCE_DIR}/include/FunctionalFactory.f90
  ${PROJECT_SOURCE_DIR}/src/FunctionalFactoryImpl.f90
  ${PROJECT_SOURCE_DIR}/include/TimeIntegratorFactory.f90
  ${PROJECT_SOURCE_DIR}/src/TimeIntegratorFactoryImpl.f90
  ${PROJECT_SOURCE_DIR}/include/ReverseMigratorFactory.f90
  ${PROJECT_SOURCE_DIR}/src/ReverseMigratorFactoryImpl.f90
  ${PROJECT_SOURCE_DIR}/include/MappingFunctionFactory.f90
  ${PROJECT_SOURCE_DIR}/src/MappingFunctionFactoryImpl.f90

  # Derived types.
  ${PROJECT_SOURCE_DIR}/include/SimulationFlags.f90
  ${PROJECT_SOURCE_DIR}/src/SimulationFlagsImpl.f90
  ${PROJECT_SOURCE_DIR}/include/SolverOptions.f90
  ${PROJECT_SOURCE_DIR}/src/SolverOptionsImpl.f90
  ${PROJECT_SOURCE_DIR}/include/StencilOperator.f90
  ${PROJECT_SOURCE_DIR}/src/StencilOperatorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/Grid.f90
  ${PROJECT_SOURCE_DIR}/src/GridImpl.f90
  ${PROJECT_SOURCE_DIR}/include/Particles.f90
  ${PROJECT_SOURCE_DIR}/src/ParticlesImpl.f90
  ${PROJECT_SOURCE_DIR}/include/State.f90
  ${PROJECT_SOURCE_DIR}/src/StateImpl.f90
  ${PROJECT_SOURCE_DIR}/include/PatchDescriptor.f90
  ${PROJECT_SOURCE_DIR}/src/PatchDescriptorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/Region.f90
  ${PROJECT_SOURCE_DIR}/src/RegionImpl.f90
  ${PROJECT_SOURCE_DIR}/include/AcousticSource.f90
  ${PROJECT_SOURCE_DIR}/src/AcousticSourceImpl.f90
  ${PROJECT_SOURCE_DIR}/include/IgnitionSource.f90
  ${PROJECT_SOURCE_DIR}/src/IgnitionSourceImpl.f90
  ${PROJECT_SOURCE_DIR}/include/FuelSource.f90
  ${PROJECT_SOURCE_DIR}/src/FuelSourceImpl.f90
  ${PROJECT_SOURCE_DIR}/include/Combustion.f90
  ${PROJECT_SOURCE_DIR}/src/CombustionImpl.f90
  ${PROJECT_SOURCE_DIR}/include/SolenoidalExcitationPatch.f90
  ${PROJECT_SOURCE_DIR}/src/SolenoidalExcitationPatchImpl.f90
  ${PROJECT_SOURCE_DIR}/include/JetExcitationPatch.f90
  ${PROJECT_SOURCE_DIR}/src/JetExcitationPatchImpl.f90

  # Abstract types.
  ${PROJECT_SOURCE_DIR}/include/Patch.f90
  ${PROJECT_SOURCE_DIR}/src/PatchImpl.f90
  ${PROJECT_SOURCE_DIR}/include/Controller.f90
  ${PROJECT_SOURCE_DIR}/src/ControllerImpl.f90
  ${PROJECT_SOURCE_DIR}/include/Functional.f90
  ${PROJECT_SOURCE_DIR}/src/FunctionalImpl.f90
  ${PROJECT_SOURCE_DIR}/include/TimeIntegrator.f90
  ${PROJECT_SOURCE_DIR}/src/TimeIntegratorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/ReverseMigrator.f90
  ${PROJECT_SOURCE_DIR}/src/ReverseMigratorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/MappingFunction.f90

  # Extended types.
  ${PROJECT_SOURCE_DIR}/include/RK4Integrator.f90
  ${PROJECT_SOURCE_DIR}/src/RK4IntegratorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/JamesonRK3Integrator.f90
  ${PROJECT_SOURCE_DIR}/src/JamesonRK3IntegratorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/SpongePatch.f90
  ${PROJECT_SOURCE_DIR}/src/SpongePatchImpl.f90
  ${PROJECT_SOURCE_DIR}/include/FarFieldPatch.f90
  ${PROJECT_SOURCE_DIR}/src/FarFieldPatchImpl.f90
  ${PROJECT_SOURCE_DIR}/include/ImpenetrableWall.f90
  ${PROJECT_SOURCE_DIR}/src/ImpenetrableWallImpl.f90
  ${PROJECT_SOURCE_DIR}/include/IsothermalWall.f90
  ${PROJECT_SOURCE_DIR}/src/IsothermalWallImpl.f90
  ${PROJECT_SOURCE_DIR}/include/AdiabaticWall.f90
  ${PROJECT_SOURCE_DIR}/src/AdiabaticWallImpl.f90
  ${PROJECT_SOURCE_DIR}/include/CostTargetPatch.f90
  ${PROJECT_SOURCE_DIR}/src/CostTargetPatchImpl.f90
  ${PROJECT_SOURCE_DIR}/include/ActuatorPatch.f90
  ${PROJECT_SOURCE_DIR}/src/ActuatorPatchImpl.f90
  ${PROJECT_SOURCE_DIR}/include/BlockInterfacePatch.f90
  ${PROJECT_SOURCE_DIR}/src/BlockInterfacePatchImpl.f90
  ${PROJECT_SOURCE_DIR}/include/UniformCheckpointer.f90
  ${PROJECT_SOURCE_DIR}/src/UniformCheckpointerImpl.f90
  ${PROJECT_SOURCE_DIR}/include/ThermalActuator.f90
  ${PROJECT_SOURCE_DIR}/src/ThermalActuatorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/FuelActuator.f90
  ${PROJECT_SOURCE_DIR}/src/FuelActuatorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/IgnitionActuator.f90
  ${PROJECT_SOURCE_DIR}/src/IgnitionActuatorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/MomentumActuator.f90
  ${PROJECT_SOURCE_DIR}/src/MomentumActuatorImpl.f90
  ${PROJECT_SOURCE_DIR}/include/AcousticNoise.f90
  ${PROJECT_SOURCE_DIR}/src/AcousticNoiseImpl.f90
  ${PROJECT_SOURCE_DIR}/include/PressureDrag.f90
  ${PROJECT_SOURCE_DIR}/src/PressureDragImpl.f90
  ${PROJECT_SOURCE_DIR}/include/DragForce.f90
  ${PROJECT_SOURCE_DIR}/src/DragForceImpl.f90
  ${PROJECT_SOURCE_DIR}/include/ReynoldsStress.f90
  ${PROJECT_SOURCE_DIR}/src/ReynoldsStressImpl.f90
  ${PROJECT_SOURCE_DIR}/include/FlameTemperature.f90
  ${PROJECT_SOURCE_DIR}/src/FlameTemperatureImpl.f90
  ${PROJECT_SOURCE_DIR}/include/ReactantDepletion.f90
  ${PROJECT_SOURCE_DIR}/src/ReactantDepletionImpl.f90
  ${PROJECT_SOURCE_DIR}/include/UniformMap.f90
  ${PROJECT_SOURCE_DIR}/src/UniformMapImpl.f90

  # C source files.
  ${PROJECT_SOURCE_DIR}/src/PLOT3DFormat.c

  )

add_library(magudiObj OBJECT ${magudiObj_SOURCES})
add_executable(magudi $<TARGET_OBJECTS:magudiObj> ${PROJECT_SOURCE_DIR}/src/main.f90)
enable_testing()
add_subdirectory(test)
add_subdirectory(utils)

set(magudi_BUILD_FILES

  # CMake-generated files
  ${PROJECT_BINARY_DIR}/CMakeCache.txt
  ${PROJECT_BINARY_DIR}/CMakeFiles
  ${PROJECT_BINARY_DIR}/cmake_install.cmake
  ${PROJECT_BINARY_DIR}/CTestTestfile.cmake
  ${PROJECT_BINARY_DIR}/Makefile

  # Compiled Python sources
  ${PROJECT_BINARY_DIR}/PLOT3D.pyc
  ${PROJECT_BINARY_DIR}/MatplotlibHelper.pyc
  ${PROJECT_BINARY_DIR}/SummationByParts.pyc

  # Build directories
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_BINARY_DIR}/modules
  ${PROJECT_BINARY_DIR}/test
  ${PROJECT_BINARY_DIR}/Testing
  ${PROJECT_BINARY_DIR}/utils

  )

set(magudi_BUILD_FILES_STRING "")
foreach(file ${magudi_BUILD_FILES})
  file(TO_NATIVE_PATH ${file} file)
  set(magudi_BUILD_FILES_STRING "${magudi_BUILD_FILES_STRING} ${file}")
endforeach()
add_custom_target(purge make clean COMMAND echo ${magudi_BUILD_FILES_STRING} | xargs rm -rf VERBATIM)
