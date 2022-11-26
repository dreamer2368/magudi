function(add_mpi_test name numProcs)
  add_executable(${name} $<TARGET_OBJECTS:magudiObj> ${name}.f90)
  foreach(np ${numProcs})
    set(test_parameters -np ${np} "./${name}")
    add_test(NAME ${name}_on_${np}_tasks COMMAND "mpirun" ${test_parameters})
  endforeach()
endfunction(add_mpi_test)

function(add_serial_test name)
  add_executable(${name} $<TARGET_OBJECTS:magudiObj> ${name}.f90)
  add_test(NAME ${name} COMMAND ${name})
endfunction(add_serial_test)

function(add_executable_simple name)
  add_executable(${name} $<TARGET_OBJECTS:magudiObj> ${name}.f90)
endfunction(add_executable_simple)

function(add_script_test name)
    add_executable(${name} $<TARGET_OBJECTS:magudiObj> ${name}.f90)

    add_test(NAME ${name}
             COMMAND sh ${name}.sh ${CMAKE_CURRENT_BINARY_DIR}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()

# Serial tests
add_script_test(SAT_block_interface_linearized)
add_script_test(full_rhs_linearized)

# MPI tests
set(numProcs 1 2)
# add_mpi_test(sponge_linearized "${numProcs}")
add_serial_test(sponge_linearized)
#add_mpi_test(SAT_farfield_linearized "${numProcs}")
add_serial_test(SAT_farfield_linearized)
#add_mpi_test(SAT_isothermal_linearized "${numProcs}")
add_serial_test(SAT_isothermal_linearized)
#add_mpi_test(artificial_dissipation_linearized "${numProcs}")
add_serial_test(artificial_dissipation_linearized)
#add_mpi_test(inviscid_flux_linearized "${numProcs}")
add_serial_test(inviscid_flux_linearized)
#add_mpi_test(viscous_flux_linearized "${numProcs}")
add_serial_test(viscous_flux_linearized)
#add_mpi_test(viscous_flux_C "${numProcs}")
add_serial_test(viscous_flux_C)
#add_mpi_test(body_force_linearized "${numProcs}")
add_serial_test(body_force_linearized)
