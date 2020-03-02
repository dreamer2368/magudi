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

# Serial tests
# add_serial_test(SAT_block_interface)
# add_serial_test(control_functional)
# add_serial_test(initial_condition)

# MPI tests
set(numProcs 1 2 4 8)
# add_mpi_test(time_integrator "${numProcs}")
# add_mpi_test(full_rhs "${numProcs}")
# add_mpi_test(sponge "${numProcs}")
# add_mpi_test(SAT_farfield "${numProcs}")
# add_mpi_test(artificial_dissipation "${numProcs}")
add_mpi_test(inviscid_flux_linearized "${numProcs}")
# add_mpi_test(viscous_flux "${numProcs}")
# add_mpi_test(SBP_operators "${numProcs}")
# add_mpi_test(body_force "${numProcs}")
