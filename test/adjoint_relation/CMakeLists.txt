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
#add_serial_test(inviscid_flux)

# MPI tests
set(numProcs 1 2 4 8)
add_mpi_test(artificial_dissipation "${numProcs}")
add_mpi_test(inviscid_flux "${numProcs}")
add_mpi_test(viscous_flux "${numProcs}")
add_mpi_test(SBP_operators "${numProcs}")
