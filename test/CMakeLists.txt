function(add_mpi_test name numprocs)
  add_executable(${name} $<TARGET_OBJECTS:magudiObj> ${name}.f90)
  foreach(np ${numprocs})
    set(test_parameters -np ${np} "./${name}")
    add_test(NAME ${name}.np${np} COMMAND "mpirun" ${test_parameters})
  endforeach()
endfunction(add_mpi_test)

# Serial tests
set(numprocs 1)
add_mpi_test(stencil_coefficients "${numprocs}")
add_mpi_test(options_parser "${numprocs}")

# MPI tests
set(numprocs 1 2 4 8)
add_mpi_test(operator_convergence "${numprocs}")
