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
add_mpi_test(inviscid_flux_jacobian "${numprocs}")
add_mpi_test(impenetrable_wall_SAT "${numprocs}")
add_mpi_test(isothermal_wall_SAT "${numprocs}")

find_package(LAPACK)
if (LAPACK_FOUND)
  add_mpi_test(incoming_inviscid_flux_jacobian "${numprocs}")
  target_link_libraries(incoming_inviscid_flux_jacobian ${LAPACK_LIBRARIES})
endif()

# MPI tests
set(numprocs 1 2 4 8)
add_mpi_test(operator_convergence "${numprocs}")
add_mpi_test(adjoint_relation "${numprocs}")
add_mpi_test(dissipation_self_adjoint "${numprocs}")
