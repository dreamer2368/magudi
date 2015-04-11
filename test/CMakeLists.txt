function(add_mpi_test name numProcs)
  add_executable(${name} $<TARGET_OBJECTS:magudiObj> ${name}.f90)
  foreach(np ${numProcs})
    set(test_parameters -np ${np} "./${name}")
    add_test(NAME ${name}_on_${np}_tasks COMMAND "mpirun" ${test_parameters})
  endforeach()
endfunction(add_mpi_test)

# Serial tests
set(numProcs 1)
add_mpi_test(stencil_coefficients "${numProcs}")
add_mpi_test(options_parser "${numProcs}")
add_mpi_test(inviscid_flux_jacobian "${numProcs}")
add_mpi_test(impenetrable_wall_SAT "${numProcs}")
add_mpi_test(isothermal_wall_SAT "${numProcs}")
add_mpi_test(drag_adjoint_forcing "${numProcs}")

find_package(LAPACK)
if (LAPACK_FOUND)
  add_mpi_test(incoming_inviscid_flux_jacobian "${numProcs}")
  target_link_libraries(incoming_inviscid_flux_jacobian ${LAPACK_LIBRARIES})
endif()

# MPI tests
set(numProcs 1 2 4 8)
add_mpi_test(operator_convergence "${numProcs}")
add_mpi_test(adjoint_relation "${numProcs}")
add_mpi_test(dissipation_self_adjoint "${numProcs}")
add_mpi_test(SBP_property "${numProcs}")
