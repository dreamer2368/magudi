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
add_serial_test(stencil_coefficients)
add_serial_test(options_parser)
add_serial_test(inviscid_flux_jacobian)
add_serial_test(impenetrable_wall_SAT)
add_serial_test(isothermal_wall_SAT)
add_serial_test(drag_adjoint_forcing)

find_package(LAPACK)
if (LAPACK_FOUND)
  add_mpi_test(incoming_inviscid_flux_jacobian "${numProcs}")
  target_link_libraries(incoming_inviscid_flux_jacobian ${LAPACK_LIBRARIES})
endif()

# MPI tests
set(numProcs 1 2 4 8)
add_mpi_test(operator_convergence "${numProcs}")
##add_mpi_test(adjoint_relation "${numProcs}")
add_mpi_test(composite_dissipation_sanity "${numProcs}")
add_mpi_test(dissipation_self_adjoint "${numProcs}")
add_mpi_test(SBP_property "${numProcs}")
add_mpi_test(patch_collectives "${numProcs}")
add_mpi_test(boundary_operator "${numProcs}")
add_mpi_test(testZAXPY "${numProcs}")
add_mpi_test(block_interface_sign "${numProcs}")

# Subdirectories
add_subdirectory(adjoint_relation)
add_subdirectory(linearized_relation)
