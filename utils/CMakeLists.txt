add_executable(vorticity_dilatation $<TARGET_OBJECTS:magudiObj> vorticity_dilatation.f90)
add_executable(jet_initial_condition $<TARGET_OBJECTS:magudiObj> jet_initial_condition.f90)

add_subdirectory(python)