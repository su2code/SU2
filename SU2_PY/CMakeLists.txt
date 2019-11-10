project(SU2_PY)

set(SU2_PY_SCRIPTS
    continuous_adjoint.py
    compute_uncertainty.py
    finite_differences.py
    mesh_deformation.py
    parallel_computation.py
    parallel_computation_fsi.py
    package_tests.py
    shape_optimization.py
    merge_solution.py
    set_ffd_design_var.py
    compute_polar.py
    compute_multipoint.py
    discrete_adjoint.py
    direct_differentiation.py
    fsi_computation.py
    SU2_CFD.py)

add_library(SU2_PY INTERFACE)
target_include_directories(SU2_PY INTERFACE ${CMAKE_CURRENT_LIST_DIR})
install(FILES ${SU2_PY_SCRIPTS} DESTINATION bin)
install(DIRECTORY SU2 FSI DESTINATION bin)
