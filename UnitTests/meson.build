# Add any new test files here
# -------------------------------------------------------------------------
# Begin unit test listings
# -------------------------------------------------------------------------

# Direct-mode tests:
su2_cfd_tests = files(['Common/geometry/primal_grid/CPrimalGrid_tests.cpp',
                       'Common/geometry/dual_grid/CDualGrid_tests.cpp',
                       'Common/geometry/CGeometry_test.cpp',
                       'Common/toolboxes/CQuasiNewtonInvLeastSquares_tests.cpp',
                       'Common/toolboxes/C1DInterpolation_tests.cpp',
                       'Common/vectorization.cpp',
                       'Common/toolboxes/ndflattener_tests.cpp',
                       'Common/containers/CLookupTable_tests.cpp',
                       'Common/toolboxes/multilayer_perceptron/CLookUp_ANN_tests.cpp',
                       'SU2_CFD/numerics/CNumerics_tests.cpp',
                       'SU2_CFD/gradients.cpp',
                       'SU2_CFD/windowing.cpp'])

# Reverse-mode (algorithmic differentiation) tests:
su2_cfd_tests_ad = files(['Common/simple_ad_test.cpp'])

# Forward-mode (direct differentiation) tests:
su2_cfd_tests_dd = files(['Common/simple_directdiff_test.cpp'])

# -------------------------------------------------------------------------
# End of unit test listings
# -------------------------------------------------------------------------

if get_option('enable-tests')
  if get_option('enable-normal')
    unit_test_files = su2_cfd_tests + files(['test_driver.cpp'])
    test_driver = executable(
        'test_driver',
        unit_test_files,
        install : true,
        dependencies : [su2_cfd_dep, common_dep, su2_deps, catch2_dep],
        cpp_args: ['-fPIC', default_warning_flags, su2_cpp_args]
    )
    test('Catch2 test driver', test_driver)
  endif

  if get_option('enable-autodiff')
    unit_test_files_ad = su2_cfd_tests_ad + files(['test_driver.cpp'])
    test_driver_AD = executable(
        'test_driver_AD',
        unit_test_files_ad,
        install : true,
        dependencies : [su2_cfd_dep_ad, commonAD_dep, su2_deps, codi_dep, catch2_dep],
        cpp_args: ['-fPIC', default_warning_flags, su2_cpp_args, codi_rev_args]
    )
    test('Catch2 test driver (AD)', test_driver_AD)
  endif

  if get_option('enable-directdiff')
    unit_test_files_dd = su2_cfd_tests_dd + files(['test_driver.cpp'])
    test_driver_DD = executable(
        'test_driver_DD',
        unit_test_files_dd,
        install : true,
        dependencies : [su2_cfd_dep_dd, commonDD_dep, su2_deps, codi_dep, catch2_dep],
        cpp_args: ['-fPIC', default_warning_flags, su2_cpp_args, codi_for_args]
    )
    test('Catch2 test driver (DD)', test_driver_DD)
  endif
endif
