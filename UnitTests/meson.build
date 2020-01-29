# Add any new test files here
su2_cfd_tests = files(['Common/geometry/primal_grid/CPrimalGrid_tests.cpp',
                       'Common/geometry/dual_grid/CDualGrid_tests.cpp',
                       'SU2_CFD/numerics/CNumerics_tests.cpp'])

# Add the test driver (the main() definition)
unit_test_files = su2_cfd_tests + files(['test_driver.cpp'])

# Ideally, I would like to add `build_by_default : false` to this, and
# only compile the tests when `ninja test` is run.  But the minute you
# add an executable to a test, it overrides `build_by_default` and will
# always compile the tests. Changing the behavior is an an open issue /
# pull request on the meson Github.
test_driver = executable(
    'test_driver',
    unit_test_files,
    install : false,
    dependencies : [su2_cfd_dep, common_dep, su2_deps, catch2_dep],
    cpp_args: ['-fPIC', default_warning_flags, su2_cpp_args]

)

test('Catch2 test driver', test_driver)
