su2_sol_src = ['SU2_SOL.cpp']
if get_option('enable-normal')

  su2_sol = executable('SU2_SOL',
                      su2_sol_src,
                      install: true,
                      dependencies: [su2_deps, common_dep, su2_cfd_dep],
                      cpp_args :[default_warning_flags, su2_cpp_args])

endif
