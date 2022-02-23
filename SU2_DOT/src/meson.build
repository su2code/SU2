su2_dot_src = 'SU2_DOT.cpp'

if get_option('enable-normal')

  su2_dot = executable('SU2_DOT',
                      su2_dot_src,
                      install: true,
                      dependencies: [su2_deps, common_dep, su2_cfd_dep, su2_def_dep],
                      cpp_args :[default_warning_flags, su2_cpp_args])

endif

if get_option('enable-autodiff')

  su2_dot_ad = executable('SU2_DOT_AD',
                          su2_dot_src,
                          install: true,
                          dependencies: [su2_deps, codi_dep, commonAD_dep, su2_cfd_dep_ad, su2_def_dep_ad],
		                      cpp_args : [default_warning_flags, su2_cpp_args, codi_rev_args])

endif
