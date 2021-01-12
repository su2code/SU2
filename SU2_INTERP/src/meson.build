su2_interp_src = ['SU2_INTERP.cpp']
if get_option('enable-normal')
  su2_interp = executable('SU2_INTERP',
                      su2_interp_src, 
                      install: true,
                      dependencies: [su2_deps, common_dep, su2_cfd_dep], 
                      cpp_args :[default_warning_flags, su2_cpp_args])


endif