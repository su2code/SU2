su2_msh_src = ['SU2_MSH.cpp']

if get_option('enable-normal')
  su2_msh = executable('SU2_MSH',
                       su2_msh_src, 
                       install: true,
		       dependencies: [su2_deps, common_dep], 
  		       cpp_args :[default_warning_flags, su2_cpp_args])
endif
