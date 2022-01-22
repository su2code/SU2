su2_def_src = ['SU2_DEF.cpp']
if get_option('enable-normal')

  su2_def = executable('SU2_DEF',
                        su2_def_src,
                        install: true,
                        dependencies: [su2_deps, common_dep, su2_cfd_dep],
                        cpp_args :[default_warning_flags, su2_cpp_args])
endif
