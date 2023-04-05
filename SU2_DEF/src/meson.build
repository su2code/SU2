su2_def_include = include_directories('./')
su2_def_src = files([
  'SU2_DEF.cpp',
  'drivers/CDeformationDriver.cpp',
  'drivers/CDiscAdjDeformationDriver.cpp'
])
if get_option('enable-normal')

  su2_def = executable('SU2_DEF',
                        su2_def_src,
                        install: true,
                        dependencies: [su2_deps, common_dep, su2_cfd_dep],
                        cpp_args :[default_warning_flags, su2_cpp_args])

  su2_def_lib = static_library('SU2core_DEF',
                               su2_def_src,
                               install : false,
                               dependencies : [su2_deps, common_dep, su2_cfd_dep],
                               cpp_args:  [default_warning_flags, su2_cpp_args])

  su2_def_dep = declare_dependency(link_with: su2_def_lib,
                                   include_directories: su2_def_include)

endif

if get_option('enable-autodiff')

  su2_def_lib_ad = static_library('SU2core_DEF_AD',
                                  su2_def_src,
                                  install : false,
                                  dependencies : [su2_deps, codi_dep, commonAD_dep, su2_cfd_dep_ad],
                                  cpp_args:  [default_warning_flags, su2_cpp_args, codi_rev_args])

  su2_def_dep_ad = declare_dependency(link_with: su2_def_lib_ad,
                                      include_directories: su2_def_include)

endif
