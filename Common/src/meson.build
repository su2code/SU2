common_include = include_directories('./')
common_src =files(['graph_coloring_structure.cpp',
           'CConfig.cpp',
           'basic_types/ad_structure.cpp',
           'wall_model.cpp',
           '../include/parallelization/mpi_structure.cpp',
           '../include/parallelization/omp_structure.cpp'])

common_cuda_cpp_args = []
foreach arg : su2_cpp_args
  if arg.startswith('-D') or arg.startswith('-U')
    common_cuda_cpp_args += arg
  endif
endforeach

common_cuda_rev_args = common_cuda_cpp_args
foreach arg : codi_rev_args
  if arg.startswith('-D') or arg.startswith('-U')
    common_cuda_rev_args += arg
  endif
endforeach

common_cuda_for_args = common_cuda_cpp_args
foreach arg : codi_for_args
  if arg.startswith('-D') or arg.startswith('-U')
    common_cuda_for_args += arg
  endif
endforeach

subdir('linear_algebra')
subdir('toolboxes')
subdir('geometry')
subdir('geometry/elements')
subdir('geometry/dual_grid')
subdir('geometry/primal_grid')
subdir('geometry/meshreader')
subdir('containers')
subdir('interface_interpolation')
subdir('fem')
subdir('grid_movement')
subdir('adt')

if get_option('enable-normal')
  common = static_library('SU2Common',
                          common_src,
                          install : false,
                          dependencies : su2_deps,
                          cpp_args: [default_warning_flags, su2_cpp_args],
                          cuda_args: common_cuda_cpp_args)

  common_dep = declare_dependency(link_with: common,
                                  include_directories : common_include)
endif

if get_option('enable-autodiff')

  commonAD = static_library('SU2CommonAD',
                            common_src,
                            install : false,
                            dependencies : [su2_deps, codi_dep],
                            cpp_args: [default_warning_flags, su2_cpp_args, codi_rev_args],
                            cuda_args: common_cuda_rev_args)

  commonAD_dep = declare_dependency(link_with: commonAD,
                                    include_directories : common_include)

endif

if get_option('enable-directdiff')

  commonDD = static_library('SU2CommonDD',
                            common_src,
                            install : false,
                            dependencies : [su2_deps, codi_dep],
                            cpp_args: [default_warning_flags, su2_cpp_args, codi_for_args],
                            cuda_args: common_cuda_for_args)

  commonDD_dep = declare_dependency(link_with: commonDD,
                                    include_directories : common_include)

endif
