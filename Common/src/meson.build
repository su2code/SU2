common_include = include_directories('./')
common_src =files(['graph_coloring_structure.cpp',
           'CConfig.cpp',
           'basic_types/ad_structure.cpp',
           'wall_model.cpp',
           '../include/parallelization/mpi_structure.cpp',
           '../include/parallelization/omp_structure.cpp'])

# nvcc cannot parse CoDiPack's tape machinery in the device pass, so the CUDA sources are
# compiled without the CODI_REVERSE_TYPE/CODI_FORWARD_TYPE defines. That gives them a
# different su2double (and hence a different layout for CConfig, CGeometry, ...) than the
# AD libraries, so they are collected separately and only linked into the primal library.
# The AD builds compile the device dispatch out entirely, see SU2_ENABLE_CUDA_KERNELS.
common_cuda_src = []
common_cuda_cpp_args = []
foreach arg : su2_cpp_args
  if arg.startswith('-D') or arg.startswith('-U')
    common_cuda_cpp_args += arg
  endif
endforeach

# Note that the OpenMP flags do reach nvcc (through omp_dep), so the .cu objects agree
# with the rest of the library about HAVE_OMP. They must still not instantiate anything
# containing OpenMP directives, because nvcc's host pass rewrites some of them; the
# device path is synchronized in the headers (see SU2_DEVICE_REGION in CSysVector.hpp),
# which only the .cpp sources instantiate.

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
                          common_src, common_cuda_src,
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
                            cpp_args: [default_warning_flags, su2_cpp_args, codi_rev_args])

  commonAD_dep = declare_dependency(link_with: commonAD,
                                    include_directories : common_include)

endif

if get_option('enable-directdiff')

  commonDD = static_library('SU2CommonDD',
                            common_src,
                            install : false,
                            dependencies : [su2_deps, codi_dep],
                            cpp_args: [default_warning_flags, su2_cpp_args, codi_for_args])

  commonDD_dep = declare_dependency(link_with: commonDD,
                                    include_directories : common_include)

endif
