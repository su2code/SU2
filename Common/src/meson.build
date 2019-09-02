common_include = include_directories('./')
common_src =files(['geometry_structure_fem_part.cpp',
           'primal_grid_structure.cpp',
           'graph_coloring_structure.cpp',
           'fem_geometry_structure.cpp',
           'fem_integration_rules.cpp',
           'config_structure.cpp',
           'dual_grid_structure.cpp',
           'element_structure.cpp',
           'geometry_structure.cpp',
           'blas_structure.cpp',
           'ad_structure.cpp',
           'grid_movement_structure.cpp',
           'fem_work_estimate_metis.cpp',
           'fem_standard_element.cpp',
           'fem_wall_distance.cpp',
           'grid_adaptation_structure.cpp',
           'element_linear.cpp',
           'fem_gauss_jacobi_quadrature.cpp',
           'gauss_structure.cpp',
           'wall_model.cpp',
           'adt_structure.cpp',
           'interpolation_structure.cpp',
           'mpi_structure.cpp',
           'fem_cgns_elements.cpp',
           'CMultiGridQueue.cpp',
           'CMeshReaderFVM.cpp',
           'CSU2ASCIIMeshReaderFVM.cpp',
           'CBoxMeshReaderFVM.cpp',
           'CCGNSMeshReaderFVM.cpp',
           'CRectangularMeshReaderFVM.cpp'])

subdir('linear_algebra')
subdir('toolboxes')
if get_option('enable-normal')
  common = static_library('SU2Common', 
                          common_src,
                          install : false, 
	                  dependencies : su2_deps, 
	 		  cpp_args: [default_warning_flags, su2_cpp_args])
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
