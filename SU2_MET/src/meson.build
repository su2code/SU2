if get_option('enable-autodiff')

  su2_met_src = ['SU2_MET.cpp']

  su2_cfd_obj_ad = su2_cfd_ad.extract_objects(['error_estimation_structure.cpp',
                                               'solver_structure.cpp',
                                               'solver_adjoint_discrete.cpp',
                                               'solver_direct_mean.cpp',
                            		               'output_structure.cpp',
                                               'output_cgns.cpp',
                                               'output_fieldview.cpp',
                                               'output_inria.cpp',
                                               'output_paraview.cpp',
                                               'output_physics.cpp',
                                               'output_su2.cpp',
                                               'output_tecplot.cpp',
        				                               'variables/CBaselineVariable.cpp',
        				                               'variables/CVariable.cpp'])
  su2_met = executable('SU2_MET',
                       su2_met_src, 
                       install: true,
		     dependencies: [su2_deps, codi_dep, commonAD_dep, inria_dep], 
		     objects : su2_cfd_obj_ad,
		     cpp_args : [default_warning_flags, su2_cpp_args, codi_rev_args, inria_args])

endif
