su2_dot_src = ['SU2_DOT.cpp']
if get_option('enable-normal')
  su2_cfd_obj = su2_cfd.extract_objects(['solver_structure.cpp',
                    		         'output_structure.cpp',
                                         'output_cgns.cpp',
                                         'output_fieldview.cpp',
                                         'output_paraview.cpp',
                                         'output_physics.cpp',
                                         'output_su2.cpp',
                                         'output_tecplot.cpp',
				         'variables/CBaselineVariable.cpp',
				         'variables/CVariable.cpp'])
  su2_dot = executable('SU2_DOT',
                       su2_dot_src, 
                       install: true,
  		       dependencies: [su2_deps, common_dep], 
		       objects : su2_cfd_obj,
                       cpp_args :[default_warning_flags, su2_cpp_args])
endif

if get_option('enable-autodiff')
  su2_dot_src_ad = ['SU2_DOT.cpp']
  su2_cfd_obj_ad = su2_cfd_ad.extract_objects(['solver_structure.cpp',
                    		       'output_structure.cpp',
                                       'output_cgns.cpp',
                                       'output_fieldview.cpp',
                                       'output_paraview.cpp',
                                       'output_physics.cpp',
                                       'output_su2.cpp',
                                       'output_tecplot.cpp',
				       'variables/CBaselineVariable.cpp',
				       'variables/CVariable.cpp'])
  su2_dot_ad = executable('SU2_DOT_AD',
                     su2_dot_src_ad, 
                     install: true,
		     dependencies: [su2_deps, codi_dep, commonAD_dep], 
		     objects : su2_cfd_obj_ad,
		     cpp_args : [default_warning_flags, su2_cpp_args, codi_rev_args])

endif
