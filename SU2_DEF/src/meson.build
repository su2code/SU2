su2_def_src = ['SU2_DEF.cpp']
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
  su2_def = executable('SU2_DEF',
                       su2_def_src, 
                       install: true,
		       dependencies: [su2_deps, common_dep], 
		       objects : su2_cfd_obj,
  		       cpp_args :[default_warning_flags, su2_cpp_args])
endif
