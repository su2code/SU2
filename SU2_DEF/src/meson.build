su2_def_src = ['SU2_DEF.cpp']
su2_cfd_obj = su2_cfd.extract_objects(['solver_structure.cpp',
                                       'output/COutput.cpp',
                                       'output/CMeshOutput.cpp',
                                       'output/output_structure_legacy.cpp',
                                       'output/output_csv.cpp',
				       'output/output_cgns.cpp',
				       'output/output_tecplot.cpp',
				       'output/output_fieldview.cpp',
				       'output/output_su2.cpp',
				       'output/output_paraview.cpp',
				       'variable_structure.cpp'])


su2_def = executable('SU2_DEF',
                     su2_def_src, 
                     install: true,
		     dependencies: [su2_deps, common_dep], 
		     objects : su2_cfd_obj,
		     cpp_args : su2_cpp_args)

