su2_dot_src = ['SU2_DOT.cpp']
su2_cfd_obj = su2_cfd.extract_objects(['solver_structure.cpp',
                                       'output/COutput.cpp',
				       'output/CBaselineOutput.cpp',
				       'output/filewriter/CCSVFileWriter.cpp',
				       'output/output_structure_legacy.cpp',
				       'output/filewriter/CTecplotFileWriter.cpp',
				       'output/filewriter/CTecplotBinaryFileWriter.cpp',
				       'output/filewriter/CSU2FileWriter.cpp',
				       'output/filewriter/CSU2BinaryFileWriter.cpp',
				       'variable_structure.cpp'])


su2_dot = executable('SU2_DOT',
                     su2_dot_src, 
                     install: true,
		     dependencies: [su2_deps, common_dep], 
		     objects : su2_cfd_obj,
		     cpp_args : su2_cpp_args)

