su2_def_src = ['SU2_DEF.cpp']
if get_option('enable-normal')

  su2_cfd_obj = su2_cfd.extract_objects(['solver_structure.cpp',
                                        'CMarkerProfileReaderFVM.cpp',
                                        'output/COutput.cpp',
					'output/tools/CWindowingTools.cpp',
                                        'output/CMeshOutput.cpp',
                                        'output/output_structure_legacy.cpp',
                                        'variables/CBaselineVariable.cpp',
                                        'variables/CVariable.cpp',
                                        'output/filewriter/CParallelDataSorter.cpp',
                                        'output/filewriter/CFVMDataSorter.cpp',
                                        'output/filewriter/CFEMDataSorter.cpp',
                                        'output/filewriter/CSurfaceFEMDataSorter.cpp',
                                        'output/filewriter/CSurfaceFVMDataSorter.cpp',
                                        'output/filewriter/CParallelFileWriter.cpp',
                                        'output/filewriter/CParaviewFileWriter.cpp',
                                        'output/filewriter/CParaviewBinaryFileWriter.cpp',
                                        'output/filewriter/CTecplotFileWriter.cpp',
                                        'output/filewriter/CTecplotBinaryFileWriter.cpp',
                                        'output/filewriter/CCSVFileWriter.cpp',
                                        'output/filewriter/CSU2FileWriter.cpp',
                                        'output/filewriter/CSU2BinaryFileWriter.cpp',
                                        'output/filewriter/CSU2MeshFileWriter.cpp'])

  su2_def = executable('SU2_DEF',
                        su2_def_src, 
                        install: true,
                        dependencies: [su2_deps, common_dep], 
                        objects : su2_cfd_obj,
                        cpp_args :[default_warning_flags, su2_cpp_args])
endif
