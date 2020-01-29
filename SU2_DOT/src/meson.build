su2_dot_src = ['SU2_DOT.cpp']
if get_option('enable-normal')
  su2_cfd_obj = su2_cfd.extract_objects(['solvers/CSolver.cpp',
                                         'solvers/CBaselineSolver.cpp',
                                         'CMarkerProfileReaderFVM.cpp',
                                         'output/COutput.cpp',
                                         'output/tools/CWindowingTools.cpp',
                                         'output/output_structure_legacy.cpp',
                                         'output/CBaselineOutput.cpp',
                                         'output/filewriter/CParallelDataSorter.cpp',
                                         'output/filewriter/CParallelFileWriter.cpp',
                                         'output/filewriter/CFEMDataSorter.cpp',
                                         'output/filewriter/CSurfaceFEMDataSorter.cpp',
                                         'output/filewriter/CFVMDataSorter.cpp',
                                         'output/filewriter/CSurfaceFVMDataSorter.cpp',
                                         'output/filewriter/CCSVFileWriter.cpp',
                                         'output/filewriter/CSTLFileWriter.cpp',
                                         'output/filewriter/CTecplotFileWriter.cpp',
                                         'output/filewriter/CTecplotBinaryFileWriter.cpp',
                                         'output/filewriter/CParaviewFileWriter.cpp',
                                         'output/filewriter/CParaviewBinaryFileWriter.cpp',
                                         'output/filewriter/CSU2FileWriter.cpp',
                                         'output/filewriter/CSU2BinaryFileWriter.cpp',
                                         'output/filewriter/CSU2MeshFileWriter.cpp',
                                         'output/filewriter/CParaviewXMLFileWriter.cpp',
                                         'output/filewriter/CParaviewVTMFileWriter.cpp',
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
  su2_cfd_obj_ad = su2_cfd_ad.extract_objects(['solvers/CSolver.cpp',
                                               'solvers/CBaselineSolver.cpp',
                                               'CMarkerProfileReaderFVM.cpp',
                                               'output/COutput.cpp',
                                               'output/tools/CWindowingTools.cpp',
                                               'output/output_structure_legacy.cpp',
                                               'output/CBaselineOutput.cpp',
                                               'output/filewriter/CParallelDataSorter.cpp',
                                               'output/filewriter/CParallelFileWriter.cpp',
                                               'output/filewriter/CFEMDataSorter.cpp',
                                               'output/filewriter/CSurfaceFEMDataSorter.cpp',
                                               'output/filewriter/CFVMDataSorter.cpp',
                                               'output/filewriter/CSurfaceFVMDataSorter.cpp',
                                               'output/filewriter/CCSVFileWriter.cpp',
                                               'output/filewriter/CSTLFileWriter.cpp',
                                               'output/filewriter/CTecplotFileWriter.cpp',
                                               'output/filewriter/CTecplotBinaryFileWriter.cpp',
                                               'output/filewriter/CParaviewFileWriter.cpp',
                                               'output/filewriter/CParaviewBinaryFileWriter.cpp',
                                               'output/filewriter/CSU2FileWriter.cpp',
                                               'output/filewriter/CSU2BinaryFileWriter.cpp',
                                               'output/filewriter/CSU2MeshFileWriter.cpp',
                                               'output/filewriter/CParaviewXMLFileWriter.cpp',
                                               'output/filewriter/CParaviewVTMFileWriter.cpp',
                                               'variables/CBaselineVariable.cpp',
                                               'variables/CVariable.cpp'])

  su2_dot_ad = executable('SU2_DOT_AD',
                          su2_dot_src_ad, 
                          install: true,
                          dependencies: [su2_deps, codi_dep, commonAD_dep], 
		                      objects : su2_cfd_obj_ad,
		                      cpp_args : [default_warning_flags, su2_cpp_args, codi_rev_args])

endif
