su2_def_include = include_directories('./')
su2_def_src = files([
  'SU2_DEF.cpp',
  'drivers/CDeformationDriver.cpp',
  'drivers/CDiscAdjDeformationDriver.cpp'
])
if get_option('enable-normal')

su2_cfd_obj = su2_cfd_lib.extract_objects(['solvers/CSolver.cpp',
                                      'solvers/CBaselineSolver.cpp',
                                      'solvers/CMeshSolver.cpp',
                                      'solvers/CFEASolver.cpp',
                                      'numerics/CNumerics.cpp',
                                      'numerics/elasticity/CFEAElasticity.cpp',
                                      'numerics/elasticity/CFEALinearElasticity.cpp',
                                      'CMarkerProfileReaderFVM.cpp',
                                      'output/COutput.cpp',
                                      'output/tools/CWindowingTools.cpp',
                                      'output/CMeshOutput.cpp',
                                      'output/output_structure_legacy.cpp',
                                      'output/CBaselineOutput.cpp',
                                      'variables/CBaselineVariable.cpp',
                                      'variables/CVariable.cpp',
                                      'variables/CFEAVariable.cpp',
                                      'variables/CFEABoundVariable.cpp',
                                      'variables/CMeshElement.cpp',
                                      'variables/CMeshVariable.cpp',
                                      'variables/CMeshBoundVariable.cpp',
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
                                      'output/filewriter/CSTLFileWriter.cpp',
                                      'output/filewriter/CSU2FileWriter.cpp',
                                      'output/filewriter/CSU2BinaryFileWriter.cpp',
                                      'output/filewriter/CParaviewXMLFileWriter.cpp',
                                      'output/filewriter/CParaviewVTMFileWriter.cpp',
                                      'output/filewriter/CSU2MeshFileWriter.cpp',
                                      'limiters/CLimiterDetails.cpp'])

  su2_def_lib = static_library('SU2core_DEF',
                               su2_def_src,
                               install : false,
                               dependencies : [su2_deps, common_dep],
                               cpp_args:  [default_warning_flags, su2_cpp_args])

  su2_def_dep = declare_dependency(link_with: su2_def_lib,
                                   include_directories: su2_def_include)

  su2_def = executable('SU2_DEF',
                        su2_def_src,
                        install: true,
                        dependencies: [su2_deps, common_dep],
                        objects : su2_cfd_obj,
                        cpp_args :[default_warning_flags, su2_cpp_args])
endif

if get_option('enable-autodiff')
su2_cfd_obj_ad = su2_cfd_lib_ad.extract_objects(['solvers/CSolver.cpp',
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
                                                 'variables/CVariable.cpp',
                                                 'limiters/CLimiterDetails.cpp'])

  su2_def_lib_ad = static_library('SU2core_DEF_AD',
                                  su2_def_src,
                                  install : false,
                                  dependencies : [su2_deps, codi_dep, commonAD_dep],
                                  cpp_args:  [default_warning_flags, su2_cpp_args, codi_rev_args])

  su2_def_dep_ad = declare_dependency(link_with: su2_def_lib_ad,
                                     include_directories: su2_def_include)

endif
