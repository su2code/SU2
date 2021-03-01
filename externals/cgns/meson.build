if build_machine.system() == 'windows' or meson.get_compiler('cpp').get_id() == 'intel'
  cgns_default_warnings = []
else
  cgns_default_warnings = ['-Wno-unused-result',
                           '-Wno-unused-parameter',
                           '-Wno-unused-variable',
                           '-Wno-unused-but-set-variable',
                           '-Wno-sign-compare',
                           '-Wno-pedantic']
endif

cgns_include = include_directories('adf', './')

lib_cgns = static_library('cgns',
                          'cgns_internals.c',
                          'cgns_error.c',
                          'cgnslib.c', 
                          'cgns_io.c',
                          'adf/ADF_interface.c',
                          'adf/ADF_internals.c',
                          install : false,
                          c_args: cgns_default_warnings)

cgns_dep = declare_dependency(link_with: lib_cgns, 
                              include_directories: cgns_include)

