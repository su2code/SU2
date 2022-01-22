if build_machine.system() == 'windows' or meson.get_compiler('cpp').get_id() == 'intel'
  cgns_default_warnings = []
else
  cgns_default_warnings = ['-Wno-unused-result',
                           '-Wno-unused-parameter',
                           '-Wno-unused-variable',
                           '-Wno-unused-but-set-variable',
                           '-Wno-sign-compare',
                           '-Wno-error=unused-function',
                           '-Wno-pedantic']
endif

subdir('hdf5')


cgns_conf_data = configuration_data()
cgns_conf_data.set('SIZEOF_LONG', cc.sizeof('long'))
cgns_conf_data.set('SIZEOF_VOID_P', cc.sizeof('void *'))


cgns_sources = ['cgns_internals.c',
'cgns_error.c',
'cgnslib.c',
'cgns_io.c',
'cg_hashmap.c',
'adf/ADF_interface.c',
'adf/ADF_internals.c',
'adfh/ADFH.c']

if mpi
  cgns_conf_data.set('BUILDPARALLEL', 1)
  cgns_conf_data.set('HDF5_HAVE_COLL_METADATA', 1)
  cgns_sources += ['pcgnslib.c']
else
  cgns_conf_data.set('BUILDPARALLEL', 0)
  cgns_conf_data.set('HDF5_HAVE_COLL_METADATA', 0)
endif


configure_file(input: 'cgnstypes.h.in',
  output: 'cgnstypes.h',
  configuration: cgns_conf_data)

configure_file(input: 'cg_hash_types.h.in',
  output: 'cg_hash_types.h',
  configuration: cgns_conf_data)

configure_file(input: 'cgnsconfig.h.in',
  output: 'cgnsconfig.h',
  configuration: cgns_conf_data)


cgns_include = include_directories('adf', 'adfh', './')

lib_cgns = static_library('cgns',
                          cgns_sources,
                          install : false,
                          dependencies : hdf5_dep,
                          c_args: cgns_default_warnings)

cgns_dep = declare_dependency(link_with: lib_cgns,
                              dependencies: hdf5_dep,
                              include_directories: cgns_include)

