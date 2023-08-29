cgns_default_warnings = []

if build_machine.system() != 'windows' and meson.get_compiler('c').get_id() != 'intel'
  desired_warnings = ['-Wno-unused-result',
                      '-Wno-unused-parameter',
                      '-Wno-unused-variable',
                      '-Wno-unused-but-set-variable',
                      '-Wno-sign-compare',
                      '-Wno-error=unused-function',
                      '-Wno-pedantic',
                      '-Wno-error=stringop-truncation',
                      '-Wno-stringop-truncation',
                      '-Wno-error=absolute-value',
                      '-Wno-error=class-memaccess']
  foreach flag : desired_warnings
    if meson.get_compiler('c').has_argument(flag)
      cgns_default_warnings += flag
    endif
  endforeach
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

