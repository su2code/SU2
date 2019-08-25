
cgns_include = include_directories('adf', './')

lib_cgns = static_library('cgns',
                          'cgns_internals.c',
			  'cgns_error.c',
			  'cgnslib.c', 
			  'cgns_io.c',
			  'adf/ADF_interface.c',
			  'adf/ADF_internals.c',
                          install : false)

cgns_dep = declare_dependency(link_with: lib_cgns, 
                              include_directories: cgns_include)

