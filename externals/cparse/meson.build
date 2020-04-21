cparse_include = include_directories('./')

lib_cparse_src = files(['shunting-yard.cpp',
                        'packToken.cpp',
			'functions.cpp',
                        'containers.cpp'])

if get_option('enable-normal')

  lib_cparse = static_library('cparse',
                              lib_cparse_src,
                              install : false)

  cparse_dep = declare_dependency(link_with: lib_cparse,
                                  include_directories: cparse_include)

endif

if get_option('enable-autodiff')

  lib_cparseAD = static_library('cparseAD',
                                lib_cparse_src,
                                install : false,
                                dependencies : codi_dep,
                                cpp_args: codi_rev_args)

  cparse_depAD = declare_dependency(link_with: lib_cparseAD,
                                    include_directories: cparse_include)

endif

if get_option('enable-directdiff')

  lib_cparseDD = static_library('cparseDD',
                                lib_cparse_src,
                                install : false,
                                dependencies : codi_dep,
                                cpp_args: codi_for_args)

  cparse_depDD = declare_dependency(link_with: lib_cparseDD,
                                    include_directories: cparse_include)

endif
