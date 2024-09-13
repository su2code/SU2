su2_geo_src = ['SU2_GEO.cpp']

if get_option('enable-normal')
  su2_geo = executable('SU2_GEO',
                       su2_geo_src,
                       install: true,
		       dependencies: [su2_deps, common_dep],
		       cpp_args : [default_warning_flags, su2_cpp_args])
endif
