if build_machine.system() == 'windows' or meson.get_compiler('cpp').get_id() == 'intel'
  rapidjson_default_warnings = []
else
  rapidjson_default_warnings = ['-Wno-unused-result',
                           '-Wno-unused-parameter',
                           '-Wno-unused-variable',
                           '-Wno-unused-but-set-variable',
                           '-Wno-sign-compare',
                           '-Wno-pedantic']
endif

rapidjson_include = include_directories('error', 'internal', './')

rapidjson_dep = declare_dependency(include_directories: rapidjson_include)

