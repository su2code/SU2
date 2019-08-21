
loguru_lib = static_library('loguru', ['loguru.cpp'], dependencies: dependency('threads'), cpp_args: '-DLOGURU_WITH_STREAMS')

loguru_dep = declare_dependency(link_with: loguru_lib, include_directories: '.', link_args: '-ldl')


