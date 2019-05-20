subdir('libmetis')
subdir('GKlib')

metis_include = include_directories('include')
metis_lib     = static_library('metis', 
                               link_with: [libmetis, gklib])

metis_dep = declare_dependency(link_with : metis_lib, 
                               include_directories : metis_include)

