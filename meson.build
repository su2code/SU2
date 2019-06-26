project('SU2', 'c', 'cpp', 
        version: '6.2.0', 
	license: 'LGPL2',
	default_options: ['buildtype=release', 'warning_level=0'])


pymod = import('python')
python = pymod.find_installation('python3')

su2_cpp_args = []
su2_deps     = []

# external dependencies 
mpi_dep = [dependency('mpi', language:'c', required: false), 
           dependency('mpi', language:'cpp', required: false)]

codi_dep = [declare_dependency(include_directories: 'externals/codi/include')]
codi_rev_args = '-DCODI_REVERSE_TYPE'
codi_for_args = '-DCODI_FORWARD_TYPE'

# Add cgns library
if get_option('enable-cgns')
  subdir('externals/cgns')
  su2_deps     += cgns_dep
  su2_cpp_args += '-DHAVE_CGNS'
endif

# check if MPI dependencies are found and add them
if mpi_dep[0].found() and mpi_dep[1].found() and get_option('enable-mpi')

  # add MPI dependency 
  su2_deps += mpi_dep
  su2_cpp_args += '-DHAVE_MPI'
  
  # compile metis
  subdir('externals/metis')
  
  # add metis dependency
  su2_deps     += metis_dep
  su2_cpp_args += '-DHAVE_METIS' 
  
  # compile parmetis
  subdir('externals/parmetis')
 
  subdir('externals/tecio')

  # add parmetis dependency
  su2_deps     += parmetis_dep
  su2_cpp_args += '-DHAVE_PARMETIS' 
  if get_option('enable-tecio')
    subdir('externals/tecio/teciompisrc')
    su2_deps += teciompi_dep
    su2_cpp_args += ['-DHAVE_TECIO', '-DTECIOMPI']
  endif
  
  # add medi dependency
  if get_option('enable-autodiff')
    codi_dep += declare_dependency(include_directories: ['externals/medi/include', 'externals/medi/src'])
  endif

else
  if get_option('enable-tecio')
    subdir('externals/tecio/teciosrc')
    pthread_dep = dependency('threads')
    su2_deps += [tecio_dep, pthread_dep]
    su2_cpp_args += '-DHAVE_TECIO'
  endif
endif

# compile common library
subdir('Common/src')
# compile SU2_CFD executable
subdir('SU2_CFD/src')
# compile SU2_DOT executable
subdir('SU2_DOT/src')
# compile SU2_DEF executable
subdir('SU2_DEF/src')
# compile SU2_GEO executable
subdir('SU2_GEO/src')
# compile SU2_MSH executable
subdir('SU2_MSH/src')
# compile SU2_SOL executable
subdir('SU2_SOL/src')
# install python scripts
subdir('SU2_PY')

if get_option('enable-pywrapper')
  subdir('SU2_PY/pySU2')
endif

ninja = find_program('ninja', required : false)
if ninja.found() == false
  message('Bootstrapping Ninja ...')
  run_config = run_command('externals/ninja/configure.py','--bootstrap')
endif
