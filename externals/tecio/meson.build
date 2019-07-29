check_dir = run_command('test', '-d', 'boost')
if check_dir.returncode() != 0
  message('Extracting boost ...')
  extract_boost = run_command('tar', '-xzf', 'boost.tar.gz', '-C', './')
else
  message('Boost sources found.')
endif

if get_option('enable-mpi')
  subdir('teciompisrc')
  su2_deps += teciompi_dep
  su2_cpp_args += ['-DHAVE_TECIO', '-DTECIOMPI']
else
  subdir('teciosrc')
  pthread_dep = dependency('threads')
  su2_deps += [tecio_dep, pthread_dep]
  su2_cpp_args += '-DHAVE_TECIO'
endif
