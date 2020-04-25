check_dir = run_command(python, 
                        script_path / 'check_dir.py',
                        'boost')
if check_dir.returncode() != 0
  message('Extracting boost ...')
  extract_boost = run_command(python, 
                              script_path / 'extract_file.py',
                              'boost.tar.gz', 
                              meson.current_source_dir(), check: true)
else
  message('Boost sources found.')
endif

if mpi
  subdir('teciompisrc')
  su2_deps += teciompi_dep
  su2_cpp_args += ['-DHAVE_TECIO', '-DTECIOMPI']
else
  subdir('teciosrc')
  pthread_dep = dependency('threads')
  su2_deps += [tecio_dep, pthread_dep]
  su2_cpp_args += '-DHAVE_TECIO'
endif
