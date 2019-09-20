# check for swig
swig = find_program('swig')

# add mpi4py include
if mpi_dep[0].found() and mpi_dep[1].found()
  if get_option('python_path') != ''
    python_dir = get_option('python_path')
  else
    python_dir = run_command(python, '-c', 'from distutils.sysconfig import get_python_lib; print(get_python_lib())').stdout().strip()
    message ('Looking for mpi4py in ' + python_dir)
  endif
  mpi4py_include = python_dir + '/mpi4py/include' 
else
  mpi4py_include = ''
endif

swig_gen = generator(
    swig,
    output: ['@BASENAME@.cxx'],
    arguments: su2_cpp_args +
    [ '-c++', '-python', '-I'+mpi4py_include, '-outdir', meson.current_build_dir(), '-o', './@OUTPUT@', '@INPUT@'],
    depfile: '@BASENAME@.d',
)

wrapper_deps = [
    dependency('python3'), 
]

if get_option('enable-normal')
  cpp_source = swig_gen.process('pySU2.i')
  swig_lib = shared_library(
      '_pysu2',
      cpp_source,
      dependencies: [wrapper_deps, common_dep, su2_deps],
      objects: su2_cfd.extract_all_objects(),
      install: true,
      include_directories : mpi4py_include,
      cpp_args : [default_warning_flags,su2_cpp_args],
      name_prefix : '',
      install_dir: 'bin'
  )
meson.add_install_script('install.sh', meson.current_build_dir(), 'pysu2.py')
endif

if get_option('enable-autodiff')
  cpp_source = swig_gen.process('pySU2ad.i')
  swig_lib = shared_library(
      '_pysu2ad',
      cpp_source,
      dependencies: [wrapper_deps, commonAD_dep, su2_deps, codi_dep],
      objects: su2_cfd_ad.extract_all_objects(),
      install: true,
      include_directories : mpi4py_include,
      cpp_args : [default_warning_flags, su2_cpp_args, codi_rev_args],
      name_prefix : '',
      install_dir: 'bin'
  )
meson.add_install_script('install.sh', meson.current_build_dir(), 'pysu2ad.py')
endif
