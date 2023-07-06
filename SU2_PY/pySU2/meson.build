# check for swig
swig = find_program('swig')
# check version, manually because meson tries --version but swig uses -version.
swig_version = run_command('swig', '-version').stdout().split('\n')
foreach line : swig_version
    if line.startswith('SWIG Version ')
        if line.replace('SWIG Version ','').version_compare('<3.0.11')
            error('The SU2 Python wrapper requires SWIG version >= 3.0.11')
        endif
    endif
endforeach

# add mpi4py include
if mpi
    mpi4py_include = run_command(python, '-c', 'import mpi4py; print(mpi4py.get_include())').stdout().strip()
    assert(not mpi4py_include.contains('Traceback'), 'python does not have mpi4py module')
    message('Using mpi4py from ' + mpi4py_include)
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
    python.dependency(embed: true),
]

if get_option('enable-normal')
  cpp_source = swig_gen.process('pySU2.i')
  swig_lib = shared_library(
      '_pysu2',
      cpp_source,
      dependencies: [wrapper_deps, common_dep, su2_deps],
      objects: [
        su2_cfd_lib.extract_all_objects(),
        su2_def_lib.extract_objects('drivers/CDeformationDriver.cpp')
      ],
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
      objects: [
        su2_cfd_lib_ad.extract_all_objects(),
        su2_def_lib_ad.extract_objects('drivers/CDiscAdjDeformationDriver.cpp')
      ],
      install: true,
      include_directories : mpi4py_include,
      cpp_args : [default_warning_flags, su2_cpp_args, codi_rev_args],
      name_prefix : '',
      install_dir: 'bin'
  )
meson.add_install_script('install.sh', meson.current_build_dir(), 'pysu2ad.py')
endif
