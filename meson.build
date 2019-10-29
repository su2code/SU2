project('SU2', 'c', 'cpp', 
        version: '6.2.0', 
        license: 'LGPL2',
        default_options: ['buildtype=release', 
                          'warning_level=0', 
                          'cpp_std=c++11'])


pymod = import('python')
python = pymod.find_installation()

su2_cpp_args = []
su2_deps     = [declare_dependency(include_directories: 'externals/CLI11')]

default_warning_flags = ['-Wno-unused-parameter',
                         '-Wno-empty-body',
                         '-Wno-format-security',
                         '-Wno-deprecated-declarations',
                         '-Wno-non-virtual-dtor']

# external dependencies 
mpi_dep = [dependency('mpi', language:'c', required : get_option('with-mpi')), 
           dependency('mpi', language:'cpp', required : get_option('with-mpi'))]

if get_option('enable-autodiff') or get_option('enable-directdiff')
  codi_dep = [declare_dependency(include_directories: 'externals/codi/include')]
  codi_rev_args = '-DCODI_REVERSE_TYPE'
  codi_for_args = '-DCODI_FORWARD_TYPE'
endif

# add cgns library
if get_option('enable-cgns')
  subdir('externals/cgns')
  su2_deps     += cgns_dep
  su2_cpp_args += '-DHAVE_CGNS'
endif

# Check for non-debug build
if get_option('buildtype')!='debug'
  su2_cpp_args += '-DNDEBUG'
endif

# check if MPI dependencies are found and add them
if mpi_dep[0].found() and mpi_dep[1].found() 

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
 
  # add parmetis dependency
  su2_deps     += parmetis_dep
  su2_cpp_args += '-DHAVE_PARMETIS' 
  
  # add medi dependency
  if get_option('enable-autodiff') or get_option('enable-directdiff')
    codi_dep += declare_dependency(include_directories: ['externals/medi/include', 'externals/medi/src'])
  endif
endif

if get_option('enable-tecio')
  subdir('externals/tecio')
endif

# PaStiX
if get_option('enable-pastix')
  assert(mpi_dep[0].found() and mpi_dep[1].found(),
         'PaStiX support requires MPI')
  assert(get_option('enable-mkl') or get_option('enable-openblas'),
         'PaStiX support requires either mkl or openblas')

  su2_cpp_args += '-DHAVE_PASTIX'

  pastix_root = get_option('pastix_root')+'/install'
  scotch_root = get_option('scotch_root')+'/lib'
  pastix_dep  = declare_dependency(include_directories: pastix_root,
                                   link_args: [ '-L../'+pastix_root, '-L../'+scotch_root,
                                   '-lpastix', '-lscotch', '-lptscotch', '-lptscotcherr',
                                   '-lm', '-lrt', '-lpthread'])
  su2_deps += pastix_dep
endif

# blas-type dependencies
if get_option('enable-mkl')

  su2_cpp_args += '-DHAVE_MKL'
  mkl_root = get_option('mkl_root')
  mkl_dep = declare_dependency(include_directories: mkl_root+'/include', compile_args: '-m64',
    link_args: ['-L'+mkl_root+'/lib/intel64','-lmkl_intel_lp64','-lmkl_sequential','-lmkl_core','-lpthread','-lm','-ldl'])
  su2_deps += mkl_dep

elif get_option('enable-openblas')

  su2_cpp_args += '-DHAVE_BLAS'
  su2_cpp_args += '-DHAVE_LAPACK'
  blas_dep = dependency(get_option('blas-name'))
  su2_deps += blas_dep

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

message('''-------------------------------------------------------------------------
         |    ___ _   _ ___                                                      | 
         |   / __| | | |_  )   Release 6.2.0 'Falcon'                            | 
         |   \__ \ |_| |/ /                                                      | 
         |   |___/\___//___|   Configuration Summary                             | 
         |                                                                       | 
         ------------------------------------------------------------------------- 
         
         Option          Value
         ---------------------
         TecIO:          @2@
         CGNS:           @3@
         AD (reverse):   @4@
         AD (forward):   @5@
         Python Wrapper: @6@
         Intel-MKL:      @7@
         OpenBlas:       @8@
         PaStiX:         @9@

         Please be sure to add the $SU2_HOME and $SU2_RUN environment variables,
         and update your $PATH (and $PYTHONPATH if applicable) with $SU2_RUN

         Based on the input to this configuration, add these lines to your .bashrc file:
        
         export SU2_RUN=@0@
         export SU2_HOME=@1@
         export PATH=$PATH:$SU2_RUN
         export PYTHONPATH=$PYTHONPATH:$SU2_RUN

         Use './ninja -C @10@ install' to compile and install SU2
'''.format(get_option('prefix')+'/bin', meson.source_root(), get_option('enable-tecio'), get_option('enable-cgns'), 
           get_option('enable-autodiff'), get_option('enable-directdiff'), get_option('enable-pywrapper'), get_option('enable-mkl'),
           get_option('enable-openblas'), get_option('enable-pastix'), meson.build_root().split('/')[-1]))
        
