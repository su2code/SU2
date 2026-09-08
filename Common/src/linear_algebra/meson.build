common_src += files(['CSysSolve_b.cpp',
                     'CSysSolve.cpp',
                     'CSysVector.cpp',
                     'CSysMatrix.cpp',
                     'CPastixWrapper.cpp',
                     'blas_structure.cpp'])

if get_option('enable-cuda')
  # Kept apart from common_src: these are compiled without the CoDiPack defines and so
  # must only go into the primal library, see common_cuda_src in Common/src/meson.build.
  common_cuda_src += files(['CSysMatrixGPU.cu', 'CSysVectorGPU.cu'])
endif
