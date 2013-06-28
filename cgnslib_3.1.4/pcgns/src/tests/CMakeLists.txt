#########
# Tests #
#########

if (NOT MPIEXEC)
  message(SEND_ERROR "need to define MPIEXEC to run tests")
endif (NOT MPIEXEC)

# Link all the executables to cgns and hdf5
if (CGNS_BUILD_SHARED AND CGNS_USE_SHARED)
  link_libraries(pcgns_shared)
else (CGNS_BUILD_SHARED AND CGNS_USE_SHARED)
  link_libraries(pcgns_static)
endif (CGNS_BUILD_SHARED AND CGNS_USE_SHARED)
link_libraries(${HDF5_LIBRARY})
link_libraries(${MPI_LIBRARY})
if (MPI_EXTRA_LIBRARY)
  link_libraries(${MPI_EXTRA_LIBRARY})
endif (MPI_EXTRA_LIBRARY)
if (HDF5_NEED_ZLIB)
  link_libraries(${ZLIB_LIBRARY}) 
endif (HDF5_NEED_ZLIB)
if (HDF5_NEED_SZIP)
  link_libraries(${SZIP_LIBRARY}) 
endif (HDF5_NEED_SZIP)
if (NOT WIN32)
  link_libraries(m)
endif (NOT WIN32)

set(benchmark_FILES
	benchmark.c)
add_executable(benchmark ${benchmark_FILES})
add_test(benchmark ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} "2"
         ${MPIEXEC_PREFLAGS} benchmark ${MPIEXEC_POSTFLAGS})

set(open_close_FILES
	open_close.c)
add_executable(open_close ${open_close_FILES})
add_test(open_close ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} "2"
         ${MPIEXEC_PREFLAGS} open_close ${MPIEXEC_POSTFLAGS})

set(test_base_FILES
	test_base.c)
add_executable(test_base ${test_base_FILES})
add_test(test_base ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} "2"
         ${MPIEXEC_PREFLAGS} test_base ${MPIEXEC_POSTFLAGS})

set(test_zone_FILES
	test_zone.c)
add_executable(test_zone ${test_zone_FILES})
add_test(test_zone ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} "2"
         ${MPIEXEC_PREFLAGS} test_zone ${MPIEXEC_POSTFLAGS})

set(test_unstructured_FILES
	test_unstructured.c)
add_executable(test_unstructured ${test_unstructured_FILES})
add_test(test_unstructured ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} "2"
         ${MPIEXEC_PREFLAGS} test_unstructured ${MPIEXEC_POSTFLAGS})

set(test_queue_FILES
	test_queue.c)
add_executable(test_queue ${test_queue_FILES})
add_test(test_queue ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} "2"
         ${MPIEXEC_PREFLAGS} test_queue ${MPIEXEC_POSTFLAGS})

set(thesis_benchmark_FILES
	thesis_benchmark.c)
add_executable(thesis_benchmark ${thesis_benchmark_FILES})
add_test(thesis_benchmark ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} "2"
         ${MPIEXEC_PREFLAGS} thesis_benchmark ${MPIEXEC_POSTFLAGS}
         -ds 1.0 -zpp 2 -ppz 2)

