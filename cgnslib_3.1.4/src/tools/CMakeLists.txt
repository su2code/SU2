#########
# Tools #
#########

# Link all the executables to cgns and hdf5
link_directories(.)
if (CGNS_USE_SHARED)
  link_libraries(cgns_shared)
else (CGNS_USE_SHARED)
  link_libraries(cgns_static)
endif (CGNS_USE_SHARED)

if (ENABLE_HDF5 AND HDF5_LIBRARY)
  link_libraries(${HDF5_LIBRARY})
  if(HDF5_NEED_ZLIB AND ZLIB_LIBRARY)
    link_libraries(${ZLIB_LIBRARY})
  endif(HDF5_NEED_ZLIB AND ZLIB_LIBRARY)
  if(HDF5_NEED_SZIP AND SZIP_LIBRARY)
    link_libraries(${SZIP_LIBRARY})
  endif(HDF5_NEED_SZIP AND SZIP_LIBRARY)
  if(HDF5_NEED_MPI AND MPI_LIBS)
    link_libraries(${MPI_LIBS})
  endif(HDF5_NEED_MPI AND MPI_LIBS)
endif (ENABLE_HDF5 AND HDF5_LIBRARY)

if (NOT WIN32)
  link_libraries(m)
endif (NOT WIN32)

# Set the files needed by each tool
set(cgnsnames_FILES
	cgnsnames.c
	cgnames.c)
set(cgnscheck_FILES
	cgnscheck.c
	getargs.c
	hash.c
	cgnames.c)
set(cgnsconvert_FILES
	cgnsconvert.c
	getargs.c)
set(cgnsdiff_FILES
	cgnsdiff.c
	getargs.c)
set(cgnslist_FILES
	cgnslist.c
	getargs.c)
set(cgnscompress_FILES
	cgnscompress.c)
set(cgnsversion_FILES
	cgnsversion.c
	getargs.c)

# Build each tool
add_executable( cgnsnames    ${cgnsnames_FILES}     )
add_executable( cgnscheck    ${cgnscheck_FILES}   )
add_executable( cgnsconvert  ${cgnsconvert_FILES} )
add_executable( cgnsdiff     ${cgnsdiff_FILES}    )
add_executable( cgnslist     ${cgnslist_FILES}    )
add_executable( cgnscompress ${cgnscompress_FILES} )
add_executable( cgnsversion  ${cgnsversion_FILES} )

install(TARGETS
	cgnsnames
	cgnscheck
	cgnsconvert
	cgnsdiff
	cgnslist
        cgnscompress
	cgnsversion
	RUNTIME DESTINATION bin)

if (WIN32)
  install(PROGRAMS
        cgnsupdate.bat
	DESTINATION bin)
  # lets also put the dll in here so they run
  if (CGNS_USE_SHARED)
    install(PROGRAMS
        ${CMAKE_BINARY_DIR}/src/cgnsdll.dll
	DESTINATION bin)
  endif (CGNS_USE_SHARED)
  if (ENABLE_HDF5)
    install(PROGRAMS
        adf2hdf.bat
        hdf2adf.bat
	DESTINATION bin)
  endif (ENABLE_HDF5)
else (WIN32)
  install(PROGRAMS
        cgnsupdate
	DESTINATION bin)
  if (ENABLE_HDF5)
    install(PROGRAMS
        adf2hdf
        hdf2adf
	DESTINATION bin)
  endif (ENABLE_HDF5)
endif (WIN32)

