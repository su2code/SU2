##########
# Global #
##########

# Include directories
include_directories(. ${HDF5_INCLUDE_PATH} ${MPI_INCLUDE_PATH})

###########
# Library #
###########

# All the files needed by the cgns library
set(pcgns_FILES
	pcgnslib.c
	pcgns_util.c)

# Build a static version of the library
add_library(pcgns_static STATIC ${pcgns_FILES})
set_target_properties(pcgns_static PROPERTIES OUTPUT_NAME pcgns)
set_target_properties(pcgns_static PROPERTIES CLEAN_DIRECT_OUTPUT 1)
set_target_properties(pcgns_static PROPERTIES VERSION "0.2")
set_target_properties(pcgns_static PROPERTIES SOVERSION "0.2")

# Build a shared version of the library
if(CGNS_BUILD_SHARED)
  mark_as_advanced(CLEAR CGNS_USE_SHARED)
  add_library(pcgns_shared SHARED ${pcgns_FILES})
  # under windows we need to link with the HDF5 and MPI libraries
  # to build a dll and change the name of the shared library
  # for both static and shared version to exist
  if (WIN32 OR CYGWIN)
    target_link_libraries(pcgns_shared ${HDF5_LIBRARY})
    if(HDF5_NEED_ZLIB AND ZLIB_LIBRARY)
      target_link_libraries(pcgns_shared ${ZLIB_LIBRARY})
    endif(HDF5_NEED_ZLIB AND ZLIB_LIBRARY)
    if(HDF5_NEED_SZIP AND SZIP_LIBRARY)
      target_link_libraries(pcgns_shared ${SZIP_LIBRARY})
    endif(HDF5_NEED_SZIP AND SZIP_LIBRARY)
    if(MPI_LIBRARY)
      target_link_libraries(pcgns_shared ${MPI_LIBRARY})
    endif(MPI_LIBRARY)
    set_target_properties(pcgns_shared PROPERTIES OUTPUT_NAME pcgnsdll)
  else (WIN32 OR CYGWIN)
    set_target_properties(pcgns_shared PROPERTIES OUTPUT_NAME pcgns)
  endif (WIN32 OR CYGWIN)
  set_target_properties(pcgns_shared PROPERTIES CLEAN_DIRECT_OUTPUT 1)
  set_target_properties(pcgns_shared PROPERTIES VERSION "0.2")
  set_target_properties(pcgns_shared PROPERTIES SOVERSION "0.2")
else(CGNS_BUILD_SHARED)
  set(CGNS_USE_SHARED "OFF")
  mark_as_advanced(FORCE CGNS_USE_SHARED)
endif(CGNS_BUILD_SHARED)

# Set the install path of the static library
install(TARGETS pcgns_static ARCHIVE DESTINATION lib)

# Set the install path of the shared library
if(CGNS_BUILD_SHARED)
  # for windows, need to install both cgnsdll.dll and cgnsdll.lib
  if (WIN32 OR CYGWIN)
    install(TARGETS pcgns_shared ARCHIVE DESTINATION lib)
    install(TARGETS pcgns_shared RUNTIME DESTINATION lib)
  else (WIN32 OR CYGWIN)
    install(TARGETS pcgns_shared LIBRARY DESTINATION lib)
  endif (WIN32 OR CYGWIN)
endif(CGNS_BUILD_SHARED)

# Set the install path of the header files
install(FILES
	pcgnslib.h
	DESTINATION include)

#########
# Tests #
#########

if (ENABLE_TESTS)
  add_subdirectory(tests)
endif (ENABLE_TESTS)

