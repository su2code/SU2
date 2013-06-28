cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  if (${CMAKE_MAJOR_VERSION} GREATER 2 OR ${CMAKE_MINOR_VERSION} GREATER 7)
    cmake_policy(SET CMP0015 OLD)
  endif (${CMAKE_MAJOR_VERSION} GREATER 2 OR ${CMAKE_MINOR_VERSION} GREATER 7)
endif(COMMAND cmake_policy)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING
  "one of: Release, Debug, RelWithDebInfo or MinSizeRel")
if (${CMAKE_MAJOR_VERSION} GREATER 2 OR ${CMAKE_MINOR_VERSION} GREATER 7)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release;Debug;RelWithDebInfo;MinSizeRel")
endif (${CMAKE_MAJOR_VERSION} GREATER 2 OR ${CMAKE_MINOR_VERSION} GREATER 7)

project("cgns" C)
set(CGNS_VERSION "3.1")

# Allow for building a package
set(CPACK_PACKAGE_VERSION "${CGNS_VERSION}-develop")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "cgns-${CGNS_VERSION}")
set(CPACK_GENERATOR "STGZ;TGZ;ZIP")
set(CPACK_SOURCE_GENERATOR "STGZ;TGZ;ZIP")
include(CPack)

CONFIGURE_FILE(
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
	IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
	"${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

# if this is cygwin, turn off win32 flag
if (CYGWIN)
  set(WIN32 "NO")
endif (CYGWIN)

##############
# build mode #
##############

# get some type sizes
include(CheckTypeSize)
CHECK_TYPE_SIZE(long SIZE_OF_LONG)
CHECK_TYPE_SIZE(off_t SIZE_OF_OFF_T)

if ("${SIZE_OF_OFF_T}" LESS 8)
  option(ENABLE_LFS "Enable or disable large file support" "OFF")
endif ("${SIZE_OF_OFF_T}" LESS 8)

if ("${CMAKE_SIZEOF_VOID_P}" GREATER 4)
  option(ENABLE_64BIT "Enable or disable 64-bit code building" "OFF")
else ("${CMAKE_SIZEOF_VOID_P}" GREATER 4)
  set(ENABLE_64BIT "OFF")
endif ("${CMAKE_SIZEOF_VOID_P}" GREATER 4)

option(ENABLE_LEGACY "Enable or disable building legacy code (3.0 compatible)" "OFF")
option(ENABLE_SCOPING "Enable or disable scoping of enumeration values" "OFF")

if (ENABLE_LEGACY)
  set(ENABLE_64BIT "OFF")
  mark_as_advanced(FORCE ENABLE_64BIT)
else (ENABLE_LEGACY)
  mark_as_advanced(CLEAR ENABLE_64BIT)
endif (ENABLE_LEGACY)

if (ENABLE_64BIT)
  set(ENABLE_LEGACY "OFF")
  mark_as_advanced(FORCE ENABLE_LEGACY)
else (ENABLE_64BIT)
  mark_as_advanced(CLEAR ENABLE_LEGACY)
endif (ENABLE_64BIT)

if (ENABLE_LFS)
  if (WIN32)
    add_definitions(-DHAVE_LSEEK64 -DHAVE_LSEEK64)
  else (WIN32)
    include(CheckFunctionExists)
    CHECK_FUNCTION_EXISTS(open64 HAVE_OPEN64)
    CHECK_FUNCTION_EXISTS(lseek64 HAVE_LSEEK64)
    if (HAVE_OPEN64)
      add_definitions(-DHAVE_LSEEK64)
    endif (HAVE_OPEN64)
    if (HAVE_OPEN64)
      add_definitions(-DHAVE_LSEEK64)
    endif (HAVE_LSEEK64)
  endif (WIN32)
else (ENABLE_LFS)
  remove_definitions(-DHAVE_OPEN64 -DHAVE_LSEEK64)
endif (ENABLE_LFS)

option(ENABLE_TESTS "Enable test programs" "OFF")
if (ENABLE_TESTS)
  enable_testing()
endif (ENABLE_TESTS)

###########
# Fortran #
###########

# Allow fortran support to be toggled
option(ENABLE_FORTRAN "Enable or disable the use of Fortran" "OFF")
if (NOT DEFINED FORTRAN_NAMING)
  set(FORTRAN_NAMING "UNKNOWN")
endif (NOT DEFINED FORTRAN_NAMING)
set(F2CLIST  "LOWERCASE" "LOWERCASE_" "LOWERCASE__" "UPPERCASE" "UPPERCASE_" "UPPERCASE__")
set(FORTRAN_NAMING_HELP
"Configures how to link the Fortran components into the C library.
Options are:
	LOWERCASE
	LOWERCASE_
	LOWERCASE__
	UPPERCASE
	UPPERRCASE_
	UPPERRCASE__
LOWERCASE_ is known to work with gfortran.  If this setting
does not work with your compiler, try the others until one does.
")

# If we are using fortran, enable it and configure the C-Fortran interface
# It would be nice to automatically detect this setting

if(ENABLE_FORTRAN)
  mark_as_advanced(CLEAR FORTRAN_NAMING)
  #only do this once
  if (FORTRAN_NAMING STREQUAL "UNKNOWN")
    set(LIST  "LOWERCASE" "LOWERCASE_" "LOWERCASE__" "UPPERCASE" "UPPERCASE_" "UPPERCASE__")
    enable_language("Fortran" OPTIONAL)
    # if no Fortran compiler, skip this
    if (CMAKE_Fortran_COMPILER_WORKS)
      foreach(nc ${F2CLIST})
        try_compile(NAMING_TEST_RESULT
	    "${CMAKE_BINARY_DIR}/fortran_test"
	    "${CMAKE_SOURCE_DIR}/fortran_test"
	    projectName ${nc})
        if(${NAMING_TEST_RESULT} STREQUAL "TRUE")
          set(FORTRAN_NAMING "${nc}" CACHE STRING ${FORTRAN_NAMING_HELP})
          if (${CMAKE_MAJOR_VERSION} GREATER 2 OR ${CMAKE_MINOR_VERSION} GREATER 7)
            set_property(CACHE FORTRAN_NAMING PROPERTY STRINGS ${F2CLIST})
          endif (${CMAKE_MAJOR_VERSION} GREATER 2 OR ${CMAKE_MINOR_VERSION} GREATER 7)
          message(STATUS "Trying Fortran naming convention ${FORTRAN_NAMING}: success")
          message(STATUS "Setting fortran naming convention to ${FORTRAN_NAMING}")
          break()
        else(${NAMING_TEST_RESULT} STREQUAL "TRUE")
          message(STATUS "Trying Fortran naming convention ${nc}: failure")
        endif(${NAMING_TEST_RESULT} STREQUAL "TRUE")
      endforeach(nc ${LIST})
    endif (CMAKE_Fortran_COMPILER_WORKS)
  endif (FORTRAN_NAMING STREQUAL "UNKNOWN")
else(ENABLE_FORTRAN)
  mark_as_advanced(FORCE FORTRAN_NAMING)
endif(ENABLE_FORTRAN)

# check that FORTRAN_NAMING is valid
if (ENABLE_FORTRAN AND FORTRAN_NAMING)
  if (FORTRAN_NAMING STREQUAL "UNKNOWN")
    set(FORTRAN_NAMING "" CACHE STRING ${FORTRAN_NAMING_HELP})
    if (${CMAKE_MAJOR_VERSION} GREATER 2 OR ${CMAKE_MINOR_VERSION} GREATER 7)
      set_property(CACHE FORTRAN_NAMING PROPERTY STRINGS ${F2CLIST})
    endif (${CMAKE_MAJOR_VERSION} GREATER 2 OR ${CMAKE_MINOR_VERSION} GREATER 7)
  else (FORTRAN_NAMING STREQUAL "UNKNOWN")
    list(FIND F2CLIST ${FORTRAN_NAMING} ListIndex)
    if (ListIndex LESS 0)
      message(SEND_ERROR "invalid FORTRAN_NAMING value. Must be empty or one of:
        LOWERCASE
        LOWERCASE_
        LOWERCASE__
        UPPERCASE
        UPPERRCASE_
        UPPERRCASE__")
    endif (ListIndex LESS 0)
  endif (FORTRAN_NAMING STREQUAL "UNKNOWN")
endif (ENABLE_FORTRAN AND FORTRAN_NAMING)

# this forces cmake to set up the required Fortran variables
if (ENABLE_FORTRAN AND CMAKE_Fortran_COMPILER_WORKS)
  enable_language("Fortran")
endif (ENABLE_FORTRAN AND CMAKE_Fortran_COMPILER_WORKS)

########
# HDF5 #
########

option(ENABLE_HDF5 "Enable or disable HDF5 interface" "OFF")
if (ENABLE_HDF5)
  find_path(HDF5_INCLUDE_PATH hdf5.h)
  find_library(HDF5_LIBRARY hdf5)

  set(HDF5_NEED_ZLIB "OFF" CACHE BOOL "Does the HDF5 library require linking to zlib?")
  if(HDF5_NEED_ZLIB)
    find_library(ZLIB_LIBRARY z)
    mark_as_advanced(CLEAR ZLIB_LIBRARY)
  else (HDF5_NEED_ZLIB)
    mark_as_advanced(FORCE ZLIB_LIBRARY)
  endif(HDF5_NEED_ZLIB)

  set(HDF5_NEED_SZIP "OFF" CACHE BOOL "Does the HDF5 library require linking to zlib?")
  if (HDF5_NEED_SZIP)
    find_library(SZIP_LIBRARY szip)
    mark_as_advanced(CLEAR SZIP_LIBRARY)
  else (HDF5_NEED_SZIP)
    mark_as_advanced(FORCE SZIP_LIBRARY)
  endif (HDF5_NEED_SZIP)

  set(HDF5_NEED_MPI "OFF" CACHE BOOL "Does the HDF5 library require linking to mpi? (Only true if using parallel HDF5)")
  set(MPI_INC)
  set(MPI_LIBS)
  if (HDF5_NEED_MPI)
    if (NOT MPI_FOUND)
      include(FindMPI)
    endif (NOT MPI_FOUND)
    # MPI variable names changed with cmake 2.8.5
    if (DEFINED MPI_C_COMPILER)
      mark_as_advanced(CLEAR MPI_C_COMPILER MPI_C_INCLUDE_PATH MPI_C_LIBRARIES)
      mark_as_advanced(FORCE MPI_LIBRARY MPI_EXTRA_LIBRARY)
      set(MPI_INC ${MPI_C_INCLUDE_PATH})
      set(MPI_LIBS ${MPI_C_LIBRARIES})
    else (DEFINED MPI_C_COMPILER)
      mark_as_advanced(CLEAR MPI_COMPILER MPI_INCLUDE_PATH MPI_LIBRARY MPI_EXTRA_LIBRARY)
      set(MPI_INC ${MPI_INCLUDE_PATH})
      if (MPI_EXTRA_LIBRARY)
        set(MPI_LIBS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
      else (MPI_EXTRA_LIBRARY)
        set(MPI_LIBS ${MPI_LIBRARY})
      endif (MPI_EXTRA_LIBRARY)
    endif (DEFINED MPI_C_COMPILER)
    if (ENABLE_TESTS AND ENABLE_PARALLEL)
      mark_as_advanced(CLEAR MPIEXEC)
    else (ENABLE_TESTS AND ENABLE_PARALLEL)
      mark_as_advanced(FORCE MPIEXEC)
    endif (ENABLE_TESTS AND ENABLE_PARALLEL)
  else (HDF5_NEED_MPI)
    mark_as_advanced(FORCE MPI_C_COMPILER MPI_C_INCLUDE_PATH MPI_C_LIBRARIES)
    mark_as_advanced(FORCE MPI_COMPILER MPI_INCLUDE_PATH MPI_LIBRARY MPI_EXTRA_LIBRARY MPIEXEC)
  endif (HDF5_NEED_MPI)
  mark_as_advanced(CLEAR HDF5_NEED_ZLIB HDF5_NEED_SZIP HDF5_NEED_MPI)
else (ENABLE_HDF5)
  mark_as_advanced(FORCE HDF5_NEED_ZLIB HDF5_NEED_SZIP HDF5_NEED_MPI)
  mark_as_advanced(FORCE HDF5_INCLUDE_PATH HDF5_LIBRARY)
  mark_as_advanced(FORCE ZLIB_LIBRARY SZIP_LIBRARY)
  mark_as_advanced(FORCE MPI_C_COMPILER MPI_C_INCLUDE_PATH MPI_C_LIBRARIES)
  mark_as_advanced(FORCE MPI_COMPILER MPI_INCLUDE_PATH MPI_LIBRARY MPI_EXTRA_LIBRARY MPIEXEC)
endif (ENABLE_HDF5)

#################
# parallel cgns #
#################

if (ENABLE_HDF5 AND HDF5_NEED_MPI)
  set(BUILD_PCGNS "OFF" CACHE BOOL "Build parallel CGNS ?")
  mark_as_advanced(CLEAR BUILD_PCGNS)
else (ENABLE_HDF5 AND HDF5_NEED_MPI)
  set(BUILD_PCGNS "OFF")
  mark_as_advanced(FORCE BUILD_PCGNS)
endif (ENABLE_HDF5 AND HDF5_NEED_MPI)

if (ENABLE_TESTS AND BUILD_PCGNS)
  mark_as_advanced(CLEAR MPIEXEC)
else (ENABLE_TESTS AND BUILD_PCGNS)
  mark_as_advanced(FORCE MPIEXEC)
endif (ENABLE_TESTS AND BUILD_PCGNS)

####################
# RPATH Management #
####################

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

########
# CGNS #
########

# Include the src directory
add_subdirectory(src)

#########
# PCGNS #
#########

if (BUILD_PCGNS)
  add_subdirectory(pcgns/src)
endif (BUILD_PCGNS)

