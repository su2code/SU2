cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  if (${CMAKE_MAJOR_VERSION} GREATER 2 OR ${CMAKE_MINOR_VERSION} GREATER 7)
    cmake_policy(SET CMP0015 OLD)
  endif (${CMAKE_MAJOR_VERSION} GREATER 2 OR ${CMAKE_MINOR_VERSION} GREATER 7)
endif(COMMAND cmake_policy)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING
  "one of: Release, Debug, RelWithDebInfo or MinSizeRel")

project("pcgns")

# Allow for building a package
set(CPACK_PACKAGE_VERSION "0.2.0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "pcgns")
set(CPACK_GENERATOR "STGZ;TGZ")
set(CPACK_SOURCE_GENERATOR "STGZ;TGZ")
include(CPack)

set(CGNS_BUILD_SHARED "ON" CACHE BOOL "Build a shared version of the library")
option(CGNS_USE_SHARED "Link programs to the CGNS shared library" "ON")

#########
# Tests #
#########

option(ENABLE_TESTS "Enable test programs" "OFF")
if (ENABLE_TESTS)
  enable_testing()
endif (ENABLE_TESTS)

########
# HDF5 #
########

# Setup the variables for HDF5
find_path(HDF5_INCLUDE_PATH hdf5.h)
find_library(HDF5_LIBRARY hdf5)

set(HDF5_NEED_ZLIB "OFF" CACHE BOOL "Does the HDF5 library require linking to zlib?")
if(HDF5_NEED_ZLIB)
  find_library(ZLIB_LIBRARY z)
endif(HDF5_NEED_ZLIB)

set(HDF5_NEED_SZIP "OFF" CACHE BOOL "Does the HDF5 library require linking to zlib?")
if (HDF5_NEED_SZIP)
  find_library(SZIP_LIBRARY szip)
endif (HDF5_NEED_SZIP)

#######
# MPI #
#######

include(FindMPI)
if (NOT MPI_FOUND)
  mark_as_advanced(CLEAR MPI_COMPILER)
endif (NOT MPI_FOUND)

if (ENABLE_TESTS AND NOT MPIEXEC)
  mark_as_advanced(CLEAR MPIEXEC)
endif (ENABLE_TESTS AND NOT MPIEXEC)

########
# CGNS #
########

# Include the src directory
add_subdirectory(src)
