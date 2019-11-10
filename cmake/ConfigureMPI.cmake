include(${CMAKE_CURRENT_LIST_DIR}/functions.cmake)

set(SU2_ENABLE_MPI
    OFF
    CACHE BOOL "Build SU2 with MPI support")

if(SU2_ENABLE_MPI)
  message(STATUS "<<< Configuring library with MPI support >>>")
  find_package(MPI REQUIRED)

  su2_add_dependency(MPI TARGETS MPI::MPI_CXX DEFINITIONS HAVE_MPI)
endif()
