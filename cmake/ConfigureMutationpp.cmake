include(${CMAKE_CURRENT_LIST_DIR}/functions.cmake)

set(SU2_ENABLE_MUTATIONPP
    OFF
    CACHE BOOL "Build SU2 with Mutationpp support")

if(SU2_ENABLE_MUTATIONPP)
  find_package(mutation++ REQUIRED)

  message(STATUS "<<< Configuring library with Mutationpp support >>>")
  su2_add_dependency(MUTATIONPP TARGETS mutation++::mutation++ DEFINITIONS HAVE_MUTATIONPP)
endif()
