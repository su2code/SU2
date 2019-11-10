include(${CMAKE_CURRENT_LIST_DIR}/functions.cmake)

set(SU2_ENABLE_METIS
    OFF
    CACHE BOOL "Build SU2 with Metis support")

if(SU2_ENABLE_METIS)
  message(STATUS "<<< Configuring library with Metis support >>>")
  add_subdirectory(${CMAKE_SOURCE_DIR}/externals/metis)

  su2_add_dependency(METIS TARGETS Metis::Metis DEFINITIONS HAVE_METIS=1)
endif()
