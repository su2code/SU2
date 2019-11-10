include(${CMAKE_CURRENT_LIST_DIR}/functions.cmake)

set(SU2_ENABLE_CGNS
    OFF
    CACHE BOOL "Build SU2 with CGNS support")

if(SU2_ENABLE_CGNS)
  message(STATUS "<<< Configuring library with CGNS support >>>")
  add_subdirectory(${CMAKE_SOURCE_DIR}/externals/cgns)
  su2_add_dependency(CGNS TARGETS CGNS::CGNS DEFINITIONS HAVE_CGNS=1)
endif()
