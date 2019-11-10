include(${CMAKE_CURRENT_LIST_DIR}/functions.cmake)

set(SU2_ENABLE_TECIO
    OFF
    CACHE BOOL "Build SU2 with TecIO support")

if(SU2_ENABLE_TECIO)
  message(STATUS "<<< Configuring library with Tecplot TecIO support >>>")
  add_subdirectory(${CMAKE_SOURCE_DIR}/externals/tecio)
  su2_add_dependency(
    TECIO
    TARGETS
    Tecio::Tecio
    DEFINITIONS
    HAVE_TECIO=1
    HAVE_TECPLOT_API=1
    HAVE_TECPLOT_API_112=1)
endif()
