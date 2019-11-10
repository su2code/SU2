project(tecio C CXX)

set(SU2_TECIO_CPPFLAGS
    "-fpermissive;USEENUM;THREED;MAKEARCHIVE;NO_ASSERTS;NO_THIRD_PARTY_LIBS"
    CACHE STRING "Specific Tecio compile definitions to use")

if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/boost")
  message(STATUS "Extracting boost source ...")
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xzf ${CMAKE_CURRENT_LIST_DIR}/boost.tar.gz
                          WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR} RESULT_VARIABLE RESULT)
  if(RESULT)
    message(
      FATAL_ERROR "Extraction of boost sources to ${CMAKE_CURRENT_LIST_DIR}/boost using 'tar' failed ...")
  endif()
endif()

set(TECIO_CPPFLAGS)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(TECIO_CPPFLAGS LINUX)
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    list(APPEND TECIO_CPPFLAGS LINUX64)
  endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(TECIO_CPPFLAGS DARWIN MAC64)
else()
  message(
    FATAL_ERROR
      ">>> Unrecognized TecIO platform ${CMAKE_SYSTEM_NAME}, see externals/tecio/Runmake for hints on how to extend <<<"
  )
endif()

if(SU2_ENABLE_MPI)
  add_subdirectory(teciompisrc)
else()
  add_subdirectory(teciosrc)
endif()
