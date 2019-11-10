# don't check if already done
if(TLS_CHECKED)
  return()
endif()

# persist through reconfigures
set(TLS_CHECKED
    ON
    CACHE INTERNAL "")

# Custom check for TLS.
if(MSVC)
  add_definitions(-D__thread=__declspec\(thread\))
else()
  # This if checks if that value is cached or not.
  try_compile(HAVE_THREADLOCALSTORAGE ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_LIST_DIR}/check_thread_storage.c)
  if(HAVE_THREADLOCALSTORAGE)
    message(STATUS "checking for thread-local storage - found")
    include(${CMAKE_CURRENT_LIST_DIR}/functions.cmake)
    add_definitions(-DTLS)
  else()
    message(STATUS "checking for thread-local storage - not found")
  endif()

  if(NOT HAVE_THREADLOCALSTORAGE)
    add_definitions(-D__thread=)
  endif()
endif()
