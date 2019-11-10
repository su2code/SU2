include(${CMAKE_CURRENT_LIST_DIR}/functions.cmake)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)

if(Threads_FOUND)
  set(LIB_THREAD Threads::Threads)
else()
  set(LIB_THREAD)
endif()

su2_add_dependency(THREADS TARGETS "${LIB_THREAD}" EXPLICIT)
