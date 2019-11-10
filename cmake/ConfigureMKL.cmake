include(${CMAKE_CURRENT_LIST_DIR}/functions.cmake)

set(SU2_ENABLE_MKL
    OFF
    CACHE BOOL "Build SU2 with Intel MKL support")

if(SU2_ENABLE_MKL)
  set(REQUIRED_MKL_VERSION 20190000)
  find_package(MKL REQUIRED)
  if(MKL_VERSION LESS REQUIRED_MKL_VERSION)
    message(FATAL_ERROR "Unsupported MKL version ${MKL_VERSION}, ${REQUIRED_MKL_VERSION} or greater required")
  endif()

  message(STATUS "<<< Configuring library with Intel MKL support >>>")
  su2_add_dependency(MKL TARGETS ${MKL_LIBRARIES} DEFINITIONS HAVE_MKL MKL_DIRECT_CALL_SEQ)
endif()
