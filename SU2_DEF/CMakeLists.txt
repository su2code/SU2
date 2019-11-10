project(SU2_DEF CXX)

if(NOT SU2_BUILD_CFD)
  message(FATAL_ERROR "SU2_DEF requires SU2_CFD to be built")
endif()

include(${CMAKE_SOURCE_DIR}/cmake/functions.cmake)


set(SU2_DEF_SOURCES include/SU2_DEF.hpp src/SU2_DEF.cpp)

add_executable(SU2_DEF ${SU2_DEF_SOURCES})
su2_setup_target(SU2_DEF LIBS SU2::SU2Core SU2::SU2 WITH_THREADS)
