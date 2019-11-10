project(SU2_SOL CXX)

if(NOT SU2_BUILD_CFD)
  message(FATAL_ERROR "SU2_SOL requires SU2_CFD to be built")
endif()

include(${CMAKE_SOURCE_DIR}/cmake/functions.cmake)

set(SU2_SOL_SOURCES include/SU2_SOL.hpp src/SU2_SOL.cpp)

add_executable(SU2_SOL ${SU2_SOL_SOURCES})
su2_setup_target(SU2_SOL LIBS SU2::SU2Core SU2::SU2 WITH_THREADS)
