project(SU2_DOT CXX)

if(NOT SU2_BUILD_CFD)
  message(FATAL_ERROR "SU2_DOT requires SU2_CFD to be built")
endif()

include(${CMAKE_SOURCE_DIR}/cmake/functions.cmake)

set(SU2_DOT_SOURCES include/SU2_DOT.hpp src/SU2_DOT.cpp)

if(SU2_BUILD_NORMAL)
  set(SU2_DOT_EXECUTABLE SU2_DOT)
elseif(SU2_BUILD_REVERSE)
  set(SU2_DOT_EXECUTABLE SU2_DOT_AD)
else()
  message(FATAL_ERROR "SU2_DOT can only be built in normal or reverse mode")
endif()

add_executable(${SU2_DOT_EXECUTABLE} ${SU2_DOT_SOURCES})

if(SU2_BUILD_NORMAL)
  su2_setup_target(${SU2_DOT_EXECUTABLE} LIBS SU2::SU2Core SU2::SU2 WITH_THREADS)
else()
  su2_setup_target(${SU2_DOT_EXECUTABLE} LIBS SU2::SU2Core SU2::SU2 WITH_THREADS WITH_CODI)
endif()
