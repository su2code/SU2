project(SU2_MSH CXX)

include(${CMAKE_SOURCE_DIR}/cmake/functions.cmake)

set(SU2_MSH_SOURCES include/SU2_MSH.hpp src/SU2_MSH.cpp)

add_executable(SU2_MSH ${SU2_MSH_SOURCES})
su2_setup_target(SU2_MSH LIBS SU2::SU2)
