project(SU2_GEO CXX)

include(${CMAKE_SOURCE_DIR}/cmake/functions.cmake)


set(SU2_GEO_SOURCES include/SU2_GEO.hpp src/SU2_GEO.cpp)

add_executable(SU2_GEO ${SU2_GEO_SOURCES})
su2_setup_target(SU2_GEO LIBS SU2::SU2)
