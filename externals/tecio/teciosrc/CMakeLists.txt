cmake_minimum_required(VERSION 2.8.12) # requiring 2.8.12 sets the CMP0022 policy to new (non-transitive linking, i.e. just the interfaces) by default

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING "CMake build type (Debug|Release)" FORCE)
endif ()

project (tecio C CXX)

FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ELSE()
    Set (Boost_INCLUDE_DIR /3rdpartylibs/boost/1.59.0/build/linux64-centos6.5/include)
    message ("Warning:  Cannot find boost include directories.  Hardcoding to ${Boost_INCLUDE_DIR}")
    INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIR})
ENDIF()
include_directories(.)

IF (WIN32)
    set(BaseFlags "/EHsc /MP /wd\"4996\" /D MSWIN /D TP_PROJECT_USES_BOOST /D BOOST_ALL_NO_LIB /D MAKEARCHIVE /D NO_THIRD_PARTY_LIBS /D NO_ASSERTS")
ELSE ()
	set(BaseFlags "-DLINUX -DLINUX64 -DTP_PROJECT_USES_BOOST -DBOOST_ALL_NO_LIB -DMAKEARCHIVE -DNO_THIRD_PARTY_LIBS -DNO_ASSERTS -fmessage-length=0 -fPIC -fvisibility=hidden -w")
ENDIF ()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BaseFlags}")
                    
FILE(GLOB mainFiles "*.cpp" "*.h")
LIST(REMOVE_ITEM mainFiles "${CMAKE_CURRENT_SOURCE_DIR}/szcombine.cpp")

add_library(tecio STATIC ${mainFiles})

add_executable(szcombine "szcombine.cpp")
if (WIN32)
    target_link_libraries(szcombine tecio)
else ()
    target_link_libraries(szcombine tecio pthread)
endif ()

# Had to add this so building tecio care package within loci-chem works...
set_target_properties(tecio PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(szcombine PROPERTIES LINKER_LANGUAGE CXX)
