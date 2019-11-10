project(pySU2)

cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)

include(${CMAKE_SOURCE_DIR}/cmake/functions.cmake)

find_package(Python REQUIRED COMPONENTS Development)
find_package(SWIG REQUIRED)
include(UseSWIG)

if(SU2_ENABLE_MPI)
  find_package(MPI4PY REQUIRED)
else()
  set(MPI4PY_INCLUDE_DIR "")
endif()

set(LINK_LIBS SU2::SU2Core SU2::SU2 Python::Python)
if(SU2_BUILD_NORMAL)
  set(SWIG_SRC ${CMAKE_CURRENT_LIST_DIR}/pySU2.i)
  set(pysu2_LIB pysu2)
elseif(SU2_BUILD_REVERSE)
  set(SWIG_SRC ${CMAKE_CURRENT_LIST_DIR}/pySU2ad.i)
  set(pysu2_LIB pysu2ad)
  list(APPEND LINK_LIBS WITH_CODI)
else()
  message(FATAL_ERROR "pySU2 can only be built in normal or reverse mode")
endif()

# compile SWIG wrapper with C++
set_property(SOURCE ${SWIG_SRC} PROPERTY CPLUSPLUS ON)
set_property(SOURCE ${SWIG_SRC} PROPERTY SWIG_MODULE_NAME ${pysu2_LIB})

# need to pull include directories and compile definitions manually
# since this is not a real cmake target
su2_get_compile_options(
  INCLUDES I
  DEFINITIONS D
  LIBS SU2::SU2Core SU2::SU2)
set_property(SOURCE ${SWIG_SRC} PROPERTY COMPILE_DEFINITIONS ${D})
set_property(SOURCE ${SWIG_SRC} PROPERTY INCLUDE_DIRECTORIES ${I} ${MPI4PY_INCLUDE_DIR})

if(SU2_BUILD_REVERSE)
  # need this to make 'using namespace ...;' work https://github.com/swig/swig/issues/1068
  # though compile times seem to increase by a lot
  set_property(SOURCE ${SWIG_SRC} PROPERTY COMPILE_OPTIONS -includeall)
endif()

swig_add_library(
  ${pysu2_LIB}
  TYPE SHARED
  LANGUAGE python
  SOURCES ${SWIG_SRC})

if(SU2_ENABLE_MPI)
  list(APPEND LINK_LIBS mpi4py::mpi4py)
endif()

su2_setup_target(${pysu2_LIB} LIBS ${LINK_LIBS})
add_library(SU2::pySU2 ALIAS ${pysu2_LIB})

install(FILES "$<TARGET_FILE_DIR:${pysu2_LIB}>/${pysu2_LIB}.py" DESTINATION bin)
