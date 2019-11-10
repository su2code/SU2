include(${CMAKE_CURRENT_LIST_DIR}/functions.cmake)

set(PARMETIS_DOCSTRING "Build SU2 with Parmetis support")
su2_option(SU2_ENABLE_PARMETIS SETUP OFF OFF BOOL ${PARMETIS_DOCSTRING})

# hide option if MPI is not used
if(NOT SU2_ENABLE_MPI)
  su2_option(SU2_ENABLE_PARMETIS HIDE)
endif()

if(SU2_ENABLE_PARMETIS)
  message(STATUS "<<< Configuring library with Parmetis support >>>")
  add_subdirectory(${CMAKE_SOURCE_DIR}/externals/parmetis)
  su2_add_dependency(PARMETIS TARGETS Parmetis::Parmetis DEFINITIONS HAVE_PARMETIS=1)
endif()
