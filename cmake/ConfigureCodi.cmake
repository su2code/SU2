include(${CMAKE_CURRENT_LIST_DIR}/functions.cmake)

set(SU2_ENABLE_CODI
    "no"
    CACHE STRING "Build executables with codi datatype")
set_property(CACHE SU2_ENABLE_CODI PROPERTY STRINGS "no" "reverse" "forward")
set(SU2_CODI_CPPFLAGS
    ""
    CACHE STRING "Specific Codi C Preprocessor flags to use")

if(SU2_ENABLE_CODI STREQUAL "forward")
  set(SU2_BUILD_NORMAL
      OFF
      CACHE INTERNAL "")
  set(SU2_BUILD_DIRECTDIFF
      ON
      CACHE INTERNAL "")
  set(SU2_BUILD_REVERSE
      OFF
      CACHE INTERNAL "")
elseif(SU2_ENABLE_CODI STREQUAL "reverse")
  set(SU2_BUILD_NORMAL
      OFF
      CACHE INTERNAL "")
  set(SU2_BUILD_DIRECTDIFF
      OFF
      CACHE INTERNAL "")
  set(SU2_BUILD_REVERSE
      ON
      CACHE INTERNAL "")
elseif(SU2_ENABLE_CODI STREQUAL "no")
  set(SU2_BUILD_NORMAL
      ON
      CACHE INTERNAL "")
  set(SU2_BUILD_DIRECTDIFF
      OFF
      CACHE INTERNAL "")
  set(SU2_BUILD_REVERSE
      OFF
      CACHE INTERNAL "")
else()
  message(FATAL_ERROR "Invalid SU2_ENABLE_CODI value ${SU2_ENABLE_CODI}, expected one of no, reverse, forward")
endif()

if(SU2_BUILD_DIRECTDIFF OR SU2_BUILD_REVERSE)
  message(STATUS "<<< Configuring library with Codi ${SU2_ENABLE_CODI} support >>>")

  set(CODIheader ${CMAKE_SOURCE_DIR}/externals/codi/include/codi.hpp)
  set(AMPIheader ${CMAKE_SOURCE_DIR}/externals/medi/include/medi/medi.hpp)

  set(sha_version_codi 501dcf0305df147481630f20ce37c2e624fb351f)
  set(github_repo_codi https://github.com/scicompkl/CoDiPack)
  set(sha_version_medi a95a23ce7585905c3a731b28c1bb512028fc02bb)
  set(github_repo_medi https://github.com/SciCompKL/MeDiPack)

  set(medi_name MeDiPack)
  set(codi_name CoDiPack)

  set(alt_name_medi externals/medi)
  set(alt_name_codi externals/codi)

  if(NOT EXISTS ${CODIheader})
    su2_download_module(${codi_name} ${alt_name_codi} ${github_repo_codi} ${sha_version_codi})
  endif()
  set(INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/externals/codi/include")

  if(SU2_ENABLE_MPI)
    if(NOT EXISTS ${AMPIheader})
      su2_download_module(${medi_name} ${alt_name_medi} ${github_repo_medi} ${sha_version_medi})
    endif()
    list(APPEND INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/externals/medi/include" "${CMAKE_SOURCE_DIR}/externals/medi/src")
  endif()

  if(SU2_BUILD_DIRECTDIFF)
    list(APPEND CODI_DEFINITIONS "CODI_FORWARD_TYPE")
  else()
    list(APPEND CODI_DEFINITIONS "CODI_REVERSE_TYPE")
  endif()

  # create a header-only target that can be linked against since compiling codi is not required
  add_library(codi INTERFACE)
  target_include_directories(codi INTERFACE ${INCLUDE_DIRS})
  target_compile_definitions(codi INTERFACE ${CODI_DEFINITIONS})
  add_library(Codi::Codi ALIAS codi)

  su2_add_dependency(CODI TARGETS Codi::Codi DEFINITIONS ${CODI_DEFINITIONS} EXPLICIT)
endif()
