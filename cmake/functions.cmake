# include guard
if(INCLUDED)
  return()
endif()

set(INCLUDED 1)

#[[
output additional messages if DEBUG is defined and TRUE
]]
function(su2_debug)
  if(DEBUG)
    message(STATUS "--Debug: ${ARGN}")
  endif()
endfunction()

#[[
merge list B into A without duplicates
]]
function(su2_merge_lists A B)
  if(NOT ${A})
    # list A is empty, return B
    set(${A} ${${B}} PARENT_SCOPE)
  elseif(${B})
    # if B is not empty, remove duplicates
    list(REMOVE_ITEM ${B} ${${A}})
    list(APPEND ${A} ${${B}}) # append
    set(${A} ${${A}} PARENT_SCOPE) # return
  endif()
endfunction()

#[[
su2_add_dependency(<NAME> [EXPLICIT]
                  [TARGETS ...]
                  [DEFINITIONS ...]
                  [INCLUDES ...])

setup dependency to SU2 that will be used when building SU2 components
NAME: name of the dependency to be used in cmake variables
EXPLICIT: whether using this dependency is explicit, i.e. requires su2_setup_target called with LIBS WITH_<NAME>
TARGETS: list of targets to link against
DEFINITIONS: list of compile definitions to use for the build target
INCLUDES: list of include directories to use for the build target

This function defines:
SU2_DEPENDENCY_<NAME>_INCLUDES: include directories merged from INCLUDES and INTERFACE_INCLUDE_DIRECTORIES of TARGETS
SU2_DEPENDENCY_<NAME>_DEFINITIONS: compile definitions merged from DEFINITIONS and INTERFACE_COMPILE_DEFINITIONS of TARGETS
SU2_DEPENDENCY_<NAME>_TARGETS: TARGETS
SU2_DEPENDENCY_<NAME>_EXPLICIT: EXPLICIT
]]
function(su2_add_dependency NAME)
  set(options EXPLICIT)
  set(oneValueArgs)
  set(multiValueArgs TARGETS DEFINITIONS INCLUDES)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT ARG_TARGETS)
    set(ARG_TARGETS "")
  endif()

  if(NOT ARG_DEFINITIONS)
    set(ARG_DEFINITIONS "")
  endif()

  set(includes)
  set(definitions ${ARG_DEFINITIONS})
  foreach(library IN LISTS ARG_TARGETS)
    if(TARGET ${library})
      get_target_property(INCLUDES ${library} INTERFACE_INCLUDE_DIRECTORIES)
      get_target_property(DEFINITIONS ${library} INTERFACE_COMPILE_DEFINITIONS)
      if(INCLUDES)
        su2_merge_lists(includes INCLUDES)
      endif()

      if(DEFINITIONS)
        su2_merge_lists(definitions DEFINITIONS)
      endif()
    endif()
  endforeach()

  if(ARG_INCLUDES)
    su2_merge_lists(includes ARG_INCLUDES)
  endif()

  set(dependencies "${SU2_DEPENDENCIES};${NAME}")
  set(SU2_DEPENDENCIES ${dependencies} CACHE INTERNAL "")
  set(SU2_DEPENDENCY_${NAME}_INCLUDES ${includes} CACHE INTERNAL "")
  set(SU2_DEPENDENCY_${NAME}_DEFINITIONS ${definitions} CACHE INTERNAL "")
  set(SU2_DEPENDENCY_${NAME}_TARGETS ${ARG_TARGETS} CACHE INTERNAL "")
  set(SU2_DEPENDENCY_${NAME}_EXPLICIT ${ARG_EXPLICIT} CACHE INTERNAL "")

endfunction()

#[[
su2_get_dependency_options(<NAME> [LINK_LIBRARIES <VAR>]
                            [DEFINITIONS <VAR>]
                            [INCLUDES <VAR>]
                            [EXPLICIT <VAR>])

Use this function to get options to a dependency that was added by su2_add_dependency(<NAME>)
LINK_LIBRARIES: sets <VAR> to SU2_DEPENDENCY_<NAME>_TARGETS
DEFINITIONS: sets <VAR> to SU2_DEPENDENCY_<NAME>_DEFINITIONS
INCLUDES: sets <VAR> to SU2_DEPENDENCY_<NAME>_INCLUDES
EXPLICIT: sets <VAR> to SU2_DEPENDENCY_<NAME>_EXPLICIT
]]
function(su2_get_dependency_options NAME)
  set(options)
  set(oneValueArgs LINK_LIBRARIES DEFINITIONS INCLUDES EXPLICIT)
  set(multiValueArgs)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT ${NAME} IN_LIST SU2_DEPENDENCIES)
    message(FATAL_ERROR "Invalid dependency specified: ${NAME}")
  endif()

  if(ARG_INCLUDES)
    set(${ARG_INCLUDES} ${SU2_DEPENDENCY_${NAME}_INCLUDES} PARENT_SCOPE)
  endif()
  if(ARG_DEFINITIONS)
    set(${ARG_DEFINITIONS} ${SU2_DEPENDENCY_${NAME}_DEFINITIONS} PARENT_SCOPE)
  endif()
  if(ARG_LINK_LIBRARIES)
    set(${ARG_LINK_LIBRARIES} ${SU2_DEPENDENCY_${NAME}_TARGETS} PARENT_SCOPE)
  endif()
  if(ARG_EXPLICIT)
    set(${ARG_EXPLICIT} ${SU2_DEPENDENCY_${NAME}_EXPLICIT} PARENT_SCOPE)
  endif()
endfunction()

#[[
su2_get_target_options(<NAME> [LINK_LIBRARIES <VAR>]
                        [DEFINITIONS <VAR>]
                        [INCLUDES <VAR>])

su2_get_dependency_options variant for valid CMake targets
LINK_LIBRARIES: sets <VAR> INTERFACE_LINK_LIBRARIES of target <NAME>
DEFINITIONS: sets <VAR> to INTERFACE_COMPILE_DEFINITIONS of target <NAME>
INCLUDES: sets <VAR> to INTERFACE_INCLUDE_DIRECTORIES of target <NAME>
]]
function(su2_get_target_options NAME)
  set(options)
  set(oneValueArgs LINK_LIBRARIES DEFINITIONS INCLUDES)
  set(multiValueArgs)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT TARGET ${NAME})
    message(FATAL_ERROR "${NAME} is not a valid CMake target")
  endif()

  if(ARG_INCLUDES)
    get_target_property(OUT ${NAME} INTERFACE_INCLUDE_DIRECTORIES)
    if(OUT)
      set(${ARG_INCLUDES} ${OUT} PARENT_SCOPE)
    endif()
  endif()
  if(ARG_DEFINITIONS)
    get_target_property(OUT ${NAME} INTERFACE_COMPILE_DEFINITIONS)
    if(OUT)
      set(${ARG_DEFINITIONS} ${OUT} PARENT_SCOPE)
    endif()
  endif()
  if(ARG_LINK_LIBRARIES)
    get_target_property(OUT ${NAME} INTERFACE_LINK_LIBRARIES)
    if(OUT)
      set(${ARG_LINK_LIBRARIES} ${OUT} PARENT_SCOPE)
    endif()
  endif()
endfunction()

#[[
su2_get_compile_options([LINK_LIBRARIES <VAR>] [DEFINITIONS <VAR>]
                        [INCLUDES <VAR>] [LIBS ...])

Get all compile options that should be used for the build target
LIBS: list of additional CMake targets to link against or WITH_<NAME> dependencies added with su2_add_dependency
LINK_LIBRARIES: sets <VAR> to libraries to be linked from implicit dependencies and LIBS
DEFINITIONS: sets <VAR> to compile definitions of implicit dependencies and LIBS
INCLUDES: sets <VAR> to include directories of implicit dependencies and LIBS
]]
function(su2_get_compile_options)
  set(options)
  set(oneValueArgs LINK_LIBRARIES DEFINITIONS INCLUDES)
  set(multiValueArgs LIBS)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(includes)
  set(definitions)
  set(links)

  foreach(VAR IN LISTS ARG_LIBS SU2_DEPENDENCIES)
    if(TARGET ${VAR})
      su2_get_target_options(${VAR}
        INCLUDES I
        DEFINITIONS D
        LINK_LIBRARIES L)
      list(APPEND L ${VAR})
    elseif(${VAR} MATCHES "WITH_(.*)")
      set(lib ${CMAKE_MATCH_1})
      if(NOT ${lib} IN_LIST SU2_DEPENDENCIES)
        message(FATAL_ERROR "Unknown dependency ${lib}")
      endif()
      su2_get_dependency_options(${lib}
        INCLUDES I
        DEFINITIONS D
        LINK_LIBRARIES L)
    elseif(${VAR} IN_LIST SU2_DEPENDENCIES)
      if(NOT SU2_DEPENDENCY_${VAR}_EXPLICIT)
        su2_get_dependency_options(${VAR}
          INCLUDES I
          DEFINITIONS D
          LINK_LIBRARIES L)
      else()
        set(I)
        set(D)
        set(L)
      endif()
    else()
      message(FATAL_ERROR "Invalid option ${VAR}")
    endif()

    su2_merge_lists(includes I)
    su2_merge_lists(definitions D)
    su2_merge_lists(links L)
  endforeach()

  if(ARG_INCLUDES)
    set(${ARG_INCLUDES} ${includes} PARENT_SCOPE)
  endif()
  if(ARG_DEFINITIONS)
    set(${ARG_DEFINITIONS} ${definitions} PARENT_SCOPE)
  endif()
  if(ARG_LINK_LIBRARIES)
    set(${ARG_LINK_LIBRARIES} ${links} PARENT_SCOPE)
  endif()
endfunction()

#[[
su2_setup_target(<TARGET> [NO_PIC]
                  [NO_INSTALL] [INSTALL_DIR <dir>]
                  [LIBS ...])

convenience function to setup provided target to compile with currenly enabled external libraries and
setup additional common compilation options/install targets

TARGET: valid CMake target to setup
NO_PIC: disable compilation with -fPIC (only valid for static library targets)
NO_INSTALL: disable installing <TARGET> to <dir> (only valid for targets that are not static libraries)
INSTALL_DIR: <TARGET> install directory, default is bin
LIBS: additional libraries to link against, either valid CMake targets or one of SU2 dependencies in the form of WITH_<NAME>
]]
function(su2_setup_target target)
  set(options NO_PIC NO_INSTALL)
  set(oneValueArgs INSTALL_DIR)
  set(multiValueArgs LIBS)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # set default install dir
  if(NOT ARG_INSTALL_DIR)
    set(ARG_INSTALL_DIR bin)
  endif()

  get_target_property(target_type ${target} TYPE)
  if(target_type STREQUAL "STATIC_LIBRARY")
    # compile static libraries with -fPIC in not disabled
    if(NOT ARG_NO_PIC)
      set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif()
  else()
    # setup install target if not disabled
    if(NOT ARG_NO_INSTALL)
      install(FILES "$<TARGET_FILE:${target}>" DESTINATION ${ARG_INSTALL_DIR})
    endif()
  endif()

  su2_get_compile_options(
    INCLUDES I
    DEFINITIONS D
    LINK_LIBRARIES L
    LIBS ${ARG_LIBS})
  target_compile_definitions(${target} PUBLIC ${D})
  target_include_directories(${target} PUBLIC ${I})
  target_link_libraries(${target} PUBLIC ${L})

  # debug
  su2_debug(
    "${target} setup with:
    includes: ${I}
    definitions: ${D}
    link libraries: ${L}")
endfunction()

#[[
su2_download_module(<NAME> <ALT_NAME> <GIT_REPO> <COMMIT_SHA>)

Download and install git repositories
NAME: name of the repository to use
ALT_NAME: directory in SU2 root to install repository at
GIT_REPO: URL of the git repository
COMMIT_SHA: SHA of the commit to download
]]
function(su2_download_module
  NAME ALT_NAME GIT_REPO COMMIT_SHA)
  message(
    STATUS
      "Initializing ${NAME} '${COMMIT_SHA}'
=====================================================================
Downloading module from ${GIT_REPO}")

set(LOCAL_FILE "${CMAKE_BINARY_DIR}/${COMMIT_SHA}.zip")

  file(
    DOWNLOAD "${GIT_REPO}/archive/${COMMIT_SHA}.zip" "${LOCAL_FILE}"
    LOG LOG
    STATUS STATUS)
  list(GET STATUS 0 RESULT)

  if(RESULT)
    list(GET STATUS 1 ERR_MSG)
    file(WRITE ${CMAKE_BINARY_DIR}/configure.log ${LOG})
    file(WRITE ${CMAKE_BINARY_DIR}/configure.err ${ERR_MSG})
    message(
      FATAL_ERROR
        "Download of module ${NAME} failed.
To download it manually, perform the following steps:
  - Download the zip at \"${GIT_REPO}/archive/${COMMIT_SHA}.zip\"
  - Extract the archive to ${CMAKE_SOURCE_DIR}/${ALT_NAME}
  - Execute command 'touch ${CMAKE_SOURCE_DIR}/${ALT_NAME}/${COMMIT_SHA}'
  - Run CMake again")
  endif()

  message(STATUS "Extracting archive...")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar -xzf ${LOCAL_FILE}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OUTPUT_VARIABLE LOG
    ERROR_VARIABLE ERR
    RESULT_VARIABLE RESULT)

  if(RESULT)
    message(FATAL_ERROR "Extraction of module ${NAME} failed: ${RESULT}
-- Command: ${CMAKE_COMMAND} -E tar -xzf ${LOCAL_FILE}
-- Output: ${LOG}
-- Error: ${ERR}")
  endif()

  message(STATUS "Creating identifier...")
  file(RENAME "${CMAKE_BINARY_DIR}/${NAME}-${COMMIT_SHA}" "${CMAKE_SOURCE_DIR}/${ALT_NAME}")
  file(TOUCH "${CMAKE_SOURCE_DIR}/${ALT_NAME}/${COMMIT_SHA}")
  file(REMOVE "${LOCAL_FILE}")

endfunction()

#[[
su2_option(<NAME> [SETUP <VALUE> <HIDE_VALUE> <TYPE> <DOCSTRING> [HIDE])

Use this to create on option that can be hidden from the cmake-gui if some conditions are not met
NAME: name of the option
SETUP: initialize <NAME> option with <VALUE>, <TYPE> and <DOCSTRING> are used the same way as setting cache variables
HIDE: hide <NAME> option
]]
function(su2_option NAME)
  set(options HIDE)
  set(oneValueArgs)
  set(multiValueArgs SETUP)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # set default value
  if(NOT DEFINED SU2_OPTION_${NAME}_HIDDEN)
    set(SU2_OPTION_${NAME}_HIDDEN
        FALSE
        CACHE INTERNAL "")
  endif()

  if(ARG_SETUP)
    list(GET ARG_SETUP 0 VALUE)
    list(GET ARG_SETUP 1 HIDE_VALUE)
    list(GET ARG_SETUP 2 TYPE)
    list(SUBLIST ARG_SETUP 3 -1 DOCSTRING)

    if(SU2_OPTION_${NAME}_HIDDEN)
      # if hidden, set to last cached value and make it visible until HIDE is called again
      set(${NAME}
          ${SU2_OPTION_${NAME}_VALUE}
          CACHE ${TYPE} "${DOCSTRING}" FORCE)
    else()
      # not hidden so use default cache variable setup
      set(${NAME}
          ${VALUE}
          CACHE ${TYPE} "${DOCSTRING}")
    endif()

    # clear hidden flag
    set(SU2_OPTION_${NAME}_HIDDEN
        FALSE
        CACHE INTERNAL "")
    set(SU2_OPTION_${NAME}_HIDE_VALUE
        ${HIDE_VALUE}
        CACHE INTERNAL "")

    su2_debug("Setup ${NAME} with ${VALUE} \"${DOCSTRING}\", value will be ${HIDE_VALUE} when hidden")
  endif()

  if(ARG_HIDE)
    # cache current value and set hidden flag before setting <NAME> to internal variable
    set(SU2_OPTION_${NAME}_VALUE
        ${${NAME}}
        CACHE INTERNAL "")
    set(SU2_OPTION_${NAME}_HIDDEN
        TRUE
        CACHE INTERNAL "")
    set(${NAME}
        ${SU2_OPTION_${NAME}_HIDE_VALUE}
        CACHE INTERNAL "")
  endif()

endfunction()
