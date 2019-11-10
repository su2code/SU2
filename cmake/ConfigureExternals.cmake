set(SU2_DEPENDENCIES
    ""
    CACHE INTERNAL "")

include(${CMAKE_CURRENT_LIST_DIR}/functions.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/TLS.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ConfigureThreads.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ConfigureMPI.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ConfigureCGNS.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ConfigureMetis.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ConfigureParmetis.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ConfigureTecio.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ConfigureCodi.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ConfigureMutationpp.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ConfigureMKL.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ConfigureCLI11.cmake)

su2_debug("SU2 dependencies: ${SU2_DEPENDENCIES}")
foreach(dep IN LISTS SU2_DEPENDENCIES)
  su2_debug(
    "  ${dep}:
        TARGETS: ${SU2_DEPENDENCY_${dep}_TARGETS}
        INCLUDES: ${SU2_DEPENDENCY_${dep}_INCLUDES}
        DEFINITIONS: ${SU2_DEPENDENCY_${dep}_DEFINITIONS}
        EXPLICIT: ${SU2_DEPENDENCY_${dep}_EXPLICIT}")
endforeach()
