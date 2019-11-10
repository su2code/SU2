project(parmetis C)

if(NOT SU2_ENABLE_METIS)
	message(FATAL_ERROR "Parmetis requires Metis to be built")
endif()

set(SU2_PARMETIS_CPPFLAGS "" CACHE STRING "Specific PARMETIS C Preprocessor flags to use")

find_package(MPI REQUIRED)

set(PARMETIS_SOURCES
    include/parmetis.h
	libparmetis/akwayfm.c
	libparmetis/ametis.c
	libparmetis/balancemylink.c
	libparmetis/comm.c
	libparmetis/csrmatch.c
	libparmetis/ctrl.c
	libparmetis/debug.c
	libparmetis/defs.h
	libparmetis/diffutil.c
	libparmetis/frename.c
	libparmetis/gkmetis.c
	libparmetis/gkmpi.c
	libparmetis/graph.c
	libparmetis/initbalance.c
	libparmetis/initmsection.c
	libparmetis/initpart.c
	libparmetis/kmetis.c
	libparmetis/kwayrefine.c
	libparmetis/macros.h
	libparmetis/match.c
	libparmetis/mdiffusion.c
	libparmetis/mesh.c
	libparmetis/mmetis.c
	libparmetis/move.c
	libparmetis/msetup.c
	libparmetis/node_refine.c
	libparmetis/ometis.c
	libparmetis/parmetislib.h
	libparmetis/proto.h
	libparmetis/pspases.c
	libparmetis/redomylink.c
	libparmetis/remap.c
	libparmetis/rename.h
	libparmetis/renumber.c
	libparmetis/rmetis.c
	libparmetis/selectq.c
	libparmetis/serial.c
	libparmetis/stat.c
	libparmetis/struct.h
	libparmetis/timer.c
	libparmetis/util.c
	libparmetis/wave.c
	libparmetis/weird.c
	libparmetis/wspace.c
	libparmetis/xyzpart.c)

add_library(parmetis STATIC ${PARMETIS_SOURCES})
target_compile_definitions(parmetis PRIVATE ${SU2_METIS_CPPFLAGS} ${SU2_PARMETIS_CPPFLAGS})
set_target_properties(parmetis PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(parmetis PUBLIC MPI::MPI_C Metis::Metis)

target_include_directories(
  parmetis
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/include ${CMAKE_CURRENT_LIST_DIR}/libparmetis ${CMAKE_SOURCE_DIR}/externals/metis/GKlib
    ${CMAKE_SOURCE_DIR}/externals/metis/include ${CMAKE_CURRENT_LIST_DIR})

add_library(Parmetis::Parmetis ALIAS parmetis)
