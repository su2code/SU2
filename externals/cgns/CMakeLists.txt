project(CGNS C)

set(SU2_CGNS_CPPFLAGS "" CACHE STRING "Specific CGNS C Preprocessor flags to use")

set(CGNS_SOURCES
    cgnslib.h
    cgns_header.h
    cgns_io.h
    cgnstypes.h
    fortran_macros.h
    adf/ADF.h
    adf/ADF_fbind.h
    adf/ADF_internals.h
    cgns_error.c
    cgns_internals.c
    cgns_io.c
    cgnslib.c
    adf/ADF_interface.c
    adf/ADF_internals.c)

add_library(cgns STATIC ${CGNS_SOURCES})
target_compile_definitions(cgns PRIVATE ${SU2_CGNS_CPPFLAGS})
set_target_properties(cgns PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(cgns PUBLIC "${CMAKE_CURRENT_LIST_DIR}" "${CMAKE_CURRENT_LIST_DIR}/adf")
add_library(CGNS::CGNS ALIAS cgns)
