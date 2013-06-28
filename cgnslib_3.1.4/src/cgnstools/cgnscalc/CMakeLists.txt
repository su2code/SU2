############
# cgnscalc #
############

include_directories(${TCL_INCLUDE_PATH} ${TK_INCLUDE_PATH})
include_directories(../calclib)
link_directories(.)
if (CGNS_USE_SHARED)
  link_libraries(cgns_shared)
else (CGNS_USE_SHARED)
  link_libraries(cgns_static)
endif (CGNS_USE_SHARED)
link_libraries(calclib ${TCL_LIBRARY} ${TK_LIBRARY})

if (ENABLE_HDF5 AND HDF5_LIBRARY)
  link_libraries(${HDF5_LIBRARY})
  if(HDF5_NEED_ZLIB AND ZLIB_LIBRARY)
    link_libraries(${ZLIB_LIBRARY})
  endif(HDF5_NEED_ZLIB AND ZLIB_LIBRARY)
  if(HDF5_NEED_SZIP AND SZIP_LIBRARY)
    link_libraries(${SZIP_LIBRARY})
  endif(HDF5_NEED_SZIP AND SZIP_LIBRARY)
  if(HDF5_NEED_MPI AND MPI_LIBS)
    link_libraries(${MPI_LIBS})
  endif(HDF5_NEED_MPI AND MPI_LIBS)
endif (ENABLE_HDF5 AND HDF5_LIBRARY)

if (WIN32)
  set(cgnscalc_FILES
	calctcl.c
	winmain.c)
else (WIN32)
  include_directories(${X11_X11_INCLUDE_PATH})
  link_libraries(${X11_X11_LIB} m)
  set(cgnscalc_FILES
	calctcl.c
	calcwish.c)
endif (WIN32)

add_executable(calcwish ${cgnscalc_FILES})
add_dependencies(calcwish calclib cgns_static)

if (WIN32)
  set_target_properties(calcwish PROPERTIES
  	LINK_FLAGS /subsystem:windows)
endif (WIN32)

if (WIN32)
  install(TARGETS
	calcwish
	RUNTIME DESTINATION bin)
  install(PROGRAMS
	cgnscalc.bat
	unitconv.bat
	DESTINATION bin)
  install(FILES
	cgnscalc.tcl
	unitconv.tcl
	cgnscalc.ico
	cgnscalc.png
	cgnscalc-icon.xbm
	cgnscalc-mask.xbm
	unitconv.ico
	unitconv.png
	unitconv.xbm
	DESTINATION share)
else (WIN32)
  install(TARGETS
	calcwish
	RUNTIME DESTINATION bin/cgnstools)
  install(PROGRAMS
	cgnscalc.sh
	DESTINATION bin
	RENAME cgnscalc)
  install(PROGRAMS
	unitconv.sh
	DESTINATION bin
	RENAME unitconv)
  install(FILES
	cgnscalc.tcl
	unitconv.tcl
	cgnscalc.ico
	cgnscalc.png
	cgnscalc-icon.xbm
	cgnscalc-mask.xbm
	unitconv.ico
	unitconv.png
	unitconv.xbm
	DESTINATION share/cgnstools)
endif (WIN32)

