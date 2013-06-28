#############
# cgnstools #
#############

# get X11 in case it is not in a standard location
if (NOT WIN32)
  include(FindX11)
  # Xmu doesn't appear to be included in FindX11
  find_library(X11_Xmu_LIB Xmu)
  if (X11_Xmu_LIB)
    mark_as_advanced(FORCE X11_Xmu_LIB)
  endif (X11_Xmu_LIB)
endif (NOT WIN32)

include(FindTCL)
include(FindOpenGL)

set(Missing "")
if (TCL_FOUND)
  mark_as_advanced(FORCE TCL_INCLUDE_PATH TCL_LIBRARY)
else (TCL_FOUND)
  mark_as_advanced(CLEAR TCL_INCLUDE_PATH TCL_LIBRARY)
  list(APPEND Missing "TCL")
endif (TCL_FOUND)

if (TK_FOUND)
  mark_as_advanced(FORCE TK_INCLUDE_PATH TK_LIBRARY)
else (TK_FOUND)
  mark_as_advanced(CLEAR TK_INCLUDE_PATH TK_LIBRARY)
  list(APPEND Missing "TK")
endif (TK_FOUND)

if (OPENGL_FOUND)
  mark_as_advanced(FORCE OPENGL_INCLUDE_PATH OPENGL_gl_LIBRARY)
else (OPENGL_FOUND)
  mark_as_advanced(CLEAR OPENGL_INCLUDE_PATH OPENGL_gl_LIBRARY)
  list(APPEND Missing "OPENGL_gl")
endif (OPENGL_FOUND)

if (OPENGL_GLU_FOUND)
  mark_as_advanced(FORCE OPENGL_glu_LIBRARY)
else (OPENGL_GLU_FOUND)
  mark_as_advanced(CLEAR OPENGL_glu_LIBRARY)
  list(APPEND Missing "OPENGL_glu")
endif (OPENGL_GLU_FOUND)

if (Missing)
  message(SEND_ERROR "The path and library needs to be defined for:
    ${Missing}")
else (Missing)
  add_subdirectory(cgnsview)
  add_subdirectory(calclib)
  add_subdirectory(cgnscalc)
  add_subdirectory(cgnsplot)
  add_subdirectory(common)
  add_subdirectory(tkogl)
  add_subdirectory(utilities)
endif (Missing)

if (WIN32)
  file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX} WIN_INSTALL_DIR)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cgconfig.bat
"set CG_BIN_DIR=${WIN_INSTALL_DIR}\\bin
set CG_LIB_DIR=${WIN_INSTALL_DIR}\\share
")
# don't need this since dll is put in the bin directory
#  if (CGNS_USE_SHARED)
#    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cgconfig.bat
#      "set PATH=${WIN_INSTALL_DIR}\\lib;%PATH%"
#    )
#  endif (CGNS_USE_SHARED)
# add Tcl path for dll's
  if (TK_WISH)
    get_filename_component(TK_DLL ${TK_WISH} PATH)
    file(TO_NATIVE_PATH ${TK_DLL} TK_DLL_DIR)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/cgconfig.bat
      "set PATH=${TK_DLL_DIR};%PATH%"
    )
  endif (TK_WISH)

  install(PROGRAMS
	${CMAKE_CURRENT_BINARY_DIR}/cgconfig.bat
	DESTINATION bin)
else (WIN32)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cgconfig
"CG_BIN_DIR=${CMAKE_INSTALL_PREFIX}/bin; export CG_BIN_DIR
CG_LIB_DIR=${CMAKE_INSTALL_PREFIX}/share/cgnstools; export CG_LIB_DIR
")

  install(PROGRAMS
	${CMAKE_CURRENT_BINARY_DIR}/cgconfig
	DESTINATION bin)
endif (WIN32)
