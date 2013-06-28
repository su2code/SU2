#########
# tkogl #
#########

include_directories(${TCL_INCLUDE_PATH} ${TK_INCLUDE_PATH})

set(tkogl_FILES 
	tkogl.c
	get.c
	gencyl.c
	load3ds.c
	nurbs.c
	quadric.c
	tess.c
	printstr.c
	feedback.c
	tkoglparse.c)

add_library(tkogl STATIC ${tkogl_FILES})

