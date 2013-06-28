###########
# calclib #
###########

set(calclib_FILES 
	calc.c
	vec.c
	vecsym.c)

add_library(calclib STATIC ${calclib_FILES})

# Add test
