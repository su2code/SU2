# SU2 Makefile.am for common library
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated: 	3/12/13
#
# Created by: 		M. Colonno
# Updated by: 		M. Colonno
#
##########################

# AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

noinst_LIBRARIES = libSU2.a

libSU2_a_SOURCES = \
                ../include/config_structure.hpp \
		../include/config_structure.inl \
		../include/dual_grid_structure.hpp \
		../include/dual_grid_structure.inl \
		../include/geometry_structure.hpp \
		../include/geometry_structure.inl \
		../include/grid_adaptation_structure.hpp \
		../include/grid_adaptation_structure.inl \
		../include/grid_movement_structure.hpp \
		../include/grid_movement_structure.inl \
		../include/linear_solvers_structure.hpp \
		../include/linear_solvers_structure.inl \
		../include/option_structure.hpp \
		../include/primal_grid_structure.hpp \
		../include/primal_grid_structure.inl \
		../include/vector_structure.hpp \
		../include/vector_structure.inl \
	        ../include/matrix_structure.hpp \
	        ../include/matrix_structure.inl \
	    ../include/su2mpi.hpp \
		../src/config_structure.cpp \
		../src/dual_grid_structure.cpp \
		../src/geometry_structure.cpp \
		../src/grid_adaptation_structure.cpp \
		../src/grid_movement_structure.cpp \
		../src/linear_solvers_structure.cpp \
		../src/primal_grid_structure.cpp \
	        ../src/vector_structure.cpp \
		../src/matrix_structure.cpp \
		../src/su2mpi.cpp 

libSU2_a_CXXFLAGS =
libSU2_a_LIBADD =

# always link to built dependencies from ./externals
libSU2_a_CXXFLAGS += @su2_externals_INCLUDES@
libSU2_a_LIBADD += @su2_externals_LIBS@

# if BUILD_CGNS
libSU2_a_CXXFLAGS += @CGNS_CXX@
libSU2_a_LIBADD += @CGNS_LD@
# endif

# if BUILD_MUTATIONPP
libSU2_a_CXXFLAGS += @MUTATIONPP_CXX@
libSU2_a_LIBADD += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
libSU2_a_CXXFLAGS += @JSONCPP_CXX@
libSU2_a_LIBADD += @JSONCPP_LD@
# endif
