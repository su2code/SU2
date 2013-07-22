# SU^2 Makefile.am for common library
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

libSU2_a_SOURCES = ../src/config_structure.hpp \
		../src/config_structure.inl \
		../src/dual_grid_structure.hpp \
		../src/dual_grid_structure.inl \
		../src/geometry_structure.hpp \
		../src/geometry_structure.inl \
		../src/grid_adaptation_structure.hpp \
		../src/grid_adaptation_structure.inl \
		../src/grid_movement_structure.hpp \
		../src/grid_movement_structure.inl \
		../src/linear_solvers_structure.hpp \
		../src/linear_solvers_structure.inl \
		../src/option_structure.hpp \
		../src/primal_grid_structure.hpp \
		../src/primal_grid_structure.inl \
		../src/vector_structure.hpp \
		../src/vector_structure.inl \
	       ../src/matrix_structure.hpp \
	       ../src/matrix_structure.inl \
		../src/config_structure.cpp \
		../src/dual_grid_structure.cpp \
		../src/geometry_structure.cpp \
		../src/grid_adaptation_structure.cpp \
		../src/grid_movement_structure.cpp \
		../src/linear_solvers_structure.cpp \
		../src/primal_grid_structure.cpp \
	       ../src/vector_structure.cpp \
		../src/matrix_structure.cpp

libSU2_a_CXXFLAGS = 
libSU2_a_LIBADD =

# if BUILD_CGNS
libSU2_a_CXXFLAGS += @CGNS_CXX@ 
libSU2_a_LIBADD += @CGNS_LD@ 
# endif

# if BUILD_METIS
libSU2_a_CXXFLAGS += @METIS_CXX@
libSU2_a_LIBADD += @METIS_LD@
# endif

# if BUILD_TECIO
libSU2_a_CXXFLAGS += @TECIO_CXX@
libSU2_a_LIBADD += @TECIO_LD@
# endif


