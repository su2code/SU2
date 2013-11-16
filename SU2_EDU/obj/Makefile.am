# SU^2 Makefile.am for EDU
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated: 	11/09/13
#
# Created by: 		M. Colonno
# Updated by: 		F. Palacios
#
##########################

# AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = ../bin/SU2_EDU

___bin_SU2_EDU_SOURCES = ../include/definition_structure.hpp \
		../include/integration_structure.hpp \
		../include/integration_structure.inl \
		../include/iteration_structure.hpp \
		../include/numerics_structure.hpp \
		../include/numerics_structure.inl \
		../include/output_structure.hpp \
		../include/solver_structure.hpp \
		../include/solver_structure.inl \
		../include/SU2_EDU.hpp \
		../include/variable_structure.hpp \
		../include/variable_structure.inl \
		../src/definition_structure.cpp \
		../src/integration_structure.cpp \
		../src/integration_time.cpp \
		../src/iteration_structure.cpp \
		../src/numerics_direct_mean.cpp \
		../src/numerics_direct_turbulent.cpp \
		../src/numerics_structure.cpp \
		../src/output_cgns.cpp \
		../src/output_structure.cpp \
		../src/output_tecplot.cpp \
		../src/output_paraview.cpp \
		../src/solver_direct_mean.cpp \
		../src/solver_direct_turbulent.cpp \
		../src/solver_structure.cpp \
		../src/SU2_EDU.cpp \
		../src/variable_direct_mean.cpp \
		../src/variable_direct_turbulent.cpp \
		../src/variable_structure.cpp

___bin_SU2_EDU_CXXFLAGS =
___bin_SU2_EDU_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
___bin_SU2_EDU_CXXFLAGS += @CGNS_CXX@
___bin_SU2_EDU_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
___bin_SU2_EDU_CXXFLAGS += @METIS_CXX@
___bin_SU2_EDU_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
___bin_SU2_EDU_CXXFLAGS += @TECIO_CXX@
___bin_SU2_EDU_LDADD += @TECIO_LD@
# endif

# if BUILD_MUTATIONPP
___bin_SU2_EDU_CXXFLAGS += @MUTATIONPP_CXX@
___bin_SU2_EDU_LDADD += @MUTATIONPP_LD@
# endif