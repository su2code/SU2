# SU^2 Makefile.am for SOL
# Stanford University Aerospace Design Lab
#
# Created: 		3/16/13
# Last updated: 	3/16/13
#
# Created by: 	M. Colonno
# Updated by: 	M. Colonno
#
##########################

ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = ../bin/SU2_SOL

___bin_SU2_SOL_SOURCES = ../include/output_structure.hpp \
		../include/solution_structure.hpp \
		../include/SU2_SOL.hpp \
		../include/variable_structure.hpp \
		../src/output_CleanUp.cpp \
		../src/output_structure.cpp \
		../src/output_WriteCGNS.cpp \
		../src/output_WriteTecplotASCII.cpp \
		../src/output_WriteTecplotBinary.cpp \
		../src/solution_structure.cpp \
		../src/SU2_SOL.cpp \
		../src/variable_structure.cpp
		
___bin_SU2_SOL_CXXFLAGS =
___bin_SU2_SOL_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
___bin_SU2_SOL_CXXFLAGS += @CGNS_CXX@
___bin_SU2_SOL_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
___bin_SU2_SOL_CXXFLAGS += @METIS_CXX@
___bin_SU2_SOL_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
___bin_SU2_SOL_CXXFLAGS += @TECIO_CXX@
___bin_SU2_SOL_LDADD += @TECIO_LD@
# endif
