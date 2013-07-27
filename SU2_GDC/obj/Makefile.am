# SU^2 Makefile.am for GDC
# Stanford University Aerospace Design Lab
#
# Created: 		3/16/13
# Last updated: 	4/09/13
#
# Created by: 	M. Colonno
# Updated by: 	T. Economon
#
##########################

ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = ../bin/SU2_GDC

___bin_SU2_GDC_SOURCES = ../include/SU2_GDC.hpp \
		../src/SU2_GDC.cpp \
		../../SU2_CFD/include/output_structure.hpp \
		../../SU2_CFD/src/output_structure.cpp \
		../../SU2_CFD/src/output_cgns.cpp \
		../../SU2_CFD/src/output_tecplot.cpp
		
___bin_SU2_GDC_CXXFLAGS =
___bin_SU2_GDC_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
___bin_SU2_GDC_CXXFLAGS += @CGNS_CXX@
___bin_SU2_GDC_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
___bin_SU2_GDC_CXXFLAGS += @METIS_CXX@
___bin_SU2_GDC_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
___bin_SU2_GDC_CXXFLAGS += @TECIO_CXX@
___bin_SU2_GDC_LDADD += @TECIO_LD@
# endif
