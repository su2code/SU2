# SU^2 Makefile.am for GDC
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated: 	11/15/12
#
# Created by: 		M. Colonno
# Updated by: 		M. Colonno
#
##########################

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = SU2_GDC

SU2_GDC_SOURCES = ../include/SU2_GDC.hpp \
		../src/SU2_GDC.cpp

SU2_GDC_CXXFLAGS =
SU2_GDC_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
SU2_GDC_CXXFLAGS += @CGNS_CXX@
SU2_GDC_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
SU2_GDC_CXXFLAGS += @METIS_CXX@
SU2_GDC_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
SU2_GDC_CXXFLAGS += @TECIO_CXX@
SU2_GDC_LDADD += @TECIO_LD@
# endif
