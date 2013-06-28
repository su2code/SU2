# SU^2 Makefile.am for DDC
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

bin_PROGRAMS = SU2_DDC

SU2_DDC_SOURCES = ../include/SU2_DDC.hpp \
		../src/SU2_DDC.cpp

SU2_DDC_CXXFLAGS =
SU2_DDC_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
SU2_DDC_CXXFLAGS += @CGNS_CXX@
SU2_DDC_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
SU2_DDC_CXXFLAGS += @METIS_CXX@
SU2_DDC_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
SU2_DDC_CXXFLAGS += @TECIO_CXX@
SU2_DDC_LDADD += @TECIO_LD@
# endif
