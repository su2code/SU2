# SU^2 Makefile.am for MAC
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

bin_PROGRAMS = SU2_MAC

SU2_MAC_SOURCES = ../include/SU2_MAC.hpp \
		../src/SU2_MAC.cpp

SU2_MAC_CXXFLAGS =
SU2_MAC_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
SU2_MAC_CXXFLAGS += @CGNS_CXX@
SU2_MAC_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
SU2_MAC_CXXFLAGS += @METIS_CXX@
SU2_MAC_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
SU2_MAC_CXXFLAGS += @TECIO_CXX@
SU2_MAC_LDADD += @TECIO_LD@
# endif
