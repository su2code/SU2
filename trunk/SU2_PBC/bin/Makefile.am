# SU^2 Makefile.am for PBC
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

bin_PROGRAMS = SU2_PBC

SU2_PBC_SOURCES = ../include/SU2_PBC.hpp \
		../src/SU2_PBC.cpp

SU2_PBC_CXXFLAGS =
SU2_PBC_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
SU2_PBC_CXXFLAGS += @CGNS_CXX@
SU2_PBC_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
SU2_PBC_CXXFLAGS += @METIS_CXX@
SU2_PBC_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
SU2_PBC_CXXFLAGS += @TECIO_CXX@
SU2_PBC_LDADD += @TECIO_LD@
# endif
