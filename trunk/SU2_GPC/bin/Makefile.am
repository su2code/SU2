# SU^2 Makefile.am for GPC
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

bin_PROGRAMS = SU2_GPC

SU2_GPC_SOURCES = ../include/SU2_GPC.hpp \
		../src/SU2_GPC.cpp

SU2_GPC_CXXFLAGS =
SU2_GPC_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
SU2_GPC_CXXFLAGS += @CGNS_CXX@
SU2_GPC_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
SU2_GPC_CXXFLAGS += @METIS_CXX@
SU2_GPC_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
SU2_GPC_CXXFLAGS += @TECIO_CXX@
SU2_GPC_LDADD += @TECIO_LD@
# endif
