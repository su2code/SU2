# SU^2 Makefile.am for SMC
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

bin_PROGRAMS = SU2_SMC

SU2_SMC_SOURCES = ../include/SU2_SMC.hpp \
		../src/SU2_SMC.cpp

SU2_SMC_CXXFLAGS =
SU2_SMC_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
SU2_SMC_CXXFLAGS += @CGNS_CXX@
SU2_SMC_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
SU2_SMC_CXXFLAGS += @METIS_CXX@
SU2_SMC_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
SU2_SMC_CXXFLAGS += @TECIO_CXX@
SU2_SMC_LDADD += @TECIO_LD@
# endif

