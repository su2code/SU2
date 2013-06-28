# SU^2 Makefile.am for MDC
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

bin_PROGRAMS = SU2_MDC

SU2_MDC_SOURCES = ../include/SU2_MDC.hpp \
		../src/SU2_MDC.cpp

SU2_MDC_CXXFLAGS =
SU2_MDC_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
SU2_MDC_CXXFLAGS += @CGNS_CXX@
SU2_MDC_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
SU2_MDC_CXXFLAGS += @METIS_CXX@
SU2_MDC_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
SU2_MDC_CXXFLAGS += @TECIO_CXX@
SU2_MDC_LDADD += @TECIO_LD@
# endif
