# SU^2 Makefile.am for GPC
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated: 	3/12/13
#
# Created by: 		M. Colonno
# Updated by: 		M. Colonno
#
##########################

# AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = ../bin/SU2_GPC

___bin_SU2_GPC_SOURCES = ../include/SU2_GPC.hpp \
		../src/SU2_GPC.cpp

___bin_SU2_GPC_CXXFLAGS =
___bin_SU2_GPC_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
___bin_SU2_GPC_CXXFLAGS += @CGNS_CXX@
___bin_SU2_GPC_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
___bin_SU2_GPC_CXXFLAGS += @METIS_CXX@
___bin_SU2_GPC_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
___bin_SU2_GPC_CXXFLAGS += @TECIO_CXX@
___bin_SU2_GPC_LDADD += @TECIO_LD@
# endif
