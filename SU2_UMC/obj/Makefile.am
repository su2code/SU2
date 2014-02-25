# SU^2 Makefile.am for UMC
# Stanford University Aerospace Design Lab
#
# Created: 		3/16/13
# Last updated: 	10/30/13
#
# Created by: 	T. Economon
# Updated by: 	T. Economon
#
##########################

ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = ../bin/SU2_UMC

___bin_SU2_UMC_SOURCES = ../include/SU2_UMC.hpp \
		../src/SU2_UMC.cpp
		
___bin_SU2_UMC_CXXFLAGS =
___bin_SU2_UMC_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
___bin_SU2_UMC_CXXFLAGS += @CGNS_CXX@
___bin_SU2_UMC_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
___bin_SU2_UMC_CXXFLAGS += @METIS_CXX@
___bin_SU2_UMC_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
___bin_SU2_UMC_CXXFLAGS += @TECIO_CXX@
___bin_SU2_UMC_LDADD += @TECIO_LD@
# endif

# if BUILD_MUTATIONPP
___bin_SU2_UMC_CXXFLAGS += @MUTATIONPP_CXX@
___bin_SU2_UMC_LDADD += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
___bin_SU2_UMC_CXXFLAGS += @JSONCPP_CXX@
___bin_SU2_UMC_LDADD += @JSONCPP_LD@
# endif