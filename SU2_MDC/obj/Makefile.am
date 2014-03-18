# SU^2 Makefile.am for MDC
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

bin_PROGRAMS = ../bin/SU2_MDC

___bin_SU2_MDC_SOURCES = ../include/SU2_MDC.hpp \
		../src/SU2_MDC.cpp

___bin_SU2_MDC_CXXFLAGS =
___bin_SU2_MDC_LDADD = ../../Common/lib/libSU2.a

# always link to built dependencies from ./contrib
___bin_SU2_MDC_CXXFLAGS += @su2_contrib_INCLUDES@
___bin_SU2_MDC_LDADD += @su2_contrib_LIBS@

# if BUILD_CGNS
___bin_SU2_MDC_CXXFLAGS += @CGNS_CXX@
___bin_SU2_MDC_LDADD += @CGNS_LD@
# endif

# if BUILD_MUTATIONPP
___bin_SU2_MDC_CXXFLAGS += @MUTATIONPP_CXX@
___bin_SU2_MDC_LDADD += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
___bin_SU2_MDC_CXXFLAGS += @JSONCPP_CXX@
___bin_SU2_MDC_LDADD += @JSONCPP_LD@
# endif
