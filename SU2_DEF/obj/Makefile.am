# SU2 Makefile.am for MDC
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

bin_PROGRAMS = ../bin/SU2_DEF

___bin_SU2_DEF_SOURCES = ../include/SU2_DEF.hpp \
		../src/SU2_DEF.cpp

___bin_SU2_DEF_CXXFLAGS =
___bin_SU2_DEF_LDADD = ../../Common/lib/libSU2.a

# always link to built dependencies from ./externals
___bin_SU2_DEF_CXXFLAGS += @su2_externals_INCLUDES@
___bin_SU2_DEF_LDADD += @su2_externals_LIBS@

# if BUILD_CGNS
___bin_SU2_DEF_CXXFLAGS += @CGNS_CXX@
___bin_SU2_DEF_LDADD += @CGNS_LD@
# endif

# if BUILD_MUTATIONPP
___bin_SU2_DEF_CXXFLAGS += @MUTATIONPP_CXX@
___bin_SU2_DEF_LDADD += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
___bin_SU2_DEF_CXXFLAGS += @JSONCPP_CXX@
___bin_SU2_DEF_LDADD += @JSONCPP_LD@
# endif
