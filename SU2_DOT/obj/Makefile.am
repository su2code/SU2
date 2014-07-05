# SU2 Makefile.am for DOT
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

bin_PROGRAMS = ../bin/SU2_DOT

___bin_SU2_DOT_SOURCES = ../include/SU2_DOT.hpp \
		../src/SU2_DOT.cpp

___bin_SU2_DOT_CXXFLAGS =
___bin_SU2_DOT_LDADD = ../../Common/lib/libSU2.a

# always link to built dependencies from ./externals
___bin_SU2_DOT_CXXFLAGS += @su2_externals_INCLUDES@
___bin_SU2_DOT_LDADD += @su2_externals_LIBS@

# if BUILD_CGNS
___bin_SU2_DOT_CXXFLAGS += @CGNS_CXX@
___bin_SU2_DOT_LDADD += @CGNS_LD@
# endif

# if BUILD_MUTATIONPP
___bin_SU2_DOT_CXXFLAGS += @MUTATIONPP_CXX@
___bin_SU2_DOT_LDADD += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
___bin_SU2_DOT_CXXFLAGS += @JSONCPP_CXX@
___bin_SU2_DOT_LDADD += @JSONCPP_LD@
# endif
