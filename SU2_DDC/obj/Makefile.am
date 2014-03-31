# SU^2 Makefile.am for DDC
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

bin_PROGRAMS = ../bin/SU2_DDC

___bin_SU2_DDC_SOURCES = ../include/SU2_DDC.hpp \
		../src/SU2_DDC.cpp

___bin_SU2_DDC_CXXFLAGS =
___bin_SU2_DDC_LDADD = ../../Common/lib/libSU2.a

# always link to built dependencies from ./externals
___bin_SU2_DDC_CXXFLAGS += @su2_externals_INCLUDES@
___bin_SU2_DDC_LDADD += @su2_externals_LIBS@

# if BUILD_CGNS
___bin_SU2_DDC_CXXFLAGS += @CGNS_CXX@
___bin_SU2_DDC_LDADD += @CGNS_LD@
# endif

# if BUILD_MUTATIONPP
___bin_SU2_DDC_CXXFLAGS += @MUTATIONPP_CXX@
___bin_SU2_DDC_LDADD += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
___bin_SU2_DDC_CXXFLAGS += @JSONCPP_CXX@
___bin_SU2_DDC_LDADD += @JSONCPP_LD@
# endif
