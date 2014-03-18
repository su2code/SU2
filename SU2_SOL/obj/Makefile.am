# SU^2 Makefile.am for SOL
# Stanford University Aerospace Design Lab
#
# Created: 		3/16/13
# Last updated: 	4/09/13
#
# Created by: 	M. Colonno
# Updated by: 	T. Economon
#
##########################

ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = ../bin/SU2_SOL

___bin_SU2_SOL_SOURCES = ../include/SU2_SOL.hpp \
		../src/SU2_SOL.cpp \
		../../SU2_CFD/include/output_structure.hpp \
		../../SU2_CFD/src/output_structure.cpp \
		../../SU2_CFD/src/output_cgns.cpp \
		../../SU2_CFD/src/output_tecplot.cpp \
		../../SU2_CFD/src/output_paraview.cpp \
		../../SU2_CFD/include/solver_structure.hpp \
                ../../SU2_CFD/src/solver_structure.cpp \
                ../../SU2_CFD/include/variable_structure.hpp \
                ../../SU2_CFD/src/variable_structure.cpp

___bin_SU2_SOL_CXXFLAGS =
___bin_SU2_SOL_LDADD = ../../Common/lib/libSU2.a

# always link to built dependencies from ./contrib
___bin_SU2_SOL_CXXFLAGS += @su2_contrib_INCLUDES@
___bin_SU2_SOL_LDADD += @su2_contrib_LIBS@

# if BUILD_CGNS
___bin_SU2_SOL_CXXFLAGS += @CGNS_CXX@
___bin_SU2_SOL_LDADD += @CGNS_LD@
# endif

# if BUILD_MUTATIONPP
___bin_SU2_SOL_CXXFLAGS += @MUTATIONPP_CXX@
___bin_SU2_SOL_LDADD += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
___bin_SU2_SOL_CXXFLAGS += @JSONCPP_CXX@
___bin_SU2_SOL_LDADD += @JSONCPP_LD@
# endif
