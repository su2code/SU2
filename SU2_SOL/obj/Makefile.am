# SU2 Makefile.am for SU2_SOL
# SU2 v3.2.7 "eagle"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2015 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.
#
# Created: 		3/16/13
# Last updated: 	4/09/13
#
# Created by: 	M. Colonno
# Updated by: 	T. Economon
#
######################################################################################

ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = ../bin/SU2_SOL

___bin_SU2_SOL_SOURCES = ../include/SU2_SOL.hpp \
		../src/SU2_SOL.cpp \
		../../SU2_CFD/include/output_structure.hpp \
		../../SU2_CFD/src/output_structure.cpp \
		../../SU2_CFD/src/output_cgns.cpp \
		../../SU2_CFD/src/output_tecplot.cpp \
    ../../SU2_CFD/src/output_fieldview.cpp \
    ../../SU2_CFD/src/output_su2.cpp \
    ../../SU2_CFD/src/output_paraview.cpp \
		../../SU2_CFD/include/solver_structure.hpp \
    ../../SU2_CFD/src/solver_structure.cpp \
    ../../SU2_CFD/include/variable_structure.hpp \
    ../../SU2_CFD/src/variable_structure.cpp

___bin_SU2_SOL_CXXFLAGS =
___bin_SU2_SOL_LDADD = ../../Common/lib/libSU2.a

# always link to built dependencies from ./externals
___bin_SU2_SOL_CXXFLAGS += @su2_externals_INCLUDES@
___bin_SU2_SOL_LDADD += @su2_externals_LIBS@

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
