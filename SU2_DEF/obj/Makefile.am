################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_DEF
# \author M. Colonno, T. Economon, F. Palacios
# \version 4.2.0 "Cardinal"
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
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
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
################################################################################

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = ../bin/SU2_DEF

___bin_SU2_DEF_SOURCES = \
  ../include/SU2_DEF.hpp \
  ../src/SU2_DEF.cpp 
  

___bin_SU2_DEF_CXXFLAGS =
___bin_SU2_DEF_LDADD = ../../Common/lib/libSU2.a \
  ../../SU2_CFD/src/___bin_SU2_CFD-output_structure.o \
  ../../SU2_CFD/src/___bin_SU2_CFD-output_cgns.o \
  ../../SU2_CFD/src/___bin_SU2_CFD-output_tecplot.o \
  ../../SU2_CFD/src/___bin_SU2_CFD-output_fieldview.o \
  ../../SU2_CFD/src/___bin_SU2_CFD-output_su2.o \
  ../../SU2_CFD/src/___bin_SU2_CFD-output_paraview.o \
  ../../SU2_CFD/src/___bin_SU2_CFD-solver_structure.o \
  ../../SU2_CFD/src/___bin_SU2_CFD-variable_structure.o

# always link to built dependencies from ./externals
___bin_SU2_DEF_CXXFLAGS += @su2_externals_INCLUDES@
___bin_SU2_DEF_LDADD += @su2_externals_LIBS@

# if BUILD_CGNS
___bin_SU2_DEF_CXXFLAGS += @CGNS_CXX@
___bin_SU2_DEF_LDADD += @CGNS_LD@
# endif

# if BUILD_HDF5
___bin_SU2_DEF_CXXFLAGS += @HDF5_CXX@
___bin_SU2_DEF_LDADD += @HDF5_LD@
# endif

# if BUILD_SZIP
___bin_SU2_DEF_CXXFLAGS += @SZIP_CXX@
___bin_SU2_DEF_LDADD += @SZIP_LD@
# endif

# if BUILD_ZLIB
___bin_SU2_DEF_CXXFLAGS += @ZLIB_CXX@
___bin_SU2_DEF_LDADD += @ZLIB_LD@
# endif

# if BUILD_MUTATIONPP
___bin_SU2_DEF_CXXFLAGS += @MUTATIONPP_CXX@
___bin_SU2_DEF_LDADD += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
___bin_SU2_DEF_CXXFLAGS += @JSONCPP_CXX@
___bin_SU2_DEF_LDADD += @JSONCPP_LD@
# endif
