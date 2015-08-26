################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_DOT
# \author M. Colonno, T. Economon, F. Palacios
# \version 4.0.0 "Cardinal"
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
################################################################################

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = 

if BUILD_NORMAL
bin_PROGRAMS += ../bin/SU2_DOT
endif

if BUILD_REVERSE
bin_PROGRAMS += ../bin/SU2_DOT_REVERSE
endif

su2_dot_sources = ../include/SU2_DOT.hpp \
		../src/SU2_DOT.cpp

su2_dot_cxx_flags = 
su2_dot_ldadd =

# always link to built dependencies from ./externals
su2_dot_cxx_flags += @su2_externals_INCLUDES@
su2_dot_ldadd += @su2_externals_LIBS@

# if BUILD_CGNS
su2_dot_cxx_flags += @CGNS_CXX@
su2_dot_ldadd += @CGNS_LD@
# endif

# if BUILD_MUTATIONPP
su2_dot_cxx_flags += @MUTATIONPP_CXX@
su2_dot_ldadd += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
su2_dot_cxx_flags += @JSONCPP_CXX@
su2_dot_ldadd += @JSONCPP_LD@
# endif

# if BUILD_HDF5
su2_dot_cxx_flags += @HDF5_CXX@
su2_dot_ldadd += @HDF5_LD@
# endif

# if BUILD_SZIP
su2_dot_cxx_flags += @SZIP_CXX@
su2_dot_ldadd += @SZIP_LD@
# endif

# if BUILD_ZLIB
su2_dot_cxx_flags += @ZLIB_CXX@
su2_dot_ldadd += @ZLIB_LD@
# endif

if BUILD_NORMAL
___bin_SU2_DOT_SOURCES = ${su2_dot_sources}
___bin_SU2_DOT_CXXFLAGS = ${su2_dot_cxx_flags}
___bin_SU2_DOT_LDADD = ../../Common/lib/libSU2.a ${su2_dot_ldadd}
endif

if BUILD_REVERSE
___bin_SU2_DOT_REVERSE_SOURCES = ${su2_dot_sources}
___bin_SU2_DOT_REVERSE_CXXFLAGS = @REVERSE_CXX@ ${su2_dot_cxx_flags} 
___bin_SU2_DOT_REVERSE_LDADD = @REVERSE_LIBS@ ../../Common/lib/libSU2_REVERSE.a ${su2_dot_ldadd}
endif
