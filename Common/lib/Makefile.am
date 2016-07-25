################################################################################
#
# \file Makefile.am
# \brief Makefile for the SU2 common library
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

noinst_LIBRARIES =

if BUILD_NORMAL
noinst_LIBRARIES += libSU2.a
endif

if BUILD_DIRECTDIFF
noinst_LIBRARIES += libSU2_DIRECTDIFF.a
endif

if BUILD_REVERSE
noinst_LIBRARIES += libSU2_AD.a
endif

lib_sources = \
  ../include/datatypes/primitive_structure.hpp \
  ../include/datatypes/primitive_structure.inl \
  ../include/datatypes/complex_structure.hpp \
  ../include/datatypes/complex_structure.inl \
  ../include/ad_structure.hpp \
  ../include/ad_structure.inl \
  ../include/config_structure.hpp \
  ../include/config_structure.hpp \
  ../include/config_structure.inl \
  ../include/dual_grid_structure.hpp \
  ../include/dual_grid_structure.inl \
  ../include/geometry_structure.hpp \
  ../include/geometry_structure.inl \
  ../include/grid_adaptation_structure.hpp \
  ../include/grid_adaptation_structure.inl \
  ../include/grid_movement_structure.hpp \
  ../include/grid_movement_structure.inl \
  ../include/linear_solvers_structure.hpp \
  ../include/linear_solvers_structure.inl \
  ../include/linear_solvers_structure_b.hpp \
  ../include/option_structure.hpp \
  ../include/primal_grid_structure.hpp \
  ../include/primal_grid_structure.inl \
  ../include/vector_structure.hpp \
  ../include/vector_structure.inl \
  ../include/matrix_structure.hpp \
  ../include/matrix_structure.inl \
  ../include/mpi_structure.hpp \
  ../include/mpi_structure.inl \
  ../include/datatype_structure.hpp \
  ../include/datatype_structure.inl \
  ../include/interpolation_structure.hpp \
  ../include/gauss_structure.hpp \
  ../include/gauss_structure.inl \
  ../include/element_structure.hpp \
  ../include/element_structure.inl \
  ../include/adt_structure.hpp \
  ../include/adt_structure.inl \
  ../src/config_structure.cpp \
  ../src/dual_grid_structure.cpp \
  ../src/geometry_structure.cpp \
  ../src/grid_adaptation_structure.cpp \
  ../src/grid_movement_structure.cpp \
  ../src/linear_solvers_structure.cpp \
  ../src/primal_grid_structure.cpp \
  ../src/vector_structure.cpp \
  ../src/matrix_structure.cpp \
  ../src/mpi_structure.cpp \
  ../src/linear_solvers_structure_b.cpp \
  ../src/ad_structure.cpp \
  ../src/gauss_structure.cpp \
  ../src/element_structure.cpp \
  ../src/element_linear.cpp \
  ../src/interpolation_structure.cpp \
  ../src/adt_structure.cpp

lib_cxxflags =
lib_ldadd =

# always link to built dependencies from ./externals
lib_cxxflags += @su2_externals_INCLUDES@
lib_ldadd += @su2_externals_LIBS@

# if BUILD_CGNS
lib_cxxflags += @CGNS_CXX@
lib_ldadd += @CGNS_LD@
# endif

# if BUILD_HDF5
lib_cxxflags += @HDF5_CXX@
lib_ldadd += @HDF5_LD@
# endif

# if BUILD_SZIP
lib_cxxflags += @SZIP_CXX@
lib_ldadd += @SZIP_LD@
# endif

# if BUILD_ZLIB
lib_cxxflags += @ZLIB_CXX@
lib_ldadd += @ZLIB_LD@
# endif

# if BUILD_MUTATIONPP
lib_cxxflags += @MUTATIONPP_CXX@
lib_ldadd += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
lib_cxxflags += @JSONCPP_CXX@
lib_ldadd += @JSONCPP_LD@
# endif

if BUILD_NORMAL
libSU2_a_SOURCES = $(lib_sources)
libSU2_a_CXXFLAGS = ${lib_cxxflags}
libSU2_a_LIBADD = ${lib_ldadd}
endif

if BUILD_DIRECTDIFF
libSU2_DIRECTDIFF_a_SOURCES = $(lib_sources)
libSU2_DIRECTDIFF_a_CXXFLAGS = @DIRECTDIFF_CXX@ ${lib_cxxflags}
libSU2_DIRECTDIFF_a_LIBADD =  ${lib_ldadd}
endif

if BUILD_REVERSE
libSU2_AD_a_SOURCES = $(lib_sources)
libSU2_AD_a_CXXFLAGS = @REVERSE_CXX@ ${lib_cxxflags}
libSU2_AD_a_LIBADD =  ${lib_ldadd}
endif

