################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_SOL
# \author M. Colonno, T. Economon, F. Palacios
# \version 7.3.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

bin_PROGRAMS = ../bin/SU2_SOL

___bin_SU2_SOL_SOURCES = \
  ../include/SU2_SOL.hpp \
  ../src/SU2_SOL.cpp

___bin_SU2_SOL_CXXFLAGS =
___bin_SU2_SOL_LDADD = \
  ../../SU2_CFD/obj/libSU2Core.a \
  ../../Common/lib/libSU2.a

# always link to built dependencies from ./externals
___bin_SU2_SOL_CXXFLAGS += @su2_externals_INCLUDES@
___bin_SU2_SOL_LDADD += @su2_externals_LIBS@
___bin_SU2_SOL_LDADD += @su2_externals_LIBPTHREAD@

# if BUILD_MUTATIONPP
___bin_SU2_SOL_CXXFLAGS += @MUTATIONPP_CXX@
___bin_SU2_SOL_LDADD += @MUTATIONPP_LD@
# endif
