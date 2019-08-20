################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_MET
# \author M. Colonno, T. Economon, F. Palacios
# \version 6.1.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.$(OBJEXT)rg>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
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
# License along with SU2. If not, see <http://www.gnu.$(OBJEXT)rg/licenses/>.
#
################################################################################

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = 

if BUILD_REVERSE
bin_PROGRAMS += ../bin/SU2_MET

su2_met_sources = ../include/SU2_MET.hpp \
    ../src/SU2_MET.cpp

su2_met_ldadd =
su2_met_cxx_flags =

su2_met_ldadd += \
  ../../SU2_CFD/obj/libSU2Core_AD.a \
  ../../Common/lib/libSU2_AD.a
endif

# always link to built dependencies from ./externals
su2_met_cxx_flags += @su2_externals_INCLUDES@
su2_met_ldadd += @su2_externals_LIBS@
su2_met_ldadd += @su2_externals_LIBPTHREAD@

# if BUILD_MUTATIONPP
su2_met_cxx_flags += @MUTATIONPP_CXX@
su2_met_ldadd += @MUTATIONPP_LD@
# endif

___bin_SU2_MET_SOURCES = ${su2_met_sources}
___bin_SU2_MET_CXXFLAGS = @REVERSE_CXX@ ${su2_met_cxx_flags} 
___bin_SU2_MET_LDADD = @REVERSE_LIBS@ ${su2_met_ldadd}

endif
