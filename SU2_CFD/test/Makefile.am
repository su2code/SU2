################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_CFD tests
# \author T. A. Oliver
# \version 5.0.0 "Raven"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
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

# Includes
AM_CPPFLAGS =  -I$(top_srcdir)/Common/include
AM_CPPFLAGS += -I$(top_srcdir)/SU2_CFD/include
AM_CPPFLAGS += @su2_externals_INCLUDES@

LDADD =  $(top_builddir)/SU2_CFD/obj/libSU2Core.a
LDADD += $(top_builddir)/Common/lib/libSU2.a
LDADD += @su2_externals_LIBS@

# Programs to be built for testing
check_PROGRAMS =
check_PROGRAMS += viscous_2d
check_PROGRAMS += viscous_residual
check_PROGRAMS += viscous_general_vs_ideal
check_PROGRAMS += viscous_incomp
check_PROGRAMS += viscous_general_2d

# Test sources
viscous_2d_SOURCES = viscous_2d_test.cpp
viscous_residual_SOURCES = viscous_residual_test.cpp
viscous_general_vs_ideal_SOURCES = viscous_general_vs_ideal_test.cpp
viscous_incomp_SOURCES = viscous_incomp_test.cpp
viscous_general_2d_SOURCES = viscous_general_2d_test.cpp

# Tests to be run
TESTS =
TESTS += viscous_2d
TESTS += viscous_residual
TESTS += viscous_general_vs_ideal
TESTS += viscous_incomp
TESTS += viscous_general_2d
