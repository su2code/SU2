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

check_PROGRAMS  =
check_PROGRAMS += unit_test_config
check_PROGRAMS += resolution_integration_test
check_PROGRAMS += resolution_3d_integration_test
check_PROGRAMS += resolution_gradient_test

# Eventually these will likely be required
AM_CPPFLAGS = -I$(top_srcdir)/Common/include
#AM_LDFLAGS = 
#LIBS = 
LDADD = $(top_builddir)/Common/lib/libSU2.a

# always link to built dependencies from ./externals
AM_CPPFLAGS += @su2_externals_INCLUDES@
LDADD += @su2_externals_LIBS@

# Test sources
unit_test_config_SOURCES = unit_test_config.cpp
resolution_integration_test_SOURCES = resolution_integration_test.cpp
resolution_3d_integration_test_SOURCES = resolution_3d_integration_test.cpp
resolution_gradient_test_SOURCES = resolution_gradient_test.cpp

# Tests to be run
TESTS  =
TESTS += unit_test_config
TESTS += resolution_integration_test
TESTS += resolution_3d_integration_test
TESTS += resolution_gradient_test

# Any of above we expect to fail?
#XFAIL_TESTS =
