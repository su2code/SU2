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

libSU2Core_cxx_flags = -fPIC
libSU2Core_libadd = 

# always link to built dependencies from ./externals
libSU2Core_cxx_flags += @su2_externals_INCLUDES@

# if BUILD_LAPACK
libSU2Core_cxx_flags += @LAPACK_CXX@
libSU2Core_libadd += @LAPACK_LD@
# endif

libSU2Core_a_CXXFLAGS = $(libSU2Core_cxx_flags)
libSU2Core_a_LIBADD = $(libSU2Core_libadd)

check_PROGRAMS  =
check_PROGRAMS += aniso_viscous_flux_inc
check_PROGRAMS += aniso_viscous_flux
check_PROGRAMS += aniso_viscous_jacobian
check_PROGRAMS += aniso_viscous_jacobian_compare
check_PROGRAMS += aniso_viscous_residual_test
check_PROGRAMS += load_hybrid_constants_test
check_PROGRAMS += eigensolver_test

# Eventually these will likely be required
AM_CPPFLAGS =
AM_CPPFLAGS += -I$(top_srcdir)/Common/include 
AM_CPPFLAGS += -I$(top_srcdir)/SU2_CFD/include
AM_CPPFLAGS += @su2_externals_INCLUDES@ @LAPACK_CXX@

#AM_LDFLAGS = 
#LIBS = 

LDADD  =
LDADD += $(top_builddir)/Common/lib/libSU2.a
LDADD += $(top_builddir)/SU2_CFD/obj/libSU2Core.a
LDADD += @su2_externals_LIBS@ @LAPACK_LD@

LDFLAGS  = -mkl

# Test sources
aniso_viscous_flux_inc_SOURCES         = aniso_viscous_flux_inc.cpp
aniso_viscous_flux_SOURCES             = aniso_viscous_flux.cpp
aniso_viscous_jacobian_SOURCES         = aniso_viscous_jacobian.cpp
aniso_viscous_jacobian_compare_SOURCES = aniso_viscous_jacobian_compare.cpp
aniso_viscous_residual_test_SOURCES    = aniso_viscous_residual_test.cpp
load_hybrid_constants_test_SOURCES     = load_hybrid_constants_test.cpp
eigensolver_test_SOURCES               = eigensolver_test.cpp

# Tests to be run
TESTS  =
TESTS += aniso_viscous_flux_inc
TESTS += aniso_viscous_flux
TESTS += aniso_viscous_jacobian
TESTS += aniso_viscous_jacobian_compare
TESTS += aniso_viscous_residual_test
TESTS += load_hybrid_constants_test
TESTS += eigensolver_test

# Any of above we expect to fail?
#XFAIL_TESTS =
