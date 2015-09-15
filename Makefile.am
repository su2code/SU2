################################################################################
#
# \file Makefile.am
# \brief Global makefile for the SU2 project
# \author M. Colonno, T. Economon, F. Palacios
# \version 4.0.1 "Cardinal"
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
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

# Build third-party optional dependencies first
SUBDIRS = externals

# Build the SU2 C++ common library
SUBDIRS += Common/lib

# Install the SU2 Python scripts to the proper location
if BUILD_NORMAL
SUBDIRS += SU2_PY
endif

# Build each of the SU2 C++ modules
if BUILD_CFD
SUBDIRS +=SU2_CFD/obj
endif

if BUILD_FSI
SUBDIRS +=SU2_FSI/obj
endif

if BUILD_DOT
SUBDIRS +=SU2_DOT/obj
endif

if BUILD_MSH
SUBDIRS +=SU2_MSH/obj
endif

if BUILD_DEF
SUBDIRS +=SU2_DEF/obj
endif

if BUILD_SOL
SUBDIRS +=SU2_SOL/obj
endif

if BUILD_GEO
SUBDIRS +=SU2_GEO/obj
endif
