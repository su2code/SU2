# SU2 global Makefile.am
# SU2 v3.2.7 "eagle"
#
# SU2 Lead Developers: Dr. Francisco Palacios (francisco.palacios@boeing.com).
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
# Created: 		10/26/12
# Last updated: 	03/24/14
#
# Created by: 		M. Colonno
# Updated by: 		T. Economon
#
# Major Changes/Contributions:
# METIS & TECIO integration added by Ben Kirk, 3/24/2014
#
######################################################################################

# AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

# build third-party optional dependencies first
SUBDIRS = externals

SUBDIRS += Common/lib

if BUILD_NORMAL
SUBDIRS += SU2_PY
endif

if BUILD_CFD
SUBDIRS +=SU2_CFD/obj
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
