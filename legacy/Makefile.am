################################################################################
#
# \file Makefile.am
# \brief Global makefile for the SU2 project
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
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

# build third-party optional dependencies first
SUBDIRS = externals

SUBDIRS += Common/lib

if BUILD_PY
SUBDIRS += SU2_PY
endif

if BUILD_CFD
SUBDIRS +=SU2_CFD/obj
endif

if BUILD_DOT
SUBDIRS +=SU2_DOT/obj
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

if BUILD_PY_WRAPPER
SUBDIRS += SU2_PY/pySU2
endif
