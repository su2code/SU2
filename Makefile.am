# SU2 global Makefile.am
# Stanford University Aerospace Design Lab
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
##########################

# AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

# build third-party optional dependencies first
SUBDIRS = externals

SUBDIRS += Common/lib

SUBDIRS += SU2_PY

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