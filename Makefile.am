# SU^2 global Makefile.am
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

if BUILD_DDC
SUBDIRS +=SU2_DDC/obj
endif

if BUILD_GPC
SUBDIRS +=SU2_GPC/obj
endif

if BUILD_MAC
SUBDIRS +=SU2_MAC/obj
endif

if BUILD_MDC
SUBDIRS +=SU2_MDC/obj
endif

if BUILD_PBC
SUBDIRS +=SU2_PBC/obj
endif

if BUILD_SOL
SUBDIRS +=SU2_SOL/obj
endif

if BUILD_GDC
SUBDIRS +=SU2_GDC/obj
endif

if BUILD_UMC
SUBDIRS +=SU2_UMC/obj
endif
