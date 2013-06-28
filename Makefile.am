# SU^2 global Makefile.am
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated: 	3/16/13
#
# Created by: 		M. Colonno
# Updated by: 		M. Colonno
# Updated by:       T. Lukaczyk 4/4/2013
#
##########################

# AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

SUBDIRS = Common/lib 

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

if BUILD_SMC
SUBDIRS +=SU2_SMC/obj
endif

if BUILD_SOL
SUBDIRS +=SU2_SOL/obj
endif

if BUILD_GDC
SUBDIRS +=SU2_GDC/obj
endif
