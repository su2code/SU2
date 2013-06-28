# SU^2 global Makefile.am
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated: 	11/15/12
#
# Created by: 		M. Colonno
# Updated by: 		M. Colonno
#
##########################

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

SUBDIRS = trunk/Common/lib trunk/SU2Py

if BUILD_CFD
SUBDIRS +=trunk/SU2_CFD/bin
endif

if BUILD_DDC
SUBDIRS +=trunk/SU2_DDC/bin
endif

if BUILD_GPC
SUBDIRS +=trunk/SU2_GPC/bin
endif

if BUILD_MAC
SUBDIRS +=trunk/SU2_MAC/bin
endif

if BUILD_MDC
SUBDIRS +=trunk/SU2_MDC/bin
endif

if BUILD_PBC
SUBDIRS +=trunk/SU2_PBC/bin
endif

if BUILD_SMC
SUBDIRS +=trunk/SU2_SMC/bin
endif

