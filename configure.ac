# SU^2 autoconf file
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated: 4/4/2013
#
# Created by: 		M. Colonno  10/26/12
# Updated by: 		M. Colonno  3/16/13
#                   T. Lukaczyk 4/4/13
#
##########################

AC_INIT([SU^2], [2.0.6], [susquared-dev@lists.stanford.edu],
        [SU^2], [http://su2.stanford.edu])
AC_PREREQ([2.59])
AM_INIT_AUTOMAKE  # [SU^2], [2.0]

# optional bindings
AC_ARG_WITH(MPI,
    	AS_HELP_STRING([--with-MPI[=ARG]], [Build parallel version with MPI tools, ARG = MPI C++ compiler]), 
	[with_MPI=$withval], [with_MPI="no"])
AC_ARG_WITH(CGNS-lib,
    	AS_HELP_STRING([--with-CGNS-lib[=ARG]], [CGNS library directory, ARG = path to libcgns.a]), 
	[with_CGNS_lib=$withval], [with_CGNS_lib="no"])
AC_ARG_WITH(CGNS-include,
    	AS_HELP_STRING([--with-CGNS-include[=ARG]], [CGNS include directory, ARG = path to cgnslib.h]), 
	[with_CGNS_include=$withval], [with_CGNS_include="no"])
AC_ARG_WITH(Tecio-lib,
    	AS_HELP_STRING([--with-Tecio-lib[=ARG]], [Tecio library directory, ARG = path to tecio.a or tecio64.a]), 
	[with_Tecio_lib=$withval], [with_Tecio_lib="no"])
AC_ARG_WITH(Tecio-include,
    	AS_HELP_STRING([--with-Tecio-include[=ARG]], [Tecio include directory, ARG = path to TECIO.h]), 
	[with_Tecio_include=$withval], [with_Tecio_include="no"])
AC_ARG_WITH(Metis-lib,
   	AS_HELP_STRING([--with-Metis-lib[=ARG]], [Metis library directory, ARG = path to libmetis.a]), 
	[with_Metis_lib=$withval], [with_Metis_lib="no"])
AC_ARG_WITH(Metis-include,
   	AS_HELP_STRING([--with-Metis-include[=ARG]], [Metis include directory, ARG = path to metis.h]), 
	[with_Metis_include=$withval], [with_Metis_include="no"])
AC_ARG_WITH(Metis-version,
   	AS_HELP_STRING([--with-Metis-version[=ARG]], [Metis major version, ARG = 4 (default) or 5]), 
	[with_Metis_version=$withval], [with_Metis_version="4"])
AC_ARG_ENABLE(CFD,
   	AS_HELP_STRING([--disable-CFD], [build the SU2_CFD executable (default = yes)]), 
	[build_CFD=$enableval], [build_CFD="yes"])
AC_ARG_ENABLE(DDC,
   	AS_HELP_STRING([--disable-DDC], [build the SU2_DDC executable (default = yes)]), 
	[build_DDC=$enableval], [build_DDC="yes"])
AC_ARG_ENABLE(GPC,
   	AS_HELP_STRING([--disable-GPC], [build the SU2_GPC executable (default = yes)]), 
	[build_GPC=$enableval], [build_GPC="yes"])
AC_ARG_ENABLE(MAC,
   	AS_HELP_STRING([--disable-MAC], [build the SU2_MAC executable (default = yes)]), 
	[build_MAC=$enableval], [build_MAC="yes"])
AC_ARG_ENABLE(MDC,
   	AS_HELP_STRING([--disable-MDC], [build the SU2_MDC executable (default = yes)]), 
	[build_MDC=$enableval], [build_MDC="yes"])
AC_ARG_ENABLE(PBC,
   	AS_HELP_STRING([--disable-PBC], [build the SU2_PBC executable (default = yes)]), 
	[build_PBC=$enableval], [build_PBC="yes"])
AC_ARG_ENABLE(SMC,
   	AS_HELP_STRING([--disable-SMC], [build the SU2_SMC executable (default = yes)]), 
	[build_SMC=$enableval], [build_SMC="yes"])
AC_ARG_ENABLE(SOL,
   	AS_HELP_STRING([--disable-SOL], [build the SU2_SOL executable (default = yes)]), 
	[build_SOL=$enableval], [build_SOL="yes"])
AC_ARG_ENABLE(GDC,
   	AS_HELP_STRING([--disable-GDC], [build the SU2_GDC executable (default = yes)]), 
	[build_GDC=$enableval], [build_GDC="yes"])

##########################

# check for lib
CGNSlib=libcgns.a
CGNSheader=cgnslib.h
have_CGNS="no"

if test "$with_CGNS_lib" != "no" 
then
  	AC_CHECK_FILE([$with_CGNS_lib/$CGNSlib],[have_CGNS="yes"],[have_CGNS="no"])   
	if test "$have_CGNS" == "no"
	then
		AC_MSG_ERROR([CGNS requested but library file not found.])
	fi
fi

# check for header
if test "$have_CGNS" != "no"
then
	if test "$with_CGNS_include" != "no"
	then
		AC_CHECK_FILE([$with_CGNS_include/$CGNSheader],[have_CGNS='yes'],[have_CGNS='no']) 
		if test "$have_CGNS" == "no"
		then
			AC_MSG_ERROR([CGNS requested but header file not found.])
		fi
	else
		have_CGNS="no"
	fi  
fi

if test "$have_CGNS" != "no"
then
  	AM_CONDITIONAL(BUILD_CGNS,true)
	CGNS_CXX=-I"$with_CGNS_include"
	CGNS_LD="$with_CGNS_lib"/$CGNSlib
else
	AM_CONDITIONAL(BUILD_CGNS,false)
	CGNS_CXX=-DNO_CGNS
	CGNS_LD=
fi
AC_SUBST([CGNS_CXX])
AC_SUBST([CGNS_LD])

###########################

# check for lib
METISlib=libmetis.a
METISheader=metis.h
METISversion=
have_METIS="no"

if test "$with_Metis_lib" != "no" 
then
  	AC_CHECK_FILE([$with_Metis_lib/$METISlib],[have_METIS="yes"],[have_METIS="no"])   
	if test "$have_METIS" == "no"
	then
		AC_MSG_ERROR([METIS requested but library file not found.])
	fi
fi

# check for header
if test "$have_METIS" != "no"
then
	if test "$with_Metis_include" != "no"
	then
		AC_CHECK_FILE([$with_Metis_include/$METISheader],[have_METIS="yes"],[have_METIS="no"]) 
		if test "$have_METIS" == "no"
		then
			AC_MSG_ERROR([METIS requested but header file not found.])
		fi
	else
		have_METIS="no"
	fi  
fi

# METIS version check
if test "$with_Metis_version" == "5"
then
	METISversion=-DMETIS_5
fi

if test "$have_METIS" != "no"
then
  	AM_CONDITIONAL(BUILD_METIS,true)
	METIS_CXX="-I"$with_Metis_include" "$METISversion""
	METIS_LD="$with_Metis_lib"/$METISlib
else
	AM_CONDITIONAL(BUILD_METIS,false)
	METIS_CXX=-DNO_METIS
	METIS_LD=
	
fi
AC_SUBST([METIS_CXX])
AC_SUBST([METIS_LD])

########################### Tecplot

Teciolib=tecio64.a
Tecioheader=TECIO.h
have_Tecio="no"

# check for lib
if test "$with_Tecio_lib" != "no" 
then
  	AC_CHECK_FILE([$with_Tecio_lib/$Teciolib],[have_Tecio="yes"],[have_Tecio="no"]) 
	if test "$have_Tecio" == "no"
	then
		Teciolib=tecio.a
		AC_CHECK_FILE([$with_Tecio_lib/$Teciolib],[have_Tecio="yes"],[have_Tecio="no"])
	fi
	if test "$have_Tecio" == "no"
	then
		AC_MSG_ERROR([Tecplot linking requested but library file not found.])
	fi
fi

# check for header
if test "$have_Tecio" != "no"
then
	if test "$with_Tecio_include" != "no"
	then
		AC_CHECK_FILE([$with_Tecio_include/$Tecioheader],[have_Tecio='yes'],[have_Tecio='no']) 
		if test "$have_Tecio" == "no"
		then
			AC_MSG_ERROR([Tecplot linking requested but header file not found.])
		fi
	else
		have_Tecio="no"
	fi  
fi

if test "$have_Tecio" != "no"
then
  	AM_CONDITIONAL(BUILD_TECIO,true)
	TECIO_CXX=-I"$with_Tecio_include"
	TECIO_LD="$with_Tecio_lib"/$Teciolib
else
	AM_CONDITIONAL(BUILD_TECIO,false)
	TECIO_CXX=-DNO_TECIO
	TECIO_LD=
fi
AC_SUBST([TECIO_CXX])
AC_SUBST([TECIO_LD])

########################### MPI

have_MPI="no"
if test "$with_MPI" != "no"
then
	CXX="$with_MPI"	
	have_MPI="yes"
else
	CXXFLAGS=""$CXXFLAGS" -DNO_MPI"
fi

###########################

AM_CONDITIONAL([BUILD_CFD],[test $build_CFD != "no"])
AM_CONDITIONAL([BUILD_DDC],[test $build_DDC != "no"])
AM_CONDITIONAL([BUILD_GPC],[test $build_GPC != "no"])
AM_CONDITIONAL([BUILD_MAC],[test $build_MAC != "no"])
AM_CONDITIONAL([BUILD_MDC],[test $build_MDC != "no"])
AM_CONDITIONAL([BUILD_PBC],[test $build_PBC != "no"])
AM_CONDITIONAL([BUILD_SMC],[test $build_SMC != "no"])
AM_CONDITIONAL([BUILD_SOL],[test $build_SOL != "no"])
AM_CONDITIONAL([BUILD_GDC],[test $build_GDC != "no"])

###########################

# AC_CONFIG_HEADERS([config.h])
AC_PROG_CXX
# AC_PROG_CC
# AM_PROG_CC_C_O
AC_PROG_RANLIB
AC_OUTPUT(Makefile \
	Common/lib/Makefile \
	SU2_CFD/obj/Makefile \
	SU2_DDC/obj/Makefile \
	SU2_GPC/obj/Makefile \
	SU2_MAC/obj/Makefile \
	SU2_MDC/obj/Makefile \
	SU2_PBC/obj/Makefile \
	SU2_SMC/obj/Makefile \
	SU2_SOL/obj/Makefile \
	SU2_GDC/obj/Makefile \
	SU2_PY/Makefile)

abs_srcdir=`(cd $srcdir && pwd)`
abs_prefix=`(cd $prefix && pwd)`


AC_MSG_RESULT([

Build Configuration Summary:

    Source code location: ${abs_srcdir}
    Install location: 	  ${abs_prefix}
    Version:              ${VERSION}
    Compiler:             ${CXX}
    Preprocessor flags:   ${CPPFLAGS}
    Compiler flags:       ${CXXFLAGS}
    Linker flags:         ${LDFLAGS}
    CGNS support:         $have_CGNS
    MPI support:          $have_MPI
    Metis support:        $have_METIS

    Build SU2_CFD:        $build_CFD
    Build SU2_DDC:        $build_DDC
    Build SU2_GPC:        $build_GPC
    Build SU2_MAC:        $build_MAC
    Build SU2_MDC:        $build_MDC
    Build SU2_PBC:        $build_PBC
    Build SU2_SMC:        $build_SMC
    Build SU2_SOL:        $build_SOL
    Build SU2_GDC:        $build_GDC

Please be sure to add the [\$]SU2_HOME and [\$]SU2_RUN environment variables, and update your [\$]PATH with [\$]SU2_RUN

Based on the input to this configuration, add these lines to your .bashrc file -

export SU2_RUN="${abs_prefix}/bin"
export SU2_HOME="${abs_srcdir}"
export PATH=[\$]PATH:[\$]SU2_RUN


])
