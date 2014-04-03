# SU^2 autoconf file
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated:         03/24/14
#
# Created by: 		M. Colonno  10/26/12
# Updated by: 		T. Economon 	3/24/14
#
# Major Changes/Contributions:
# METIS & TECIO integration added by Ben Kirk, 3/24/2014
#
##########################

AC_INIT([SU^2], [3.0.0], [susquared-dev@lists.stanford.edu],
        [SU^2], [http://su2.stanford.edu])
AC_PREREQ([2.59])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([subdir-objects])  # [SU^2], [2.0]
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES(yes)])


# AC_CONFIG_HEADERS([config.h])
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_RANLIB
AC_PROG_MKDIR_P

# --------------------------------------------------------------
# Check for important type sizes
# --------------------------------------------------------------
AC_CHECK_SIZEOF(short int)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(long int)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(void *)



# --------------------------------------------------------------
# Build up contributed packages that we will build from source
# --------------------------------------------------------------
su2_externals_INCLUDES=""
su2_externals_LIBS=""

# Tecplot
CONFIGURE_TECIO
if (test $enabletecio = yes); then
  su2_externals_INCLUDES="$TECIO_INCLUDE $su2_externals_INCLUDES"
  su2_externals_LIBS="$TECIO_LIB $su2_externals_LIBS"
else
  su2_externals_INCLUDES="-DNO_TECIO $su2_externals_INCLUDES"
fi
AM_CONDITIONAL(BUILD_TECIO, test x$enabletecio = xyes)
AC_CONFIG_FILES([externals/tecio/Makefile])

# Metis
CONFIGURE_METIS
if (test $enablemetis = yes); then
  su2_externals_INCLUDES="$METIS_INCLUDE $su2_externals_INCLUDES"
  su2_externals_LIBS="$METIS_LIB $su2_externals_LIBS"
else
  su2_externals_INCLUDES="-DNO_METIS $su2_externals_INCLUDES"
fi
AM_CONDITIONAL(BUILD_METIS, test x$enablemetis = xyes)
AC_CONFIG_FILES([externals/metis/Makefile])


AC_SUBST([su2_externals_INCLUDES])
AC_SUBST([su2_externals_LIBS])


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
AC_ARG_WITH(Mutationpp-lib,
   	AS_HELP_STRING([--with-Mutationpp-lib[=ARG]], [Mutation++ library directory, ARG = path to libmutation++.dylib]),
	[with_Mutationpp_lib=$withval], [with_Mutationpp_lib="no"])
AC_ARG_WITH(Mutationpp-include,
   	AS_HELP_STRING([--with-Mutationpp-include[=ARG]], [Mutation++ include directory, ARG = path to mutation++.h]),
	[with_Mutationpp_include=$withval], [with_Mutationpp_include="no"])
AC_ARG_WITH(Jsoncpp-lib,
   	AS_HELP_STRING([--with-Jsoncpp-lib[=ARG]], [Jsoncpp library directory, ARG = path to jsoncpp_SYS.a]),
	[with_Jsoncpp_lib=$withval], [with_Jsoncpp_lib="no"])
AC_ARG_WITH(Jsoncpp-include,
   	AS_HELP_STRING([--with-Jsoncpp-include[=ARG]], [Jsoncpp include directory, ARG = path to json.h]),
	[with_Jsoncpp_include=$withval], [with_Jsoncpp_include="no"])
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
AC_ARG_ENABLE(SOL,
   	AS_HELP_STRING([--disable-SOL], [build the SU2_SOL executable (default = yes)]),
	[build_SOL=$enableval], [build_SOL="yes"])
AC_ARG_ENABLE(GDC,
   	AS_HELP_STRING([--disable-GDC], [build the SU2_GDC executable (default = yes)]),
	[build_GDC=$enableval], [build_GDC="yes"])
AC_ARG_ENABLE(UMC,
   	AS_HELP_STRING([--disable-UMC], [build the SU2_UMC executable (default = yes)]),
	[build_UMC=$enableval], [build_UMC="yes"])

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

##########################

# check for lib
Mutationpplib=libmutation++.dylib
Mutationppheader=mutation++.h
have_Mutationpp="no"

if test "$with_Mutationpp_lib" != "no"
then
  	AC_CHECK_FILE([$with_Mutationpp_lib/$Mutationpplib],[have_Mutationpp="yes"],[have_Mutationpp="no"])
	if test "$have_Mutationpp" == "no"
	then
		AC_MSG_ERROR([Mutation++ requested but library file not found.])
	fi
fi

# check for header
if test "$have_Mutationpp" != "no"
then
	if test "$with_Mutationpp_include" != "no"
	then
		AC_CHECK_FILE([$with_Mutationpp_include/$Mutationppheader],[have_Mutationpp='yes'],[have_Mutationpp='no'])
		if test "$have_Mutationpp" == "no"
		then
			AC_MSG_ERROR([Mutation++ requested but header file not found.])
		fi
	else
		have_Mutationpp="no"
	fi
fi

if test "$have_Mutationpp" != "no"
then
  	AM_CONDITIONAL(BUILD_MUTATIONPP,true)
	MUTATIONPP_CXX=-I"$with_Mutationpp_include"
	MUTATIONPP_LD="$with_Mutationpp_lib"/$Mutationpplib
else
	AM_CONDITIONAL(BUILD_MUTATIONPP,false)
	MUTATIONPP_CXX=-DNO_MUTATIONPP
	MUTATION_LD=
fi
AC_SUBST([MUTATIONPP_CXX])
AC_SUBST([MUTATIONPP_LD])

##########################

# check for lib
Jsoncpplib=libjson.a
Jsoncppheader=json/json.h
have_Jsoncpp="no"

if test "$with_Jsoncpp_lib" != "no"
then
  	AC_CHECK_FILE([$with_Jsoncpp_lib/$Jsoncpplib],[have_Jsoncpp="yes"],[have_Jsoncpp="no"])
	if test "$have_Jsoncpp" == "no"
	then
		AC_MSG_ERROR([Jsoncpp requested but library file not found.])
	fi
fi

# check for header
if test "$have_Jsoncpp" != "no"
then
	if test "$with_Jsoncpp_include" != "no"
	then
		AC_CHECK_FILE([$with_Jsoncpp_include/$Jsoncppheader],[have_Jsoncpp='yes'],[have_Jsoncpp='no'])
		if test "$have_Jsoncpp" == "no"
		then
			AC_MSG_ERROR([Jsoncpp requested but header file not found.])
		fi
	else
		have_Jsoncpp="no"
	fi
fi

if test "$have_Jsoncpp" != "no"
then
  	AM_CONDITIONAL(BUILD_JSONCPP,true)
	JSONCPP_CXX=-I"$with_Jsoncpp_include"
	JSONCPP_LD="$with_Jsoncpp_lib"/$Jsoncpplib
else
	AM_CONDITIONAL(BUILD_JSONCPP,false)
	JSONCPP_CXX=-DNO_JSONCPP
	JSONCPP_LD=
fi
AC_SUBST([JSONCPP_CXX])
AC_SUBST([JSONCPP_LD])


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
AM_CONDITIONAL([BUILD_SOL],[test $build_SOL != "no"])
AM_CONDITIONAL([BUILD_GDC],[test $build_GDC != "no"])
AM_CONDITIONAL([BUILD_UMC],[test $build_UMC != "no"])

###########################

AC_OUTPUT(Makefile \
        externals/Makefile \
	Common/lib/Makefile \
	SU2_CFD/obj/Makefile \
	SU2_DDC/obj/Makefile \
	SU2_GPC/obj/Makefile \
	SU2_MAC/obj/Makefile \
	SU2_MDC/obj/Makefile \
	SU2_PBC/obj/Makefile \
	SU2_SOL/obj/Makefile \
	SU2_GDC/obj/Makefile \
	SU2_UMC/obj/Makefile \
	SU2_PY/Makefile)

abs_srcdir=`(cd $srcdir && pwd)`
abs_prefix=$prefix #`(cd $prefix && pwd)`


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
    Metis support:        $enablemetis
    TecIO support:        $enabletecio
    Mutation++ support:   $have_Mutationpp
    Jsoncpp support: 	  $have_Jsoncpp

    External includes:    $su2_externals_INCLUDES
    External libs:        $su2_externals_LIBS

    Build SU2_CFD:        $build_CFD
    Build SU2_DDC:        $build_DDC
    Build SU2_GPC:        $build_GPC
    Build SU2_MAC:        $build_MAC
    Build SU2_MDC:        $build_MDC
    Build SU2_PBC:        $build_PBC
    Build SU2_SOL:        $build_SOL
    Build SU2_GDC:        $build_GDC
    Build SU2_UMC:        $build_UMC

Please be sure to add the [\$]SU2_HOME and [\$]SU2_RUN environment variables,
and update your [\$]PATH (and [\$]PYTHONPATH if applicable) with [\$]SU2_RUN.

Based on the input to this configuration, add these lines to your .bashrc file:

export SU2_RUN="${abs_prefix}/bin"
export SU2_HOME="${abs_srcdir}"
export PATH=[\$]PATH:[\$]SU2_RUN
export PYTHONPATH=[\$]PYTHONPATH:[\$]SU2_RUN

])
