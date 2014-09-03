# SU2 autoconf file
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

AC_INIT([SU2], [3.2.1], [susquared-dev@lists.stanford.edu],
        [SU2], [http://su2.stanford.edu])
AC_PREREQ([2.59])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([subdir-objects foreign])
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
  su2_externals_INCLUDES="-DHAVE_TECIO $TECIO_INCLUDE $su2_externals_INCLUDES"
  su2_externals_LIBS="$TECIO_LIB $su2_externals_LIBS"
else
  su2_externals_INCLUDES="$su2_externals_INCLUDES"
fi
AM_CONDITIONAL(BUILD_TECIO, test x$enabletecio = xyes)
AC_CONFIG_FILES([externals/tecio/Makefile])

# Metis
CONFIGURE_METIS
if (test $enablemetis = yes); then
  su2_externals_INCLUDES="-DHAVE_METIS $METIS_INCLUDE $su2_externals_INCLUDES"
  su2_externals_LIBS="$METIS_LIB $su2_externals_LIBS"
else
  su2_externals_INCLUDES="$su2_externals_INCLUDES"
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
AC_ARG_ENABLE(PRT,
   	AS_HELP_STRING([--disable-PRT], [build the SU2_PRT executable (default = yes)]),
	[build_PRT=$enableval], [build_PRT="yes"])
AC_ARG_ENABLE(DOT,
   	AS_HELP_STRING([--disable-DOT], [build the SU2_DOT executable (default = yes)]),
	[build_DOT=$enableval], [build_DOT="yes"])
AC_ARG_ENABLE(MSH,
   	AS_HELP_STRING([--disable-MSH], [build the SU2_MSH executable (default = yes)]),
	[build_MSH=$enableval], [build_MSH="yes"])
AC_ARG_ENABLE(DEF,
   	AS_HELP_STRING([--disable-DEF], [build the SU2_DEF executable (default = yes)]),
	[build_DEF=$enableval], [build_DEF="yes"])
AC_ARG_ENABLE(SOL,
   	AS_HELP_STRING([--disable-SOL], [build the SU2_SOL executable (default = yes)]),
	[build_SOL=$enableval], [build_SOL="yes"])
AC_ARG_ENABLE(GEO,
   	AS_HELP_STRING([--disable-GEO], [build the SU2_GEO executable (default = yes)]),
	[build_GEO=$enableval], [build_GEO="yes"])

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
	CGNS_CXX="-DHAVE_CGNS -I"$with_CGNS_include""
	CGNS_LD="$with_CGNS_lib"/$CGNSlib
else
	AM_CONDITIONAL(BUILD_CGNS,false)
	CGNS_CXX=
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
	MUTATIONPP_CXX="-DHAVE_MUTATIONPP -I"$with_Mutationpp_include""
	MUTATIONPP_LD="$with_Mutationpp_lib"/$Mutationpplib
else
	AM_CONDITIONAL(BUILD_MUTATIONPP,false)
	MUTATIONPP_CXX=
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
	JSONCPP_CXX="-DHAVE_JSONCPP -I"$with_Jsoncpp_include""
	JSONCPP_LD="$with_Jsoncpp_lib"/$Jsoncpplib
else
	AM_CONDITIONAL(BUILD_JSONCPP,false)
	JSONCPP_CXX=
	JSONCPP_LD=
fi
AC_SUBST([JSONCPP_CXX])
AC_SUBST([JSONCPP_LD])


########################### MPI

have_MPI="no"
if test "$with_MPI" != "no"
then
	CXX="$with_MPI"
	CXXFLAGS="-DHAVE_MPI "$CXXFLAGS""
	have_MPI="yes"
fi

###########################

AM_CONDITIONAL([BUILD_CFD],[test $build_CFD != "no"])
AM_CONDITIONAL([BUILD_PRT],[test $build_PRT != "no"])
AM_CONDITIONAL([BUILD_DOT],[test $build_DOT != "no"])
AM_CONDITIONAL([BUILD_MSH],[test $build_MSH != "no"])
AM_CONDITIONAL([BUILD_DEF],[test $build_DEF != "no"])
AM_CONDITIONAL([BUILD_SOL],[test $build_SOL != "no"])
AM_CONDITIONAL([BUILD_GEO],[test $build_GEO != "no"])

###########################

AC_OUTPUT(Makefile \
        externals/Makefile \
	Common/lib/Makefile \
	SU2_CFD/obj/Makefile \
	SU2_PRT/obj/Makefile \
	SU2_DOT/obj/Makefile \
	SU2_MSH/obj/Makefile \
	SU2_DEF/obj/Makefile \
	SU2_SOL/obj/Makefile \
	SU2_GEO/obj/Makefile \
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
    Build SU2_PRT:        $build_PRT
    Build SU2_DOT:        $build_DOT
    Build SU2_MSH:        $build_MSH
    Build SU2_DEF:        $build_DEF
    Build SU2_SOL:        $build_SOL
    Build SU2_GEO:        $build_GEO

Please be sure to add the [\$]SU2_HOME and [\$]SU2_RUN environment variables,
and update your [\$]PATH (and [\$]PYTHONPATH if applicable) with [\$]SU2_RUN.

Based on the input to this configuration, add these lines to your .bashrc file:

export SU2_RUN="${abs_prefix}/bin"
export SU2_HOME="${abs_srcdir}"
export PATH=[\$]PATH:[\$]SU2_RUN
export PYTHONPATH=[\$]PYTHONPATH:[\$]SU2_RUN

])
