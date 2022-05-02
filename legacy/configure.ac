################################################################################
#
# \file configure.ac
# \brief Main file for configuring the autoconf/automake build process
# \author M. Colonno, T. Economon, F. Palacios, B. Kirk
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
#
# Major Changes/Contributions:
# METIS & TECIO integration added by Ben Kirk, 3/24/2014
# LAPACK connection added by Teus van der Stelt, 11/13/2014
# ParMETIS integration added by Ben Kirk & Thomas D. Economon, 2/3/2015
# CGNS integration added by Thomas D. Economon, 8/17/2016
#
##########################

AC_INIT([SU2], [6.2.0], [su2code-dev@lists.stanford.edu],
        [SU2], [https://su2code.github.io])
AC_PREREQ([2.59])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([foreign])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES(yes)])

# AC_CONFIG_HEADERS([config.h])
SU2_SET_COMPILERS
AM_PROG_CC_C_O
AC_PROG_RANLIB
AC_PROG_MKDIR_P
AC_PROG_CXX

# Python configuration
AM_PATH_PYTHON
AC_ARG_VAR([PYTHON_INCLUDE], [Include flags for python, bypassing python-config])
AC_ARG_VAR([PYTHON_CONFIG], [Path to python-config])
AC_ARG_VAR([PYTHON_LIB], [Python lib flags])
AC_ARG_VAR([PYTHON_EXEC_PREFIX], [Python exec-prefix])
AS_IF([test -z "$PYTHON_INCLUDE"], [
  AS_IF([test -z "$PYTHON_CONFIG"], [
    AC_PATH_PROGS([PYTHON_CONFIG],
                  [python$PYTHON_VERSION-config python-config],
                  [no],
                  [`dirname $PYTHON`])
    AS_IF([test "$PYTHON_CONFIG" = no], [AC_MSG_ERROR([cannot find python-config for $PYTHON.])])
  ])
  AC_MSG_CHECKING([python include flags])
  PYTHON_INCLUDE=`$PYTHON_CONFIG --includes`
  PYTHON_LIB=`$PYTHON_CONFIG --libs --embed`
  PYTHON_EXEC_PREFIX=`$PYTHON_CONFIG --exec-prefix`
  AC_MSG_RESULT([$PYTHON_INCLUDE])
])

#PKG_PROG_PKG_CONFIG

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

# optional bindings for external packages

AC_ARG_WITH(Mutationpp-lib,
    AS_HELP_STRING([--with-Mutationpp-lib[=ARG]], [Mutation++ library directory, ARG = path to libmutation++.a]),
  [with_Mutationpp_lib=$withval], [with_Mutationpp_lib="no"])
AC_ARG_WITH(Mutationpp-include,
    AS_HELP_STRING([--with-Mutationpp-include[=ARG]], [Mutation++ include directory, ARG = path to mutation++.h]),
  [with_Mutationpp_include=$withval], [with_Mutationpp_include="no"])

# Check for disabling of SU2 modules

AC_ARG_ENABLE(CFD,
    AS_HELP_STRING([--disable-CFD], [build the SU2_CFD executable (default = yes)]),
  [build_CFD=$enableval], [build_CFD="yes"])
AC_ARG_ENABLE(DOT,
   	AS_HELP_STRING([--disable-DOT], [build the SU2_DOT executable (default = yes)]),
	[build_DOT=$enableval], [build_DOT="yes"])
AC_ARG_ENABLE(DEF,
    AS_HELP_STRING([--disable-DEF], [build the SU2_DEF executable (default = yes)]),
  [build_DEF=$enableval], [build_DEF="yes"])
AC_ARG_ENABLE(SOL,
    AS_HELP_STRING([--disable-SOL], [build the SU2_SOL executable (default = yes)]),
  [build_SOL=$enableval], [build_SOL="yes"])
AC_ARG_ENABLE(GEO,
    AS_HELP_STRING([--disable-GEO], [build the SU2_GEO executable (default = yes)]),
  [build_GEO=$enableval], [build_GEO="yes"])
AC_ARG_ENABLE(PY_WRAPPER,
    AS_HELP_STRING([--enable-PY_WRAPPER], [wrap the SU2 code with Python (default = no)]),
  [build_PY_WRAPPER="yes"], [build_PY_WRAPPER="no"])

# Check for data types

AC_ARG_ENABLE(normal,
    AS_HELP_STRING([--disable-normal], [build executables with normal datatype (default = yes)]),
    [build_NORMAL=$enableval], [build_NORMAL="yes"])

# Check for the old MPI option so that we can throw an error

AC_ARG_WITH(MPI,
AS_HELP_STRING(),
        [with_MPI=$withval], [with_MPI="no"])
if test "$with_MPI" != "no"
then
    AC_MSG_ERROR([unsupported option --with-MPI. For MPI support, use --enable-mpi and specify your MPI implementations with --with-cc=/path/to/mpicc and --with-cxx=/path/to/mpicxx in your configure call.])
fi

# --------------------------------------------------------------
# Build up contributed packages that we will build from source
# --------------------------------------------------------------

su2_externals_INCLUDES=""
su2_externals_LIBS=""
su2_externals_LIBPTHREAD=""

# Tecplot
CONFIGURE_TECIO
if (test $enabletecio = yes); then
  su2_externals_INCLUDES="-DHAVE_TECIO $TECIO_INCLUDE $su2_externals_INCLUDES"
  if ( test x$have_MPI = xyes ); then
    su2_externals_INCLUDES="-DTECIOMPI $su2_externals_INCLUDES"
  fi
  su2_externals_LIBS="$TECIO_LIB $su2_externals_LIBS"
  # Check if boost directory exists, otherwise extract the sources
  if ( ! test -d "$srcdir/externals/tecio/boost" ); then
    AC_MSG_RESULT([Extracting boost source ...])
    if ( ! tar -xzf $srcdir/externals/tecio/boost.tar.gz -C $srcdir/externals/tecio/ ); then
	AC_MSG_ERROR([Extraction of boost sources to $srcdir/externals/tecio/boost using 'tar' failed ...])	    
    fi
  fi
   # Attempt to get pthread to work
   AC_CHECK_LIB(pthread,pthread_create,LIBPTHREAD="-lpthread")
  su2_externals_LIBPTHREAD="$LIBPTHREAD"
else
  su2_externals_INCLUDES="$su2_externals_INCLUDES"
fi

AM_CONDITIONAL(BUILD_TECIO, test x$enabletecio = xyes)
AM_CONDITIONAL(BUILD_TECIOMPI, test x$enabletecio = xyes -a x$have_MPI = xyes)
AC_CONFIG_FILES([externals/tecio/Makefile])
AC_CONFIG_FILES([externals/tecio/teciosrc/Makefile])
AC_CONFIG_FILES([externals/tecio/teciompisrc/Makefile])

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

# Parmetis
CONFIGURE_PARMETIS
if (test $enableparmetis = yes); then
  su2_externals_INCLUDES="-DHAVE_PARMETIS $PARMETIS_INCLUDE $su2_externals_INCLUDES"
  su2_externals_LIBS="$PARMETIS_LIB $su2_externals_LIBS"
else
  su2_externals_INCLUDES="$su2_externals_INCLUDES"
fi
AM_CONDITIONAL(BUILD_PARMETIS, test x$enableparmetis = xyes)
AC_CONFIG_FILES([externals/parmetis/Makefile])

# CGNS
CONFIGURE_CGNS
if (test $enablecgns = yes); then
  su2_externals_INCLUDES="-DHAVE_CGNS $CGNS_INCLUDE $su2_externals_INCLUDES"
  su2_externals_LIBS="$CGNS_LIB $su2_externals_LIBS"
else
  su2_externals_INCLUDES="$su2_externals_INCLUDES"
fi
AM_CONDITIONAL(BUILD_CGNS, test x$enablecgns = xyes)
AC_CONFIG_FILES([externals/cgns/Makefile])

AC_SUBST([su2_externals_INCLUDES])
AC_SUBST([su2_externals_LIBS])
AC_SUBST([su2_externals_LIBPTHREAD])

su2_externals_INCLUDES="-I\$(top_srcdir)/externals/CLI11 $su2_externals_INCLUDES"

##########################

# check for lib
Mutationpplib=libmutation++.a
Mutationppheader=mutation++.h
have_MPP="no"

if test "$with_Mutationpp_lib" != "no"
then
  	AC_CHECK_FILE([$with_Mutationpp_lib/$Mutationpplib],[have_MPP="yes"],[have_MPP="no"])
	if test "$have_MPP" == "no"
	then
		AC_MSG_ERROR([Mutation++ requested but library file not found.])
	fi
fi

# check for header
if test "$have_MPP" != "no"
then
	if test "$with_Mutationpp_include" != "no"
	then
		AC_CHECK_FILE([$with_Mutationpp_include/$Mutationppheader],[have_MPP='yes'],[have_MPP='no'])
		if test "$have_MPP" == "no"
		then
			AC_MSG_ERROR([Mutation++ requested but header file not found.])
		fi
	else
		have_MPP="no"
	fi
fi

if test "$have_MPP" != "no"
then
  	AM_CONDITIONAL(BUILD_MUTATIONPP,true,test x$enablempp = xyes -a x$have_MPP = xyes)
	MUTATIONPP_CXX="-DHAVE_MPP -I"$with_Mutationpp_include""
	MUTATIONPP_LD="$with_Mutationpp_lib"/$Mutationpplib
else
	AM_CONDITIONAL(BUILD_MUTATIONPP,false)
	MUTATIONPP_CXX=
	MUTATION_LD=
fi
AC_SUBST([MUTATIONPP_CXX])
AC_SUBST([MUTATIONPP_LD])

##########################

# Check if MKL is available, enable if 2019.0 or newer.

have_MKL="no"
if test "${MKLROOT+set}" = set; then
  FOUND_MKL_VERSION=$(grep -F "INTEL_MKL_VERSION" ${MKLROOT}/include/mkl_version.h | tr -dc '0-9')
  REQUIRED_MKL_VERSION=20190000
  if [[ "$FOUND_MKL_VERSION" -ge "$REQUIRED_MKL_VERSION" ]]; then
    CPPFLAGS="-DHAVE_MKL -DMKL_DIRECT_CALL_SEQ $CPPFLAGS"
    CXXFLAGS="-I${MKLROOT}/include $CXXFLAGS"
    LIBS="-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl $LIBS" 
    have_MKL="yes"
  fi
fi

###########################
# Determine what versions of the code to build

build_PY=yes

# direct differentation type
build_DIRECTDIFF=no
DIRECTDIFF_CXX=
DIRECTDIFF_LIBS=
# reverse mode type
build_REVERSE=no
REVERSE_CXX=
FORWARD_CXX=
CONFIGURE_CODI

AC_SUBST([DIRECTDIFF_CXX])
AC_SUBST([DIRECTDIFF_LIBS])
AC_SUBST([REVERSE_CXX])
AC_SUBST([REVERSE_LIBS])
# Set the variables
AM_CONDITIONAL([BUILD_NORMAL], [test $build_NORMAL != "no"])
AM_CONDITIONAL([BUILD_DIRECTDIFF], [test $build_DIRECTDIFF != "no"])
AM_CONDITIONAL([BUILD_REVERSE], [test $build_REVERSE != "no"])
if test $build_NORMAL == "no"
then
  build_DOT=no
  build_DEF=no
  build_SOL=no
  build_GEO=no
  build_PY=no
fi

if test $build_DIRECTDIFF == "yes"
then
  build_PY_WRAPPER=no
fi

if test $build_REVERSE == "yes"
then
  build_DOT=yes
fi
###########################

AM_CONDITIONAL([BUILD_CFD],[test $build_CFD != "no"])
AM_CONDITIONAL([BUILD_DOT],[test $build_DOT != "no"])
AM_CONDITIONAL([BUILD_DEF],[test $build_DEF != "no"])
AM_CONDITIONAL([BUILD_SOL],[test $build_SOL != "no"])
AM_CONDITIONAL([BUILD_GEO],[test $build_GEO != "no"])
AM_CONDITIONAL([BUILD_PY],[test $build_PY != "no"])
AM_CONDITIONAL([BUILD_PY_WRAPPER], [test $build_PY_WRAPPER != "no"])

###########################

AC_OUTPUT(Makefile \
        externals/Makefile \
	      Common/lib/Makefile \
	      SU2_CFD/obj/Makefile \
	      SU2_DOT/obj/Makefile \
	      SU2_DEF/obj/Makefile \
	      SU2_SOL/obj/Makefile \
	      SU2_GEO/obj/Makefile \
	      SU2_PY/Makefile \
	      SU2_PY/pySU2/Makefile)

abs_srcdir=`(cd $srcdir && pwd)`
abs_prefix=$prefix #`(cd $prefix && pwd)`


AC_MSG_RESULT([

-------------------------------------------------------------------------
|    ___ _   _ ___                                                      |
|   / __| | | |_  )   Release 7.3.0 'Blackbird'                         |
|   \\__ \\ |_| |/ /                                                      |
|   |___/\\___//___|   Suite                                             |
|                                                                       |
-------------------------------------------------------------------------
| SU2 Project Website: https://su2code.github.io                        |
|                                                                       |
| The SU2 Project is maintained by the SU2 Foundation                   |
| (http://su2foundation.org)                                            |
-------------------------------------------------------------------------
| Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)                |
|                                                                       |
| SU2 is free software; you can redistribute it and/or                  |
| modify it under the terms of the GNU Lesser General Public            |
| License as published by the Free Software Foundation; either          |
| version 2.1 of the License, or (at your option) any later version.    |
|                                                                       |
| SU2 is distributed in the hope that it will be useful,                |
| but WITHOUT ANY WARRANTY; without even the implied warranty of        |
| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |
| Lesser General Public License for more details.                       |
|                                                                       |
| You should have received a copy of the GNU Lesser General Public      |
| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |
-------------------------------------------------------------------------

Build Configuration Summary:

    Source code location: ${abs_srcdir}
    Install location: 	  ${abs_prefix}
    Version:              ${VERSION}
    C++ Compiler:         ${CXX}
    C Compiler:           ${CC}
    Preprocessor flags:   ${CPPFLAGS}
    Compiler flags:       ${CXXFLAGS}
    Linker flags:         ${LDFLAGS}
    MPI support:          $have_MPI
    Metis support:        $enablemetis
    Parmetis support:     $enableparmetis
    TecIO support:        $enabletecio
    CGNS support:         $enablecgns
    Mutation++ support:   $have_MPP
    MKL support:          $have_MKL
    Datatype support:
        double            $build_NORMAL
        codi_reverse      $build_CODI_REVERSE
        codi_forward      $build_CODI_FORWARD

    External includes:    $su2_externals_INCLUDES
    External libs:        $su2_externals_LIBS
    External libpthreads: $su2_externals_LIBPTHREAD

    Build SU2_CFD:        $build_CFD
    Build SU2_DOT:        $build_DOT
    Build SU2_DEF:        $build_DEF
    Build SU2_SOL:        $build_SOL
    Build SU2_GEO:        $build_GEO
    Build Py Wrapper:     $build_PY_WRAPPER

Please be sure to add the [\$]SU2_HOME and [\$]SU2_RUN environment variables,
and update your [\$]PATH (and [\$]PYTHONPATH if applicable) with [\$]SU2_RUN.

Based on the input to this configuration, add these lines to your .bashrc file:

export SU2_RUN="${abs_prefix}/bin"
export SU2_HOME="${abs_srcdir}"
export PATH=[\$]PATH:[\$]SU2_RUN
export PYTHONPATH=[\$]PYTHONPATH:[\$]SU2_RUN

])
