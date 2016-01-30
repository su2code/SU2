AC_DEFUN([CONFIGURE_ADOLC],
[
    AC_ARG_ENABLE(adolc-reverse,
        AS_HELP_STRING([--enable-adolc-reverse], [build executables with adolc reverse datatype (default = no)]),
        [build_ADOLC_REVERSE="yes"], [build_ADOLC_REVERSE="no"])
    AC_ARG_ENABLE(adolc-forward,
        AS_HELP_STRING([--enable-adolc-forward], [build executables with adolc forward datatype (default = no)]),
        [build_ADOLC_FORWARD="yes"], [build_ADOLC_FORWARD="no"])
        ADOLC_REVERSE_CXX=
        ADOLC_FORWARD_CXX=
        ADOLC_LIBS=
        ADOLC_VERSION=

        if test "$have_MPI" == "yes" && test "$build_ADOLC_REVERSE" == "yes"
        then
          ADOLC_VERSION="adolc_ampi = 2.5.3-trunk"
        else
          ADOLC_VERSION="adolc = 2.5.3-trunk"
        fi

        if test "$build_ADOLC_REVERSE" == "yes" || test "$build_ADOLC_FORWARD" == "yes"
        then
          PKG_CHECK_MODULES([ADOLC], [${ADOLC_VERSION}], [have_ADOLC='yes'], [have_ADOLC='no'])
          if test "$have_ADOLC" == "no"
          then
            AC_MSG_ERROR([ADOLC requested but library file not found.])
          fi
        fi

        if test "$build_ADOLC_FORWARD" == "yes"
        then
           DIRECTDIFF_CXX="-std=c++11 -DADOLC_FORWARD_TYPE "$ADOLC_CFLAGS""
           DIRECTDIFF_LIBS=$ADOLC_LIBS
           build_DIRECTDIFF=yes
           build_NORMAL=no
        elif test "$build_ADOLC_REVERSE" == "yes"
        then
           REVERSE_CXX="-std=c++11 -DADOLC_REVERSE_TYPE "$ADOLC_CFLAGS""
           REVERSE_LIBS=$ADOLC_LIBS
           build_REVERSE=yes
           build_NORMAL=no
           else
              build_REVERSE=no
              build_DIRECTDIFF=no
           fi
])
