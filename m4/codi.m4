AC_DEFUN([CONFIGURE_CODI],
[
    AC_ARG_ENABLE(codi-reverse,
        AS_HELP_STRING([--enable-codi-reverse], [build executables with codi reverse datatype (default = no)]),
        [build_CODI_REVERSE="yes"], [build_CODI_REVERSE="no"])
    AC_ARG_ENABLE(codi-forward,
        AS_HELP_STRING([--enable-codi-forward], [build executables with codi forward datatype (default = no)]),
        [build_CODI_FORWARD="yes"], [build_CODI_FORWARD="no"])
    AC_ARG_WITH(codi-include,
        AS_HELP_STRING([--with-CODI-include[=ARG]], [CODI include directory, ARG = path to codi.hpp]),
        [with_CODI_include=$withval], [with_CODI_include="no"])
        REVERSE_CXX=
        FORWARD_CXX=
        CODIheader=codi.hpp

        if test "$build_CODI_REVERSE" == "yes" || test "$build_CODI_FORWARD" == "yes"
        then
          AC_CHECK_FILE([$with_CODI_include/$CODIheader],[have_CODI='yes'],[have_CODI='no'])
          if test "$have_CODI" == "no"
          then
            AC_MSG_ERROR([CODI requested but header file not found.])
          fi
        fi

        if test "$build_CODI_FORWARD" == "yes"
        then
           DIRECTDIFF_CXX="-std=c++11 -DCODI_FORWARD_TYPE -I"$with_CODI_include""
           build_DIRECTDIFF=yes
           build_REVERSE=no
           build_NORMAL=no
        elif test "$build_CODI_REVERSE" == "yes"
        then
           REVERSE_CXX="-std=c++11 -DCODI_REVERSE_TYPE -I"$with_CODI_include""
           build_REVERSE=yes
           build_DIRECTDIFF=no
           build_NORMAL=no
        else
           build_REVERSE=no
           build_DIRECTDIFF=no
        fi
])
