AC_DEFUN([CONFIGURE_CODI],
[
    AC_ARG_ENABLE(codi-reverse,
        AS_HELP_STRING([--enable-codi-reverse], [build executables with codi reverse datatype (default = no)]),
        [build_CODI_REVERSE="yes"], [build_CODI_REVERSE="no"])
    AC_ARG_ENABLE(codi-forward,
        AS_HELP_STRING([--enable-codi-forward], [build executables with codi forward datatype (default = no)]),
        [build_CODI_FORWARD="yes"], [build_CODI_FORWARD="no"])

        CODIheader=${srcdir}/externals/codi/include/codi.hpp
        AMPIlib=${srcdir}/externals/adjointmpi/libAMPI.a
        AMPIheader=${srcdir}/externals/adjointmpi/include/ampi.h

        if test "$build_CODI_REVERSE" == "yes" || test "$build_CODI_FORWARD" == "yes"
        then
          AC_CHECK_FILE([$CODIheader],[have_CODI='yes'],[have_CODI='no'])
          if test "$have_CODI" == "no"
          then
            AC_MSG_ERROR([CODI header was not found in externals/CoDi/include. Use 'preconfigure.py --autodiff' to download it.])
          fi
        fi

        if test "$build_CODI_FORWARD" == "yes"
        then
           DIRECTDIFF_CXX="-std=c++0x -DCODI_FORWARD_TYPE -I\$(top_srcdir)/externals/codi/include"
           build_DIRECTDIFF=yes
           build_REVERSE=no
           build_NORMAL=no
        elif test "$build_CODI_REVERSE" == "yes"
        then
           REVERSE_CXX="-std=c++0x -DCODI_REVERSE_TYPE -I\$(top_srcdir)/externals/codi/include"
           if test "$enablempi" == "yes"
           then
              AC_CHECK_FILE([$AMPIlib], [have_AMPIlib='yes'], [have_AMPIlib='no'])
              AC_CHECK_FILE([$AMPIheader], [have_AMPIheader='yes'], [have_AMPIheader='no'])
              if test "$have_AMPIlib" == "no"
              then
                AC_MSG_ERROR([AMPI library not found in externals/adjointmpi. Use 'preconfigure.py --autodiff --enable-mpi' to download it.])
              fi
              if test "$have_AMPIheader" == "no"
              then
                AC_MSG_ERROR([AMPI header not found in externals/adjointmpi/include.  Use 'preconfigure.py --autodiff --enable-mpi' to download it.])
              fi
              REVERSE_CXX=$REVERSE_CXX" -I\$(top_srcdir)/externals/adjointmpi/include"
              REVERSE_LIBS="\$(top_srcdir)/externals/adjointmpi/libAMPI.a"
           fi
           build_REVERSE=yes
           build_DIRECTDIFF=no
           build_NORMAL=no
        fi
])
