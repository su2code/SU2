# -------------------------------------------------------------
# amgio
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_AMG],
[
  AC_ARG_ENABLE(amg,
                AC_HELP_STRING(),
		[case "${enableval}" in
		  yes)  enableamg=yes ;;
		   no)  enableamg=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-amg) ;;
		 esac],
		 [enableamg=yes])

  # The amgio API is distributed with SU2, so we don't have to guess
  # where it might be installed...
  if (test $enableamg = yes); then

    # amg platform-specific compiler flags
    AMGIO_CPPFLAGS=""
    case "${host_os}" in
      *linux*)
        AMGIO_CPPFLAGS="-DLINUX $AMGIO_CPPFLAGS"
        AC_CHECK_SIZEOF([void *])
        if (test $ac_cv_sizeof_void_p = 8); then
          AMGIO_CPPFLAGS="-DLINUX64 $AMGIO_CPPFLAGS"
        fi
        ;;

      *darwin*)
        AMGIO_CPPFLAGS="-DDARWIN -DLONGIS64 $AMGIO_CPPFLAGS"
        ;;

      *)
        AC_MSG_RESULT([>>> Unrecognized AMG platform <<<])
        ;;
    esac

     AMGIO_INCLUDE="-I\$(top_srcdir)/externals/amgio"
     AMGIO_LIB="\$(top_builddir)/externals/amgio/libMeshb.a"
     AC_DEFINE(HAVE_AMG, 1, [Flag indicating whether the library will be compiled with GMF support])
     AC_MSG_RESULT(<<< Configuring library with GMF support >>>)

     # look for thread-local storage
     #AX_TLS
 else
     AMGIO_INCLUDE=""
     AMGIO_LIB=""
     enableamg=no
  fi

  AC_SUBST(AMGIO_INCLUDE)
  AC_SUBST(AMGIO_LIB)
  AC_SUBST(AMGIO_CPPFLAGS)
])
