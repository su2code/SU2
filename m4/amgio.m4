# -------------------------------------------------------------
# Inria
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_INRIA],
[
  AC_ARG_ENABLE(inria,
                AC_HELP_STRING(),
		[case "${enableval}" in
		  yes)  enableinria=yes ;;
		   no)  enableinria=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-inria) ;;
		 esac],
		 [enableinria=yes])

  # The amgio API is distributed with SU2, so we don't have to guess
  # where it might be installed...
  if (test $enableinria = yes); then

    # inria platform-specific compiler flags
    INRIA_CPPFLAGS=""
    case "${host_os}" in
      *linux*)
        INRIA_CPPFLAGS="-DLINUX $INRIA_CPPFLAGS"
        AC_CHECK_SIZEOF([void *])
        if (test $ac_cv_sizeof_void_p = 8); then
          INRIA_CPPFLAGS="-DLINUX64 $INRIA_CPPFLAGS"
        fi
        ;;

      *darwin*)
        INRIA_CPPFLAGS="-DDARWIN -DLONGIS64 $INRIA_CPPFLAGS"
        ;;

      *)
        AC_MSG_RESULT([>>> Unrecognized Inria platform <<<])
        ;;
    esac

     INRIA_INCLUDE="-I\$(top_srcdir)/externals/amgio/sources"
     INRIA_LIB="\$(top_builddir)/externals/amgio/libMeshb.a"
     AC_DEFINE(HAVE_INRIA, 1, [Flag indicating whether the library will be compiled with GMF support])
     AC_MSG_RESULT(<<< Configuring library with GMF support >>>)

     # look for thread-local storage
     #AX_TLS
 else
     INRIA_INCLUDE=""
     INRIA_LIB=""
     enableinria=no
  fi

  AC_SUBST(INRIA_INCLUDE)
  AC_SUBST(INRIA_LIB)
  AC_SUBST(INRIA_CPPFLAGS)
])
