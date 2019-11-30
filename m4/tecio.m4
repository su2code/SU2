dnl -------------------------------------------------------------
dnl tecio
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_TECIO],
[
  AC_ARG_ENABLE(tecio,
                AC_HELP_STRING([--enable-tecio],
                               [build with Tecplot TecIO API support (from source)]),
		[case "${enableval}" in
		  yes)  enabletecio=yes ;;
		   no)  enabletecio=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-tecio) ;;
		 esac],
		 [enabletecio=yes])


  # The TECIO API is distributed with libmesh, so we don't have to guess
  # where it might be installed...
  if (test $enabletecio = yes); then

    # tecio platform-specific compiler flags
    TECIO_CPPFLAGS=""
    case "${host_os}" in
      *linux*)
	TECIO_CPPFLAGS="-DLINUX $TECIO_CPPFLAGS"
	AC_CHECK_SIZEOF([void *])
	if (test $ac_cv_sizeof_void_p = 8); then
	  TECIO_CPPFLAGS="-DLINUX64 $TECIO_CPPFLAGS"
	fi
	;;

      *darwin*)
	TECIO_CPPFLAGS="-DDARWIN -DMAC64 $TECIO_CPPFLAGS"
        ;;

        *)
	AC_MSG_RESULT([>>> Unrecognized TecIO platform, see externals/tecio/Runmake for hints on how to extend <<<])
	;;
    esac


     if (test $have_MPI = yes); then
         TECIO_INCLUDE="-I\$(top_srcdir)/externals/tecio/teciompisrc"
         TECIO_LIB="\$(top_builddir)/externals/tecio/teciompisrc/libteciompi.a"
     else
         TECIO_INCLUDE="-I\$(top_srcdir)/externals/tecio/teciosrc"
         TECIO_LIB="\$(top_builddir)/externals/tecio/teciosrc/libtecio.a"
     fi
     AC_DEFINE(HAVE_TECPLOT_API, 1, [Flag indicating whether the library will be compiled with Tecplot TecIO API support])
     AC_DEFINE(HAVE_TECPLOT_API_112, 1, [Flag indicating tecplot API understands newer features])
     AC_MSG_RESULT(<<< Configuring library with Tecplot TecIO support >>>)
     have_tecio=yes
  else
     TECIO_INCLUDE=""
     enabletecio=no
     have_tecio=no
  fi

  AC_SUBST(TECIO_INCLUDE)
  AC_SUBST(TECIO_CPPFLAGS)
])
