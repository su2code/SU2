dnl -------------------------------------------------------------
dnl Metis
dnl -------------------------------------------------------------
AC_DEFUN([CONFIGURE_METIS],
[
  AC_ARG_ENABLE(metis,
                AC_HELP_STRING([--enable-metis],
                               [build with Metis graph partitioning suppport]),
		[case "${enableval}" in
		  yes)  enablemetis=yes ;;
		   no)  enablemetis=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-metis) ;;
		 esac],
		 [enablemetis=yes])

  dnl The METIS API is distributed with SU2, so we don't have to guess
  dnl where it might be installed...
  if (test $enablemetis = yes); then
     METIS_INCLUDE="-DMETIS_5 -I\$(top_srcdir)/contrib/metis/include"
     METIS_LIB="\$(top_builddir)/contrib/metis/libmetis.a"
     AC_DEFINE(HAVE_METIS, 1, [Flag indicating whether the library will be compiled with Metis support])
     AC_MSG_RESULT(<<< Configuring library with Metis support >>>)

     dnl look for thread-local storage
     AX_TLS
 else
     METIS_INCLUDE=""
     METIS_LIB=""
     enablemetis=no
  fi

  AC_SUBST(METIS_INCLUDE)
  AC_SUBST(METIS_LIB)
])
