# -------------------------------------------------------------
# Parmetis
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_PARMETIS],
[
  AC_ARG_ENABLE(parmetis,
                AS_HELP_STRING([--disable-parmetis],
                               [build without Parmetis parallel graph partitioning suppport]),
		[case "${enableval}" in
		  yes)  enableparmetis=yes ;;
		   no)  enableparmetis=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-parmetis) ;;
		 esac],
		 [enableparmetis=yes])

  # Trump --enable-parmetis with --disable-mpi
  if (test "x$enablempi" = xno); then
    enableparmetis=no
  fi


  # The PARMETIS API is distributed with SU2, so we don't have to guess
  # where it might be installed...
  if (test $enableparmetis = yes); then
     PARMETIS_INCLUDE="-I\$(top_srcdir)/externals/parmetis/include"
     PARMETIS_LIB="\$(top_builddir)/externals/parmetis/libparmetis.a"
     AC_DEFINE(HAVE_PARMETIS, 1, [Flag indicating whether the library will be compiled with Parmetis support])
     AC_MSG_RESULT(<<< Configuring library with Parmetis support >>>)
  else
     PARMETIS_INCLUDE=""
     PARMETIS_LIB=""
     enableparmetis=no
  fi

  AC_SUBST(PARMETIS_INCLUDE)
  AC_SUBST(PARMETIS_LIB)
])
