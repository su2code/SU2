# -------------------------------------------------------------
# Metis
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_METIS],
[
  AC_ARG_ENABLE(metis,
                AC_HELP_STRING(),
		[case "${enableval}" in
		  yes)  enablemetis=yes ;;
		   no)  enablemetis=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-metis) ;;
		 esac],
		 [enablemetis=yes])

  # Trump --enable-parmetis with --disable-mpi
  if (test "x$enablempi" = xno); then
    enablemetis=no
  fi

  # The METIS API is distributed with SU2, so we don't have to guess
  # where it might be installed...
  if (test $enablemetis = yes); then

    # look for METIS build cppflags by honoring the --with-metis-cppflags="..." flag,
    # defaulting to what we know works
    AC_ARG_WITH([metis-cppflags],
                 AC_HELP_STRING([--with-metis-cppflags="-D_FILE_OFFSET_BITS=64 -DNDEBUG -DNDEBUG2 -DHAVE_EXECINFO_H -DHAVE_GETLINE"],
                                [Specific METIS C Preprocessor flags to use]),
                 [SU2_METIS_CPPFLAGS="$withval"],
                 [SU2_METIS_CPPFLAGS="-D_FILE_OFFSET_BITS=64 -DNDEBUG -DNDEBUG2 -DHAVE_EXECINFO_H -DHAVE_GETLINE"])


     METIS_INCLUDE="-I\$(top_srcdir)/externals/metis/include"
     METIS_LIB="\$(top_builddir)/externals/metis/libmetis.a"
     AC_DEFINE(HAVE_METIS, 1, [Flag indicating whether the library will be compiled with Metis support])
     AC_MSG_RESULT(<<< Configuring library with Metis support >>>)

     # look for thread-local storage
     #AX_TLS
 else
     METIS_INCLUDE=""
     METIS_LIB=""
     SU2_METIS_CPPFLAGS=""
     enablemetis=no
  fi

  AC_SUBST(METIS_INCLUDE)
  AC_SUBST(METIS_LIB)
  AC_SUBST(SU2_METIS_CPPFLAGS)
])
