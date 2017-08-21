# -------------------------------------------------------------
# CGNS
# -------------------------------------------------------------
AC_DEFUN([CONFIGURE_CGNS],
[
  AC_ARG_ENABLE(cgns,
                AC_HELP_STRING([--enable-cgns],
                               [build with CFD General Notation System (CGNS) standard support (ADF only)]),
		[case "${enableval}" in
		  yes)  enablecgns=yes ;;
		   no)  enablecgns=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-cgns) ;;
		 esac],
		 [enablecgns=yes])

  # The CGNS lib is distributed with SU2, so we don't have to guess
  # where it might be installed...
  if (test $enablecgns = yes); then

    # look for CGNS build cppflags by honoring the --with-cgns-cppflags="..." flag,
    # defaulting to what we know works
    AC_ARG_WITH([cgns-cppflags],
                 AC_HELP_STRING([--with-cgns-cppflags="-fPIC"],
                                [Specific CGNS C Preprocessor flags to use]),
                 [SU2_CGNS_CPPFLAGS="$withval"],
                 [SU2_CGNS_CPPFLAGS=""])


     CGNS_INCLUDE="-I\$(top_srcdir)/externals/cgns -I\$(top_srcdir)/externals/cgns/adf"
     CGNS_LIB="\$(top_builddir)/externals/cgns/libcgns.a"
     AC_DEFINE(HAVE_CGNS, 1, [Flag indicating whether the library will be compiled with CGNS support])
     AC_MSG_RESULT(<<< Configuring library with CGNS support >>>)

     # look for thread-local storage
     #AX_TLS
 else
     CGNS_INCLUDE=""
     CGNS_LIB=""
     SU2_CGNS_CPPFLAGS=""
     enablecgns=no
  fi

  AC_SUBST(CGNS_INCLUDE)
  AC_SUBST(CGNS_LIB)
  AC_SUBST(SU2_CGNS_CPPFLAGS)
])
