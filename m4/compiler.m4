# -------------------------------------------------------------
# -------------------------------------------------------------
AC_DEFUN([SU2_SET_COMPILERS],
[
  # --------------------------------------------------------------
  # look for a decent C++ compiler or honor --with-cxx=...
  CXX_TRY_LIST="g++ icpc icc pgCC c++"

     # -------------------------------------------------------------------
     # MPI -- disabled by default.  Check for it now so we can be somewhat
     #                             smart about which compilers to look for
     # -------------------------------------------------------------------
     AC_ARG_ENABLE(mpi,
                   AS_HELP_STRING([--enable-mpi],
                                  [build with MPI message passing support]),
   		   [case "${enableval}" in
   		     yes)  enablempi=yes ;;
   		      no)  enablempi=no ;;
    		       *)  AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
   		    esac],
   		    [enablempi=no])

  have_MPI="no"

  if  (test "$enablempi" != no) ; then
    have_MPI="yes"
    CPPFLAGS="-DHAVE_MPI $CPPFLAGS"
    CXX_TRY_LIST="mpicxx mpiCC mpicc $CXX_TRY_LIST"
  else
    AC_MSG_RESULT(>>> MPI support disabled by default <<<)
  fi

  AC_ARG_WITH([cxx],
  	    AS_HELP_STRING([--with-cxx=CXX],
                             [C++ compiler to use]),
  	    [CXX="$withval"],
  	    [])

  # --------------------------------------------------------------
  # Determines a C++ compiler to use.  First checks if the variable CXX is
  # already set.  If not, then searches under g++, c++, and other names.
  # --------------------------------------------------------------
  AC_PROG_CXX([$CXX_TRY_LIST])
  # --------------------------------------------------------------



  # --------------------------------------------------------------
  # look for a decent C compiler or honor --with-cc=...
  CC_TRY_LIST="gcc icc pgcc cc"
  if  (test "$enablempi" != no) ; then
    CC_TRY_LIST="mpicc $CC_TRY_LIST"
  fi
  AC_ARG_WITH([cc],
  	    AS_HELP_STRING([--with-cc=CC],
                             [C compiler to use]),
  	    [CC="$withval"],
  	    [])

  # --------------------------------------------------------------
  # Determine a C compiler to use.  If CC is not already set, checks for
  # gcc, cc, and other C compilers.  Then sets the CC variable to the result.
  # --------------------------------------------------------------
  AC_PROG_CC([$CC_TRY_LIST])
  # --------------------------------------------------------------
])
