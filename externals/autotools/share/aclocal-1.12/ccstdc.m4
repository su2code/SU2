## ----------------------------------------- ##             -*- Autoconf -*-
## ANSIfy the C compiler whenever possible.  ##
## From Franc,ois Pinard                     ##
## ----------------------------------------- ##

# Copyright (C) 1996-2012 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This was merged into AC_PROG_CC in Autoconf.

AU_DEFUN([AM_PROG_CC_STDC],
[AC_PROG_CC
AC_DIAGNOSE([obsolete], [$0:
	your code should no longer depend upon 'am_cv_prog_cc_stdc', but upon
	'ac_cv_prog_cc_stdc'.  Remove this warning and the assignment when
	you adjust the code.  You can also remove the above call to
	AC_PROG_CC if you already called it elsewhere.])
am_cv_prog_cc_stdc=$ac_cv_prog_cc_stdc
])
AU_DEFUN([fp_PROG_CC_STDC])
