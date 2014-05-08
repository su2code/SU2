# This file is part of Autoconf.                       -*- Autoconf -*-
# Go language support.
# Copyright (C) 2011-2012 Free Software Foundation, Inc.

# This file is part of Autoconf.  This program is free
# software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Under Section 7 of GPL version 3, you are granted additional
# permissions described in the Autoconf Configure Script Exception,
# version 3.0, as published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# and a copy of the Autoconf Configure Script Exception along with
# this program; see the files COPYINGv3 and COPYING.EXCEPTION
# respectively.  If not, see <http://www.gnu.org/licenses/>.

# Go support contributed by Ian Lance Taylor.

# This currently only supports gccgo, not 6g/8g/5g.

# ------------------- #
# Language selection.
# ------------------- #

# AC_LANG(Go)
# -----------
AC_LANG_DEFINE([Go], [go], [GO], [GOC], [],
[ac_ext=go
ac_compile='$GOC -c $GOFLAGS conftest.$ac_ext >&AS_MESSAGE_LOG_FD'
ac_link='$GOC -o conftest$ac_exeext $GOFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&AS_MESSAGE_LOG_FD'
ac_compiler_gnu=yes
])

# AC_LANG_GO
# ----------
AU_DEFUN([AC_LANG_GO], [AC_LANG(Go)])

# ------------------- #
# Producing programs.
# ------------------- #

# AC_LANG_PROGRAM(Go)([PROLOGUE], [BODY])
# ---------------------------------------
m4_define([AC_LANG_PROGRAM(Go)],
[package main
$1
func main() {
$2
}])

# _AC_LANG_IO_PROGRAM(Go)
# -----------------------
# Produce source that performs I/O.
m4_define([_AC_LANG_IO_PROGRAM(Go)],
[AC_LANG_PROGRAM([import ( "fmt"; "os" )],
[f, err := os.Open("conftest.out", os.O_CREATE|os.O_WRONLY, 0777)
 if err != nil {
	fmt.Println(err)
	os.Exit(1)
 }
 if err = f.Close(); err != nil {
	fmt.Println(err)
	os.Exit(1)
 }
 os.Exit(0)
])])

# AC_LANG_CALL(Go)(PROLOGUE, FUNCTION)
# ------------------------------------
# Avoid conflicting decl of main.
m4_define([AC_LANG_CALL(Go)],
[AC_LANG_PROGRAM([$1
m4_if([$2], [main], ,
[func $2()])],[$2()])])

# AC_LANG_FUNC_LINK_TRY(Go)(FUNCTION)
# -----------------------------------
# Try to link a program which calls FUNCTION.
m4_define([AC_LANG_FUNC_LINK_TRY(Go)],
[AC_LANG_PROGRAM(
[func $1() int
var f = $1
], [return f()])])

# AC_LANG_BOOL_COMPILE_TRY(Go)(PROLOGUE, EXPRESSION)
# --------------------------------------------------
# Return a program which is valid if EXPRESSION is nonzero.
m4_define([AC_LANG_BOOL_COMPILE_TRY(Go)],
[AC_LANG_PROGRAM([$1], [var test_array @<:@1 - 2 * !($2)@:>@int
test_array @<:@0@:>@ = 0
])])

# AC_LANG_INT_SAVE(Go)(PROLOGUE, EXPRESSION)
# ------------------------------------------
m4_define([AC_LANG_INT_SAVE(Go)],
[AC_LANG_PROGRAM([$1
import (
	"fmt"
	"os"
)
],
[f, err := os.Open("conftest.val", os.O_CREATE|os.O_WRONLY, 0777)
 if err != nil {
	os.Exit(1)
 }
 if $2 < 0 {
	int64 i = int64($2)
	if i != $2 {
		os.Exit(1)
	}
	if _, err := fmt.Print(f, i); err != nil {
		os.Exit(1)
	}
 } else {
	uint64 i = uint64($2)
	if i != $2 {
		os.Exit(1)
	}
	if _, err := fmt.Print(f, i); err != nil {
		os.Exit(1)
	}
 }
 if err = f.Close(); err != nil {
	os.Exit(1)
 }
 os.Exit(0);
])])

# ---------------------- #
# Looking for compilers. #
# ---------------------- #

# AC_LANG_COMPILER(Go)
# --------------------
AC_DEFUN([AC_LANG_COMPILER(Go)],
[AC_REQUIRE([AC_PROG_GO])])

# AC_PROG_GO
# ----------
AN_MAKEVAR([GOC], [AC_PROG_GO])
AN_PROGRAM([gccgo], [AC_PROG_GO])
AC_DEFUN([AC_PROG_GO],
[AC_LANG_PUSH(Go)dnl
AC_ARG_VAR([GOC],   [Go compiler command])dnl
AC_ARG_VAR([GOFLAGS], [Go compiler flags])dnl
_AC_ARG_VAR_LDFLAGS()dnl
m4_ifval([$1],
      [AC_CHECK_TOOLS(GOC, [$1])],
[AC_CHECK_TOOL(GOC, gccgo)
if test -z "$GOC"; then
  if test -n "$ac_tool_prefix"; then
    AC_CHECK_PROG(GOC, [${ac_tool_prefix}gccgo], [$ac_tool_prefix}gccgo])
  fi
fi
if test -z "$GOC"; then
  AC_CHECK_PROG(GOC, gccgo, gccgo, , , false)
fi
])

# Provide some information about the compiler.
_AS_ECHO_LOG([checking for _AC_LANG compiler version])
set X $ac_compile
ac_compiler=$[2]
_AC_DO_LIMIT([$ac_compiler --version >&AS_MESSAGE_LOG_FD])
m4_expand_once([_AC_COMPILER_EXEEXT])[]dnl
m4_expand_once([_AC_COMPILER_OBJEXT])[]dnl
GOFLAGS="-g -O2"
AC_LANG_POP(Go)dnl
])# AC_PROG_GO
