#!/bin/sh

# sh script to launch CGNS calculator

dir=`dirname $0`

# source the setup script

for d in $dir $dir/cgnstools $dir/.. ; do
  if test -f $d/cgconfig ; then
    . $d/cgconfig
    break
  fi
done

# get the calcwish executable

calcwish=""
for d in $CG_BIN_DIR $dir $dir/cgnstools $dir/cgnscalc \
         /usr/local/bin /usr/local/bin/cgnstools ; do
  if test -x $d/calcwish ; then
    calcwish=$d/calcwish
    break
  fi
  if test -x $d/cgnswish/calcwish ; then
    calcwish=$d/cgnswish/calcwish
    break
  fi
done
if test -z "$calcwish" ; then
  echo "Error: calcwish executable not found"
  exit 1
fi

# find the cgnscalc tcl script

cgnscalc=""
for d in $CG_LIB_DIR $dir $dir/cgnstools $dir/cgnscalc \
         /usr/local/share /usr/local/share/cgnstools ; do
  if test -f $d/cgnscalc.tcl ; then
    cgnscalc=$d/cgnscalc.tcl
    break
  fi
done
if test -z "$cgnscalc" ; then
  echo "Error: cgnscalc.tcl script not found"
  exit 1
fi

# check that display is set

#if test -z "$DISPLAY" ; then
#  echo "Error: DISPLAY environment variable not set"
#  exit 1
#fi

# execute

if test "$#" = 0 ; then
  exec $calcwish $cgnscalc
else
  exec $calcwish $cgnscalc "$@"
fi

