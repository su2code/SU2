#!/bin/sh

# sh script to launch CGNS File viewer/editor

dir=`dirname $0`

# source the setup script

for d in $dir $dir/cgnstools $dir/.. ; do
  if test -f $d/cgconfig ; then
    . $d/cgconfig
    break
  fi
done

# get the cgiowish executable

cgiowish=""
for d in $CG_BIN_DIR $dir $dir/cgnstools $dir/cgnsview \
         /usr/local/bin /usr/local/bin/cgnstools ; do
  if test -x $d/cgiowish ; then
    cgiowish=$d/cgiowish
    break
  fi
  if test -x $d/cgnswish/cgiowish ; then
    cgiowish=$d/cgnswish/cgiowish
    break
  fi
done
if test -z "$cgiowish" ; then
  echo "Error: cgiowish executable not found"
  exit 1
fi

# find the cgnsview tcl script

cgnsview=""
for d in $CG_LIB_DIR $dir $dir/cgnstools $dir/cgnsview \
         /usr/local/share /usr/local/share/cgnstools ; do
  if test -f $d/cgnsview.tcl ; then
    cgnsview=$d/cgnsview.tcl
    break
  fi
done
if test -z "$cgnsview" ; then
  echo "Error: cgnsview.tcl script not found"
  exit 1
fi

# check that display is set

#if test -z "$DISPLAY" ; then
#  echo "Error: DISPLAY environment variable not set"
#  exit 1
#fi

# execute

if test "$#" = 0 ; then
  exec $cgiowish $cgnsview
else
  exec $cgiowish $cgnsview "$@"
fi
