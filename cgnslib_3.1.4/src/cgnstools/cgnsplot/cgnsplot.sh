#!/bin/sh

# sh script to launch CGNS plot

dir=`dirname $0`

# source the setup script

for d in $dir $dir/cgnstools $dir/.. ; do
  if test -f $d/cgconfig ; then
    . $d/cgconfig
    break
  fi
done

# get the plotwish executable

plotwish=""
for d in $CG_BIN_DIR $dir $dir/cgnstools $dir/cgnsplot \
         /usr/local/bin /usr/local/bin/cgnstools ; do
  if test -x $d/plotwish ; then
    plotwish=$d/plotwish
    break
  fi
  if test -x $d/cgnswish/plotwish ; then
    plotwish=$d/cgnswish/plotwish
    break
  fi
done
if test -z "$plotwish" ; then
  echo "Error: plotwish executable not found"
  exit 1
fi

# find the cgnsplot tcl script

cgnsplot=""
for d in $CG_LIB_DIR $dir $dir/cgnstools $dir/cgnsplot \
         /usr/local/share /usr/local/share/cgnstools ; do
  if test -f $d/cgnsplot.tcl ; then
    cgnsplot=$d/cgnsplot.tcl
    break
  fi
done
if test -z "$cgnsplot" ; then
  echo "Error: cgnsplot.tcl script not found"
  exit 1
fi

# check that display is set

#if test -z "$DISPLAY" ; then
#  echo "Error: DISPLAY environment variable not set"
#  exit 1
#fi

# execute

if test "$#" = 0 ; then
  exec $plotwish $cgnsplot
else
  exec $plotwish $cgnsplot "$@"
fi
