#!/bin/sh

# sh script to launch unit convertor

dir=`dirname $0`

# source the setup script

for d in $dir $dir/cgnstools $dir/.. ; do
  if test -f $d/cgconfig ; then
    . $d/cgconfig
    break
  fi
done

# The normal wish will work here, but calcwish should
# be available, and may also be used

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

# find the unitconv tcl script

unitconv=""
for d in $CG_LIB_DIR $dir $dir/cgnstools $dir/cgnscalc \
         /usr/local/share /usr/local/share/cgnstools ; do
  if test -f $d/unitconv.tcl ; then
    unitconv=$d/unitconv.tcl
    break
  fi
done

if test -z "$unitconv" ; then
  echo "Error: unitconv.tcl script not found"
  exit 1
fi

# check that display is set

if test -z "$DISPLAY" ; then
  echo "Error: DISPLAY environment variable not set"
  exit 1
fi

# execute

exec $calcwish $unitconv

