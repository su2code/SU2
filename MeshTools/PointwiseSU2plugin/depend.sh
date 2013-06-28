#!/bin/sh

case "$machine" in
"linux_x86_64")
  DEP="g++ -MM -MG -m64 -mtune=opteron"
  ;;
"linux" )
  DEP="g++ -MM -MG"
  ;;
"macosx" )
  DEP="g++ -MM -MG"
  ;;
*)
  echo "Unsupported machine type: $machine"
  exit 1
  ;;
esac

DIR="$1"
shift 1
case "$DIR" in
"" | ".")
  $DEP "$@" | sed -e 's@^\(.*\)\.o:@\1.d \1.o:@'
  ;;
*)
  $DEP "$@" | sed -e "s@^\(.*\)\.o:@$DIR/\1.d $DIR/\1.o:@"
  ;;
esac
