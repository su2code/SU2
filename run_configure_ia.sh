#!/bin/sh
#$1 = 0/1 toggle between AVX2/AVX512
#$2 = build type BASE/MKL/LIBXSMM

source /opt/intel/compilers_and_libraries_2016/linux/bin/compilervars.sh intel64
source /opt/intel/impi_latest/bin64/mpivars.sh

LIBXSMM_include_dir=/home/dmudiger/tools/libxsmm_github/include
LIBXSMM_lib_dir=/home/dmudiger/tools/libxsmm_github/lib

#single-node
SU2_dir="$PWD"
config="_HOM_$2"
INSTALL_DIR=$SU2_dir/$config
HOME_dir="$HOME"

openmp=1
sim=0
autovec=0

CPU=1
  
INCFLAGS="-I$SU2_dir/externals"


FLAGS="-DTIME -DPROFILE -DIPM_STATS -D_SU2_ -Wall -Wextra -DNDEBUG -g -O2"
if [ $1 -eq 0 ]; then
  FLAGS+=" -march=core-avx2"
else
  FLAGS+=" -xCOMMON-AVX512"
fi
FLAGS+=" -qopenmp -mkl"

if [ "$2" == "MKL" ]; then
  FLAGS+=" -DHAVE_MKL"
fi 
if [ "$2" == "LIBXSMM" ]; then
  FLAGS+=" -DHAVE_LIBXSMM"
  INCFLAGS+=" -I$LIBXSMM_include_dir"
  LDFLAGS+=" -L$LIBXSMM_lib_dir -lxsmm"
  LIBFLAGS+=" -lxsmm"
fi

CFLAGS="$FLAGS" 
CXXFLAGS="$FLAGS -std=c++11"

cd $SU2_dir

build_cmd="./configure --enable-mpi --with-cc=mpiicc --with-cxx=mpiicpc CXXFLAGS=\"$CXXFLAGS $INCFLAGS\" CFLAGS=\"$CFLAGS $INCFLAGS\" LDFLAGS=\"$LDFLAGS\" LIBS=\"$LIBFLAGS\" --with-metis-cppflags=\"-O3\" --with-parmetis-cppflags=\"-O3\" --prefix=$INSTALL_DIR --exec-prefix=$INSTALL_DIR"

echo "$build_cmd"
read -p "Continue ?"

eval "$build_cmd"

#make clean
make -j 
make install
