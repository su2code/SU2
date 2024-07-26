#!bin/bash

rm -r build

export PKG_CONFIG_PATH="/home/areen/Programs/OpenBLAS_install/lib:$PKG_CONFIG_PATH"
export PKG_CONFIG_PATH="/usr/local/cuda/lib64:$PKG_CONFIG_PATH"
export CUDA_PATH="/usr/local/cuda:$CUDA_PATH"
./meson.py build --buildtype=debug -Denable-pywrapper=true -Dwith-mpi=disabled -Denable-openblas=true -Denable-cuda=true --prefix=/home/areen/Programs/SU2_Install_GPU
./ninja -C build install