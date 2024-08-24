#!bin/bash

rm -r build

export PKG_CONFIG_PATH="/home/areen/Programs/OpenBLAS_install/lib:$PKG_CONFIG_PATH"
./meson.py build --buildtype=debug -Denable-pywrapper=true -Dwith-mpi=disabled -Denable-openblas=true  --prefix=/home/areen/Programs/SU2_Install_CPU
./ninja -C build install