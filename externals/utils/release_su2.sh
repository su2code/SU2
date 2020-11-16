#!/bin/bash

# Shell script to create SU2 binaries on Linux/Mac

# First check for input args.
if [ $# != 4 ]; then
    printf "\n\n Incorrect input args. Must specify tag, OS, compiler, and build cores."
    printf "\n Example: '$ ./release_su2.sh 6.0.0 macos10.13 llvm 4'"
    printf "\n Now exiting...\n\n"
    exit 1
fi

# Set up some directory variables
CURDIR=$(pwd)
BUILDDIR=${CURDIR}/tmp_su2
SU2DIR=${BUILDDIR}/SU2
RELEASEDIR=${BUILDDIR}/su2-$1
PKGNAME=su2-$1-$2-$3.tgz

printf "\n\n -- Creating SU2 v$1 Binaries --"

# Create a temp directory and download the code (master branch)
printf "\n\n Downloading the current master branch...\n\n"
if [ ! -d ${BUILDDIR} ]
then
 mkdir ${BUILDDIR}
fi
cd ${BUILDDIR}
git clone https://github.com/su2code/SU2.git

# Check out the requested release tag
printf "\n\n Checking out the requested tag for v$1...\n\n"
cd ${SU2DIR}
git checkout tags/v$1

# Build the vanilla version of SU2
printf "\n\n Configuring and building the vanilla version...\n\n"
./bootstrap
./configure --prefix=${SU2DIR} CXXFLAGS='-O3' LDFLAGS='-static-libstdc++'
make clean
make -j$4 install
cd ${BUILDDIR}

# Create an appropriately named folder for the release
printf "\n\n Packing up the binaries and quick start..."
if [ ! -d ${RELEASEDIR} ]
then
 mkdir ${RELEASEDIR}
 mkdir ${RELEASEDIR}/bin
 mkdir ${RELEASEDIR}/QuickStart
fi
cp -R ${SU2DIR}/bin/* ${RELEASEDIR}/bin/
cp ${SU2DIR}/QuickStart/inv_NACA0012.cfg ${RELEASEDIR}/QuickStart/
cp ${SU2DIR}/QuickStart/mesh_NACA0012_inv.su2 ${RELEASEDIR}/QuickStart/
rm -fr ${SU2DIR}

# Tar up the finished build directory and clean up
printf "\n\n Cleaning up..."
tar -czf ${PKGNAME} -C ${RELEASEDIR}/.. .
mv ${PKGNAME} ${CURDIR}
cd ${CURDIR}
rm -fr ${BUILDDIR}
printf "\n\n Finished build. Generated su2-$1-$2-$3.tgz for distribution.\n\n"
