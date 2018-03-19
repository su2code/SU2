#!/bin/sh

# Temporarily fixes Travis CI issue with paths for Python packages
export PATH=/usr/bin:$PATH

# Install the necessary packages using apt-get with sudo
sudo apt-get update -qq
sudo apt-get install -qq build-essential libopenmpi-dev

# Install Python dependencies
# http://conda.pydata.org/docs/travis.html#the-travis-yml-file
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda install -q python=$TRAVIS_PYTHON_VERSION numpy scipy mpi4py swig

# to avoid interference with MPI
test -n $CC  && unset CC
test -n $CXX && unset CXX