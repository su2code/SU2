#!/bin/sh

# Configure, make, and install SU2
- echo $TRAVIS_BUILD_DIR
- echo $CONFIGURE_COMMAND
#- $CONFIGURE_COMMAND
make -j 4
make install

# Add environmental variables according to the configure step
export SU2_RUN=$TRAVIS_BUILD_DIR/bin
export SU2_HOME=$TRAVIS_BUILD_DIR
export PATH=$PATH:$SU2_RUN
export PYTHONPATH=$PYTHONPATH:$SU2_RUN

