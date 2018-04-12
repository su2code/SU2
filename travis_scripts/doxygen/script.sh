#!/bin/sh

export GH_REPO_NAME=Doxygen
export DOXYFILE=$TRAVIS_BUILD_DIR/doc/Doxyfile
export GH_REPO_REF=github.com/su2code/Doxygen.git
sh ../../doc/generateDocumentationAndDeploy.sh

