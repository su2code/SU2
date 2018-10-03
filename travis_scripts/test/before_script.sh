#!/bin/sh

# Get the test cases
git clone -b develop https://github.com/su2code/TestCases.git ./TestData
cp -R ./TestData/* ./TestCases/

# Get the tutorial cases
git clone -b develop https://github.com/su2code/su2code.github.io ./Tutorials

