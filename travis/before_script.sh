#!/bin/sh

# Get the test cases
git clone -b develop https://github.com/su2code/TestCases.git ./TestData
cp -R ./TestData/* ./TestCases/

# Get the test cases
git clone -b develop https://github.com/su2code/Tutorials.git ./Tutorials

# Enter the SU2/TestCases/ directory, which is now ready to run
cd TestCases/