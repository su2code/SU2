#!/bin/sh

# Get the test cases
git clone -b develop https://github.com/su2code/TestCases.git ./TestData
cp -R ./TestData/* ./TestCases/

# Get the tutorial cases
git clone -b master https://github.com/su2code/su2code.github.io ./Tutorials

# Enter the SU2/TestCases/ directory, which is now ready to run
cd TestCases/

