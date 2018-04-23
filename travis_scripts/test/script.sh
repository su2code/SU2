#!/bin/sh

# Run the tests via the Python scripts
ls
echo $PWD
travis_wait 90 python $TEST_SCRIPT


