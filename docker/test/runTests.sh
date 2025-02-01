#!/bin/sh -l
# List of arguments to this script
# $1 : Tag/SHA of su2code/SU2
# $2 : Tag/SHA of su2code/TestCases
# $3 : Tag/SHA of su2code/Tutorials
# $4 : Path of the installation directory
# $5 : Test script to execute

usage="$(basename "$0") [-h] [-t tutorial_branch] [-b su2_branch] [-c testcases_branch] [-s test_script]
where:
    -h  show this help text
    -t  branch of su2code/Tutorials repo
    (if not provided, it is assumed that it is mounted at /src/SU2)
    -b  branch of su2code/SU2 repo.
    (if not provided, it is assumed that it is mounted at /src/Tutorials)
    -c  branch of su2code/TestCases repo.
    (if not provided, it is assumed that it is mounted at /src/TestData)
    -s  name of the test script to execute (default: parallel_regression.py).

Compiled binaries must be mounted at /install/ !

Note: If you specify a working directory using the --workdir option for docker,
      append this directory to all paths above (e.g. use --workdir=/tmp if running in user mode)."

su2branch=""
testbranch=""
tutorialbranch=""
script="parallel_regression.py"

while [ "`echo $1 | cut -c1`" = "-" ]
do
    case "$1" in
        -t)
                tutorialbranch=$2
                shift 2
            ;;
        -b)
                su2branch=$2
                shift 2
            ;;
        -c)
                testbranch=$2
                shift 2
            ;;
        -s)
                script=$2
                shift 2
            ;;
        *)
                echo "$usage" >&2
                exit 1
            ;;
esac
done


if [ ! -d "tests" ]; then
  mkdir tests
fi
if [ ! -d "src" ]; then
  mkdir "src"
fi

if [ ! -z "$su2branch" ]; then
  name="SU2_$(echo $su2branch | sed 's/\//_/g')"
  echo "Branch provided. Cloning to $PWD/src/$name"
  cd "src"

  # Clone su2code/SU2, su2code/TestCases and su2code/Tutorials
  git clone -b master https://github.com/su2code/SU2 $name
  cd $name
  git config --add remote.origin.fetch '+refs/pull/*/merge:refs/remotes/origin/refs/pull/*/merge'
  git config --add remote.origin.fetch '+refs/heads/*:refs/remotes/origin/refs/heads/*'
  git fetch origin
  git checkout $su2branch
  git submodule update
  cd ..
  cd ..
  cp -r src/$name/TestCases tests/.
else
  if [ ! -d "src/SU2" ]; then
    echo "SU2 source directory not found. Make sure to mount existing SU2 at directory at /src/SU2. Otherwise use -b to provide a branch."
    exit 1
  fi
  cp -r src/SU2/TestCases tests/.
fi
if [ ! -z "$testbranch" ]; then
  git clone --depth=1 -b $testbranch https://github.com/su2code/TestCases.git ./TestData
else
  if [ ! -d "src/TestData" ]; then
    echo "$PWD/src/TestData not found. Make sure to mount existing su2code/TestCases repo or use -c to provide a branch to clone."
    exit 1
  fi
fi
cp -R ./TestData/* tests/TestCases/
if [ ! -z "$tutorialbranch" ]; then
  git clone --depth=1 -b $tutorialbranch https://github.com/su2code/su2code.github.io ./Tutorials
else
  if [ ! -d "src/Tutorials" ]; then
    echo "$PWD/src/Tutorials not found. Make sure to mount existing su2code/su2code.github.io repo or use -t to provide a branch to clone."
    exit 1
  fi
fi
cp -R ./Tutorials/ tests/.

# Set the environment variables
export SU2_RUN=$PWD/install/bin
export PATH=$SU2_RUN:$PATH
export PYTHONPATH=$SU2_RUN:$PYTHONPATH
export OMPI_MCA_btl_vader_single_copy_mechanism=none
export SU2_MPI_COMMAND='mpirun --allow-run-as-root -n %i %s'
alias mpirun='mpirun --allow-run-as-root'

# Run Test Script
cd tests/TestCases
python $script

if [ $? -eq 0 ]; then
    echo "Tests passed"
    exit 0
else
    echo "Tests failed"
    exit 1
fi
