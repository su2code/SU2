#!/bin/sh -l

echo "SU2 v7 Docker Compilation Container"
usage="$(basename "$0") [-h] [-f build_flags] [-b branch_name] 
where:
    -h  show this help text
    -f  meson build flags (ignored if a directory \"docker_build\" is found at /src/SU2).
    -b  branch name (if not given, existing SU2 directory must be mounted in /src/SU2).

Compiled binaries can be found at /install/. Mount that directory for access.
Note: If you specify a working directory using the --workdir option for docker, 
      append this directory to all paths above (e.g. use --workdir=/tmp if running in user mode)."

flags=""
branch=""
workdir=$PWD

export CCACHE_DIR=$workdir/ccache

if [ "$#" -ne 0 ]; then
  while [ "$(echo $1 | cut -c1)" = "-" ]
    do
        case "$1" in
            -f)
                    flags=$2
                    shift 2
                ;;  
            -b)
                    branch=$2
                    shift 2
                ;;
            *)
                    echo "$usage" >&2
                    exit 1
                ;;
    esac
    done
fi

if [ ! -z "$branch" ]; then
  name="SU2_$(echo $branch | sed 's/\//_/g')"
  echo "Branch provided. Cloning to $PWD/src/$name"
  if [ ! -d "src" ]; then
    mkdir "src"
  fi
  cd "src"
  git clone --recursive https://github.com/su2code/SU2 $name
  cd $name
  git config --add remote.origin.fetch '+refs/pull/*/merge:refs/remotes/origin/refs/pull/*/merge'
  git config --add remote.origin.fetch '+refs/heads/*:refs/remotes/origin/refs/heads/*'
  git fetch origin
  git checkout $branch
  git submodule update
else
  if [ ! -d "src/SU2" ]; then
    echo "SU2 source directory not found. Make sure to mount existing SU2 at directory at /src/SU2. Otherwise use -b to provide a branch."
    exit 1
  fi
  cd src/SU2
fi

if [ ! -d "docker_build" ]; then
  ./meson.py docker_build --prefix=$workdir/install $flags
else
  echo "Build Directory found. Ignoring build flags. To set new flags, remove docker_build directory."
  ./meson.py docker_build --reconfigure $flags
fi

./ninja -C docker_build install


