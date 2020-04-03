<p align="center">
<img width="250" height="154" src="Common/doc/logoSU2small.png">
</p>


# SU2 (ver. 7.0.1 "Blackbird"): The Open-Source CFD Code

Computational analysis tools have revolutionized the way we design engineering systems, but most established codes are proprietary, unavailable, or prohibitively expensive for many users. The SU2 team is changing this, making multiphysics analysis and design optimization freely available as open-source software and involving everyone in its creation and development. 

For an overview of the technical details in SU2, please see the following AIAA Journal article:

"SU2: An open-source suite for multiphysics simulation and design," AIAA Journal, 54(3):828-846, 2016. http://arc.aiaa.org/doi/10.2514/1.J053813

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

Continuous Integration:<br/>
[![Regression Testing](https://github.com/su2code/SU2/workflows/Regression%20Testing/badge.svg?branch=develop)](https://github.com/su2code/SU2/actions)
[![Release](https://github.com/su2code/SU2/workflows/Release%20Management/badge.svg?branch=develop)](https://github.com/su2code/SU2/actions)

# SU2 Introduction

SU2 is a suite of open-source software tools written in C++ for the numerical solution of partial differential equations (PDE) and performing PDE constrained optimization. 

The primary applications are computational fluid dynamics and aerodynamic shape optimization, but has been extended to treat more general equations such as electrodynamics and chemically reacting flows. 

You will find more information and the latest news in:
   - SU2 Home Page: https://su2code.github.io
   - GitHub repository: https://github.com/su2code
   - CFD Online: http://www.cfd-online.com/Forums/su2/
   - Twitter: https://twitter.com/su2code
   - Facebook: https://www.facebook.com/su2code


# SU2 Installation

## Precompiled binaries for Linux, MacOS, Windows

You can find precompiled binaries of the latest version on our [download page](https://su2code.github.io/download/) or under [releases](https://github.com/su2code/SU2/releases).

## Build SU2
The build system of SU2 is based on a combination of [meson](http://mesonbuild.com/) (as the front-end) and [ninja](https://ninja-build.org/) (as the back-end). Meson is an open source build system meant to be both extremely fast, and, even more importantly, as user friendly as possible. Ninja is a small low-level build system with a focus on speed. 

Short summary of the minimal requirements:

- C/C++ compiler
- Python 3

**Note:** all other necessary build tools and dependencies are shipped with the source code or are downloaded automatically.

If you have these tools installed, you can create a configuration using the `meson.py` found in the root source code folder:
```
./meson.py build
```
Use `ninja` to compile and install the code

```
./ninja -C build install
```

For more information on how to install and build SU2 on Linux, MacOS or Windows, have a look at the [documentation](https://su2code.github.io/docs_v7/).

##  SU2 Path setup

When installation is complete, please be sure to add the `$SU2_HOME` and `$SU2_RUN` environment variables, and update your `$PATH` with `$SU2_RUN`. 

For example, add these lines to your `.bashrc` file:
```
export SU2_RUN="your_prefix/bin"
export SU2_HOME="/path/to/SU2vX.X.X/"
export PATH=$PATH:$SU2_RUN
export PYTHONPATH=$SU2_RUN:$PYTHONPATH
```

`$SU2_RUN` should point to the folder where all binaries and python scripts were installed. This is the prefix you set with the --prefix option to meson. Note that the bin/ directory is automatically added to your prefix path.

`$SU2_HOME` should point to the root directory of the source code distribution, i.e., `/path/to/SU2vX.X.X/`.

Thanks for building, and happy optimizing!

- The SU2 Development Team


# SU2 Developers


We follow the popular "GitFlow" branching model for scalable development. In the SU2 repository, the master branch represents the latest stable major or minor release (7.0, 6.2.0, etc.), it should only be modified during version releases. Work that is staged for release is put into the develop branch via Pull Requests on GitHub from various "feature" branches where folks do their day-to-day work on the code. At release time, the work that has been merged into the develop branch is pushed to the master branch and tagged as a release.

SU2 is being developed by individuals and organized teams all around the world. 

A list of current contributors can be found in the AUTHORS.md file.


# Instructions to configure SU2 for implicitDG discretization

If you are on sherlock, load the modules necessary. Last option superlu module is optional, but I would recommend u to configure that since it will be easier to run comparison studies later on.
```
module load gcc eigen metis parmetis openmpi superlu
```
Configuring SU2. SU2 offers the conventional Makefile build and meson/ninja build. Personally I’ve been using the old fashion makefile build as linking/unlinking library on sherlock using meson gave me a lot of headache.

Initiating configurtion
```
./bootstrap
```

Installing necessary 3rd party libraries, in particular, CodiPack (automatic differentiation tool) and MediPack (MPI wrapper for CoDiPack)
```
./preconfigure.py --enable-direct-diff --enable-mpi
```

Configure SU2 setup. A couple of notes:
1. SU2 requires metis and parmetis and it actually comes with the two libraries by default, which means you wouldnt have to load metis/parmetis on sherlock. However, the Superlu solver on sherlock loads metis/parmetis modules on sherlock, which then gives a version conflict error if the SU2 metis/parmetis are present. Hence, my solution at the moment is to disable the installation of metis/parmetis in SU2, but using metis/parmetis on sherlock server. Two flags `-DHAVE_METIS -DHAVE_PARMETIS` and two links '-lmetis -lparmetis' are used to link to the sherlock metis/parmetis libraries. I hope this doesnt give you any problem when you use your spaND solver.
2. You are welcome to change the optimization flag if you start to run performance tests.
3. Please install SU2 into a SCRATCH folder. Installing it into a HOME directory will easily lead to memory issues when you start to run it locally for quick tests. SCRATCH directory seems to be less stringent and allow you to run some quick tests without the need to submit jobs to the queue.
```
./configure --prefix=$SCRATCH/<your directory> --enable-mpi --enable-codi-forward --disable-metis --disable-parmetis CXXFLAGS="-O0 -g -DHAVE_METIS -DHAVE_PARMETIS -lmetis -lparmetis -lsuperlu_dist"
```

Once you are done with the above commend, the screen output towards the end will have 4 lines of path variables that you need to copy and add to the bash file. Typically they look something like this (this is my version, yours will be different depending on the installing path you specified earlier with the prefix option)
```
export SU2_RUN="/scratch/users/zanxu/SU2/bin"
export SU2_HOME="/scratch/users/zanxu/SU2"
export PATH=$PATH:$SU2_RUN
export PYTHONPATH=$PYTHONPATH:$SU2_RUN
```
(Later on you will need to copy and paste these 4 lines and add to the end of the bash file. You will need to reload the terminal or refresh the bashfile to run SU2. but this can be done last.)
```
cd
vi .bashrc
(copy and paste)
```

In the original SU2 directory, start installing the objects and executables. On sherlock, I would highly recommend you to request multiple N nodes via `-j N` to parallelize the installation for speedup. It takes quite some time. If the parallelization fails at some point, just rerun it. I have experienced some minor problems before they usually can be fixed with a re-run.
```
make install -j 8
```

This should be the end of installation. Once you are done with installation and reloading/refreshing bash file, you should be able to run a quick test. you can go to QuickStart directory and run SU2 to test its completion.
```
cd QuickStart
```
Usually the SU2 executable is `SU2_CFD`, but for our case, we use automatic differentiation for computing Jacobian, so our SU2 is built with automatic differentiation support and the executable is `SU2_CFD_DIRECTDIFF`. The argument required is a single `.cfg` configuration file. You shouldnt worry about the configuration file. Mostly they are used to specify the flow solver, not really the linear solver. I will make sure you have the right config file when we start to run tests together.
To run the QuickStart test case sequentially
```
SU2_CFD_DIRECTDIFF inv_NACA0012.cfg
```
To run in parallel
```
mpirun -n 8 SU2_CFD_DIRECTDIFF inv_NACA0012.cfg
```