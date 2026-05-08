%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SU2 configuration file                                                 %
% PaStiX support build instructions.                                     %
% Institution: Imperial College London                                   %
% File Version 8.5.0 "Harrier"                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1 - Build and install dependencies
%  - A BLAS library with LAPACK and LAPACKE support (e.g. OpenBLAS),
%    sequential versions are required (threading is managed by the solver).
%  - Scotch and PT-Scotch.
%  - HWLOC (recommended for PaStiX).
% If you use your OS package manager for these, the PaStiX and SU2 build
% systems should be able to find them automatically.
% If you build them from source and install in a custom location, you will
% need to use CPATH, LIBRARY_PATH, and LD_LIBRARY_PATH for the build
% systems to find them.
%
% 2 - Build PaStiX
% Get PaStiX 6.X.X from https://gitlab.inria.fr/solverstack/pastix.git
% Follow their build instructions, make sure to compile it with support
% for MPI and PT-Scotch (see the cmake flags).
% You may need to adjust the integer type from 64 to 32 bits depending
% on how your version of Scotch was built.
%
% 4 - Build SU2
% Follow the normal meson build instructions, add -Denable-pastix=true,
% this requires you to compile with MKL (-Denable-mkl=true) or OpenBLAS
% (-Denable-openblas=true) support in your call to meson.py.
% The reliable way for meson to find PaStiX is use PKG_CONFIG_PATH to
% make the pastix.pc.in file discoverable. For the other dependencies,
% see #1. You should build with OpenMP to get the best performance,
% using only MPI with PaStiX is not efficient.
%
% 5 - Running SU2
% When running with MPI + OpenMP you need to pay attention to core and
% thread binding to get the best performance. For example, to use 2 ranks
% with 8 threads each (effectively) you may need something like:
% mpirun -n 2 --map-by NUMA:PE=8 SU2_CFD -t 8 --thread_multiple config.cfg
% otherwise the MPI implementation may not let each rank use the requested
% threads, or PaStiX/HWLOC may bind the threads of different ranks to the
% same core. Different MPI versions may have other arguments.
%
% 6 - Tested platforms
% - Ubuntu 24.04 (gcc 13), ompi 5.0.6, default OpenBLAS, Scotch, PT-Scotch,
%   and HWLOC (i.e. from apt install ...).
%
