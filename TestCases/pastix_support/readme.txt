%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SU2 configuration file                                                 %
% PaStiX support build instructions.                                     %
% Institution: Imperial College London                                   %
% File Version 8.0.0 "Harrier"                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1 - Download
% Get PaStiX 5.2.3 from https://gforge.inria.fr/frs/?group_id=186
% Get Scotch 6.0.6 from https://gforge.inria.fr/projects/scotch
% Note: These two versions were tested on a number of platforms, some
% issues were encountered with more recent version of Scotch, and PaStiX
% 6.0.X is not compatible with SU2 as it does not support MPI yet.
%
% 2 - Build Scotch
% Extract the tarball downloaded in 1 into "externals"
% Rename the directory as "scotch"
% cd scotch/src && cp Make.inc/Makefile.inc.x86-64_pc_linux2.XXXX Makefile.inc
% (choose the XXXX that matches your compiler)
% Edit Makefile.inc and delete the cflag -DSCOTCH_PTHREAD (see why in 3-ii)
% "make ptscotch"
%
% Note: If you build SU2 (7.0.1+) with OpenMP support (-Dwith-omp=true),
%       AND your system supports MPI_THREAD_MULTIPLE, you do not need to
%       delete the -DSCOTCH_PTHREAD flag (but doing so does no harm).
%
% 3 - Build PaStiX
% Extract the tarball downloaded in 1 into "externals"
% Rename the directory as "pastix"
% cd pastix/src && cp config/LINUX-XXXX.in config.in
% (again choose the XXXX that matches your compiler)
% Edit config.in
%  i   - Uncomment the lines for "VERSIONINT  = _int32"
%  ii  - Uncomment the lines for "VERSIONSMP  = _nosmp",
%        SU2 does not currently support MPI+Threads.
%  iii - Set SCOTCH_HOME as SCOTCH_HOME ?= ${PWD}/../../scotch/
%  iv  - Comment out the lines for "Hardware Locality", this may only be
%        important for an MPI+Threads build (usually it is not).
%  v   - Optionally look at the BLAS section (only required by "make examples")
% "make all"
%
% Note: If you build SU2 (7.0.1+) with OpenMP support (-Dwith-omp=true),
%       skip 3-ii, note however that this may not work well with SU2_CFD_AD.
%       If you do use MPI+Threads, it is important for good performance that your
%       system supports MPI_THREAD_MULTIPLE (SU2_CFD --thread_multiple ...)
%       Furthermore, if MPI_THREAD_MULTIPLE is NOT supported, you need to
%       uncomment the line with "-DPASTIX_FUNNELED" in config.in.
%       Finally, if you just use threads (no MPI) this is not important.
%
% 4 - Build SU2
% Follow the normal meson build instructions, add -Denable-pastix=true,
% this requires you to compile with MKL (-Denable-mkl=true) or OpenBLAS
% (-Denable-openblas=true) support in your call to meson.py.
% If you did not build PaStiX and Scotch in the externals folders you must
% use -Dpastix_root="some path" and -Dscotch_root="another path" to
% indicate where they are RELATIVE to the SU2 directory.
% You need sequential versions of BLAS. But when using MPI+Threads beware that
% OpenBLAS needs to have parallel support otherwise the solver will crash, if
% you get poor performance see 5.4 below.
%
% 5 - Common problems and known issues
% - OpenMPI 4 does not work with PaStiX 5, downgrade to 3.1.4.
% - Very early versions of OpenMPI 3 may have problems with MPI types.
% - OpenBLAS build fails when linking executables. Old versions (e.g.
%   0.2.18) did not provide LAPACK support, rebuild or upgrade.
% - Very bad performance with OpenBLAS on some systems (observed on Ubuntu
%   16.04) try "export OMP_NUM_THREADS=1" before running SU2, check that
%   you only see N SU2 processes running at 100% (mpirun -n N SU2_XXX).
% - Cannot find BLAS dependency:
%   i   - On some OS the package has a different name (e.g. Ubuntu 16.04
%     blas-openblas instead of openblas), use -Dblas-name="right name" in
%     call to meson.py
%   ii  - The name is right but meson cannot find it. Set env variable
%     PKG_CONFIG_PATH=$PKG_CONFIG_PATH:"directory with someblas.pc file"
% - MKL is not in its standard directory (/opt/intel/mkl), use option
%   -Dmkl_root="non standard directory" in call to meson.py (headers are
%   expected in "include" and libraries in "lib/intel64").
%
% 6 - Tested platforms
% - Ubuntu 18.04, gcc 7.4, ompi 3.1.4, mkl 2017, openblas 0.2.20 and 0.3.2.dev
% - Ubuntu 16.04, gcc 5.4, ompi 3.1.4, mkl 2017 and 2019
% - CentOS 7.6.1810, gcc 5.4, ompi 3.1.4, mkl 2017
% - CentOS 7.6.1810, gcc 5.4, impi 2018, mkl 2019
% - CentOS 7.6.1810, gcc 8.2, impi 2018, mkl 2019
%
