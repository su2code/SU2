%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SU2 configuration file                                                 %
% PaStiX support build instructions.                                     %
% Institution: Imperial College London                                   %
% File Version 6.2.0 "Falcon"                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1 - Download
% Get PaStiX 5.2.3 from https://gforge.inria.fr/frs/?group_id=186
% Get Scotch 6.X.X from https://gforge.inria.fr/projects/scotch
%
% 2 - Build Scotch
% Extract the tarball downloaded in 1 into "externals"
% Rename the directory as "scotch"
% cd scotch/src && cp Make.inc/Makefile.inc.x86-64_pc_linux2.XXXX Makefile.inc
% (choose the XXXX that matches your compiler)
% Edit Makefile.inc and delete the cflag -DSCOTCH_PTHREAD (see why in 3-ii)
% make ptscotch
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
%  iii - Set SCOTCH_HOME as SCOTCH_HOME ?= $(PWD)/../../scotch/
%  iv  - Optionally look at the BLAS section (required by "make examples")
% make all
%
% 4 - Build SU2
% Follow the normal meson build instructions, add -Denable-pastix=true,
% this requires you to compile with MKL (-Denable-mkl=true) or OpenBLAS
% (-Denable-openblas=true) support in your call to meson.py.
% If you did not build PaStiX and Scotch in the externals folders you must
% use -Dpastix_root="some path" and -Dscotch_root="another path" to
% indicate where they are RELATIVE to the meson "build" directory.
%
