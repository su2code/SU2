################################################################################
#
# \file Makefile.am
# \brief Makefile for the SU2 Python framework
# \author M. Colonno, T. Economon, F. Palacios, T. Lukaczyk
# \version 7.0.4 "Blackbird"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2020, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

bin_SCRIPTS = \
    continuous_adjoint.py \
    compute_uncertainty.py \
    finite_differences.py \
    mesh_deformation.py \
    parallel_computation.py \
    parallel_computation_fsi.py \
    package_tests.py \
    shape_optimization.py \
    merge_solution.py \
    set_ffd_design_var.py \
    compute_polar.py \
    compute_multipoint.py \
    discrete_adjoint.py \
    direct_differentiation.py \
    fsi_computation.py \
    SU2_CFD.py


mypkgdir = $(prefix)/bin
nobase_dist_mypkg_DATA= \
    SU2/__init__.py \
    SU2/eval/design.py \
    SU2/eval/functions.py \
    SU2/eval/gradients.py \
    SU2/eval/__init__.py \
    SU2/io/config.py \
    SU2/io/config_options.py \
    SU2/io/data.py \
    SU2/io/filelock.py \
    SU2/io/redirect.py \
    SU2/io/state.py \
    SU2/io/tools.py \
    SU2/io/historyMap.py \
    SU2/io/__init__.py \
    SU2/mesh/adapt.py \
    SU2/mesh/tools.py \
    SU2/mesh/__init__.py \
    SU2/opt/project.py \
    SU2/opt/scipy_tools.py \
    SU2/opt/__init__.py \
    SU2/run/adaptation.py \
    SU2/run/adjoint.py \
    SU2/run/deform.py \
    SU2/run/direct.py \
    SU2/run/interface.py \
    SU2/run/merge.py \
    SU2/run/geometry.py \
    SU2/run/projection.py \
    SU2/run/__init__.py \
    SU2/util/bunch.py \
    SU2/util/filter_adjoint.py \
    SU2/util/lhc_unif.py \
    SU2/util/misc.py \
    SU2/util/mp_eval.py \
    SU2/util/ordered_bunch.py \
    SU2/util/ordered_dict.py \
    SU2/util/plot.py \
    SU2/util/polarSweepLib.py \
    SU2/util/switch.py \
    SU2/util/which.py \
    SU2/util/__init__.py \
    FSI/__init__.py \
    FSI/FSIInterface.py \
    FSI/PitchPlungeAirfoilStructuralTester.py \
    FSI/io/__init__.py \
    FSI/io/FSI_config.py \
    FSI/util/__init__.py \
    FSI/util/switch.py

EXTRA_DIST = $(bin_SCRIPTS)
