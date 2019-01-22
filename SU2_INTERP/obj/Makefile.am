################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_SOL
# \author M. Colonno, T. Economon, F. Palacios
# \version 6.1.0 "Falcon"
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
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = ../bin/SU2_INTERP

___bin_SU2_INTERP_SOURCES = \
  ../include/SU2_INTERP.hpp \
  ../src/SU2_INTERP.cpp

___bin_SU2_INTERP_CXXFLAGS =
___bin_SU2_INTERP_LDADD = \
  ../../SU2_CFD/src/libSU2Core_a-data_manufactured_solutions.o \
  ../../SU2_CFD/src/libSU2Core_a-definition_structure.o \
  ../../SU2_CFD/src/libSU2Core_a-fluid_model.o \
  ../../SU2_CFD/src/libSU2Core_a-fluid_model_pig.o \
  ../../SU2_CFD/src/libSU2Core_a-fluid_model_pvdw.o \
  ../../SU2_CFD/src/libSU2Core_a-fluid_model_ppr.o \
  ../../SU2_CFD/src/libSU2Core_a-fluid_model_inc.o \
  ../../SU2_CFD/src/libSU2Core_a-integration_structure.o \
  ../../SU2_CFD/src/libSU2Core_a-integration_time.o \
  ../../SU2_CFD/src/libSU2Core_a-driver_direct_multizone.o \
  ../../SU2_CFD/src/libSU2Core_a-driver_direct_singlezone.o \
  ../../SU2_CFD/src/libSU2Core_a-driver_structure.o \
  ../../SU2_CFD/src/libSU2Core_a-iteration_structure.o \
  ../../SU2_CFD/src/libSU2Core_a-numerics_adjoint_mean.o \
  ../../SU2_CFD/src/libSU2Core_a-numerics_adjoint_turbulent.o \
  ../../SU2_CFD/src/libSU2Core_a-numerics_direct_heat.o \
  ../../SU2_CFD/src/libSU2Core_a-numerics_direct_mean.o \
  ../../SU2_CFD/src/libSU2Core_a-numerics_direct_mean_inc.o \
  ../../SU2_CFD/src/libSU2Core_a-numerics_direct_transition.o \
  ../../SU2_CFD/src/libSU2Core_a-numerics_direct_turbulent.o \
  ../../SU2_CFD/src/libSU2Core_a-numerics_direct_elasticity_nonlinear.o \
  ../../SU2_CFD/src/libSU2Core_a-numerics_direct_elasticity_linear.o \
  ../../SU2_CFD/src/libSU2Core_a-numerics_direct_elasticity.o \
  ../../SU2_CFD/src/libSU2Core_a-numerics_structure.o \
  ../../SU2_CFD/src/libSU2Core_a-numerics_template.o \
  ../../SU2_CFD/src/libSU2Core_a-output_cgns.o \
  ../../SU2_CFD/src/libSU2Core_a-output_structure.o \
  ../../SU2_CFD/src/libSU2Core_a-output_physics.o \
  ../../SU2_CFD/src/libSU2Core_a-output_tecplot.o \
  ../../SU2_CFD/src/libSU2Core_a-output_fieldview.o \
  ../../SU2_CFD/src/libSU2Core_a-output_su2.o \
  ../../SU2_CFD/src/libSU2Core_a-output_paraview.o \
  ../../SU2_CFD/src/libSU2Core_a-python_wrapper_structure.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_adjoint_mean.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_adjoint_turbulent.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_adjoint_discrete.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_adjoint_elasticity.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_direct_heat.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_direct_mean.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_direct_mean_fem.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_direct_mean_inc.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_direct_transition.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_direct_turbulent.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_direct_elasticity.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_structure.o \
  ../../SU2_CFD/src/libSU2Core_a-solver_template.o \
  ../../SU2_CFD/src/libSU2Core_a-transfer_physics.o \
  ../../SU2_CFD/src/libSU2Core_a-transfer_structure.o \
  ../../SU2_CFD/src/libSU2Core_a-transport_model.o \
  ../../SU2_CFD/src/libSU2Core_a-variable_adjoint_mean.o \
  ../../SU2_CFD/src/libSU2Core_a-variable_adjoint_turbulent.o \
  ../../SU2_CFD/src/libSU2Core_a-variable_adjoint_discrete.o \
  ../../SU2_CFD/src/libSU2Core_a-variable_adjoint_elasticity.o \
  ../../SU2_CFD/src/libSU2Core_a-variable_direct_heat.o \
  ../../SU2_CFD/src/libSU2Core_a-variable_direct_mean.o \
  ../../SU2_CFD/src/libSU2Core_a-variable_direct_mean_inc.o \
  ../../SU2_CFD/src/libSU2Core_a-variable_direct_transition.o \
  ../../SU2_CFD/src/libSU2Core_a-variable_direct_turbulent.o \
  ../../SU2_CFD/src/libSU2Core_a-variable_direct_elasticity.o \
  ../../SU2_CFD/src/libSU2Core_a-variable_structure.o \
  ../../SU2_CFD/src/libSU2Core_a-variable_template.o \
  ../../Common/lib/libSU2.a

# always link to built dependencies from ./externals
___bin_SU2_INTERP_CXXFLAGS += @su2_externals_INCLUDES@
___bin_SU2_INTERP_LDADD += @su2_externals_LIBS@
___bin_SU2_INTERP_LDADD += @su2_externals_LIBPTHREAD@

# if BUILD_MUTATIONPP
___bin_SU2_INTERP_CXXFLAGS += @MUTATIONPP_CXX@
___bin_SU2_INTERP_LDADD += @MUTATIONPP_LD@
# endif
