################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_MET
# \author M. Colonno, T. Economon, F. Palacios
# \version 6.1.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.$(OBJEXT)rg>
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
# License along with SU2. If not, see <http://www.gnu.$(OBJEXT)rg/licenses/>.
#
################################################################################

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = 

___bin_SU2_MET_SOURCES = \
  ../include/SU2_MET.hpp \
  ../src/SU2_MET.cpp
___bin_SU2_MET_CXXFLAGS =
___bin_SU2_MET_LDADD = 

if BUILD_REVERSE

bin_PROGRAMS += ../bin/SU2_MET

___bin_SU2_MET_LDADD += \
  ../../SU2_CFD/src/libSU2Core_AD_a-data_manufactured_solutions.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-definition_structure.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-error_estimation_structure.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-fluid_model.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-fluid_model_pig.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-fluid_model_pvdw.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-fluid_model_ppr.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-fluid_model_inc.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-integration_structure.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-integration_time.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-driver_direct_multizone.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-driver_direct_singlezone.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-driver_structure.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-iteration_structure.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-numerics_adjoint_mean.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-numerics_adjoint_turbulent.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-numerics_direct_heat.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-numerics_direct_mean.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-numerics_direct_mean_inc.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-numerics_direct_transition.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-numerics_direct_turbulent.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-numerics_direct_elasticity_nonlinear.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-numerics_direct_elasticity_linear.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-numerics_direct_elasticity.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-numerics_structure.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-numerics_template.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-output_cgns.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-output_structure.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-output_physics.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-output_tecplot.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-output_fieldview.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-output_inria.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-output_su2.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-output_paraview.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-python_wrapper_structure.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_adjoint_mean.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_adjoint_turbulent.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_adjoint_discrete.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_adjoint_elasticity.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_direct_heat.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_direct_mean.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_direct_mean_fem.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_direct_mean_inc.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_direct_transition.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_direct_turbulent.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_direct_elasticity.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_structure.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-solver_template.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-transfer_physics.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-transfer_structure.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-transport_model.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-variable_adjoint_mean.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-variable_adjoint_turbulent.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-variable_adjoint_discrete.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-variable_adjoint_elasticity.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-variable_direct_heat.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-variable_direct_mean.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-variable_direct_mean_inc.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-variable_direct_transition.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-variable_direct_turbulent.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-variable_direct_elasticity.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-variable_structure.$(OBJEXT) \
  ../../SU2_CFD/src/libSU2Core_AD_a-variable_template.$(OBJEXT) \
  ../../Common/lib/libSU2_AD.a

# always link to built dependencies from ./externals
___bin_SU2_MET_CXXFLAGS += @su2_externals_INCLUDES@
___bin_SU2_MET_LDADD += @su2_externals_LIBS@
___bin_SU2_MET_LDADD += @su2_externals_LIBPTHREAD@

# add reverse flags
___bin_SU2_MET_CXXFLAGS += @REVERSE_CXX@
___bin_SU2_MET_LDADD += @REVERSE_LIBS@ 

# if BUILD_MUTATIONPP
___bin_SU2_MET_CXXFLAGS += @MUTATIONPP_CXX@
___bin_SU2_MET_LDADD += @MUTATIONPP_LD@
# endif

endif

