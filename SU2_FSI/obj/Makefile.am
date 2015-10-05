################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_FSI
# \author M. Colonno, T. Economon, F. Palacios, R. Sanchez
# \version 4.0.0 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

bin_PROGRAMS = ../bin/SU2_FSI

su2_fsi_sources = ../include/SU2_FSI.hpp \
		../src/SU2_FSI.cpp \
		../include/iteration_structure_fsi.hpp \
		../src/iteration_structure_fsi.cpp \
		../../SU2_CFD/include/definition_structure.hpp \
		../../SU2_CFD/src/definition_structure.cpp

su2_fsi_cxx_flags =
su2_fsi_ldadd =

su2_fsi_ldadd += \
    ../../SU2_CFD/src/___bin_SU2_CFD-fluid_model.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-fluid_model_pig.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-fluid_model_pvdw.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-fluid_model_ppr.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-gauss_structure.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-element_linear.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-element_structure.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-integration_structure.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-integration_time.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-iteration_structure.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_structure.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_adjoint_levelset.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_adjoint_mean.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_adjoint_tne2.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_adjoint_turbulent.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_direct_poisson.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_direct_elasticity.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_direct_heat.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_direct_mean.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_direct_tne2.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_direct_transition.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_direct_turbulent.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_direct_wave.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_machine_learning.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_linearized_mean.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_linearized_turbulent.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_fem_nonlinear_elasticity.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_fem_linear_elasticity.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_fem_elasticity.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_machine_learning_turbulent.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-numerics_template.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-output_cgns.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-output_structure.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-output_tecplot.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-output_fieldview.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-output_su2.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-output_paraview.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_structure.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_adjoint_levelset.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_adjoint_mean.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_adjoint_tne2.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_adjoint_turbulent.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_adjoint_discrete.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_direct_poisson.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_direct_elasticity.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_direct_heat.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_direct_mean.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_direct_tne2.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_direct_transition.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_direct_turbulent.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_direct_wave.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_linearized_mean.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_linearized_turbulent.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_fem_elasticity.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-solver_template.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-transfer_physics.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-transfer_structure.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-transport_model.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_structure.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_adjoint_levelset.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_adjoint_mean.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_adjoint_tne2.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_adjoint_turbulent.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_adjoint_discrete.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_direct_poisson.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_direct_elasticity.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_direct_heat.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_direct_mean.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_direct_tne2.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_direct_transition.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_direct_turbulent.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_direct_wave.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_linearized_mean.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_linearized_turbulent.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_fem_elasticity.o \
    ../../SU2_CFD/src/___bin_SU2_CFD-variable_template.o

# always link to built dependencies from ./externals
su2_fsi_cxx_flags += @su2_externals_INCLUDES@
su2_fsi_ldadd += @su2_externals_LIBS@

# if BUILD_CGNS
su2_fsi_cxx_flags += @CGNS_CXX@
su2_fsi_ldadd += @CGNS_LD@
# endif

# if BUILD_MUTATIONPP
su2_fsi_cxx_flags += @MUTATIONPP_CXX@
su2_fsi_ldadd += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
su2_fsi_cxx_flags += @JSONCPP_CXX@
su2_fsi_ldadd += @JSONCPP_LD@
# endif

# if BUILD_HDF5
su2_fsi_cxx_flags += @HDF5_CXX@
su2_fsi_ldadd += @HDF5_LD@
# endif

# if BUILD_SZIP
su2_fsi_cxx_flags += @SZIP_CXX@
su2_fsi_ldadd += @SZIP_LD@
# endif

# if BUILD_ZLIB
su2_fsi_cxx_flags += @ZLIB_CXX@
su2_fsi_ldadd += @ZLIB_LD@
# endif

if BUILD_NORMAL
___bin_SU2_FSI_SOURCES = $(su2_fsi_sources)
___bin_SU2_FSI_CXXFLAGS = ${su2_fsi_cxx_flags}
___bin_SU2_FSI_LDADD = ../../Common/lib/libSU2.a ${su2_fsi_ldadd}
endif
