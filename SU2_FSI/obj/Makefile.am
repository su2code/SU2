# SU2 Makefile.am for SU2_FSI
# SU2 v3.2.9 "eagle"
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
# Created: 		05/18/15
# Last updated: 	05/19/15
#
# Created by: 		R. Sanchez
# Updated by: 		R. Sanchez
#
######################################################################################

# AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = ../bin/SU2_FSI

___bin_SU2_FSI_SOURCES = ../include/SU2_FSI.hpp \
		../src/SU2_FSI.cpp \
	../include/iteration_structure_fsi.hpp \
		../src/iteration_structure_fsi.cpp \
    ../../SU2_CFD/include/fluid_model.hpp \
    ../../SU2_CFD/include/fluid_model.inl \
		../../SU2_CFD/include/integration_structure.hpp \
		../../SU2_CFD/include/integration_structure.inl \
		../../SU2_CFD/include/iteration_structure.hpp \
		../../SU2_CFD/include/numerics_structure.hpp \
   	../../SU2_CFD/include/numerics_machine_learning.hpp \
		../../SU2_CFD/include/numerics_machine_learning.inl \
		../../SU2_CFD/include/numerics_machine_learning_turbulent.hpp \
		../../SU2_CFD/include/numerics_structure.inl \
		../../SU2_CFD/include/output_structure.hpp \
		../../SU2_CFD/include/solver_structure.hpp \
		../../SU2_CFD/include/solver_structure.inl \
    ../../SU2_CFD/include/transport_model.hpp \
    ../../SU2_CFD/include/transport_model.inl \
		../../SU2_CFD/include/variable_structure.hpp \
		../../SU2_CFD/include/variable_structure.inl \
		../../SU2_CFD/src/definition_structure.cpp \
		../../SU2_CFD/src/fluid_model.cpp \
    ../../SU2_CFD/src/fluid_model_pig.cpp \
    ../../SU2_CFD/src/fluid_model_pvdw.cpp \
    ../../SU2_CFD/src/fluid_model_ppr.cpp \
		../../SU2_CFD/src/integration_structure.cpp \
		../../SU2_CFD/src/integration_time.cpp \
		../../SU2_CFD/src/iteration_structure.cpp \
		../../SU2_CFD/src/numerics_adjoint_levelset.cpp \
		../../SU2_CFD/src/numerics_adjoint_mean.cpp \
		../../SU2_CFD/src/numerics_adjoint_tne2.cpp \
		../../SU2_CFD/src/numerics_adjoint_turbulent.cpp \
		../../SU2_CFD/src/numerics_direct_poisson.cpp \
		../../SU2_CFD/src/numerics_direct_elasticity.cpp \
		../../SU2_CFD/src/numerics_direct_heat.cpp \
		../../SU2_CFD/src/numerics_direct_mean.cpp \
		../../SU2_CFD/src/numerics_direct_tne2.cpp \
		../../SU2_CFD/src/numerics_direct_transition.cpp \
		../../SU2_CFD/src/numerics_direct_turbulent.cpp \
		../../SU2_CFD/src/numerics_direct_wave.cpp \
		../../SU2_CFD/src/numerics_machine_learning.cpp \
		../../SU2_CFD/src/numerics_linearized_mean.cpp \
		../../SU2_CFD/src/numerics_linearized_turbulent.cpp \
		../../SU2_CFD/src/numerics_structure.cpp \
		../../SU2_CFD/src/numerics_machine_learning_turbulent.cpp \
		../../SU2_CFD/src/numerics_template.cpp \
		../../SU2_CFD/src/output_cgns.cpp \
		../../SU2_CFD/src/output_structure.cpp \
		../../SU2_CFD/src/output_tecplot.cpp \
    ../../SU2_CFD/src/output_fieldview.cpp \
    ../../SU2_CFD/src/output_su2.cpp \
    ../../SU2_CFD/src/output_paraview.cpp \
		../../SU2_CFD/src/solver_adjoint_levelset.cpp \
		../../SU2_CFD/src/solver_adjoint_mean.cpp \
		../../SU2_CFD/src/solver_adjoint_tne2.cpp \
		../../SU2_CFD/src/solver_adjoint_turbulent.cpp \
		../../SU2_CFD/src/solver_direct_poisson.cpp \
		../../SU2_CFD/src/solver_direct_elasticity.cpp \
		../../SU2_CFD/src/solver_direct_heat.cpp \
		../../SU2_CFD/src/solver_direct_mean.cpp \
		../../SU2_CFD/src/solver_direct_tne2.cpp \
		../../SU2_CFD/src/solver_direct_transition.cpp \
		../../SU2_CFD/src/solver_direct_turbulent.cpp \
		../../SU2_CFD/src/solver_direct_wave.cpp \
		../../SU2_CFD/src/solver_linearized_mean.cpp \
		../../SU2_CFD/src/solver_linearized_turbulent.cpp \
		../../SU2_CFD/src/solver_structure.cpp \
		../../SU2_CFD/src/solver_template.cpp \
    ../../SU2_CFD/src/transport_model.cpp \
		../../SU2_CFD/src/variable_adjoint_levelset.cpp \
		../../SU2_CFD/src/variable_adjoint_mean.cpp \
		../../SU2_CFD/src/variable_adjoint_tne2.cpp \
		../../SU2_CFD/src/variable_adjoint_turbulent.cpp \
		../../SU2_CFD/src/variable_direct_poisson.cpp \
		../../SU2_CFD/src/variable_direct_elasticity.cpp \
		../../SU2_CFD/src/variable_direct_heat.cpp \
		../../SU2_CFD/src/variable_direct_mean.cpp \
		../../SU2_CFD/src/variable_direct_tne2.cpp \
		../../SU2_CFD/src/variable_direct_transition.cpp \
		../../SU2_CFD/src/variable_direct_turbulent.cpp \
		../../SU2_CFD/src/variable_direct_wave.cpp \
		../../SU2_CFD/src/variable_linearized_mean.cpp \
		../../SU2_CFD/src/variable_linearized_turbulent.cpp \
		../../SU2_CFD/src/variable_structure.cpp \
		../../SU2_CFD/src/variable_template.cpp

___bin_SU2_FSI_CXXFLAGS =
___bin_SU2_FSI_LDADD = ../../Common/lib/libSU2.a

# always link to built dependencies from ./externals
___bin_SU2_FSI_CXXFLAGS += @su2_externals_INCLUDES@
___bin_SU2_FSI_LDADD += @su2_externals_LIBS@

# if BUILD_CGNS
___bin_SU2_FSI_CXXFLAGS += @CGNS_CXX@
___bin_SU2_FSI_LDADD += @CGNS_LD@
# endif

# if BUILD_MUTATIONPP
___bin_SU2_FSI_CXXFLAGS += @MUTATIONPP_CXX@
___bin_SU2_FSI_LDADD += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
___bin_SU2_FSI_CXXFLAGS += @JSONCPP_CXX@
___bin_SU2_FSI_LDADD += @JSONCPP_LD@
# endif
