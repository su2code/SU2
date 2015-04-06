# SU2 Makefile.am for FSI
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated: 	3/11/13
#
# Created by: 		A. Variyar
# Updated by: 		A. Variyar
#
##########################

# AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = ../bin/pySU2

___bin_pySU2_SOURCES = ../include/SU2_FSI.hpp \
		../include/SU2_interface.hpp \
../src/SU2_FSI.cpp \
../src/SU2_interface.cpp \
../../SU2_CFD/include/definition_structure.hpp \
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
../../SU2_CFD/src/output_paraview.cpp \
../../SU2_CFD/src/output_fieldview.cpp \
../../SU2_CFD/src/output_su2.cpp \
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


___bin_pySU2_CXXFLAGS =
___bin_pySU2_LDADD = ../../Common/lib/libSU2.a

# always link to built dependencies from ./externals
___bin_pySU2_CXXFLAGS += @su2_externals_INCLUDES@
___bin_pySU2_LDADD += @su2_externals_LIBS@

# if BUILD_CGNS
___bin_pySU2_CXXFLAGS += @CGNS_CXX@
___bin_pySU2_LDADD += @CGNS_LD@
# endif

# if BUILD_MUTATIONPP
___bin_pySU2_CXXFLAGS += @MUTATIONPP_CXX@
___bin_pySU2_LDADD += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
___bin_pySU2_CXXFLAGS += @JSONCPP_CXX@
___bin_pySU2_LDADD += @JSONCPP_LD@
# endif
