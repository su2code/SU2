# SU^2 Makefile.am for CFD
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated: 	3/11/13
#
# Created by: 		M. Colonno
# Updated by: 		M. Colonno
#
##########################

# AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = ../bin/SU2_CFD

___bin_SU2_CFD_SOURCES = ../include/definition_structure.hpp \
		../include/integration_structure.hpp \
		../include/integration_structure.inl \
		../include/iteration_structure.hpp \
		../include/math_ad.hpp \
		../include/numerics_structure.hpp \
		../include/numerics_structure.inl \
		../include/output_structure.hpp \
		../include/solver_structure.hpp \
		../include/solver_structure.inl \
		../include/SU2_CFD.hpp \
		../include/variable_structure.hpp \
		../include/variable_structure.inl \
		../src/definition_structure.cpp \
		../src/integration_notime.cpp \
		../src/integration_structure.cpp \
		../src/integration_time.cpp \
		../src/iteration_structure.cpp \
		../src/math_ad.cpp \
		../src/numerics_adjoint_electric.cpp \
		../src/numerics_adjoint_levelset.cpp \
		../src/numerics_adjoint_mean.cpp \
		../src/numerics_adjoint_plasma.cpp \
		../src/numerics_adjoint_tne2.cpp \
		../src/numerics_adjoint_turbulent.cpp \
		../src/numerics_direct_electric.cpp \
		../src/numerics_direct_fea.cpp \
		../src/numerics_direct_heat.cpp \
		../src/numerics_direct_levelset.cpp \
		../src/numerics_direct_mean.cpp \
		../src/numerics_direct_plasma.cpp \
		../src/numerics_direct_tne2.cpp \
		../src/numerics_direct_transition.cpp \
		../src/numerics_direct_turbulent.cpp \
		../src/numerics_direct_wave.cpp \
		../src/numerics_linearized_mean.cpp \
		../src/numerics_linearized_turbulent.cpp \
		../src/numerics_structure.cpp \
		../src/numerics_template.cpp \
		../src/output_cgns.cpp \
		../src/output_structure.cpp \
		../src/output_tecplot.cpp \
		../src/solver_adjoint_levelset.cpp \
		../src/solver_adjoint_mean.cpp \
		../src/solver_adjoint_plasma.cpp \
		../src/solver_adjoint_tne2.cpp \
		../src/solver_adjoint_turbulent.cpp \
		../src/solver_direct_electric.cpp \
		../src/solver_direct_fea.cpp \
		../src/solver_direct_heat.cpp \
		../src/solver_direct_levelset.cpp \
		../src/solver_direct_mean.cpp \
		../src/solver_direct_plasma.cpp \
		../src/solver_direct_tne2.cpp \
		../src/solver_direct_transition.cpp \
		../src/solver_direct_turbulent.cpp \
		../src/solver_direct_wave.cpp \
		../src/solver_linearized_mean.cpp \
		../src/solver_linearized_turbulent.cpp \
		../src/solver_structure.cpp \
		../src/solver_template.cpp \
		../src/SU2_CFD.cpp \
		../src/variable_adjoint_electric.cpp \
		../src/variable_adjoint_levelset.cpp \
		../src/variable_adjoint_mean.cpp \
		../src/variable_adjoint_plasma.cpp \
		../src/variable_adjoint_tne2.cpp \
		../src/variable_adjoint_turbulent.cpp \
		../src/variable_direct_electric.cpp \
		../src/variable_direct_fea.cpp \
		../src/variable_direct_heat.cpp \
		../src/variable_direct_levelset.cpp \
		../src/variable_direct_mean.cpp \
		../src/variable_direct_plasma.cpp \
		../src/variable_direct_tne2.cpp \
		../src/variable_direct_transition.cpp \
		../src/variable_direct_turbulent.cpp \
		../src/variable_direct_wave.cpp \
		../src/variable_linearized_mean.cpp \
		../src/variable_linearized_turbulent.cpp \
		../src/variable_structure.cpp \
		../src/variable_template.cpp

___bin_SU2_CFD_CXXFLAGS =
___bin_SU2_CFD_LDADD = ../../Common/lib/libSU2.a

# if BUILD_CGNS
___bin_SU2_CFD_CXXFLAGS += @CGNS_CXX@
___bin_SU2_CFD_LDADD += @CGNS_LD@
# endif

# if BUILD_METIS
___bin_SU2_CFD_CXXFLAGS += @METIS_CXX@
___bin_SU2_CFD_LDADD += @METIS_LD@
# endif

# if BUILD_TECIO
___bin_SU2_CFD_CXXFLAGS += @TECIO_CXX@
___bin_SU2_CFD_LDADD += @TECIO_LD@
# endif
