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
		../include/solution_structure.hpp \
		../include/solution_structure.inl \
		../include/SU2_CFD.hpp \
		../include/variable_structure.hpp \
		../include/variable_structure.inl \
		../src/definition_structure.cpp \
		../src/integration_notime.cpp \
		../src/integration_structure.cpp \
		../src/integration_time.cpp \
		../src/iteration_structure.cpp \
		../src/math_ad.cpp \
		../src/numerics_convective_ad.cpp \
		../src/numerics_convective.cpp \
		../src/numerics_source_ad.cpp \
		../src/numerics_source.cpp \
		../src/numerics_structure.cpp \
		../src/numerics_viscous.cpp \
		../src/numerics_viscous_ad.cpp \
		../src/output_structure.cpp \
		../src/output_cgns.cpp \
		../src/output_tecplot.cpp \
		../src/solution_adjoint_levelset.cpp \
		../src/solution_adjoint_mean.cpp \
		../src/solution_adjoint_plasma.cpp \
		../src/solution_adjoint_turbulent.cpp \
		../src/solution_direct_electric.cpp \
		../src/solution_direct_fea.cpp \
		../src/solution_direct_levelset.cpp \
		../src/solution_direct_mean.cpp \
		../src/solution_direct_plasma.cpp \
		../src/solution_direct_transition.cpp \
		../src/solution_direct_turbulent.cpp \
		../src/solution_direct_turbulent_ad.cpp \
		../src/solution_direct_wave.cpp \
		../src/solution_linearized_mean.cpp \
		../src/solution_linearized_turbulent.cpp \
		../src/solution_structure.cpp \
		../src/solution_template.cpp \
		../src/SU2_CFD.cpp \
		../src/variable_adjoint.cpp \
		../src/variable_direct.cpp \
		../src/variable_linearized.cpp \
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
