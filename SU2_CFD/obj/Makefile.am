################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_CFD
# \author M. Colonno, T. Economon, F. Palacios
# \version 5.0.0 "Raven"
#
# SU2 Original Developers: Dr. Francisco D. Palacios.
#                          Dr. Thomas D. Economon.
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

bin_PROGRAMS =

if BUILD_NORMAL
bin_PROGRAMS += ../bin/SU2_CFD
endif

if BUILD_DIRECTDIFF
bin_PROGRAMS += ../bin/SU2_CFD_DIRECTDIFF
endif

if BUILD_REVERSE
bin_PROGRAMS += ../bin/SU2_CFD_AD
endif

noinst_LIBRARIES= 

if BUILD_NORMAL
noinst_LIBRARIES+= libSU2Core.a
endif

libSU2Core_sources = \
  ../include/definition_structure.hpp \
  ../include/fluid_model.hpp \
  ../include/fluid_model.inl \
  ../include/integration_structure.hpp \
  ../include/integration_structure.inl \
  ../include/driver_structure.hpp \
  ../include/iteration_structure.hpp \
  ../include/numerics_structure.hpp \
  ../include/numerics_structure.inl \
  ../include/output_structure.hpp \
  ../include/output_structure.inl \
  ../include/solver_structure.hpp \
  ../include/solver_structure.inl \
  ../include/SU2_CFD.hpp \
  ../include/transport_model.hpp \
  ../include/transport_model.inl \
  ../include/variable_structure.hpp \
  ../include/variable_structure.inl \
  ../include/transfer_structure.hpp \
  ../include/transfer_structure.inl \
  ../src/definition_structure.cpp \
  ../src/fluid_model.cpp \
  ../src/fluid_model_pig.cpp \
  ../src/fluid_model_pvdw.cpp \
  ../src/fluid_model_ppr.cpp \
  ../src/integration_structure.cpp \
  ../src/integration_time.cpp \
  ../src/driver_structure.cpp \
  ../src/iteration_structure.cpp \
  ../src/numerics_adjoint_mean.cpp \
  ../src/numerics_adjoint_mean_inc.cpp \
  ../src/numerics_adjoint_turbulent.cpp \
  ../src/numerics_direct_poisson.cpp \
  ../src/numerics_direct_heat.cpp \
  ../src/numerics_direct_mean.cpp \
  ../src/numerics_direct_mean_inc.cpp \
  ../src/numerics_direct_transition.cpp \
  ../src/numerics_direct_turbulent.cpp \
  ../src/numerics_direct_wave.cpp \
  ../src/numerics_direct_elasticity_nonlinear.cpp \
  ../src/numerics_direct_elasticity_linear.cpp \
  ../src/numerics_direct_elasticity.cpp \
  ../src/numerics_structure.cpp \
  ../src/numerics_template.cpp \
  ../src/output_cgns.cpp \
  ../src/output_structure.cpp \
  ../src/output_physics.cpp \
  ../src/output_tecplot.cpp \
  ../src/output_fieldview.cpp \
  ../src/output_su2.cpp \
  ../src/output_paraview.cpp \
  ../src/solver_adjoint_mean.cpp \
  ../src/solver_adjoint_mean_inc.cpp \
  ../src/solver_adjoint_turbulent.cpp \
  ../src/solver_adjoint_discrete.cpp \
  ../src/solver_adjoint_elasticity.cpp \
  ../src/solver_direct_poisson.cpp \
  ../src/solver_direct_heat.cpp \
  ../src/solver_direct_mean.cpp \
  ../src/solver_direct_mean_inc.cpp \
  ../src/solver_direct_transition.cpp \
  ../src/solver_direct_turbulent.cpp \
  ../src/solver_direct_wave.cpp \
  ../src/solver_direct_elasticity.cpp \
  ../src/solver_structure.cpp \
  ../src/solver_template.cpp \
  ../src/transfer_physics.cpp \
  ../src/transfer_structure.cpp \
  ../src/transport_model.cpp \
  ../src/variable_adjoint_mean.cpp \
  ../src/variable_adjoint_mean_inc.cpp \
  ../src/variable_adjoint_turbulent.cpp \
  ../src/variable_adjoint_discrete.cpp \
  ../src/variable_adjoint_elasticity.cpp \
  ../src/variable_direct_poisson.cpp \
  ../src/variable_direct_heat.cpp \
  ../src/variable_direct_mean.cpp \
  ../src/variable_direct_mean_inc.cpp \
  ../src/variable_direct_transition.cpp \
  ../src/variable_direct_turbulent.cpp \
  ../src/variable_direct_wave.cpp \
  ../src/variable_direct_elasticity.cpp \
  ../src/variable_structure.cpp \
  ../src/variable_template.cpp

su2_cfd_sources = \
  ../include/SU2_CFD.hpp \
  ../src/SU2_CFD.cpp

libSU2Core_cxx_flags = -fPIC
libSU2Core_libadd = 

su2_cfd_cxx_flags =
su2_cfd_ldadd =

# always link to built dependencies from ./externals
su2_cfd_cxx_flags += @su2_externals_INCLUDES@
su2_cfd_ldadd += @su2_externals_LIBS@
libSU2Core_cxx_flags += @su2_externals_INCLUDES@

# if BUILD_MUTATIONPP
su2_cfd_cxx_flags += @MUTATIONPP_CXX@
su2_cfd_ldadd += @MUTATIONPP_LD@
# endif

# if BUILD_JSONCPP
su2_cfd_cxx_flags += @JSONCPP_CXX@
su2_cfd_ldadd += @JSONCPP_LD@
# endif

# if BUILD_HDF5
su2_cfd_cxx_flags += @HDF5_CXX@
su2_cfd_ldadd += @HDF5_LD@
# endif

# if BUILD_SZIP
su2_cfd_cxx_flags += @SZIP_CXX@
su2_cfd_ldadd += @SZIP_LD@
# endif

# if BUILD_ZLIB
su2_cfd_cxx_flags += @ZLIB_CXX@
su2_cfd_ldadd += @ZLIB_LD@
# endif

if BUILD_NORMAL
libSU2Core_a_SOURCES = $(libSU2Core_sources)
libSU2Core_a_CXXFLAGS = $(libSU2Core_cxx_flags)
libSU2Core_a_LIBADD = $(libSU2Core_libadd)
___bin_SU2_CFD_SOURCES = $(su2_cfd_sources)
___bin_SU2_CFD_CXXFLAGS = ${su2_cfd_cxx_flags}
___bin_SU2_CFD_LDADD = libSU2Core.a ../../Common/lib/libSU2.a ${su2_cfd_ldadd}
endif

if BUILD_DIRECTDIFF
___bin_SU2_CFD_DIRECTDIFF_SOURCES = $(su2_cfd_sources) ${libSU2Core_sources}
___bin_SU2_CFD_DIRECTDIFF_CXXFLAGS = @DIRECTDIFF_CXX@ ${su2_cfd_cxx_flags}
___bin_SU2_CFD_DIRECTDIFF_LDADD = @DIRECTDIFF_LIBS@ ../../Common/lib/libSU2_DIRECTDIFF.a ${su2_cfd_ldadd}
endif

if BUILD_REVERSE
___bin_SU2_CFD_AD_SOURCES = $(su2_cfd_sources) ${libSU2Core_sources}
___bin_SU2_CFD_AD_CXXFLAGS = @REVERSE_CXX@ ${su2_cfd_cxx_flags}
___bin_SU2_CFD_AD_LDADD = @REVERSE_LIBS@ ../../Common/lib/libSU2_AD.a ${su2_cfd_ldadd}
endif
