################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_CFD
# \author M. Colonno, T. Economon, F. Palacios
# \version 7.0.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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
noinst_LIBRARIES =

if BUILD_NORMAL
bin_PROGRAMS += ../bin/SU2_CFD
noinst_LIBRARIES+= libSU2Core.a
endif

if BUILD_DIRECTDIFF
bin_PROGRAMS += ../bin/SU2_CFD_DIRECTDIFF
noinst_LIBRARIES+= libSU2Core_DIRECTDIFF.a
endif

if BUILD_REVERSE
bin_PROGRAMS += ../bin/SU2_CFD_AD
noinst_LIBRARIES+= libSU2Core_AD.a
endif

libSU2Core_sources = ../src/definition_structure.cpp \
  ../src/fluid_model.cpp \
  ../src/fluid_model_pig.cpp \
  ../src/fluid_model_pvdw.cpp \
  ../src/fluid_model_ppr.cpp \
  ../src/fluid_model_inc.cpp \
  ../src/integration_structure.cpp \
  ../src/integration_time.cpp \
  ../src/drivers/CMultizoneDriver.cpp \
  ../src/drivers/CSinglezoneDriver.cpp \
  ../src/drivers/CDiscAdjSinglezoneDriver.cpp \
  ../src/drivers/CDiscAdjMultizoneDriver.cpp \
  ../src/drivers/CDriver.cpp \
  ../src/drivers/CDummyDriver.cpp \
  ../src/iteration_structure.cpp \
  ../src/numerics_adjoint_mean.cpp \
  ../src/numerics_adjoint_turbulent.cpp \
  ../src/numerics_direct_heat.cpp \
  ../src/numerics_direct_mean.cpp \
  ../src/numerics_direct_mean_inc.cpp \
  ../src/numerics_direct_transition.cpp \
  ../src/numerics_direct_turbulent.cpp \
  ../src/numerics/elasticity/CFEAElasticity.cpp \
  ../src/numerics/elasticity/CFEALinearElasticity.cpp \
  ../src/numerics/elasticity/CFEAMeshElasticity.cpp \
  ../src/numerics/elasticity/CFEANonlinearElasticity.cpp \
  ../src/numerics/elasticity/CFEM_NeoHookean_Comp.cpp \
  ../src/numerics/elasticity/CFEM_Knowles_NearInc.cpp \
  ../src/numerics/elasticity/CFEM_IdealDE.cpp \
  ../src/numerics/elasticity/CFEM_DielectricElastomer.cpp \
  ../src/numerics_structure.cpp \
  ../src/numerics_template.cpp \
  ../src/output/filewriter/CCSVFileWriter.cpp \
  ../src/output/filewriter/CFEMDataSorter.cpp \
  ../src/output/filewriter/CFVMDataSorter.cpp \
  ../src/output/filewriter/CParallelDataSorter.cpp \
  ../src/output/filewriter/CParallelFileWriter.cpp \
  ../src/output/filewriter/CParaviewBinaryFileWriter.cpp \
  ../src/output/filewriter/CParaviewFileWriter.cpp \
  ../src/output/filewriter/CSurfaceFEMDataSorter.cpp \
  ../src/output/filewriter/CSurfaceFVMDataSorter.cpp \
  ../src/output/filewriter/CSU2BinaryFileWriter.cpp \
  ../src/output/filewriter/CSU2FileWriter.cpp \
  ../src/output/filewriter/CSU2MeshFileWriter.cpp \
  ../src/output/filewriter/CTecplotFileWriter.cpp \
  ../src/output/filewriter/CTecplotBinaryFileWriter.cpp \
  ../src/output/COutput.cpp \
  ../src/output/output_physics.cpp \
  ../src/output/CMeshOutput.cpp \
  ../src/output/CElasticityOutput.cpp \
  ../src/output/CFlowOutput.cpp \
  ../src/output/CFlowCompOutput.cpp \
  ../src/output/CFlowCompFEMOutput.cpp \
  ../src/output/CFlowIncOutput.cpp \
  ../src/output/CHeatOutput.cpp \
  ../src/output/CBaselineOutput.cpp \
  ../src/output/CAdjElasticityOutput.cpp \
  ../src/output/CAdjHeatOutput.cpp \
  ../src/output/CAdjFlowCompOutput.cpp \
  ../src/output/CAdjFlowIncOutput.cpp \
  ../src/output/CMultizoneOutput.cpp \
  ../src/output/output_structure_legacy.cpp \
  ../src/python_wrapper_structure.cpp \
  ../src/solver_adjoint_mean.cpp \
  ../src/solver_adjoint_turbulent.cpp \
  ../src/solver_adjoint_discrete.cpp \
  ../src/solver_adjoint_elasticity.cpp \
  ../src/solvers/CDiscAdjMeshSolver.cpp \
  ../src/solver_direct_heat.cpp \
  ../src/solver_direct_mean.cpp \
  ../src/solver_direct_mean_fem.cpp \
  ../src/solver_direct_mean_inc.cpp \
  ../src/solver_direct_transition.cpp \
  ../src/solver_direct_turbulent.cpp \
  ../src/solver_direct_elasticity.cpp \
  ../src/solvers/CMeshSolver.cpp \
  ../src/solver_structure.cpp \
  ../src/solver_template.cpp \
  ../src/CMarkerProfileReaderFVM.cpp \
  ../src/interfaces/CInterface.cpp \
  ../src/interfaces/cfd/CConservativeVarsInterface.cpp \
  ../src/interfaces/cfd/CMixingPlaneInterface.cpp \
  ../src/interfaces/cfd/CSlidingInterface.cpp \
  ../src/interfaces/cht/CConjugateHeatInterface.cpp \
  ../src/interfaces/fsi/CDisplacementsInterface.cpp \
  ../src/interfaces/fsi/CFlowTractionInterface.cpp \
  ../src/interfaces/fsi/CDiscAdjFlowTractionInterface.cpp \
  ../src/interfaces/fsi/CDisplacementsInterfaceLegacy.cpp \
  ../src/interfaces/fsi/CDiscAdjDisplacementsInterfaceLegacy.cpp \
  ../src/transport_model.cpp \
  ../src/variables/CFEABoundVariable.cpp \
  ../src/variables/CDiscAdjMeshBoundVariable.cpp \
  ../src/variables/CMeshBoundVariable.cpp \
  ../src/variables/CMeshElement.cpp \
  ../src/variables/CMeshVariable.cpp \
  ../src/variables/CHeatFVMVariable.cpp \
  ../src/variables/CVariable.cpp \
  ../src/variables/CAdjNSVariable.cpp \
  ../src/variables/CTurbSSTVariable.cpp \
  ../src/variables/CAdjTurbVariable.cpp \
  ../src/variables/CTransLMVariable.cpp \
  ../src/variables/CDiscAdjFEABoundVariable.cpp \
  ../src/variables/CDiscAdjFEAVariable.cpp \
  ../src/variables/CIncEulerVariable.cpp \
  ../src/variables/CTurbVariable.cpp \
  ../src/variables/CNSVariable.cpp \
  ../src/variables/CBaselineVariable.cpp \
  ../src/variables/CTurbSAVariable.cpp \
  ../src/variables/CFEAVariable.cpp \
  ../src/variables/CAdjEulerVariable.cpp \
  ../src/variables/CDiscAdjVariable.cpp \
  ../src/variables/CIncNSVariable.cpp \
  ../src/variables/CEulerVariable.cpp \
  ../src/limiters/computeLimiters.cpp

su2_cfd_sources = \
  ../include/SU2_CFD.hpp \
  ../src/SU2_CFD.cpp

libSU2Core_cxx_flags = -fPIC -std=c++11
libSU2Core_libadd = 

su2_cfd_cxx_flags = -fPIC -std=c++11
su2_cfd_ldadd =

# always link to built dependencies from ./externals
su2_cfd_cxx_flags += @su2_externals_INCLUDES@
su2_cfd_ldadd += @su2_externals_LIBS@
su2_cfd_ldadd += @su2_externals_LIBPTHREAD@
libSU2Core_cxx_flags += @su2_externals_INCLUDES@

# if BUILD_MUTATIONPP
su2_cfd_cxx_flags += @MUTATIONPP_CXX@
su2_cfd_ldadd += @MUTATIONPP_LD@
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
libSU2Core_DIRECTDIFF_a_SOURCES = $(libSU2Core_sources)
libSU2Core_DIRECTDIFF_a_CXXFLAGS = @DIRECTDIFF_CXX@ $(libSU2Core_cxx_flags)
libSU2Core_DIRECTDIFF_a_LIBADD = @DIRECTDIFF_LIBS@ $(libSU2Core_libadd)
___bin_SU2_CFD_DIRECTDIFF_SOURCES = $(su2_cfd_sources) 
___bin_SU2_CFD_DIRECTDIFF_CXXFLAGS = @DIRECTDIFF_CXX@ ${su2_cfd_cxx_flags}
___bin_SU2_CFD_DIRECTDIFF_LDADD = libSU2Core_DIRECTDIFF.a @DIRECTDIFF_LIBS@ ../../Common/lib/libSU2_DIRECTDIFF.a ${su2_cfd_ldadd}
endif

if BUILD_REVERSE
libSU2Core_AD_a_SOURCES = $(libSU2Core_sources)
libSU2Core_AD_a_CXXFLAGS = @REVERSE_CXX@ $(libSU2Core_cxx_flags)
libSU2Core_AD_a_LIBADD = @REVERSE_LIBS@ $(libSU2Core_libadd)
___bin_SU2_CFD_AD_SOURCES = $(su2_cfd_sources) 
___bin_SU2_CFD_AD_CXXFLAGS = @REVERSE_CXX@ ${su2_cfd_cxx_flags}
___bin_SU2_CFD_AD_LDADD = libSU2Core_AD.a @REVERSE_LIBS@ ../../Common/lib/libSU2_AD.a ${su2_cfd_ldadd}
endif
