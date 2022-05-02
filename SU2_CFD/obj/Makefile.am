################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_CFD
# \author M. Colonno, T. Economon, F. Palacios
# \version 7.3.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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
  ../src/fluid/CFluidModel.cpp \
  ../src/fluid/CIdealGas.cpp \
  ../src/fluid/CPengRobinson.cpp \
  ../src/fluid/CVanDerWaalsGas.cpp \
  ../src/fluid/CNEMOGas.cpp \
  ../src/fluid/CSU2TCLib.cpp \
  ../src/fluid/CMutationTCLib.cpp \
  ../src/integration/CIntegration.cpp \
  ../src/integration/CSingleGridIntegration.cpp \
  ../src/integration/CMultiGridIntegration.cpp \
  ../src/integration/CNewtonIntegration.cpp \
  ../src/integration/CStructuralIntegration.cpp \
  ../src/integration/CFEM_DG_Integration.cpp \
  ../src/integration/CIntegrationFactory.cpp \
  ../src/drivers/CMultizoneDriver.cpp \
  ../src/drivers/CSinglezoneDriver.cpp \
  ../src/drivers/CDiscAdjSinglezoneDriver.cpp \
  ../src/drivers/CDiscAdjMultizoneDriver.cpp \
  ../src/drivers/CDriver.cpp \
  ../src/drivers/CDummyDriver.cpp \
  ../src/iteration/CIteration.cpp \
  ../src/iteration/CIterationFactory.cpp \
  ../src/iteration/CAdjFluidIteration.cpp \
  ../src/iteration/CDiscAdjFluidIteration.cpp \
  ../src/iteration/CDiscAdjHeatIteration.cpp \
  ../src/iteration/CDiscAdjFEAIteration.cpp \
  ../src/iteration/CFEAIteration.cpp \
  ../src/iteration/CFEMFluidIteration.cpp \
  ../src/iteration/CFluidIteration.cpp \
  ../src/iteration/CHeatIteration.cpp \
  ../src/iteration/CTurboIteration.cpp \
  ../src/numerics/CNumerics.cpp \
  ../src/numerics/template.cpp \
  ../src/numerics/transition.cpp \
  ../src/numerics/heat.cpp \
  ../src/numerics/radiation.cpp \
  ../src/numerics/flow/convection/roe.cpp \
  ../src/numerics/flow/convection/fds.cpp \
  ../src/numerics/flow/convection/fvs.cpp \
  ../src/numerics/flow/convection/cusp.cpp \
  ../src/numerics/flow/convection/hllc.cpp \
  ../src/numerics/flow/convection/ausm_slau.cpp \
  ../src/numerics/flow/convection/centered.cpp \
  ../src/numerics/flow/flow_diffusion.cpp \
  ../src/numerics/flow/flow_sources.cpp \
  ../src/numerics/NEMO/convection/roe.cpp \
  ../src/numerics/NEMO/convection/ausm.cpp \
  ../src/numerics/NEMO/convection/ausmplusup2.cpp \
  ../src/numerics/NEMO/convection/ausmpwplus.cpp \
  ../src/numerics/NEMO/convection/lax.cpp \
  ../src/numerics/NEMO/convection/msw.cpp \
  ../src/numerics/NEMO/CNEMONumerics.cpp \
  ../src/numerics/NEMO/NEMO_diffusion.cpp \
  ../src/numerics/NEMO/NEMO_sources.cpp \
  ../src/numerics/continuous_adjoint/adj_convection.cpp \
  ../src/numerics/continuous_adjoint/adj_diffusion.cpp \
  ../src/numerics/continuous_adjoint/adj_sources.cpp \
  ../src/numerics/scalar/scalar_sources.cpp \
  ../src/numerics/species/species_sources.cpp \
  ../src/numerics/elasticity/CFEAElasticity.cpp \
  ../src/numerics/elasticity/CFEALinearElasticity.cpp \
  ../src/numerics/elasticity/CFEANonlinearElasticity.cpp \
  ../src/numerics/elasticity/nonlinear_models.cpp \
  ../include/numerics_simd/CNumericsSIMD.cpp \
  ../src/output/filewriter/CCSVFileWriter.cpp \
  ../src/output/filewriter/CSTLFileWriter.cpp \
  ../src/output/filewriter/CFEMDataSorter.cpp \
  ../src/output/filewriter/CFVMDataSorter.cpp \
  ../src/output/filewriter/CParallelDataSorter.cpp \
  ../src/output/filewriter/CParallelFileWriter.cpp \
  ../src/output/filewriter/CParaviewBinaryFileWriter.cpp \
  ../src/output/filewriter/CParaviewXMLFileWriter.cpp \
  ../src/output/filewriter/CParaviewVTMFileWriter.cpp \
  ../src/output/filewriter/CParaviewFileWriter.cpp \
  ../src/output/filewriter/CSurfaceFEMDataSorter.cpp \
  ../src/output/filewriter/CSurfaceFVMDataSorter.cpp \
  ../src/output/filewriter/CSU2BinaryFileWriter.cpp \
  ../src/output/filewriter/CSU2FileWriter.cpp \
  ../src/output/filewriter/CSU2MeshFileWriter.cpp \
  ../src/output/filewriter/CTecplotFileWriter.cpp \
  ../src/output/filewriter/CTecplotBinaryFileWriter.cpp \
  ../src/output/filewriter/CCGNSFileWriter.cpp \  
  ../src/output/tools/CWindowingTools.cpp \
  ../src/output/COutput.cpp \
  ../src/output/output_physics.cpp \
  ../src/output/CMeshOutput.cpp \
  ../src/output/CElasticityOutput.cpp \
  ../src/output/CFVMOutput.cpp \
  ../src/output/CFlowOutput.cpp \
  ../src/output/CFlowCompOutput.cpp \
  ../src/output/CFlowCompFEMOutput.cpp \
  ../src/output/CFlowIncOutput.cpp \
  ../src/output/CHeatOutput.cpp \
  ../src/output/CBaselineOutput.cpp \
  ../src/output/CAdjElasticityOutput.cpp \
  ../src/output/CAdjHeatOutput.cpp \
  ../src/output/CAdjFlowOutput.cpp \
  ../src/output/CAdjFlowCompOutput.cpp \
  ../src/output/CAdjFlowIncOutput.cpp \
  ../src/output/CMultizoneOutput.cpp \
  ../src/output/CNEMOCompOutput.cpp \
  ../src/output/COutputFactory.cpp \
  ../src/output/output_structure_legacy.cpp \
  ../src/python_wrapper_structure.cpp \
  ../src/solvers/CAdjEulerSolver.cpp \
  ../src/solvers/CAdjNSSolver.cpp \
  ../src/solvers/CAdjTurbSolver.cpp \
  ../src/solvers/CBaselineSolver.cpp \
  ../src/solvers/CBaselineSolver_FEM.cpp \
  ../src/solvers/CDiscAdjFEASolver.cpp \
  ../src/solvers/CDiscAdjMeshSolver.cpp \
  ../src/solvers/CDiscAdjSolver.cpp \
  ../src/solvers/CEulerSolver.cpp \
  ../src/solvers/CFEASolver.cpp \
  ../src/solvers/CFEM_DG_EulerSolver.cpp \
  ../src/solvers/CFEM_DG_NSSolver.cpp \
  ../src/solvers/CHeatSolver.cpp \
  ../src/solvers/CIncEulerSolver.cpp \
  ../src/solvers/CIncNSSolver.cpp \
  ../src/solvers/CMeshSolver.cpp \
  ../src/solvers/CRadSolver.cpp \
  ../src/solvers/CRadP1Solver.cpp \
  ../src/solvers/CNEMOEulerSolver.cpp \
  ../src/solvers/CNEMONSSolver.cpp \
  ../src/solvers/CNSSolver.cpp \
  ../src/solvers/CSolver.cpp \
  ../src/solvers/CTemplateSolver.cpp \
  ../src/solvers/CTransLMSolver.cpp \
  ../src/solvers/CSpeciesSolver.cpp \
  ../src/solvers/CTurbSolver.cpp \
  ../src/solvers/CTurbSASolver.cpp \
  ../src/solvers/CTurbSSTSolver.cpp \
  ../src/solvers/CSolverFactory.cpp \
  ../src/limiters/CLimiterDetails.cpp \
  ../src/CMarkerProfileReaderFVM.cpp \
  ../src/interfaces/CInterface.cpp \
  ../src/interfaces/cfd/CConservativeVarsInterface.cpp \
  ../src/interfaces/cfd/CMixingPlaneInterface.cpp \
  ../src/interfaces/cfd/CSlidingInterface.cpp \
  ../src/interfaces/cht/CConjugateHeatInterface.cpp \
  ../src/interfaces/fsi/CDisplacementsInterface.cpp \
  ../src/interfaces/fsi/CFlowTractionInterface.cpp \
  ../src/interfaces/fsi/CDiscAdjFlowTractionInterface.cpp \
  ../src/variables/CFEABoundVariable.cpp \
  ../src/variables/CRadVariable.cpp \
  ../src/variables/CRadP1Variable.cpp \
  ../src/variables/CDiscAdjMeshBoundVariable.cpp \
  ../src/variables/CMeshBoundVariable.cpp \
  ../src/variables/CMeshElement.cpp \
  ../src/variables/CMeshVariable.cpp \
  ../src/variables/CHeatVariable.cpp \
  ../src/variables/CVariable.cpp \
  ../src/variables/CAdjNSVariable.cpp \
  ../src/variables/CTurbSSTVariable.cpp \
  ../src/variables/CAdjTurbVariable.cpp \
  ../src/variables/CTransLMVariable.cpp \
  ../src/variables/CDiscAdjFEABoundVariable.cpp \
  ../src/variables/CIncEulerVariable.cpp \
  ../src/variables/CScalarVariable.cpp \
  ../src/variables/CSpeciesVariable.cpp \
  ../src/variables/CTurbVariable.cpp \
  ../src/variables/CNSVariable.cpp \
  ../src/variables/CNEMOEulerVariable.cpp \
  ../src/variables/CNEMONSVariable.cpp \
  ../src/variables/CBaselineVariable.cpp \
  ../src/variables/CTurbSAVariable.cpp \
  ../src/variables/CFEAVariable.cpp \
  ../src/variables/CAdjEulerVariable.cpp \
  ../src/variables/CDiscAdjVariable.cpp \
  ../src/variables/CIncNSVariable.cpp \
  ../src/variables/CEulerVariable.cpp \
  ../src/variables/CFlowVariable.cpp

su2_cfd_sources = \
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
