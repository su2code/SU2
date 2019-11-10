project(SU2_CFD CXX)

include(${CMAKE_SOURCE_DIR}/cmake/functions.cmake)

set(SU2_CORE_SOURCES
    src/definition_structure.cpp
    src/fluid_model.cpp
    src/fluid_model_pig.cpp
    src/fluid_model_pvdw.cpp
    src/fluid_model_ppr.cpp
    src/fluid_model_inc.cpp
    src/integration_structure.cpp
    src/integration_time.cpp
    src/drivers/CMultizoneDriver.cpp
    src/drivers/CSinglezoneDriver.cpp
    src/drivers/CDiscAdjSinglezoneDriver.cpp
    src/drivers/CDiscAdjMultizoneDriver.cpp
    src/drivers/CDriver.cpp
    src/drivers/CDummyDriver.cpp
    src/iteration_structure.cpp
    src/numerics_adjoint_mean.cpp
    src/numerics_adjoint_turbulent.cpp
    src/numerics_direct_heat.cpp
    src/numerics_direct_mean.cpp
    src/numerics_direct_mean_inc.cpp
    src/numerics_direct_transition.cpp
    src/numerics_direct_turbulent.cpp
    src/numerics_direct_elasticity_nonlinear.cpp
    src/numerics_direct_elasticity_linear.cpp
    src/numerics_direct_elasticity.cpp
    src/numerics/CFEAMeshElasticity.cpp
    src/numerics_structure.cpp
    src/numerics_template.cpp
    src/output/filewriter/CCSVFileWriter.cpp
    src/output/filewriter/CFEMDataSorter.cpp
    src/output/filewriter/CFVMDataSorter.cpp
    src/output/filewriter/CParallelDataSorter.cpp
    src/output/filewriter/CParallelFileWriter.cpp
    src/output/filewriter/CParaviewBinaryFileWriter.cpp
    src/output/filewriter/CParaviewFileWriter.cpp
    src/output/filewriter/CSurfaceFEMDataSorter.cpp
    src/output/filewriter/CSurfaceFVMDataSorter.cpp
    src/output/filewriter/CSU2BinaryFileWriter.cpp
    src/output/filewriter/CSU2FileWriter.cpp
    src/output/filewriter/CSU2MeshFileWriter.cpp
    src/output/filewriter/CTecplotFileWriter.cpp
    src/output/filewriter/CTecplotBinaryFileWriter.cpp
    src/output/COutput.cpp
    src/output/output_physics.cpp
    src/output/CMeshOutput.cpp
    src/output/CElasticityOutput.cpp
    src/output/CFlowOutput.cpp
    src/output/CFlowCompOutput.cpp
    src/output/CFlowCompFEMOutput.cpp
    src/output/CFlowIncOutput.cpp
    src/output/CHeatOutput.cpp
    src/output/CBaselineOutput.cpp
    src/output/CAdjElasticityOutput.cpp
    src/output/CAdjHeatOutput.cpp
    src/output/CAdjFlowCompOutput.cpp
    src/output/CAdjFlowIncOutput.cpp
    src/output/CMultizoneOutput.cpp
    src/output/output_structure_legacy.cpp
    src/python_wrapper_structure.cpp
    src/solver_adjoint_mean.cpp
    src/solver_adjoint_turbulent.cpp
    src/solver_adjoint_discrete.cpp
    src/solver_adjoint_elasticity.cpp
    src/solvers/CDiscAdjMeshSolver.cpp
    src/solver_direct_heat.cpp
    src/solver_direct_mean.cpp
    src/solver_direct_mean_fem.cpp
    src/solver_direct_mean_inc.cpp
    src/solver_direct_transition.cpp
    src/solver_direct_turbulent.cpp
    src/solver_direct_elasticity.cpp
    src/solvers/CMeshSolver.cpp
    src/solver_structure.cpp
    src/solver_template.cpp
    src/interfaces/CInterface.cpp
    src/interfaces/cfd/CConservativeVarsInterface.cpp
    src/interfaces/cfd/CMixingPlaneInterface.cpp
    src/interfaces/cfd/CSlidingInterface.cpp
    src/interfaces/cht/CConjugateHeatInterface.cpp
    src/interfaces/fsi/CDisplacementsInterface.cpp
    src/interfaces/fsi/CFlowTractionInterface.cpp
    src/interfaces/fsi/CDiscAdjFlowTractionInterface.cpp
    src/interfaces/fsi/CDisplacementsInterfaceLegacy.cpp
    src/interfaces/fsi/CDiscAdjDisplacementsInterfaceLegacy.cpp
    src/transport_model.cpp
    src/variables/CFEABoundVariable.cpp
    src/variables/CDiscAdjMeshBoundVariable.cpp
    src/variables/CMeshBoundVariable.cpp
    src/variables/CMeshElement.cpp
    src/variables/CMeshVariable.cpp
    src/variables/CHeatFVMVariable.cpp
    src/variables/CVariable.cpp
    src/variables/CAdjNSVariable.cpp
    src/variables/CTurbSSTVariable.cpp
    src/variables/CAdjTurbVariable.cpp
    src/variables/CTransLMVariable.cpp
    src/variables/CDiscAdjFEABoundVariable.cpp
    src/variables/CDiscAdjFEAVariable.cpp
    src/variables/CIncEulerVariable.cpp
    src/variables/CTurbVariable.cpp
    src/variables/CNSVariable.cpp
    src/variables/CBaselineVariable.cpp
    src/variables/CTurbSAVariable.cpp
    src/variables/CFEAVariable.cpp
    src/variables/CAdjEulerVariable.cpp
    src/variables/CDiscAdjVariable.cpp
    src/variables/CIncNSVariable.cpp
    src/variables/CEulerVariable.cpp)

set(SU2_CFD_SOURCES include/SU2_CFD.hpp src/SU2_CFD.cpp)

if(SU2_BUILD_NORMAL)
  set(SU2_CORE_LIB SU2Core)
  set(SU2_CFD_EXE SU2_CFD)
elseif(SU2_BUILD_DIRECTDIFF)
  set(SU2_CORE_LIB SU2Core_DIRECTDIFF)
  set(SU2_CFD_EXE SU2_CFD_DIRECTDIFF)
elseif(SU2_BUILD_REVERSE)
  set(SU2_CORE_LIB SU2Core_AD)
  set(SU2_CFD_EXE SU2_CFD_AD)
else()
  message(
    FATAL_ERROR "Invalid SU2 build type, neither normal, directdiff or reverse")
endif()

add_library(${SU2_CORE_LIB} STATIC ${SU2_CORE_SOURCES})
target_include_directories(${SU2_CORE_LIB} PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
add_library(SU2::SU2Core ALIAS ${SU2_CORE_LIB})

add_executable(${SU2_CFD_EXE} ${SU2_CFD_SOURCES})

if (SU2_BUILD_NORMAL)
  su2_setup_target(${SU2_CORE_LIB})
  su2_setup_target(${SU2_CFD_EXE} LIBS SU2::SU2Core SU2::SU2 WITH_THREADS)
else()
  su2_setup_target(${SU2_CORE_LIB} LIBS WITH_CODI)
  su2_setup_target(${SU2_CFD_EXE} LIBS SU2::SU2Core SU2::SU2 WITH_THREADS WITH_CODI)
endif()
