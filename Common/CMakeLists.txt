project(SU2_common CXX)

include(${CMAKE_SOURCE_DIR}/cmake/functions.cmake)

set(SOURCES
    include/datatypes/primitive_structure.hpp
    include/datatypes/primitive_structure.inl
    # include/datatypes/complex_structure.hpp
    # include/datatypes/complex_structure.inl
    include/ad_structure.hpp
    include/ad_structure.inl
    include/fem_cgns_elements.hpp
    include/config_structure.hpp
    include/config_structure.hpp
    include/config_structure.inl
    include/blas_structure.hpp
    include/dual_grid_structure.hpp
    include/dual_grid_structure.inl
    include/fem_geometry_structure.hpp
    include/fem_geometry_structure.inl
    include/fem_standard_element.hpp
    include/fem_standard_element.inl
    include/geometry_structure.hpp
    include/geometry_structure.inl
    include/CMultiGridQueue.hpp
    include/CMeshReaderFVM.hpp
    include/CSU2ASCIIMeshReaderFVM.hpp
    include/CCGNSMeshReaderFVM.hpp
    include/CRectangularMeshReaderFVM.hpp
    include/CBoxMeshReaderFVM.hpp
    include/graph_coloring_structure.hpp
    include/grid_adaptation_structure.hpp
    include/grid_adaptation_structure.inl
    include/grid_movement_structure.hpp
    include/grid_movement_structure.inl
    include/option_structure.hpp
    include/primal_grid_structure.hpp
    include/primal_grid_structure.inl
    include/mpi_structure.hpp
    include/mpi_structure.inl
    include/datatype_structure.hpp
    include/datatype_structure.inl
    include/interpolation_structure.hpp
    include/fem_gauss_jacobi_quadrature.hpp
    include/fem_gauss_jacobi_quadrature.inl
    include/gauss_structure.hpp
    include/gauss_structure.inl
    include/element_structure.hpp
    include/element_structure.inl
    include/adt_structure.hpp
    include/adt_structure.inl
    include/wall_model.hpp
    include/wall_model.inl
    include/toolboxes/printing_toolbox.hpp
    include/toolboxes/CLinearPartitioner.hpp
    include/toolboxes/MMS/CVerificationSolution.hpp
    include/toolboxes/MMS/CVerificationSolution.inl
    include/toolboxes/MMS/CIncTGVSolution.hpp
    include/toolboxes/MMS/CInviscidVortexSolution.hpp
    include/toolboxes/MMS/CMMSIncEulerSolution.hpp
    include/toolboxes/MMS/CMMSIncNSSolution.hpp
    include/toolboxes/MMS/CMMSNSTwoHalfCirclesSolution.hpp
    include/toolboxes/MMS/CMMSNSTwoHalfSpheresSolution.hpp
    include/toolboxes/MMS/CMMSNSUnitQuadSolution.hpp
    include/toolboxes/MMS/CMMSNSUnitQuadSolutionWallBC.hpp
    include/toolboxes/MMS/CNSUnitQuadSolution.hpp
    include/toolboxes/MMS/CRinglebSolution.hpp
    include/toolboxes/MMS/CTGVSolution.hpp
    include/toolboxes/MMS/CUserDefinedSolution.hpp
    include/linear_algebra/CSysVector.hpp
    include/linear_algebra/CSysMatrix.hpp
    include/linear_algebra/CSysMatrix.inl
    include/linear_algebra/CMatrixVectorProduct.hpp
    include/linear_algebra/CPreconditioner.hpp
    include/linear_algebra/CSysSolve.hpp
    include/linear_algebra/CSysSolve_b.hpp
    src/fem_cgns_elements.cpp
    src/config_structure.cpp
    src/blas_structure.cpp
    src/dual_grid_structure.cpp
    src/fem_geometry_structure.cpp
    src/fem_integration_rules.cpp
    src/fem_standard_element.cpp
    src/fem_wall_distance.cpp
    src/fem_work_estimate_metis.cpp
    src/geometry_structure.cpp
    src/CMultiGridQueue.cpp
    src/CMeshReaderFVM.cpp
    src/CSU2ASCIIMeshReaderFVM.cpp
    src/CCGNSMeshReaderFVM.cpp
    src/CRectangularMeshReaderFVM.cpp
    src/CBoxMeshReaderFVM.cpp
    src/geometry_structure_fem_part.cpp
    src/graph_coloring_structure.cpp
    src/grid_adaptation_structure.cpp
    src/grid_movement_structure.cpp
    src/primal_grid_structure.cpp
    src/mpi_structure.cpp
    src/ad_structure.cpp
    src/fem_gauss_jacobi_quadrature.cpp
    src/gauss_structure.cpp
    src/element_structure.cpp
    src/element_linear.cpp
    src/interpolation_structure.cpp
    src/adt_structure.cpp
    src/wall_model.cpp
    src/toolboxes/printing_toolbox.cpp
    src/toolboxes/signal_processing_toolbox.cpp
    src/toolboxes/CLinearPartitioner.cpp
    src/toolboxes/MMS/CVerificationSolution.cpp
    src/toolboxes/MMS/CIncTGVSolution.cpp
    src/toolboxes/MMS/CInviscidVortexSolution.cpp
    src/toolboxes/MMS/CMMSIncEulerSolution.cpp
    src/toolboxes/MMS/CMMSIncNSSolution.cpp
    src/toolboxes/MMS/CMMSNSTwoHalfCirclesSolution.cpp
    src/toolboxes/MMS/CMMSNSTwoHalfSpheresSolution.cpp
    src/toolboxes/MMS/CMMSNSUnitQuadSolution.cpp
    src/toolboxes/MMS/CMMSNSUnitQuadSolutionWallBC.cpp
    src/toolboxes/MMS/CNSUnitQuadSolution.cpp
    src/toolboxes/MMS/CRinglebSolution.cpp
    src/toolboxes/MMS/CTGVSolution.cpp
    src/toolboxes/MMS/CUserDefinedSolution.cpp
    src/linear_algebra/CSysVector.cpp
    src/linear_algebra/CSysMatrix.cpp
    src/linear_algebra/CSysSolve.cpp
    src/linear_algebra/CSysSolve_b.cpp
    src/linear_algebra/CPastixWrapper.cpp)

if(SU2_BUILD_NORMAL)
  set(TARGET_NAME SU2)
elseif(SU2_BUILD_DIRECTDIFF)
  set(TARGET_NAME SU2_DIRECTDIFF)
elseif(SU2_BUILD_REVERSE)
  set(TARGET_NAME SU2_AD)
else()
  message(
    FATAL_ERROR "Invalid SU2 build type, neither normal, directdiff or reverse")
endif()

add_library(${TARGET_NAME} STATIC ${SOURCES})
target_include_directories(${TARGET_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)

if(SU2_BUILD_NORMAL)
  su2_setup_target(${TARGET_NAME})
else()
  su2_setup_target(${TARGET_NAME} LIBS WITH_CODI)
endif()

add_library(SU2::SU2 ALIAS ${TARGET_NAME})
