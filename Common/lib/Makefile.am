################################################################################
#
# \file Makefile.am
# \brief Makefile for the SU2 common library
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

noinst_LIBRARIES =

if BUILD_NORMAL
noinst_LIBRARIES += libSU2.a
endif

if BUILD_DIRECTDIFF
noinst_LIBRARIES += libSU2_DIRECTDIFF.a
endif

if BUILD_REVERSE
noinst_LIBRARIES += libSU2_AD.a
endif

lib_sources = \
  ../src/fem/fem_cgns_elements.cpp \
  ../src/CConfig.cpp \
  ../src/linear_algebra/blas_structure.cpp \
  ../src/fem/fem_geometry_structure.cpp \
  ../src/fem/fem_integration_rules.cpp \
  ../src/fem/fem_standard_element.cpp \
  ../src/fem/fem_wall_distance.cpp \
  ../src/fem/fem_work_estimate_metis.cpp \
  ../src/fem/geometry_structure_fem_part.cpp \
  ../src/graph_coloring_structure.cpp \
  ../src/grid_movement/CGridMovement.cpp \
  ../src/grid_movement/CFreeFormBlending.cpp \
  ../src/grid_movement/CBSplineBlending.cpp \
  ../src/grid_movement/CBezierBlending.cpp \
  ../src/grid_movement/CFreeFormDefBox.cpp \
  ../src/grid_movement/CVolumetricMovement.cpp \
  ../src/grid_movement/CSurfaceMovement.cpp \
  ../include/parallelization/mpi_structure.cpp \
  ../include/parallelization/omp_structure.cpp \
  ../src/basic_types/ad_structure.cpp \
  ../src/fem/fem_gauss_jacobi_quadrature.cpp \
  ../src/geometry/CGeometry.cpp \
  ../src/geometry/CPhysicalGeometry.cpp \
  ../src/geometry/CMultiGridGeometry.cpp \
  ../src/geometry/CMultiGridQueue.cpp \
  ../src/geometry/CDummyGeometry.cpp \
  ../src/geometry/elements/CElement.cpp \
  ../src/geometry/elements/CTRIA1.cpp \
  ../src/geometry/elements/CQUAD4.cpp \
  ../src/geometry/elements/CTETRA1.cpp \
  ../src/geometry/elements/CPYRAM5.cpp \
  ../src/geometry/elements/CPRISM6.cpp \
  ../src/geometry/elements/CHEXA8.cpp \
  ../src/geometry/meshreader/CMeshReaderFVM.cpp \
  ../src/geometry/meshreader/CSU2ASCIIMeshReaderFVM.cpp \
  ../src/geometry/meshreader/CCGNSMeshReaderFVM.cpp \
  ../src/geometry/meshreader/CRectangularMeshReaderFVM.cpp \
  ../src/geometry/meshreader/CBoxMeshReaderFVM.cpp \
  ../src/geometry/dual_grid/CDualGrid.cpp \
  ../src/geometry/dual_grid/CEdge.cpp \
  ../src/geometry/dual_grid/CPoint.cpp \
  ../src/geometry/dual_grid/CVertex.cpp \
  ../src/geometry/dual_grid/CTurboVertex.cpp \
  ../src/geometry/primal_grid/CPrimalGrid.cpp \
  ../src/geometry/primal_grid/CPrimalGridFEM.cpp \
  ../src/geometry/primal_grid/CPrimalGridBoundFEM.cpp \
  ../src/geometry/primal_grid/CLine.cpp \
  ../src/geometry/primal_grid/CTriangle.cpp \
  ../src/geometry/primal_grid/CPrism.cpp \
  ../src/geometry/primal_grid/CPyramid.cpp \
  ../src/geometry/primal_grid/CHexahedron.cpp \
  ../src/geometry/primal_grid/CTetrahedron.cpp \
  ../src/geometry/primal_grid/CQuadrilateral.cpp \
  ../src/geometry/primal_grid/CVertexMPI.cpp \
  ../src/interface_interpolation/CInterpolatorFactory.cpp \
  ../src/interface_interpolation/CInterpolator.cpp \
  ../src/interface_interpolation/CMirror.cpp \
  ../src/interface_interpolation/CSlidingMesh.cpp \
  ../src/interface_interpolation/CIsoparametric.cpp \
  ../src/interface_interpolation/CNearestNeighbor.cpp \
  ../src/interface_interpolation/CRadialBasisFunction.cpp \
  ../src/adt/CADTBaseClass.cpp \
  ../src/adt/CADTPointsOnlyClass.cpp \
  ../src/adt/CADTElemClass.cpp \
  ../src/wall_model.cpp \
  ../src/toolboxes/printing_toolbox.cpp \
  ../src/toolboxes/CLinearPartitioner.cpp \
  ../src/toolboxes/C1DInterpolation.cpp \
  ../src/toolboxes/CSymmetricMatrix.cpp \
  ../src/toolboxes/CSquareMatrixCM.cpp \
  ../src/toolboxes/MMS/CVerificationSolution.cpp \
  ../src/toolboxes/MMS/CIncTGVSolution.cpp \
  ../src/toolboxes/MMS/CInviscidVortexSolution.cpp \
  ../src/toolboxes/MMS/CMMSIncEulerSolution.cpp \
  ../src/toolboxes/MMS/CMMSIncNSSolution.cpp \
  ../src/toolboxes/MMS/CMMSNSTwoHalfCirclesSolution.cpp \
  ../src/toolboxes/MMS/CMMSNSTwoHalfSpheresSolution.cpp \
  ../src/toolboxes/MMS/CMMSNSUnitQuadSolution.cpp \
  ../src/toolboxes/MMS/CMMSNSUnitQuadSolutionWallBC.cpp \
  ../src/toolboxes/MMS/CNSUnitQuadSolution.cpp \
  ../src/toolboxes/MMS/CRinglebSolution.cpp \
  ../src/toolboxes/MMS/CTGVSolution.cpp \
  ../src/toolboxes/MMS/CUserDefinedSolution.cpp \
  ../src/linear_algebra/CSysVector.cpp \
  ../src/linear_algebra/CSysMatrix.cpp \
  ../src/linear_algebra/CSysSolve.cpp \
  ../src/linear_algebra/CSysSolve_b.cpp \
  ../src/linear_algebra/CPastixWrapper.cpp

lib_cxxflags = -fPIC -std=c++11
lib_ldadd =

# always link to built dependencies from ./externals
lib_cxxflags += @su2_externals_INCLUDES@
lib_ldadd += @su2_externals_LIBS@

# if BUILD_MUTATIONPP
lib_cxxflags += @MUTATIONPP_CXX@
lib_ldadd += @MUTATIONPP_LD@
# endif

if BUILD_NORMAL
libSU2_a_SOURCES = $(lib_sources)
libSU2_a_CXXFLAGS = ${lib_cxxflags}
libSU2_a_LIBADD = ${lib_ldadd}
endif

if BUILD_DIRECTDIFF
libSU2_DIRECTDIFF_a_SOURCES = $(lib_sources)
libSU2_DIRECTDIFF_a_CXXFLAGS = @DIRECTDIFF_CXX@ ${lib_cxxflags}
libSU2_DIRECTDIFF_a_LIBADD =  ${lib_ldadd}
endif

if BUILD_REVERSE
libSU2_AD_a_SOURCES = $(lib_sources)
libSU2_AD_a_CXXFLAGS = @REVERSE_CXX@ ${lib_cxxflags}
libSU2_AD_a_LIBADD =  ${lib_ldadd}
endif

