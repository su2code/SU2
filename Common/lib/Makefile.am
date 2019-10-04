################################################################################
#
# \file Makefile.am
# \brief Makefile for the SU2 common library
# \author M. Colonno, T. Economon, F. Palacios
# \version 6.2.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
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
  ../include/datatypes/primitive_structure.hpp \
  ../include/datatypes/primitive_structure.inl \
  ../include/datatypes/complex_structure.hpp \
  ../include/datatypes/complex_structure.inl \
  ../include/ad_structure.hpp \
  ../include/ad_structure.inl \
  ../include/fem_cgns_elements.hpp \
  ../include/config_structure.hpp \
  ../include/config_structure.hpp \
  ../include/config_structure.inl \
  ../include/blas_structure.hpp \
  ../include/dual_grid_structure.hpp \
  ../include/dual_grid_structure.inl \
  ../include/fem_geometry_structure.hpp \
  ../include/fem_geometry_structure.inl \
  ../include/fem_standard_element.hpp \
  ../include/fem_standard_element.inl \
  ../include/geometry_structure.hpp \
  ../include/geometry_structure.inl \
  ../include/CMultiGridQueue.hpp \
  ../include/CMeshReaderFVM.hpp \
  ../include/CSU2ASCIIMeshReaderFVM.hpp \
  ../include/CCGNSMeshReaderFVM.hpp \
  ../include/CRectangularMeshReaderFVM.hpp \
  ../include/CBoxMeshReaderFVM.hpp \
  ../include/graph_coloring_structure.hpp \
  ../include/grid_adaptation_structure.hpp \
  ../include/grid_adaptation_structure.inl \
  ../include/grid_movement_structure.hpp \
  ../include/grid_movement_structure.inl \
  ../include/option_structure.hpp \
  ../include/primal_grid_structure.hpp \
  ../include/primal_grid_structure.inl \
  ../include/mpi_structure.hpp \
  ../include/mpi_structure.inl \
  ../include/datatype_structure.hpp \
  ../include/datatype_structure.inl \
  ../include/interpolation_structure.hpp \
  ../include/fem_gauss_jacobi_quadrature.hpp \
  ../include/fem_gauss_jacobi_quadrature.inl \
  ../include/gauss_structure.hpp \
  ../include/gauss_structure.inl \
  ../include/element_structure.hpp \
  ../include/element_structure.inl \
  ../include/adt_structure.hpp \
  ../include/adt_structure.inl \
  ../include/wall_model.hpp \
  ../include/wall_model.inl \
  ../include/toolboxes/printing_toolbox.hpp \
  ../include/toolboxes/CLinearPartitioner.hpp \
  ../include/toolboxes/MMS/CVerificationSolution.hpp \
  ../include/toolboxes/MMS/CVerificationSolution.inl \
  ../include/toolboxes/MMS/CIncTGVSolution.hpp \
  ../include/toolboxes/MMS/CInviscidVortexSolution.hpp \
  ../include/toolboxes/MMS/CMMSIncEulerSolution.hpp \
  ../include/toolboxes/MMS/CMMSIncNSSolution.hpp \
  ../include/toolboxes/MMS/CMMSNSTwoHalfCirclesSolution.hpp \
  ../include/toolboxes/MMS/CMMSNSTwoHalfSpheresSolution.hpp \
  ../include/toolboxes/MMS/CMMSNSUnitQuadSolution.hpp \
  ../include/toolboxes/MMS/CMMSNSUnitQuadSolutionWallBC.hpp \
  ../include/toolboxes/MMS/CNSUnitQuadSolution.hpp \
  ../include/toolboxes/MMS/CRinglebSolution.hpp \
  ../include/toolboxes/MMS/CTGVSolution.hpp \
  ../include/toolboxes/MMS/CUserDefinedSolution.hpp \
  ../include/linear_algebra/CSysVector.hpp \
  ../include/linear_algebra/CSysMatrix.hpp \
  ../include/linear_algebra/CSysMatrix.inl \
  ../include/linear_algebra/CMatrixVectorProduct.hpp \
  ../include/linear_algebra/CPreconditioner.inl \
  ../include/linear_algebra/CSysSolve.hpp \
  ../include/linear_algebra/CSysSolve_b.hpp \
  ../src/fem_cgns_elements.cpp \
  ../src/config_structure.cpp \
  ../src/blas_structure.cpp \
  ../src/dual_grid_structure.cpp \
  ../src/fem_geometry_structure.cpp \
  ../src/fem_integration_rules.cpp \
  ../src/fem_standard_element.cpp \
  ../src/fem_wall_distance.cpp \
  ../src/fem_work_estimate_metis.cpp \
  ../src/geometry_structure.cpp \
  ../src/CMultiGridQueue.cpp \
  ../src/CMeshReaderFVM.cpp \
  ../src/CSU2ASCIIMeshReaderFVM.cpp \
  ../src/CCGNSMeshReaderFVM.cpp \
  ../src/CRectangularMeshReaderFVM.cpp \
  ../src/CBoxMeshReaderFVM.cpp \
  ../src/geometry_structure_fem_part.cpp \
  ../src/graph_coloring_structure.cpp \
  ../src/grid_adaptation_structure.cpp \
  ../src/grid_movement_structure.cpp \
  ../src/primal_grid_structure.cpp \
  ../src/mpi_structure.cpp \
  ../src/ad_structure.cpp \
  ../src/fem_gauss_jacobi_quadrature.cpp \
  ../src/gauss_structure.cpp \
  ../src/element_structure.cpp \
  ../src/element_linear.cpp \
  ../src/interpolation_structure.cpp \
  ../src/adt_structure.cpp \
  ../src/wall_model.cpp \
  ../src/toolboxes/printing_toolbox.cpp \
  ../src/toolboxes/signal_processing_toolbox.cpp \
  ../src/toolboxes/CLinearPartitioner.cpp \
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

