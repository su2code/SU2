/*!
 * \file CPoint.cpp
 * \brief Main classes for defining the points of the dual grid
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../../include/geometry/dual_grid/CPoint.hpp"
#include "../../../include/CConfig.hpp"
#include "../../../include/parallelization/omp_structure.hpp"

CPoint::CPoint(unsigned long npoint, unsigned long ndim) : nDim(ndim) { MinimalAllocation(npoint); }

void CPoint::MinimalAllocation(unsigned long npoint) {
  /*--- Global index a parallel simulation. ---*/
  GlobalIndex.resize(npoint) = 0;

  /*--- Set the color for mesh partitioning. ---*/
  Color.resize(npoint) = 0;

  /*--- Coordinates. ---*/
  Coord.resize(npoint, nDim) = su2double(0.0);
}

CPoint::CPoint(unsigned long npoint, unsigned long ndim, unsigned short imesh, const CConfig* config) : nDim(ndim) {
  MinimalAllocation(npoint);

  FullAllocation(imesh, config);
}

void CPoint::FullAllocation(unsigned short imesh, const CConfig* config) {
  const auto npoint = GlobalIndex.size();

  /*--- Volumes ---*/

  Volume.resize(npoint) = su2double(0.0);
  Periodic_Volume.resize(npoint) = su2double(0.0);

  if (config->GetTime_Marching() != TIME_MARCHING::STEADY) {
    Volume_n.resize(npoint) = su2double(0.0);
    Volume_nM1.resize(npoint) = su2double(0.0);
    if (config->GetDynamic_Grid() && config->GetDiscrete_Adjoint()) {
      Volume_Old.resize(npoint) = su2double(0.0);
      Volume_n_Old.resize(npoint) = su2double(0.0);
      Volume_nM1_Old.resize(npoint) = su2double(0.0);
    }
  }

  if (config->GetDiscrete_Adjoint()) {
    AD_InputIndex.resize(npoint, nDim) = 0;
    AD_OutputIndex.resize(npoint, nDim) = 0;
  }

  /*--- Multigrid structures. ---*/
  if (config->GetnMGLevels() > 0) {
    Parent_CV.resize(npoint) = 0;
    Agglomerate.resize(npoint) = false;
    Agglomerate_Indirect.resize(npoint) = false;
    /*--- The finest grid does not have children CV's. ---*/
    if (imesh != MESH_0) {
      nChildren_CV.resize(npoint) = 0;
      Children_CV.resize(npoint);
    }
  }

  /*--- Identify boundaries, physical boundaries (not send-receive condition), detect if
   *    an element belong to the domain or it must be computed with other processor. ---*/
  Domain.resize(npoint) = true;
  Boundary.resize(npoint) = false;
  SolidBoundary.resize(npoint) = false;
  ViscousBoundary.resize(npoint) = false;
  PhysicalBoundary.resize(npoint) = false;
  PeriodicBoundary.resize(npoint) = false;

  Vertex.resize(npoint);

  /*--- For smoothing the numerical grid coordinates ---*/
  if (config->GetSmoothNumGrid()) {
    Coord_Old.resize(npoint, nDim) = su2double(0.0);
    Coord_Sum.resize(npoint, nDim) = su2double(0.0);
  }

  /*--- Storage of grid velocities for dynamic meshes. ---*/

  if (config->GetDynamic_Grid()) {
    GridVel.resize(npoint, nDim) = su2double(0.0);

    /*--- Grid velocity gradients are needed for the continuous adjoint. ---*/
    if (config->GetContinuous_Adjoint()) GridVel_Grad.resize(npoint, nDim, nDim, 0.0);

    /*--- Structures for storing old node coordinates for computing grid
     *    velocities via finite differencing with dynamically deforming meshes. ---*/
    /*--- In the case of CMeshSolver, these coordinates are stored as solutions to the mesh problem. ---*/
    if (config->GetGrid_Movement() && (config->GetTime_Marching() != TIME_MARCHING::STEADY)) {
      Coord_n.resize(npoint, nDim) = su2double(0.0);
      Coord_p1.resize(npoint, nDim) = su2double(0.0);
      Coord_n1.resize(npoint, nDim) = su2double(0.0);
      if (Coord_Old.empty()) Coord_Old.resize(npoint, nDim) = su2double(0.0);
    }
  }

  /*--- Other geometric properties of the CV's required by numerical methods. ---*/
  nNeighbor.resize(npoint) = 0;
  MaxLength.resize(npoint) = su2double(0.0);
  Curvature.resize(npoint) = su2double(0.0);

  Wall_Distance.resize(npoint) = su2double(0.0);
  ClosestWall_Rank.resize(npoint) = -1;
  ClosestWall_Zone.resize(npoint) = numeric_limits<unsigned short>::max();
  ClosestWall_Marker = ClosestWall_Zone;
  ClosestWall_Elem.resize(npoint) = numeric_limits<unsigned long>::max();

  RoughnessHeight.resize(npoint) = su2double(0.0);
  SharpEdge_Distance.resize(npoint) = su2double(0.0);
}

void CPoint::SetElems(const vector<vector<long> >& elemsMatrix) { Elem = CCompressedSparsePatternL(elemsMatrix); }

void CPoint::SetPoints(const vector<vector<unsigned long> >& pointsMatrix) {
  Point = CCompressedSparsePatternUL(pointsMatrix);
  Edge = CCompressedSparsePatternL(Point.outerPtr(), Point.outerPtr() + Point.getOuterSize() + 1, long(-1));
}

void CPoint::SetVolume_n() {
  assert(Volume_n.size() == Volume.size());
  parallelCopy(Volume.size(), Volume.data(), Volume_n.data());
}

void CPoint::SetVolume_nM1() {
  assert(Volume_nM1.size() == Volume_n.size());
  parallelCopy(Volume_n.size(), Volume_n.data(), Volume_nM1.data());
}

void CPoint::SetVolume_Old() {
  assert(Volume_Old.size() == Volume.size());
  parallelCopy(Volume.size(), Volume.data(), Volume_Old.data());
}

void CPoint::SetVolume_n_Old() {
  assert(Volume_n_Old.size() == Volume_n.size());
  parallelCopy(Volume_n.size(), Volume_n.data(), Volume_n_Old.data());
}

void CPoint::SetVolume_nM1_Old() {
  assert(Volume_nM1_Old.size() == Volume_nM1.size());
  parallelCopy(Volume_nM1.size(), Volume_nM1.data(), Volume_nM1_Old.data());
}

void CPoint::SetVolume_n_from_OldnM1() {
  assert(Volume_n.size() == Volume_nM1_Old.size());
  parallelCopy(Volume_nM1_Old.size(), Volume_nM1_Old.data(), Volume_n.data());
}

void CPoint::SetVolume_from_Oldn() {
  assert(Volume.size() == Volume_n_Old.size());
  parallelCopy(Volume_n_Old.size(), Volume_n_Old.data(), Volume.data());
}

void CPoint::SetCoord_n() {
  assert(Coord_n.size() == Coord.size());
  parallelCopy(Coord.size(), Coord.data(), Coord_n.data());
}

void CPoint::SetCoord_n1() {
  assert(Coord_n1.size() == Coord_n.size());
  parallelCopy(Coord_n.size(), Coord_n.data(), Coord_n1.data());
}

void CPoint::SetCoord_Old() {
  assert(Coord_Old.size() == Coord.size());
  parallelCopy(Coord.size(), Coord.data(), Coord_Old.data());
}

void CPoint::SetCoord_SumZero() { parallelSet(Coord_Sum.size(), 0.0, Coord_Sum.data()); }
