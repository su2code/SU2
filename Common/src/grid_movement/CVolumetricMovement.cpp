/*!
 * \file CVolumetricMovement.cpp
 * \brief Subroutines for moving mesh volume elements
 * \author F. Palacios, T. Economon, S. Padron
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/grid_movement/CVolumetricMovement.hpp"

CVolumetricMovement::CVolumetricMovement() : CGridMovement() {}

CVolumetricMovement::CVolumetricMovement(CGeometry* geometry) : CGridMovement() {
  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();

  nDim = geometry->GetnDim();
}

CVolumetricMovement::~CVolumetricMovement() = default;

void CVolumetricMovement::UpdateDualGrid(CGeometry* geometry, CConfig* config) {
  /*--- After moving all nodes, update the dual mesh. Recompute the edges and
dual mesh control volumes in the domain and on the boundaries. ---*/

  geometry->SetControlVolume(config, UPDATE);
  geometry->SetBoundControlVolume(config, UPDATE);
  geometry->SetMaxLength(config);
}

void CVolumetricMovement::UpdateMultiGrid(CGeometry** geometry, CConfig* config) {
  unsigned short iMGfine, iMGlevel, nMGlevel = config->GetnMGLevels();

  /*--- Update the multigrid structure after moving the finest grid,
including computing the grid velocities on the coarser levels. ---*/

  for (iMGlevel = 1; iMGlevel <= nMGlevel; iMGlevel++) {
    iMGfine = iMGlevel - 1;
    geometry[iMGlevel]->SetControlVolume(geometry[iMGfine], UPDATE);
    geometry[iMGlevel]->SetBoundControlVolume(geometry[iMGfine], config, UPDATE);
    geometry[iMGlevel]->SetCoord(geometry[iMGfine]);
    if (config->GetGrid_Movement()) geometry[iMGlevel]->SetRestricted_GridVelocity(geometry[iMGfine]);
  }
}

void CVolumetricMovement::ComputeDeforming_Element_Volume(CGeometry* geometry, su2double& MinVolume,
                                                          su2double& MaxVolume, bool Screen_Output) {
  unsigned long iElem, ElemCounter = 0, PointCorners[8];
  su2double Volume = 0.0, CoordCorners[8][3];
  unsigned short nNodes = 0, iNodes, iDim;
  bool RightVol = true;

  if (rank == MASTER_NODE && Screen_Output) cout << "Computing volumes of the grid elements." << endl;

  MaxVolume = -1E22;
  MinVolume = 1E22;

  /*--- Load up each triangle and tetrahedron to check for negative volumes. ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) nNodes = 3;
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID) nNodes = 5;
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM) nNodes = 6;
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) nNodes = 8;

    for (iNodes = 0; iNodes < nNodes; iNodes++) {
      PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
      for (iDim = 0; iDim < nDim; iDim++) {
        CoordCorners[iNodes][iDim] = geometry->nodes->GetCoord(PointCorners[iNodes], iDim);
      }
    }

    /*--- 2D elements ---*/

    if (nDim == 2) {
      if (nNodes == 3) Volume = GetTriangle_Area(CoordCorners);
      if (nNodes == 4) Volume = GetQuadrilateral_Area(CoordCorners);
    }

    /*--- 3D Elementes ---*/

    if (nDim == 3) {
      if (nNodes == 4) Volume = GetTetra_Volume(CoordCorners);
      if (nNodes == 5) Volume = GetPyram_Volume(CoordCorners);
      if (nNodes == 6) Volume = GetPrism_Volume(CoordCorners);
      if (nNodes == 8) Volume = GetHexa_Volume(CoordCorners);
    }

    RightVol = true;
    if (Volume < 0.0) RightVol = false;

    MaxVolume = max(MaxVolume, Volume);
    MinVolume = min(MinVolume, Volume);
    geometry->elem[iElem]->SetVolume(Volume);

    if (!RightVol) ElemCounter++;
  }

#ifdef HAVE_MPI
  unsigned long ElemCounter_Local = ElemCounter;
  ElemCounter = 0;
  su2double MaxVolume_Local = MaxVolume;
  MaxVolume = 0.0;
  su2double MinVolume_Local = MinVolume;
  MinVolume = 0.0;
  SU2_MPI::Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MaxVolume_Local, &MaxVolume, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MinVolume_Local, &MinVolume, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
#endif

  /*--- Volume from  0 to 1 ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    Volume = geometry->elem[iElem]->GetVolume() / MaxVolume;
    geometry->elem[iElem]->SetVolume(Volume);
  }

  if ((ElemCounter != 0) && (rank == MASTER_NODE) && (Screen_Output))
    cout << "There are " << ElemCounter << " elements with negative volume.\n" << endl;
}

void CVolumetricMovement::ComputenNonconvexElements(CGeometry* geometry, bool Screen_Output) {
  unsigned long iElem;
  unsigned short iDim;
  unsigned long nNonconvexElements = 0;

  /*--- Load up each tetrahedron to check for convex properties. ---*/
  if (nDim == 2) {
    for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      su2double minCrossProduct = 1.e6, maxCrossProduct = -1.e6;

      const auto nNodes = geometry->elem[iElem]->GetnNodes();

      /*--- Get coordinates of corner points ---*/
      unsigned short iNodes;
      unsigned long PointCorners[8];
      const su2double* CoordCorners[8];

      for (iNodes = 0; iNodes < nNodes; iNodes++) {
        PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
        CoordCorners[iNodes] = geometry->nodes->GetCoord(PointCorners[iNodes]);
      }

      /*--- Determine whether element is convex ---*/
      for (iNodes = 0; iNodes < nNodes; iNodes++) {
        /*--- Calculate minimum and maximum angle between edge vectors adjacent to each node ---*/
        su2double edgeVector_i[3], edgeVector_j[3];

        for (iDim = 0; iDim < nDim; iDim++) {
          if (iNodes == 0) {
            edgeVector_i[iDim] = CoordCorners[nNodes - 1][iDim] - CoordCorners[iNodes][iDim];
          } else {
            edgeVector_i[iDim] = CoordCorners[iNodes - 1][iDim] - CoordCorners[iNodes][iDim];
          }

          if (iNodes == nNodes - 1) {
            edgeVector_j[iDim] = CoordCorners[0][iDim] - CoordCorners[iNodes][iDim];
          } else {
            edgeVector_j[iDim] = CoordCorners[iNodes + 1][iDim] - CoordCorners[iNodes][iDim];
          }
        }

        /*--- Calculate cross product of edge vectors ---*/
        su2double crossProduct;
        crossProduct = edgeVector_i[1] * edgeVector_j[0] - edgeVector_i[0] * edgeVector_j[1];

        if (crossProduct < minCrossProduct) minCrossProduct = crossProduct;
        if (crossProduct > maxCrossProduct) maxCrossProduct = crossProduct;
      }

      /*--- Element is nonconvex if cross product of at least one set of adjacent edges is negative ---*/
      if (minCrossProduct < 0 && maxCrossProduct > 0) {
        nNonconvexElements++;
      }
    }
  } else if (rank == MASTER_NODE) {
    cout << "\nWARNING: Convexity is not checked for 3D elements (issue #1171).\n" << endl;
  }

  unsigned long nNonconvexElements_Local = nNonconvexElements;
  nNonconvexElements = 0;
  SU2_MPI::Allreduce(&nNonconvexElements_Local, &nNonconvexElements, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  /*--- Set number of nonconvex elements in geometry ---*/
  geometry->SetnNonconvexElements(nNonconvexElements);
}

su2double CVolumetricMovement::GetTriangle_Area(su2double CoordCorners[8][3]) const {
  unsigned short iDim;
  su2double a[3] = {0.0, 0.0, 0.0}, b[3] = {0.0, 0.0, 0.0};
  su2double *Coord_0, *Coord_1, *Coord_2, Area;

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];

  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = Coord_0[iDim] - Coord_2[iDim];
    b[iDim] = Coord_1[iDim] - Coord_2[iDim];
  }

  Area = 0.5 * fabs(a[0] * b[1] - a[1] * b[0]);

  return Area;
}

su2double CVolumetricMovement::GetQuadrilateral_Area(su2double CoordCorners[8][3]) const {
  unsigned short iDim;
  su2double a[3] = {0.0, 0.0, 0.0}, b[3] = {0.0, 0.0, 0.0};
  su2double *Coord_0, *Coord_1, *Coord_2, Area;

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];

  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = Coord_0[iDim] - Coord_2[iDim];
    b[iDim] = Coord_1[iDim] - Coord_2[iDim];
  }

  Area = 0.5 * fabs(a[0] * b[1] - a[1] * b[0]);

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[3];

  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = Coord_0[iDim] - Coord_2[iDim];
    b[iDim] = Coord_1[iDim] - Coord_2[iDim];
  }

  Area += 0.5 * fabs(a[0] * b[1] - a[1] * b[0]);

  return Area;
}

su2double CVolumetricMovement::GetTetra_Volume(su2double CoordCorners[8][3]) const {
  unsigned short iDim;
  su2double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
  su2double r1[3] = {0.0, 0.0, 0.0}, r2[3] = {0.0, 0.0, 0.0}, r3[3] = {0.0, 0.0, 0.0},
            CrossProduct[3] = {0.0, 0.0, 0.0}, Volume;

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  Coord_3 = CoordCorners[3];

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  return Volume;
}

su2double CVolumetricMovement::GetPyram_Volume(su2double CoordCorners[8][3]) const {
  unsigned short iDim;
  su2double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
  su2double r1[3] = {0.0, 0.0, 0.0}, r2[3] = {0.0, 0.0, 0.0}, r3[3] = {0.0, 0.0, 0.0},
            CrossProduct[3] = {0.0, 0.0, 0.0}, Volume;

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  Coord_3 = CoordCorners[4];

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[3];
  Coord_3 = CoordCorners[4];

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  return Volume;
}

su2double CVolumetricMovement::GetPrism_Volume(su2double CoordCorners[8][3]) const {
  unsigned short iDim;
  su2double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
  su2double r1[3] = {0.0, 0.0, 0.0}, r2[3] = {0.0, 0.0, 0.0}, r3[3] = {0.0, 0.0, 0.0},
            CrossProduct[3] = {0.0, 0.0, 0.0}, Volume;

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[1];
  Coord_3 = CoordCorners[5];

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[5];
  Coord_2 = CoordCorners[1];
  Coord_3 = CoordCorners[4];

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[5];
  Coord_2 = CoordCorners[4];
  Coord_3 = CoordCorners[3];

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  return Volume;
}

su2double CVolumetricMovement::GetHexa_Volume(su2double CoordCorners[8][3]) const {
  unsigned short iDim;
  su2double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
  su2double r1[3] = {0.0, 0.0, 0.0}, r2[3] = {0.0, 0.0, 0.0}, r3[3] = {0.0, 0.0, 0.0},
            CrossProduct[3] = {0.0, 0.0, 0.0}, Volume;

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  Coord_3 = CoordCorners[5];

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[7];
  Coord_3 = CoordCorners[5];

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[3];
  Coord_3 = CoordCorners[7];

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[5];
  Coord_2 = CoordCorners[7];
  Coord_3 = CoordCorners[4];

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  Coord_0 = CoordCorners[2];
  Coord_1 = CoordCorners[7];
  Coord_2 = CoordCorners[5];
  Coord_3 = CoordCorners[6];

  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  return Volume;
}

void CVolumetricMovement::Rigid_Rotation(CGeometry* geometry, CConfig* config, unsigned short iZone,
                                         unsigned long iter) {
  /*--- Local variables ---*/
  unsigned short iDim, nDim;
  unsigned long iPoint;
  su2double r[3] = {0.0, 0.0, 0.0}, rotCoord[3] = {0.0, 0.0, 0.0}, *Coord;
  su2double Center[3] = {0.0, 0.0, 0.0}, Omega[3] = {0.0, 0.0, 0.0}, Lref;
  su2double dt, Center_Moment[3] = {0.0, 0.0, 0.0};
  su2double *GridVel, newGridVel[3] = {0.0, 0.0, 0.0};
  su2double rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  su2double dtheta, dphi, dpsi, cosTheta, sinTheta;
  su2double cosPhi, sinPhi, cosPsi, sinPsi;
  bool harmonic_balance = (config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE);
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());

  /*--- Problem dimension and physical time step ---*/
  nDim = geometry->GetnDim();
  dt = config->GetDelta_UnstTimeND();
  Lref = config->GetLength_Ref();

  /*--- For the unsteady adjoint, use reverse time ---*/
  if (adjoint) {
    /*--- Set the first adjoint mesh position to the final direct one ---*/
    if (iter == 0) dt = ((su2double)config->GetnTime_Iter() - 1) * dt;
    /*--- Reverse the rotation direction for the adjoint ---*/
    else
      dt = -1.0 * dt;
  } else {
    /*--- No rotation at all for the first direct solution ---*/
    if (iter == 0) dt = 0;
  }

  /*--- Center of rotation & angular velocity vector from config ---*/

  for (iDim = 0; iDim < 3; iDim++) {
    Center[iDim] = config->GetMotion_Origin(iDim);
    Omega[iDim] = config->GetRotation_Rate(iDim) / config->GetOmega_Ref();
  }

  /*-- Set dt for harmonic balance cases ---*/
  if (harmonic_balance) {
    /*--- period of oscillation & compute time interval using nTimeInstances ---*/
    su2double period = config->GetHarmonicBalance_Period();
    period /= config->GetTime_Ref();
    dt = period * (su2double)iter / (su2double)(config->GetnTimeInstances());
  }

  /*--- Compute delta change in the angle about the x, y, & z axes. ---*/

  dtheta = Omega[0] * dt;
  dphi = Omega[1] * dt;
  dpsi = Omega[2] * dt;

  if (rank == MASTER_NODE && iter == 0) {
    cout << " Angular velocity: (" << Omega[0] << ", " << Omega[1];
    cout << ", " << Omega[2] << ") rad/s." << endl;
  }

  /*--- Store angles separately for clarity. Compute sines/cosines. ---*/

  cosTheta = cos(dtheta);
  cosPhi = cos(dphi);
  cosPsi = cos(dpsi);
  sinTheta = sin(dtheta);
  sinPhi = sin(dphi);
  sinPsi = sin(dpsi);

  /*--- Compute the rotation matrix. Note that the implicit
   ordering is rotation about the x-axis, y-axis, then z-axis. ---*/

  rotMatrix[0][0] = cosPhi * cosPsi;
  rotMatrix[1][0] = cosPhi * sinPsi;
  rotMatrix[2][0] = -sinPhi;

  rotMatrix[0][1] = sinTheta * sinPhi * cosPsi - cosTheta * sinPsi;
  rotMatrix[1][1] = sinTheta * sinPhi * sinPsi + cosTheta * cosPsi;
  rotMatrix[2][1] = sinTheta * cosPhi;

  rotMatrix[0][2] = cosTheta * sinPhi * cosPsi + sinTheta * sinPsi;
  rotMatrix[1][2] = cosTheta * sinPhi * sinPsi - sinTheta * cosPsi;
  rotMatrix[2][2] = cosTheta * cosPhi;

  /*--- Loop over and rotate each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    /*--- Coordinates of the current point ---*/
    Coord = geometry->nodes->GetCoord(iPoint);
    GridVel = geometry->nodes->GetGridVel(iPoint);

    /*--- Calculate non-dim. position from rotation center ---*/
    r[0] = (Coord[0] - Center[0]) / Lref;
    r[1] = (Coord[1] - Center[1]) / Lref;
    if (nDim == 3) r[2] = (Coord[2] - Center[2]) / Lref;

    /*--- Compute transformed point coordinates ---*/
    rotCoord[0] = rotMatrix[0][0] * r[0] + rotMatrix[0][1] * r[1] + rotMatrix[0][2] * r[2];

    rotCoord[1] = rotMatrix[1][0] * r[0] + rotMatrix[1][1] * r[1] + rotMatrix[1][2] * r[2];

    rotCoord[2] = rotMatrix[2][0] * r[0] + rotMatrix[2][1] * r[1] + rotMatrix[2][2] * r[2];

    /*--- Cross Product of angular velocity and distance from center.
     Note that we have assumed the grid velocities have been set to
     an initial value in the plunging routine. ---*/

    newGridVel[0] = GridVel[0] + Omega[1] * rotCoord[2] - Omega[2] * rotCoord[1];
    newGridVel[1] = GridVel[1] + Omega[2] * rotCoord[0] - Omega[0] * rotCoord[2];
    if (nDim == 3) newGridVel[2] = GridVel[2] + Omega[0] * rotCoord[1] - Omega[1] * rotCoord[0];

    /*--- Store new node location & grid velocity. Add center.
     Do not store the grid velocity if this is an adjoint calculation.---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->nodes->SetCoord(iPoint, iDim, rotCoord[iDim] + Center[iDim]);
      if (!adjoint) geometry->nodes->SetGridVel(iPoint, iDim, newGridVel[iDim]);
    }
  }

  /*--- Set the moment computation center to the new location after
   incrementing the position with the rotation. ---*/

  for (unsigned short jMarker = 0; jMarker < config->GetnMarker_Monitoring(); jMarker++) {
    Center_Moment[0] = config->GetRefOriginMoment_X(jMarker);
    Center_Moment[1] = config->GetRefOriginMoment_Y(jMarker);
    Center_Moment[2] = config->GetRefOriginMoment_Z(jMarker);

    /*--- Calculate non-dim. position from rotation center ---*/

    for (iDim = 0; iDim < nDim; iDim++) r[iDim] = (Center_Moment[iDim] - Center[iDim]) / Lref;
    if (nDim == 2) r[nDim] = 0.0;

    /*--- Compute transformed point coordinates ---*/

    rotCoord[0] = rotMatrix[0][0] * r[0] + rotMatrix[0][1] * r[1] + rotMatrix[0][2] * r[2];

    rotCoord[1] = rotMatrix[1][0] * r[0] + rotMatrix[1][1] * r[1] + rotMatrix[1][2] * r[2];

    rotCoord[2] = rotMatrix[2][0] * r[0] + rotMatrix[2][1] * r[1] + rotMatrix[2][2] * r[2];

    config->SetRefOriginMoment_X(jMarker, Center[0] + rotCoord[0]);
    config->SetRefOriginMoment_Y(jMarker, Center[1] + rotCoord[1]);
    config->SetRefOriginMoment_Z(jMarker, Center[2] + rotCoord[2]);
  }

  /*--- After moving all nodes, update geometry class ---*/

  UpdateDualGrid(geometry, config);
}

void CVolumetricMovement::Rigid_Pitching(CGeometry* geometry, CConfig* config, unsigned short iZone,
                                         unsigned long iter) {
  /*--- Local variables ---*/
  su2double r[3] = {0.0, 0.0, 0.0}, rotCoord[3] = {0.0, 0.0, 0.0}, *Coord, Center[3] = {0.0, 0.0, 0.0},
            Omega[3] = {0.0, 0.0, 0.0}, Ampl[3] = {0.0, 0.0, 0.0}, Phase[3] = {0.0, 0.0, 0.0};
  su2double Lref, deltaT, alphaDot[3], *GridVel, newGridVel[3] = {0.0, 0.0, 0.0};
  su2double rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  su2double dtheta, dphi, dpsi, cosTheta, sinTheta;
  su2double cosPhi, sinPhi, cosPsi, sinPsi;
  su2double time_new, time_old;
  su2double DEG2RAD = PI_NUMBER / 180.0;
  unsigned short iDim;
  unsigned short nDim = geometry->GetnDim();
  unsigned long iPoint;
  bool harmonic_balance = (config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE);
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());

  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND();
  Lref = config->GetLength_Ref();

  /*--- Pitching origin, frequency, and amplitude from config. ---*/

  for (iDim = 0; iDim < 3; iDim++) {
    Center[iDim] = config->GetMotion_Origin(iDim);
    Omega[iDim] = config->GetPitching_Omega(iDim) / config->GetOmega_Ref();
    Ampl[iDim] = config->GetPitching_Ampl(iDim) * DEG2RAD;
    Phase[iDim] = config->GetPitching_Phase(iDim) * DEG2RAD;
  }

  if (harmonic_balance) {
    /*--- period of oscillation & compute time interval using nTimeInstances ---*/
    su2double period = config->GetHarmonicBalance_Period();
    period /= config->GetTime_Ref();
    deltaT = period / (su2double)(config->GetnTimeInstances());
  }

  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/
    unsigned long nFlowIter = config->GetnTime_Iter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<su2double>(directIter) * deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<su2double>(directIter) + 1.0) * deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<su2double>(iter) * deltaT;
    if (harmonic_balance) {
      /*--- For harmonic balance, begin movement from the zero position ---*/
      time_old = 0.0;
    } else {
      time_old = time_new;
      if (iter != 0) time_old = (static_cast<su2double>(iter) - 1.0) * deltaT;
    }
  }

  /*--- Compute delta change in the angle about the x, y, & z axes. ---*/

  dtheta = -Ampl[0] * (sin(Omega[0] * time_new + Phase[0]) - sin(Omega[0] * time_old + Phase[0]));
  dphi = -Ampl[1] * (sin(Omega[1] * time_new + Phase[1]) - sin(Omega[1] * time_old + Phase[1]));
  dpsi = -Ampl[2] * (sin(Omega[2] * time_new + Phase[2]) - sin(Omega[2] * time_old + Phase[2]));

  /*--- Angular velocity at the new time ---*/

  alphaDot[0] = -Omega[0] * Ampl[0] * cos(Omega[0] * time_new + Phase[0]);
  alphaDot[1] = -Omega[1] * Ampl[1] * cos(Omega[1] * time_new + Phase[1]);
  alphaDot[2] = -Omega[2] * Ampl[2] * cos(Omega[2] * time_new + Phase[2]);

  if (rank == MASTER_NODE && iter == 0) {
    cout << " Pitching frequency: (" << Omega[0] << ", " << Omega[1];
    cout << ", " << Omega[2] << ") rad/s." << endl;
    cout << " Pitching amplitude: (" << Ampl[0] / DEG2RAD << ", ";
    cout << Ampl[1] / DEG2RAD << ", " << Ampl[2] / DEG2RAD;
    cout << ") degrees." << endl;
    cout << " Pitching phase lag: (" << Phase[0] / DEG2RAD << ", ";
    cout << Phase[1] / DEG2RAD << ", " << Phase[2] / DEG2RAD;
    cout << ") degrees." << endl;
  }

  /*--- Store angles separately for clarity. Compute sines/cosines. ---*/

  cosTheta = cos(dtheta);
  cosPhi = cos(dphi);
  cosPsi = cos(dpsi);
  sinTheta = sin(dtheta);
  sinPhi = sin(dphi);
  sinPsi = sin(dpsi);

  /*--- Compute the rotation matrix. Note that the implicit
   ordering is rotation about the x-axis, y-axis, then z-axis. ---*/

  rotMatrix[0][0] = cosPhi * cosPsi;
  rotMatrix[1][0] = cosPhi * sinPsi;
  rotMatrix[2][0] = -sinPhi;

  rotMatrix[0][1] = sinTheta * sinPhi * cosPsi - cosTheta * sinPsi;
  rotMatrix[1][1] = sinTheta * sinPhi * sinPsi + cosTheta * cosPsi;
  rotMatrix[2][1] = sinTheta * cosPhi;

  rotMatrix[0][2] = cosTheta * sinPhi * cosPsi + sinTheta * sinPsi;
  rotMatrix[1][2] = cosTheta * sinPhi * sinPsi - sinTheta * cosPsi;
  rotMatrix[2][2] = cosTheta * cosPhi;

  /*--- Loop over and rotate each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    /*--- Coordinates of the current point ---*/
    Coord = geometry->nodes->GetCoord(iPoint);
    GridVel = geometry->nodes->GetGridVel(iPoint);

    /*--- Calculate non-dim. position from rotation center ---*/
    for (iDim = 0; iDim < nDim; iDim++) r[iDim] = (Coord[iDim] - Center[iDim]) / Lref;
    if (nDim == 2) r[nDim] = 0.0;

    /*--- Compute transformed point coordinates ---*/
    rotCoord[0] = rotMatrix[0][0] * r[0] + rotMatrix[0][1] * r[1] + rotMatrix[0][2] * r[2];

    rotCoord[1] = rotMatrix[1][0] * r[0] + rotMatrix[1][1] * r[1] + rotMatrix[1][2] * r[2];

    rotCoord[2] = rotMatrix[2][0] * r[0] + rotMatrix[2][1] * r[1] + rotMatrix[2][2] * r[2];

    /*--- Cross Product of angular velocity and distance from center.
     Note that we have assumed the grid velocities have been set to
     an initial value in the plunging routine. ---*/

    newGridVel[0] = GridVel[0] + alphaDot[1] * rotCoord[2] - alphaDot[2] * rotCoord[1];
    newGridVel[1] = GridVel[1] + alphaDot[2] * rotCoord[0] - alphaDot[0] * rotCoord[2];
    if (nDim == 3) newGridVel[2] = GridVel[2] + alphaDot[0] * rotCoord[1] - alphaDot[1] * rotCoord[0];

    /*--- Store new node location & grid velocity. Add center location.
     Do not store the grid velocity if this is an adjoint calculation.---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->nodes->SetCoord(iPoint, iDim, rotCoord[iDim] + Center[iDim]);
      if (!adjoint) geometry->nodes->SetGridVel(iPoint, iDim, newGridVel[iDim]);
    }
  }

  /*--- For pitching we don't update the motion origin and moment reference origin. ---*/

  /*--- After moving all nodes, update geometry class ---*/

  UpdateDualGrid(geometry, config);
}

void CVolumetricMovement::Rigid_Plunging(CGeometry* geometry, CConfig* config, unsigned short iZone,
                                         unsigned long iter) {
  /*--- Local variables ---*/
  su2double deltaX[3], newCoord[3] = {0.0, 0.0, 0.0}, Center[3], *Coord, Omega[3], Ampl[3], Lref;
  su2double *GridVel, newGridVel[3] = {0.0, 0.0, 0.0}, xDot[3];
  su2double deltaT, time_new, time_old;
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint;
  bool harmonic_balance = (config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE);
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());

  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND();
  Lref = config->GetLength_Ref();

  for (iDim = 0; iDim < 3; iDim++) {
    Center[iDim] = config->GetMotion_Origin(iDim);
    Omega[iDim] = config->GetPlunging_Omega(iDim) / config->GetOmega_Ref();
    Ampl[iDim] = config->GetPlunging_Ampl(iDim) / Lref;
  }

  /*--- Plunging frequency and amplitude from config. ---*/

  if (harmonic_balance) {
    /*--- period of oscillation & time interval using nTimeInstances ---*/
    su2double period = config->GetHarmonicBalance_Period();
    period /= config->GetTime_Ref();
    deltaT = period / (su2double)(config->GetnTimeInstances());
  }

  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/
    unsigned long nFlowIter = config->GetnTime_Iter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<su2double>(directIter) * deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<su2double>(directIter) + 1.0) * deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<su2double>(iter) * deltaT;
    if (harmonic_balance) {
      /*--- For harmonic balance, begin movement from the zero position ---*/
      time_old = 0.0;
    } else {
      time_old = time_new;
      if (iter != 0) time_old = (static_cast<su2double>(iter) - 1.0) * deltaT;
    }
  }

  /*--- Compute delta change in the position in the x, y, & z directions. ---*/
  deltaX[0] = -Ampl[0] * (sin(Omega[0] * time_new) - sin(Omega[0] * time_old));
  deltaX[1] = -Ampl[1] * (sin(Omega[1] * time_new) - sin(Omega[1] * time_old));
  deltaX[2] = -Ampl[2] * (sin(Omega[2] * time_new) - sin(Omega[2] * time_old));

  /*--- Compute grid velocity due to plunge in the x, y, & z directions. ---*/
  xDot[0] = -Ampl[0] * Omega[0] * (cos(Omega[0] * time_new));
  xDot[1] = -Ampl[1] * Omega[1] * (cos(Omega[1] * time_new));
  xDot[2] = -Ampl[2] * Omega[2] * (cos(Omega[2] * time_new));

  if (rank == MASTER_NODE && iter == 0) {
    cout << " Plunging frequency: (" << Omega[0] << ", " << Omega[1];
    cout << ", " << Omega[2] << ") rad/s." << endl;
    cout << " Plunging amplitude: (" << Ampl[0] << ", ";
    cout << Ampl[1] << ", " << Ampl[2] << ") m." << endl;
  }

  /*--- Loop over and move each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    /*--- Coordinates of the current point ---*/
    Coord = geometry->nodes->GetCoord(iPoint);
    GridVel = geometry->nodes->GetGridVel(iPoint);

    /*--- Increment the node position using the delta values. ---*/
    for (iDim = 0; iDim < nDim; iDim++) newCoord[iDim] = Coord[iDim] + deltaX[iDim];

    /*--- Cross Product of angular velocity and distance from center.
     Note that we have assumed the grid velocities have been set to
     an initial value in the plunging routine. ---*/

    newGridVel[0] = GridVel[0] + xDot[0];
    newGridVel[1] = GridVel[1] + xDot[1];
    if (nDim == 3) newGridVel[2] = GridVel[2] + xDot[2];

    /*--- Store new node location & grid velocity. Do not store the grid
     velocity if this is an adjoint calculation. ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->nodes->SetCoord(iPoint, iDim, newCoord[iDim]);
      if (!adjoint) geometry->nodes->SetGridVel(iPoint, iDim, newGridVel[iDim]);
    }
  }

  /*--- Set the mesh motion center to the new location after
   incrementing the position with the rigid translation. This
   new location will be used for subsequent pitching/rotation.---*/

  for (iDim = 0; iDim < 3; iDim++) {
    Center[iDim] = config->GetMotion_Origin(iDim) + deltaX[iDim];
  }
  config->SetMotion_Origin(Center);

  /*--- As the body origin may have moved, print it to the console ---*/

  //  if (rank == MASTER_NODE) {
  //    cout << " Body origin: (" << Center[0]+deltaX[0];
  //    cout << ", " << Center[1]+deltaX[1] << ", " << Center[2]+deltaX[2];
  //    cout << ")." << endl;
  //  }

  /*--- Set the moment computation center to the new location after
   incrementing the position with the plunging. ---*/

  for (unsigned short jMarker = 0; jMarker < config->GetnMarker_Monitoring(); jMarker++) {
    Center[0] = config->GetRefOriginMoment_X(jMarker) + deltaX[0];
    Center[1] = config->GetRefOriginMoment_Y(jMarker) + deltaX[1];
    Center[2] = config->GetRefOriginMoment_Z(jMarker) + deltaX[2];
    config->SetRefOriginMoment_X(jMarker, Center[0]);
    config->SetRefOriginMoment_Y(jMarker, Center[1]);
    config->SetRefOriginMoment_Z(jMarker, Center[2]);
  }

  /*--- After moving all nodes, update geometry class ---*/

  UpdateDualGrid(geometry, config);
}

void CVolumetricMovement::Rigid_Translation(CGeometry* geometry, CConfig* config, unsigned short iZone,
                                            unsigned long iter) {
  /*--- Local variables ---*/
  su2double deltaX[3], newCoord[3] = {0.0, 0.0, 0.0}, Center[3], *Coord;
  su2double xDot[3];
  su2double deltaT, time_new, time_old;
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint;
  bool harmonic_balance = (config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE);
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());

  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND();

  /*--- Get motion center and translation rates from config ---*/

  for (iDim = 0; iDim < 3; iDim++) {
    Center[iDim] = config->GetMotion_Origin(iDim);
    xDot[iDim] = config->GetTranslation_Rate(iDim);
  }

  if (harmonic_balance) {
    /*--- period of oscillation & time interval using nTimeInstances ---*/
    su2double period = config->GetHarmonicBalance_Period();
    period /= config->GetTime_Ref();
    deltaT = period / (su2double)(config->GetnTimeInstances());
  }

  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/
    unsigned long nFlowIter = config->GetnTime_Iter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<su2double>(directIter) * deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<su2double>(directIter) + 1.0) * deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<su2double>(iter) * deltaT;
    if (harmonic_balance) {
      /*--- For harmonic balance, begin movement from the zero position ---*/
      time_old = 0.0;
    } else {
      time_old = time_new;
      if (iter != 0) time_old = (static_cast<su2double>(iter) - 1.0) * deltaT;
    }
  }

  /*--- Compute delta change in the position in the x, y, & z directions. ---*/
  deltaX[0] = xDot[0] * (time_new - time_old);
  deltaX[1] = xDot[1] * (time_new - time_old);
  deltaX[2] = xDot[2] * (time_new - time_old);

  if (rank == MASTER_NODE) {
    cout << " New physical time: " << time_new << " seconds." << endl;
    if (iter == 0) {
      cout << " Translational velocity: (" << xDot[0] * config->GetVelocity_Ref() << ", "
           << xDot[1] * config->GetVelocity_Ref();
      cout << ", " << xDot[2] * config->GetVelocity_Ref();
      if (config->GetSystemMeasurements() == SI)
        cout << ") m/s." << endl;
      else
        cout << ") ft/s." << endl;
    }
  }

  /*--- Loop over and move each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    /*--- Coordinates of the current point ---*/
    Coord = geometry->nodes->GetCoord(iPoint);

    /*--- Increment the node position using the delta values. ---*/
    for (iDim = 0; iDim < nDim; iDim++) newCoord[iDim] = Coord[iDim] + deltaX[iDim];

    /*--- Store new node location & grid velocity. Do not store the grid
     velocity if this is an adjoint calculation. ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->nodes->SetCoord(iPoint, iDim, newCoord[iDim]);
      if (!adjoint) geometry->nodes->SetGridVel(iPoint, iDim, xDot[iDim]);
    }
  }

  /*--- Set the mesh motion center to the new location after
   incrementing the position with the rigid translation. This
   new location will be used for subsequent pitching/rotation.---*/

  for (iDim = 0; iDim < 3; iDim++) {
    Center[iDim] = config->GetMotion_Origin(iDim) + deltaX[iDim];
  }
  config->SetMotion_Origin(Center);

  /*--- Set the moment computation center to the new location after
   incrementing the position with the translation. ---*/

  for (unsigned short jMarker = 0; jMarker < config->GetnMarker_Monitoring(); jMarker++) {
    Center[0] = config->GetRefOriginMoment_X(jMarker) + deltaX[0];
    Center[1] = config->GetRefOriginMoment_Y(jMarker) + deltaX[1];
    Center[2] = config->GetRefOriginMoment_Z(jMarker) + deltaX[2];
    config->SetRefOriginMoment_X(jMarker, Center[0]);
    config->SetRefOriginMoment_Y(jMarker, Center[1]);
    config->SetRefOriginMoment_Z(jMarker, Center[2]);
  }

  /*--- After moving all nodes, update geometry class ---*/

  UpdateDualGrid(geometry, config);
}

void CVolumetricMovement::SetVolume_Scaling(CGeometry* geometry, CConfig* config, bool UpdateGeo) {
  unsigned short iDim;
  unsigned long iPoint;
  su2double newCoord[3] = {0.0, 0.0, 0.0}, *Coord;

  /*--- The scaling factor is the only input to this option. Currently,
   the mesh must be scaled the same amount in all three directions. ---*/

  su2double Scale = config->GetDV_Value(0) * config->GetOpt_RelaxFactor();

  if (rank == MASTER_NODE) {
    cout << "Scaling the mesh by a constant factor of " << Scale << "." << endl;
  }

  /*--- Loop over and move each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    /*--- Coordinates of the current point ---*/
    Coord = geometry->nodes->GetCoord(iPoint);

    /*--- Scale the node position by the specified factor. ---*/
    for (iDim = 0; iDim < nDim; iDim++) newCoord[iDim] = Scale * Coord[iDim];

    /*--- Store the new node location. ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->nodes->SetCoord(iPoint, iDim, newCoord[iDim]);
    }
  }

  /*--- After moving all nodes, update geometry class ---*/
  if (UpdateGeo) UpdateDualGrid(geometry, config);
}

void CVolumetricMovement::SetVolume_Translation(CGeometry* geometry, CConfig* config, bool UpdateGeo) {
  unsigned short iDim;
  unsigned long iPoint;
  su2double *Coord, deltaX[3] = {0.0, 0.0, 0.0}, newCoord[3] = {0.0, 0.0, 0.0};
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Get the unit vector and magnitude of displacement. Note that we
   assume this is the first DV entry since it is for mesh translation.
   Create the displacement vector from the magnitude and direction. ---*/

  su2double Ampl = config->GetDV_Value(0) * Scale;
  su2double length = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    deltaX[iDim] = config->GetParamDV(0, iDim);
    length += deltaX[iDim] * deltaX[iDim];
  }
  length = sqrt(length);
  for (iDim = 0; iDim < nDim; iDim++) deltaX[iDim] = Ampl * deltaX[iDim] / length;
  if (rank == MASTER_NODE) {
    cout << "Translational displacement: (" << deltaX[0] << ", ";
    cout << deltaX[1] << ", " << deltaX[2] << ")." << endl;
  }

  /*--- Loop over and move each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    /*--- Coordinates of the current point ---*/
    Coord = geometry->nodes->GetCoord(iPoint);

    /*--- Increment the node position using the delta values. ---*/
    for (iDim = 0; iDim < nDim; iDim++) newCoord[iDim] = Coord[iDim] + deltaX[iDim];

    /*--- Store new node location. ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->nodes->SetCoord(iPoint, iDim, newCoord[iDim]);
    }
  }

  /*--- After moving all nodes, update geometry class ---*/
  if (UpdateGeo) UpdateDualGrid(geometry, config);
}

void CVolumetricMovement::SetVolume_Rotation(CGeometry* geometry, CConfig* config, bool UpdateGeo) {
  unsigned short iDim;
  unsigned long iPoint;
  su2double x, y, z;
  su2double *Coord, deltaX[3] = {0.0, 0.0, 0.0}, newCoord[3] = {0.0, 0.0, 0.0};
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- xyz-coordinates of a point on the line of rotation. */
  su2double a = config->GetParamDV(0, 0);
  su2double b = config->GetParamDV(0, 1);
  su2double c = 0.0;
  if (geometry->GetnDim() == 3) c = config->GetParamDV(0, 2);

  /*--- xyz-coordinate of the line's direction vector. ---*/
  su2double u = config->GetParamDV(0, 3) - config->GetParamDV(0, 0);
  su2double v = config->GetParamDV(0, 4) - config->GetParamDV(0, 1);
  su2double w = 1.0;
  if (geometry->GetnDim() == 3) w = config->GetParamDV(0, 5) - config->GetParamDV(0, 2);

  /*--- The angle of rotation. ---*/
  su2double theta = config->GetDV_Value(0) * Scale * PI_NUMBER / 180.0;

  /*--- Print to the console. ---*/
  if (rank == MASTER_NODE) {
    cout << "Rotation axis vector: (" << u << ", ";
    cout << v << ", " << w << ")." << endl;
    cout << "Angle of rotation: " << config->GetDV_Value(0) * Scale;
    cout << " degrees." << endl;
  }

  /*--- Intermediate values used in computations. ---*/
  su2double u2 = u * u;
  su2double v2 = v * v;
  su2double w2 = w * w;
  su2double cosT = cos(theta);
  su2double sinT = sin(theta);
  su2double l2 = u2 + v2 + w2;
  su2double l = sqrt(l2);

  /*--- Loop over and move each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    /*--- Coordinates of the current point ---*/
    Coord = geometry->nodes->GetCoord(iPoint);

    /*--- Displacement for this point due to the rotation. ---*/
    x = Coord[0];
    y = Coord[1];
    z = 0.0;
    if (geometry->GetnDim() == 3) z = Coord[2];

    deltaX[0] = a * (v2 + w2) + u * (-b * v - c * w + u * x + v * y + w * z) +
                (-a * (v2 + w2) + u * (b * v + c * w - v * y - w * z) + (v2 + w2) * x) * cosT +
                l * (-c * v + b * w - w * y + v * z) * sinT;
    deltaX[0] = deltaX[0] / l2 - x;

    deltaX[1] = b * (u2 + w2) + v * (-a * u - c * w + u * x + v * y + w * z) +
                (-b * (u2 + w2) + v * (a * u + c * w - u * x - w * z) + (u2 + w2) * y) * cosT +
                l * (c * u - a * w + w * x - u * z) * sinT;
    deltaX[1] = deltaX[1] / l2 - y;

    deltaX[2] = c * (u2 + v2) + w * (-a * u - b * v + u * x + v * y + w * z) +
                (-c * (u2 + v2) + w * (a * u + b * v - u * x - v * y) + (u2 + v2) * z) * cosT +
                l * (-b * u + a * v - v * x + u * y) * sinT;
    if (geometry->GetnDim() == 3)
      deltaX[2] = deltaX[2] / l2 - z;
    else
      deltaX[2] = 0.0;

    /*--- Increment the node position using the delta values. ---*/
    for (iDim = 0; iDim < nDim; iDim++) newCoord[iDim] = Coord[iDim] + deltaX[iDim];

    /*--- Store new node location. ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->nodes->SetCoord(iPoint, iDim, newCoord[iDim]);
    }
  }

  /*--- After moving all nodes, update geometry class ---*/
  if (UpdateGeo) UpdateDualGrid(geometry, config);
}
