/*!
 * \file CVolumetricMovement.cpp
 * \brief Subroutines for moving mesh volume elements
 * \author F. Palacios, T. Economon, S. Padron
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

#include "../../include/grid_movement/CVolumetricMovement.hpp"
#include "../../include/adt/CADTPointsOnlyClass.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"

CVolumetricMovement::CVolumetricMovement() : CGridMovement(), System(LINEAR_SOLVER_MODE::MESH_DEFORM) {}

CVolumetricMovement::CVolumetricMovement(CGeometry* geometry, CConfig* config)
    : CGridMovement(), System(LINEAR_SOLVER_MODE::MESH_DEFORM) {
  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();

  /*--- Initialize the number of spatial dimensions, length of the state
   vector (same as spatial dimensions for grid deformation), and grid nodes. ---*/

  nDim = geometry->GetnDim();
  nVar = geometry->GetnDim();
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  nIterMesh = 0;

  /*--- Initialize matrix, solution, and r.h.s. structures for the linear solver. ---*/
  if (config->GetVolumetric_Movement() || config->GetSmoothGradient()) {
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    StiffMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
  }
}

CVolumetricMovement::~CVolumetricMovement() = default;

void CVolumetricMovement::UpdateGridCoord(CGeometry* geometry, CConfig* config) {
  unsigned short iDim;
  unsigned long iPoint, total_index;
  su2double new_coord;

  /*--- Update the grid coordinates using the solution of the linear system
   after grid deformation (LinSysSol contains the x, y, z displacements). ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iDim = 0; iDim < nDim; iDim++) {
      total_index = iPoint * nDim + iDim;
      new_coord = geometry->nodes->GetCoord(iPoint, iDim) + LinSysSol[total_index];
      if (fabs(new_coord) < EPS * EPS) new_coord = 0.0;
      geometry->nodes->SetCoord(iPoint, iDim, new_coord);
    }

  /*--- LinSysSol contains the non-transformed displacements in the periodic halo cells.
   * Hence we still need a communication of the transformed coordinates, otherwise periodicity
   * is not maintained. ---*/

  geometry->InitiateComms(geometry, config, COORDINATES);
  geometry->CompleteComms(geometry, config, COORDINATES);
}

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
    geometry[iMGlevel]->SetBoundControlVolume(geometry[iMGfine], UPDATE);
    geometry[iMGlevel]->SetCoord(geometry[iMGfine]);
    if (config->GetGrid_Movement()) geometry[iMGlevel]->SetRestricted_GridVelocity(geometry[iMGfine]);
  }
}

void CVolumetricMovement::SetVolume_Deformation(CGeometry* geometry, CConfig* config, bool UpdateGeo, bool Derivative,
                                                bool ForwardProjectionDerivative) {
  unsigned long Tot_Iter = 0;
  su2double MinVolume, MaxVolume;

  /*--- Retrieve number or iterations, tol, output, etc. from config ---*/

  auto Screen_Output = config->GetDeform_Output();
  auto Nonlinear_Iter = config->GetGridDef_Nonlinear_Iter();

  /*--- Disable the screen output if we're running SU2_CFD ---*/

  if (config->GetKind_SU2() == SU2_COMPONENT::SU2_CFD && !Derivative) Screen_Output = false;
  if (config->GetSmoothGradient()) Screen_Output = true;

  /*--- Set the number of nonlinear iterations to 1 if Derivative computation is enabled ---*/

  if (Derivative) Nonlinear_Iter = 1;

  /*--- Loop over the total number of grid deformation iterations. The surface
   deformation can be divided into increments to help with stability. In
   particular, the linear elasticity equations hold only for small deformations. ---*/

  for (auto iNonlinear_Iter = 0ul; iNonlinear_Iter < Nonlinear_Iter; iNonlinear_Iter++) {
    /*--- Initialize vector and sparse matrix ---*/

    LinSysSol.SetValZero();
    LinSysRes.SetValZero();
    StiffMatrix.SetValZero();

    /*--- Compute the stiffness matrix entries for all nodes/elements in the
     mesh. FEA uses a finite element method discretization of the linear
     elasticity equations (transfers element stiffnesses to point-to-point). ---*/

    MinVolume = SetFEAMethodContributions_Elem(geometry, config);

    /*--- Set the boundary and volume displacements (as prescribed by the
     design variable perturbations controlling the surface shape)
     as a Dirichlet BC. ---*/

    SetBoundaryDisplacements(geometry, config);

    /*--- Fix the location of any points in the domain, if requested. ---*/

    SetDomainDisplacements(geometry, config);

    /*--- Set the boundary derivatives (overrides the actual displacements) ---*/

    if (Derivative) {
      SetBoundaryDerivatives(geometry, config, ForwardProjectionDerivative);
    }

    /*--- Communicate any prescribed boundary displacements via MPI,
     so that all nodes have the same solution and r.h.s. entries
     across all partitions. ---*/

    CSysMatrixComms::Initiate(LinSysSol, geometry, config);
    CSysMatrixComms::Complete(LinSysSol, geometry, config);

    CSysMatrixComms::Initiate(LinSysRes, geometry, config);
    CSysMatrixComms::Complete(LinSysRes, geometry, config);

    /*--- Definition of the preconditioner matrix vector multiplication, and linear solver ---*/

    /*--- To keep legacy behavior ---*/
    System.SetToleranceType(LinearToleranceType::RELATIVE);

    /*--- If we want no derivatives or the direct derivatives, we solve the system using the
     * normal matrix vector product and preconditioner. For the mesh sensitivities using
     * the discrete adjoint method we solve the system using the transposed matrix. ---*/
    if (!Derivative || ((config->GetKind_SU2() == SU2_COMPONENT::SU2_CFD) && Derivative) ||
        (config->GetSmoothGradient() && ForwardProjectionDerivative)) {
      Tot_Iter = System.Solve(StiffMatrix, LinSysRes, LinSysSol, geometry, config);

    } else if (Derivative && (config->GetKind_SU2() == SU2_COMPONENT::SU2_DOT)) {
      Tot_Iter = System.Solve_b(StiffMatrix, LinSysRes, LinSysSol, geometry, config);
    }
    su2double Residual = System.GetResidual();

    /*--- Update the grid coordinates and cell volumes using the solution
     of the linear system (usol contains the x, y, z displacements). ---*/

    if (!Derivative) {
      UpdateGridCoord(geometry, config);
    } else {
      UpdateGridCoord_Derivatives(geometry, config, ForwardProjectionDerivative);
    }
    if (UpdateGeo) {
      UpdateDualGrid(geometry, config);
    }

    if (!Derivative) {
      /*--- Check for failed deformation (negative volumes). ---*/

      ComputeDeforming_Element_Volume(geometry, MinVolume, MaxVolume, Screen_Output);

      /*--- Calculate amount of nonconvex elements ---*/

      ComputenNonconvexElements(geometry, Screen_Output);
    }

    /*--- Set number of iterations in the mesh update. ---*/

    Set_nIterMesh(Tot_Iter);

    if (rank == MASTER_NODE && Screen_Output) {
      cout << "Non-linear iter.: " << iNonlinear_Iter + 1 << "/" << Nonlinear_Iter << ". Linear iter.: " << Tot_Iter
           << ". ";
      if (nDim == 2)
        cout << "Min. area: " << MinVolume << ". Error: " << Residual << "." << endl;
      else
        cout << "Min. volume: " << MinVolume << ". Error: " << Residual << "." << endl;
    }
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

void CVolumetricMovement::ComputeSolid_Wall_Distance(CGeometry* geometry, CConfig* config, su2double& MinDistance,
                                                     su2double& MaxDistance) const {
  unsigned long nVertex_SolidWall, ii, jj, iVertex, iPoint, pointID;
  unsigned short iMarker, iDim;
  su2double dist, MaxDistance_Local, MinDistance_Local;
  int rankID;

  /*--- Initialize min and max distance ---*/

  MaxDistance = -1E22;
  MinDistance = 1E22;

  /*--- Compute the total number of nodes on no-slip boundaries ---*/

  nVertex_SolidWall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); ++iMarker) {
    if (config->GetSolid_Wall(iMarker)) nVertex_SolidWall += geometry->GetnVertex(iMarker);
  }

  /*--- Allocate the vectors to hold boundary node coordinates
   and its local ID. ---*/

  vector<su2double> Coord_bound(nDim * nVertex_SolidWall);
  vector<unsigned long> PointIDs(nVertex_SolidWall);

  /*--- Retrieve and store the coordinates of the no-slip boundary nodes
   and their local point IDs. ---*/

  ii = 0;
  jj = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); ++iMarker) {
    if (config->GetSolid_Wall(iMarker)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); ++iVertex) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        PointIDs[jj++] = iPoint;
        for (iDim = 0; iDim < nDim; ++iDim) Coord_bound[ii++] = geometry->nodes->GetCoord(iPoint, iDim);
      }
    }
  }

  /*--- Build the ADT of the boundary nodes. ---*/

  CADTPointsOnlyClass WallADT(nDim, nVertex_SolidWall, Coord_bound.data(), PointIDs.data(), true);

  /*--- Loop over all interior mesh nodes and compute the distances to each
   of the no-slip boundary nodes. Store the minimum distance to the wall
   for each interior mesh node. ---*/

  if (WallADT.IsEmpty()) {
    /*--- No solid wall boundary nodes in the entire mesh.
     Set the wall distance to zero for all nodes. ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) geometry->nodes->SetWall_Distance(iPoint, 0.0);
  } else {
    /*--- Solid wall boundary nodes are present. Compute the wall
     distance for all nodes. ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {
      WallADT.DetermineNearestNode(geometry->nodes->GetCoord(iPoint), dist, pointID, rankID);
      geometry->nodes->SetWall_Distance(iPoint, dist);

      MaxDistance = max(MaxDistance, dist);

      /*--- To discard points on the surface we use > EPS ---*/

      if (sqrt(dist) > EPS) MinDistance = min(MinDistance, dist);
    }

    MaxDistance_Local = MaxDistance;
    MaxDistance = 0.0;
    MinDistance_Local = MinDistance;
    MinDistance = 0.0;

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&MaxDistance_Local, &MaxDistance, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&MinDistance_Local, &MinDistance, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
#else
    MaxDistance = MaxDistance_Local;
    MinDistance = MinDistance_Local;
#endif
  }
}

su2double CVolumetricMovement::SetFEAMethodContributions_Elem(CGeometry* geometry, CConfig* config) {
  unsigned short iVar, iDim, nNodes = 0, iNodes, StiffMatrix_nElem = 0;
  unsigned long iElem, PointCorners[8];
  su2double **StiffMatrix_Elem = nullptr, CoordCorners[8][3];
  su2double MinVolume = 0.0, MaxVolume = 0.0, MinDistance = 0.0, MaxDistance = 0.0, ElemVolume = 0.0,
            ElemDistance = 0.0;

  bool Screen_Output = config->GetDeform_Output();

  /*--- Allocate maximum size (quadrilateral and hexahedron) ---*/

  if (nDim == 2)
    StiffMatrix_nElem = 8;
  else
    StiffMatrix_nElem = 24;

  StiffMatrix_Elem = new su2double*[StiffMatrix_nElem];
  for (iVar = 0; iVar < StiffMatrix_nElem; iVar++) StiffMatrix_Elem[iVar] = new su2double[StiffMatrix_nElem];

  /*--- Compute min volume in the entire mesh. ---*/

  ComputeDeforming_Element_Volume(geometry, MinVolume, MaxVolume, Screen_Output);
  if (rank == MASTER_NODE && Screen_Output)
    cout << "Min. volume: " << MinVolume << ", max. volume: " << MaxVolume << "." << endl;

  /*--- Compute the distance to the nearest surface if needed
   as part of the stiffness calculation.. ---*/

  if ((config->GetDeform_Stiffness_Type() == SOLID_WALL_DISTANCE) || (config->GetDeform_Limit() < 1E6)) {
    ComputeSolid_Wall_Distance(geometry, config, MinDistance, MaxDistance);
    if (rank == MASTER_NODE && Screen_Output)
      cout << "Min. distance: " << MinDistance << ", max. distance: " << MaxDistance << "." << endl;
  }

  /*--- Compute contributions from each element by forming the stiffness matrix (FEA) ---*/

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

    /*--- Extract Element volume and distance to compute the stiffness ---*/

    ElemVolume = geometry->elem[iElem]->GetVolume();

    if ((config->GetDeform_Stiffness_Type() == SOLID_WALL_DISTANCE)) {
      ElemDistance = 0.0;
      for (iNodes = 0; iNodes < nNodes; iNodes++)
        ElemDistance += geometry->nodes->GetWall_Distance(PointCorners[iNodes]);
      ElemDistance = ElemDistance / (su2double)nNodes;
    }

    if (nDim == 2)
      SetFEA_StiffMatrix2D(geometry, config, StiffMatrix_Elem, PointCorners, CoordCorners, nNodes, ElemVolume,
                           ElemDistance);
    if (nDim == 3)
      SetFEA_StiffMatrix3D(geometry, config, StiffMatrix_Elem, PointCorners, CoordCorners, nNodes, ElemVolume,
                           ElemDistance);

    AddFEA_StiffMatrix(geometry, StiffMatrix_Elem, PointCorners, nNodes);
  }

  /*--- Deallocate memory and exit ---*/

  for (iVar = 0; iVar < StiffMatrix_nElem; iVar++) delete[] StiffMatrix_Elem[iVar];
  delete[] StiffMatrix_Elem;

  return MinVolume;
}

su2double CVolumetricMovement::ShapeFunc_Triangle(su2double Xi, su2double Eta, su2double CoordCorners[8][3],
                                                  su2double DShapeFunction[8][4]) {
  int i, j, k;
  su2double c0, c1, xsj;
  su2double xs[3][3], ad[3][3];

  /*--- Shape functions ---*/

  DShapeFunction[0][3] = Xi;
  DShapeFunction[1][3] = Eta;
  DShapeFunction[2][3] = 1 - Xi - Eta;

  /*--- dN/d xi, dN/d eta ---*/

  DShapeFunction[0][0] = 1.0;
  DShapeFunction[0][1] = 0.0;
  DShapeFunction[1][0] = 0.0;
  DShapeFunction[1][1] = 1.0;
  DShapeFunction[2][0] = -1.0;
  DShapeFunction[2][1] = -1.0;

  /*--- Jacobian transformation ---*/

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 3; k++) {
        xs[i][j] = xs[i][j] + CoordCorners[k][j] * DShapeFunction[k][i];
      }
    }
  }

  /*--- Adjoint to Jacobian ---*/

  ad[0][0] = xs[1][1];
  ad[0][1] = -xs[0][1];
  ad[1][0] = -xs[1][0];
  ad[1][1] = xs[0][0];

  /*--- Determinant of Jacobian ---*/

  xsj = ad[0][0] * ad[1][1] - ad[0][1] * ad[1][0];

  /*--- Jacobian inverse ---*/

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = ad[i][j] / xsj;
    }
  }

  /*--- Derivatives with repect to global coordinates ---*/

  for (k = 0; k < 3; k++) {
    c0 = xs[0][0] * DShapeFunction[k][0] + xs[0][1] * DShapeFunction[k][1];  // dN/dx
    c1 = xs[1][0] * DShapeFunction[k][0] + xs[1][1] * DShapeFunction[k][1];  // dN/dy
    DShapeFunction[k][0] = c0;                                               // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1;                                               // store dN/dy instead of dN/d eta
  }

  return xsj;
}

su2double CVolumetricMovement::ShapeFunc_Quadrilateral(su2double Xi, su2double Eta, su2double CoordCorners[8][3],
                                                       su2double DShapeFunction[8][4]) {
  int i, j, k;
  su2double c0, c1, xsj;
  su2double xs[3][3], ad[3][3];

  /*--- Shape functions ---*/

  DShapeFunction[0][3] = 0.25 * (1.0 - Xi) * (1.0 - Eta);
  DShapeFunction[1][3] = 0.25 * (1.0 + Xi) * (1.0 - Eta);
  DShapeFunction[2][3] = 0.25 * (1.0 + Xi) * (1.0 + Eta);
  DShapeFunction[3][3] = 0.25 * (1.0 - Xi) * (1.0 + Eta);

  /*--- dN/d xi, dN/d eta ---*/

  DShapeFunction[0][0] = -0.25 * (1.0 - Eta);
  DShapeFunction[0][1] = -0.25 * (1.0 - Xi);
  DShapeFunction[1][0] = 0.25 * (1.0 - Eta);
  DShapeFunction[1][1] = -0.25 * (1.0 + Xi);
  DShapeFunction[2][0] = 0.25 * (1.0 + Eta);
  DShapeFunction[2][1] = 0.25 * (1.0 + Xi);
  DShapeFunction[3][0] = -0.25 * (1.0 + Eta);
  DShapeFunction[3][1] = 0.25 * (1.0 - Xi);

  /*--- Jacobian transformation ---*/

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 4; k++) {
        xs[i][j] = xs[i][j] + CoordCorners[k][j] * DShapeFunction[k][i];
      }
    }
  }

  /*--- Adjoint to Jacobian ---*/

  ad[0][0] = xs[1][1];
  ad[0][1] = -xs[0][1];
  ad[1][0] = -xs[1][0];
  ad[1][1] = xs[0][0];

  /*--- Determinant of Jacobian ---*/

  xsj = ad[0][0] * ad[1][1] - ad[0][1] * ad[1][0];

  /*--- Jacobian inverse ---*/

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = ad[i][j] / xsj;
    }
  }

  /*--- Derivatives with repect to global coordinates ---*/

  for (k = 0; k < 4; k++) {
    c0 = xs[0][0] * DShapeFunction[k][0] + xs[0][1] * DShapeFunction[k][1];  // dN/dx
    c1 = xs[1][0] * DShapeFunction[k][0] + xs[1][1] * DShapeFunction[k][1];  // dN/dy
    DShapeFunction[k][0] = c0;                                               // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1;                                               // store dN/dy instead of dN/d eta
  }

  return xsj;
}

su2double CVolumetricMovement::ShapeFunc_Tetra(su2double Xi, su2double Eta, su2double Zeta,
                                               su2double CoordCorners[8][3], su2double DShapeFunction[8][4]) {
  int i, j, k;
  su2double c0, c1, c2, xsj;
  su2double xs[3][3], ad[3][3];

  /*--- Shape functions ---*/

  DShapeFunction[0][3] = Xi;
  DShapeFunction[1][3] = Zeta;
  DShapeFunction[2][3] = 1.0 - Xi - Eta - Zeta;
  DShapeFunction[3][3] = Eta;

  /*--- dN/d xi, dN/d eta, dN/d zeta ---*/

  DShapeFunction[0][0] = 1.0;
  DShapeFunction[0][1] = 0.0;
  DShapeFunction[0][2] = 0.0;
  DShapeFunction[1][0] = 0.0;
  DShapeFunction[1][1] = 0.0;
  DShapeFunction[1][2] = 1.0;
  DShapeFunction[2][0] = -1.0;
  DShapeFunction[2][1] = -1.0;
  DShapeFunction[2][2] = -1.0;
  DShapeFunction[3][0] = 0.0;
  DShapeFunction[3][1] = 1.0;
  DShapeFunction[3][2] = 0.0;

  /*--- Jacobian transformation ---*/

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 4; k++) {
        xs[i][j] = xs[i][j] + CoordCorners[k][j] * DShapeFunction[k][i];
      }
    }
  }

  /*--- Adjoint to Jacobian ---*/

  ad[0][0] = xs[1][1] * xs[2][2] - xs[1][2] * xs[2][1];
  ad[0][1] = xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2];
  ad[0][2] = xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1];
  ad[1][0] = xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2];
  ad[1][1] = xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0];
  ad[1][2] = xs[0][2] * xs[1][0] - xs[0][0] * xs[1][2];
  ad[2][0] = xs[1][0] * xs[2][1] - xs[1][1] * xs[2][0];
  ad[2][1] = xs[0][1] * xs[2][0] - xs[0][0] * xs[2][1];
  ad[2][2] = xs[0][0] * xs[1][1] - xs[0][1] * xs[1][0];

  /*--- Determinant of Jacobian ---*/

  xsj = xs[0][0] * ad[0][0] + xs[0][1] * ad[1][0] + xs[0][2] * ad[2][0];

  /*--- Jacobian inverse ---*/

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j] / xsj;
    }
  }

  /*--- Derivatives with repect to global coordinates ---*/

  for (k = 0; k < 4; k++) {
    c0 = xs[0][0] * DShapeFunction[k][0] + xs[0][1] * DShapeFunction[k][1] + xs[0][2] * DShapeFunction[k][2];  // dN/dx
    c1 = xs[1][0] * DShapeFunction[k][0] + xs[1][1] * DShapeFunction[k][1] + xs[1][2] * DShapeFunction[k][2];  // dN/dy
    c2 = xs[2][0] * DShapeFunction[k][0] + xs[2][1] * DShapeFunction[k][1] + xs[2][2] * DShapeFunction[k][2];  // dN/dz
    DShapeFunction[k][0] = c0;  // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1;  // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2;  // store dN/dz instead of dN/d zeta
  }

  return xsj;
}

su2double CVolumetricMovement::ShapeFunc_Pyram(su2double Xi, su2double Eta, su2double Zeta,
                                               su2double CoordCorners[8][3], su2double DShapeFunction[8][4]) {
  int i, j, k;
  su2double c0, c1, c2, xsj;
  su2double xs[3][3], ad[3][3];

  /*--- Shape functions ---*/

  DShapeFunction[0][3] = 0.25 * (-Xi + Eta + Zeta - 1.0) * (-Xi - Eta + Zeta - 1.0) / (1.0 - Zeta);
  DShapeFunction[1][3] = 0.25 * (-Xi - Eta + Zeta - 1.0) * (Xi - Eta + Zeta - 1.0) / (1.0 - Zeta);
  DShapeFunction[2][3] = 0.25 * (Xi + Eta + Zeta - 1.0) * (Xi - Eta + Zeta - 1.0) / (1.0 - Zeta);
  DShapeFunction[3][3] = 0.25 * (Xi + Eta + Zeta - 1.0) * (-Xi + Eta + Zeta - 1.0) / (1.0 - Zeta);
  DShapeFunction[4][3] = Zeta;

  /*--- dN/d xi ---*/

  DShapeFunction[0][0] = 0.5 * (Zeta - Xi - 1.0) / (Zeta - 1.0);
  DShapeFunction[1][0] = 0.5 * Xi / (Zeta - 1.0);
  DShapeFunction[2][0] = 0.5 * (1.0 - Zeta - Xi) / (Zeta - 1.0);
  DShapeFunction[3][0] = DShapeFunction[1][0];
  DShapeFunction[4][0] = 0.0;

  /*--- dN/d eta ---*/

  DShapeFunction[0][1] = 0.5 * Eta / (Zeta - 1.0);
  DShapeFunction[1][1] = 0.5 * (Zeta - Eta - 1.0) / (Zeta - 1.0);
  DShapeFunction[2][1] = DShapeFunction[0][1];
  DShapeFunction[3][1] = 0.5 * (1.0 - Zeta - Eta) / (Zeta - 1.0);
  DShapeFunction[4][1] = 0.0;

  /*--- dN/d zeta ---*/

  DShapeFunction[0][2] = 0.25 * (-1.0 + 2.0 * Zeta - Zeta * Zeta - Eta * Eta + Xi * Xi) / ((1.0 - Zeta) * (1.0 - Zeta));
  DShapeFunction[1][2] = 0.25 * (-1.0 + 2.0 * Zeta - Zeta * Zeta + Eta * Eta - Xi * Xi) / ((1.0 - Zeta) * (1.0 - Zeta));
  DShapeFunction[2][2] = DShapeFunction[0][2];
  DShapeFunction[3][2] = DShapeFunction[1][2];
  DShapeFunction[4][2] = 1.0;

  /*--- Jacobian transformation ---*/

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 5; k++) {
        xs[i][j] = xs[i][j] + CoordCorners[k][j] * DShapeFunction[k][i];
      }
    }
  }

  /*--- Adjoint to Jacobian ---*/

  ad[0][0] = xs[1][1] * xs[2][2] - xs[1][2] * xs[2][1];
  ad[0][1] = xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2];
  ad[0][2] = xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1];
  ad[1][0] = xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2];
  ad[1][1] = xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0];
  ad[1][2] = xs[0][2] * xs[1][0] - xs[0][0] * xs[1][2];
  ad[2][0] = xs[1][0] * xs[2][1] - xs[1][1] * xs[2][0];
  ad[2][1] = xs[0][1] * xs[2][0] - xs[0][0] * xs[2][1];
  ad[2][2] = xs[0][0] * xs[1][1] - xs[0][1] * xs[1][0];

  /*--- Determinant of Jacobian ---*/

  xsj = xs[0][0] * ad[0][0] + xs[0][1] * ad[1][0] + xs[0][2] * ad[2][0];

  /*--- Jacobian inverse ---*/

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j] / xsj;
    }
  }

  /*--- Derivatives with repect to global coordinates ---*/

  for (k = 0; k < 5; k++) {
    c0 = xs[0][0] * DShapeFunction[k][0] + xs[0][1] * DShapeFunction[k][1] + xs[0][2] * DShapeFunction[k][2];  // dN/dx
    c1 = xs[1][0] * DShapeFunction[k][0] + xs[1][1] * DShapeFunction[k][1] + xs[1][2] * DShapeFunction[k][2];  // dN/dy
    c2 = xs[2][0] * DShapeFunction[k][0] + xs[2][1] * DShapeFunction[k][1] + xs[2][2] * DShapeFunction[k][2];  // dN/dz
    DShapeFunction[k][0] = c0;  // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1;  // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2;  // store dN/dz instead of dN/d zeta
  }

  return xsj;
}

su2double CVolumetricMovement::ShapeFunc_Prism(su2double Xi, su2double Eta, su2double Zeta,
                                               su2double CoordCorners[8][3], su2double DShapeFunction[8][4]) {
  int i, j, k;
  su2double c0, c1, c2, xsj;
  su2double xs[3][3], ad[3][3];

  /*--- Shape functions ---*/

  DShapeFunction[0][3] = 0.5 * Eta * (1.0 - Xi);
  DShapeFunction[1][3] = 0.5 * Zeta * (1.0 - Xi);
  DShapeFunction[2][3] = 0.5 * (1.0 - Eta - Zeta) * (1.0 - Xi);
  DShapeFunction[3][3] = 0.5 * Eta * (Xi + 1.0);
  DShapeFunction[4][3] = 0.5 * Zeta * (Xi + 1.0);
  DShapeFunction[5][3] = 0.5 * (1.0 - Eta - Zeta) * (Xi + 1.0);

  /*--- dN/d Xi, dN/d Eta, dN/d Zeta ---*/

  DShapeFunction[0][0] = -0.5 * Eta;
  DShapeFunction[0][1] = 0.5 * (1.0 - Xi);
  DShapeFunction[0][2] = 0.0;
  DShapeFunction[1][0] = -0.5 * Zeta;
  DShapeFunction[1][1] = 0.0;
  DShapeFunction[1][2] = 0.5 * (1.0 - Xi);
  DShapeFunction[2][0] = -0.5 * (1.0 - Eta - Zeta);
  DShapeFunction[2][1] = -0.5 * (1.0 - Xi);
  DShapeFunction[2][2] = -0.5 * (1.0 - Xi);
  DShapeFunction[3][0] = 0.5 * Eta;
  DShapeFunction[3][1] = 0.5 * (Xi + 1.0);
  DShapeFunction[3][2] = 0.0;
  DShapeFunction[4][0] = 0.5 * Zeta;
  DShapeFunction[4][1] = 0.0;
  DShapeFunction[4][2] = 0.5 * (Xi + 1.0);
  DShapeFunction[5][0] = 0.5 * (1.0 - Eta - Zeta);
  DShapeFunction[5][1] = -0.5 * (Xi + 1.0);
  DShapeFunction[5][2] = -0.5 * (Xi + 1.0);

  /*--- Jacobian transformation ---*/

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 6; k++) {
        xs[i][j] = xs[i][j] + CoordCorners[k][j] * DShapeFunction[k][i];
      }
    }
  }

  /*--- Adjoint to Jacobian ---*/

  ad[0][0] = xs[1][1] * xs[2][2] - xs[1][2] * xs[2][1];
  ad[0][1] = xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2];
  ad[0][2] = xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1];
  ad[1][0] = xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2];
  ad[1][1] = xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0];
  ad[1][2] = xs[0][2] * xs[1][0] - xs[0][0] * xs[1][2];
  ad[2][0] = xs[1][0] * xs[2][1] - xs[1][1] * xs[2][0];
  ad[2][1] = xs[0][1] * xs[2][0] - xs[0][0] * xs[2][1];
  ad[2][2] = xs[0][0] * xs[1][1] - xs[0][1] * xs[1][0];

  /*--- Determinant of Jacobian ---*/

  xsj = xs[0][0] * ad[0][0] + xs[0][1] * ad[1][0] + xs[0][2] * ad[2][0];

  /*--- Jacobian inverse ---*/

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j] / xsj;
    }
  }

  /*--- Derivatives with repect to global coordinates ---*/

  for (k = 0; k < 6; k++) {
    c0 = xs[0][0] * DShapeFunction[k][0] + xs[0][1] * DShapeFunction[k][1] + xs[0][2] * DShapeFunction[k][2];  // dN/dx
    c1 = xs[1][0] * DShapeFunction[k][0] + xs[1][1] * DShapeFunction[k][1] + xs[1][2] * DShapeFunction[k][2];  // dN/dy
    c2 = xs[2][0] * DShapeFunction[k][0] + xs[2][1] * DShapeFunction[k][1] + xs[2][2] * DShapeFunction[k][2];  // dN/dz
    DShapeFunction[k][0] = c0;  // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1;  // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2;  // store dN/dz instead of dN/d zeta
  }

  return xsj;
}

su2double CVolumetricMovement::ShapeFunc_Hexa(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3],
                                              su2double DShapeFunction[8][4]) {
  int i, j, k;
  su2double c0, c1, c2, xsj;
  su2double xs[3][3], ad[3][3];

  /*--- Shape functions ---*/

  DShapeFunction[0][3] = 0.125 * (1.0 - Xi) * (1.0 - Eta) * (1.0 - Zeta);
  DShapeFunction[1][3] = 0.125 * (1.0 + Xi) * (1.0 - Eta) * (1.0 - Zeta);
  DShapeFunction[2][3] = 0.125 * (1.0 + Xi) * (1.0 + Eta) * (1.0 - Zeta);
  DShapeFunction[3][3] = 0.125 * (1.0 - Xi) * (1.0 + Eta) * (1.0 - Zeta);
  DShapeFunction[4][3] = 0.125 * (1.0 - Xi) * (1.0 - Eta) * (1.0 + Zeta);
  DShapeFunction[5][3] = 0.125 * (1.0 + Xi) * (1.0 - Eta) * (1.0 + Zeta);
  DShapeFunction[6][3] = 0.125 * (1.0 + Xi) * (1.0 + Eta) * (1.0 + Zeta);
  DShapeFunction[7][3] = 0.125 * (1.0 - Xi) * (1.0 + Eta) * (1.0 + Zeta);

  /*--- dN/d xi ---*/

  DShapeFunction[0][0] = -0.125 * (1.0 - Eta) * (1.0 - Zeta);
  DShapeFunction[1][0] = 0.125 * (1.0 - Eta) * (1.0 - Zeta);
  DShapeFunction[2][0] = 0.125 * (1.0 + Eta) * (1.0 - Zeta);
  DShapeFunction[3][0] = -0.125 * (1.0 + Eta) * (1.0 - Zeta);
  DShapeFunction[4][0] = -0.125 * (1.0 - Eta) * (1.0 + Zeta);
  DShapeFunction[5][0] = 0.125 * (1.0 - Eta) * (1.0 + Zeta);
  DShapeFunction[6][0] = 0.125 * (1.0 + Eta) * (1.0 + Zeta);
  DShapeFunction[7][0] = -0.125 * (1.0 + Eta) * (1.0 + Zeta);

  /*--- dN/d eta ---*/

  DShapeFunction[0][1] = -0.125 * (1.0 - Xi) * (1.0 - Zeta);
  DShapeFunction[1][1] = -0.125 * (1.0 + Xi) * (1.0 - Zeta);
  DShapeFunction[2][1] = 0.125 * (1.0 + Xi) * (1.0 - Zeta);
  DShapeFunction[3][1] = 0.125 * (1.0 - Xi) * (1.0 - Zeta);
  DShapeFunction[4][1] = -0.125 * (1.0 - Xi) * (1.0 + Zeta);
  DShapeFunction[5][1] = -0.125 * (1.0 + Xi) * (1.0 + Zeta);
  DShapeFunction[6][1] = 0.125 * (1.0 + Xi) * (1.0 + Zeta);
  DShapeFunction[7][1] = 0.125 * (1.0 - Xi) * (1.0 + Zeta);

  /*--- dN/d zeta ---*/

  DShapeFunction[0][2] = -0.125 * (1.0 - Xi) * (1.0 - Eta);
  DShapeFunction[1][2] = -0.125 * (1.0 + Xi) * (1.0 - Eta);
  DShapeFunction[2][2] = -0.125 * (1.0 + Xi) * (1.0 + Eta);
  DShapeFunction[3][2] = -0.125 * (1.0 - Xi) * (1.0 + Eta);
  DShapeFunction[4][2] = 0.125 * (1.0 - Xi) * (1.0 - Eta);
  DShapeFunction[5][2] = 0.125 * (1.0 + Xi) * (1.0 - Eta);
  DShapeFunction[6][2] = 0.125 * (1.0 + Xi) * (1.0 + Eta);
  DShapeFunction[7][2] = 0.125 * (1.0 - Xi) * (1.0 + Eta);

  /*--- Jacobian transformation ---*/

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 8; k++) {
        xs[i][j] = xs[i][j] + CoordCorners[k][j] * DShapeFunction[k][i];
      }
    }
  }

  /*--- Adjoint to Jacobian ---*/

  ad[0][0] = xs[1][1] * xs[2][2] - xs[1][2] * xs[2][1];
  ad[0][1] = xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2];
  ad[0][2] = xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1];
  ad[1][0] = xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2];
  ad[1][1] = xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0];
  ad[1][2] = xs[0][2] * xs[1][0] - xs[0][0] * xs[1][2];
  ad[2][0] = xs[1][0] * xs[2][1] - xs[1][1] * xs[2][0];
  ad[2][1] = xs[0][1] * xs[2][0] - xs[0][0] * xs[2][1];
  ad[2][2] = xs[0][0] * xs[1][1] - xs[0][1] * xs[1][0];

  /*--- Determinant of Jacobian ---*/

  xsj = xs[0][0] * ad[0][0] + xs[0][1] * ad[1][0] + xs[0][2] * ad[2][0];

  /*--- Jacobian inverse ---*/
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j] / xsj;
    }
  }

  /*--- Derivatives with repect to global coordinates ---*/

  for (k = 0; k < 8; k++) {
    c0 = xs[0][0] * DShapeFunction[k][0] + xs[0][1] * DShapeFunction[k][1] + xs[0][2] * DShapeFunction[k][2];  // dN/dx
    c1 = xs[1][0] * DShapeFunction[k][0] + xs[1][1] * DShapeFunction[k][1] + xs[1][2] * DShapeFunction[k][2];  // dN/dy
    c2 = xs[2][0] * DShapeFunction[k][0] + xs[2][1] * DShapeFunction[k][1] + xs[2][2] * DShapeFunction[k][2];  // dN/dz
    DShapeFunction[k][0] = c0;  // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1;  // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2;  // store dN/dz instead of dN/d zeta
  }

  return xsj;
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

void CVolumetricMovement::SetFEA_StiffMatrix2D(CGeometry* geometry, CConfig* config, su2double** StiffMatrix_Elem,
                                               unsigned long PointCorners[8], su2double CoordCorners[8][3],
                                               unsigned short nNodes, su2double ElemVolume, su2double ElemDistance) {
  su2double B_Matrix[3][8], D_Matrix[3][3], Aux_Matrix[8][3];
  su2double Xi = 0.0, Eta = 0.0, Det = 0.0, E = 1 / EPS, Lambda = 0.0, Mu = 0.0, Nu = 0.0;
  unsigned short iNode, iVar, jVar, kVar, iGauss, nGauss = 0;
  su2double DShapeFunction[8][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0},
                                    {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0},
                                    {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
  su2double Location[4][3], Weight[4];
  unsigned short nVar = geometry->GetnDim();

  for (iVar = 0; iVar < nNodes * nVar; iVar++) {
    for (jVar = 0; jVar < nNodes * nVar; jVar++) {
      StiffMatrix_Elem[iVar][jVar] = 0.0;
    }
  }

  /*--- Integration formulae from "Shape functions and points of
   integration of the Résumé" by Josselin DELMAS (2013) ---*/

  /*--- Triangle. Nodes of numerical integration at 1 point (order 1). ---*/

  if (nNodes == 3) {
    nGauss = 1;
    Location[0][0] = 0.333333333333333;
    Location[0][1] = 0.333333333333333;
    Weight[0] = 0.5;
  }

  /*--- Quadrilateral. Nodes of numerical integration at 4 points (order 2). ---*/

  if (nNodes == 4) {
    nGauss = 4;
    Location[0][0] = -0.577350269189626;
    Location[0][1] = -0.577350269189626;
    Weight[0] = 1.0;
    Location[1][0] = 0.577350269189626;
    Location[1][1] = -0.577350269189626;
    Weight[1] = 1.0;
    Location[2][0] = 0.577350269189626;
    Location[2][1] = 0.577350269189626;
    Weight[2] = 1.0;
    Location[3][0] = -0.577350269189626;
    Location[3][1] = 0.577350269189626;
    Weight[3] = 1.0;
  }

  for (iGauss = 0; iGauss < nGauss; iGauss++) {
    Xi = Location[iGauss][0];
    Eta = Location[iGauss][1];

    if (nNodes == 3) Det = ShapeFunc_Triangle(Xi, Eta, CoordCorners, DShapeFunction);
    if (nNodes == 4) Det = ShapeFunc_Quadrilateral(Xi, Eta, CoordCorners, DShapeFunction);

    /*--- Compute the B Matrix ---*/

    for (iVar = 0; iVar < 3; iVar++)
      for (jVar = 0; jVar < nNodes * nVar; jVar++) B_Matrix[iVar][jVar] = 0.0;

    for (iNode = 0; iNode < nNodes; iNode++) {
      B_Matrix[0][0 + iNode * nVar] = DShapeFunction[iNode][0];
      B_Matrix[1][1 + iNode * nVar] = DShapeFunction[iNode][1];

      B_Matrix[2][0 + iNode * nVar] = DShapeFunction[iNode][1];
      B_Matrix[2][1 + iNode * nVar] = DShapeFunction[iNode][0];
    }

    /*--- Impose a type of stiffness for each element ---*/

    switch (config->GetDeform_Stiffness_Type()) {
      case INVERSE_VOLUME:
        E = 1.0 / ElemVolume;
        break;
      case SOLID_WALL_DISTANCE:
        E = 1.0 / ElemDistance;
        break;
      case CONSTANT_STIFFNESS:
        E = 1.0 / EPS;
        break;
    }

    Nu = config->GetDeform_Coeff();
    Mu = E / (2.0 * (1.0 + Nu));
    Lambda = Nu * E / ((1.0 + Nu) * (1.0 - 2.0 * Nu));

    /*--- Compute the D Matrix (for plane strain and 3-D)---*/

    D_Matrix[0][0] = Lambda + 2.0 * Mu;
    D_Matrix[0][1] = Lambda;
    D_Matrix[0][2] = 0.0;
    D_Matrix[1][0] = Lambda;
    D_Matrix[1][1] = Lambda + 2.0 * Mu;
    D_Matrix[1][2] = 0.0;
    D_Matrix[2][0] = 0.0;
    D_Matrix[2][1] = 0.0;
    D_Matrix[2][2] = Mu;

    /*--- Compute the BT.D Matrix ---*/

    for (iVar = 0; iVar < nNodes * nVar; iVar++) {
      for (jVar = 0; jVar < 3; jVar++) {
        Aux_Matrix[iVar][jVar] = 0.0;
        for (kVar = 0; kVar < 3; kVar++) Aux_Matrix[iVar][jVar] += B_Matrix[kVar][iVar] * D_Matrix[kVar][jVar];
      }
    }

    /*--- Compute the BT.D.B Matrix (stiffness matrix), and add to the original
     matrix using Gauss integration ---*/

    for (iVar = 0; iVar < nNodes * nVar; iVar++) {
      for (jVar = 0; jVar < nNodes * nVar; jVar++) {
        for (kVar = 0; kVar < 3; kVar++) {
          StiffMatrix_Elem[iVar][jVar] += Weight[iGauss] * Aux_Matrix[iVar][kVar] * B_Matrix[kVar][jVar] * fabs(Det);
        }
      }
    }
  }
}

void CVolumetricMovement::SetFEA_StiffMatrix3D(CGeometry* geometry, CConfig* config, su2double** StiffMatrix_Elem,
                                               unsigned long PointCorners[8], su2double CoordCorners[8][3],
                                               unsigned short nNodes, su2double ElemVolume, su2double ElemDistance) {
  su2double B_Matrix[6][24],
      D_Matrix[6][6] = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      Aux_Matrix[24][6];
  su2double Xi = 0.0, Eta = 0.0, Zeta = 0.0, Det = 0.0, Mu = 0.0, E = 0.0, Lambda = 0.0, Nu = 0.0;
  unsigned short iNode, iVar, jVar, kVar, iGauss, nGauss = 0;
  su2double DShapeFunction[8][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0},
                                    {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0},
                                    {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
  su2double Location[8][3], Weight[8];
  unsigned short nVar = geometry->GetnDim();

  for (iVar = 0; iVar < nNodes * nVar; iVar++) {
    for (jVar = 0; jVar < nNodes * nVar; jVar++) {
      StiffMatrix_Elem[iVar][jVar] = 0.0;
    }
  }

  /*--- Integration formulae from "Shape functions and points of
   integration of the Résumé" by Josselin Delmas (2013) ---*/

  /*--- Tetrahedrons. Nodes of numerical integration at 1 point (order 1). ---*/

  if (nNodes == 4) {
    nGauss = 1;
    Location[0][0] = 0.25;
    Location[0][1] = 0.25;
    Location[0][2] = 0.25;
    Weight[0] = 0.166666666666666;
  }

  /*--- Pyramids. Nodes numerical integration at 5 points. ---*/

  if (nNodes == 5) {
    nGauss = 5;
    Location[0][0] = 0.5;
    Location[0][1] = 0.0;
    Location[0][2] = 0.1531754163448146;
    Weight[0] = 0.133333333333333;
    Location[1][0] = 0.0;
    Location[1][1] = 0.5;
    Location[1][2] = 0.1531754163448146;
    Weight[1] = 0.133333333333333;
    Location[2][0] = -0.5;
    Location[2][1] = 0.0;
    Location[2][2] = 0.1531754163448146;
    Weight[2] = 0.133333333333333;
    Location[3][0] = 0.0;
    Location[3][1] = -0.5;
    Location[3][2] = 0.1531754163448146;
    Weight[3] = 0.133333333333333;
    Location[4][0] = 0.0;
    Location[4][1] = 0.0;
    Location[4][2] = 0.6372983346207416;
    Weight[4] = 0.133333333333333;
  }

  /*--- Prism. Nodes of numerical integration at 6 points (order 3 in Xi, order 2 in Eta and Mu ). ---*/

  if (nNodes == 6) {
    nGauss = 6;
    Location[0][0] = -0.577350269189626;
    Location[0][1] = 0.166666666666667;
    Location[0][2] = 0.166666666666667;
    Weight[0] = 0.166666666666667;
    Location[1][0] = -0.577350269189626;
    Location[1][1] = 0.666666666666667;
    Location[1][2] = 0.166666666666667;
    Weight[1] = 0.166666666666667;
    Location[2][0] = -0.577350269189626;
    Location[2][1] = 0.166666666666667;
    Location[2][2] = 0.666666666666667;
    Weight[2] = 0.166666666666667;
    Location[3][0] = 0.577350269189626;
    Location[3][1] = 0.166666666666667;
    Location[3][2] = 0.166666666666667;
    Weight[3] = 0.166666666666667;
    Location[4][0] = 0.577350269189626;
    Location[4][1] = 0.666666666666667;
    Location[4][2] = 0.166666666666667;
    Weight[4] = 0.166666666666667;
    Location[5][0] = 0.577350269189626;
    Location[5][1] = 0.166666666666667;
    Location[5][2] = 0.666666666666667;
    Weight[5] = 0.166666666666667;
  }

  /*--- Hexahedrons. Nodes of numerical integration at 6 points (order 3). ---*/

  if (nNodes == 8) {
    nGauss = 8;
    Location[0][0] = -0.577350269189626;
    Location[0][1] = -0.577350269189626;
    Location[0][2] = -0.577350269189626;
    Weight[0] = 1.0;
    Location[1][0] = -0.577350269189626;
    Location[1][1] = -0.577350269189626;
    Location[1][2] = 0.577350269189626;
    Weight[1] = 1.0;
    Location[2][0] = -0.577350269189626;
    Location[2][1] = 0.577350269189626;
    Location[2][2] = -0.577350269189626;
    Weight[2] = 1.0;
    Location[3][0] = -0.577350269189626;
    Location[3][1] = 0.577350269189626;
    Location[3][2] = 0.577350269189626;
    Weight[3] = 1.0;
    Location[4][0] = 0.577350269189626;
    Location[4][1] = -0.577350269189626;
    Location[4][2] = -0.577350269189626;
    Weight[4] = 1.0;
    Location[5][0] = 0.577350269189626;
    Location[5][1] = -0.577350269189626;
    Location[5][2] = 0.577350269189626;
    Weight[5] = 1.0;
    Location[6][0] = 0.577350269189626;
    Location[6][1] = 0.577350269189626;
    Location[6][2] = -0.577350269189626;
    Weight[6] = 1.0;
    Location[7][0] = 0.577350269189626;
    Location[7][1] = 0.577350269189626;
    Location[7][2] = 0.577350269189626;
    Weight[7] = 1.0;
  }

  for (iGauss = 0; iGauss < nGauss; iGauss++) {
    Xi = Location[iGauss][0];
    Eta = Location[iGauss][1];
    Zeta = Location[iGauss][2];

    if (nNodes == 4) Det = ShapeFunc_Tetra(Xi, Eta, Zeta, CoordCorners, DShapeFunction);
    if (nNodes == 5) Det = ShapeFunc_Pyram(Xi, Eta, Zeta, CoordCorners, DShapeFunction);
    if (nNodes == 6) Det = ShapeFunc_Prism(Xi, Eta, Zeta, CoordCorners, DShapeFunction);
    if (nNodes == 8) Det = ShapeFunc_Hexa(Xi, Eta, Zeta, CoordCorners, DShapeFunction);

    /*--- Compute the B Matrix ---*/

    for (iVar = 0; iVar < 6; iVar++)
      for (jVar = 0; jVar < nNodes * nVar; jVar++) B_Matrix[iVar][jVar] = 0.0;

    for (iNode = 0; iNode < nNodes; iNode++) {
      B_Matrix[0][0 + iNode * nVar] = DShapeFunction[iNode][0];
      B_Matrix[1][1 + iNode * nVar] = DShapeFunction[iNode][1];
      B_Matrix[2][2 + iNode * nVar] = DShapeFunction[iNode][2];

      B_Matrix[3][0 + iNode * nVar] = DShapeFunction[iNode][1];
      B_Matrix[3][1 + iNode * nVar] = DShapeFunction[iNode][0];

      B_Matrix[4][1 + iNode * nVar] = DShapeFunction[iNode][2];
      B_Matrix[4][2 + iNode * nVar] = DShapeFunction[iNode][1];

      B_Matrix[5][0 + iNode * nVar] = DShapeFunction[iNode][2];
      B_Matrix[5][2 + iNode * nVar] = DShapeFunction[iNode][0];
    }

    /*--- Impose a type of stiffness for each element ---*/

    switch (config->GetDeform_Stiffness_Type()) {
      case INVERSE_VOLUME:
        E = 1.0 / ElemVolume;
        break;
      case SOLID_WALL_DISTANCE:
        E = 1.0 / ElemDistance;
        break;
      case CONSTANT_STIFFNESS:
        E = 1.0 / EPS;
        break;
    }

    Nu = config->GetDeform_Coeff();
    Mu = E / (2.0 * (1.0 + Nu));
    Lambda = Nu * E / ((1.0 + Nu) * (1.0 - 2.0 * Nu));

    /*--- Compute the D Matrix (for plane strain and 3-D)---*/

    D_Matrix[0][0] = Lambda + 2.0 * Mu;
    D_Matrix[0][1] = Lambda;
    D_Matrix[0][2] = Lambda;
    D_Matrix[1][0] = Lambda;
    D_Matrix[1][1] = Lambda + 2.0 * Mu;
    D_Matrix[1][2] = Lambda;
    D_Matrix[2][0] = Lambda;
    D_Matrix[2][1] = Lambda;
    D_Matrix[2][2] = Lambda + 2.0 * Mu;
    D_Matrix[3][3] = Mu;
    D_Matrix[4][4] = Mu;
    D_Matrix[5][5] = Mu;

    /*--- Compute the BT.D Matrix ---*/

    for (iVar = 0; iVar < nNodes * nVar; iVar++) {
      for (jVar = 0; jVar < 6; jVar++) {
        Aux_Matrix[iVar][jVar] = 0.0;
        for (kVar = 0; kVar < 6; kVar++) Aux_Matrix[iVar][jVar] += B_Matrix[kVar][iVar] * D_Matrix[kVar][jVar];
      }
    }

    /*--- Compute the BT.D.B Matrix (stiffness matrix), and add to the original
     matrix using Gauss integration ---*/

    for (iVar = 0; iVar < nNodes * nVar; iVar++) {
      for (jVar = 0; jVar < nNodes * nVar; jVar++) {
        for (kVar = 0; kVar < 6; kVar++) {
          StiffMatrix_Elem[iVar][jVar] += Weight[iGauss] * Aux_Matrix[iVar][kVar] * B_Matrix[kVar][jVar] * fabs(Det);
        }
      }
    }
  }
}

void CVolumetricMovement::AddFEA_StiffMatrix(CGeometry* geometry, su2double** StiffMatrix_Elem,
                                             unsigned long PointCorners[8], unsigned short nNodes) {
  unsigned short iVar, jVar, iDim, jDim;

  unsigned short nVar = geometry->GetnDim();

  su2double** StiffMatrix_Node;
  StiffMatrix_Node = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) StiffMatrix_Node[iVar] = new su2double[nVar];

  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++) StiffMatrix_Node[iVar][jVar] = 0.0;

  /*--- Transform the stiffness matrix for the hexahedral element into the
   contributions for the individual nodes relative to each other. ---*/

  for (iVar = 0; iVar < nNodes; iVar++) {
    for (jVar = 0; jVar < nNodes; jVar++) {
      for (iDim = 0; iDim < nVar; iDim++) {
        for (jDim = 0; jDim < nVar; jDim++) {
          StiffMatrix_Node[iDim][jDim] = StiffMatrix_Elem[(iVar * nVar) + iDim][(jVar * nVar) + jDim];
        }
      }

      StiffMatrix.AddBlock(PointCorners[iVar], PointCorners[jVar], StiffMatrix_Node);
    }
  }

  /*--- Deallocate memory and exit ---*/

  for (iVar = 0; iVar < nVar; iVar++) delete[] StiffMatrix_Node[iVar];
  delete[] StiffMatrix_Node;
}

void CVolumetricMovement::SetBoundaryDisplacements(CGeometry* geometry, CConfig* config) {
  unsigned short iDim, nDim = geometry->GetnDim(), iMarker, axis = 0;
  unsigned long iPoint, total_index, iVertex;
  su2double *VarCoord, MeanCoord[3] = {0.0, 0.0, 0.0}, VarIncrement = 1.0;

  /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
   deforming meshes (MARKER_MOVING), while SU2_DEF will use it for deforming
   meshes after imposing design variable surface deformations (DV_MARKER). ---*/

  SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();

  /*--- If requested (no by default) impose the surface deflections in
   increments and solve the grid deformation equations iteratively with
   successive small deformations. ---*/

  VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());

  /*--- As initialization, set to zero displacements of all the surfaces except the symmetry
   plane (which is treated specially, see below), internal and the send-receive boundaries ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (((config->GetMarker_All_KindBC(iMarker) != SYMMETRY_PLANE) &&
         (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) &&
         (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY))) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        for (iDim = 0; iDim < nDim; iDim++) {
          total_index = iPoint * nDim + iDim;
          LinSysRes[total_index] = 0.0;
          LinSysSol[total_index] = 0.0;
          StiffMatrix.DeleteValsRowi(total_index);
        }
      }
    }
  }

  /*--- Set the known displacements, note that some points of the moving surfaces
   could be on on the symmetry plane, we should specify DeleteValsRowi again (just in case) ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) ||
        ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DEF)) ||
        ((config->GetDirectDiff() == D_DESIGN) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD) &&
         (config->GetMarker_All_DV(iMarker) == YES)) ||
        ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DOT))) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
        for (iDim = 0; iDim < nDim; iDim++) {
          total_index = iPoint * nDim + iDim;
          LinSysRes[total_index] = SU2_TYPE::GetValue(VarCoord[iDim] * VarIncrement);
          LinSysSol[total_index] = SU2_TYPE::GetValue(VarCoord[iDim] * VarIncrement);
          StiffMatrix.DeleteValsRowi(total_index);
        }
      }
    }
  }

  /*--- Set to zero displacements of the normal component for the symmetry plane condition ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE)) {
      su2double* Coord_0 = nullptr;

      for (iDim = 0; iDim < nDim; iDim++) MeanCoord[iDim] = 0.0;

      /*--- Store the coord of the first point to help identify the axis. ---*/

      iPoint = geometry->vertex[iMarker][0]->GetNode();
      Coord_0 = geometry->nodes->GetCoord(iPoint);

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        VarCoord = geometry->nodes->GetCoord(iPoint);
        for (iDim = 0; iDim < nDim; iDim++)
          MeanCoord[iDim] += (VarCoord[iDim] - Coord_0[iDim]) * (VarCoord[iDim] - Coord_0[iDim]);
      }
      for (iDim = 0; iDim < nDim; iDim++) MeanCoord[iDim] = sqrt(MeanCoord[iDim]);
      if (nDim == 3) {
        if ((MeanCoord[0] <= MeanCoord[1]) && (MeanCoord[0] <= MeanCoord[2])) axis = 0;
        if ((MeanCoord[1] <= MeanCoord[0]) && (MeanCoord[1] <= MeanCoord[2])) axis = 1;
        if ((MeanCoord[2] <= MeanCoord[0]) && (MeanCoord[2] <= MeanCoord[1])) axis = 2;
      } else {
        if ((MeanCoord[0] <= MeanCoord[1])) axis = 0;
        if ((MeanCoord[1] <= MeanCoord[0])) axis = 1;
      }

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        total_index = iPoint * nDim + axis;
        LinSysRes[total_index] = 0.0;
        LinSysSol[total_index] = 0.0;
        StiffMatrix.DeleteValsRowi(total_index);
      }
    }
  }

  /*--- Don't move the nearfield plane ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        for (iDim = 0; iDim < nDim; iDim++) {
          total_index = iPoint * nDim + iDim;
          LinSysRes[total_index] = 0.0;
          LinSysSol[total_index] = 0.0;
          StiffMatrix.DeleteValsRowi(total_index);
        }
      }
    }
  }

  /*--- Move the FSI interfaces ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_ZoneInterface(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
        for (iDim = 0; iDim < nDim; iDim++) {
          total_index = iPoint * nDim + iDim;
          LinSysRes[total_index] = SU2_TYPE::GetValue(VarCoord[iDim] * VarIncrement);
          LinSysSol[total_index] = SU2_TYPE::GetValue(VarCoord[iDim] * VarIncrement);
          StiffMatrix.DeleteValsRowi(total_index);
        }
      }
    }
  }
}

void CVolumetricMovement::SetBoundaryDerivatives(CGeometry* geometry, CConfig* config,
                                                 bool ForwardProjectionDerivative) {
  unsigned short iDim, iMarker;
  unsigned long iPoint, total_index, iVertex;

  su2double* VarCoord;
  SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();
  if ((config->GetDirectDiff() == D_DESIGN) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_DV(iMarker) == YES)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
          for (iDim = 0; iDim < nDim; iDim++) {
            total_index = iPoint * nDim + iDim;
            LinSysRes[total_index] = SU2_TYPE::GetDerivative(VarCoord[iDim]);
            LinSysSol[total_index] = SU2_TYPE::GetDerivative(VarCoord[iDim]);
          }
        }
      }
    }
    if (LinSysRes.norm() == 0.0) cout << "Warning: Derivatives are zero!" << endl;
  } else if ((Kind_SU2 == SU2_COMPONENT::SU2_DOT) && !ForwardProjectionDerivative) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        total_index = iPoint * nDim + iDim;
        LinSysRes[total_index] = SU2_TYPE::GetValue(geometry->GetSensitivity(iPoint, iDim));
        LinSysSol[total_index] = SU2_TYPE::GetValue(geometry->GetSensitivity(iPoint, iDim));
      }
    }
  } else if (config->GetSmoothGradient() && ForwardProjectionDerivative) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_DV(iMarker) == YES)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          for (iDim = 0; iDim < nDim; iDim++) {
            total_index = iPoint * nDim + iDim;
            LinSysRes[total_index] = SU2_TYPE::GetValue(geometry->GetSensitivity(iPoint, iDim));
            LinSysSol[total_index] = SU2_TYPE::GetValue(geometry->GetSensitivity(iPoint, iDim));
          }
        }
      }
    }
    if (LinSysRes.norm() == 0.0) cout << "Warning: Derivatives are zero!" << endl;
  }
}

void CVolumetricMovement::UpdateGridCoord_Derivatives(CGeometry* geometry, CConfig* config,
                                                      bool ForwardProjectionDerivative) {
  unsigned short iDim, iMarker;
  unsigned long iPoint, total_index, iVertex;
  auto* new_coord = new su2double[3];

  SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();

  /*--- Update derivatives of the grid coordinates using the solution of the linear system
     after grid deformation (LinSysSol contains the derivatives of the x, y, z displacements). ---*/
  if ((config->GetDirectDiff() == D_DESIGN) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) {
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      new_coord[0] = 0.0;
      new_coord[1] = 0.0;
      new_coord[2] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        total_index = iPoint * nDim + iDim;
        new_coord[iDim] = geometry->nodes->GetCoord(iPoint, iDim);
        SU2_TYPE::SetDerivative(new_coord[iDim], SU2_TYPE::GetValue(LinSysSol[total_index]));
      }
      geometry->nodes->SetCoord(iPoint, new_coord);
    }
  } else if ((Kind_SU2 == SU2_COMPONENT::SU2_DOT) && !ForwardProjectionDerivative) {
    // need to reset here, since we read out the whole vector, but are only interested in boundary derivatives.
    if (config->GetSmoothGradient()) {
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          total_index = iPoint * nDim + iDim;
          geometry->SetSensitivity(iPoint, iDim, 0.0);
        }
      }
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetSolid_Wall(iMarker) || (config->GetMarker_All_DV(iMarker) == YES)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (geometry->nodes->GetDomain(iPoint)) {
            for (iDim = 0; iDim < nDim; iDim++) {
              total_index = iPoint * nDim + iDim;
              geometry->SetSensitivity(iPoint, iDim, LinSysSol[total_index]);
            }
          }
        }
      }
    }
  } else if (config->GetSmoothGradient() && ForwardProjectionDerivative) {
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        total_index = iPoint * nDim + iDim;
        geometry->SetSensitivity(iPoint, iDim, LinSysSol[total_index]);
      }
    }
  }

  delete[] new_coord;
}

void CVolumetricMovement::SetDomainDisplacements(CGeometry* geometry, CConfig* config) {
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, total_index;

  if (config->GetHold_GridFixed()) {
    auto MinCoordValues = config->GetHold_GridFixed_Coord();
    auto MaxCoordValues = &config->GetHold_GridFixed_Coord()[3];

    /*--- Set to zero displacements of all the points that are not going to be moved
     except the surfaces ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      auto Coord = geometry->nodes->GetCoord(iPoint);
      for (iDim = 0; iDim < nDim; iDim++) {
        if ((Coord[iDim] < MinCoordValues[iDim]) || (Coord[iDim] > MaxCoordValues[iDim])) {
          total_index = iPoint * nDim + iDim;
          LinSysRes[total_index] = 0.0;
          LinSysSol[total_index] = 0.0;
          StiffMatrix.DeleteValsRowi(total_index);
        }
      }
    }
  }

  /*--- Don't move the volume grid outside the limits based
   on the distance to the solid surface ---*/

  if (config->GetDeform_Limit() < 1E6) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      if (geometry->nodes->GetWall_Distance(iPoint) >= config->GetDeform_Limit()) {
        for (iDim = 0; iDim < nDim; iDim++) {
          total_index = iPoint * nDim + iDim;
          LinSysRes[total_index] = 0.0;
          LinSysSol[total_index] = 0.0;
          StiffMatrix.DeleteValsRowi(total_index);
        }
      }
    }
  }
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
