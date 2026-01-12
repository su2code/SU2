/*!
 * \file CLinearElasticity.cpp
 * \brief Subroutines for moving mesh volume elements using the linear elasticity analogy
 * \author F. Palacios, T. Economon, S. Padron
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/grid_movement/CLinearElasticity.hpp"
#include "../../include/adt/CADTPointsOnlyClass.hpp"

CLinearElasticity::CLinearElasticity(CGeometry* geometry, CConfig* config)
    : CVolumetricMovement(geometry), System(LINEAR_SOLVER_MODE::MESH_DEFORM) {
  /*--- Initialize the number of spatial dimensions, length of the state
   vector (same as spatial dimensions for grid deformation), and grid nodes. ---*/

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

CLinearElasticity::~CLinearElasticity(void) = default;

void CLinearElasticity::SetVolume_Deformation(CGeometry* geometry, CConfig* config, bool UpdateGeo, bool Derivative,
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

    nIterMesh = Tot_Iter;

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

void CLinearElasticity::UpdateGridCoord(CGeometry* geometry, CConfig* config) {
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

  geometry->InitiateComms(geometry, config, MPI_QUANTITIES::COORDINATES);
  geometry->CompleteComms(geometry, config, MPI_QUANTITIES::COORDINATES);
}

void CLinearElasticity::UpdateGridCoord_Derivatives(CGeometry* geometry, CConfig* config,
                                                    bool ForwardProjectionDerivative) {
  unsigned short iDim, iMarker;
  unsigned long iPoint, total_index, iVertex;

  SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();

  /*--- Update derivatives of the grid coordinates using the solution of the linear system
     after grid deformation (LinSysSol contains the derivatives of the x, y, z displacements). ---*/
  if ((config->GetDirectDiff() == D_DESIGN) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) {
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      su2double new_coord[3] = {};
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
}

void CLinearElasticity::ComputeSolid_Wall_Distance(CGeometry* geometry, CConfig* config, su2double& MinDistance,
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

su2double CLinearElasticity::SetFEAMethodContributions_Elem(CGeometry* geometry, CConfig* config) {
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

void CLinearElasticity::SetFEA_StiffMatrix2D(CGeometry* geometry, CConfig* config, su2double** StiffMatrix_Elem,
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

  for (iVar = 0; iVar < static_cast<unsigned short>(nNodes * nVar); iVar++) {
    for (jVar = 0; jVar < static_cast<unsigned short>(nNodes * nVar); jVar++) {
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
      for (jVar = 0; jVar < static_cast<unsigned short>(nNodes * nVar); jVar++) B_Matrix[iVar][jVar] = 0.0;

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

    for (iVar = 0; iVar < static_cast<unsigned short>(nNodes * nVar); iVar++) {
      for (jVar = 0; jVar < 3; jVar++) {
        Aux_Matrix[iVar][jVar] = 0.0;
        for (kVar = 0; kVar < 3; kVar++) Aux_Matrix[iVar][jVar] += B_Matrix[kVar][iVar] * D_Matrix[kVar][jVar];
      }
    }

    /*--- Compute the BT.D.B Matrix (stiffness matrix), and add to the original
     matrix using Gauss integration ---*/

    for (iVar = 0; iVar < static_cast<unsigned short>(nNodes * nVar); iVar++) {
      for (jVar = 0; jVar < static_cast<unsigned short>(nNodes * nVar); jVar++) {
        for (kVar = 0; kVar < 3; kVar++) {
          StiffMatrix_Elem[iVar][jVar] += Weight[iGauss] * Aux_Matrix[iVar][kVar] * B_Matrix[kVar][jVar] * fabs(Det);
        }
      }
    }
  }
}

void CLinearElasticity::SetFEA_StiffMatrix3D(CGeometry* geometry, CConfig* config, su2double** StiffMatrix_Elem,
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

  for (iVar = 0; iVar < static_cast<unsigned short>(nNodes * nVar); iVar++) {
    for (jVar = 0; jVar < static_cast<unsigned short>(nNodes * nVar); jVar++) {
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
      for (jVar = 0; jVar < static_cast<unsigned short>(nNodes * nVar); jVar++) B_Matrix[iVar][jVar] = 0.0;

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

    for (iVar = 0; iVar < static_cast<unsigned short>(nNodes * nVar); iVar++) {
      for (jVar = 0; jVar < 6; jVar++) {
        Aux_Matrix[iVar][jVar] = 0.0;
        for (kVar = 0; kVar < 6; kVar++) Aux_Matrix[iVar][jVar] += B_Matrix[kVar][iVar] * D_Matrix[kVar][jVar];
      }
    }

    /*--- Compute the BT.D.B Matrix (stiffness matrix), and add to the original
     matrix using Gauss integration ---*/

    for (iVar = 0; iVar < static_cast<unsigned short>(nNodes * nVar); iVar++) {
      for (jVar = 0; jVar < static_cast<unsigned short>(nNodes * nVar); jVar++) {
        for (kVar = 0; kVar < 6; kVar++) {
          StiffMatrix_Elem[iVar][jVar] += Weight[iGauss] * Aux_Matrix[iVar][kVar] * B_Matrix[kVar][jVar] * fabs(Det);
        }
      }
    }
  }
}

void CLinearElasticity::AddFEA_StiffMatrix(CGeometry* geometry, su2double** StiffMatrix_Elem,
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

su2double CLinearElasticity::ShapeFunc_Triangle(su2double Xi, su2double Eta, su2double CoordCorners[8][3],
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

su2double CLinearElasticity::ShapeFunc_Quadrilateral(su2double Xi, su2double Eta, su2double CoordCorners[8][3],
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

su2double CLinearElasticity::ShapeFunc_Tetra(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3],
                                             su2double DShapeFunction[8][4]) {
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

su2double CLinearElasticity::ShapeFunc_Pyram(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3],
                                             su2double DShapeFunction[8][4]) {
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

su2double CLinearElasticity::ShapeFunc_Prism(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3],
                                             su2double DShapeFunction[8][4]) {
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

su2double CLinearElasticity::ShapeFunc_Hexa(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3],
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

void CLinearElasticity::SetDomainDisplacements(CGeometry* geometry, CConfig* config) {
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

void CLinearElasticity::SetBoundaryDisplacements(CGeometry* geometry, CConfig* config) {
  unsigned short iDim, nDim = geometry->GetnDim(), iMarker, axis = 0;
  unsigned long iPoint, total_index, iVertex;
  su2double *VarCoord, MeanCoord[3] = {0.0, 0.0, 0.0}, VarIncrement = 1.0;

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
    if (IsDeformationMarker(config, iMarker)) {
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
    if ((config->GetMarker_All_ZoneInterface(iMarker) == YES) && (config->GetKind_SU2() == SU2_COMPONENT::SU2_CFD)) {
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

void CLinearElasticity::SetBoundaryDerivatives(CGeometry* geometry, CConfig* config, bool ForwardProjectionDerivative) {
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
