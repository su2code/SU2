/*!
 * \file grid_movement_structure.cpp
 * \brief Subroutines for doing the grid movement using different strategies
 * \author F. Palacios, T. Economon, S. Padron
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../include/grid_movement_structure.hpp"
#include "../include/adt_structure.hpp"
#include <list>

#include "../include/linear_algebra/CMatrixVectorProduct.hpp"
#include "../include/linear_algebra/CPreconditioner.hpp"

using namespace std;

CGridMovement::CGridMovement(void) { }

CGridMovement::~CGridMovement(void) { }

CVolumetricMovement::CVolumetricMovement(void) : CGridMovement() {



}

CVolumetricMovement::CVolumetricMovement(CGeometry *geometry, CConfig *config) : CGridMovement() {
  
  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();
  
    /*--- Initialize the number of spatial dimensions, length of the state
     vector (same as spatial dimensions for grid deformation), and grid nodes. ---*/

    nDim   = geometry->GetnDim();
    nVar   = geometry->GetnDim();
    nPoint = geometry->GetnPoint();
    nPointDomain = geometry->GetnPointDomain();

    nIterMesh = 0;

    /*--- Initialize matrix, solution, and r.h.s. structures for the linear solver. ---*/
    if (config->GetVolumetric_Movement()){
      LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
      LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
      StiffMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
    }
}

CVolumetricMovement::~CVolumetricMovement(void) { }

void CVolumetricMovement::UpdateGridCoord(CGeometry *geometry, CConfig *config) {
  
  unsigned short iDim;
  unsigned long iPoint, total_index;
  su2double new_coord;
  
  /*--- Update the grid coordinates using the solution of the linear system
   after grid deformation (LinSysSol contains the x, y, z displacements). ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iDim = 0; iDim < nDim; iDim++) {
      total_index = iPoint*nDim + iDim;
      new_coord = geometry->node[iPoint]->GetCoord(iDim)+LinSysSol[total_index];
      if (fabs(new_coord) < EPS*EPS) new_coord = 0.0;
      geometry->node[iPoint]->SetCoord(iDim, new_coord);
    }

  /*--- LinSysSol contains the non-transformed displacements in the periodic halo cells.
   * Hence we still need a communication of the transformed coordinates, otherwise periodicity
   * is not maintained. ---*/

  geometry->InitiateComms(geometry, config, COORDINATES);
  geometry->CompleteComms(geometry, config, COORDINATES);
  
}

void CVolumetricMovement::UpdateDualGrid(CGeometry *geometry, CConfig *config) {
  
  /*--- After moving all nodes, update the dual mesh. Recompute the edges and
   dual mesh control volumes in the domain and on the boundaries. ---*/
  
  geometry->SetCoord_CG();
  geometry->SetControlVolume(config, UPDATE);
  geometry->SetBoundControlVolume(config, UPDATE);
  geometry->SetMaxLength(config);
  
}

void CVolumetricMovement::UpdateMultiGrid(CGeometry **geometry, CConfig *config) {
  
  unsigned short iMGfine, iMGlevel, nMGlevel = config->GetnMGLevels();
  
  /*--- Update the multigrid structure after moving the finest grid,
   including computing the grid velocities on the coarser levels. ---*/
  
  for (iMGlevel = 1; iMGlevel <= nMGlevel; iMGlevel++) {
    iMGfine = iMGlevel-1;
    geometry[iMGlevel]->SetControlVolume(config, geometry[iMGfine], UPDATE);
    geometry[iMGlevel]->SetBoundControlVolume(config, geometry[iMGfine],UPDATE);
    geometry[iMGlevel]->SetCoord(geometry[iMGfine]);
    if (config->GetGrid_Movement())
      geometry[iMGlevel]->SetRestricted_GridVelocity(geometry[iMGfine], config);
  }
 
}

void CVolumetricMovement::SetVolume_Deformation(CGeometry *geometry, CConfig *config, bool UpdateGeo, bool Derivative) {
  
  unsigned long IterLinSol = 0, Smoothing_Iter, iNonlinear_Iter, MaxIter = 0, RestartIter = 50, Tot_Iter = 0, Nonlinear_Iter = 0;
  su2double MinVolume, MaxVolume, NumError, Residual = 0.0, Residual_Init = 0.0;
  bool Screen_Output;

  /*--- Retrieve number or iterations, tol, output, etc. from config ---*/
  
  Smoothing_Iter = config->GetDeform_Linear_Solver_Iter();
  Screen_Output  = config->GetDeform_Output();
  NumError       = config->GetDeform_Linear_Solver_Error();
  Nonlinear_Iter = config->GetGridDef_Nonlinear_Iter();
  
  /*--- Disable the screen output if we're running SU2_CFD ---*/
  
  if (config->GetKind_SU2() == SU2_CFD && !Derivative) Screen_Output = false;

  /*--- Set the number of nonlinear iterations to 1 if Derivative computation is enabled ---*/

  if (Derivative) Nonlinear_Iter = 1;
  
  /*--- Loop over the total number of grid deformation iterations. The surface
   deformation can be divided into increments to help with stability. In
   particular, the linear elasticity equations hold only for small deformations. ---*/
  
  for (iNonlinear_Iter = 0; iNonlinear_Iter < Nonlinear_Iter; iNonlinear_Iter++) {
    
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

    if (Derivative) { SetBoundaryDerivatives(geometry, config); }
    
    CMatrixVectorProduct<su2double>* mat_vec = NULL;
    CPreconditioner<su2double>* precond = NULL;

    /*--- Communicate any prescribed boundary displacements via MPI,
     so that all nodes have the same solution and r.h.s. entries
     across all partitions. ---*/
    
    StiffMatrix.InitiateComms(LinSysSol, geometry, config, SOLUTION_MATRIX);
    StiffMatrix.CompleteComms(LinSysSol, geometry, config, SOLUTION_MATRIX);
    
    StiffMatrix.InitiateComms(LinSysRes, geometry, config, SOLUTION_MATRIX);
    StiffMatrix.CompleteComms(LinSysRes, geometry, config, SOLUTION_MATRIX);

    /*--- Definition of the preconditioner matrix vector multiplication, and linear solver ---*/

    /*--- If we want no derivatives or the direct derivatives,
     * we solve the system using the normal matrix vector product and preconditioner.
     * For the mesh sensitivities using the discrete adjoint method we solve the system using the transposed matrix,
     * hence we need the corresponding matrix vector product and the preconditioner.  ---*/
    if (!Derivative || ((config->GetKind_SU2() == SU2_CFD) && Derivative)) {

      if (config->GetKind_Deform_Linear_Solver_Prec() == LU_SGS) {
        if ((rank == MASTER_NODE) && Screen_Output) cout << "\n# LU_SGS preconditioner." << endl;
    		mat_vec = new CSysMatrixVectorProduct<su2double>(StiffMatrix, geometry, config);
    		precond = new CLU_SGSPreconditioner<su2double>(StiffMatrix, geometry, config);
    	}
    	if (config->GetKind_Deform_Linear_Solver_Prec() == ILU) {
        if ((rank == MASTER_NODE) && Screen_Output) cout << "\n# ILU preconditioner." << endl;
    		StiffMatrix.BuildILUPreconditioner();
    		mat_vec = new CSysMatrixVectorProduct<su2double>(StiffMatrix, geometry, config);
    		precond = new CILUPreconditioner<su2double>(StiffMatrix, geometry, config, false);
    	}
    	if (config->GetKind_Deform_Linear_Solver_Prec() == JACOBI) {
        if ((rank == MASTER_NODE) && Screen_Output) cout << "\n# Jacobi preconditioner." << endl;
    		StiffMatrix.BuildJacobiPreconditioner();
    		mat_vec = new CSysMatrixVectorProduct<su2double>(StiffMatrix, geometry, config);
    		precond = new CJacobiPreconditioner<su2double>(StiffMatrix, geometry, config, false);
    	}

    } else if (Derivative && (config->GetKind_SU2() == SU2_DOT)) {

      /*--- Build the ILU or Jacobi preconditioner for the transposed system ---*/

      if ((config->GetKind_Deform_Linear_Solver_Prec() == ILU) ||
          (config->GetKind_Deform_Linear_Solver_Prec() == LU_SGS)) {
        if ((rank == MASTER_NODE) && Screen_Output) cout << "\n# ILU preconditioner." << endl;
    		StiffMatrix.BuildILUPreconditioner(true);
    		mat_vec = new CSysMatrixVectorProductTransposed<su2double>(StiffMatrix, geometry, config);
    		precond = new CILUPreconditioner<su2double>(StiffMatrix, geometry, config, true);
    	}
    	if (config->GetKind_Deform_Linear_Solver_Prec() == JACOBI) {
        if ((rank == MASTER_NODE) && Screen_Output) cout << "\n# Jacobi preconditioner." << endl;
    		StiffMatrix.BuildJacobiPreconditioner(true);
    		mat_vec = new CSysMatrixVectorProductTransposed<su2double>(StiffMatrix, geometry, config);
    		precond = new CJacobiPreconditioner<su2double>(StiffMatrix, geometry, config, true);
    	}

    }
    
    if (LinSysRes.norm() != 0.0){
      switch (config->GetKind_Deform_Linear_Solver()) {
        
        /*--- Solve the linear system (GMRES with restart) ---*/
        
        case RESTARTED_FGMRES:

          Tot_Iter = 0; MaxIter = RestartIter;

          System.FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, NumError, 1, Residual_Init, false, config);

          if ((rank == MASTER_NODE) && Screen_Output) {
            cout << "\n# FGMRES (with restart) residual history" << endl;
            cout << "# Residual tolerance target = " << NumError << endl;
            cout << "# Initial residual norm     = " << Residual_Init << endl;
          }

          if (rank == MASTER_NODE) { cout << "     " << Tot_Iter << "     " << Residual_Init/Residual_Init << endl; }

          while (Tot_Iter < Smoothing_Iter) {

            if (IterLinSol + RestartIter > Smoothing_Iter)
              MaxIter = Smoothing_Iter - IterLinSol;

            IterLinSol = System.FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, NumError, MaxIter, Residual, false, config);
            Tot_Iter += IterLinSol;

            if ((rank == MASTER_NODE) && Screen_Output) { cout << "     " << Tot_Iter << "     " << Residual/Residual_Init << endl; }

            if (Residual < Residual_Init*NumError) { break; }

          }

          if ((rank == MASTER_NODE) && Screen_Output) {
            cout << "# FGMRES (with restart) final (true) residual:" << endl;
            cout << "# Iteration = " << Tot_Iter << ": |res|/|res0| = " << Residual/Residual_Init << ".\n" << endl;
          }

          break;

          /*--- Solve the linear system (GMRES) ---*/

        case FGMRES:

          Tot_Iter = System.FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, NumError, Smoothing_Iter, Residual, Screen_Output, config);

          break;

          /*--- Solve the linear system (BCGSTAB) ---*/

        case BCGSTAB:

          Tot_Iter = System.BCGSTAB_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, NumError, Smoothing_Iter, Residual, Screen_Output, config);

          break;


        case CONJUGATE_GRADIENT:

          Tot_Iter = System.CG_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, NumError, Smoothing_Iter, Residual, Screen_Output, config);

          break;

      }
    }
    
    /*--- Deallocate memory needed by the Krylov linear solver ---*/
    
    delete mat_vec;
    delete precond;
    
    /*--- Update the grid coordinates and cell volumes using the solution
     of the linear system (usol contains the x, y, z displacements). ---*/

    if (!Derivative) { UpdateGridCoord(geometry, config); }
    else { UpdateGridCoord_Derivatives(geometry, config); }
    if (UpdateGeo) { UpdateDualGrid(geometry, config); }
    
    /*--- Check for failed deformation (negative volumes). ---*/
    
    ComputeDeforming_Element_Volume(geometry, MinVolume, MaxVolume, Screen_Output);
    
    /*--- Set number of iterations in the mesh update. ---*/

    Set_nIterMesh(Tot_Iter);

    if (rank == MASTER_NODE && Screen_Output) {
      cout << "Non-linear iter.: " << iNonlinear_Iter+1 << "/" << Nonlinear_Iter  << ". Linear iter.: " << Tot_Iter << ". ";
      if (nDim == 2) cout << "Min. area: " << MinVolume << ". Error: " << Residual << "." << endl;
      else cout << "Min. volume: " << MinVolume << ". Error: " << Residual << "." << endl;
    }
    
  }
  

}

void CVolumetricMovement::ComputeDeforming_Element_Volume(CGeometry *geometry, su2double &MinVolume, su2double &MaxVolume, bool Screen_Output) {
  
  unsigned long iElem, ElemCounter = 0, PointCorners[8];
  su2double Volume = 0.0, CoordCorners[8][3];
  unsigned short nNodes = 0, iNodes, iDim;
  bool RightVol = true;
  
  if (rank == MASTER_NODE && Screen_Output)
    cout << "Computing volumes of the grid elements." << endl;
  
  MaxVolume = -1E22; MinVolume = 1E22;
  
  /*--- Load up each triangle and tetrahedron to check for negative volumes. ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     nNodes = 3;
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL)    nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      nNodes = 5;
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        nNodes = 6;
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   nNodes = 8;
    
    for (iNodes = 0; iNodes < nNodes; iNodes++) {
      PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
      for (iDim = 0; iDim < nDim; iDim++) {
        CoordCorners[iNodes][iDim] = geometry->node[PointCorners[iNodes]]->GetCoord(iDim);
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
  unsigned long ElemCounter_Local = ElemCounter; ElemCounter = 0;
  su2double MaxVolume_Local = MaxVolume; MaxVolume = 0.0;
  su2double MinVolume_Local = MinVolume; MinVolume = 0.0;
  SU2_MPI::Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MaxVolume_Local, &MaxVolume, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MinVolume_Local, &MinVolume, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
  
  /*--- Volume from  0 to 1 ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    Volume = geometry->elem[iElem]->GetVolume()/MaxVolume;
    geometry->elem[iElem]->SetVolume(Volume);
  }
  
  if ((ElemCounter != 0) && (rank == MASTER_NODE) && (Screen_Output))
    cout <<"There are " << ElemCounter << " elements with negative volume.\n" << endl;
  
}

  
  
void CVolumetricMovement::ComputeSolid_Wall_Distance(CGeometry *geometry, CConfig *config, su2double &MinDistance, su2double &MaxDistance) {
  
  unsigned long nVertex_SolidWall, ii, jj, iVertex, iPoint, pointID;
  unsigned short iMarker, iDim;
  su2double dist, MaxDistance_Local, MinDistance_Local;
  int rankID;

  /*--- Initialize min and max distance ---*/

  MaxDistance = -1E22; MinDistance = 1E22;
  
  /*--- Compute the total number of nodes on no-slip boundaries ---*/
  
  nVertex_SolidWall = 0;
  for(iMarker=0; iMarker<config->GetnMarker_All(); ++iMarker) {
    if( (config->GetMarker_All_KindBC(iMarker) == EULER_WALL ||
         config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)  ||
       (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ||
        (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE)) {
      nVertex_SolidWall += geometry->GetnVertex(iMarker);
    }
  }
  
  /*--- Allocate the vectors to hold boundary node coordinates
   and its local ID. ---*/
  
  vector<su2double>     Coord_bound(nDim*nVertex_SolidWall);
  vector<unsigned long> PointIDs(nVertex_SolidWall);
  
  /*--- Retrieve and store the coordinates of the no-slip boundary nodes
   and their local point IDs. ---*/
  
  ii = 0; jj = 0;
  for (iMarker=0; iMarker<config->GetnMarker_All(); ++iMarker) {
    if ( (config->GetMarker_All_KindBC(iMarker) == EULER_WALL ||
         config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)  ||
       (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)  ||
         (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE)) {
      for (iVertex=0; iVertex<geometry->GetnVertex(iMarker); ++iVertex) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        PointIDs[jj++] = iPoint;
        for (iDim=0; iDim<nDim; ++iDim)
          Coord_bound[ii++] = geometry->node[iPoint]->GetCoord(iDim);
      }
    }
  }
  
  /*--- Build the ADT of the boundary nodes. ---*/
  
  CADTPointsOnlyClass WallADT(nDim, nVertex_SolidWall, Coord_bound.data(),
                              PointIDs.data(), true);
  
  /*--- Loop over all interior mesh nodes and compute the distances to each
   of the no-slip boundary nodes. Store the minimum distance to the wall
   for each interior mesh node. ---*/
  
  if( WallADT.IsEmpty() ) {
    
    /*--- No solid wall boundary nodes in the entire mesh.
     Set the wall distance to zero for all nodes. ---*/
    
    for (iPoint=0; iPoint<geometry->GetnPoint(); ++iPoint)
      geometry->node[iPoint]->SetWall_Distance(0.0);
  }
  else {
    
    /*--- Solid wall boundary nodes are present. Compute the wall
     distance for all nodes. ---*/
    
    for(iPoint=0; iPoint<geometry->GetnPoint(); ++iPoint) {
      
      WallADT.DetermineNearestNode(geometry->node[iPoint]->GetCoord(), dist,
                                   pointID, rankID);
      geometry->node[iPoint]->SetWall_Distance(dist);
      
      MaxDistance = max(MaxDistance, dist);
      
      /*--- To discard points on the surface we use > EPS ---*/
       
      if (sqrt(dist) > EPS)  MinDistance = min(MinDistance, dist);
      
    }
    
    MaxDistance_Local = MaxDistance; MaxDistance = 0.0;
    MinDistance_Local = MinDistance; MinDistance = 0.0;
    
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&MaxDistance_Local, &MaxDistance, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MinDistance_Local, &MinDistance, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else
    MaxDistance = MaxDistance_Local;
    MinDistance = MinDistance_Local;
#endif
    
  }
  
}

su2double CVolumetricMovement::SetFEAMethodContributions_Elem(CGeometry *geometry, CConfig *config) {
  
  unsigned short iVar, iDim, nNodes = 0, iNodes, StiffMatrix_nElem = 0;
  unsigned long iElem, PointCorners[8];
  su2double **StiffMatrix_Elem = NULL, CoordCorners[8][3];
  su2double MinVolume = 0.0, MaxVolume = 0.0, MinDistance = 0.0, MaxDistance = 0.0, ElemVolume = 0.0, ElemDistance = 0.0;
  
  bool Screen_Output  = config->GetDeform_Output();
  
  /*--- Allocate maximum size (quadrilateral and hexahedron) ---*/
  
  if (nDim == 2) StiffMatrix_nElem = 8;
  else StiffMatrix_nElem = 24;
    
  StiffMatrix_Elem = new su2double* [StiffMatrix_nElem];
  for (iVar = 0; iVar < StiffMatrix_nElem; iVar++)
    StiffMatrix_Elem[iVar] = new su2double [StiffMatrix_nElem];
  
  /*--- Compute min volume in the entire mesh. ---*/
  
  ComputeDeforming_Element_Volume(geometry, MinVolume, MaxVolume, Screen_Output);
  if (rank == MASTER_NODE && Screen_Output) cout <<"Min. volume: "<< MinVolume <<", max. volume: "<< MaxVolume <<"." << endl;
  
  /*--- Compute the distance to the nearest surface if needed
   as part of the stiffness calculation.. ---*/

  if ((config->GetDeform_Stiffness_Type() == SOLID_WALL_DISTANCE) ||
      (config->GetDeform_Limit() < 1E6)) {
    ComputeSolid_Wall_Distance(geometry, config, MinDistance, MaxDistance);
    if (rank == MASTER_NODE && Screen_Output) cout <<"Min. distance: "<< MinDistance <<", max. distance: "<< MaxDistance <<"." << endl;
  }
  
  /*--- Compute contributions from each element by forming the stiffness matrix (FEA) ---*/
  
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      nNodes = 3;
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       nNodes = 5;
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         nNodes = 6;
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    nNodes = 8;
    
    for (iNodes = 0; iNodes < nNodes; iNodes++) {
      PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
      for (iDim = 0; iDim < nDim; iDim++) {
        CoordCorners[iNodes][iDim] = geometry->node[PointCorners[iNodes]]->GetCoord(iDim);
      }
    }
    
    /*--- Extract Element volume and distance to compute the stiffness ---*/
    
    ElemVolume = geometry->elem[iElem]->GetVolume();
    
    if ((config->GetDeform_Stiffness_Type() == SOLID_WALL_DISTANCE)) {
      ElemDistance = 0.0;
      for (iNodes = 0; iNodes < nNodes; iNodes++)
        ElemDistance += geometry->node[PointCorners[iNodes]]->GetWall_Distance();
      ElemDistance = ElemDistance/(su2double)nNodes;
    }
    
    if (nDim == 2) SetFEA_StiffMatrix2D(geometry, config, StiffMatrix_Elem, PointCorners, CoordCorners, nNodes, ElemVolume, ElemDistance);
    if (nDim == 3) SetFEA_StiffMatrix3D(geometry, config, StiffMatrix_Elem, PointCorners, CoordCorners, nNodes, ElemVolume, ElemDistance);
    
    AddFEA_StiffMatrix(geometry, StiffMatrix_Elem, PointCorners, nNodes);
    
  }
  
  /*--- Deallocate memory and exit ---*/
  
  for (iVar = 0; iVar < StiffMatrix_nElem; iVar++)
    delete [] StiffMatrix_Elem[iVar];
  delete [] StiffMatrix_Elem;
  
  return MinVolume;

}

su2double CVolumetricMovement::ShapeFunc_Triangle(su2double Xi, su2double Eta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]) {
  
  int i, j, k;
  su2double c0, c1, xsj;
  su2double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = Xi;
  DShapeFunction[1][3] = Eta;
  DShapeFunction[2][3] = 1-Xi-Eta;
  
  /*--- dN/d xi, dN/d eta ---*/
  
  DShapeFunction[0][0] = 1.0;  DShapeFunction[0][1] = 0.0;
  DShapeFunction[1][0] = 0.0;  DShapeFunction[1][1] = 1.0;
  DShapeFunction[2][0] = -1.0; DShapeFunction[2][1] = -1.0;
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 3; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to Jacobian ---*/
  
  ad[0][0] = xs[1][1];
  ad[0][1] = -xs[0][1];
  ad[1][0] = -xs[1][0];
  ad[1][1] = xs[0][0];
  
  /*--- Determinant of Jacobian ---*/
  
  xsj = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 3; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]; // dN/dy
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
  }
  
  return xsj;
  
}

su2double CVolumetricMovement::ShapeFunc_Quadrilateral(su2double Xi, su2double Eta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]) {
  
  int i, j, k;
  su2double c0, c1, xsj;
  su2double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = 0.25*(1.0-Xi)*(1.0-Eta);
  DShapeFunction[1][3] = 0.25*(1.0+Xi)*(1.0-Eta);
  DShapeFunction[2][3] = 0.25*(1.0+Xi)*(1.0+Eta);
  DShapeFunction[3][3] = 0.25*(1.0-Xi)*(1.0+Eta);
  
  /*--- dN/d xi, dN/d eta ---*/
  
  DShapeFunction[0][0] = -0.25*(1.0-Eta); DShapeFunction[0][1] = -0.25*(1.0-Xi);
  DShapeFunction[1][0] =  0.25*(1.0-Eta); DShapeFunction[1][1] = -0.25*(1.0+Xi);
  DShapeFunction[2][0] =  0.25*(1.0+Eta); DShapeFunction[2][1] =  0.25*(1.0+Xi);
  DShapeFunction[3][0] = -0.25*(1.0+Eta); DShapeFunction[3][1] =  0.25*(1.0-Xi);
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 4; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to Jacobian ---*/
  
  ad[0][0] = xs[1][1];
  ad[0][1] = -xs[0][1];
  ad[1][0] = -xs[1][0];
  ad[1][1] = xs[0][0];
  
  /*--- Determinant of Jacobian ---*/
  
  xsj = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 4; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]; // dN/dy
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
  }
  
  return xsj;
  
}

su2double CVolumetricMovement::ShapeFunc_Tetra(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]) {
  
  int i, j, k;
  su2double c0, c1, c2, xsj;
  su2double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = Xi;
  DShapeFunction[1][3] = Zeta;
  DShapeFunction[2][3] = 1.0 - Xi - Eta - Zeta;
  DShapeFunction[3][3] = Eta;
  
  /*--- dN/d xi, dN/d eta, dN/d zeta ---*/
  
  DShapeFunction[0][0] = 1.0;  DShapeFunction[0][1] = 0.0;  DShapeFunction[0][2] = 0.0;
  DShapeFunction[1][0] = 0.0;  DShapeFunction[1][1] = 0.0;  DShapeFunction[1][2] = 1.0;
  DShapeFunction[2][0] = -1.0; DShapeFunction[2][1] = -1.0; DShapeFunction[2][2] = -1.0;
  DShapeFunction[3][0] = 0.0;  DShapeFunction[3][1] = 1.0;  DShapeFunction[3][2] = 0.0;
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 4; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to Jacobian ---*/
  
  ad[0][0] = xs[1][1]*xs[2][2]-xs[1][2]*xs[2][1];
  ad[0][1] = xs[0][2]*xs[2][1]-xs[0][1]*xs[2][2];
  ad[0][2] = xs[0][1]*xs[1][2]-xs[0][2]*xs[1][1];
  ad[1][0] = xs[1][2]*xs[2][0]-xs[1][0]*xs[2][2];
  ad[1][1] = xs[0][0]*xs[2][2]-xs[0][2]*xs[2][0];
  ad[1][2] = xs[0][2]*xs[1][0]-xs[0][0]*xs[1][2];
  ad[2][0] = xs[1][0]*xs[2][1]-xs[1][1]*xs[2][0];
  ad[2][1] = xs[0][1]*xs[2][0]-xs[0][0]*xs[2][1];
  ad[2][2] = xs[0][0]*xs[1][1]-xs[0][1]*xs[1][0];
  
  /*--- Determinant of Jacobian ---*/
  
  xsj = xs[0][0]*ad[0][0]+xs[0][1]*ad[1][0]+xs[0][2]*ad[2][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 4; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]+xs[0][2]*DShapeFunction[k][2]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]+xs[1][2]*DShapeFunction[k][2]; // dN/dy
    c2 = xs[2][0]*DShapeFunction[k][0]+xs[2][1]*DShapeFunction[k][1]+xs[2][2]*DShapeFunction[k][2]; // dN/dz
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d zeta
  }
  
  return xsj;
  
}

su2double CVolumetricMovement::ShapeFunc_Pyram(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]) {
  
  int i, j, k;
  su2double c0, c1, c2, xsj;
  su2double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = 0.25*(-Xi+Eta+Zeta-1.0)*(-Xi-Eta+Zeta-1.0)/(1.0-Zeta);
  DShapeFunction[1][3] = 0.25*(-Xi-Eta+Zeta-1.0)*( Xi-Eta+Zeta-1.0)/(1.0-Zeta);
  DShapeFunction[2][3] = 0.25*( Xi+Eta+Zeta-1.0)*( Xi-Eta+Zeta-1.0)/(1.0-Zeta);
  DShapeFunction[3][3] = 0.25*( Xi+Eta+Zeta-1.0)*(-Xi+Eta+Zeta-1.0)/(1.0-Zeta);
  DShapeFunction[4][3] = Zeta;
  
  /*--- dN/d xi ---*/
  
  DShapeFunction[0][0] = 0.5*(Zeta-Xi-1.0)/(Zeta-1.0);
  DShapeFunction[1][0] = 0.5*Xi/(Zeta-1.0);
  DShapeFunction[2][0] = 0.5*(1.0-Zeta-Xi)/(Zeta-1.0);
  DShapeFunction[3][0] = DShapeFunction[1][0];
  DShapeFunction[4][0] = 0.0;
  
  /*--- dN/d eta ---*/
  
  DShapeFunction[0][1] = 0.5*Eta/(Zeta-1.0);
  DShapeFunction[1][1] = 0.5*(Zeta-Eta-1.0)/(Zeta-1.0);
  DShapeFunction[2][1] = DShapeFunction[0][1];
  DShapeFunction[3][1] = 0.5*(1.0-Zeta-Eta)/(Zeta-1.0);
  DShapeFunction[4][1] = 0.0;
  
  /*--- dN/d zeta ---*/
  
  DShapeFunction[0][2] = 0.25*(-1.0 + 2.0*Zeta - Zeta*Zeta - Eta*Eta + Xi*Xi)/((1.0-Zeta)*(1.0-Zeta));
  DShapeFunction[1][2] = 0.25*(-1.0 + 2.0*Zeta - Zeta*Zeta + Eta*Eta - Xi*Xi)/((1.0-Zeta)*(1.0-Zeta));
  DShapeFunction[2][2] = DShapeFunction[0][2];
  DShapeFunction[3][2] = DShapeFunction[1][2];
  DShapeFunction[4][2] = 1.0;
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 5; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to Jacobian ---*/
  
  ad[0][0] = xs[1][1]*xs[2][2]-xs[1][2]*xs[2][1];
  ad[0][1] = xs[0][2]*xs[2][1]-xs[0][1]*xs[2][2];
  ad[0][2] = xs[0][1]*xs[1][2]-xs[0][2]*xs[1][1];
  ad[1][0] = xs[1][2]*xs[2][0]-xs[1][0]*xs[2][2];
  ad[1][1] = xs[0][0]*xs[2][2]-xs[0][2]*xs[2][0];
  ad[1][2] = xs[0][2]*xs[1][0]-xs[0][0]*xs[1][2];
  ad[2][0] = xs[1][0]*xs[2][1]-xs[1][1]*xs[2][0];
  ad[2][1] = xs[0][1]*xs[2][0]-xs[0][0]*xs[2][1];
  ad[2][2] = xs[0][0]*xs[1][1]-xs[0][1]*xs[1][0];
  
  /*--- Determinant of Jacobian ---*/
  
  xsj = xs[0][0]*ad[0][0]+xs[0][1]*ad[1][0]+xs[0][2]*ad[2][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 5; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]+xs[0][2]*DShapeFunction[k][2]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]+xs[1][2]*DShapeFunction[k][2]; // dN/dy
    c2 = xs[2][0]*DShapeFunction[k][0]+xs[2][1]*DShapeFunction[k][1]+xs[2][2]*DShapeFunction[k][2]; // dN/dz
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d zeta
  }
  
  return xsj;
  
}

su2double CVolumetricMovement::ShapeFunc_Prism(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]) {
  
  int i, j, k;
  su2double c0, c1, c2, xsj;
  su2double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = 0.5*Eta*(1.0-Xi);
  DShapeFunction[1][3] = 0.5*Zeta*(1.0-Xi);
  DShapeFunction[2][3] = 0.5*(1.0-Eta-Zeta)*(1.0-Xi);
  DShapeFunction[3][3] = 0.5*Eta*(Xi+1.0);
  DShapeFunction[4][3] = 0.5*Zeta*(Xi+1.0);
  DShapeFunction[5][3] = 0.5*(1.0-Eta-Zeta)*(Xi+1.0);
  
  /*--- dN/d Xi, dN/d Eta, dN/d Zeta ---*/
  
  DShapeFunction[0][0] = -0.5*Eta;            DShapeFunction[0][1] = 0.5*(1.0-Xi);      DShapeFunction[0][2] = 0.0;
  DShapeFunction[1][0] = -0.5*Zeta;           DShapeFunction[1][1] = 0.0;               DShapeFunction[1][2] = 0.5*(1.0-Xi);
  DShapeFunction[2][0] = -0.5*(1.0-Eta-Zeta); DShapeFunction[2][1] = -0.5*(1.0-Xi);     DShapeFunction[2][2] = -0.5*(1.0-Xi);
  DShapeFunction[3][0] = 0.5*Eta;             DShapeFunction[3][1] = 0.5*(Xi+1.0);      DShapeFunction[3][2] = 0.0;
  DShapeFunction[4][0] = 0.5*Zeta;            DShapeFunction[4][1] = 0.0;               DShapeFunction[4][2] = 0.5*(Xi+1.0);
  DShapeFunction[5][0] = 0.5*(1.0-Eta-Zeta);  DShapeFunction[5][1] = -0.5*(Xi+1.0);     DShapeFunction[5][2] = -0.5*(Xi+1.0);
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 6; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to Jacobian ---*/
  
  ad[0][0] = xs[1][1]*xs[2][2]-xs[1][2]*xs[2][1];
  ad[0][1] = xs[0][2]*xs[2][1]-xs[0][1]*xs[2][2];
  ad[0][2] = xs[0][1]*xs[1][2]-xs[0][2]*xs[1][1];
  ad[1][0] = xs[1][2]*xs[2][0]-xs[1][0]*xs[2][2];
  ad[1][1] = xs[0][0]*xs[2][2]-xs[0][2]*xs[2][0];
  ad[1][2] = xs[0][2]*xs[1][0]-xs[0][0]*xs[1][2];
  ad[2][0] = xs[1][0]*xs[2][1]-xs[1][1]*xs[2][0];
  ad[2][1] = xs[0][1]*xs[2][0]-xs[0][0]*xs[2][1];
  ad[2][2] = xs[0][0]*xs[1][1]-xs[0][1]*xs[1][0];
  
  /*--- Determinant of Jacobian ---*/
  
  xsj = xs[0][0]*ad[0][0]+xs[0][1]*ad[1][0]+xs[0][2]*ad[2][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 6; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]+xs[0][2]*DShapeFunction[k][2]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]+xs[1][2]*DShapeFunction[k][2]; // dN/dy
    c2 = xs[2][0]*DShapeFunction[k][0]+xs[2][1]*DShapeFunction[k][1]+xs[2][2]*DShapeFunction[k][2]; // dN/dz
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d zeta
  }
  
  return xsj;
  
}

su2double CVolumetricMovement::ShapeFunc_Hexa(su2double Xi, su2double Eta, su2double Zeta, su2double CoordCorners[8][3], su2double DShapeFunction[8][4]) {
  
  int i, j, k;
  su2double c0, c1, c2, xsj;
  su2double xs[3][3], ad[3][3];
  
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0-Zeta);
  DShapeFunction[1][3] = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0-Zeta);
  DShapeFunction[2][3] = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0-Zeta);
  DShapeFunction[3][3] = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0-Zeta);
  DShapeFunction[4][3] = 0.125*(1.0-Xi)*(1.0-Eta)*(1.0+Zeta);
  DShapeFunction[5][3] = 0.125*(1.0+Xi)*(1.0-Eta)*(1.0+Zeta);
  DShapeFunction[6][3] = 0.125*(1.0+Xi)*(1.0+Eta)*(1.0+Zeta);
  DShapeFunction[7][3] = 0.125*(1.0-Xi)*(1.0+Eta)*(1.0+Zeta);
  
  /*--- dN/d xi ---*/
  
  DShapeFunction[0][0] = -0.125*(1.0-Eta)*(1.0-Zeta);
  DShapeFunction[1][0] = 0.125*(1.0-Eta)*(1.0-Zeta);
  DShapeFunction[2][0] = 0.125*(1.0+Eta)*(1.0-Zeta);
  DShapeFunction[3][0] = -0.125*(1.0+Eta)*(1.0-Zeta);
  DShapeFunction[4][0] = -0.125*(1.0-Eta)*(1.0+Zeta);
  DShapeFunction[5][0] = 0.125*(1.0-Eta)*(1.0+Zeta);
  DShapeFunction[6][0] = 0.125*(1.0+Eta)*(1.0+Zeta);
  DShapeFunction[7][0] = -0.125*(1.0+Eta)*(1.0+Zeta);
  
  /*--- dN/d eta ---*/
  
  DShapeFunction[0][1] = -0.125*(1.0-Xi)*(1.0-Zeta);
  DShapeFunction[1][1] = -0.125*(1.0+Xi)*(1.0-Zeta);
  DShapeFunction[2][1] = 0.125*(1.0+Xi)*(1.0-Zeta);
  DShapeFunction[3][1] = 0.125*(1.0-Xi)*(1.0-Zeta);
  DShapeFunction[4][1] = -0.125*(1.0-Xi)*(1.0+Zeta);
  DShapeFunction[5][1] = -0.125*(1.0+Xi)*(1.0+Zeta);
  DShapeFunction[6][1] = 0.125*(1.0+Xi)*(1.0+Zeta);
  DShapeFunction[7][1] = 0.125*(1.0-Xi)*(1.0+Zeta);
  
  /*--- dN/d zeta ---*/
  
  DShapeFunction[0][2] = -0.125*(1.0-Xi)*(1.0-Eta);
  DShapeFunction[1][2] = -0.125*(1.0+Xi)*(1.0-Eta);
  DShapeFunction[2][2] = -0.125*(1.0+Xi)*(1.0+Eta);
  DShapeFunction[3][2] = -0.125*(1.0-Xi)*(1.0+Eta);
  DShapeFunction[4][2] = 0.125*(1.0-Xi)*(1.0-Eta);
  DShapeFunction[5][2] = 0.125*(1.0+Xi)*(1.0-Eta);
  DShapeFunction[6][2] = 0.125*(1.0+Xi)*(1.0+Eta);
  DShapeFunction[7][2] = 0.125*(1.0-Xi)*(1.0+Eta);
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 8; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to Jacobian ---*/
  
  ad[0][0] = xs[1][1]*xs[2][2]-xs[1][2]*xs[2][1];
  ad[0][1] = xs[0][2]*xs[2][1]-xs[0][1]*xs[2][2];
  ad[0][2] = xs[0][1]*xs[1][2]-xs[0][2]*xs[1][1];
  ad[1][0] = xs[1][2]*xs[2][0]-xs[1][0]*xs[2][2];
  ad[1][1] = xs[0][0]*xs[2][2]-xs[0][2]*xs[2][0];
  ad[1][2] = xs[0][2]*xs[1][0]-xs[0][0]*xs[1][2];
  ad[2][0] = xs[1][0]*xs[2][1]-xs[1][1]*xs[2][0];
  ad[2][1] = xs[0][1]*xs[2][0]-xs[0][0]*xs[2][1];
  ad[2][2] = xs[0][0]*xs[1][1]-xs[0][1]*xs[1][0];
  
  /*--- Determinant of Jacobian ---*/
  
  xsj = xs[0][0]*ad[0][0]+xs[0][1]*ad[1][0]+xs[0][2]*ad[2][0];
  
  /*--- Jacobian inverse ---*/
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 8; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]+xs[0][2]*DShapeFunction[k][2]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]+xs[1][2]*DShapeFunction[k][2]; // dN/dy
    c2 = xs[2][0]*DShapeFunction[k][0]+xs[2][1]*DShapeFunction[k][1]+xs[2][2]*DShapeFunction[k][2]; // dN/dz
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d zeta
  }
  
  return xsj;
  
}

su2double CVolumetricMovement::GetTriangle_Area(su2double CoordCorners[8][3]) {
  
  unsigned short iDim;
  su2double a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0};
  su2double *Coord_0, *Coord_1, *Coord_2, Area;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = Coord_0[iDim]-Coord_2[iDim];
    b[iDim] = Coord_1[iDim]-Coord_2[iDim];
  }
  
  Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  
  return Area;
  
}

su2double CVolumetricMovement::GetQuadrilateral_Area(su2double CoordCorners[8][3]) {
  
  unsigned short iDim;
  su2double a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0};
  su2double *Coord_0, *Coord_1, *Coord_2, Area;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = Coord_0[iDim]-Coord_2[iDim];
    b[iDim] = Coord_1[iDim]-Coord_2[iDim];
  }
  
  Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[3];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = Coord_0[iDim]-Coord_2[iDim];
    b[iDim] = Coord_1[iDim]-Coord_2[iDim];
  }
  
  Area += 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  
  return Area;
  
}

su2double CVolumetricMovement::GetTetra_Volume(su2double CoordCorners[8][3]) {
  
  unsigned short iDim;
  su2double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0}, Volume;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  Coord_3 = CoordCorners[3];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  return Volume;
  
}

su2double CVolumetricMovement::GetPyram_Volume(su2double CoordCorners[8][3]) {
  
  unsigned short iDim;
  su2double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0}, Volume;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  Coord_3 = CoordCorners[4];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[3];
  Coord_3 = CoordCorners[4];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  return Volume;

}

su2double CVolumetricMovement::GetPrism_Volume(su2double CoordCorners[8][3]) {
  
  unsigned short iDim;
  su2double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0}, Volume;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[1];
  Coord_3 = CoordCorners[5];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
    
  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[5];
  Coord_2 = CoordCorners[1];
  Coord_3 = CoordCorners[4];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[5];
  Coord_2 = CoordCorners[4];
  Coord_3 = CoordCorners[3];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  return Volume;

}

su2double CVolumetricMovement::GetHexa_Volume(su2double CoordCorners[8][3]) {
  
  unsigned short iDim;
  su2double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
  su2double r1[3] = {0.0,0.0,0.0}, r2[3] = {0.0,0.0,0.0}, r3[3] = {0.0,0.0,0.0}, CrossProduct[3] = {0.0,0.0,0.0}, Volume;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  Coord_3 = CoordCorners[5];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[7];
  Coord_3 = CoordCorners[5];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[3];
  Coord_3 = CoordCorners[7];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[5];
  Coord_2 = CoordCorners[7];
  Coord_3 = CoordCorners[4];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[2];
  Coord_1 = CoordCorners[7];
  Coord_2 = CoordCorners[5];
  Coord_3 = CoordCorners[6];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
  CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
  CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
  CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  return Volume;

}

void CVolumetricMovement::SetFEA_StiffMatrix2D(CGeometry *geometry, CConfig *config, su2double **StiffMatrix_Elem, unsigned long PointCorners[8], su2double CoordCorners[8][3],
                                               unsigned short nNodes, su2double ElemVolume, su2double ElemDistance) {
  
  su2double B_Matrix[3][8], D_Matrix[3][3], Aux_Matrix[8][3];
  su2double Xi = 0.0, Eta = 0.0, Det = 0.0, E = 1/EPS, Lambda = 0.0, Mu = 0.0, Nu = 0.0;
  unsigned short iNode, iVar, jVar, kVar, iGauss, nGauss = 0;
  su2double DShapeFunction[8][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
  su2double Location[4][3], Weight[4];
  unsigned short nVar = geometry->GetnDim();
  
  for (iVar = 0; iVar < nNodes*nVar; iVar++) {
    for (jVar = 0; jVar < nNodes*nVar; jVar++) {
      StiffMatrix_Elem[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Integration formulae from "Shape functions and points of
   integration of the Résumé" by Josselin DELMAS (2013) ---*/
  
  /*--- Triangle. Nodes of numerical integration at 1 point (order 1). ---*/
  
  if (nNodes == 3) {
    nGauss = 1;
    Location[0][0] = 0.333333333333333;  Location[0][1] = 0.333333333333333;  Weight[0] = 0.5;
  }
  
  /*--- Quadrilateral. Nodes of numerical integration at 4 points (order 2). ---*/
  
  if (nNodes == 4) {
    nGauss = 4;
    Location[0][0] = -0.577350269189626;  Location[0][1] = -0.577350269189626;  Weight[0] = 1.0;
    Location[1][0] = 0.577350269189626;   Location[1][1] = -0.577350269189626;  Weight[1] = 1.0;
    Location[2][0] = 0.577350269189626;   Location[2][1] = 0.577350269189626;   Weight[2] = 1.0;
    Location[3][0] = -0.577350269189626;  Location[3][1] = 0.577350269189626;   Weight[3] = 1.0;
  }
  
  for (iGauss = 0; iGauss < nGauss; iGauss++) {
    
    Xi = Location[iGauss][0]; Eta = Location[iGauss][1];
    
    if (nNodes == 3) Det = ShapeFunc_Triangle(Xi, Eta, CoordCorners, DShapeFunction);
    if (nNodes == 4) Det = ShapeFunc_Quadrilateral(Xi, Eta, CoordCorners, DShapeFunction);
    
    /*--- Compute the B Matrix ---*/
    
    for (iVar = 0; iVar < 3; iVar++)
      for (jVar = 0; jVar < nNodes*nVar; jVar++)
        B_Matrix[iVar][jVar] = 0.0;
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      B_Matrix[0][0+iNode*nVar] = DShapeFunction[iNode][0];
      B_Matrix[1][1+iNode*nVar] = DShapeFunction[iNode][1];
      
      B_Matrix[2][0+iNode*nVar] = DShapeFunction[iNode][1];
      B_Matrix[2][1+iNode*nVar] = DShapeFunction[iNode][0];
    }
    
    /*--- Impose a type of stiffness for each element ---*/
    
    switch (config->GetDeform_Stiffness_Type()) {
      case INVERSE_VOLUME: E = 1.0 / ElemVolume; break;
      case SOLID_WALL_DISTANCE: E = 1.0 / ElemDistance; break;
      case CONSTANT_STIFFNESS: E = 1.0 / EPS; break;
    }
    
    Nu = config->GetDeform_Coeff();
    Mu = E / (2.0*(1.0 + Nu));
    Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));
    
    /*--- Compute the D Matrix (for plane strain and 3-D)---*/
    
    D_Matrix[0][0] = Lambda + 2.0*Mu;    D_Matrix[0][1] = Lambda;            D_Matrix[0][2] = 0.0;
    D_Matrix[1][0] = Lambda;            D_Matrix[1][1] = Lambda + 2.0*Mu;   D_Matrix[1][2] = 0.0;
    D_Matrix[2][0] = 0.0;               D_Matrix[2][1] = 0.0;               D_Matrix[2][2] = Mu;
    
    
    /*--- Compute the BT.D Matrix ---*/
    
    for (iVar = 0; iVar < nNodes*nVar; iVar++) {
      for (jVar = 0; jVar < 3; jVar++) {
        Aux_Matrix[iVar][jVar] = 0.0;
        for (kVar = 0; kVar < 3; kVar++)
          Aux_Matrix[iVar][jVar] += B_Matrix[kVar][iVar]*D_Matrix[kVar][jVar];
      }
    }
    
    /*--- Compute the BT.D.B Matrix (stiffness matrix), and add to the original
     matrix using Gauss integration ---*/
    
    for (iVar = 0; iVar < nNodes*nVar; iVar++) {
      for (jVar = 0; jVar < nNodes*nVar; jVar++) {
        for (kVar = 0; kVar < 3; kVar++) {
          StiffMatrix_Elem[iVar][jVar] += Weight[iGauss] * Aux_Matrix[iVar][kVar]*B_Matrix[kVar][jVar] * fabs(Det);
        }
      }
    }
    
  }
  
}

void CVolumetricMovement::SetFEA_StiffMatrix3D(CGeometry *geometry, CConfig *config, su2double **StiffMatrix_Elem, unsigned long PointCorners[8], su2double CoordCorners[8][3],
                                               unsigned short nNodes, su2double ElemVolume, su2double ElemDistance) {
  
  su2double B_Matrix[6][24], D_Matrix[6][6] = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}, Aux_Matrix[24][6];
  su2double Xi = 0.0, Eta = 0.0, Zeta = 0.0, Det = 0.0, Mu = 0.0, E = 0.0, Lambda = 0.0, Nu = 0.0;
  unsigned short iNode, iVar, jVar, kVar, iGauss, nGauss = 0;
  su2double DShapeFunction[8][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
  su2double Location[8][3], Weight[8];
  unsigned short nVar = geometry->GetnDim();
  
  for (iVar = 0; iVar < nNodes*nVar; iVar++) {
    for (jVar = 0; jVar < nNodes*nVar; jVar++) {
      StiffMatrix_Elem[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Integration formulae from "Shape functions and points of
   integration of the Résumé" by Josselin Delmas (2013) ---*/
  
  /*--- Tetrahedrons. Nodes of numerical integration at 1 point (order 1). ---*/
  
  if (nNodes == 4) {
    nGauss = 1;
    Location[0][0] = 0.25;  Location[0][1] = 0.25;  Location[0][2] = 0.25;  Weight[0] = 0.166666666666666;
  }
  
  /*--- Pyramids. Nodes numerical integration at 5 points. ---*/
  
  if (nNodes == 5) {
    nGauss = 5;
    Location[0][0] = 0.5;   Location[0][1] = 0.0;   Location[0][2] = 0.1531754163448146;  Weight[0] = 0.133333333333333;
    Location[1][0] = 0.0;   Location[1][1] = 0.5;   Location[1][2] = 0.1531754163448146;  Weight[1] = 0.133333333333333;
    Location[2][0] = -0.5;  Location[2][1] = 0.0;   Location[2][2] = 0.1531754163448146;  Weight[2] = 0.133333333333333;
    Location[3][0] = 0.0;   Location[3][1] = -0.5;  Location[3][2] = 0.1531754163448146;  Weight[3] = 0.133333333333333;
    Location[4][0] = 0.0;   Location[4][1] = 0.0;   Location[4][2] = 0.6372983346207416;  Weight[4] = 0.133333333333333;
  }
  
  /*--- Prism. Nodes of numerical integration at 6 points (order 3 in Xi, order 2 in Eta and Mu ). ---*/
  
  if (nNodes == 6) {
    nGauss = 6;
    Location[0][0] = -0.577350269189626;  Location[0][1] = 0.166666666666667;  Location[0][2] = 0.166666666666667;  Weight[0] = 0.166666666666667;
    Location[1][0] = -0.577350269189626;  Location[1][1] = 0.666666666666667;  Location[1][2] = 0.166666666666667;  Weight[1] = 0.166666666666667;
    Location[2][0] = -0.577350269189626;  Location[2][1] = 0.166666666666667;  Location[2][2] = 0.666666666666667;  Weight[2] = 0.166666666666667;
    Location[3][0] =  0.577350269189626;  Location[3][1] = 0.166666666666667;  Location[3][2] = 0.166666666666667;  Weight[3] = 0.166666666666667;
    Location[4][0] =  0.577350269189626;  Location[4][1] = 0.666666666666667;  Location[4][2] = 0.166666666666667;  Weight[4] = 0.166666666666667;
    Location[5][0] =  0.577350269189626;  Location[5][1] = 0.166666666666667;  Location[5][2] = 0.666666666666667;  Weight[5] = 0.166666666666667;
  }
  
  /*--- Hexahedrons. Nodes of numerical integration at 6 points (order 3). ---*/
  
  if (nNodes == 8) {
    nGauss = 8;
    Location[0][0] = -0.577350269189626;  Location[0][1] = -0.577350269189626;  Location[0][2] = -0.577350269189626;  Weight[0] = 1.0;
    Location[1][0] = -0.577350269189626;  Location[1][1] = -0.577350269189626;  Location[1][2] = 0.577350269189626;   Weight[1] = 1.0;
    Location[2][0] = -0.577350269189626;  Location[2][1] = 0.577350269189626;   Location[2][2] = -0.577350269189626;  Weight[2] = 1.0;
    Location[3][0] = -0.577350269189626;  Location[3][1] = 0.577350269189626;   Location[3][2] = 0.577350269189626;   Weight[3] = 1.0;
    Location[4][0] = 0.577350269189626;   Location[4][1] = -0.577350269189626;  Location[4][2] = -0.577350269189626;  Weight[4] = 1.0;
    Location[5][0] = 0.577350269189626;   Location[5][1] = -0.577350269189626;  Location[5][2] = 0.577350269189626;   Weight[5] = 1.0;
    Location[6][0] = 0.577350269189626;   Location[6][1] = 0.577350269189626;   Location[6][2] = -0.577350269189626;  Weight[6] = 1.0;
    Location[7][0] = 0.577350269189626;   Location[7][1] = 0.577350269189626;   Location[7][2] = 0.577350269189626;   Weight[7] = 1.0;
  }
  
  for (iGauss = 0; iGauss < nGauss; iGauss++) {
    
    Xi = Location[iGauss][0]; Eta = Location[iGauss][1];  Zeta = Location[iGauss][2];
    
    if (nNodes == 4) Det = ShapeFunc_Tetra(Xi, Eta, Zeta, CoordCorners, DShapeFunction);
    if (nNodes == 5) Det = ShapeFunc_Pyram(Xi, Eta, Zeta, CoordCorners, DShapeFunction);
    if (nNodes == 6) Det = ShapeFunc_Prism(Xi, Eta, Zeta, CoordCorners, DShapeFunction);
    if (nNodes == 8) Det = ShapeFunc_Hexa(Xi, Eta, Zeta, CoordCorners, DShapeFunction);
    
    /*--- Compute the B Matrix ---*/
    
    for (iVar = 0; iVar < 6; iVar++)
      for (jVar = 0; jVar < nNodes*nVar; jVar++)
        B_Matrix[iVar][jVar] = 0.0;
    
    for (iNode = 0; iNode < nNodes; iNode++) {
      B_Matrix[0][0+iNode*nVar] = DShapeFunction[iNode][0];
      B_Matrix[1][1+iNode*nVar] = DShapeFunction[iNode][1];
      B_Matrix[2][2+iNode*nVar] = DShapeFunction[iNode][2];
      
      B_Matrix[3][0+iNode*nVar] = DShapeFunction[iNode][1];
      B_Matrix[3][1+iNode*nVar] = DShapeFunction[iNode][0];
      
      B_Matrix[4][1+iNode*nVar] = DShapeFunction[iNode][2];
      B_Matrix[4][2+iNode*nVar] = DShapeFunction[iNode][1];
      
      B_Matrix[5][0+iNode*nVar] = DShapeFunction[iNode][2];
      B_Matrix[5][2+iNode*nVar] = DShapeFunction[iNode][0];
    }
    
    /*--- Impose a type of stiffness for each element ---*/
    
    switch (config->GetDeform_Stiffness_Type()) {
      case INVERSE_VOLUME: E = 1.0 / ElemVolume; break;
      case SOLID_WALL_DISTANCE: E = 1.0 / ElemDistance; break;
      case CONSTANT_STIFFNESS: E = 1.0 / EPS; break;
    }
    
    Nu = config->GetDeform_Coeff();
    Mu = E / (2.0*(1.0 + Nu));
    Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));
    
    /*--- Compute the D Matrix (for plane strain and 3-D)---*/
    
    D_Matrix[0][0] = Lambda + 2.0*Mu;  D_Matrix[0][1] = Lambda;          D_Matrix[0][2] = Lambda;
    D_Matrix[1][0] = Lambda;          D_Matrix[1][1] = Lambda + 2.0*Mu;  D_Matrix[1][2] = Lambda;
    D_Matrix[2][0] = Lambda;          D_Matrix[2][1] = Lambda;          D_Matrix[2][2] = Lambda + 2.0*Mu;
    D_Matrix[3][3] = Mu;
    D_Matrix[4][4] = Mu;
    D_Matrix[5][5] = Mu;
    
    
    /*--- Compute the BT.D Matrix ---*/
    
    for (iVar = 0; iVar < nNodes*nVar; iVar++) {
      for (jVar = 0; jVar < 6; jVar++) {
        Aux_Matrix[iVar][jVar] = 0.0;
        for (kVar = 0; kVar < 6; kVar++)
          Aux_Matrix[iVar][jVar] += B_Matrix[kVar][iVar]*D_Matrix[kVar][jVar];
      }
    }
    
    /*--- Compute the BT.D.B Matrix (stiffness matrix), and add to the original
     matrix using Gauss integration ---*/
    
    for (iVar = 0; iVar < nNodes*nVar; iVar++) {
      for (jVar = 0; jVar < nNodes*nVar; jVar++) {
        for (kVar = 0; kVar < 6; kVar++) {
          StiffMatrix_Elem[iVar][jVar] += Weight[iGauss] * Aux_Matrix[iVar][kVar]*B_Matrix[kVar][jVar] * fabs(Det);
        }
      }
    }
    
  }
  
}

void CVolumetricMovement::AddFEA_StiffMatrix(CGeometry *geometry, su2double **StiffMatrix_Elem, unsigned long PointCorners[8], unsigned short nNodes) {
  
  unsigned short iVar, jVar, iDim, jDim;
  
  unsigned short nVar = geometry->GetnDim();

  su2double **StiffMatrix_Node;
  StiffMatrix_Node = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    StiffMatrix_Node[iVar] = new su2double [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      StiffMatrix_Node[iVar][jVar] = 0.0;
  
  /*--- Transform the stiffness matrix for the hexahedral element into the
   contributions for the individual nodes relative to each other. ---*/
  
  for (iVar = 0; iVar < nNodes; iVar++) {
    for (jVar = 0; jVar < nNodes; jVar++) {
      
      for (iDim = 0; iDim < nVar; iDim++) {
        for (jDim = 0; jDim < nVar; jDim++) {
          StiffMatrix_Node[iDim][jDim] = StiffMatrix_Elem[(iVar*nVar)+iDim][(jVar*nVar)+jDim];
        }
      }

      StiffMatrix.AddBlock(PointCorners[iVar], PointCorners[jVar], StiffMatrix_Node);
      
    }
  }
  
  /*--- Deallocate memory and exit ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] StiffMatrix_Node[iVar];
  delete [] StiffMatrix_Node;
  
}

void CVolumetricMovement::SetBoundaryDisplacements(CGeometry *geometry, CConfig *config) {
  
  unsigned short iDim, nDim = geometry->GetnDim(), iMarker, axis = 0;
  unsigned long iPoint, total_index, iVertex;
  su2double *VarCoord, MeanCoord[3] = {0.0,0.0,0.0}, VarIncrement = 1.0;
  
  /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
   deforming meshes (MARKER_MOVING), while SU2_DEF will use it for deforming
   meshes after imposing design variable surface deformations (DV_MARKER). ---*/
  
  unsigned short Kind_SU2 = config->GetKind_SU2();
  
  /*--- If requested (no by default) impose the surface deflections in
   increments and solve the grid deformation equations iteratively with
   successive small deformations. ---*/
  
  VarIncrement = 1.0/((su2double)config->GetGridDef_Nonlinear_Iter());
  
  /*--- As initialization, set to zero displacements of all the surfaces except the symmetry
   plane, internal and periodic bc the receive boundaries and periodic boundaries. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (((config->GetMarker_All_KindBC(iMarker) != SYMMETRY_PLANE) &&
         (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) &&
         (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
         (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY))) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        for (iDim = 0; iDim < nDim; iDim++) {
          total_index = iPoint*nDim + iDim;
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
    if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_CFD)) ||
        ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_DEF)) ||
        ((config->GetDirectDiff() == D_DESIGN) && (Kind_SU2 == SU2_CFD) && (config->GetMarker_All_DV(iMarker) == YES)) ||
        ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_DOT))) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
        for (iDim = 0; iDim < nDim; iDim++) {
          total_index = iPoint*nDim + iDim;
          LinSysRes[total_index] = SU2_TYPE::GetValue(VarCoord[iDim] * VarIncrement);
          LinSysSol[total_index] = SU2_TYPE::GetValue(VarCoord[iDim] * VarIncrement);
          StiffMatrix.DeleteValsRowi(total_index);
        }
      }
    }
  }
  
  /*--- Set to zero displacements of the normal component for the symmetry plane condition ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) ) {
      
      su2double *Coord_0 = NULL;
      for (iDim = 0; iDim < nDim; iDim++) MeanCoord[iDim] = 0.0;
      
      /*--- Store the coord of the first point to help identify the axis. ---*/
      
      iPoint  = geometry->vertex[iMarker][0]->GetNode();
      Coord_0 = geometry->node[iPoint]->GetCoord();
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        VarCoord = geometry->node[iPoint]->GetCoord();
        for (iDim = 0; iDim < nDim; iDim++)
          MeanCoord[iDim] += (VarCoord[iDim]-Coord_0[iDim])*(VarCoord[iDim]-Coord_0[iDim]);
      }
      for (iDim = 0; iDim < nDim; iDim++) MeanCoord[iDim] = sqrt(MeanCoord[iDim]);
      if (nDim==3) {
        if ((MeanCoord[0] <= MeanCoord[1]) && (MeanCoord[0] <= MeanCoord[2])) axis = 0;
        if ((MeanCoord[1] <= MeanCoord[0]) && (MeanCoord[1] <= MeanCoord[2])) axis = 1;
        if ((MeanCoord[2] <= MeanCoord[0]) && (MeanCoord[2] <= MeanCoord[1])) axis = 2;
      }
      else {
        if ((MeanCoord[0] <= MeanCoord[1]) ) axis = 0;
        if ((MeanCoord[1] <= MeanCoord[0]) ) axis = 1;
      }
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        total_index = iPoint*nDim + axis;
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
          total_index = iPoint*nDim + iDim;
          LinSysRes[total_index] = 0.0;
          LinSysSol[total_index] = 0.0;
          StiffMatrix.DeleteValsRowi(total_index);
        }
      }
    }
  }

  /*--- Move the FSI interfaces ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_ZoneInterface(iMarker) != 0) && (Kind_SU2 == SU2_CFD)) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
        for (iDim = 0; iDim < nDim; iDim++) {
          total_index = iPoint*nDim + iDim;
          LinSysRes[total_index] = SU2_TYPE::GetValue(VarCoord[iDim] * VarIncrement);
          LinSysSol[total_index] = SU2_TYPE::GetValue(VarCoord[iDim] * VarIncrement);
          StiffMatrix.DeleteValsRowi(total_index);
        }
      }
    }
  }

}

void CVolumetricMovement::SetBoundaryDerivatives(CGeometry *geometry, CConfig *config) {
  unsigned short iDim, iMarker;
  unsigned long iPoint, total_index, iVertex;

  su2double * VarCoord;
  unsigned short Kind_SU2 = config->GetKind_SU2();
  if ((config->GetDirectDiff() == D_DESIGN) && (Kind_SU2 == SU2_CFD)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_DV(iMarker) == YES)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
          for (iDim = 0; iDim < nDim; iDim++) {
            total_index = iPoint*nDim + iDim;
            LinSysRes[total_index] = SU2_TYPE::GetDerivative(VarCoord[iDim]);
            LinSysSol[total_index] = SU2_TYPE::GetDerivative(VarCoord[iDim]);
          }
        }
      }
    }
    if (LinSysRes.norm() == 0.0) cout << "Warning: Derivatives are zero!" << endl;
  } else if (Kind_SU2 == SU2_DOT) {

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        total_index = iPoint*nDim + iDim;
        LinSysRes[total_index] = SU2_TYPE::GetValue(geometry->GetSensitivity(iPoint, iDim));
        LinSysSol[total_index] = SU2_TYPE::GetValue(geometry->GetSensitivity(iPoint, iDim));
      }
    }
  }
}

void CVolumetricMovement::UpdateGridCoord_Derivatives(CGeometry *geometry, CConfig *config) {
  unsigned short iDim, iMarker;
  unsigned long iPoint, total_index, iVertex;
  su2double *new_coord = new su2double[3];

  unsigned short Kind_SU2 = config->GetKind_SU2();

  /*--- Update derivatives of the grid coordinates using the solution of the linear system
     after grid deformation (LinSysSol contains the derivatives of the x, y, z displacements). ---*/
  if ((config->GetDirectDiff() == D_DESIGN) && (Kind_SU2 == SU2_CFD)) {
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      new_coord[0] = 0.0; new_coord[1] = 0.0; new_coord[2] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        total_index = iPoint*nDim + iDim;
        new_coord[iDim] = geometry->node[iPoint]->GetCoord(iDim);
        SU2_TYPE::SetDerivative(new_coord[iDim], SU2_TYPE::GetValue(LinSysSol[total_index]));
      }
      geometry->node[iPoint]->SetCoord(new_coord);
    }
  } else if (Kind_SU2 == SU2_DOT) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX ) ||
         (config->GetMarker_All_KindBC(iMarker) == EULER_WALL ) ||
         (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL ) ||
         (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE) ||
         (config->GetMarker_All_DV(iMarker) == YES)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (geometry->node[iPoint]->GetDomain()) {
            for (iDim = 0; iDim < nDim; iDim++) {
              total_index = iPoint*nDim + iDim;
              geometry->SetSensitivity(iPoint,iDim, LinSysSol[total_index]);
            }
          }
        }
      }
    }
  }
  
  delete [] new_coord;
}

void CVolumetricMovement::SetDomainDisplacements(CGeometry *geometry, CConfig *config) {
  
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, total_index;
  su2double *Coord, *MinCoordValues, *MaxCoordValues, *Hold_GridFixed_Coord;
  
  if (config->GetHold_GridFixed()) {
    
    MinCoordValues = new su2double [nDim];
    MaxCoordValues = new su2double [nDim];
    
    for (iDim = 0; iDim < nDim; iDim++) {
      MinCoordValues[iDim] = 0.0;
      MaxCoordValues[iDim] = 0.0;
    }
    
    Hold_GridFixed_Coord = config->GetHold_GridFixed_Coord();
    
    MinCoordValues[0] = Hold_GridFixed_Coord[0];
    MinCoordValues[1] = Hold_GridFixed_Coord[1];
    MinCoordValues[2] = Hold_GridFixed_Coord[2];
    MaxCoordValues[0] = Hold_GridFixed_Coord[3];
    MaxCoordValues[1] = Hold_GridFixed_Coord[4];
    MaxCoordValues[2] = Hold_GridFixed_Coord[5];
    
    /*--- Set to zero displacements of all the points that are not going to be moved
     except the surfaces ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      Coord = geometry->node[iPoint]->GetCoord();
      for (iDim = 0; iDim < nDim; iDim++) {
        if ((Coord[iDim] < MinCoordValues[iDim]) || (Coord[iDim] > MaxCoordValues[iDim])) {
          total_index = iPoint*nDim + iDim;
          LinSysRes[total_index] = 0.0;
          LinSysSol[total_index] = 0.0;
          StiffMatrix.DeleteValsRowi(total_index);
        }
      }
    }
    
    delete [] MinCoordValues;
    delete [] MaxCoordValues;
    
  }
  
  /*--- Don't move the volume grid outside the limits based 
   on the distance to the solid surface ---*/
  
  if (config->GetDeform_Limit() < 1E6) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      if (geometry->node[iPoint]->GetWall_Distance() >= config->GetDeform_Limit()) {
        for (iDim = 0; iDim < nDim; iDim++) {
          total_index = iPoint*nDim + iDim;
          LinSysRes[total_index] = 0.0;
          LinSysSol[total_index] = 0.0;
          StiffMatrix.DeleteValsRowi(total_index);
        }
      }
    }
  }
  
}

void CVolumetricMovement::Rigid_Rotation(CGeometry *geometry, CConfig *config,
                                         unsigned short iZone, unsigned long iter) {
  
  /*--- Local variables ---*/
  unsigned short iDim, nDim; 
  unsigned long iPoint;
  su2double r[3] = {0.0,0.0,0.0}, rotCoord[3] = {0.0,0.0,0.0}, *Coord;
  su2double Center[3] = {0.0,0.0,0.0}, Omega[3] = {0.0,0.0,0.0}, Lref;
  su2double dt, Center_Moment[3] = {0.0,0.0,0.0};
  su2double *GridVel, newGridVel[3] = {0.0,0.0,0.0};
  su2double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  su2double dtheta, dphi, dpsi, cosTheta, sinTheta;
  su2double cosPhi, sinPhi, cosPsi, sinPsi;
  bool harmonic_balance = (config->GetTime_Marching() == HARMONIC_BALANCE);
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());

  /*--- Problem dimension and physical time step ---*/
  nDim = geometry->GetnDim();
  dt   = config->GetDelta_UnstTimeND();
  Lref = config->GetLength_Ref();

  /*--- For the unsteady adjoint, use reverse time ---*/
  if (adjoint) {
    /*--- Set the first adjoint mesh position to the final direct one ---*/
    if (iter == 0) dt = ((su2double)config->GetnTime_Iter()-1)*dt;
    /*--- Reverse the rotation direction for the adjoint ---*/
    else dt = -1.0*dt;
  } else {
    /*--- No rotation at all for the first direct solution ---*/
    if (iter == 0) dt = 0;
  }
  
  /*--- Center of rotation & angular velocity vector from config ---*/
  
  for (iDim = 0; iDim < 3; iDim++){
    Center[iDim] = config->GetMotion_Origin(iDim);
    Omega[iDim]  = config->GetRotation_Rate(iDim)/config->GetOmega_Ref();
  }

  /*-- Set dt for harmonic balance cases ---*/
  if (harmonic_balance) {
    /*--- period of oscillation & compute time interval using nTimeInstances ---*/
    su2double period = config->GetHarmonicBalance_Period();
    period /= config->GetTime_Ref();
    dt = period * (su2double)iter/(su2double)(config->GetnTimeInstances());
  }
  
  /*--- Compute delta change in the angle about the x, y, & z axes. ---*/

  dtheta = Omega[0]*dt;
  dphi   = Omega[1]*dt;
  dpsi   = Omega[2]*dt;

  if (rank == MASTER_NODE && iter == 0) {
    cout << " Angular velocity: (" << Omega[0] << ", " << Omega[1];
    cout << ", " << Omega[2] << ") rad/s." << endl;
  }
  
  /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
  
  cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
  sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
  
  /*--- Compute the rotation matrix. Note that the implicit
   ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
  
  rotMatrix[0][0] = cosPhi*cosPsi;
  rotMatrix[1][0] = cosPhi*sinPsi;
  rotMatrix[2][0] = -sinPhi;
  
  rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
  rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
  rotMatrix[2][1] = sinTheta*cosPhi;
  
  rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
  rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
  rotMatrix[2][2] = cosTheta*cosPhi;
  
  /*--- Loop over and rotate each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord   = geometry->node[iPoint]->GetCoord();
    GridVel = geometry->node[iPoint]->GetGridVel();
    
    /*--- Calculate non-dim. position from rotation center ---*/
    r[0] = (Coord[0]-Center[0])/Lref;
    r[1] = (Coord[1]-Center[1])/Lref;
    if (nDim == 3) r[2] = (Coord[2]-Center[2])/Lref;
    
    /*--- Compute transformed point coordinates ---*/
    rotCoord[0] = rotMatrix[0][0]*r[0] 
                + rotMatrix[0][1]*r[1] 
                + rotMatrix[0][2]*r[2];
    
    rotCoord[1] = rotMatrix[1][0]*r[0] 
                + rotMatrix[1][1]*r[1] 
                + rotMatrix[1][2]*r[2];
    
    rotCoord[2] = rotMatrix[2][0]*r[0] 
                + rotMatrix[2][1]*r[1] 
                + rotMatrix[2][2]*r[2];
    
    /*--- Cross Product of angular velocity and distance from center.
     Note that we have assumed the grid velocities have been set to
     an initial value in the plunging routine. ---*/
    
    newGridVel[0] = GridVel[0] + Omega[1]*rotCoord[2] - Omega[2]*rotCoord[1];
    newGridVel[1] = GridVel[1] + Omega[2]*rotCoord[0] - Omega[0]*rotCoord[2];
    if (nDim == 3) newGridVel[2] = GridVel[2] + Omega[0]*rotCoord[1] - Omega[1]*rotCoord[0];
    
    /*--- Store new node location & grid velocity. Add center. 
     Do not store the grid velocity if this is an adjoint calculation.---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim, rotCoord[iDim] + Center[iDim]);
      if (!adjoint) geometry->node[iPoint]->SetGridVel(iDim, newGridVel[iDim]);
      
    }
  }
  
  /*--- Set the moment computation center to the new location after
   incrementing the position with the rotation. ---*/
  
  for (unsigned short jMarker=0; jMarker<config->GetnMarker_Monitoring(); jMarker++) {
    
    Center_Moment[0] = config->GetRefOriginMoment_X(jMarker);
    Center_Moment[1] = config->GetRefOriginMoment_Y(jMarker);
    Center_Moment[2] = config->GetRefOriginMoment_Z(jMarker);
    
    /*--- Calculate non-dim. position from rotation center ---*/
    
    for (iDim = 0; iDim < nDim; iDim++)
      r[iDim] = (Center_Moment[iDim]-Center[iDim])/Lref;
    if (nDim == 2) r[nDim] = 0.0;
    
    /*--- Compute transformed point coordinates ---*/
    
    rotCoord[0] = rotMatrix[0][0]*r[0]
    + rotMatrix[0][1]*r[1]
    + rotMatrix[0][2]*r[2];
    
    rotCoord[1] = rotMatrix[1][0]*r[0]
    + rotMatrix[1][1]*r[1]
    + rotMatrix[1][2]*r[2];
    
    rotCoord[2] = rotMatrix[2][0]*r[0]
    + rotMatrix[2][1]*r[1]
    + rotMatrix[2][2]*r[2];
    
    config->SetRefOriginMoment_X(jMarker, Center[0]+rotCoord[0]);
    config->SetRefOriginMoment_Y(jMarker, Center[1]+rotCoord[1]);
    config->SetRefOriginMoment_Z(jMarker, Center[2]+rotCoord[2]);
  }
  
  /*--- After moving all nodes, update geometry class ---*/
  
  UpdateDualGrid(geometry, config);

}

void CVolumetricMovement::Rigid_Pitching(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter) {
  
  /*--- Local variables ---*/
  su2double r[3] = {0.0,0.0,0.0}, rotCoord[3] = {0.0,0.0,0.0}, *Coord, Center[3] = {0.0,0.0,0.0},
  Omega[3] = {0.0,0.0,0.0}, Ampl[3] = {0.0,0.0,0.0}, Phase[3] = {0.0,0.0,0.0};
  su2double Lref, deltaT, alphaDot[3], *GridVel, newGridVel[3] = {0.0,0.0,0.0};
  su2double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  su2double dtheta, dphi, dpsi, cosTheta, sinTheta;
  su2double cosPhi, sinPhi, cosPsi, sinPsi;
  su2double time_new, time_old;
  su2double DEG2RAD = PI_NUMBER/180.0;
  unsigned short iDim;
  unsigned short nDim = geometry->GetnDim();
  unsigned long iPoint;
  bool harmonic_balance = (config->GetTime_Marching() == HARMONIC_BALANCE);
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());
  
  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND(); 
  Lref   = config->GetLength_Ref();

  /*--- Pitching origin, frequency, and amplitude from config. ---*/	
  
  for (iDim = 0; iDim < 3; iDim++){
    Center[iDim] = config->GetMotion_Origin(iDim);
    Omega[iDim]  = config->GetPitching_Omega(iDim)/config->GetOmega_Ref();
    Ampl[iDim]   = config->GetPitching_Ampl(iDim)*DEG2RAD;
    Phase[iDim]  = config->GetPitching_Phase(iDim)*DEG2RAD;
  }


  if (harmonic_balance) {    
    /*--- period of oscillation & compute time interval using nTimeInstances ---*/
    su2double period = config->GetHarmonicBalance_Period();
    period /= config->GetTime_Ref();
    deltaT = period/(su2double)(config->GetnTimeInstances());
  }

  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/ 
    unsigned long nFlowIter  = config->GetnTime_Iter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<su2double>(directIter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<su2double>(directIter)+1.0)*deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<su2double>(iter)*deltaT;
    if (harmonic_balance) {
      /*--- For harmonic balance, begin movement from the zero position ---*/
      time_old = 0.0;
    } else {
      time_old = time_new;
      if (iter != 0) time_old = (static_cast<su2double>(iter)-1.0)*deltaT;
    }
  }
  
  /*--- Compute delta change in the angle about the x, y, & z axes. ---*/
  
  dtheta = -Ampl[0]*(sin(Omega[0]*time_new + Phase[0]) - sin(Omega[0]*time_old + Phase[0]));
  dphi   = -Ampl[1]*(sin(Omega[1]*time_new + Phase[1]) - sin(Omega[1]*time_old + Phase[1]));
  dpsi   = -Ampl[2]*(sin(Omega[2]*time_new + Phase[2]) - sin(Omega[2]*time_old + Phase[2]));
  
  /*--- Angular velocity at the new time ---*/
  
  alphaDot[0] = -Omega[0]*Ampl[0]*cos(Omega[0]*time_new);
  alphaDot[1] = -Omega[1]*Ampl[1]*cos(Omega[1]*time_new);
  alphaDot[2] = -Omega[2]*Ampl[2]*cos(Omega[2]*time_new);

  if (rank == MASTER_NODE && iter == 0) {
      cout << " Pitching frequency: (" << Omega[0] << ", " << Omega[1];
      cout << ", " << Omega[2] << ") rad/s." << endl;
      cout << " Pitching amplitude: (" << Ampl[0]/DEG2RAD << ", ";
      cout << Ampl[1]/DEG2RAD << ", " << Ampl[2]/DEG2RAD;
      cout << ") degrees."<< endl;
      cout << " Pitching phase lag: (" << Phase[0]/DEG2RAD << ", ";
      cout << Phase[1]/DEG2RAD <<", "<< Phase[2]/DEG2RAD;
      cout << ") degrees."<< endl;
  }
  
  /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
  
  cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
  sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
  
  /*--- Compute the rotation matrix. Note that the implicit
   ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
  
  rotMatrix[0][0] = cosPhi*cosPsi;
  rotMatrix[1][0] = cosPhi*sinPsi;
  rotMatrix[2][0] = -sinPhi;
  
  rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
  rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
  rotMatrix[2][1] = sinTheta*cosPhi;
  
  rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
  rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
  rotMatrix[2][2] = cosTheta*cosPhi;
  
  /*--- Loop over and rotate each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord   = geometry->node[iPoint]->GetCoord();
    GridVel = geometry->node[iPoint]->GetGridVel();
    
    /*--- Calculate non-dim. position from rotation center ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
    if (nDim == 2) r[nDim] = 0.0;
    
    /*--- Compute transformed point coordinates ---*/
    rotCoord[0] = rotMatrix[0][0]*r[0] 
                + rotMatrix[0][1]*r[1] 
                + rotMatrix[0][2]*r[2];
    
    rotCoord[1] = rotMatrix[1][0]*r[0] 
                + rotMatrix[1][1]*r[1] 
                + rotMatrix[1][2]*r[2];
    
    rotCoord[2] = rotMatrix[2][0]*r[0] 
                + rotMatrix[2][1]*r[1] 
                + rotMatrix[2][2]*r[2];
    
    /*--- Cross Product of angular velocity and distance from center.
     Note that we have assumed the grid velocities have been set to 
     an initial value in the plunging routine. ---*/
    
    newGridVel[0] = GridVel[0] + alphaDot[1]*rotCoord[2] - alphaDot[2]*rotCoord[1];
    newGridVel[1] = GridVel[1] + alphaDot[2]*rotCoord[0] - alphaDot[0]*rotCoord[2];
    if (nDim == 3) newGridVel[2] = GridVel[2] + alphaDot[0]*rotCoord[1] - alphaDot[1]*rotCoord[0];
    
    /*--- Store new node location & grid velocity. Add center location.
     Do not store the grid velocity if this is an adjoint calculation.---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim, rotCoord[iDim]+Center[iDim]);
      if (!adjoint) geometry->node[iPoint]->SetGridVel(iDim, newGridVel[iDim]);
    }
  }
  
  /*--- For pitching we don't update the motion origin and moment reference origin. ---*/

  /*--- After moving all nodes, update geometry class ---*/
  
  UpdateDualGrid(geometry, config);
  
}

void CVolumetricMovement::Rigid_Plunging(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter) {
  
  /*--- Local variables ---*/
  su2double deltaX[3], newCoord[3] = {0.0, 0.0, 0.0}, Center[3], *Coord, Omega[3], Ampl[3], Lref;
  su2double *GridVel, newGridVel[3] = {0.0, 0.0, 0.0}, xDot[3];
  su2double deltaT, time_new, time_old;
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint;
  bool harmonic_balance = (config->GetTime_Marching() == HARMONIC_BALANCE);
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());
  
  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND();
  Lref   = config->GetLength_Ref();
  
  for (iDim = 0; iDim < 3; iDim++){
    Center[iDim] = config->GetMotion_Origin(iDim);
    Omega[iDim]  = config->GetPlunging_Omega(iDim)/config->GetOmega_Ref();
    Ampl[iDim]   = config->GetPlunging_Ampl(iDim)/Lref;
  }
  
  /*--- Plunging frequency and amplitude from config. ---*/
  
  if (harmonic_balance) {
    /*--- period of oscillation & time interval using nTimeInstances ---*/
    su2double period = config->GetHarmonicBalance_Period();
    period /= config->GetTime_Ref();
    deltaT = period/(su2double)(config->GetnTimeInstances());
  }
  
  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/
    unsigned long nFlowIter  = config->GetnTime_Iter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<su2double>(directIter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<su2double>(directIter)+1.0)*deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<su2double>(iter)*deltaT;
    if (harmonic_balance) {
      /*--- For harmonic balance, begin movement from the zero position ---*/
      time_old = 0.0;
    } else {
      time_old = time_new;
      if (iter != 0) time_old = (static_cast<su2double>(iter)-1.0)*deltaT;
    }
  }
  
  /*--- Compute delta change in the position in the x, y, & z directions. ---*/
  deltaX[0] = -Ampl[0]*(sin(Omega[0]*time_new) - sin(Omega[0]*time_old));
  deltaX[1] = -Ampl[1]*(sin(Omega[1]*time_new) - sin(Omega[1]*time_old));
  deltaX[2] = -Ampl[2]*(sin(Omega[2]*time_new) - sin(Omega[2]*time_old));
  
  /*--- Compute grid velocity due to plunge in the x, y, & z directions. ---*/
  xDot[0] = -Ampl[0]*Omega[0]*(cos(Omega[0]*time_new));
  xDot[1] = -Ampl[1]*Omega[1]*(cos(Omega[1]*time_new));
  xDot[2] = -Ampl[2]*Omega[2]*(cos(Omega[2]*time_new));
  
  if (rank == MASTER_NODE && iter == 0) {
    cout << " Plunging frequency: (" << Omega[0] << ", " << Omega[1];
    cout << ", " << Omega[2] << ") rad/s." << endl;
    cout << " Plunging amplitude: (" << Ampl[0] << ", ";
    cout << Ampl[1] << ", " << Ampl[2] <<  ") m."<< endl;
  }
  
  /*--- Loop over and move each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord   = geometry->node[iPoint]->GetCoord();
    GridVel = geometry->node[iPoint]->GetGridVel();
    
    /*--- Increment the node position using the delta values. ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      newCoord[iDim] = Coord[iDim] + deltaX[iDim];
    
    /*--- Cross Product of angular velocity and distance from center.
     Note that we have assumed the grid velocities have been set to
     an initial value in the plunging routine. ---*/
    
    newGridVel[0] = GridVel[0] + xDot[0];
    newGridVel[1] = GridVel[1] + xDot[1];
   if (nDim == 3) newGridVel[2] = GridVel[2] + xDot[2];
    
    /*--- Store new node location & grid velocity. Do not store the grid
     velocity if this is an adjoint calculation. ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim, newCoord[iDim]);
      if (!adjoint) geometry->node[iPoint]->SetGridVel(iDim, newGridVel[iDim]);
    }
  }
  
  /*--- Set the mesh motion center to the new location after
   incrementing the position with the rigid translation. This
   new location will be used for subsequent pitching/rotation.---*/
  
  for (iDim = 0; iDim < 3; iDim++){
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
  
  for (unsigned short jMarker=0; jMarker<config->GetnMarker_Monitoring(); jMarker++) {
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

void CVolumetricMovement::Rigid_Translation(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter) {
  
  /*--- Local variables ---*/
  su2double deltaX[3], newCoord[3], Center[3], *Coord;
  su2double xDot[3];
  su2double deltaT, time_new, time_old;
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint;
  bool harmonic_balance = (config->GetTime_Marching() == HARMONIC_BALANCE);
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());
  
  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND();
  
  /*--- Get motion center and translation rates from config ---*/
  
  for (iDim = 0; iDim < 3; iDim++){
    Center[iDim] = config->GetMotion_Origin(iDim);
    xDot[iDim]   = config->GetTranslation_Rate(iDim);
  }
  
  if (harmonic_balance) {
    /*--- period of oscillation & time interval using nTimeInstances ---*/
    su2double period = config->GetHarmonicBalance_Period();
    period /= config->GetTime_Ref();
    deltaT = period/(su2double)(config->GetnTimeInstances());
  }
  
  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/
    unsigned long nFlowIter  = config->GetnTime_Iter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<su2double>(directIter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<su2double>(directIter)+1.0)*deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<su2double>(iter)*deltaT;
    if (harmonic_balance) {
      /*--- For harmonic balance, begin movement from the zero position ---*/
      time_old = 0.0;
    } else {
      time_old = time_new;
      if (iter != 0) time_old = (static_cast<su2double>(iter)-1.0)*deltaT;
    }
  }
  
  /*--- Compute delta change in the position in the x, y, & z directions. ---*/
  deltaX[0] = xDot[0]*(time_new-time_old);
  deltaX[1] = xDot[1]*(time_new-time_old);
  deltaX[2] = xDot[2]*(time_new-time_old);

  if (rank == MASTER_NODE) {
    cout << " New physical time: " << time_new << " seconds." << endl;
    if (iter == 0) {
    cout << " Translational velocity: (" << xDot[0]*config->GetVelocity_Ref() << ", " << xDot[1]*config->GetVelocity_Ref();
      cout << ", " << xDot[2]*config->GetVelocity_Ref();
      if (config->GetSystemMeasurements() == SI) cout << ") m/s." << endl;
      else cout << ") ft/s." << endl;
    }
  }
  
  /*--- Loop over and move each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord = geometry->node[iPoint]->GetCoord();
    
    /*--- Increment the node position using the delta values. ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      newCoord[iDim] = Coord[iDim] + deltaX[iDim];
    
    /*--- Store new node location & grid velocity. Do not store the grid
     velocity if this is an adjoint calculation. ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim, newCoord[iDim]);
      if (!adjoint) geometry->node[iPoint]->SetGridVel(iDim,xDot[iDim]);
    }
  }
  
  /*--- Set the mesh motion center to the new location after
   incrementing the position with the rigid translation. This
   new location will be used for subsequent pitching/rotation.---*/
  
  for (iDim = 0; iDim < 3; iDim++){
    Center[iDim] = config->GetMotion_Origin(iDim) + deltaX[iDim];
  } 
  config->SetMotion_Origin(Center);

  
  /*--- Set the moment computation center to the new location after
   incrementing the position with the translation. ---*/
  
  for (unsigned short jMarker=0; jMarker<config->GetnMarker_Monitoring(); jMarker++) {
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

void CVolumetricMovement::SetVolume_Scaling(CGeometry *geometry, CConfig *config, bool UpdateGeo) {

  unsigned short iDim;
  unsigned long iPoint;
  su2double newCoord[3] = {0.0,0.0,0.0}, *Coord;

  /*--- The scaling factor is the only input to this option. Currently, 
   the mesh must be scaled the same amount in all three directions. ---*/
  
  su2double Scale = config->GetDV_Value(0)*config->GetOpt_RelaxFactor();
  
  if (rank == MASTER_NODE) {
    cout << "Scaling the mesh by a constant factor of " << Scale << "." << endl;
  }
  
  /*--- Loop over and move each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord = geometry->node[iPoint]->GetCoord();
    
    /*--- Scale the node position by the specified factor. ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      newCoord[iDim] = Scale*Coord[iDim];
    
    /*--- Store the new node location. ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim, newCoord[iDim]);
    }
  }

  /*--- After moving all nodes, update geometry class ---*/
  if (UpdateGeo) UpdateDualGrid(geometry, config);
  
}

void CVolumetricMovement::SetVolume_Translation(CGeometry *geometry, CConfig *config, bool UpdateGeo)  {

  unsigned short iDim;
  unsigned long iPoint;
  su2double *Coord, deltaX[3] = {0.0,0.0,0.0}, newCoord[3] = {0.0,0.0,0.0};
  su2double Scale = config->GetOpt_RelaxFactor();
  
  /*--- Get the unit vector and magnitude of displacement. Note that we
   assume this is the first DV entry since it is for mesh translation.
   Create the displacement vector from the magnitude and direction. ---*/
  
  su2double Ampl = config->GetDV_Value(0)*Scale;
  su2double length = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    deltaX[iDim] = config->GetParamDV(0, iDim);
    length += deltaX[iDim]*deltaX[iDim];
  }
  length = sqrt(length);
  for (iDim = 0; iDim < nDim; iDim++)
    deltaX[iDim] = Ampl*deltaX[iDim]/length;
  if (rank == MASTER_NODE) {
    cout << "Translational displacement: (" << deltaX[0] << ", ";
    cout  << deltaX[1] << ", " << deltaX[2] << ")." << endl;
  }
  
  /*--- Loop over and move each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord = geometry->node[iPoint]->GetCoord();
    
    /*--- Increment the node position using the delta values. ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      newCoord[iDim] = Coord[iDim] + deltaX[iDim];
    
    /*--- Store new node location. ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim, newCoord[iDim]);
    }
  }
  
  /*--- After moving all nodes, update geometry class ---*/
  if (UpdateGeo) UpdateDualGrid(geometry, config);
  
}

void CVolumetricMovement::SetVolume_Rotation(CGeometry *geometry, CConfig *config, bool UpdateGeo) {
  
  unsigned short iDim;
  unsigned long iPoint;
  su2double x, y, z;
  su2double *Coord, deltaX[3] = {0.0,0.0,0.0}, newCoord[3] = {0.0,0.0,0.0};
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- xyz-coordinates of a point on the line of rotation. */
  su2double a = config->GetParamDV(0, 0);
  su2double b = config->GetParamDV(0, 1);
  su2double c = 0.0;
  if (geometry->GetnDim() == 3) c = config->GetParamDV(0,2);
  
  /*--- xyz-coordinate of the line's direction vector. ---*/
  su2double u = config->GetParamDV(0, 3)-config->GetParamDV(0, 0);
  su2double v = config->GetParamDV(0, 4)-config->GetParamDV(0, 1);
  su2double w = 1.0;
  if (geometry->GetnDim() == 3)
    w = config->GetParamDV(0, 5)-config->GetParamDV(0, 2);
  
  /*--- The angle of rotation. ---*/
  su2double theta = config->GetDV_Value(0)*Scale*PI_NUMBER/180.0;
  
  /*--- Print to the console. ---*/
  if (rank == MASTER_NODE) {
    cout << "Rotation axis vector: (" << u << ", ";
    cout << v << ", " << w << ")." << endl;
    cout << "Angle of rotation: " << config->GetDV_Value(0)*Scale;
    cout << " degrees." << endl;
  }
  
  /*--- Intermediate values used in computations. ---*/
  su2double u2=u*u; su2double v2=v*v; su2double w2=w*w;
  su2double cosT = cos(theta); su2double sinT = sin(theta);
  su2double l2 = u2 + v2 + w2; su2double l = sqrt(l2);
  
  /*--- Loop over and move each node in the volume mesh ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord = geometry->node[iPoint]->GetCoord();
    
    /*--- Displacement for this point due to the rotation. ---*/
    x = Coord[0]; y = Coord[1]; z = 0.0;
    if (geometry->GetnDim() == 3) z = Coord[2];
    
    deltaX[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
    + (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
    + l*(-c*v + b*w - w*y + v*z)*sinT;
    deltaX[0] = deltaX[0]/l2 - x;
    
    deltaX[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z)
    + (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
    + l*(c*u - a*w + w*x - u*z)*sinT;
    deltaX[1] = deltaX[1]/l2 - y;
    
    deltaX[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z)
    + (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
    + l*(-b*u + a*v - v*x + u*y)*sinT;
    if (geometry->GetnDim() == 3) deltaX[2] = deltaX[2]/l2 - z;
    else deltaX[2] = 0.0;
    
    /*--- Increment the node position using the delta values. ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      newCoord[iDim] = Coord[iDim] + deltaX[iDim];
    
    /*--- Store new node location. ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim, newCoord[iDim]);
    }
  }
 
  /*--- After moving all nodes, update geometry class ---*/
  if (UpdateGeo) UpdateDualGrid(geometry, config);
  
}

CSurfaceMovement::CSurfaceMovement(void) : CGridMovement() {
  
  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();
  
  nFFDBox = 0;
  nLevel = 0;
  FFDBoxDefinition = false;
}

CSurfaceMovement::~CSurfaceMovement(void) {}

void CSurfaceMovement::SetSurface_Deformation(CGeometry *geometry, CConfig *config) {
  
  unsigned short iFFDBox, iDV, iLevel, iChild, iParent, jFFDBox, iMarker;
  unsigned short Degree_Unitary [] = {1,1,1}, BSpline_Unitary [] = {2,2,2};
  su2double MaxDiff, Current_Scale, Ratio, New_Scale;
  string FFDBoxTag;
   bool allmoving;
  
  bool cylindrical = (config->GetFFD_CoordSystem() == CYLINDRICAL);
  bool spherical   = (config->GetFFD_CoordSystem() == SPHERICAL);
  bool polar       = (config->GetFFD_CoordSystem() == POLAR);
  bool cartesian   = (config->GetFFD_CoordSystem() == CARTESIAN);
  su2double BoundLimit = config->GetOpt_LineSearch_Bound();

  /*--- Setting the Free Form Deformation ---*/
  
  if (config->GetDesign_Variable(0) == FFD_SETTING) {
    
    /*--- Definition of the FFD deformation class ---*/
    
    FFDBox = new CFreeFormDefBox*[MAX_NUMBER_FFD];
    
    /*--- Read the FFD information from the config file ---*/
    
    ReadFFDInfo(geometry, config, FFDBox);
    
    /*--- If there is a FFDBox in the input file ---*/
    
    if (nFFDBox != 0) {
      
      /*--- if polar coordinates, trnasform the corner to polar ---*/
      
      if (cylindrical) {
        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
          FFDBox[iFFDBox]->SetCart2Cyl_CornerPoints(config);
        }
      }
      else if (spherical || polar) {
        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
          FFDBox[iFFDBox]->SetCart2Sphe_CornerPoints(config);
        }
      }
      
      /*--- If the FFDBox was not defined in the input file ---*/
      
      if ((rank == MASTER_NODE) && (GetnFFDBox() != 0)) {
        if (cartesian) cout << endl <<"----------------- FFD technique (cartesian -> parametric) ---------------" << endl;
        else if (cylindrical) cout << endl <<"----------------- FFD technique (cylinder -> parametric) ---------------" << endl;
        else if (spherical) cout << endl <<"----------------- FFD technique (spherical -> parametric) ---------------" << endl;
        else if (polar) cout << endl <<"----------------- FFD technique (polar -> parametric) ---------------" << endl;
      }
      
      /*--- Create a unitary FFDBox as baseline for other FFDBoxes shapes ---*/
      
      CFreeFormDefBox FFDBox_unitary(Degree_Unitary, BSpline_Unitary, BEZIER);
      FFDBox_unitary.SetUnitCornerPoints();

      /*--- Compute the control points of the unitary box, in this case the degree is 1 and the order is 2 ---*/
      
      FFDBox_unitary.SetControlPoints_Parallelepiped();
      
      for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
        
        /*--- Compute the support control points for the final FFD using the unitary box ---*/
        
        FFDBox_unitary.SetSupportCP(FFDBox[iFFDBox]);
        
        /*--- Compute control points in the support box ---*/
        
        FFDBox_unitary.SetSupportCPChange(FFDBox[iFFDBox]);
        
        /*--- Compute the parametric coordinates, it also find the points in
         the FFDBox using the parametrics coordinates ---*/
        
        SetParametricCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);
        
        
        /*--- If polar coordinates, transform the corners and control points to cartesians ---*/
        
        if (cylindrical) {
          FFDBox[iFFDBox]->SetCyl2Cart_CornerPoints(config);
          FFDBox[iFFDBox]->SetCyl2Cart_ControlPoints(config);
        }
        else if (spherical || polar) {
          FFDBox[iFFDBox]->SetSphe2Cart_CornerPoints(config);
          FFDBox[iFFDBox]->SetSphe2Cart_ControlPoints(config);
        }
      }
      /*--- Output original FFD FFDBox ---*/
      
      if (rank == MASTER_NODE) {
        for (unsigned short iFile = 0; iFile < config->GetnVolumeOutputFiles(); iFile++){
          unsigned short *FileFormat = config->GetVolumeOutputFiles();
          if (FileFormat[iFile] == PARAVIEW || FileFormat[iFile] == PARAVIEW_BINARY) {
            cout << "Writing a Paraview file of the FFD boxes." << endl;
            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetParaview(geometry, iFFDBox, true);
            }
          } else if (FileFormat[iFile] == TECPLOT || FileFormat[iFile] == TECPLOT_BINARY) {
            cout << "Writing a Tecplot file of the FFD boxes." << endl;
            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetTecplot(geometry, iFFDBox, true);
            }
          }
          else if (FileFormat[iFile] == CGNS)  {
            cout << "Writing a CGNS file of the FFD boxes." << endl;
            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetCGNS(geometry, iFFDBox, true);
            }
          }
        }  
      }
    }
    
    else {
      SU2_MPI::Error("There are not FFD boxes in the mesh file!!", CURRENT_FUNCTION);
    }
    
  }
  
  /*--- Free Form deformation based ---*/
  
  if ((config->GetDesign_Variable(0) == FFD_CONTROL_POINT_2D) ||
      (config->GetDesign_Variable(0) == FFD_CAMBER_2D) ||
      (config->GetDesign_Variable(0) == FFD_THICKNESS_2D) ||
      (config->GetDesign_Variable(0) == FFD_TWIST_2D) ||
      (config->GetDesign_Variable(0) == FFD_CONTROL_POINT) ||
      (config->GetDesign_Variable(0) == FFD_NACELLE) ||
      (config->GetDesign_Variable(0) == FFD_GULL) ||
      (config->GetDesign_Variable(0) == FFD_TWIST) ||
      (config->GetDesign_Variable(0) == FFD_ROTATION) ||
      (config->GetDesign_Variable(0) == FFD_CONTROL_SURFACE) ||
      (config->GetDesign_Variable(0) == FFD_CAMBER) ||
      (config->GetDesign_Variable(0) == FFD_THICKNESS) ||
      (config->GetDesign_Variable(0) == FFD_ANGLE_OF_ATTACK)) {
    
    /*--- Definition of the FFD deformation class ---*/
    
    FFDBox = new CFreeFormDefBox*[MAX_NUMBER_FFD];
    
    /*--- Read the FFD information from the grid file ---*/
    
    ReadFFDInfo(geometry, config, FFDBox, config->GetMesh_FileName());
    
    /*--- If there is a FFDBox in the input file ---*/
    
    if (nFFDBox != 0) {
      
      /*--- If the FFDBox was not defined in the input file ---*/
      
      if (!GetFFDBoxDefinition()) {
        SU2_MPI::Error(string("There is not FFD box definition in the mesh file,\n") +
                       string("run DV_KIND=FFD_SETTING first !!"), CURRENT_FUNCTION);
      }
      
      /* --- Check if the FFD boxes referenced in the design variable definition can be found --- */
      
      for (iDV = 0; iDV < config->GetnDV(); iDV++) {
        if (!CheckFFDBoxDefinition(config, iDV)) {
          SU2_MPI::Error(string("There is no FFD box with tag \"") + config->GetFFDTag(iDV) + string("\" defined in the mesh file.\n") +
                         string("Check the definition of the design variables and/or the FFD settings !!"), CURRENT_FUNCTION);
        }
      }
      
      /*--- Check that the user has specified a non-zero number of surfaces to move with DV_MARKER. ---*/
      
      if (config->GetnMarker_DV() == 0) {
        SU2_MPI::Error(string("No markers are specified in DV_MARKER, so no deformation will occur.\n") +
                       string("List markers to be deformed in DV_MARKER."), CURRENT_FUNCTION);
      }
    
      /*--- Output original FFD FFDBox ---*/
      
      if ((rank == MASTER_NODE) && (config->GetKind_SU2() != SU2_DOT)) {
        
        for (unsigned short iFile = 0; iFile < config->GetnVolumeOutputFiles(); iFile++){
          unsigned short *FileFormat = config->GetVolumeOutputFiles();
          
          if (FileFormat[iFile] == PARAVIEW || FileFormat[iFile] == PARAVIEW_BINARY) {
            cout << "Writing a Paraview file of the FFD boxes." << endl;
            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetParaview(geometry, iFFDBox, true);
            }
          } else if (FileFormat[iFile] == TECPLOT || FileFormat[iFile] == TECPLOT_BINARY) {
            cout << "Writing a Tecplot file of the FFD boxes." << endl;
            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetTecplot(geometry, iFFDBox, true);
            }
          }
          else if (FileFormat[iFile] == CGNS)  {
            cout << "Writing a CGNS file of the FFD boxes." << endl;
            for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
              FFDBox[iFFDBox]->SetCGNS(geometry, iFFDBox, true);
            }
          }
        }
      }
      
      /*--- If polar FFD, change the coordinates system ---*/
      
      if (cylindrical) {
        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
          FFDBox[iFFDBox]->SetCart2Cyl_CornerPoints(config);
          FFDBox[iFFDBox]->SetCart2Cyl_ControlPoints(config);
        }
      }
      else if (spherical || polar) {
        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
          FFDBox[iFFDBox]->SetCart2Sphe_CornerPoints(config);
          FFDBox[iFFDBox]->SetCart2Sphe_ControlPoints(config);
        }
      }
      
      /*--- Apply the deformation to the orifinal FFD box ---*/
      
      if ((rank == MASTER_NODE) && (GetnFFDBox() != 0))
        cout << endl <<"----------------- FFD technique (parametric -> cartesian) ---------------" << endl;
      
      /*--- Loop over all the FFD boxes levels ---*/
      
      for (iLevel = 0; iLevel < GetnLevel(); iLevel++) {
        
        /*--- Loop over all FFD FFDBoxes ---*/
        
        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
          
          /*--- Check the level of the FFD box ---*/
          
          if (FFDBox[iFFDBox]->GetLevel() == iLevel) {
            
            /*--- Check the dimension of the FFD compared with the design variables ---*/
            
            if (rank == MASTER_NODE) cout << "Checking FFD box dimension." << endl;
            CheckFFDDimension(geometry, config, FFDBox[iFFDBox], iFFDBox);
            
            /*--- Compute intersections of the FFD box with the surface to eliminate design
             variables and satisfy surface continuity ---*/
            
            if (rank == MASTER_NODE) cout << "Checking FFD box intersections with the solid surfaces." << endl;
            CheckFFDIntersections(geometry, config, FFDBox[iFFDBox], iFFDBox);
            
            /*--- Compute the parametric coordinates of the child box
             control points (using the parent FFDBox)  ---*/
            
            for (iChild = 0; iChild < FFDBox[iFFDBox]->GetnChildFFDBox(); iChild++) {
              FFDBoxTag = FFDBox[iFFDBox]->GetChildFFDBoxTag(iChild);
              for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
                if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
              SetParametricCoordCP(geometry, config, FFDBox[iFFDBox], FFDBox[jFFDBox]);
            }
            
            /*--- Update the parametric coordinates if it is a child FFDBox ---*/
            
            if (iLevel > 0) UpdateParametricCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);
            
            /*--- Apply the design variables to the control point position ---*/
            
            for (iDV = 0; iDV < config->GetnDV(); iDV++) {
              switch ( config->GetDesign_Variable(iDV) ) {
                case FFD_CONTROL_POINT_2D : SetFFDCPChange_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                case FFD_CAMBER_2D :        SetFFDCamber_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                case FFD_THICKNESS_2D :     SetFFDThickness_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                case FFD_TWIST_2D :         SetFFDTwist_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                case FFD_CONTROL_POINT :    SetFFDCPChange(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                case FFD_NACELLE :          SetFFDNacelle(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                case FFD_GULL :             SetFFDGull(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                case FFD_TWIST :            SetFFDTwist(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                case FFD_ROTATION :         SetFFDRotation(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                case FFD_CONTROL_SURFACE :  SetFFDControl_Surface(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                case FFD_CAMBER :           SetFFDCamber(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                case FFD_THICKNESS :        SetFFDThickness(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                case FFD_ANGLE_OF_ATTACK :  SetFFDAngleOfAttack(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
              }
            }
            
            /*--- Recompute cartesian coordinates using the new control point location ---*/
            
            MaxDiff = SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, false);
            
            if ((MaxDiff > BoundLimit) && (config->GetKind_SU2() == SU2_DEF)) {
              
              if (rank == MASTER_NODE) cout << "Out-of-bounds, re-adjusting scale factor to safisfy line search limit." << endl;
              
              Current_Scale = config->GetOpt_RelaxFactor();
              Ratio = (BoundLimit/MaxDiff);
              New_Scale = Current_Scale *(Ratio-1.0);
              config->SetOpt_RelaxFactor(New_Scale);
              
              /*--- Apply the design variables to the control point position ---*/
              
              for (iDV = 0; iDV < config->GetnDV(); iDV++) {
                switch ( config->GetDesign_Variable(iDV) ) {
                  case FFD_CONTROL_POINT_2D : SetFFDCPChange_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                  case FFD_CAMBER_2D :        SetFFDCamber_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                  case FFD_THICKNESS_2D :     SetFFDThickness_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                  case FFD_TWIST_2D :         SetFFDTwist_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                  case FFD_CONTROL_POINT :    SetFFDCPChange(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                  case FFD_NACELLE :          SetFFDNacelle(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                  case FFD_GULL :             SetFFDGull(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                  case FFD_TWIST :            SetFFDTwist(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                  case FFD_ROTATION :         SetFFDRotation(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                  case FFD_CONTROL_SURFACE :  SetFFDControl_Surface(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                  case FFD_CAMBER :           SetFFDCamber(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                  case FFD_THICKNESS :        SetFFDThickness(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                  case FFD_ANGLE_OF_ATTACK :  SetFFDAngleOfAttack(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, false); break;
                }
              }
              
              /*--- Recompute cartesian coordinates using the new control point location ---*/
              
              MaxDiff = SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, false);
              
            }
            
            /*--- Reparametrization of the parent FFD box ---*/
            
            for (iParent = 0; iParent < FFDBox[iFFDBox]->GetnParentFFDBox(); iParent++) {
              FFDBoxTag = FFDBox[iFFDBox]->GetParentFFDBoxTag(iParent);
              for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
                if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
              UpdateParametricCoord(geometry, config, FFDBox[jFFDBox], jFFDBox);
            }
            
            /*--- Compute the new location of the control points of the child boxes
             (using the parent FFDBox) ---*/
            
            for (iChild = 0; iChild < FFDBox[iFFDBox]->GetnChildFFDBox(); iChild++) {
              FFDBoxTag = FFDBox[iFFDBox]->GetChildFFDBoxTag(iChild);
              for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
                if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
              GetCartesianCoordCP(geometry, config, FFDBox[iFFDBox], FFDBox[jFFDBox]);
            }
          }
        }
        
        /*--- If polar, compute the cartesians coordinates ---*/
        
        if (cylindrical) {
          for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
            FFDBox[iFFDBox]->SetCyl2Cart_CornerPoints(config);
            FFDBox[iFFDBox]->SetCyl2Cart_ControlPoints(config);
          }
        }
        else if (spherical || polar) {
          for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
            FFDBox[iFFDBox]->SetSphe2Cart_CornerPoints(config);
            FFDBox[iFFDBox]->SetSphe2Cart_ControlPoints(config);
          }
        }
        
        /*--- Output the deformed FFD Boxes ---*/
        
        if ((rank == MASTER_NODE) && (config->GetKind_SU2() != SU2_DOT)) {
          
          for (unsigned short iFile = 0; iFile < config->GetnVolumeOutputFiles(); iFile++){
            unsigned short *FileFormat = config->GetVolumeOutputFiles();
            
            if (FileFormat[iFile] == PARAVIEW || FileFormat[iFile] == PARAVIEW_BINARY) {
              cout << "Writing a Paraview file of the FFD boxes." << endl;
              for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
                FFDBox[iFFDBox]->SetParaview(geometry, iFFDBox, false);
              }
            } else if (FileFormat[iFile] == TECPLOT || FileFormat[iFile] == TECPLOT_BINARY) {
              cout << "Writing a Tecplot file of the FFD boxes." << endl;
              for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
                FFDBox[iFFDBox]->SetTecplot(geometry, iFFDBox, false);
              }
            }
            else if (FileFormat[iFile] == CGNS)  {
              cout << "Writing a CGNS file of the FFD boxes." << endl;
              for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
                FFDBox[iFFDBox]->SetCGNS(geometry, iFFDBox, false);
              }
            }
          }
        }
      }
    }
    
    else {
      SU2_MPI::Error("There are no FFD Boxes in the mesh file!!", CURRENT_FUNCTION);
    }
    
  }
  
  /*--- External surface file based ---*/
  
  else if (config->GetDesign_Variable(0) == SURFACE_FILE) {
    
    /*--- Check whether a surface file exists for input ---*/
    ofstream Surface_File;
    string filename = config->GetDV_Filename();
    Surface_File.open(filename.c_str(), ios::in);
    
    /*--- A surface file does not exist, so write a new one for the
     markers that are specified as part of the motion. ---*/
    if (Surface_File.fail()) {
      
      if (rank == MASTER_NODE && size == SINGLE_NODE) {
        cout << "No surface positions file found. Writing a template file: " << filename << "." << endl;
        
        Surface_File.open(filename.c_str(), ios::out);
        Surface_File.precision(15);
        unsigned long iMarker, jPoint, GlobalIndex, iVertex; su2double *Coords;
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
          if (config->GetMarker_All_DV(iMarker) == YES) {
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
              jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              GlobalIndex = geometry->node[jPoint]->GetGlobalIndex();
              Coords = geometry->node[jPoint]->GetCoord();
              Surface_File << GlobalIndex << "\t" << Coords[0] << "\t" << Coords[1];
              if (geometry->GetnDim() == 2) Surface_File << endl;
              else Surface_File << "\t" << Coords[2] << endl;
            }
          }
        }
        Surface_File.close();
        
      } else {
        SU2_MPI::Error("No surface positions file found and template writing not yet supported in parallel.\n To generate a template surface positions file, run SU2_DEF again in serial.", CURRENT_FUNCTION);
      }
    }
    
    else {
      /*--- A surface file exists, so read in the coordinates ---*/
      Surface_File.close();
      if (rank == MASTER_NODE) cout << "Updating the surface coordinates from the input file." << endl;
      SetExternal_Deformation(geometry, config, ZONE_0, 0);
    }
    
  }
    
  else if ((config->GetDesign_Variable(0) == ROTATION) ||
           (config->GetDesign_Variable(0) == TRANSLATION) ||
           (config->GetDesign_Variable(0) == SCALE) ||
           (config->GetDesign_Variable(0) == HICKS_HENNE) ||
           (config->GetDesign_Variable(0) == SURFACE_BUMP) ||
           (config->GetDesign_Variable(0) == ANGLE_OF_ATTACK)) {
    
    /*--- Apply rotation, displacement and stretching design variables (this
     should be done before the bump function design variables) ---*/
    
    for (iDV = 0; iDV < config->GetnDV(); iDV++) {
      switch ( config->GetDesign_Variable(iDV) ) {
        case SCALE :  SetScale(geometry, config, iDV, false); break;
        case TRANSLATION :  SetTranslation(geometry, config, iDV, false); break;
        case ROTATION :     SetRotation(geometry, config, iDV, false); break;
      }
    }
    
    /*--- Apply the design variables to the control point position ---*/
    
    for (iDV = 0; iDV < config->GetnDV(); iDV++) {
      switch ( config->GetDesign_Variable(iDV) ) {
        case HICKS_HENNE :  SetHicksHenne(geometry, config, iDV, false); break;
      }
    }
    
    /*--- Apply the design variables to the control point position ---*/

    for (iDV = 0; iDV < config->GetnDV(); iDV++) {
      switch ( config->GetDesign_Variable(iDV) ) {
        case SURFACE_BUMP :  SetSurface_Bump(geometry, config, iDV, false); break;
      }
    }

    /*--- Apply the angle of attack design variable ---*/
    
    for (iDV = 0; iDV < config->GetnDV(); iDV++) {
      switch ( config->GetDesign_Variable(iDV) ) {
        case ANGLE_OF_ATTACK :  SetAngleOfAttack(geometry, config, iDV, false); break;
      }
    }
    
  }
  
  /*--- NACA_4Digits design variable ---*/
  
  else if (config->GetDesign_Variable(0) == NACA_4DIGITS) { SetNACA_4Digits(geometry, config); }
  
  /*--- Parabolic airfoil design variable ---*/
  
  else if (config->GetDesign_Variable(0) == PARABOLIC) { SetParabolic(geometry, config); }
  
  /*--- Airfoil from file design variable ---*/
  
  else if (config->GetDesign_Variable(0) == AIRFOIL) { SetAirfoil(geometry, config); }
  
  /*--- FFD setting ---*/
  
  else if (config->GetDesign_Variable(0) == FFD_SETTING) {
    if (rank == MASTER_NODE)
      cout << "No surface deformation (setting FFD)." << endl;
  }
  
  /*--- Scale, Translate, and Rotate will be done with rigid mesh transforms. ---*/
  
  else if ((config->GetDesign_Variable(0) == ROTATION) ||
           (config->GetDesign_Variable(0) == TRANSLATION) ||
           (config->GetDesign_Variable(0) == SCALE)) {
    
    /*--- If all markers are deforming, use volume method.
     If only some are deforming, use surface method ---*/
    
    /*--- iDV was uninitialized, so hard-coding to one. Check intended
     behavior (might want to loop over all iDV in case we have trans & rotate. ---*/
    iDV = 0;
    allmoving = true;
    
    /*--- Loop over markers ---*/
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_DV(iMarker) == NO)
        allmoving = false;
    }
    
    if (!allmoving) {
      /*---Only some markers are moving, use the surface method ---*/
      if (config->GetDesign_Variable(0) == ROTATION)
        SetRotation(geometry, config, iDV, false);
      if (config->GetDesign_Variable(0) == SCALE)
        SetScale(geometry, config, iDV, false);
      if (config->GetDesign_Variable(0) == TRANSLATION)
        SetTranslation(geometry, config, iDV, false);
    }
    else {
      if (rank == MASTER_NODE)
        cout << "No surface deformation (scaling, rotation, or translation)." << endl;
    }
  }
  
  /*--- Design variable not implement ---*/
  
  else {
    if (rank == MASTER_NODE)
      cout << "Design Variable not implemented yet" << endl;
  }
  
}


void CSurfaceMovement::SetSurface_Derivative(CGeometry *geometry, CConfig *config) {

  su2double DV_Value = 0.0;

  unsigned short iDV = 0, iDV_Value = 0;

  for (iDV = 0; iDV < config->GetnDV(); iDV++) {
    for (iDV_Value = 0; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++) {

      DV_Value = config->GetDV_Value(iDV, iDV_Value);

      /*--- If value of the design variable is not 0.0 we apply the differentation.
     *     Note if multiple variables are non-zero, we end up with the sum of all the derivatives. ---*/

      if (DV_Value != 0.0) {

        DV_Value = 0.0;

        SU2_TYPE::SetDerivative(DV_Value, 1.0);

        config->SetDV_Value(iDV, iDV_Value, DV_Value);
      }
    }
  }

  /*--- Run the surface deformation with DV_Value = 0.0 (no deformation at all) ---*/

  SetSurface_Deformation(geometry, config);
}

void CSurfaceMovement::CopyBoundary(CGeometry *geometry, CConfig *config) {
  
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  su2double *Coord;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Coord = geometry->node[iPoint]->GetCoord();
      geometry->vertex[iMarker][iVertex]->SetCoord(Coord);
    }
  }
  
}

void CSurfaceMovement::SetParametricCoord(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox) {
  
  unsigned short iMarker, iDim, iOrder, jOrder, kOrder, lOrder, mOrder, nOrder;
  unsigned long iVertex, iPoint, TotalVertex = 0;
  su2double *CartCoordNew, *ParamCoord, CartCoord[3], ParamCoordGuess[3], MaxDiff, my_MaxDiff = 0.0, Diff, *Coord;
  unsigned short nDim = geometry->GetnDim();
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;

  unsigned short BoxFFD = true;
  bool cylindrical = (config->GetFFD_CoordSystem() == CYLINDRICAL);
  bool spherical = (config->GetFFD_CoordSystem() == SPHERICAL);
  bool polar = (config->GetFFD_CoordSystem() == POLAR);
  
  /*--- Change order and control points reduce the
   complexity of the point inversion (this only works with boxes,
 in case of Bezier curves, and we maintain an internal copy)---*/
  
  if (BoxFFD && (config->GetFFD_Blending() == BEZIER)) {
    
    for (iOrder = 0; iOrder < 2; iOrder++) {
      for (jOrder = 0; jOrder < 2; jOrder++) {
        for (kOrder = 0; kOrder < 2; kOrder++) {
          
          lOrder = 0; mOrder = 0; nOrder = 0;
          if (iOrder == 1) {lOrder = FFDBox->GetlOrder()-1;}
          if (jOrder == 1) {mOrder = FFDBox->GetmOrder()-1;}
          if (kOrder == 1) {nOrder = FFDBox->GetnOrder()-1;}
          
          Coord = FFDBox->GetCoordControlPoints(lOrder, mOrder, nOrder);
          
          FFDBox->SetCoordControlPoints(Coord, iOrder, jOrder, kOrder);
          
        }
      }
    }
    
    FFDBox->SetlOrder(2); FFDBox->SetmOrder(2); FFDBox->SetnOrder(2);
    FFDBox->SetnControlPoints();
    FFDBox->BlendingFunction[0]->SetOrder(2, 2);
    FFDBox->BlendingFunction[1]->SetOrder(2, 2);
    FFDBox->BlendingFunction[2]->SetOrder(2, 2);
  }
  /*--- Point inversion algorithm with a basic box ---*/
  
  ParamCoordGuess[0]  = 0.5; ParamCoordGuess[1] = 0.5; ParamCoordGuess[2] = 0.5;
  CartCoord[0]        = 0.0; CartCoord[1]       = 0.0; CartCoord[2]       = 0.0;
  
  /*--- Count the number of vertices ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_DV(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
        TotalVertex++;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_DV(iMarker) == YES) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        /*--- Get the cartesian coordinates ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          CartCoord[iDim] = geometry->vertex[iMarker][iVertex]->GetCoord(iDim);
        
        /*--- Transform the cartesian into polar ---*/
        
        if (cylindrical) {
          X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
          
          Xbar =  CartCoord[0] - X_0; Ybar =  CartCoord[1] - Y_0; Zbar =  CartCoord[2] - Z_0;
          
          CartCoord[0] = sqrt(Ybar*Ybar + Zbar*Zbar);
          CartCoord[1] = atan2(Zbar, Ybar); if (CartCoord[1] > PI_NUMBER/2.0) CartCoord[1] -= 2.0*PI_NUMBER;
          CartCoord[2] = Xbar;
        }
        else if (spherical || polar) {
          X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
          
          Xbar =  CartCoord[0] - X_0; Ybar =  CartCoord[1] - Y_0; Zbar =  CartCoord[2] - Z_0;
          
          CartCoord[0] = sqrt(Xbar*Xbar + Ybar*Ybar + Zbar*Zbar);
          CartCoord[1] = atan2(Zbar, Ybar);  if (CartCoord[1] > PI_NUMBER/2.0) CartCoord[1] -= 2.0*PI_NUMBER;
          CartCoord[2] = acos(Xbar/CartCoord[0]);
        }
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- If the point is inside the FFD, compute the value of the parametric coordinate ---*/
        
        if (FFDBox->GetPointFFD(geometry, config, iPoint)) {
          
          /*--- Find the parametric coordinate ---*/
          
          ParamCoord = FFDBox->GetParametricCoord_Iterative(iPoint, CartCoord, ParamCoordGuess, config);
          
          /*--- Compute the cartesian coordinates using the parametric coordinates
           to check that everything is correct ---*/
          
          CartCoordNew = FFDBox->EvalCartesianCoord(ParamCoord);
          
          /*--- Compute max difference between original value and the recomputed value ---*/
          
          Diff = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Diff += (CartCoordNew[iDim]-CartCoord[iDim])*(CartCoordNew[iDim]-CartCoord[iDim]);
          Diff = sqrt(Diff);
          my_MaxDiff = max(my_MaxDiff, Diff);
          
          /*--- If the parametric coordinates are in (0,1) the point belongs to the FFDBox, using the input tolerance  ---*/
          
          if (((ParamCoord[0] >= - config->GetFFD_Tol()) && (ParamCoord[0] <= 1.0 + config->GetFFD_Tol())) &&
              ((ParamCoord[1] >= - config->GetFFD_Tol()) && (ParamCoord[1] <= 1.0 + config->GetFFD_Tol())) &&
              ((ParamCoord[2] >= - config->GetFFD_Tol()) && (ParamCoord[2] <= 1.0 + config->GetFFD_Tol()))) {
            
            
            /*--- Rectification of the initial tolerance (we have detected situations
             where 0.0 and 1.0 doesn't work properly ---*/
            
            su2double lower_limit = config->GetFFD_Tol();
            su2double upper_limit = 1.0-config->GetFFD_Tol();
            
            if (ParamCoord[0] < lower_limit) ParamCoord[0] = lower_limit;
            if (ParamCoord[1] < lower_limit) ParamCoord[1] = lower_limit;
            if (ParamCoord[2] < lower_limit) ParamCoord[2] = lower_limit;
            if (ParamCoord[0] > upper_limit) ParamCoord[0] = upper_limit;
            if (ParamCoord[1] > upper_limit) ParamCoord[1] = upper_limit;
            if (ParamCoord[2] > upper_limit) ParamCoord[2] = upper_limit;
            
            /*--- Set the value of the parametric coordinate ---*/
            
            FFDBox->Set_MarkerIndex(iMarker);
            FFDBox->Set_VertexIndex(iVertex);
            FFDBox->Set_PointIndex(iPoint);
            FFDBox->Set_ParametricCoord(ParamCoord);
            FFDBox->Set_CartesianCoord(CartCoord);
            
            ParamCoordGuess[0] = ParamCoord[0]; ParamCoordGuess[1] = ParamCoord[1]; ParamCoordGuess[2] = ParamCoord[2];
            
            if (Diff >= config->GetFFD_Tol()) {
              cout << "Please check this point: Local (" << ParamCoord[0] <<" "<< ParamCoord[1] <<" "<< ParamCoord[2] <<") <-> Global ("
              << CartCoord[0] <<" "<< CartCoord[1] <<" "<< CartCoord[2] <<") <-> Error "<< Diff <<" vs "<< config->GetFFD_Tol() <<"." << endl;
            }
            
          }
          else {
            
            if (Diff >= config->GetFFD_Tol()) {
              cout << "Please check this point: Local (" << ParamCoord[0] <<" "<< ParamCoord[1] <<" "<< ParamCoord[2] <<") <-> Global ("
              << CartCoord[0] <<" "<< CartCoord[1] <<" "<< CartCoord[2] <<") <-> Error "<< Diff <<" vs "<< config->GetFFD_Tol() <<"." << endl;
            }
            
          }
          
        }
      }
    }
  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&my_MaxDiff, &MaxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  MaxDiff = my_MaxDiff;
#endif
  
  if (rank == MASTER_NODE)
    cout << "Compute parametric coord      | FFD box: " << FFDBox->GetTag() << ". Max Diff: " << MaxDiff <<"."<< endl;
  
  
  /*--- After the point inversion, copy the original 
   information back (this only works with boxes,
   and we maintain an internal copy) ---*/
  
  if (BoxFFD) {
    FFDBox->SetOriginalControlPoints();
    if (config->GetFFD_Blending() == BEZIER){
      FFDBox->BlendingFunction[0]->SetOrder(FFDBox->GetlOrder(), FFDBox->GetlOrder());
      FFDBox->BlendingFunction[1]->SetOrder(FFDBox->GetmOrder(), FFDBox->GetmOrder());
      FFDBox->BlendingFunction[2]->SetOrder(FFDBox->GetnOrder(), FFDBox->GetnOrder());
    }
  }
}

void CSurfaceMovement::SetParametricCoordCP(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBoxParent, CFreeFormDefBox *FFDBoxChild) {
  unsigned short iOrder, jOrder, kOrder;
  su2double *CartCoord, *ParamCoord, ParamCoordGuess[3];
  
  for (iOrder = 0; iOrder < FFDBoxChild->GetlOrder(); iOrder++)
    for (jOrder = 0; jOrder < FFDBoxChild->GetmOrder(); jOrder++)
      for (kOrder = 0; kOrder < FFDBoxChild->GetnOrder(); kOrder++) {
        CartCoord = FFDBoxChild->GetCoordControlPoints(iOrder, jOrder, kOrder);
        ParamCoord = FFDBoxParent->GetParametricCoord_Iterative(0, CartCoord, ParamCoordGuess, config);
        FFDBoxChild->SetParCoordControlPoints(ParamCoord, iOrder, jOrder, kOrder);
      }

  if (rank == MASTER_NODE)
    cout << "Compute parametric coord (CP) | FFD parent box: " << FFDBoxParent->GetTag() << ". FFD child box: " << FFDBoxChild->GetTag() <<"."<< endl;


}

void CSurfaceMovement::GetCartesianCoordCP(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBoxParent, CFreeFormDefBox *FFDBoxChild) {
  unsigned short iOrder, jOrder, kOrder, iDim;
  su2double *CartCoord, *ParamCoord;
    
  for (iOrder = 0; iOrder < FFDBoxChild->GetlOrder(); iOrder++)
    for (jOrder = 0; jOrder < FFDBoxChild->GetmOrder(); jOrder++)
      for (kOrder = 0; kOrder < FFDBoxChild->GetnOrder(); kOrder++) {
        ParamCoord = FFDBoxChild->GetParCoordControlPoints(iOrder, jOrder, kOrder);
        
        /*--- Clip the value of the parametric coordinates (just in case)  ---*/
        for (iDim = 0; iDim < 3; iDim++) {
          if (ParamCoord[iDim] >= 1.0) ParamCoord[iDim] = 1.0;
          if (ParamCoord[iDim] <= 0.0) ParamCoord[iDim] = 0.0;
        }

        CartCoord = FFDBoxParent->EvalCartesianCoord(ParamCoord);
        FFDBoxChild->SetCoordControlPoints(CartCoord, iOrder, jOrder, kOrder);
        FFDBoxChild->SetCoordControlPoints_Copy(CartCoord, iOrder, jOrder, kOrder);
        
      }
  
  if (rank == MASTER_NODE)
    cout << "Update cartesian coord (CP)   | FFD parent box: " << FFDBoxParent->GetTag() << ". FFD child box: " << FFDBoxChild->GetTag() <<"."<< endl;

}

void CSurfaceMovement::CheckFFDDimension(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox) {
  
  unsigned short iIndex, jIndex, kIndex, lDegree, mDegree, nDegree, iDV;
  bool OutOffLimits;
  bool polar = (config->GetFFD_CoordSystem() == POLAR);
  
  lDegree = FFDBox->GetlOrder()-1;
  mDegree = FFDBox->GetmOrder()-1;
  nDegree = FFDBox->GetnOrder()-1;
  
  OutOffLimits = false;
  for (iDV = 0; iDV < config->GetnDV(); iDV++) {
    if (config->GetFFDTag(iDV)== FFDBox->GetTag()){
      switch ( config->GetDesign_Variable(iDV) ) {
        case FFD_CONTROL_POINT_2D :
          if (polar) {
            iIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 1)));
            kIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 2)));
            if ((iIndex > lDegree) || (kIndex > nDegree)) OutOffLimits = true;
          }
          else {
            iIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 1)));
            jIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 2)));
            if ((iIndex > lDegree) || (jIndex > mDegree)) OutOffLimits = true;
          }
          break;
        case FFD_CAMBER :  case FFD_THICKNESS :
            iIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 1)));
            jIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 2)));
            if ((iIndex > lDegree) || (jIndex > mDegree)) OutOffLimits = true;
          break;
        case FFD_CAMBER_2D :  case FFD_THICKNESS_2D :
          iIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 1)));
          if (iIndex > lDegree) OutOffLimits = true;
          break;
        case FFD_CONTROL_POINT :  case FFD_NACELLE :
          iIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 1)));
          jIndex= SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 2)));
          kIndex = SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 3)));
          if ((iIndex > lDegree) || (jIndex > mDegree) || (kIndex > nDegree)) OutOffLimits = true;
          break;
        case FFD_GULL :  case FFD_TWIST :
          jIndex= SU2_TYPE::Int(fabs(config->GetParamDV(iDV, 1)));
          if (jIndex > mDegree) OutOffLimits = true;
          break;
      }
    }
  }
  
  if (rank == MASTER_NODE) {
    if (OutOffLimits) {
      char buf1[100], buf2[100];
      SPRINTF(buf1, "Design variables out off FFD limits (%u, %u, %u).\n", lDegree, mDegree, nDegree);
      SPRINTF(buf2, "Please check the ijk indices of the design variables.");
      SU2_MPI::Error(string(buf1) + string(buf2), CURRENT_FUNCTION);
    }
  }
  
  /*--- This barrier is important to guaranty that we will stop the software in a clean way ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
}

void CSurfaceMovement::CheckFFDIntersections(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox) {
  
  su2double Coord_0[] = {0,0,0}, Coord_1[] = {0,0,0};
  unsigned short index, iMarker, iNode, jNode, lDegree, mDegree, nDegree, iDim;
  unsigned long iElem, iPoint, jPoint;
  bool IPlane_Intersect_A = false, IPlane_Intersect_B = false;
  bool JPlane_Intersect_A = false, JPlane_Intersect_B = false;
  bool KPlane_Intersect_A = false, KPlane_Intersect_B = false;
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;

  unsigned short Kind_SU2 = config->GetKind_SU2();
  bool FFD_Symmetry_Plane = config->GetFFD_Symmetry_Plane();
  bool cylindrical = (config->GetFFD_CoordSystem() == CYLINDRICAL);
  bool spherical = (config->GetFFD_CoordSystem() == SPHERICAL);
  bool polar = (config->GetFFD_CoordSystem() == POLAR);
  bool cartesian = (config->GetFFD_CoordSystem() == CARTESIAN);
  
  lDegree = FFDBox->GetlOrder()-1;
  mDegree = FFDBox->GetmOrder()-1;
  nDegree = FFDBox->GetnOrder()-1;
  
  if (config->GetFFD_Continuity() != USER_INPUT) {
    
    /*--- Check intersection with plane i=0 ---*/
    
    su2double *IPlane_Coord_0_A = FFDBox->GetCoordControlPoints(0, 0, 0);
    su2double *IPlane_Coord_1_A = FFDBox->GetCoordControlPoints(0, 0, nDegree);
    su2double *IPlane_Coord_2_A = FFDBox->GetCoordControlPoints(0, mDegree, 0);
    
    su2double *IPlane_Coord_0_A_ = FFDBox->GetCoordControlPoints(0, mDegree, nDegree);
    su2double *IPlane_Coord_1_A_ = FFDBox->GetCoordControlPoints(0, mDegree, 0);
    su2double *IPlane_Coord_2_A_ = FFDBox->GetCoordControlPoints(0, 0, nDegree);
    
    /*--- Check intersection with plane i=lDegree ---*/
    
    su2double *IPlane_Coord_0_B = FFDBox->GetCoordControlPoints(lDegree, 0, 0);
    su2double *IPlane_Coord_1_B = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);
    su2double *IPlane_Coord_2_B = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);
    
    su2double *IPlane_Coord_0_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, nDegree);
    su2double *IPlane_Coord_1_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);
    su2double *IPlane_Coord_2_B_ = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);
    
    /*--- Check intersection with plane j=0 ---*/
    
    su2double *JPlane_Coord_0_A = FFDBox->GetCoordControlPoints(0,      0, 0);
    su2double *JPlane_Coord_1_A = FFDBox->GetCoordControlPoints(0,      0, nDegree);
    su2double *JPlane_Coord_2_A = FFDBox->GetCoordControlPoints(lDegree, 0, 0);
    
    su2double *JPlane_Coord_0_A_ = FFDBox->GetCoordControlPoints(lDegree, 0, nDegree);
    su2double *JPlane_Coord_1_A_ = FFDBox->GetCoordControlPoints(lDegree, 0, 0);
    su2double *JPlane_Coord_2_A_ = FFDBox->GetCoordControlPoints(0,      0, nDegree);
    
    /*--- Check intersection with plane j=mDegree ---*/
    
    su2double *JPlane_Coord_0_B = FFDBox->GetCoordControlPoints(0,      mDegree, 0);
    su2double *JPlane_Coord_1_B = FFDBox->GetCoordControlPoints(0,      mDegree, nDegree);
    su2double *JPlane_Coord_2_B = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);
    
    su2double *JPlane_Coord_0_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, nDegree);
    su2double *JPlane_Coord_1_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);
    su2double *JPlane_Coord_2_B_ = FFDBox->GetCoordControlPoints(0,      mDegree, nDegree);
    
    /*--- Check intersection with plane k=0 ---*/
    
    su2double *KPlane_Coord_0_A = FFDBox->GetCoordControlPoints(0,      0,      0);
    su2double *KPlane_Coord_1_A = FFDBox->GetCoordControlPoints(0,      mDegree, 0);
    su2double *KPlane_Coord_2_A = FFDBox->GetCoordControlPoints(lDegree, 0,      0);
    
    su2double *KPlane_Coord_0_A_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, 0);
    su2double *KPlane_Coord_1_A_ = FFDBox->GetCoordControlPoints(lDegree, 0,      0);
    su2double *KPlane_Coord_2_A_ = FFDBox->GetCoordControlPoints(0,      mDegree, 0);
    
    /*--- Check intersection with plane k=nDegree ---*/
    
    su2double *KPlane_Coord_0_B = FFDBox->GetCoordControlPoints(0,      0,      nDegree);
    su2double *KPlane_Coord_1_B = FFDBox->GetCoordControlPoints(0,      mDegree, nDegree);
    su2double *KPlane_Coord_2_B = FFDBox->GetCoordControlPoints(lDegree, 0,      nDegree);
    
    su2double *KPlane_Coord_0_B_ = FFDBox->GetCoordControlPoints(lDegree, mDegree, nDegree);
    su2double *KPlane_Coord_1_B_ = FFDBox->GetCoordControlPoints(lDegree, 0,      nDegree);
    su2double *KPlane_Coord_2_B_ = FFDBox->GetCoordControlPoints(0,      mDegree, nDegree);
    
    /*--- Loop over all the grid triangles ---*/
    
    IPlane_Intersect_A = false; IPlane_Intersect_B = false;
    JPlane_Intersect_A = false; JPlane_Intersect_B = false;
    KPlane_Intersect_A = false; KPlane_Intersect_B = false;
    
    /*--- Only the markers in the moving list ---*/
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      
      if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_CFD)) ||
          ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_DEF)) ||
          ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_GEO)) ||
          ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_DOT)) ||
          ((config->GetMarker_All_DV(iMarker) == YES) && (config->GetDirectDiff() == D_DESIGN))) {
        
        for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
          
          for (iNode = 0; iNode < geometry->bound[iMarker][iElem]->GetnNodes(); iNode++) {
            iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
            
            for (jNode = 0; jNode < geometry->bound[iMarker][iElem]->GetnNodes(); jNode++) {
              jPoint = geometry->bound[iMarker][iElem]->GetNode(jNode);
              
              if (jPoint > iPoint) {
                
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                  Coord_0[iDim] = geometry->node[iPoint]->GetCoord()[iDim];
                  Coord_1[iDim] = geometry->node[jPoint]->GetCoord()[iDim];
                }
                
                /*--- Write the coordinates in the right parametric system ---*/
                
                if (cylindrical) {
                  
                  X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
                  
                  Xbar =  Coord_0[0] - X_0; Ybar =  Coord_0[1] - Y_0; Zbar =  Coord_0[2] - Z_0;
                  
                  Coord_0[0] = sqrt(Ybar*Ybar + Zbar*Zbar);
                  Coord_0[1] = atan2(Zbar, Ybar); if (Coord_0[1] > PI_NUMBER/2.0) Coord_0[1] -= 2.0*PI_NUMBER;
                  Coord_0[2] = Xbar;
                  
                  Xbar =  Coord_1[0] - X_0; Ybar =  Coord_1[1] - Y_0; Zbar =  Coord_1[2] - Z_0;
                  
                  Coord_1[0] = sqrt(Ybar*Ybar + Zbar*Zbar);
                  Coord_1[1] = atan2(Zbar, Ybar); if (Coord_1[1] > PI_NUMBER/2.0) Coord_1[1] -= 2.0*PI_NUMBER;
                  Coord_1[2] = Xbar;
                  
                }
                
                else if (spherical || polar) {
                  
                  X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
                  
                  Xbar =  Coord_0[0] - X_0; Ybar =  Coord_0[1] - Y_0; Zbar =  Coord_0[2] - Z_0;
                  
                  Coord_0[0] = sqrt(Xbar*Xbar + Ybar*Ybar + Zbar*Zbar);
                  Coord_0[1] = atan2(Zbar, Ybar);  if (Coord_0[1] > PI_NUMBER/2.0) Coord_0[1] -= 2.0*PI_NUMBER;
                  Coord_0[2] = acos (Xbar/Coord_0[0]);
                  
                  Xbar =  Coord_1[0] - X_0; Ybar =  Coord_1[1] - Y_0; Zbar =  Coord_1[2] - Z_0;
                  
                  Coord_1[0] = sqrt(Xbar*Xbar + Ybar*Ybar + Zbar*Zbar);
                  Coord_1[1] = atan2(Zbar, Ybar);  if (Coord_1[1] > PI_NUMBER/2.0) Coord_1[1] -= 2.0*PI_NUMBER;
                  Coord_1[2] = acos(Xbar/Coord_1[0]);
                  
                }
                
                if (geometry->GetnDim() == 3) {
                  
                  if (!IPlane_Intersect_A) {
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, IPlane_Coord_0_A, IPlane_Coord_1_A, IPlane_Coord_2_A)) { IPlane_Intersect_A = true; }
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, IPlane_Coord_0_A_, IPlane_Coord_1_A_, IPlane_Coord_2_A_)) { IPlane_Intersect_A = true; }
                  }
                  
                  if (!IPlane_Intersect_B) {
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, IPlane_Coord_0_B, IPlane_Coord_1_B, IPlane_Coord_2_B)) { IPlane_Intersect_B = true; }
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, IPlane_Coord_0_B_, IPlane_Coord_1_B_, IPlane_Coord_2_B_)) { IPlane_Intersect_B = true; }
                  }
                  
                  if ((!JPlane_Intersect_A) && (!FFD_Symmetry_Plane)) {
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_A, JPlane_Coord_1_A, JPlane_Coord_2_A)) { JPlane_Intersect_A = true; }
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_A_, JPlane_Coord_1_A_, JPlane_Coord_2_A_)) { JPlane_Intersect_A = true; }
                  }
                  
                  if (cartesian) {
                    if ((!JPlane_Intersect_B) && (!FFD_Symmetry_Plane)) {
                      if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_B, JPlane_Coord_1_B, JPlane_Coord_2_B)) { JPlane_Intersect_B = true; }
                      if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_B_, JPlane_Coord_1_B_, JPlane_Coord_2_B_)) { JPlane_Intersect_B = true; }
                    }
                  }
                  else {
                    if (!JPlane_Intersect_B) {
                      if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_B, JPlane_Coord_1_B, JPlane_Coord_2_B)) { JPlane_Intersect_B = true; }
                      if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, JPlane_Coord_0_B_, JPlane_Coord_1_B_, JPlane_Coord_2_B_)) { JPlane_Intersect_B = true; }
                    }
                  }
                  
                  if (!KPlane_Intersect_A) {
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, KPlane_Coord_0_A, KPlane_Coord_1_A, KPlane_Coord_2_A)) { KPlane_Intersect_A = true; }
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, KPlane_Coord_0_A_, KPlane_Coord_1_A_, KPlane_Coord_2_A_)) { KPlane_Intersect_A = true; }
                  }
                  
                  if (!KPlane_Intersect_B) {
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, KPlane_Coord_0_B, KPlane_Coord_1_B, KPlane_Coord_2_B)) { KPlane_Intersect_B = true; }
                    if (geometry->SegmentIntersectsTriangle(Coord_0, Coord_1, KPlane_Coord_0_B_, KPlane_Coord_1_B_, KPlane_Coord_2_B_)) { KPlane_Intersect_B = true; }
                  }
                  
                } else {
                  
                  if (!IPlane_Intersect_A) {
                    if (geometry->SegmentIntersectsLine(Coord_0, Coord_1, IPlane_Coord_0_A, IPlane_Coord_2_A)) { IPlane_Intersect_A = true;}
                  }
                  if (!IPlane_Intersect_B) {
                    if (geometry->SegmentIntersectsLine(Coord_0, Coord_1, IPlane_Coord_0_B, IPlane_Coord_2_B)) { IPlane_Intersect_B = true;}
                  }
                  if (!JPlane_Intersect_A) {
                    if (geometry->SegmentIntersectsLine(Coord_0, Coord_1, JPlane_Coord_0_A, JPlane_Coord_2_A)) { JPlane_Intersect_A = true;}
                  }
                  if (!JPlane_Intersect_B) {
                    if (geometry->SegmentIntersectsLine(Coord_0, Coord_1, JPlane_Coord_0_B, JPlane_Coord_2_B)) { JPlane_Intersect_B = true;}
                  }
                }
              }
            }
          }
        }
      }
    }
    
    /*--- Comunicate the planes that interesect the surface ---*/
    
    unsigned short MyCode[6] = {0,0,0,0,0,0}, Code[6] = {0,0,0,0,0,0};
    
    if (IPlane_Intersect_A) MyCode[0] = 1;
    if (IPlane_Intersect_B) MyCode[1] = 1;
    if (JPlane_Intersect_A) MyCode[2] = 1;
    if (JPlane_Intersect_B) MyCode[3] = 1;
    if (KPlane_Intersect_A) MyCode[4] = 1;
    if (KPlane_Intersect_B) MyCode[5] = 1;
    
#ifdef HAVE_MPI
    
    /*--- Add SU2_MPI::Allreduce information using all the nodes ---*/
    
    SU2_MPI::Allreduce(&MyCode, &Code, 6, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
    
#else
    
    Code[0] = MyCode[0]; Code[1] = MyCode[1]; Code[2] = MyCode[2];
    Code[3] = MyCode[3]; Code[4] = MyCode[4]; Code[5] = MyCode[5];
    
#endif
    
    if (Code[0] != 0) IPlane_Intersect_A = true; else IPlane_Intersect_A = false;
    if (Code[1] != 0) IPlane_Intersect_B = true; else IPlane_Intersect_B = false;
    if (Code[2] != 0) JPlane_Intersect_A = true; else JPlane_Intersect_A = false;
    if (Code[3] != 0) JPlane_Intersect_B = true; else JPlane_Intersect_B = false;
    if (Code[4] != 0) KPlane_Intersect_A = true; else KPlane_Intersect_A = false;
    if (Code[5] != 0) KPlane_Intersect_B = true; else KPlane_Intersect_B = false;
    
    /*--- Screen output ---*/
    
    if (rank == MASTER_NODE) {
      
      if (IPlane_Intersect_A || IPlane_Intersect_B ||
          JPlane_Intersect_A || JPlane_Intersect_B ||
          KPlane_Intersect_A || KPlane_Intersect_B ) {
        
        cout << "The FFD planes ";
        
        if (cartesian) {
          if (IPlane_Intersect_A) cout << "i=0, ";
          if (IPlane_Intersect_B) cout << "i="<< lDegree << ", ";
          if (JPlane_Intersect_A) cout << "j=0, ";
          if (JPlane_Intersect_B) cout << "j="<< mDegree << ", ";
          if (KPlane_Intersect_A) cout << "k=0, ";
          if (KPlane_Intersect_B) cout << "k="<< nDegree << ", ";
        }
        else if (cylindrical) {
          if (IPlane_Intersect_A) cout << "r=0, ";
          if (IPlane_Intersect_B) cout << "r="<< lDegree << ", ";
          if (JPlane_Intersect_A) cout << "theta=0, ";
          if (JPlane_Intersect_B) cout << "theta="<< mDegree << ", ";
          if (KPlane_Intersect_A) cout << "z=0, ";
          if (KPlane_Intersect_B) cout << "z="<< nDegree << ", ";
        }
        else if (spherical) {
          if (IPlane_Intersect_A) cout << "r=0, ";
          if (IPlane_Intersect_B) cout << "r="<< lDegree << ", ";
          if (JPlane_Intersect_A) cout << "theta=0, ";
          if (JPlane_Intersect_B) cout << "theta="<< mDegree << ", ";
          if (KPlane_Intersect_A) cout << "phi=0, ";
          if (KPlane_Intersect_B) cout << "phi="<< nDegree << ", ";
        }
        else if (polar) {
          if (IPlane_Intersect_A) cout << "r=0, ";
          if (IPlane_Intersect_B) cout << "r="<< lDegree << ", ";
          if (KPlane_Intersect_A) cout << "theta=0, ";
          if (KPlane_Intersect_B) cout << "theta="<< nDegree << ", ";
        }
        
        cout << "intersect solid surfaces." << endl;
      }
      
    }
    
  }
  
  /*--- Fix the FFD planes based on the intersections with solid surfaces,
   and the continuity level, check that we have enough degree for the continuity
   that we are looking for ---*/
  
  if (config->GetFFD_Continuity() == USER_INPUT) {
    if (rank == MASTER_NODE)
      cout << "SU2 is fixing user's input planes." << endl;
    
    for (index = 0; index < config->GetnFFD_Fix_IDir(); index++)
      if ((config->GetFFD_Fix_IDir(index) <= lDegree) && (config->GetFFD_Fix_IDir(index) >= 0))
        FFDBox->Set_Fix_IPlane(config->GetFFD_Fix_IDir(index));
    for (index = 0; index < config->GetnFFD_Fix_JDir(); index++)
      if ((config->GetFFD_Fix_JDir(index) <= mDegree) && (config->GetFFD_Fix_JDir(index) >= 0))
        FFDBox->Set_Fix_JPlane(config->GetFFD_Fix_JDir(index));
    for (index = 0; index < config->GetnFFD_Fix_KDir(); index++)
      if ((config->GetFFD_Fix_KDir(index) <= nDegree) && (config->GetFFD_Fix_KDir(index) >= 0))
        FFDBox->Set_Fix_KPlane(config->GetFFD_Fix_KDir(index));
    
  }
  
  if (config->GetFFD_Continuity() == DERIVATIVE_NONE) {
    if (rank == MASTER_NODE)
      cout << "SU2 is fixing the planes to maintain a continuous surface." << endl;
    
    if (IPlane_Intersect_A) { FFDBox->Set_Fix_IPlane(0); }
    if (IPlane_Intersect_B) { FFDBox->Set_Fix_IPlane(lDegree); }
    if (JPlane_Intersect_A) { FFDBox->Set_Fix_JPlane(0); }
    if (JPlane_Intersect_B) { FFDBox->Set_Fix_JPlane(mDegree); }
    if (KPlane_Intersect_A) { FFDBox->Set_Fix_KPlane(0); }
    if (KPlane_Intersect_B) { FFDBox->Set_Fix_KPlane(nDegree); }
    
  }
  
  if (config->GetFFD_Continuity() == DERIVATIVE_1ST) {
    if (rank == MASTER_NODE)
      cout << "SU2 is fixing the planes to maintain a continuous 1st order derivative." << endl;
    
    if (IPlane_Intersect_A) { FFDBox->Set_Fix_IPlane(0); FFDBox->Set_Fix_IPlane(1); }
    if (IPlane_Intersect_B) { FFDBox->Set_Fix_IPlane(lDegree); FFDBox->Set_Fix_IPlane(lDegree-1); }
    if (JPlane_Intersect_A) { FFDBox->Set_Fix_JPlane(0); FFDBox->Set_Fix_JPlane(1); }
    if (JPlane_Intersect_B) { FFDBox->Set_Fix_JPlane(mDegree); FFDBox->Set_Fix_JPlane(mDegree-1); }
    if (KPlane_Intersect_A) { FFDBox->Set_Fix_KPlane(0); FFDBox->Set_Fix_KPlane(1); }
    if (KPlane_Intersect_B) { FFDBox->Set_Fix_KPlane(nDegree); FFDBox->Set_Fix_KPlane(nDegree-1); }
    
  }
  
  if (config->GetFFD_Continuity() == DERIVATIVE_2ND) {
    if (rank == MASTER_NODE)
      cout << "SU2 is fixing the planes to maintain a continuous 2nd order derivative." << endl;
    
    if ((IPlane_Intersect_A) && (lDegree > 1)) { FFDBox->Set_Fix_IPlane(0); FFDBox->Set_Fix_IPlane(1); FFDBox->Set_Fix_IPlane(2); }
    if ((IPlane_Intersect_B) && (lDegree > 1)) { FFDBox->Set_Fix_IPlane(lDegree); FFDBox->Set_Fix_IPlane(lDegree-1); FFDBox->Set_Fix_IPlane(lDegree-2); }
    if ((JPlane_Intersect_A) && (mDegree > 1)) { FFDBox->Set_Fix_JPlane(0); FFDBox->Set_Fix_JPlane(1); FFDBox->Set_Fix_JPlane(2); }
    if ((JPlane_Intersect_B) && (mDegree > 1)) { FFDBox->Set_Fix_JPlane(mDegree); FFDBox->Set_Fix_JPlane(mDegree-1); FFDBox->Set_Fix_JPlane(mDegree-2); }
    if ((KPlane_Intersect_A) && (nDegree > 1)) { FFDBox->Set_Fix_KPlane(0); FFDBox->Set_Fix_KPlane(1);FFDBox->Set_Fix_KPlane(2); }
    if ((KPlane_Intersect_B) && (nDegree > 1)) { FFDBox->Set_Fix_KPlane(nDegree); FFDBox->Set_Fix_KPlane(nDegree-1); FFDBox->Set_Fix_KPlane(nDegree-2); }
    
  }
  
}

void CSurfaceMovement::UpdateParametricCoord(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox) {
  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint, iSurfacePoints;
  su2double CartCoord[3] = {0.0,0.0,0.0}, *CartCoordNew, *CartCoordOld;
  su2double *ParamCoord, *var_coord, ParamCoordGuess[3] = {0.0,0.0,0.0};
  su2double MaxDiff, my_MaxDiff = 0.0, Diff;
      
  /*--- Recompute the parametric coordinates ---*/
  
  for (iSurfacePoints = 0; iSurfacePoints < FFDBox->GetnSurfacePoint(); iSurfacePoints++) {
    
    /*--- Get the marker of the surface point ---*/
    
    iMarker = FFDBox->Get_MarkerIndex(iSurfacePoints);
    
    if (config->GetMarker_All_DV(iMarker) == YES) {
      
      /*--- Get the vertex of the surface point ---*/
      
      iVertex = FFDBox->Get_VertexIndex(iSurfacePoints);
      iPoint = FFDBox->Get_PointIndex(iSurfacePoints);
  
      /*--- Get the parametric and cartesians coordinates of the 
       surface point (they don't mach) ---*/
      
      ParamCoord = FFDBox->Get_ParametricCoord(iSurfacePoints);
      
      /*--- Compute and set the cartesian coord using the variation computed 
       with the previous deformation ---*/
      
      var_coord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
      CartCoordOld = geometry->node[iPoint]->GetCoord();
      for (iDim = 0; iDim < 3; iDim++)
        CartCoord[iDim] = CartCoordOld[iDim] + var_coord[iDim];
      FFDBox->Set_CartesianCoord(CartCoord, iSurfacePoints);

      /*--- Find the parametric coordinate using as ParamCoordGuess the previous value ---*/
      
      ParamCoordGuess[0] = ParamCoord[0]; ParamCoordGuess[1] = ParamCoord[1]; ParamCoordGuess[2] = ParamCoord[2];
      ParamCoord = FFDBox->GetParametricCoord_Iterative(iPoint, CartCoord, ParamCoordGuess, config);
          
      /*--- Set the new value of the parametric coordinates ---*/
      
      FFDBox->Set_ParametricCoord(ParamCoord, iSurfacePoints);
      
      /*--- Compute the cartesian coordinates using the parametric coordinates 
       to check that everything is correct ---*/
      
      CartCoordNew = FFDBox->EvalCartesianCoord(ParamCoord);
      
      /*--- Compute max difference between original value and the recomputed value ---*/
      
      Diff = 0.0;
      for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
        Diff += (CartCoordNew[iDim]-CartCoord[iDim])*(CartCoordNew[iDim]-CartCoord[iDim]);
      Diff = sqrt(Diff);
      my_MaxDiff = max(my_MaxDiff, Diff);
        
    }
  }
    
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&my_MaxDiff, &MaxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  MaxDiff = my_MaxDiff;
#endif
  
  if (rank == MASTER_NODE) 
    cout << "Update parametric coord       | FFD box: " << FFDBox->GetTag() << ". Max Diff: " << MaxDiff <<"."<< endl;
  
}

su2double CSurfaceMovement::SetCartesianCoord(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, bool ResetDef) {
  
  su2double *CartCoordNew, Diff, my_MaxDiff = 0.0, MaxDiff,
  *ParamCoord, VarCoord[3] = {0.0, 0.0, 0.0}, CartCoordOld[3] = {0.0, 0.0, 0.0};
  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint, iSurfacePoints;
  
  bool cylindrical = (config->GetFFD_CoordSystem() == CYLINDRICAL);
  bool spherical = (config->GetFFD_CoordSystem() == SPHERICAL);
  bool polar = (config->GetFFD_CoordSystem() == POLAR);
  unsigned short nDim = geometry->GetnDim();
  
  /*--- Set to zero all the porints in VarCoord, this is important when we are dealing with different boxes
    because a loop over GetnSurfacePoint is no sufficient ---*/
  
  if (ResetDef) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
    }
  }
  
  /*--- Recompute the cartesians coordinates ---*/
  
  for (iSurfacePoints = 0; iSurfacePoints < FFDBox->GetnSurfacePoint(); iSurfacePoints++) {
    
    /*--- Get the marker of the surface point ---*/
    
    iMarker = FFDBox->Get_MarkerIndex(iSurfacePoints);
    
    if (config->GetMarker_All_DV(iMarker) == YES) {
      
      /*--- Get the vertex of the surface point ---*/
      
      iVertex = FFDBox->Get_VertexIndex(iSurfacePoints);
      iPoint = FFDBox->Get_PointIndex(iSurfacePoints);
      
      /*--- Set to zero the variation of the coordinates ---*/
      
      geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      
      /*--- Get the parametric coordinate of the surface point ---*/
      
      ParamCoord = FFDBox->Get_ParametricCoord(iSurfacePoints);
      
      /*--- Compute the new cartesian coordinate, and set the value in
       the FFDBox structure ---*/
      
      CartCoordNew = FFDBox->EvalCartesianCoord(ParamCoord);
      
      /*--- If polar coordinates, compute the cartesians from the polar value ---*/
      
      if (cylindrical) {
        
        su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
        X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
        
        Xbar = CartCoordNew[2];
        Ybar = CartCoordNew[0] * cos(CartCoordNew[1]);
        Zbar = CartCoordNew[0] * sin(CartCoordNew[1]);
        
        CartCoordNew[0] =  Xbar + X_0;  CartCoordNew[1] = Ybar + Y_0; CartCoordNew[2] = Zbar + Z_0;
        
      }
      else if (spherical || polar) {
        
        su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
        X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
        
        Xbar = CartCoordNew[0] * cos(CartCoordNew[2]);
        Ybar = CartCoordNew[0] * cos(CartCoordNew[1]) * sin(CartCoordNew[2]);
        Zbar = CartCoordNew[0] * sin(CartCoordNew[1]) * sin(CartCoordNew[2]);
        
        CartCoordNew[0] =  Xbar + X_0;  CartCoordNew[1] = Ybar + Y_0; CartCoordNew[2] = Zbar + Z_0;
        
      }
      
      FFDBox->Set_CartesianCoord(CartCoordNew, iSurfacePoints);
      
      /*--- Get the original cartesian coordinates of the surface point ---*/
      
      for (iDim = 0; iDim < nDim; iDim++) {
        CartCoordOld[iDim] = geometry->node[iPoint]->GetCoord(iDim);
      }
      
      /*--- Set the value of the variation of the coordinates ---*/
      
      Diff = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        VarCoord[iDim] = CartCoordNew[iDim] - CartCoordOld[iDim];
        if ((fabs(VarCoord[iDim]) <= EPS) && (config->GetDirectDiff() != D_DESIGN) && (!config->GetAD_Mode()))
          VarCoord[iDim] = 0.0;
        Diff += (VarCoord[iDim]*VarCoord[iDim]);
      }
      Diff = sqrt(Diff);
      
      my_MaxDiff = max(my_MaxDiff, Diff);
      
      /*--- Set the variation of the coordinates ---*/
      
      geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      
    }
  }
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&my_MaxDiff, &MaxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  MaxDiff = my_MaxDiff;
#endif
  
  if (rank == MASTER_NODE)
    cout << "Update cartesian coord        | FFD box: " << FFDBox->GetTag() << ". Max Diff: " << MaxDiff <<"."<< endl;
  
  return MaxDiff;

}


bool CSurfaceMovement::SetFFDCPChange_2D(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
    unsigned short iDV, bool ResetDef) {
  
  su2double movement[3] = {0.0,0.0,0.0}, Ampl;
  unsigned short index[3], i, j, iFFDBox, iPlane;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();
  bool polar = (config->GetFFD_CoordSystem() == POLAR);

  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/
  
  if (ResetDef == true) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
      ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }
  
  design_FFDBox = config->GetFFDTag(iDV);
  
  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    
    /*--- Compute deformation ---*/
    
    /*--- If we have only design value, than this value is the amplitude,
     * otherwise we have a general movement. ---*/
    
    if (config->GetnDV_Value(iDV) == 1) {
      
      Ampl = config->GetDV_Value(iDV)*Scale;
      
      if (polar){
        movement[0] = config->GetParamDV(iDV, 3)*Ampl;
        movement[1] = 0.0;
        movement[2] = config->GetParamDV(iDV, 4)*Ampl;
      }
      else {
        movement[0] = config->GetParamDV(iDV, 3)*Ampl;
        movement[1] = config->GetParamDV(iDV, 4)*Ampl;
        movement[2] = 0.0;
      }
      
    } else {
      if (polar){
        movement[0] = config->GetDV_Value(iDV, 0);
        movement[1] = 0.0;
        movement[2] = config->GetDV_Value(iDV, 1);
      }
      else {
      movement[0] = config->GetDV_Value(iDV, 0);
      movement[1] = config->GetDV_Value(iDV, 1);
      movement[2] = 0.0;
      }
      
    }
    
    if (polar){
       index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
       index[1] = 0;
      index[2] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
    }
    else {
       index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
       index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
      index[2] = 0;
    }
    
    /*--- Check that it is possible to move the control point ---*/
    
    for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
      if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) return false;
    }
    
    for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
      if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return false;
    }
    
    for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
      if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) return false;
    }
    
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        FFDBox->SetControlPoints(index, movement);
      }
    }
    
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1)) {
      for (j = 0; j < FFDBox->GetmOrder(); j++) {
        index[1] = j;
        FFDBox->SetControlPoints(index, movement);
      }
    }
    
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        for (j = 0; j < FFDBox->GetmOrder(); j++) {
          index[1] = j;
          FFDBox->SetControlPoints(index, movement);
        }
      }
    }
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1)) {
      
      FFDBox->SetControlPoints(index, movement);
    }
    
    /*--- Upper surface ---*/
    
    if (polar) index[1] = 1;
    else index[2] = 1;

    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        FFDBox->SetControlPoints(index, movement);
      }
    }
    
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1)) {
      for (j = 0; j < FFDBox->GetmOrder(); j++) {
        index[1] = j;
        FFDBox->SetControlPoints(index, movement);
      }
    }
    
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        for (j = 0; j < FFDBox->GetmOrder(); j++) {
          index[1] = j;
          FFDBox->SetControlPoints(index, movement);
        }
      }
    }
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1)) {
      
      FFDBox->SetControlPoints(index, movement);
    }
  }
  else {
    return false;
  }
  
  return true;
  
}

bool CSurfaceMovement::SetFFDCPChange(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
                                      unsigned short iDV, bool ResetDef) {
  
  su2double movement[3] = {0.0,0.0,0.0}, Ampl;
  unsigned short index[3], i, j, k, iPlane, iFFDBox;
  bool CheckIndex;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();
  
  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/
  
  if (ResetDef == true) {
    FFDBox->SetOriginalControlPoints();
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
      ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }
  
  design_FFDBox = config->GetFFDTag(iDV);
  
  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    
    /*--- Compute deformation ---*/
    
    /*--- If we have only design value, than this value is the amplitude,
     * otherwise we have a general movement. ---*/
    
    if (config->GetnDV_Value(iDV) == 1) {
      
      Ampl = config->GetDV_Value(iDV)*Scale;
      
      movement[0] = config->GetParamDV(iDV, 4)*Ampl;
      movement[1] = config->GetParamDV(iDV, 5)*Ampl;
      movement[2] = config->GetParamDV(iDV, 6)*Ampl;
      
    } else {
      
      movement[0] = config->GetDV_Value(iDV, 0);
      movement[1] = config->GetDV_Value(iDV, 1);
      movement[2] = config->GetDV_Value(iDV, 2);
      
    }
    
    index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
    index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
    index[2] = SU2_TYPE::Int(config->GetParamDV(iDV, 3));
    
    /*--- Check that it is possible to move the control point ---*/
    
    for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
      if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) return false;
    }
    
    for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
      if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return false;
    }
    
    for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
      if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) return false;
    }
    
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        
        CheckIndex = true;
        for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
          if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) CheckIndex = false;
        }
        
        if (CheckIndex) FFDBox->SetControlPoints(index, movement);
        
      }
    }
    
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
      for (j = 0; j < FFDBox->GetmOrder(); j++) {
        index[1] = j;
        
        CheckIndex = true;
        for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
          if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) CheckIndex = false;
        }
        
        if (CheckIndex) FFDBox->SetControlPoints(index, movement);
        
      }
    }
    
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) == -1)) {
      for (k = 0; k < FFDBox->GetnOrder(); k++) {
        index[2] = k;
        
        CheckIndex = true;
        for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
          if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) CheckIndex = false;
        }
        
        if (CheckIndex) FFDBox->SetControlPoints(index, movement);
        
      }
    }
    
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        for (j = 0; j < FFDBox->GetmOrder(); j++) {
          index[1] = j;
          FFDBox->SetControlPoints(index, movement);
        }
      }
    }
    
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) == -1)) {
      for (j = 0; j < FFDBox->GetmOrder(); j++) {
        index[1] = j;
        for (k = 0; k < FFDBox->GetnOrder(); k++) {
          index[2] = k;
          FFDBox->SetControlPoints(index, movement);
        }
      }
    }
    
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) == -1)) {
      for (i = 0; i < FFDBox->GetlOrder(); i++) {
        index[0] = i;
        for (k = 0; k < FFDBox->GetnOrder(); k++) {
          index[2] = k;
          FFDBox->SetControlPoints(index, movement);
        }
      }
    }
    
    if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
        (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
      FFDBox->SetControlPoints(index, movement);
    }
    
  }
  else {
    return false;
  }
  
  return true;
  
}

bool CSurfaceMovement::SetFFDGull(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
                                  unsigned short iDV, bool ResetDef) {
  
  su2double movement[3] = {0.0,0.0,0.0}, Ampl;
  unsigned short index[3], i, k, iPlane, iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();
  
  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/
  
  if (ResetDef == true) {
    FFDBox->SetOriginalControlPoints();
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
      ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }
  
  design_FFDBox = config->GetFFDTag(iDV);
  
  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    
    /*--- Compute deformation ---*/
    
    Ampl = config->GetDV_Value(iDV)*Scale;
    
    movement[0] = 0.0;
    movement[1] = 0.0;
    movement[2] = Ampl;
    
    /*--- Change the control points ---*/
    
    index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
    
    /*--- Check that it is possible to move the control point ---*/
    
    for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
      if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return false;
    }
    
    for (i = 0; i < FFDBox->GetlOrder(); i++) {
      index[0] = i;
      for (k = 0; k < FFDBox->GetnOrder(); k++) {
        index[2] = k;
        FFDBox->SetControlPoints(index, movement);
      }
    }
    
  }
  else {
    return false;
  }
  
  return true;
  
}

bool CSurfaceMovement::SetFFDNacelle(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
                                     unsigned short iDV, bool ResetDef) {
  
  su2double movement[3] = {0.0,0.0,0.0}, Ampl;
  unsigned short index[3], i, j, k, iPlane, iFFDBox, Theta, ThetaMax;
  string design_FFDBox;
  bool SameCP = false;
  su2double Scale = config->GetOpt_RelaxFactor();
  
  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/
  
  if (ResetDef == true) {
    FFDBox->SetOriginalControlPoints();
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
      ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }
  
  design_FFDBox = config->GetFFDTag(iDV);
  
  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    
    /*--- Compute deformation ---*/
    
    Ampl = config->GetDV_Value(iDV)*Scale;
    
    movement[0] = config->GetParamDV(iDV, 4)*Ampl;
    movement[1] = 0.0;
    movement[2] = config->GetParamDV(iDV, 5)*Ampl;
    
    index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
    index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
    index[2] = SU2_TYPE::Int(config->GetParamDV(iDV, 3));
    if (index[1] == SU2_TYPE::Int(FFDBox->GetmOrder()) - index[1] -1) SameCP = true;
    
    ThetaMax = 2;
    if (SameCP) ThetaMax = 1;
    
    for (Theta = 0; Theta < ThetaMax; Theta++) {
      
      if (Theta == 1) index[1] = SU2_TYPE::Int(FFDBox->GetmOrder()) - index[1] -1;
      
      /*--- Check that it is possible to move the control point ---*/
      
      for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
        if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) return false;
      }
      
      for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
        if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return false;
      }
      
      for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
        if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) return false;
      }
      
      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
        for (i = 0; i < FFDBox->GetlOrder(); i++) {
          index[0] = i;
          FFDBox->SetControlPoints(index, movement);
        }
      }
      
      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
        for (j = 0; j < FFDBox->GetmOrder(); j++) {
          index[1] = j;
          FFDBox->SetControlPoints(index, movement);
        }
      }
      
      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) == -1)) {
        for (k = 0; k < FFDBox->GetnOrder(); k++) {
          index[2] = k;
          FFDBox->SetControlPoints(index, movement);
        }
      }
      
      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
        for (i = 0; i < FFDBox->GetlOrder(); i++) {
          index[0] = i;
          for (j = 0; j < FFDBox->GetmOrder(); j++) {
            index[1] = j;
            FFDBox->SetControlPoints(index, movement);
          }
        }
      }
      
      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) == -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) == -1)) {
        for (j = 0; j < FFDBox->GetmOrder(); j++) {
          index[1] = j;
          for (k = 0; k < FFDBox->GetnOrder(); k++) {
            index[2] = k;
            FFDBox->SetControlPoints(index, movement);
          }
        }
      }
      
      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) == -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) == -1)) {
        for (i = 0; i < FFDBox->GetlOrder(); i++) {
          index[0] = i;
          for (k = 0; k < FFDBox->GetnOrder(); k++) {
            index[2] = k;
            FFDBox->SetControlPoints(index, movement);
          }
        }
      }
      
      if ((SU2_TYPE::Int(config->GetParamDV(iDV, 1)) != -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 2)) != -1) &&
          (SU2_TYPE::Int(config->GetParamDV(iDV, 3)) != -1)) {
        FFDBox->SetControlPoints(index, movement);
      }
    }
    
  }
  else {
    return false;
  }
  
  return true;
  
}

bool CSurfaceMovement::SetFFDCamber_2D(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
                                       unsigned short iDV, bool ResetDef) {
  
  su2double Ampl, movement[3] = {0.0,0.0,0.0};
  unsigned short index[3], kIndex, iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();
  
  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/
  
  if (ResetDef == true) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
      ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }
  
  design_FFDBox = config->GetFFDTag(iDV);
  
  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    
    for (kIndex = 0; kIndex < 2; kIndex++) {
      
      Ampl = config->GetDV_Value(iDV)*Scale;
            
      movement[0] = 0.0;
      if (kIndex == 0) movement[1] = Ampl;
      else movement[1] = Ampl;
      movement[2] = 0.0;
      
      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1)); index[1] = kIndex; index[2] = 0;
      FFDBox->SetControlPoints(index, movement);
      
      index[2] = 1;
      FFDBox->SetControlPoints(index, movement);
      
    }
    
  }
  else {
    return false;
  }
  
  return true;
  
}

bool CSurfaceMovement::SetFFDThickness_2D(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
                                          unsigned short iDV, bool ResetDef) {
  
  su2double Ampl, movement[3]= {0.0,0.0,0.0};
  unsigned short index[3], kIndex, iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();
  
  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/
  
  if (ResetDef == true) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
      ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }
  
  design_FFDBox = config->GetFFDTag(iDV);
  
  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
        
    for (kIndex = 0; kIndex < 2; kIndex++) {
      
      Ampl = config->GetDV_Value(iDV)*Scale;
      
      movement[0] = 0.0;
      if (kIndex == 0) movement[1] = -Ampl;
      else movement[1] = Ampl;
      movement[2] = 0.0;
      
      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1)); index[1] = kIndex; index[2] = 0;
      FFDBox->SetControlPoints(index, movement);
      
      index[2] = 1;
      FFDBox->SetControlPoints(index, movement);
      
    }
    
  }
  else {
    return false;
  }
  
  return true;
  
}

bool CSurfaceMovement::SetFFDTwist_2D(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
                                      unsigned short iDV, bool ResetDef) {
  
  return true;
  
}

bool CSurfaceMovement::SetFFDCamber(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
                                    unsigned short iDV, bool ResetDef) {
  
  su2double Ampl, movement[3] = {0.0,0.0,0.0};
  unsigned short index[3], kIndex, iPlane, iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();
  
  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/
  
  if (ResetDef == true) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
      ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }
  
  design_FFDBox = config->GetFFDTag(iDV);
  
  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    
    /*--- Check that it is possible to move the control point ---*/
    
    for (kIndex = 0; kIndex < 2; kIndex++) {
      
      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
      index[2] = kIndex;
      
      for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
        if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) return false;
      }
      
      for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
        if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return false;
      }
      
      for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
        if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) return false;
      }
      
    }
    
    for (kIndex = 0; kIndex < 2; kIndex++) {
            
      Ampl = config->GetDV_Value(iDV)*Scale;
            
      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2)); 
      index[2] = kIndex;
      
      movement[0] = 0.0; movement[1] = 0.0; 
      if (kIndex == 0) movement[2] = Ampl;
      else movement[2] = Ampl;
      
      FFDBox->SetControlPoints(index, movement);
      
    }
    
  }
  else {
    return false;
  }
  
  return true;
  
}

void CSurfaceMovement::SetFFDAngleOfAttack(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
                                           unsigned short iDV, bool ResetDef) {
  
  su2double Scale = config->GetOpt_RelaxFactor();

  su2double Ampl = config->GetDV_Value(iDV)*Scale;
  
  config->SetAoA_Offset(Ampl);
  
}

bool CSurfaceMovement::SetFFDThickness(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
                                       unsigned short iDV, bool ResetDef) {
  
  su2double Ampl, movement[3] = {0.0,0.0,0.0};
  unsigned short index[3], kIndex, iPlane, iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();
  
  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/
  
  if (ResetDef == true) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
      ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }
  
  design_FFDBox = config->GetFFDTag(iDV);
  
  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    
    /*--- Check that it is possible to move the control point ---*/
    
    for (kIndex = 0; kIndex < 2; kIndex++) {
      
      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
      index[2] = kIndex;
      
      for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
        if (index[0] == FFDBox->Get_Fix_IPlane(iPlane)) return false;
      }
      
      for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
        if (index[1] == FFDBox->Get_Fix_JPlane(iPlane)) return false;
      }
      
      for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
        if (index[2] == FFDBox->Get_Fix_KPlane(iPlane)) return false;
      }
      
    }
    
    
    for (kIndex = 0; kIndex < 2; kIndex++) {
      
      Ampl = config->GetDV_Value(iDV)*Scale;
      
      index[0] = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      index[1] = SU2_TYPE::Int(config->GetParamDV(iDV, 2));
      index[2] = kIndex;
      
      movement[0] = 0.0; movement[1] = 0.0;
      if (kIndex == 0) movement[2] = -Ampl;
      else movement[2] = Ampl;
      
      FFDBox->SetControlPoints(index, movement);
      
    }
    
  }
  else {
    return false;
  }
  
  return true;
  
}

bool CSurfaceMovement::SetFFDTwist(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
                                   unsigned short iDV, bool ResetDef) {
  
  unsigned short iOrder, jOrder, kOrder;
  su2double  x, y, z, movement[3], Segment_P0[3], Segment_P1[3], Plane_P0[3], Plane_Normal[3],
  Variable_P0, Variable_P1, Intersection[3], Variable_Interp;
  unsigned short index[3], iPlane, iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();
  
  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/
  
  if (ResetDef == true) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
      ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }
  
  design_FFDBox = config->GetFFDTag(iDV);
  
  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    
    /*--- Check that it is possible to move the control point ---*/
    
    jOrder = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
    for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_JPlane(); iPlane++) {
      if (jOrder == FFDBox->Get_Fix_JPlane(iPlane)) return false;
    }
    
    /*--- Line plane intersection to find the origin of rotation ---*/
    
    Segment_P0[0] = config->GetParamDV(iDV, 2);
    Segment_P0[1] = config->GetParamDV(iDV, 3);
    Segment_P0[2] = config->GetParamDV(iDV, 4);
    
    Segment_P1[0] = config->GetParamDV(iDV, 5);
    Segment_P1[1] = config->GetParamDV(iDV, 6);
    Segment_P1[2] = config->GetParamDV(iDV, 7);
    
    iOrder = 0;
    jOrder = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
    kOrder = 0;
    su2double *coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
    Plane_P0[0] = coord[0]; Plane_P0[1] = coord[1]; Plane_P0[2] = coord[2];
    Plane_Normal[0] = 0.0; Plane_Normal[1] = 1.0; Plane_Normal[2] = 0.0;
    
    Variable_P0 = 0.0; Variable_P1 = 0.0;
    
    Intersection[0] = 0.0; Intersection[1] = 0.0;  Intersection[2] = 0.0;
    
    bool result = geometry->SegmentIntersectsPlane(Segment_P0, Segment_P1, Variable_P0, Variable_P1,
                                                   Plane_P0, Plane_Normal, Intersection, Variable_Interp);
    
    if (result) {
      
      /*--- xyz-coordinates of a point on the line of rotation. ---*/
      
      su2double a = Intersection[0];
      su2double b = Intersection[1];
      su2double c = Intersection[2];
      
      /*--- xyz-coordinate of the line's direction vector. ---*/
      
      su2double u = Plane_Normal[0];
      su2double v = Plane_Normal[1];
      su2double w = Plane_Normal[2];
      
      /*--- The angle of rotation is computed based on a characteristic length of the wing,
       otherwise it is difficult to compare with other length based design variables. ---*/
      
      su2double RefLength = config->GetRefLength();
      su2double theta = atan(config->GetDV_Value(iDV)*Scale/RefLength);
      
      /*--- An intermediate value used in computations. ---*/
      
      su2double u2=u*u; su2double v2=v*v; su2double w2=w*w;
      su2double l2 = u2 + v2 + w2; su2double l = sqrt(l2);
      su2double cosT; su2double sinT;
      
      /*--- Change the value of the control point if move is true ---*/
      
      jOrder = SU2_TYPE::Int(config->GetParamDV(iDV, 1));
      for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++)
        for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
          index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
          su2double *coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
          x = coord[0]; y = coord[1]; z = coord[2];
          
          cosT = cos(theta);
          sinT = sin(theta);
          
          movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
          + (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
          + l*(-c*v + b*w - w*y + v*z)*sinT;
          movement[0] = movement[0]/l2 - x;
          
          movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z)
          + (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
          + l*(c*u - a*w + w*x - u*z)*sinT;
          movement[1] = movement[1]/l2 - y;
          
          movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z)
          + (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
          + l*(-b*u + a*v - v*x + u*y)*sinT;
          movement[2] = movement[2]/l2 - z;
          
          /*--- Check that it is possible to move the control point ---*/
          
          for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_IPlane(); iPlane++) {
            if (iOrder == FFDBox->Get_Fix_IPlane(iPlane)) {
              movement[0] = 0.0; movement[1] = 0.0; movement[2] = 0.0;
            }
          }
          
          for (iPlane = 0 ; iPlane < FFDBox->Get_nFix_KPlane(); iPlane++) {
            if (kOrder == FFDBox->Get_Fix_KPlane(iPlane)) {
              movement[0] = 0.0; movement[1] = 0.0; movement[2] = 0.0;
            }
          }

          FFDBox->SetControlPoints(index, movement);
          
        }
      
    }
    
  }
  else {
    return false;
  }
  
  return true;
  
}

bool CSurfaceMovement::SetFFDRotation(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
                                      unsigned short iDV, bool ResetDef) {
  
  unsigned short iOrder, jOrder, kOrder;
  su2double movement[3] = {0.0,0.0,0.0}, x, y, z;
  unsigned short index[3], iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();
  
  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/
  
  if (ResetDef == true) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
      ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }
  
  design_FFDBox = config->GetFFDTag(iDV);
  
  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    
    /*--- xyz-coordinates of a point on the line of rotation. ---*/
    
    su2double a = config->GetParamDV(iDV, 1);
    su2double b = config->GetParamDV(iDV, 2);
    su2double c = config->GetParamDV(iDV, 3);
    
    /*--- xyz-coordinate of the line's direction vector. ---*/
    
    su2double u = config->GetParamDV(iDV, 4)-config->GetParamDV(iDV, 1);
    su2double v = config->GetParamDV(iDV, 5)-config->GetParamDV(iDV, 2);
    su2double w = config->GetParamDV(iDV, 6)-config->GetParamDV(iDV, 3);
    
    /*--- The angle of rotation. ---*/
    
    su2double theta = config->GetDV_Value(iDV)*Scale*PI_NUMBER/180.0;
    
    /*--- An intermediate value used in computations. ---*/
    
    su2double u2=u*u; su2double v2=v*v; su2double w2=w*w;
    su2double cosT = cos(theta); su2double sinT = sin(theta);
    su2double l2 = u2 + v2 + w2; su2double l = sqrt(l2);
    
    /*--- Change the value of the control point if move is true ---*/
    
    for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++)
      for (jOrder = 0; jOrder < FFDBox->GetmOrder(); jOrder++)
        for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
          index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
          su2double *coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
          x = coord[0]; y = coord[1]; z = coord[2];
          movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
          + (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
          + l*(-c*v + b*w - w*y + v*z)*sinT;
          movement[0] = movement[0]/l2 - x;
          
          movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z)
          + (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
          + l*(c*u - a*w + w*x - u*z)*sinT;
          movement[1] = movement[1]/l2 - y;
          
          movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z)
          + (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
          + l*(-b*u + a*v - v*x + u*y)*sinT;
          movement[2] = movement[2]/l2 - z;
          
          FFDBox->SetControlPoints(index, movement);
          
        }
  }
  else {
    return false;
  }
  
  return true;
  
}

bool CSurfaceMovement::SetFFDControl_Surface(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, CFreeFormDefBox **ResetFFDBox,
                                             unsigned short iDV, bool ResetDef) {
  
  unsigned short iOrder, jOrder, kOrder;
  su2double movement[3] = {0.0,0.0,0.0}, x, y, z;
  unsigned short index[3], iFFDBox;
  string design_FFDBox;
  su2double Scale = config->GetOpt_RelaxFactor();
  
  /*--- Set control points to its original value (even if the
   design variable is not in this box) ---*/
  
  if (ResetDef == true) {
    for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
      ResetFFDBox[iFFDBox]->SetOriginalControlPoints();
  }
  
  design_FFDBox = config->GetFFDTag(iDV);
  
  if (design_FFDBox.compare(FFDBox->GetTag()) == 0) {
    
    /*--- xyz-coordinates of a point on the line of rotation. ---*/
    
    su2double a = config->GetParamDV(iDV, 1);
    su2double b = config->GetParamDV(iDV, 2);
    su2double c = config->GetParamDV(iDV, 3);
    
    /*--- xyz-coordinate of the line's direction vector. ---*/
    
    su2double u = config->GetParamDV(iDV, 4)-config->GetParamDV(iDV, 1);
    su2double v = config->GetParamDV(iDV, 5)-config->GetParamDV(iDV, 2);
    su2double w = config->GetParamDV(iDV, 6)-config->GetParamDV(iDV, 3);
    
    /*--- The angle of rotation. ---*/
    
    su2double theta = -config->GetDV_Value(iDV)*Scale*PI_NUMBER/180.0;
    
    /*--- An intermediate value used in computations. ---*/
    
    su2double u2=u*u; su2double v2=v*v; su2double w2=w*w;
    su2double cosT = cos(theta); su2double sinT = sin(theta);
    su2double l2 = u2 + v2 + w2; su2double l = sqrt(l2);
    
    /*--- Change the value of the control point if move is true ---*/
    
    for (iOrder = 0; iOrder < FFDBox->GetlOrder()-2; iOrder++)
      for (jOrder = 2; jOrder < FFDBox->GetmOrder()-2; jOrder++)
        for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
          index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
          su2double *coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
          x = coord[0]; y = coord[1]; z = coord[2];
          movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
          + (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
          + l*(-c*v + b*w - w*y + v*z)*sinT;
          movement[0] = movement[0]/l2 - x;
          
          movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z)
          + (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
          + l*(c*u - a*w + w*x - u*z)*sinT;
          movement[1] = movement[1]/l2 - y;
          
          movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z)
          + (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
          + l*(-b*u + a*v - v*x + u*y)*sinT;
          movement[2] = movement[2]/l2 - z;
          
          FFDBox->SetControlPoints(index, movement);
          
        }
  }
  else {
    return false;
  }
  
  return true;
  
}

void CSurfaceMovement::SetAngleOfAttack(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
  
  su2double Scale = config->GetOpt_RelaxFactor();
  su2double Ampl = config->GetDV_Value(iDV)*Scale;
  config->SetAoA_Offset(Ampl);
  
}

void CSurfaceMovement::SetHicksHenne(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3] = {0.0,0.0,0.0}, VarCoord_[3] = {0.0,0.0,0.0}, *Coord_, *Normal_, ek, fk,
      Coord[3] = {0.0,0.0,0.0}, Normal[3] = {0.0,0.0,0.0},
  TPCoord[2] = {0.0, 0.0}, LPCoord[2] = {0.0, 0.0}, Distance, Chord, AoA, ValCos, ValSin;
  
  bool upper = true;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/
  
  if ((iDV == 0) || (ResetDef == true)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
  }
  
  /*--- Compute the angle of attack to apply the deformation ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      Coord_ = boundary->vertex[iMarker][0]->GetCoord();
      TPCoord[0] = Coord_[0]; TPCoord[1] = Coord_[1];
      for (iVertex = 1; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        if (Coord_[0] > TPCoord[0]) { TPCoord[0] = Coord_[0]; TPCoord[1] = Coord_[1]; }
      }
    }
  }
  
#ifdef HAVE_MPI

  int iProcessor, nProcessor = size;
  su2double *Buffer_Send_Coord, *Buffer_Receive_Coord;
  
  Buffer_Receive_Coord = new su2double [nProcessor*2];
  Buffer_Send_Coord = new su2double [2];
  
  Buffer_Send_Coord[0] = TPCoord[0]; Buffer_Send_Coord[1] = TPCoord[1];

  SU2_MPI::Allgather(Buffer_Send_Coord, 2, MPI_DOUBLE, Buffer_Receive_Coord, 2, MPI_DOUBLE, MPI_COMM_WORLD);

  TPCoord[0] = Buffer_Receive_Coord[0]; TPCoord[1] = Buffer_Receive_Coord[1];
  for (iProcessor = 1; iProcessor < nProcessor; iProcessor++) {
    Coord[0] = Buffer_Receive_Coord[iProcessor*2 + 0];
    Coord[1] = Buffer_Receive_Coord[iProcessor*2 + 1];
    if (Coord[0] > TPCoord[0]) { TPCoord[0] = Coord[0]; TPCoord[1] = Coord[1]; }
  }
  
  delete[] Buffer_Send_Coord;   delete[] Buffer_Receive_Coord;
  
#endif


  Chord = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        Distance = sqrt(pow(Coord_[0] - TPCoord[0], 2.0) + pow(Coord_[1] - TPCoord[1], 2.0));
        if (Chord < Distance) { Chord = Distance; LPCoord[0] = Coord_[0]; LPCoord[1] = Coord_[1]; }
      }
    }
  }
  
#ifdef HAVE_MPI
   
  Buffer_Receive_Coord = new su2double [nProcessor*2];
  Buffer_Send_Coord = new su2double [2];
  
  Buffer_Send_Coord[0] = LPCoord[0]; Buffer_Send_Coord[1] = LPCoord[1];

  SU2_MPI::Allgather(Buffer_Send_Coord, 2, MPI_DOUBLE, Buffer_Receive_Coord, 2, MPI_DOUBLE, MPI_COMM_WORLD);
  
  Chord = 0.0;
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    Coord[0] = Buffer_Receive_Coord[iProcessor*2 + 0];
    Coord[1] = Buffer_Receive_Coord[iProcessor*2 + 1];
    Distance = sqrt(pow(Coord[0] - TPCoord[0], 2.0) + pow(Coord[1] - TPCoord[1], 2.0));
    if (Chord < Distance) { Chord = Distance; LPCoord[0] = Coord[0]; LPCoord[1] = Coord[1]; }
  }
  
  delete[] Buffer_Send_Coord;   delete[] Buffer_Receive_Coord;
  
#endif
  
  AoA = atan((LPCoord[1] - TPCoord[1]) / (TPCoord[0] - LPCoord[0]))*180/PI_NUMBER;
  
  /*--- WARNING: AoA currently overwritten to zero. ---*/
  AoA = 0.0;

  /*--- Perform multiple airfoil deformation ---*/
  
  su2double Ampl = config->GetDV_Value(iDV)*Scale;
  su2double xk = config->GetParamDV(iDV, 1);
  const su2double t2 = 3.0;
  
  if (config->GetParamDV(iDV, 0) == NO) { upper = false; }
  if (config->GetParamDV(iDV, 0) == YES) { upper = true; }
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
      
      if (config->GetMarker_All_DV(iMarker) == YES) {
        
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        Normal_ = boundary->vertex[iMarker][iVertex]->GetNormal();
        
        /*--- The Hicks Henne bump functions should be applied to a basic airfoil without AoA,
         and unitary chord, a tranformation is required ---*/
        
        ValCos = cos(AoA*PI_NUMBER/180.0);
        ValSin = sin(AoA*PI_NUMBER/180.0);
        
        Coord[0] = Coord_[0]*ValCos - Coord_[1]*ValSin;
        Coord[0] = max(0.0, Coord[0]); // Coord x should be always positive
        Coord[1] = Coord_[1]*ValCos + Coord_[0]*ValSin;

        Normal[0] = Normal_[0]*ValCos - Normal_[1]*ValSin;
        Normal[1] = Normal_[1]*ValCos + Normal_[0]*ValSin;

        /*--- Bump computation ---*/

        ek = log10(0.5)/log10(xk);
        if (Coord[0] > 10*EPS) fk = pow( sin( PI_NUMBER * pow(Coord[0], ek) ), t2);
        else fk = 0.0;

        /*--- Upper and lower surface ---*/

        if (( upper) && (Normal[1] > 0)) { VarCoord[1] =  Ampl*fk; }
        if ((!upper) && (Normal[1] < 0)) { VarCoord[1] = -Ampl*fk; }

      }
      
      /*--- Apply the transformation to the coordinate variation ---*/
      
      ValCos = cos(-AoA*PI_NUMBER/180.0);
      ValSin = sin(-AoA*PI_NUMBER/180.0);
      
      VarCoord_[0] = VarCoord[0]*ValCos - VarCoord[1]*ValSin;
      VarCoord_[1] = VarCoord[1]*ValCos + VarCoord[0]*ValSin;

      boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord_);
      
    }
  }
  
}

void CSurfaceMovement::SetSurface_Bump(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3] = {0.0,0.0,0.0}, ek, fk, *Coord, xCoord;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/

  if ((iDV == 0) || (ResetDef == true)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
  }

  /*--- Perform multiple airfoil deformation ---*/

  su2double Ampl = config->GetDV_Value(iDV)*Scale;
  su2double x_start = config->GetParamDV(iDV, 0);
  su2double x_end = config->GetParamDV(iDV, 1);
  su2double BumpSize = x_end - x_start;
  su2double BumpLoc = x_start;
  su2double xk = config->GetParamDV(iDV, 2);
  const su2double t2 = 3.0;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;

      if (config->GetMarker_All_DV(iMarker) == YES) {

        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();

        xCoord = (Coord[0] - BumpLoc);
        ek = log10(0.5)/log10((xk-BumpLoc+EPS)/BumpSize);
        if (xCoord > 0.0) fk = pow( sin( PI_NUMBER * pow((xCoord+EPS)/BumpSize, ek)), t2);
        else fk = 0.0;

        if ((xCoord <= 0.0) || (xCoord >= BumpSize)) VarCoord[1] =  0.0;
        else { VarCoord[1] =  Ampl*fk; }

      }

      boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);

    }
  }

}

void CSurfaceMovement::SetCST(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3] = {0.0,0.0,0.0}, VarCoord_[3] = {0.0,0.0,0.0}, *Coord_, *Normal_, fk,
    Coord[3] = {0.0,0.0,0.0}, Normal[3] = {0.0,0.0,0.0},
    TPCoord[2] = {0.0, 0.0}, LPCoord[2] = {0.0, 0.0}, Distance, Chord, AoA, ValCos, ValSin;
  
  bool upper = true;
  su2double Scale = config->GetOpt_RelaxFactor();

  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/

  if ((iDV == 0) || (ResetDef == true)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
  }
  
    /*--- Compute the angle of attack to apply the deformation ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      Coord_ = boundary->vertex[iMarker][0]->GetCoord();
      TPCoord[0] = Coord_[0]; TPCoord[1] = Coord_[1];
      for (iVertex = 1; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        if (Coord_[0] > TPCoord[0]) { TPCoord[0] = Coord_[0]; TPCoord[1] = Coord_[1]; }
      }
    }
  }
  
#ifdef HAVE_MPI

  int iProcessor, nProcessor = size;
  su2double *Buffer_Send_Coord, *Buffer_Receive_Coord;
  
  Buffer_Receive_Coord = new su2double [nProcessor*2];
  Buffer_Send_Coord = new su2double [2];
  
  Buffer_Send_Coord[0] = TPCoord[0]; Buffer_Send_Coord[1] = TPCoord[1];

  SU2_MPI::Allgather(Buffer_Send_Coord, 2, MPI_DOUBLE, Buffer_Receive_Coord, 2, MPI_DOUBLE, MPI_COMM_WORLD);

  TPCoord[0] = Buffer_Receive_Coord[0]; TPCoord[1] = Buffer_Receive_Coord[1];
  for (iProcessor = 1; iProcessor < nProcessor; iProcessor++) {
    Coord[0] = Buffer_Receive_Coord[iProcessor*2 + 0];
    Coord[1] = Buffer_Receive_Coord[iProcessor*2 + 1];
    if (Coord[0] > TPCoord[0]) { TPCoord[0] = Coord[0]; TPCoord[1] = Coord[1]; }
  }
  
  delete[] Buffer_Send_Coord;   delete[] Buffer_Receive_Coord;
  
#endif


  Chord = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        Distance = sqrt(pow(Coord_[0] - TPCoord[0], 2.0) + pow(Coord_[1] - TPCoord[1], 2.0));
        if (Chord < Distance) { Chord = Distance; LPCoord[0] = Coord_[0]; LPCoord[1] = Coord_[1]; }
      }
    }
  }
  
#ifdef HAVE_MPI
   
  Buffer_Receive_Coord = new su2double [nProcessor*2];
  Buffer_Send_Coord = new su2double [2];
  
  Buffer_Send_Coord[0] = LPCoord[0]; Buffer_Send_Coord[1] = LPCoord[1];

  SU2_MPI::Allgather(Buffer_Send_Coord, 2, MPI_DOUBLE, Buffer_Receive_Coord, 2, MPI_DOUBLE, MPI_COMM_WORLD);
  
  Chord = 0.0;
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    Coord[0] = Buffer_Receive_Coord[iProcessor*2 + 0];
    Coord[1] = Buffer_Receive_Coord[iProcessor*2 + 1];
    Distance = sqrt(pow(Coord[0] - TPCoord[0], 2.0) + pow(Coord[1] - TPCoord[1], 2.0));
    if (Chord < Distance) { Chord = Distance; LPCoord[0] = Coord[0]; LPCoord[1] = Coord[1]; }
  }
  
  delete[] Buffer_Send_Coord;   delete[] Buffer_Receive_Coord;
  
#endif
  
  AoA = atan((LPCoord[1] - TPCoord[1]) / (TPCoord[0] - LPCoord[0]))*180/PI_NUMBER;
  
  /*--- WARNING: AoA currently overwritten to zero. ---*/
  AoA = 0.0;

  /*--- Perform multiple airfoil deformation ---*/
    
  su2double Ampl = config->GetDV_Value(iDV)*Scale;
  su2double KulfanNum = config->GetParamDV(iDV, 1) - 1.0;
  su2double maxKulfanNum = config->GetParamDV(iDV, 2) - 1.0;
  if (KulfanNum < 0) {
    std::cout << "Warning: Kulfan number should be greater than 1." << std::endl;
  }
  if (KulfanNum > maxKulfanNum) {
    std::cout << "Warning: Kulfan number should be less than provided maximum." << std::endl;
  }

  if (config->GetParamDV(iDV, 0) == NO) { upper = false;}
  if (config->GetParamDV(iDV, 0) == YES) { upper = true;}
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
      
      if (config->GetMarker_All_DV(iMarker) == YES) {
        
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        Normal_ = boundary->vertex[iMarker][iVertex]->GetNormal();
        
        /*--- The CST functions should be applied to a basic airfoil without AoA,
         and unitary chord, a tranformation is required ---*/

        ValCos = cos(AoA*PI_NUMBER/180.0);
        ValSin = sin(AoA*PI_NUMBER/180.0);
        
        Coord[0] = Coord_[0]*ValCos - Coord_[1]*ValSin;
        Coord[0] = max(0.0, Coord[0]); // Coord x should be always positive
        Coord[1] = Coord_[1]*ValCos + Coord_[0]*ValSin;
        
        Normal[0] = Normal_[0]*ValCos - Normal_[1]*ValSin;
        Normal[1] = Normal_[1]*ValCos + Normal_[0]*ValSin;
  
        /*--- CST computation ---*/
        su2double fact_n = 1;
  su2double fact_cst = 1;
        su2double fact_cst_n = 1;
  
  for (int i = 1; i <= maxKulfanNum; i++) {
    fact_n = fact_n * i;
  }
  for (int i = 1; i <= KulfanNum; i++) {
    fact_cst = fact_cst * i;
  }
  for (int i = 1; i <= maxKulfanNum - KulfanNum; i++) {
    fact_cst_n = fact_cst_n * i;
  } 
  
  // CST method only for 2D NACA type airfoils  
  su2double N1, N2;       
  N1 = 0.5;
  N2 = 1.0;
 
  /*--- Upper and lower surface change in coordinates based on CST equations by Kulfan et. al (www.brendakulfan.com/docs/CST3.pdf)  ---*/
        fk = pow(Coord[0],N1)*pow((1-Coord[0]), N2) * fact_n/(fact_cst*(fact_cst_n)) * pow(Coord[0], KulfanNum) * pow((1-Coord[0]), (maxKulfanNum-(KulfanNum)));

  if (( upper) && (Normal[1] > 0)) { VarCoord[1] =  Ampl*fk; }

        if ((!upper) && (Normal[1] < 0)) { VarCoord[1] =  Ampl*fk; }

  
  }
      
      /*--- Apply the transformation to the coordinate variation ---*/

      ValCos = cos(-AoA*PI_NUMBER/180.0);
      ValSin = sin(-AoA*PI_NUMBER/180.0);

      VarCoord_[0] = VarCoord[0]*ValCos - VarCoord[1]*ValSin;
      VarCoord_[1] = VarCoord[1]*ValCos + VarCoord[0]*ValSin;

            boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord_);
    }
  }
}

void CSurfaceMovement::SetRotation(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3] = {0.0,0.0,0.0}, *Coord;
  su2double movement[3] = {0.0,0.0,0.0}, x, y, z;
  su2double Scale = config->GetOpt_RelaxFactor();
  
  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/
  
  if ((iDV == 0) || (ResetDef == true)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
  }
  
  /*--- xyz-coordinates of a point on the line of rotation. */
  
  su2double a = config->GetParamDV(iDV, 0);
  su2double b = config->GetParamDV(iDV, 1);
  su2double c = 0.0;
  if (boundary->GetnDim() == 3) c = config->GetParamDV(0,2);
  
  /*--- xyz-coordinate of the line's direction vector. ---*/
  
  su2double u = config->GetParamDV(iDV, 3)-config->GetParamDV(iDV, 0);
  su2double v = config->GetParamDV(iDV, 4)-config->GetParamDV(iDV, 1);
  su2double w = 1.0;
  if (boundary->GetnDim() == 3) w = config->GetParamDV(iDV, 5)-config->GetParamDV(iDV, 2);
  
  /*--- The angle of rotation. ---*/
  
  su2double theta = config->GetDV_Value(iDV)*Scale*PI_NUMBER/180.0;
  
  /*--- An intermediate value used in computations. ---*/
  
  su2double u2=u*u; su2double v2=v*v; su2double w2=w*w;
  su2double cosT = cos(theta); su2double sinT = sin(theta);
  su2double l2 = u2 + v2 + w2; su2double l = sqrt(l2);
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
      if (config->GetMarker_All_DV(iMarker) == YES) {
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
        x = Coord[0]; y = Coord[1]; z = Coord[2];
        
        movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
        + (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
        + l*(-c*v + b*w - w*y + v*z)*sinT;
        movement[0] = movement[0]/l2 - x;
        
        movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z)
        + (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
        + l*(c*u - a*w + w*x - u*z)*sinT;
        movement[1] = movement[1]/l2 - y;
        
        movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z)
        + (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
        + l*(-b*u + a*v - v*x + u*y)*sinT;
        if (boundary->GetnDim() == 3) movement[2] = movement[2]/l2 - z;
        else movement[2] = 0.0;
        
        VarCoord[0] = movement[0];
        VarCoord[1] = movement[1];
        if (boundary->GetnDim() == 3) VarCoord[2] = movement[2];
        
      }
      boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
    }
}

void CSurfaceMovement::SetTranslation(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3] = {0.0,0.0,0.0};
  su2double Scale = config->GetOpt_RelaxFactor();
  su2double Ampl = config->GetDV_Value(iDV)*Scale;
  
  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/
  
  if ((iDV == 0) || (ResetDef == true)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
  }
  
  su2double xDispl = config->GetParamDV(iDV, 0);
  su2double yDispl = config->GetParamDV(iDV, 1);
  su2double zDispl = 0;
  if (boundary->GetnDim() == 3) zDispl = config->GetParamDV(iDV, 2);
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
      if (config->GetMarker_All_DV(iMarker) == YES) {
        VarCoord[0] = Ampl*xDispl;
        VarCoord[1] = Ampl*yDispl;
        if (boundary->GetnDim() == 3) VarCoord[2] = Ampl*zDispl;
      }
      boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
    }
  
}

void CSurfaceMovement::SetScale(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3] = {0.0,0.0,0.0}, x, y, z, *Coord;
  su2double Scale = config->GetOpt_RelaxFactor();
  su2double Ampl = config->GetDV_Value(iDV)*Scale;
  
  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/
  
  if ((iDV == 0) || (ResetDef == true)) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
        boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
  }
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
      if (config->GetMarker_All_DV(iMarker) == YES) {
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
        x = Coord[0]; y = Coord[1]; z = Coord[2];
        VarCoord[0] = (Ampl-1.0)*x;
        VarCoord[1] = (Ampl-1.0)*y;
        if (boundary->GetnDim() == 3) VarCoord[2] = (Ampl-1.0)*z;
      }
      boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
    }
  
}

void CSurfaceMovement::Moving_Walls(CGeometry *geometry, CConfig *config,
                                    unsigned short iZone, unsigned long iter) {
  
  /*--- Local variables ---*/
  unsigned short iMarker, jMarker, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iVertex;
  su2double xDot[3] = {0.0,0.0,0.0}, *Coord, Center[3] = {0.0,0.0,0.0}, Omega[3] = {0.0,0.0,0.0}, r[3] = {0.0,0.0,0.0}, GridVel[3] = {0.0,0.0,0.0};
  su2double L_Ref     = config->GetLength_Ref();
  su2double Omega_Ref = config->GetOmega_Ref();
  su2double Vel_Ref   = config->GetVelocity_Ref();
  string Marker_Tag;
  
  /*--- Store grid velocity for each node on the moving surface(s).
   Sum and store the x, y, & z velocities due to translation and rotation. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Moving(iMarker) == YES) {
      
      /*--- Identify iMarker from the list of those under MARKER_MOVING ---*/
      
      Marker_Tag = config->GetMarker_All_TagBound(iMarker);
      jMarker    = config->GetMarker_Moving(Marker_Tag);
      
      /*--- Get prescribed wall speed from config for this marker ---*/
      
      for (iDim = 0; iDim < 3; iDim++){
        Center[iDim] = config->GetMarkerMotion_Origin(jMarker, iDim);
        Omega[iDim]  = config->GetMarkerRotationRate(jMarker, iDim)/Omega_Ref;
        xDot[iDim]   = config->GetMarkerTranslationRate(jMarker, iDim)/Vel_Ref;
      }
      
      
      if (rank == MASTER_NODE && iter == 0) {
        cout << " Storing grid velocity for marker: ";
        cout << Marker_Tag << "." << endl;
        cout << " Translational velocity: (" << xDot[0]*config->GetVelocity_Ref() << ", " << xDot[1]*config->GetVelocity_Ref();
        cout << ", " << xDot[2]*config->GetVelocity_Ref();
        if (config->GetSystemMeasurements() == SI) cout << ") m/s." << endl;
        else cout << ") ft/s." << endl;
        cout << " Angular velocity: (" << Omega[0] << ", " << Omega[1];
        cout << ", " << Omega[2] << ") rad/s about origin: (" << Center[0];
        cout << ", " << Center[1] << ", " << Center[2] << ")." << endl;
      }
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        /*--- Get the index and coordinates of the current point ---*/
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Coord  = geometry->node[iPoint]->GetCoord();
        
        /*--- Calculate non-dim. position from rotation center ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          r[iDim] = (Coord[iDim]-Center[iDim])/L_Ref;
        if (nDim == 2) r[nDim] = 0.0;
        
        /*--- Cross Product of angular velocity and distance from center to
         get the rotational velocity. Note that we are adding on the velocity
         due to pure translation as well. ---*/
        
        GridVel[0] = xDot[0] + Omega[1]*r[2] - Omega[2]*r[1];
        GridVel[1] = xDot[1] + Omega[2]*r[0] - Omega[0]*r[2];
        GridVel[2] = xDot[2] + Omega[0]*r[1] - Omega[1]*r[0];
        
        /*--- Store the moving wall velocity for this node ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          geometry->node[iPoint]->SetGridVel(iDim, GridVel[iDim]);
  
      }
    }
  }
}

void CSurfaceMovement::Surface_Translating(CGeometry *geometry, CConfig *config,
                                        unsigned long iter, unsigned short iZone) {
  
  su2double deltaT, time_new, time_old;
  su2double Center[3] = {0.0,0.0,0.0}, VarCoord[3] = {0.0,0.0,0.0};
  su2double xDot[3] = {0.0,0.0,0.0};
  unsigned short iMarker, jMarker, Moving;
  unsigned long iVertex;
  string Marker_Tag, Moving_Tag;
  unsigned short iDim;
  
  /*--- Initialize the delta variation in coordinates ---*/
  VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
  
  /*--- Retrieve values from the config file ---*/
  
  deltaT = config->GetDelta_UnstTimeND();
  
  /*--- Compute delta time based on physical time step ---*/
  time_new = static_cast<su2double>(iter)*deltaT;
  if (iter == 0) {
    time_old = time_new;
  } else {
    time_old = static_cast<su2double>(iter-1)*deltaT;
  }
  
  /*--- Store displacement of each node on the translating surface ---*/
    /*--- Loop over markers and find the particular marker(s) (surface) to translate ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Moving = config->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config->GetnMarker_Moving(); jMarker++) {
        
        Moving_Tag = config->GetMarker_Moving_TagBound(jMarker);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        
        if (Marker_Tag == Moving_Tag && (config->GetKind_SurfaceMovement(jMarker) == DEFORMING)) {

          for (iDim = 0; iDim < 3; iDim++){
            xDot[iDim]   = config->GetMarkerTranslationRate(jMarker, iDim);
            Center[iDim] = config->GetMarkerMotion_Origin(jMarker, iDim);
          }
          
          /*--- Print some information to the console. Be verbose at the first
           iteration only (mostly for debugging purposes). ---*/
          // Note that the MASTER_NODE might not contain all the markers being moved.
          
          if (rank == MASTER_NODE) {
            cout << " Storing translating displacement for marker: ";
            cout << Marker_Tag << "." << endl;
            if (iter == 0) {
              cout << " Translational velocity: (" << xDot[0]*config->GetVelocity_Ref() << ", " << xDot[1]*config->GetVelocity_Ref();
              cout << ", " << xDot[2]*config->GetVelocity_Ref();
              if (config->GetSystemMeasurements() == SI) cout << ") m/s." << endl;
              else cout << ") ft/s." << endl;
            }
          }
          
          /*--- Compute delta change in the position in the x, y, & z directions. ---*/
          
          VarCoord[0] = xDot[0]*(time_new-time_old);
          VarCoord[1] = xDot[1]*(time_new-time_old);
          VarCoord[2] = xDot[2]*(time_new-time_old);
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            /*--- Set node displacement for volume deformation ---*/
            geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
            
          }
        }
      }
    }
  }
  
  /*--- When updating the origins it is assumed that all markers have the
        same translational velocity, because we use the last VarCoord set ---*/
  
  /*--- Set the mesh motion center to the new location after
   incrementing the position with the translation. This new
   location will be used for subsequent mesh motion for the given marker.---*/
  
  for (jMarker=0; jMarker<config->GetnMarker_Moving(); jMarker++) {
    
    /*-- Check if we want to update the motion origin for the given marker ---*/
    
    if (config->GetMoveMotion_Origin(jMarker) == YES) {
      for (iDim = 0; iDim < 3; iDim++){
        Center[iDim] += VarCoord[iDim];
      }
      config->SetMarkerMotion_Origin(Center, jMarker);      
    }
  }
  
  /*--- Set the moment computation center to the new location after
   incrementing the position with the translation. ---*/
  
  for (jMarker=0; jMarker<config->GetnMarker_Monitoring(); jMarker++) {
    Center[0] = config->GetRefOriginMoment_X(jMarker) + VarCoord[0];
    Center[1] = config->GetRefOriginMoment_Y(jMarker) + VarCoord[1];
    Center[2] = config->GetRefOriginMoment_Z(jMarker) + VarCoord[2];
    config->SetRefOriginMoment_X(jMarker, Center[0]);
    config->SetRefOriginMoment_Y(jMarker, Center[1]);
    config->SetRefOriginMoment_Z(jMarker, Center[2]);
  }
}

void CSurfaceMovement::Surface_Plunging(CGeometry *geometry, CConfig *config,
                                           unsigned long iter, unsigned short iZone) {
  
  su2double deltaT, time_new, time_old, Lref;
  su2double Center[3] = {0.0, 0.0, 0.0}, VarCoord[3], Omega[3], Ampl[3];
  su2double DEG2RAD = PI_NUMBER/180.0;
  unsigned short iMarker, jMarker, Moving;
  unsigned long iVertex;
  string Marker_Tag, Moving_Tag;
  unsigned short iDim;
  
  /*--- Initialize the delta variation in coordinates ---*/
  VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
  
  /*--- Retrieve values from the config file ---*/
  
  deltaT = config->GetDelta_UnstTimeND();
  Lref   = config->GetLength_Ref();
  
  /*--- Compute delta time based on physical time step ---*/
  time_new = static_cast<su2double>(iter)*deltaT;
  if (iter == 0) {
    time_old = time_new;
  } else {
    time_old = static_cast<su2double>(iter-1)*deltaT;
  }
  
  /*--- Store displacement of each node on the plunging surface ---*/
    /*--- Loop over markers and find the particular marker(s) (surface) to plunge ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Moving = config->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config->GetnMarker_Moving(); jMarker++) {
        
        Moving_Tag = config->GetMarker_Moving_TagBound(jMarker);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        
        if (Marker_Tag == Moving_Tag && (config->GetKind_SurfaceMovement(jMarker) == DEFORMING)) {
          
          /*--- Plunging frequency and amplitude from config. ---*/
          
          for (iDim = 0; iDim < 3; iDim++){
            Ampl[iDim]   = config->GetMarkerPlunging_Ampl(jMarker, iDim)/Lref;
            Omega[iDim]  = config->GetMarkerPlunging_Omega(jMarker, iDim)/config->GetOmega_Ref();
            Center[iDim] = config->GetMarkerMotion_Origin(jMarker, iDim);
          }
          /*--- Print some information to the console. Be verbose at the first
           iteration only (mostly for debugging purposes). ---*/
          // Note that the MASTER_NODE might not contain all the markers being moved.

          if (rank == MASTER_NODE) {
            cout << " Storing plunging displacement for marker: ";
            cout << Marker_Tag << "." << endl;
            if (iter == 0) {
              cout << " Plunging frequency: (" << Omega[0] << ", " << Omega[1];
              cout << ", " << Omega[2] << ") rad/s." << endl;
              cout << " Plunging amplitude: (" << Ampl[0]/DEG2RAD;
              cout << ", " << Ampl[1]/DEG2RAD << ", " << Ampl[2]/DEG2RAD;
              cout << ") degrees."<< endl;
            }
          }
          
          /*--- Compute delta change in the position in the x, y, & z directions. ---*/
          
          VarCoord[0] = -Ampl[0]*(sin(Omega[0]*time_new) - sin(Omega[0]*time_old));
          VarCoord[1] = -Ampl[1]*(sin(Omega[1]*time_new) - sin(Omega[1]*time_old));
          VarCoord[2] = -Ampl[2]*(sin(Omega[2]*time_new) - sin(Omega[2]*time_old));
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            /*--- Set node displacement for volume deformation ---*/
            geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
            
          }
        }
      }
    }
  }
  
  /*--- When updating the origins it is assumed that all markers have the
   same plunging movement, because we use the last VarCoord set ---*/
  
  /*--- Set the mesh motion center to the new location after
   incrementing the position with the translation. This new
   location will be used for subsequent mesh motion for the given marker.---*/
  
  for (jMarker=0; jMarker<config->GetnMarker_Moving(); jMarker++) {
    
    /*-- Check if we want to update the motion origin for the given marker ---*/
    
    if (config->GetMoveMotion_Origin(jMarker) == YES) {
      for (iDim = 0; iDim < 3; iDim++){
        Center[iDim] += VarCoord[iDim];
      }
      config->SetMarkerMotion_Origin(Center, jMarker);      
    }
  }
  
  /*--- Set the moment computation center to the new location after
   incrementing the position with the plunging. ---*/
  
  for (jMarker=0; jMarker<config->GetnMarker_Monitoring(); jMarker++) {
    Center[0] = config->GetRefOriginMoment_X(jMarker) + VarCoord[0];
    Center[1] = config->GetRefOriginMoment_Y(jMarker) + VarCoord[1];
    Center[2] = config->GetRefOriginMoment_Z(jMarker) + VarCoord[2];
    config->SetRefOriginMoment_X(jMarker, Center[0]);
    config->SetRefOriginMoment_Y(jMarker, Center[1]);
    config->SetRefOriginMoment_Z(jMarker, Center[2]);
  }
}

void CSurfaceMovement::Surface_Pitching(CGeometry *geometry, CConfig *config,
                                        unsigned long iter, unsigned short iZone) {
  
  su2double deltaT, time_new, time_old, Lref, *Coord;
  su2double Center[3], VarCoord[3], Omega[3], Ampl[3], Phase[3];
  su2double rotCoord[3], r[3] = {0.0,0.0,0.0};
  su2double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  su2double dtheta, dphi, dpsi, cosTheta, sinTheta;
  su2double cosPhi, sinPhi, cosPsi, sinPsi;
  su2double DEG2RAD = PI_NUMBER/180.0;
  unsigned short iMarker, jMarker, Moving, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iVertex;
  string Marker_Tag, Moving_Tag;
  
  /*--- Initialize the delta variation in coordinates ---*/
  VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
  
  /*--- Retrieve values from the config file ---*/
  
  deltaT = config->GetDelta_UnstTimeND();
  Lref   = config->GetLength_Ref();
  
  /*--- Compute delta time based on physical time step ---*/
  time_new = static_cast<su2double>(iter)*deltaT;
  if (iter == 0) {
    time_old = time_new;
  } else {
    time_old = static_cast<su2double>(iter-1)*deltaT;
  }

  /*--- Store displacement of each node on the pitching surface ---*/
    /*--- Loop over markers and find the particular marker(s) (surface) to pitch ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Moving = config->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config->GetnMarker_Moving(); jMarker++) {
        
        Moving_Tag = config->GetMarker_Moving_TagBound(jMarker);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        
        if (Marker_Tag == Moving_Tag && (config->GetKind_SurfaceMovement(jMarker) == DEFORMING)) {
          
          /*--- Pitching origin, frequency, and amplitude from config. ---*/
          
          for (iDim = 0; iDim < 3; iDim++){
            Ampl[iDim]   = config->GetMarkerPitching_Ampl(jMarker, iDim)*DEG2RAD;
            Omega[iDim]  = config->GetMarkerPitching_Omega(jMarker, iDim)/config->GetOmega_Ref();
            Phase[iDim]  = config->GetMarkerPitching_Phase(jMarker, iDim)*DEG2RAD;
            Center[iDim] = config->GetMarkerMotion_Origin(jMarker, iDim);
          }
          /*--- Print some information to the console. Be verbose at the first
           iteration only (mostly for debugging purposes). ---*/
          // Note that the MASTER_NODE might not contain all the markers being moved.

          if (rank == MASTER_NODE) {
            cout << " Storing pitching displacement for marker: ";
            cout << Marker_Tag << "." << endl;
            if (iter == 0) {
              cout << " Pitching frequency: (" << Omega[0] << ", " << Omega[1];
              cout << ", " << Omega[2] << ") rad/s about origin: (" << Center[0];
              cout << ", " << Center[1] << ", " << Center[2] << ")." << endl;
              cout << " Pitching amplitude about origin: (" << Ampl[0]/DEG2RAD;
              cout << ", " << Ampl[1]/DEG2RAD << ", " << Ampl[2]/DEG2RAD;
              cout << ") degrees."<< endl;
              cout << " Pitching phase lag about origin: (" << Phase[0]/DEG2RAD;
              cout << ", " << Phase[1]/DEG2RAD <<", "<< Phase[2]/DEG2RAD;
              cout << ") degrees."<< endl;
            }
          }
          
          /*--- Compute delta change in the angle about the x, y, & z axes. ---*/
          
          dtheta = -Ampl[0]*(sin(Omega[0]*time_new + Phase[0])
                             - sin(Omega[0]*time_old + Phase[0]));
          dphi   = -Ampl[1]*(sin(Omega[1]*time_new + Phase[1])
                             - sin(Omega[1]*time_old + Phase[1]));
          dpsi   = -Ampl[2]*(sin(Omega[2]*time_new + Phase[2])
                             - sin(Omega[2]*time_old + Phase[2]));
          
          /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
          
          cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
          sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
          
          /*--- Compute the rotation matrix. Note that the implicit
           ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = cosPhi*sinPsi;
          rotMatrix[2][0] = -sinPhi;
          
          rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = sinTheta*cosPhi;
          
          rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
          rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
          rotMatrix[2][2] = cosTheta*cosPhi;
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            /*--- Index and coordinates of the current point ---*/
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            Coord  = geometry->node[iPoint]->GetCoord();
            
            /*--- Calculate non-dim. position from rotation center ---*/
            
            for (iDim = 0; iDim < nDim; iDim++)
              r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
            if (nDim == 2) r[nDim] = 0.0;
            
            /*--- Compute transformed point coordinates ---*/
            
            rotCoord[0] = rotMatrix[0][0]*r[0]
                        + rotMatrix[0][1]*r[1]
                        + rotMatrix[0][2]*r[2] + Center[0];
            
            rotCoord[1] = rotMatrix[1][0]*r[0]
                        + rotMatrix[1][1]*r[1]
                        + rotMatrix[1][2]*r[2] + Center[1];
            
            rotCoord[2] = rotMatrix[2][0]*r[0]
                        + rotMatrix[2][1]*r[1]
                        + rotMatrix[2][2]*r[2] + Center[2];
            
            /*--- Calculate delta change in the x, y, & z directions ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              VarCoord[iDim] = (rotCoord[iDim]-Coord[iDim])/Lref;
            if (nDim == 2) VarCoord[nDim] = 0.0;
            
            /*--- Set node displacement for volume deformation ---*/
            geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
            
          }
        }
      }
    }
  }
  /*--- For pitching we don't update the motion origin and moment reference origin. ---*/
}

void CSurfaceMovement::Surface_Rotating(CGeometry *geometry, CConfig *config,
                                        unsigned long iter, unsigned short iZone) {
  
  su2double deltaT, time_new, time_old, Lref, *Coord;
  su2double Center[3] = {0.0,0.0,0.0}, VarCoord[3] = {0.0,0.0,0.0}, Omega[3] = {0.0,0.0,0.0},
  rotCoord[3] = {0.0,0.0,0.0}, r[3] = {0.0,0.0,0.0}, Center_Aux[3] = {0.0,0.0,0.0};
  su2double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  su2double dtheta, dphi, dpsi, cosTheta, sinTheta;
  su2double cosPhi, sinPhi, cosPsi, sinPsi;
  unsigned short iMarker, jMarker, Moving, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iVertex;
  string Marker_Tag, Moving_Tag;
  
  /*--- Initialize the delta variation in coordinates ---*/
  VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
  
  /*--- Retrieve values from the config file ---*/
  
  deltaT = config->GetDelta_UnstTimeND();
  Lref   = config->GetLength_Ref();
  
  /*--- Compute delta time based on physical time step ---*/
  time_new = static_cast<su2double>(iter)*deltaT;
  if (iter == 0) {
    time_old = time_new;
  } else {
    time_old = static_cast<su2double>(iter-1)*deltaT;
  }
  
  /*--- Store displacement of each node on the rotating surface ---*/
    /*--- Loop over markers and find the particular marker(s) (surface) to rotate ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Moving = config->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config->GetnMarker_Moving(); jMarker++) {
        
        Moving_Tag = config->GetMarker_Moving_TagBound(jMarker);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        
        if (Marker_Tag == Moving_Tag && (config->GetKind_SurfaceMovement(jMarker) == DEFORMING)) {
          
          /*--- Rotation origin and angular velocity from config. ---*/
          
          for (iDim = 0; iDim < 3; iDim++){
            Omega[iDim]  = config->GetMarkerRotationRate(jMarker, iDim)/config->GetOmega_Ref();
            Center[iDim] = config->GetMarkerMotion_Origin(jMarker, iDim);
          }
          /*--- Print some information to the console. Be verbose at the first
           iteration only (mostly for debugging purposes). ---*/
          // Note that the MASTER_NODE might not contain all the markers being moved.

          if (rank == MASTER_NODE) {
            cout << " Storing rotating displacement for marker: ";
            cout << Marker_Tag << "." << endl;
            if (iter == 0) {
              cout << " Angular velocity: (" << Omega[0] << ", " << Omega[1];
              cout << ", " << Omega[2] << ") rad/s about origin: (" << Center[0];
              cout << ", " << Center[1] << ", " << Center[2] << ")." << endl;
            }
          }
          
          /*--- Compute delta change in the angle about the x, y, & z axes. ---*/
          
          dtheta = Omega[0]*(time_new-time_old);
          dphi   = Omega[1]*(time_new-time_old);
          dpsi   = Omega[2]*(time_new-time_old);
          
          /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
          
          cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
          sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
          
          /*--- Compute the rotation matrix. Note that the implicit
           ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = cosPhi*sinPsi;
          rotMatrix[2][0] = -sinPhi;
          
          rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = sinTheta*cosPhi;
          
          rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
          rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
          rotMatrix[2][2] = cosTheta*cosPhi;
          
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            /*--- Index and coordinates of the current point ---*/
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            Coord  = geometry->node[iPoint]->GetCoord();
            
            /*--- Calculate non-dim. position from rotation center ---*/
            
            for (iDim = 0; iDim < nDim; iDim++)
              r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
            if (nDim == 2) r[nDim] = 0.0;
            
            /*--- Compute transformed point coordinates ---*/
            
            rotCoord[0] = rotMatrix[0][0]*r[0]
            + rotMatrix[0][1]*r[1]
            + rotMatrix[0][2]*r[2] + Center[0];
            
            rotCoord[1] = rotMatrix[1][0]*r[0]
            + rotMatrix[1][1]*r[1]
            + rotMatrix[1][2]*r[2] + Center[1];
            
            rotCoord[2] = rotMatrix[2][0]*r[0]
            + rotMatrix[2][1]*r[1]
            + rotMatrix[2][2]*r[2] + Center[2];
            
            /*--- Calculate delta change in the x, y, & z directions ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              VarCoord[iDim] = (rotCoord[iDim]-Coord[iDim])/Lref;
            if (nDim == 2) VarCoord[nDim] = 0.0;
            
            /*--- Set node displacement for volume deformation ---*/
            geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
            
          }
        }
      }
    }
  }
  
  /*--- When updating the origins it is assumed that all markers have the
   same rotation movement, because we use the last markers rotation matrix and center ---*/
  
  /*--- Set the mesh motion center to the new location after
   incrementing the position with the rotation. This new
   location will be used for subsequent mesh motion for the given marker.---*/
  
  for (jMarker=0; jMarker<config->GetnMarker_Moving(); jMarker++) {
    
    /*-- Check if we want to update the motion origin for the given marker ---*/
    
    if (config->GetMoveMotion_Origin(jMarker) == YES) {
        
      for (iDim = 0; iDim < 3; iDim++){
        Center_Aux[iDim] = config->GetMarkerMotion_Origin(jMarker, iDim);
      }
      
      /*--- Calculate non-dim. position from rotation center ---*/
      
      for (iDim = 0; iDim < nDim; iDim++)
        r[iDim] = (Center_Aux[iDim]-Center[iDim])/Lref;
      if (nDim == 2) r[nDim] = 0.0;
      
      /*--- Compute transformed point coordinates ---*/
      
      rotCoord[0] = rotMatrix[0][0]*r[0]
      + rotMatrix[0][1]*r[1]
      + rotMatrix[0][2]*r[2] + Center[0];
      
      rotCoord[1] = rotMatrix[1][0]*r[0]
      + rotMatrix[1][1]*r[1]
      + rotMatrix[1][2]*r[2] + Center[1];
      
      rotCoord[2] = rotMatrix[2][0]*r[0]
      + rotMatrix[2][1]*r[1]
      + rotMatrix[2][2]*r[2] + Center[2];
      
      /*--- Calculate delta change in the x, y, & z directions ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        VarCoord[iDim] = (rotCoord[iDim]-Center_Aux[iDim])/Lref;
      if (nDim == 2) VarCoord[nDim] = 0.0;
      
      for (iDim = 0; iDim < 3; iDim++){
        Center_Aux[iDim] += VarCoord[iDim];
      }
      config->SetMarkerMotion_Origin(Center_Aux, jMarker);      
    }
  }

  /*--- Set the moment computation center to the new location after
   incrementing the position with the rotation. ---*/
  
  for (jMarker=0; jMarker<config->GetnMarker_Monitoring(); jMarker++) {
      
    Center_Aux[0] = config->GetRefOriginMoment_X(jMarker);
    Center_Aux[1] = config->GetRefOriginMoment_Y(jMarker);
    Center_Aux[2] = config->GetRefOriginMoment_Z(jMarker);

    /*--- Calculate non-dim. position from rotation center ---*/
    
    for (iDim = 0; iDim < nDim; iDim++)
      r[iDim] = (Center_Aux[iDim]-Center[iDim])/Lref;
    if (nDim == 2) r[nDim] = 0.0;
    
    /*--- Compute transformed point coordinates ---*/
    
    rotCoord[0] = rotMatrix[0][0]*r[0]
    + rotMatrix[0][1]*r[1]
    + rotMatrix[0][2]*r[2] + Center[0];
    
    rotCoord[1] = rotMatrix[1][0]*r[0]
    + rotMatrix[1][1]*r[1]
    + rotMatrix[1][2]*r[2] + Center[1];
    
    rotCoord[2] = rotMatrix[2][0]*r[0]
    + rotMatrix[2][1]*r[1]
    + rotMatrix[2][2]*r[2] + Center[2];
    
    /*--- Calculate delta change in the x, y, & z directions ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      VarCoord[iDim] = (rotCoord[iDim]-Center_Aux[iDim])/Lref;
    if (nDim == 2) VarCoord[nDim] = 0.0;
    
    config->SetRefOriginMoment_X(jMarker, Center_Aux[0]+VarCoord[0]);
    config->SetRefOriginMoment_Y(jMarker, Center_Aux[1]+VarCoord[1]);
    config->SetRefOriginMoment_Z(jMarker, Center_Aux[2]+VarCoord[2]);
  }
}

void CSurfaceMovement::AeroelasticDeform(CGeometry *geometry, CConfig *config, unsigned long TimeIter, unsigned short iMarker, unsigned short iMarker_Monitoring, vector<su2double>& displacements) {
  
  /* The sign conventions of these are those of the Typical Section Wing Model, below the signs are corrected */
  su2double dh = -displacements[0];           // relative plunge
  su2double dalpha = -displacements[1];       // relative pitch
  su2double dh_x, dh_y;
  su2double Center[2];
  unsigned short iDim;
  su2double Lref = config->GetLength_Ref();
  su2double *Coord;
  unsigned long iPoint, iVertex;
  su2double x_new, y_new;
  su2double VarCoord[3];
  string Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
  
  /*--- Calculate the plunge displacement for the Typical Section Wing Model taking into account rotation ---*/
  if (config->GetKind_GridMovement() == AEROELASTIC_RIGID_MOTION) {
    su2double Omega, dt, psi;
    dt = config->GetDelta_UnstTimeND();
    Omega  = (config->GetRotation_Rate(3)/config->GetOmega_Ref());
    psi = Omega*(dt*TimeIter);
    
    /*--- Correct for the airfoil starting position (This is hardcoded in here) ---*/
    if (Monitoring_Tag == "Airfoil1") {
      psi = psi + 0.0;
    }
    else if (Monitoring_Tag == "Airfoil2") {
      psi = psi + 2.0/3.0*PI_NUMBER;
    }
    else if (Monitoring_Tag == "Airfoil3") {
      psi = psi + 4.0/3.0*PI_NUMBER;
    }
    else
      cout << "WARNING: There is a marker that we are monitoring that doesn't match the values hardcoded above!" << endl;
    
    dh_x = -dh*sin(psi);
    dh_y = dh*cos(psi);
    
  } else {
    dh_x = 0;
    dh_y = dh;
  }
  
  /*--- Pitching origin from config. ---*/
  
  Center[0] = config->GetRefOriginMoment_X(iMarker_Monitoring);
  Center[1] = config->GetRefOriginMoment_Y(iMarker_Monitoring);
  
  for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
    /*--- Coordinates of the current point ---*/
    Coord = geometry->node[iPoint]->GetCoord();
    
    /*--- Calculate non-dim. position from rotation center ---*/
    su2double r[2] = {0,0};
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
        r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
    
    /*--- Compute delta of transformed point coordinates ---*/
    // The deltas are needed for the FEA grid deformation Method.
    // rotation contribution - previous position + plunging contribution
    x_new = cos(dalpha)*r[0] - sin(dalpha)*r[1] -r[0] + dh_x;
    y_new = sin(dalpha)*r[0] + cos(dalpha)*r[1] -r[1] + dh_y;
    
    VarCoord[0] = x_new;
    VarCoord[1] = y_new;
    VarCoord[2] = 0.0;
    
    /*--- Store new delta node locations for the surface ---*/
    geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
  }
  /*--- Set the elastic axis to the new location after incrementing the position with the plunge ---*/
  config->SetRefOriginMoment_X(iMarker_Monitoring, Center[0]+dh_x);
  config->SetRefOriginMoment_Y(iMarker_Monitoring, Center[1]+dh_y);

  
}

void CSurfaceMovement::SetBoundary_Flutter3D(CGeometry *geometry, CConfig *config,
                                             CFreeFormDefBox **FFDBox, unsigned long iter, unsigned short iZone) {
  
  su2double omega, deltaT;
  su2double alpha, alpha_new, alpha_old;
  su2double time_new, time_old;
  su2double Omega[3], Ampl[3];
  su2double DEG2RAD = PI_NUMBER/180.0;
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());
  unsigned short iDim = 0;
  
  /*--- Retrieve values from the config file ---*/
  
  deltaT = config->GetDelta_UnstTimeND();
  
  /*--- Pitching origin, frequency, and amplitude from config. ---*/
  
  for (iDim = 0; iDim < 3; iDim++){
    Omega[iDim] = config->GetPitching_Omega(iDim)/config->GetOmega_Ref();
    Ampl[iDim] = config->GetPitching_Ampl(iDim)*DEG2RAD;
  }
  
  /*--- Compute delta time based on physical time step ---*/
  
  if (adjoint) {
    
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/
    
    unsigned long nFlowIter  = config->GetnTime_Iter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<su2double>(directIter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<su2double>(directIter)+1.0)*deltaT;
  } else {
    
    /*--- Forward time for the direct problem ---*/
    
    time_new = static_cast<su2double>(iter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<su2double>(iter)-1.0)*deltaT;
  }
  
  /*--- Update the pitching angle at this time step. Flip sign for
   nose-up positive convention. ---*/
  
  omega     = Omega[2];
  alpha_new = Ampl[2]*sin(omega*time_new);
  alpha_old = Ampl[2]*sin(omega*time_old);
  alpha     = (1E-10 + (alpha_new - alpha_old))*(-PI_NUMBER/180.0);
  
  if (rank == MASTER_NODE)
    cout << "New dihedral angle (alpha): " << alpha_new/DEG2RAD << " degrees." << endl;
  
  unsigned short iOrder, jOrder, kOrder;
  short iFFDBox;
  su2double movement[3] = {0.0,0.0,0.0};
  bool *move = new bool [nFFDBox];
  unsigned short *index = new unsigned short[3];
  
  move[0] = true; move[1] = true; move[2] = true;  

  /*--- Change the value of the control point if move is true ---*/
  
  for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
    if (move[iFFDBox])
      for (iOrder = 0; iOrder < FFDBox[iFFDBox]->GetlOrder(); iOrder++)
        for (jOrder = 0; jOrder < FFDBox[iFFDBox]->GetmOrder(); jOrder++)
          for (kOrder = 0; kOrder < FFDBox[iFFDBox]->GetnOrder(); kOrder++) {
            index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
            su2double *coord = FFDBox[iFFDBox]->GetCoordControlPoints(iOrder, jOrder, kOrder);
            movement[0] = 0.0; movement[1] = 0.0; movement[2] = coord[1]*tan(alpha);
            FFDBox[iFFDBox]->SetControlPoints(index, movement);
          }
  
  /*--- Recompute cartesian coordinates using the new control points position ---*/
  
  for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
    SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, false);
  
  delete [] index;
  delete [] move;
  
}

void CSurfaceMovement::SetExternal_Deformation(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter) {
  
  /*--- Local variables ---*/
  
  unsigned short iDim, nDim; 
  unsigned long iPoint = 0, flowIter = 0;
  unsigned long jPoint, GlobalIndex;
  su2double VarCoord[3], *Coord_Old = NULL, *Coord_New = NULL, Center[3] = {0.0,0.0,0.0};
  su2double Lref   = config->GetLength_Ref();
  su2double NewCoord[3] = {0.0,0.0,0.0}, rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  su2double r[3] = {0.0,0.0,0.0}, rotCoord[3] = {0.0,0.0,0.0};
  unsigned long iVertex;
  unsigned short iMarker;
  char buffer[50];
  string DV_Filename, UnstExt, text_line;
  ifstream surface_positions;
  bool unsteady = config->GetTime_Marching();
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());
  
  /*--- Load stuff from config ---*/
  
  nDim = geometry->GetnDim();
  DV_Filename = config->GetDV_Filename();
  
  /*--- Set the extension for the correct unsteady mesh motion file ---*/
  
  if (unsteady) {
    if (adjoint) {
      /*--- For the unsteady adjoint, we integrate backwards through
       physical time, so perform mesh motion in reverse. ---*/
      unsigned long nFlowIter = config->GetnTime_Iter() - 1;
      flowIter  = nFlowIter - iter;
      unsigned short lastindex = DV_Filename.find_last_of(".");
      DV_Filename = DV_Filename.substr(0, lastindex);
      if ((SU2_TYPE::Int(flowIter) >= 0) && (SU2_TYPE::Int(flowIter) < 10)) SPRINTF (buffer, "_0000%d.dat", SU2_TYPE::Int(flowIter));
      if ((SU2_TYPE::Int(flowIter) >= 10) && (SU2_TYPE::Int(flowIter) < 100)) SPRINTF (buffer, "_000%d.dat", SU2_TYPE::Int(flowIter));
      if ((SU2_TYPE::Int(flowIter) >= 100) && (SU2_TYPE::Int(flowIter) < 1000)) SPRINTF (buffer, "_00%d.dat", SU2_TYPE::Int(flowIter));
      if ((SU2_TYPE::Int(flowIter) >= 1000) && (SU2_TYPE::Int(flowIter) < 10000)) SPRINTF (buffer, "_0%d.dat", SU2_TYPE::Int(flowIter));
      if (SU2_TYPE::Int(flowIter) >= 10000) SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(flowIter));
      UnstExt = string(buffer);
      DV_Filename.append(UnstExt);
    } else {
      /*--- Forward time for the direct problem ---*/
      flowIter = iter;
      unsigned short lastindex = DV_Filename.find_last_of(".");
      DV_Filename = DV_Filename.substr(0, lastindex);
      if ((SU2_TYPE::Int(flowIter) >= 0) && (SU2_TYPE::Int(flowIter) < 10)) SPRINTF (buffer, "_0000%d.dat", SU2_TYPE::Int(flowIter));
      if ((SU2_TYPE::Int(flowIter) >= 10) && (SU2_TYPE::Int(flowIter) < 100)) SPRINTF (buffer, "_000%d.dat", SU2_TYPE::Int(flowIter));
      if ((SU2_TYPE::Int(flowIter) >= 100) && (SU2_TYPE::Int(flowIter) < 1000)) SPRINTF (buffer, "_00%d.dat", SU2_TYPE::Int(flowIter));
      if ((SU2_TYPE::Int(flowIter) >= 1000) && (SU2_TYPE::Int(flowIter) < 10000)) SPRINTF (buffer, "_0%d.dat", SU2_TYPE::Int(flowIter));
      if (SU2_TYPE::Int(flowIter) >= 10000) SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(flowIter));
      UnstExt = string(buffer);
      DV_Filename.append(UnstExt);
    }
    
    if (rank == MASTER_NODE)
      cout << "Reading in the arbitrary mesh motion from direct iteration " << flowIter << "." << endl;
  }
  
  /*--- Open the motion file ---*/

  surface_positions.open(DV_Filename.data(), ios::in);
  
  /*--- Throw error if there is no file ---*/
  
  if (surface_positions.fail()) {
    SU2_MPI::Error(string("There is no surface positions file ") + DV_Filename, CURRENT_FUNCTION);
  }
  
  /*--- Read in and store the new mesh node locations ---*/ 
  
  while (getline(surface_positions, text_line)) {
    istringstream point_line(text_line);
    if (nDim == 2) point_line >> iPoint >> NewCoord[0] >> NewCoord[1];
    if (nDim == 3) point_line >> iPoint >> NewCoord[0] >> NewCoord[1] >> NewCoord[2];
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_DV(iMarker) == YES && config->GetKind_SU2() == SU2_DEF) ||
          (config->GetMarker_All_Moving(iMarker) == YES && config->GetKind_SU2() == SU2_CFD)) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex = geometry->node[jPoint]->GetGlobalIndex();
          if (GlobalIndex == iPoint) {
            geometry->vertex[iMarker][iVertex]->SetVarCoord(NewCoord);
            break;
          }
        }
      }
    }
  }
  
  /*--- Close the surface positions file ---*/
  
  surface_positions.close();
  
  /*--- If rotating as well, prepare the rotation matrix ---*/
  
  if (config->GetKind_GridMovement() == EXTERNAL_ROTATION) {
    
    /*--- Variables needed only for rotation ---*/
    
    su2double Omega[3], dt;
    su2double dtheta, dphi, dpsi, cosTheta, sinTheta;
    su2double cosPhi, sinPhi, cosPsi, sinPsi;
    
    /*--- Center of rotation & angular velocity vector from config ---*/
    Center[0] = config->GetMotion_Origin(0);
    Center[1] = config->GetMotion_Origin(1);
    Center[2] = config->GetMotion_Origin(2);
    
    /*--- Angular velocity vector from config ---*/
    
    dt = static_cast<su2double>(iter)*config->GetDelta_UnstTimeND();
    Omega[0]  = config->GetRotation_Rate(0);
    Omega[1]  = config->GetRotation_Rate(1);
    Omega[2]  = config->GetRotation_Rate(2);
    
    /*--- For the unsteady adjoint, use reverse time ---*/
    if (adjoint) {
      /*--- Set the first adjoint mesh position to the final direct one ---*/
      if (iter == 0) dt = ((su2double)config->GetnTime_Iter()-1) * dt;
      /*--- Reverse the rotation direction for the adjoint ---*/
      else dt = -1.0*dt;
    } else {
      /*--- No rotation at all for the first direct solution ---*/
      if (iter == 0) dt = 0;
    }
    
    /*--- Compute delta change in the angle about the x, y, & z axes. ---*/
    
    dtheta = Omega[0]*dt;   
    dphi   = Omega[1]*dt; 
    dpsi   = Omega[2]*dt;
    
    /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
    
    cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
    sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
    
    /*--- Compute the rotation matrix. Note that the implicit
     ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
    
    rotMatrix[0][0] = cosPhi*cosPsi;
    rotMatrix[1][0] = cosPhi*sinPsi;
    rotMatrix[2][0] = -sinPhi;
    
    rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
    rotMatrix[2][1] = sinTheta*cosPhi;
    
    rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
    rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
    rotMatrix[2][2] = cosTheta*cosPhi;
    
  }
  
  /*--- Loop through to find only moving surface markers ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_DV(iMarker) == YES && config->GetKind_SU2() == SU2_DEF) ||
        (config->GetMarker_All_Moving(iMarker) == YES && config->GetKind_SU2() == SU2_CFD)) {
      
      /*--- Loop over all surface points for this marker ---*/
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- Get current and new coordinates from file ---*/
        
        Coord_Old = geometry->node[iPoint]->GetCoord();
        Coord_New = geometry->vertex[iMarker][iVertex]->GetVarCoord();
        
        /*--- If we're also rotating, multiply each point by the
         rotation matrix. It is assumed that the coordinates in
         Coord_Old have already been rotated using SetRigid_Rotation(). ---*/
        
        if (config->GetKind_GridMovement() == EXTERNAL_ROTATION) {
          
          /*--- Calculate non-dim. position from rotation center ---*/
          
          for (iDim = 0; iDim < nDim; iDim++)
            r[iDim] = (Coord_New[iDim]-Center[iDim])/Lref;
          if (nDim == 2) r[nDim] = 0.0;
          
          /*--- Compute transformed point coordinates ---*/
          
          rotCoord[0] = rotMatrix[0][0]*r[0] 
                      + rotMatrix[0][1]*r[1] 
                      + rotMatrix[0][2]*r[2] + Center[0];
          
          rotCoord[1] = rotMatrix[1][0]*r[0] 
                      + rotMatrix[1][1]*r[1] 
                      + rotMatrix[1][2]*r[2] + Center[1];
          
          rotCoord[2] = rotMatrix[2][0]*r[0] 
                      + rotMatrix[2][1]*r[1] 
                      + rotMatrix[2][2]*r[2] + Center[2];
          
          /*--- Copy rotated coords back to original array for consistency ---*/
          for (iDim = 0; iDim < nDim; iDim++)
            Coord_New[iDim] = rotCoord[iDim];
        }
        
        /*--- Calculate delta change in the x, y, & z directions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          VarCoord[iDim] = (Coord_New[iDim]-Coord_Old[iDim])/Lref;
        if (nDim == 2) VarCoord[nDim] = 0.0;

        /*--- Set position changes to be applied by the spring analogy ---*/
        geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
        
      }
    }  
  }
}

void CSurfaceMovement::SetNACA_4Digits(CGeometry *boundary, CConfig *config) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3], *Coord, *Normal, Ycurv, Yesp;

  if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple deformations."; cin.get();  }

  su2double Ya = config->GetParamDV(0,0) / 100.0; /*--- Maximum camber as a fraction of the chord 
          (100 m is the first of the four digits) ---*/
  su2double Xa = config->GetParamDV(0,1) / 10.0; /*--- Location of maximum camber as a fraction of 
          the chord (10 p is the second digit in the NACA xxxx description) ---*/
  su2double t = config->GetParamDV(0,2) / 100.0; /*--- Maximum thickness as a fraction of the
            chord (so 100 t gives the last two digits in 
            the NACA 4-digit denomination) ---*/
    
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
      if (config->GetMarker_All_DV(iMarker) == YES) {
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
        Normal = boundary->vertex[iMarker][iVertex]->GetNormal();
        
        if (Coord[0] < Xa) Ycurv = (2.0*Xa*Coord[0]-pow(Coord[0],2.0))*(Ya/pow(Xa,2.0));
        else Ycurv = ((1.0-2.0*Xa)+2.0*Xa*Coord[0]-pow(Coord[0],2.0))*(Ya/pow((1.0-Xa), 2.0));
        
        Yesp = t*(1.4845*sqrt(Coord[0])-0.6300*Coord[0]-1.7580*pow(Coord[0],2.0)+
              1.4215*pow(Coord[0],3.0)-0.518*pow(Coord[0],4.0));
        
        if (Normal[1] > 0) VarCoord[1] =  (Ycurv + Yesp) - Coord[1];
        if (Normal[1] < 0) VarCoord[1] =  (Ycurv - Yesp) - Coord[1];

      }
      boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
    }
}

void CSurfaceMovement::SetParabolic(CGeometry *boundary, CConfig *config) {
  unsigned long iVertex;
  unsigned short iMarker;
  su2double VarCoord[3], *Coord, *Normal;
  
  if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple deformations."; cin.get();  }
  
  su2double c = config->GetParamDV(0,0); /*--- Center of the parabola ---*/
  su2double t = config->GetParamDV(0,1) / 100.0; /*--- Thickness of the parabola ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
      if (config->GetMarker_All_DV(iMarker) == YES) {
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
        Normal = boundary->vertex[iMarker][iVertex]->GetNormal();
        
        if (Normal[1] > 0) {
          VarCoord[1] =  t*(Coord[0]*Coord[0]-Coord[0])/(2.0*(c*c-c)) - Coord[1];
        }
        if (Normal[1] < 0) {
          VarCoord[1] =  t*(Coord[0]-Coord[0]*Coord[0])/(2.0*(c*c-c)) - Coord[1];
        }
      }
      boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
    }
}

void CSurfaceMovement::SetAirfoil(CGeometry *boundary, CConfig *config) {
  unsigned long iVertex, n_Airfoil = 0;
  unsigned short iMarker, nUpper, nLower, iUpper, iLower, iVar, iDim;
  su2double *VarCoord, *Coord, NewYCoord, NewXCoord, *Coord_i, *Coord_ip1, yp1, ypn,
  Airfoil_Coord[2]= {0.0,0.0}, factor, coeff = 10000, Upper, Lower, Arch = 0.0, TotalArch = 0.0,
  x_i, x_ip1, y_i, y_ip1;
  passivedouble AirfoilScale;
  vector<su2double> Svalue, Xcoord, Ycoord, Xcoord2, Ycoord2, Xcoord_Aux, Ycoord_Aux;
  bool AddBegin = true, AddEnd = true;
  char AirfoilFile[256], AirfoilFormat[15], MeshOrientation[15], AirfoilClose[15];
  ifstream airfoil_file;
  string text_line;
  int ierr = 0;

  unsigned short nDim = boundary->GetnDim();
  
  VarCoord = new su2double[nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    VarCoord[iDim] = 0.0;

  /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
   deforming meshes (MARKER_MOVING), while SU2_DEF will use it for deforming
   meshes after imposing design variable surface deformations (DV_MARKER). ---*/
  
  unsigned short Kind_SU2 = config->GetKind_SU2();
  
  /*--- Read the coordinates. Two main formats:
   - Selig are in an x, y format starting from trailing edge, along the upper surface to the leading
   edge and back around the lower surface to trailing edge.
   - Lednicer are upper surface points leading edge to trailing edge and then lower surface leading
   edge to trailing edge.
   ---*/
  
  /*--- Open the restart file, throw an error if this fails. ---*/
  
  cout << "Enter the name of file with the airfoil information: ";
  ierr = scanf("%255s", AirfoilFile);
  if (ierr == 0) { SU2_MPI::Error("No input read!!", CURRENT_FUNCTION); }
  airfoil_file.open(AirfoilFile, ios::in);
  if (airfoil_file.fail()) {
    SU2_MPI::Error(string("There is no airfoil file ") + string(AirfoilFile), CURRENT_FUNCTION);
  }
  cout << "Enter the format of the airfoil (Selig or Lednicer): ";
  ierr = scanf("%14s", AirfoilFormat);
  if (ierr == 0) { SU2_MPI::Error("No input read!!", CURRENT_FUNCTION); }

  cout << "Thickness scaling (1.0 means no scaling)?: ";
  ierr = scanf("%lf", &AirfoilScale);
  if (ierr == 0) { SU2_MPI::Error("No input read!!", CURRENT_FUNCTION); }

  cout << "Close the airfoil (Yes or No)?: ";
  ierr = scanf("%14s", AirfoilClose);
  if (ierr == 0) { SU2_MPI::Error("No input read!!", CURRENT_FUNCTION); }

  cout << "Surface mesh orientation (clockwise, or anticlockwise): ";
  ierr = scanf("%14s", MeshOrientation);
  if (ierr == 0) { SU2_MPI::Error("No input read!!", CURRENT_FUNCTION); }

  /*--- The first line is the header ---*/
  
  getline (airfoil_file, text_line);
  cout << "File info: " << text_line << endl;
  
  if (strcmp (AirfoilFormat,"Selig") == 0) {

    while (getline (airfoil_file, text_line)) {
      istringstream point_line(text_line);
      
      /*--- Read the x & y coordinates from this line of the file (anticlockwise) ---*/
      
      point_line >> Airfoil_Coord[0] >> Airfoil_Coord[1];
      
      /*--- Close the arifoil ---*/
      
      if (strcmp (AirfoilClose,"Yes") == 0)
        factor = -atan(coeff*(Airfoil_Coord[0]-1.0))*2.0/PI_NUMBER;
      else factor = 1.0;
      
      /*--- Store the coordinates in vectors ---*/
      
      Xcoord.push_back(Airfoil_Coord[0]);
      Ycoord.push_back(Airfoil_Coord[1]*factor*AirfoilScale);
    }
    
  }
  if (strcmp (AirfoilFormat,"Lednicer") == 0) {
    
    /*--- The second line is the number of points ---*/

    getline(airfoil_file, text_line);
    istringstream point_line(text_line);
    point_line >> Upper >> Lower;
    
    nUpper = SU2_TYPE::Int(Upper);
    nLower = SU2_TYPE::Int(Lower);
  
    Xcoord.resize(nUpper+nLower-1);
    Ycoord.resize(nUpper+nLower-1);
    
    /*--- White line ---*/

    getline (airfoil_file, text_line);

    for (iUpper = 0; iUpper < nUpper; iUpper++) {
      getline (airfoil_file, text_line);
      istringstream point_line(text_line);
      point_line >> Airfoil_Coord[0] >> Airfoil_Coord[1];
      Xcoord[nUpper-iUpper-1] = Airfoil_Coord[0];
      
      if (strcmp (AirfoilClose,"Yes") == 0)
        factor = -atan(coeff*(Airfoil_Coord[0]-1.0))*2.0/PI_NUMBER;
      else factor = 1.0;
      
      Ycoord[nUpper-iUpper-1] = Airfoil_Coord[1]*AirfoilScale*factor;
    }
    
    getline (airfoil_file, text_line);

    for (iLower = 0; iLower < nLower; iLower++) {
      getline (airfoil_file, text_line);
      istringstream point_line(text_line);
      point_line >> Airfoil_Coord[0] >> Airfoil_Coord[1];
      
      if (strcmp (AirfoilClose,"Yes") == 0)
        factor = -atan(coeff*(Airfoil_Coord[0]-1.0))*2.0/PI_NUMBER;
      else factor = 1.0;
      
      Xcoord[nUpper+iLower-1] = Airfoil_Coord[0];
      Ycoord[nUpper+iLower-1] = Airfoil_Coord[1]*AirfoilScale*factor;
    }
      
  }
  
  /*--- Check the coordinate (1,0) at the beginning and end of the file ---*/
  
  if (Xcoord[0] == 1.0) AddBegin = false;
  if (Xcoord[Xcoord.size()-1] == 1.0) AddEnd = false;
  
  if (AddBegin) { Xcoord.insert(Xcoord.begin(), 1.0);   Ycoord.insert(Ycoord.begin(), 0.0);}
  if (AddEnd) { Xcoord.push_back(1.0);                Ycoord.push_back(0.0);}
  
  /*--- Change the orientation (depend on the input file, and the mesh file) ---*/
  
  if (strcmp (MeshOrientation,"clockwise") == 0) {
    for (iVar = 0; iVar < Xcoord.size(); iVar++) {
      Xcoord_Aux.push_back(Xcoord[iVar]);
      Ycoord_Aux.push_back(Ycoord[iVar]);
    }
    
    for (iVar = 0; iVar < Xcoord.size(); iVar++) {
      Xcoord[iVar] = Xcoord_Aux[Xcoord.size()-iVar-1];
      Ycoord[iVar] = Ycoord_Aux[Xcoord.size()-iVar-1];
    }
  }
  
  /*--- Compute the total arch length ---*/
  
  Arch = 0.0; Svalue.push_back(Arch);

  for (iVar = 0; iVar < Xcoord.size()-1; iVar++) {
    x_i = Xcoord[iVar];  x_ip1 = Xcoord[iVar+1];
    y_i = Ycoord[iVar];  y_ip1 = Ycoord[iVar+1];
    Arch += sqrt((x_ip1-x_i)*(x_ip1-x_i)+(y_ip1-y_i)*(y_ip1-y_i));
    Svalue.push_back(Arch);
  }
  x_i = Xcoord[Xcoord.size()-1];  x_ip1 = Xcoord[0];
  y_i = Ycoord[Xcoord.size()-1];  y_ip1 = Ycoord[0];
  Arch += sqrt((x_ip1-x_i)*(x_ip1-x_i)+(y_ip1-y_i)*(y_ip1-y_i));
  
  /*--- Non dimensionalization ---*/
  
  for (iVar = 0; iVar < Svalue.size(); iVar++) { Svalue[iVar] /= Arch; }

  /*--- Close the restart file ---*/
  
  airfoil_file.close();
  
  /*--- Create a spline for X and Y coordiantes using the arch length ---*/
  
  n_Airfoil = Svalue.size();
  yp1 = (Xcoord[1]-Xcoord[0])/(Svalue[1]-Svalue[0]);
  ypn = (Xcoord[n_Airfoil-1]-Xcoord[n_Airfoil-2])/(Svalue[n_Airfoil-1]-Svalue[n_Airfoil-2]);
  
  Xcoord2.resize(n_Airfoil+1);
  boundary->SetSpline(Svalue, Xcoord, n_Airfoil, yp1, ypn, Xcoord2);
  
  n_Airfoil = Svalue.size();
  yp1 = (Ycoord[1]-Ycoord[0])/(Svalue[1]-Svalue[0]);
  ypn = (Ycoord[n_Airfoil-1]-Ycoord[n_Airfoil-2])/(Svalue[n_Airfoil-1]-Svalue[n_Airfoil-2]);
  
  Ycoord2.resize(n_Airfoil+1);
  boundary->SetSpline(Svalue, Ycoord, n_Airfoil, yp1, ypn, Ycoord2);
  
  TotalArch = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_CFD)) ||
        ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_DEF))) {
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]-1; iVertex++) {
        Coord_i = boundary->vertex[iMarker][iVertex]->GetCoord();
        Coord_ip1 = boundary->vertex[iMarker][iVertex+1]->GetCoord();
        
        x_i = Coord_i[0]; x_ip1 = Coord_ip1[0];
        y_i = Coord_i[1]; y_ip1 = Coord_ip1[1];
        
        TotalArch += sqrt((x_ip1-x_i)*(x_ip1-x_i)+(y_ip1-y_i)*(y_ip1-y_i));
      }
      Coord_i = boundary->vertex[iMarker][boundary->nVertex[iMarker]-1]->GetCoord();
      Coord_ip1 = boundary->vertex[iMarker][0]->GetCoord();
      x_i = Coord_i[0]; x_ip1 = Coord_ip1[0];
      y_i = Coord_i[1]; y_ip1 = Coord_ip1[1];
      TotalArch += sqrt((x_ip1-x_i)*(x_ip1-x_i)+(y_ip1-y_i)*(y_ip1-y_i));
    }
  }
  
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Arch = 0.0;
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
      if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_CFD)) ||
          ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_DEF))) {
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
        
        if (iVertex == 0) Arch = 0.0;
        else {
          Coord_i = boundary->vertex[iMarker][iVertex-1]->GetCoord();
          Coord_ip1 = boundary->vertex[iMarker][iVertex]->GetCoord();
          x_i = Coord_i[0]; x_ip1 = Coord_ip1[0];
          y_i = Coord_i[1]; y_ip1 = Coord_ip1[1];
          Arch += sqrt((x_ip1-x_i)*(x_ip1-x_i)+(y_ip1-y_i)*(y_ip1-y_i))/TotalArch;
        }
        
        NewXCoord = boundary->GetSpline(Svalue, Xcoord, Xcoord2, n_Airfoil, Arch);
        NewYCoord = boundary->GetSpline(Svalue, Ycoord, Ycoord2, n_Airfoil, Arch);
        
        /*--- Store the delta change in the x & y coordinates ---*/
        
        VarCoord[0] = NewXCoord - Coord[0];
        VarCoord[1] = NewYCoord - Coord[1];
      }

      boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      
    }
  }

  delete [] VarCoord;
  
}

void CSurfaceMovement::ReadFFDInfo(CGeometry *geometry, CConfig *config, CFreeFormDefBox **FFDBox, string val_mesh_filename) {
  
  string text_line, iTag;
  ifstream mesh_file;
  su2double CPcoord[3], coord[] = {0,0,0};
  unsigned short degree[3], iFFDBox, iCornerPoints, iControlPoints, iMarker, iDegree, jDegree, kDegree,
  iChar, LevelFFDBox, nParentFFDBox, iParentFFDBox, nChildFFDBox, iChildFFDBox, nMarker, *nCornerPoints,
  *nControlPoints;
  unsigned long iSurfacePoints, iPoint, jPoint, iVertex, nVertex, nPoint, iElem = 0,
  nElem, my_nSurfPoints, nSurfPoints, *nSurfacePoints;
  su2double XCoord, YCoord;

  bool polar = (config->GetFFD_CoordSystem() == POLAR);
  unsigned short nDim = geometry->GetnDim(), iDim;
  unsigned short SplineOrder[3];
  unsigned short Blending = 0;

  char *cstr = new char [val_mesh_filename.size()+1];
  strcpy (cstr, val_mesh_filename.c_str());
  
  mesh_file.open(cstr, ios::in);
  if (mesh_file.fail()) {
    SU2_MPI::Error("There is no geometry file (ReadFFDInfo)!!", CURRENT_FUNCTION);
  }
  
  while (getline (mesh_file, text_line)) {
    
    /*--- Read the inner elements ---*/
    
    string::size_type position = text_line.find ("NELEM=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nElem = atoi(text_line.c_str());
      for (iElem = 0; iElem < nElem; iElem++) {
        getline(mesh_file, text_line);
      }
    }
    
    /*--- Read the inner points ---*/
    
    position = text_line.find ("NPOIN=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nPoint = atoi(text_line.c_str());
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        getline(mesh_file, text_line);
      }
    }

    /*--- Read the boundaries  ---*/
    
    position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nMarker = atoi(text_line.c_str());
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        getline(mesh_file, text_line);
        getline(mesh_file, text_line);
        text_line.erase (0,13); nVertex = atoi(text_line.c_str());
        for (iVertex = 0; iVertex < nVertex; iVertex++) {
          getline(mesh_file, text_line);
        }
      }
    }
    
    /*--- Read the FFDBox information  ---*/
    
    position = text_line.find ("FFD_NBOX=",0);
    if (position != string::npos) {
      text_line.erase (0,9);
      nFFDBox = atoi(text_line.c_str());
      
      if (rank == MASTER_NODE) cout << nFFDBox << " Free Form Deformation boxes." << endl;
      
      nCornerPoints = new unsigned short[nFFDBox];
      nControlPoints = new unsigned short[nFFDBox];
      nSurfacePoints = new unsigned long[nFFDBox];
      
      getline (mesh_file, text_line);
      text_line.erase (0,11);
      nLevel = atoi(text_line.c_str());
      
      if (rank == MASTER_NODE) cout << nLevel << " Free Form Deformation nested levels." << endl;

      for (iFFDBox = 0 ; iFFDBox < nFFDBox; iFFDBox++) {
        
        /*--- Read the name of the FFD box ---*/
        
        getline (mesh_file, text_line);
        text_line.erase (0,8);
        
        /*--- Remove extra data from the FFDBox name ---*/
        
        string::size_type position;
        for (iChar = 0; iChar < 20; iChar++) {
          position = text_line.find( " ", 0 );
          if (position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\r", 0 );
          if (position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\n", 0 );
          if (position != string::npos) text_line.erase (position,1);
        }
        
        string TagFFDBox = text_line.c_str();
        
        if (rank == MASTER_NODE) cout << "FFD box tag: " << TagFFDBox <<". ";

        /*--- Read the level of the FFD box ---*/
        
        getline (mesh_file, text_line);
        text_line.erase (0,10);
        LevelFFDBox = atoi(text_line.c_str());
        
        if (rank == MASTER_NODE) cout << "FFD box level: " << LevelFFDBox <<". ";
        
        /*--- Read the degree of the FFD box ---*/
        
        
        if (nDim == 2) {
          if (polar) {
            getline (mesh_file, text_line);
            text_line.erase (0,13); degree[0] = atoi(text_line.c_str());
            degree[1] = 1;
            getline (mesh_file, text_line);
            text_line.erase (0,13); degree[2] = atoi(text_line.c_str());
          }
          else {
            getline (mesh_file, text_line);
            text_line.erase (0,13); degree[0] = atoi(text_line.c_str());
            getline (mesh_file, text_line);
            text_line.erase (0,13); degree[1] = atoi(text_line.c_str());
            degree[2] = 1;
          }
        }
        else {
          getline (mesh_file, text_line);
          text_line.erase (0,13); degree[0] = atoi(text_line.c_str());
          getline (mesh_file, text_line);
          text_line.erase (0,13); degree[1] = atoi(text_line.c_str());
          getline (mesh_file, text_line);
          text_line.erase (0,13); degree[2] = atoi(text_line.c_str());
        }
        
        if (rank == MASTER_NODE) {
          if (nDim == 2) {
            if (polar) cout << "Degrees: " << degree[0] << ", " << degree[2] << "." << endl;
            else cout << "Degrees: " << degree[0] << ", " << degree[1] << "." << endl;
          }
          else cout << "Degrees: " << degree[0] << ", " << degree[1] << ", " << degree[2] << "." << endl;
        }

        getline (mesh_file, text_line);
        if (text_line.substr(0,12) != "FFD_BLENDING"){
          SU2_MPI::Error(string("Deprecated FFD information found in mesh file.\n") +
                         string("FFD information generated with SU2 version <= 4.3 is incompatible with the current version.") +
                         string("Run SU2_DEF again with DV_KIND= FFD_SETTING."), CURRENT_FUNCTION);
        }
        text_line.erase(0,14);
        if (text_line == "BEZIER"){
          Blending = BEZIER;
        }
        if (text_line == "BSPLINE_UNIFORM"){
          Blending = BSPLINE_UNIFORM;
        }

        if (Blending == BSPLINE_UNIFORM) {
          getline (mesh_file, text_line);
          text_line.erase (0,17); SplineOrder[0] = atoi(text_line.c_str());
          getline (mesh_file, text_line);
          text_line.erase (0,17); SplineOrder[1] = atoi(text_line.c_str());
          if (nDim == 3){
            getline (mesh_file, text_line);
            text_line.erase (0,17); SplineOrder[2] = atoi(text_line.c_str());
          } else {
            SplineOrder[2] = 2;
          }
        }
        if (rank == MASTER_NODE){
          if (Blending == BSPLINE_UNIFORM){
            cout << "FFD Blending using B-Splines. ";
            cout << "Order: " << SplineOrder[0] << ", " << SplineOrder[1];
            if (nDim == 3) cout << ", " << SplineOrder[2];
            cout << ". " << endl;
          }
          if (Blending == BEZIER){
            cout << "FFD Blending using Bezier Curves." << endl;
          }
        }

        FFDBox[iFFDBox] = new CFreeFormDefBox(degree, SplineOrder, Blending);
        FFDBox[iFFDBox]->SetTag(TagFFDBox); FFDBox[iFFDBox]->SetLevel(LevelFFDBox);

        /*--- Read the number of parents boxes ---*/
        
        getline (mesh_file, text_line);
        text_line.erase (0,12);
        nParentFFDBox = atoi(text_line.c_str());
        if (rank == MASTER_NODE) cout << "Number of parent boxes: " << nParentFFDBox <<". ";
        for (iParentFFDBox = 0; iParentFFDBox < nParentFFDBox; iParentFFDBox++) {
          getline(mesh_file, text_line);
          
          /*--- Remove extra data from the FFDBox name ---*/
          
          string::size_type position;
          for (iChar = 0; iChar < 20; iChar++) {
            position = text_line.find( " ", 0 );
            if (position != string::npos) text_line.erase (position,1);
            position = text_line.find( "\r", 0 );
            if (position != string::npos) text_line.erase (position,1);
            position = text_line.find( "\n", 0 );
            if (position != string::npos) text_line.erase (position,1);
          }
          
          string ParentFFDBox = text_line.c_str();
          FFDBox[iFFDBox]->SetParentFFDBox(ParentFFDBox);
        }
        
        /*--- Read the number of children boxes ---*/
        
        getline (mesh_file, text_line);
        text_line.erase (0,13);
        nChildFFDBox = atoi(text_line.c_str());
        if (rank == MASTER_NODE) cout << "Number of child boxes: " << nChildFFDBox <<"." << endl;
        
        for (iChildFFDBox = 0; iChildFFDBox < nChildFFDBox; iChildFFDBox++) {
          getline(mesh_file, text_line);
          
          /*--- Remove extra data from the FFDBox name ---*/
          
          string::size_type position;
          for (iChar = 0; iChar < 20; iChar++) {
            position = text_line.find( " ", 0 );
            if (position != string::npos) text_line.erase (position,1);
            position = text_line.find( "\r", 0 );
            if (position != string::npos) text_line.erase (position,1);
            position = text_line.find( "\n", 0 );
            if (position != string::npos) text_line.erase (position,1);
          }
          
          string ChildFFDBox = text_line.c_str();
          FFDBox[iFFDBox]->SetChildFFDBox(ChildFFDBox);
        }

        /*--- Read the number of the corner points ---*/
        
        getline (mesh_file, text_line);
        text_line.erase (0,18); nCornerPoints[iFFDBox] = atoi(text_line.c_str());
        if (rank == MASTER_NODE) cout << "Corner points: " << nCornerPoints[iFFDBox] <<". ";
        if (nDim == 2) nCornerPoints[iFFDBox] = nCornerPoints[iFFDBox]*SU2_TYPE::Int(2);



        /*--- Read the coordinates of the corner points ---*/
        
          
        if (nDim == 2) {

          if (polar) {

            getline(mesh_file, text_line); istringstream FFDBox_line_1(text_line);
            FFDBox_line_1 >> XCoord; FFDBox_line_1 >> YCoord;

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1)*YCoord;
            CPcoord[2] = -sin(0.1)*YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 4);

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1)*YCoord;
            CPcoord[2] = sin(0.1)*YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 7);

            getline(mesh_file, text_line); istringstream FFDBox_line_2(text_line);
            FFDBox_line_2 >> XCoord; FFDBox_line_2 >> YCoord;

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1)*YCoord;
            CPcoord[2] = -sin(0.1)*YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 0);

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1)*YCoord;
            CPcoord[2] = sin(0.1)*YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 3);

            getline(mesh_file, text_line); istringstream FFDBox_line_3(text_line);
            FFDBox_line_3 >> XCoord; FFDBox_line_3 >> YCoord;

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1)*YCoord;
            CPcoord[2] = -sin(0.1)*YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 1);

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1)*YCoord;
            CPcoord[2] = sin(0.1)*YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 2);

            getline(mesh_file, text_line); istringstream FFDBox_line_4(text_line);
            FFDBox_line_4 >> XCoord; FFDBox_line_4 >> YCoord;

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1)*YCoord;
            CPcoord[2] = -sin(0.1)*YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 5);

            CPcoord[0] = XCoord;
            CPcoord[1] = cos(0.1)*YCoord;
            CPcoord[2] = sin(0.1)*YCoord;
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, 6);

          }
          else {
            for (iCornerPoints = 0; iCornerPoints < nCornerPoints[iFFDBox]; iCornerPoints++) {
              if (iCornerPoints < nCornerPoints[iFFDBox]/SU2_TYPE::Int(2)) {
                getline(mesh_file, text_line); istringstream FFDBox_line(text_line);
                FFDBox_line >> CPcoord[0]; FFDBox_line >> CPcoord[1]; CPcoord[2] = -0.5;
              }
              else {
                CPcoord[0] = FFDBox[iFFDBox]->GetCoordCornerPoints(0, iCornerPoints-nCornerPoints[iFFDBox]/SU2_TYPE::Int(2));
                CPcoord[1] = FFDBox[iFFDBox]->GetCoordCornerPoints(1, iCornerPoints-nCornerPoints[iFFDBox]/SU2_TYPE::Int(2));
                CPcoord[2] = 0.5;
              }
              FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, iCornerPoints);
            }
          }

        }
        else {
          for (iCornerPoints = 0; iCornerPoints < nCornerPoints[iFFDBox]; iCornerPoints++) {
            getline(mesh_file, text_line); istringstream FFDBox_line(text_line);
            FFDBox_line >> CPcoord[0]; FFDBox_line >> CPcoord[1]; FFDBox_line >> CPcoord[2];
            FFDBox[iFFDBox]->SetCoordCornerPoints(CPcoord, iCornerPoints);
          }
        }

        /*--- Read the number of the control points ---*/
        
        getline (mesh_file, text_line);
        text_line.erase (0,19); nControlPoints[iFFDBox] = atoi(text_line.c_str());
        
        if (rank == MASTER_NODE) cout << "Control points: " << nControlPoints[iFFDBox] <<". ";
        
        /*--- Method to identify if there is a FFDBox definition ---*/
        
        if (nControlPoints[iFFDBox] != 0) FFDBoxDefinition = true;

        /*--- Read the coordinates of the control points ---*/
        
        for (iControlPoints = 0; iControlPoints < nControlPoints[iFFDBox]; iControlPoints++) {
          getline(mesh_file, text_line); istringstream FFDBox_line(text_line);
          FFDBox_line >> iDegree; FFDBox_line >> jDegree; FFDBox_line >> kDegree; 
          FFDBox_line >> CPcoord[0]; FFDBox_line >> CPcoord[1]; FFDBox_line >> CPcoord[2];
          FFDBox[iFFDBox]->SetCoordControlPoints(CPcoord, iDegree, jDegree, kDegree);
          FFDBox[iFFDBox]->SetCoordControlPoints_Copy(CPcoord, iDegree, jDegree, kDegree);
        }
        
        getline (mesh_file, text_line);
        text_line.erase (0,19); nSurfacePoints[iFFDBox] = atoi(text_line.c_str());

        /*--- The surface points parametric coordinates, all the nodes read the FFD 
         information but they only store their part ---*/
        
        my_nSurfPoints = 0;
        for (iSurfacePoints = 0; iSurfacePoints < nSurfacePoints[iFFDBox]; iSurfacePoints++) {
          getline(mesh_file, text_line); istringstream FFDBox_line(text_line);
          FFDBox_line >> iTag; FFDBox_line >> iPoint;
          
          if (config->GetMarker_All_TagBound(iTag) != -1) {

            iMarker = config->GetMarker_All_TagBound(iTag);
            FFDBox_line >> CPcoord[0]; FFDBox_line >> CPcoord[1]; FFDBox_line >> CPcoord[2];
            
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
              jPoint =  geometry->vertex[iMarker][iVertex]->GetNode();
              if (iPoint == geometry->node[jPoint]->GetGlobalIndex()) {
                for (iDim = 0; iDim < nDim; iDim++) {
                  coord[iDim] = geometry->node[jPoint]->GetCoord()[iDim];
                }
                FFDBox[iFFDBox]->Set_MarkerIndex(iMarker);
                FFDBox[iFFDBox]->Set_VertexIndex(iVertex);
                FFDBox[iFFDBox]->Set_PointIndex(jPoint);
                FFDBox[iFFDBox]->Set_ParametricCoord(CPcoord);
                FFDBox[iFFDBox]->Set_CartesianCoord(coord);
                my_nSurfPoints++;
              }
            }

          }
          
        }
        
        nSurfacePoints[iFFDBox] = my_nSurfPoints;
        
#ifdef HAVE_MPI
        nSurfPoints = 0;
        SU2_MPI::Allreduce(&my_nSurfPoints, &nSurfPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (rank == MASTER_NODE) cout << "Surface points: " << nSurfPoints <<"."<< endl;
#else
        nSurfPoints = my_nSurfPoints;
        if (rank == MASTER_NODE) cout << "Surface points: " << nSurfPoints <<"."<< endl;
#endif
        
      }
      
      delete [] nCornerPoints;
      delete [] nControlPoints;
      delete [] nSurfacePoints;    
    }
  }
  mesh_file.close();
  
  if (nFFDBox == 0) {
    if (rank == MASTER_NODE) cout <<"There is no FFD box definition. Just in case, check the .su2 file" << endl;
  }

}

void CSurfaceMovement::ReadFFDInfo(CGeometry *geometry, CConfig *config, CFreeFormDefBox **FFDBox) {
  
  string text_line, iTag;
  ifstream mesh_file;
  su2double coord[3];
  unsigned short degree[3], iFFDBox, iCornerPoints, LevelFFDBox, nParentFFDBox,
  iParentFFDBox, nChildFFDBox, iChildFFDBox, *nCornerPoints;

  bool polar = (config->GetFFD_CoordSystem() == POLAR);
  unsigned short nDim = geometry->GetnDim(), iDim;
  unsigned short SplineOrder[3]={2,2,2};

  for (iDim = 0; iDim < 3; iDim++){
    SplineOrder[iDim] = SU2_TYPE::Short(config->GetFFD_BSplineOrder()[iDim]);
  }
  
  
  /*--- Read the FFDBox information from the config file ---*/
  
  nFFDBox = config->GetnFFDBox();
  
  if (rank == MASTER_NODE) cout << nFFDBox << " Free Form Deformation boxes." << endl;
  
  nCornerPoints = new unsigned short[nFFDBox];
  
  nLevel = 1; // Nested FFD is not active
  
  if (rank == MASTER_NODE) cout << nLevel << " Free Form Deformation nested levels." << endl;
  
  for (iFFDBox = 0 ; iFFDBox < nFFDBox; iFFDBox++) {
    
    /*--- Read the name of the FFD box ---*/
    
    string TagFFDBox = config->GetTagFFDBox(iFFDBox);
    
    if (rank == MASTER_NODE) cout << "FFD box tag: " << TagFFDBox <<". ";
    
    /*--- Read the level of the FFD box ---*/
    
    LevelFFDBox = 0; // Nested FFD is not active
    
    if (rank == MASTER_NODE) cout << "FFD box level: " << LevelFFDBox <<". ";
    
    /*--- Read the degree of the FFD box ---*/
    
    if (nDim == 2) {
      if (polar) {
        degree[0] = config->GetDegreeFFDBox(iFFDBox, 0);
        degree[1] = 1;
        degree[2] = config->GetDegreeFFDBox(iFFDBox, 1);
      }
      else {
        degree[0] = config->GetDegreeFFDBox(iFFDBox, 0);
        degree[1] = config->GetDegreeFFDBox(iFFDBox, 1);
        degree[2] = 1;
      }
    }
    else {
      degree[0] = config->GetDegreeFFDBox(iFFDBox, 0);
      degree[1] = config->GetDegreeFFDBox(iFFDBox, 1);
      degree[2] = config->GetDegreeFFDBox(iFFDBox, 2);
    }

    if (rank == MASTER_NODE) {
      if (nDim == 2) {
        if (polar) cout << "Degrees: " << degree[0] << ", " << degree[2] << "." << endl;
        else cout << "Degrees: " << degree[0] << ", " << degree[1] << "." << endl;
      }
      else cout << "Degrees: " << degree[0] << ", " << degree[1] << ", " << degree[2] << "." << endl;
    }
    
    if (rank == MASTER_NODE){
      if (config->GetFFD_Blending() == BSPLINE_UNIFORM){
        cout << "FFD Blending using B-Splines. ";
        cout << "Order: " << SplineOrder[0] << ", " << SplineOrder[1];
        if (nDim == 3) cout << ", " << SplineOrder[2];
        cout << ". " << endl;
      }
      if (config->GetFFD_Blending() == BEZIER){
        cout << "FFD Blending using Bezier Curves." << endl;
      }
    }

    FFDBox[iFFDBox] = new CFreeFormDefBox(degree, SplineOrder, config->GetFFD_Blending());
    FFDBox[iFFDBox]->SetTag(TagFFDBox); FFDBox[iFFDBox]->SetLevel(LevelFFDBox);

    /*--- Read the number of parents boxes ---*/
    
    nParentFFDBox = 0; // Nested FFD is not active
    if (rank == MASTER_NODE) cout << "Number of parent boxes: " << nParentFFDBox <<". ";
    
    for (iParentFFDBox = 0; iParentFFDBox < nParentFFDBox; iParentFFDBox++) {
      string ParentFFDBox = "NONE"; // Nested FFD is not active
      FFDBox[iFFDBox]->SetParentFFDBox(ParentFFDBox);
    }
    
    /*--- Read the number of children boxes ---*/
    
    nChildFFDBox = 0; // Nested FFD is not active
    if (rank == MASTER_NODE) cout << "Number of child boxes: " << nChildFFDBox <<"." << endl;
    
    for (iChildFFDBox = 0; iChildFFDBox < nChildFFDBox; iChildFFDBox++) {
      string ChildFFDBox = "NONE"; // Nested FFD is not active
      FFDBox[iFFDBox]->SetChildFFDBox(ChildFFDBox);
    }
    
    /*--- Read the number of the corner points ---*/
    
    nCornerPoints[iFFDBox] = 8;
    
    /*--- Read the coordinates of the corner points ---*/
    
    for (iCornerPoints = 0; iCornerPoints < nCornerPoints[iFFDBox]; iCornerPoints++) {
      
      if (nDim == 2) {

        if (polar) {

          coord[0] = config->GetCoordFFDBox(iFFDBox, 1*3);
          coord[1] = cos(0.1)*config->GetCoordFFDBox(iFFDBox, 1*3+1);
          coord[2] = -sin(0.1)*config->GetCoordFFDBox(iFFDBox, 1*3+1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 0);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 2*3);
          coord[1] = cos(0.1)*config->GetCoordFFDBox(iFFDBox, 2*3+1);
          coord[2] = -sin(0.1)*config->GetCoordFFDBox(iFFDBox, 2*3+1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 1);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 2*3);
          coord[1] = cos(0.1)*config->GetCoordFFDBox(iFFDBox, 2*3+1);
          coord[2] = sin(0.1)*config->GetCoordFFDBox(iFFDBox, 2*3+1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 2);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 1*3);
          coord[1] = cos(0.1)*config->GetCoordFFDBox(iFFDBox, 1*3+1);
          coord[2] = sin(0.1)*config->GetCoordFFDBox(iFFDBox, 1*3+1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 3);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 0*3);
          coord[1] = cos(0.1)*config->GetCoordFFDBox(iFFDBox, 0*3+1);
          coord[2] = -sin(0.1)*config->GetCoordFFDBox(iFFDBox, 0*3+1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 4);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 3*3);
          coord[1] = cos(0.1)*config->GetCoordFFDBox(iFFDBox, 3*3+1);
          coord[2] = -sin(0.1)*config->GetCoordFFDBox(iFFDBox, 3*3+1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 5);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 3*3);
          coord[1] = cos(0.1)*config->GetCoordFFDBox(iFFDBox, 3*3+1);
          coord[2] = sin(0.1)*config->GetCoordFFDBox(iFFDBox, 3*3+1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 6);

          coord[0] = config->GetCoordFFDBox(iFFDBox, 0*3);
          coord[1] = cos(0.1)*config->GetCoordFFDBox(iFFDBox, 0*3+1);
          coord[2] = sin(0.1)*config->GetCoordFFDBox(iFFDBox, 0*3+1);
          FFDBox[iFFDBox]->SetCoordCornerPoints(coord, 7);

        }

        else {
          if (iCornerPoints < nCornerPoints[iFFDBox]/SU2_TYPE::Int(2)) {
            coord[0] = config->GetCoordFFDBox(iFFDBox, iCornerPoints*3);
            coord[1] = config->GetCoordFFDBox(iFFDBox, iCornerPoints*3+1);
            coord[2] = -0.5;
          }
          else {
            coord[0] = FFDBox[iFFDBox]->GetCoordCornerPoints(0, iCornerPoints-nCornerPoints[iFFDBox]/SU2_TYPE::Int(2));
            coord[1] = FFDBox[iFFDBox]->GetCoordCornerPoints(1, iCornerPoints-nCornerPoints[iFFDBox]/SU2_TYPE::Int(2));
            coord[2] = 0.5;
          }
        }

      }
      else {
        coord[0] = config->GetCoordFFDBox(iFFDBox, iCornerPoints*3);
        coord[1] = config->GetCoordFFDBox(iFFDBox, iCornerPoints*3+1);
        coord[2] = config->GetCoordFFDBox(iFFDBox, iCornerPoints*3+2);
      }
      
      FFDBox[iFFDBox]->SetCoordCornerPoints(coord, iCornerPoints);
      
    }
    
    /*--- Method to identify if there is a FFDBox definition ---*/
    
    FFDBoxDefinition = false;
    
  }
  
  delete [] nCornerPoints;
  
  if (nFFDBox == 0) {
    SU2_MPI::Error("There is no FFD box definition. Check the config file.", CURRENT_FUNCTION);
  }
  
}

void CSurfaceMovement::MergeFFDInfo(CGeometry *geometry, CConfig *config) {
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned long iPoint;
  unsigned short iFFDBox;
  
#ifndef HAVE_MPI
  
  /*--- In serial, the single process has access to all geometry, so simply
   load the coordinates into the data structure. ---*/
  
  /*--- Total number of points in each FFD box. ---*/
  
  for (iFFDBox = 0 ; iFFDBox < nFFDBox; iFFDBox++) {
    
    /*--- Loop over the mesh to collect the coords of the local points. ---*/
    
    for (iPoint = 0; iPoint < FFDBox[iFFDBox]->GetnSurfacePoint(); iPoint++) {
      
      /*--- Retrieve the current parametric coordinates at this node. ---*/
      
      GlobalCoordX[iFFDBox].push_back(FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[0]);
      GlobalCoordY[iFFDBox].push_back(FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[1]);
      GlobalCoordZ[iFFDBox].push_back(FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[2]);
      GlobalPoint[iFFDBox].push_back(FFDBox[iFFDBox]->Get_PointIndex(iPoint));
      
      /*--- Marker of the boundary in the local domain. ---*/
      
      unsigned short MarkerIndex = FFDBox[iFFDBox]->Get_MarkerIndex(iPoint);
      string TagBound = config->GetMarker_All_TagBound(MarkerIndex);
      
      /*--- Find the Marker of the boundary in the config file. ---*/
      
      unsigned short MarkerIndex_CfgFile = config->GetMarker_CfgFile_TagBound(TagBound);
      string TagBound_CfgFile = config->GetMarker_CfgFile_TagBound(MarkerIndex_CfgFile);
      
      /*--- Set the value of the tag at this node. ---*/
      
      GlobalTag[iFFDBox].push_back(TagBound_CfgFile);
      
    }
    
  }
  
#else
  
  /*--- MPI preprocessing ---*/
  
  int iProcessor, nProcessor = size;
  
  /*--- Local variables needed for merging the geometry with MPI. ---*/
  
  unsigned long jPoint, iPointLocal;
  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long nBuffer_Scalar = 0;
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];
  
  for (iFFDBox = 0 ; iFFDBox < nFFDBox; iFFDBox++) {
    
    nLocalPoint = 0;
    for (iPoint = 0; iPoint < FFDBox[iFFDBox]->GetnSurfacePoint(); iPoint++) {
      
      iPointLocal = FFDBox[iFFDBox]->Get_PointIndex(iPoint);
      
      if (iPointLocal < geometry->GetnPointDomain()) {
        nLocalPoint++;
      }
      
    }
    Buffer_Send_nPoint[0] = nLocalPoint;

    /*--- Communicate the total number of nodes on this domain. ---*/
    
    SU2_MPI::Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG,
               Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    
    nBuffer_Scalar = MaxLocalPoint;

    /*--- Send and Recv buffers. ---*/
    
    su2double *Buffer_Send_X = new su2double[MaxLocalPoint];
    su2double *Buffer_Recv_X = NULL;
    
    su2double *Buffer_Send_Y = new su2double[MaxLocalPoint];
    su2double *Buffer_Recv_Y = NULL;
    
    su2double *Buffer_Send_Z = new su2double[MaxLocalPoint];
    su2double *Buffer_Recv_Z = NULL;
    
    unsigned long *Buffer_Send_Point = new unsigned long[MaxLocalPoint];
    unsigned long *Buffer_Recv_Point = NULL;
    
    unsigned short *Buffer_Send_MarkerIndex_CfgFile = new unsigned short[MaxLocalPoint];
    unsigned short *Buffer_Recv_MarkerIndex_CfgFile = NULL;
    
    /*--- Prepare the receive buffers in the master node only. ---*/
    
    if (rank == MASTER_NODE) {
      
      Buffer_Recv_X = new su2double[nProcessor*MaxLocalPoint];
      Buffer_Recv_Y = new su2double[nProcessor*MaxLocalPoint];
      Buffer_Recv_Z = new su2double[nProcessor*MaxLocalPoint];
      Buffer_Recv_Point = new unsigned long[nProcessor*MaxLocalPoint];
      Buffer_Recv_MarkerIndex_CfgFile = new unsigned short[nProcessor*MaxLocalPoint];
      
    }
    
    /*--- Main communication routine. Loop over each coordinate and perform
     the MPI comm. Temporary 1-D buffers are used to send the coordinates at
     all nodes on each partition to the master node. These are then unpacked
     by the master and sorted by global index in one large n-dim. array. ---*/
    
    /*--- Loop over this partition to collect the coords of the local points. ---*/
    
    jPoint = 0;
    for (iPoint = 0; iPoint < FFDBox[iFFDBox]->GetnSurfacePoint(); iPoint++) {
      
      iPointLocal = FFDBox[iFFDBox]->Get_PointIndex(iPoint);
      
      if (iPointLocal < geometry->GetnPointDomain()) {
        
        /*--- Load local coords into the temporary send buffer. ---*/
        
        Buffer_Send_X[jPoint] = FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[0];
        Buffer_Send_Y[jPoint] = FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[1];
        Buffer_Send_Z[jPoint] = FFDBox[iFFDBox]->Get_ParametricCoord(iPoint)[2];
        
        /*--- Store the global index for this local node. ---*/
        
        Buffer_Send_Point[jPoint] = geometry->node[FFDBox[iFFDBox]->Get_PointIndex(iPoint)]->GetGlobalIndex();
        
        /*--- Marker of the boundary in the local domain. ---*/
        
        unsigned short MarkerIndex = FFDBox[iFFDBox]->Get_MarkerIndex(iPoint);
        string TagBound = config->GetMarker_All_TagBound(MarkerIndex);
        
        /*--- Find the Marker of the boundary in the config file.---*/
        
        unsigned short MarkerIndex_CfgFile = config->GetMarker_CfgFile_TagBound(TagBound);
        Buffer_Send_MarkerIndex_CfgFile[jPoint] = MarkerIndex_CfgFile;
        
        jPoint++;
        
      }
      
    }
    
    /*--- Gather the coordinate data on the master node using MPI. ---*/
    
    SU2_MPI::Gather(Buffer_Send_X, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_X, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_Y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_Z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_Point, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Point, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_MarkerIndex_CfgFile, nBuffer_Scalar, MPI_UNSIGNED_SHORT, Buffer_Recv_MarkerIndex_CfgFile, nBuffer_Scalar, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);

    /*--- The master node unpacks and sorts this variable by global index ---*/
    
    if (rank == MASTER_NODE) {
      
      jPoint = 0;
      
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          
          GlobalCoordX[iFFDBox].push_back(Buffer_Recv_X[jPoint]);
          GlobalCoordY[iFFDBox].push_back(Buffer_Recv_Y[jPoint]);
          GlobalCoordZ[iFFDBox].push_back(Buffer_Recv_Z[jPoint]);
          GlobalPoint[iFFDBox].push_back(Buffer_Recv_Point[jPoint]);
          
          string TagBound_CfgFile = config->GetMarker_CfgFile_TagBound(Buffer_Recv_MarkerIndex_CfgFile[jPoint]);
          GlobalTag[iFFDBox].push_back(TagBound_CfgFile);
          jPoint++;
          
        }
        
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        
        jPoint = (iProcessor+1)*nBuffer_Scalar;
        
      }
    }
    
    /*--- Immediately release the temporary data buffers. ---*/
    
    delete [] Buffer_Send_X;
    delete [] Buffer_Send_Y;
    delete [] Buffer_Send_Z;
    delete [] Buffer_Send_Point;
    delete [] Buffer_Send_MarkerIndex_CfgFile;
    
    if (rank == MASTER_NODE) {
      delete [] Buffer_Recv_X;
      delete [] Buffer_Recv_Y;
      delete [] Buffer_Recv_Z;
      delete [] Buffer_Recv_Point;
      delete [] Buffer_Recv_MarkerIndex_CfgFile;
    }
    
  }

  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nPoint;
  }
  
#endif
  
}

void CSurfaceMovement::WriteFFDInfo(CSurfaceMovement** surface_movement, CGeometry **geometry, CConfig **config) {
  
  
  unsigned short iOrder, jOrder, kOrder, iFFDBox, iCornerPoints, iParentFFDBox, iChildFFDBox, iZone;
  unsigned long iSurfacePoints;
  char cstr[MAX_STRING_SIZE], mesh_file[MAX_STRING_SIZE];
  string str;
  ofstream output_file;
  su2double *coord;
  string text_line;
  
  bool polar = (config[ZONE_0]->GetFFD_CoordSystem() == POLAR);

  unsigned short nDim = geometry[ZONE_0]->GetnDim();
  
  for (iZone = 0; iZone < config[ZONE_0]->GetnZone(); iZone++){

    /*--- Merge the parallel FFD info ---*/
  
    surface_movement[iZone]->MergeFFDInfo(geometry[iZone], config[iZone]);

    if (iZone > 0){

      /* --- Merge the per-zone FFD info from the other zones into ZONE_0 ---*/

      for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++){

        surface_movement[ZONE_0]->GlobalCoordX[iFFDBox].insert(surface_movement[ZONE_0]->GlobalCoordX[iFFDBox].end(),
                                                               surface_movement[iZone]->GlobalCoordX[iFFDBox].begin(),
                                                               surface_movement[iZone]->GlobalCoordX[iFFDBox].end());
        surface_movement[ZONE_0]->GlobalCoordY[iFFDBox].insert(surface_movement[ZONE_0]->GlobalCoordY[iFFDBox].end(),
                                                               surface_movement[iZone]->GlobalCoordY[iFFDBox].begin(),
                                                               surface_movement[iZone]->GlobalCoordY[iFFDBox].end());
        surface_movement[ZONE_0]->GlobalCoordZ[iFFDBox].insert(surface_movement[ZONE_0]->GlobalCoordZ[iFFDBox].end(),
                                                               surface_movement[iZone]->GlobalCoordZ[iFFDBox].begin(),
                                                               surface_movement[iZone]->GlobalCoordZ[iFFDBox].end());
        surface_movement[ZONE_0]->GlobalTag[iFFDBox].insert(surface_movement[ZONE_0]->GlobalTag[iFFDBox].end(),
                                                               surface_movement[iZone]->GlobalTag[iFFDBox].begin(),
                                                               surface_movement[iZone]->GlobalTag[iFFDBox].end());
        surface_movement[ZONE_0]->GlobalPoint[iFFDBox].insert(surface_movement[ZONE_0]->GlobalPoint[iFFDBox].end(),
                                                               surface_movement[iZone]->GlobalPoint[iFFDBox].begin(),
                                                               surface_movement[iZone]->GlobalPoint[iFFDBox].end());
      }
    }
  }



  
  /*--- Attach to the mesh file the FFD information (all information is in ZONE_0) ---*/
  
  if (rank == MASTER_NODE) {
    
    /*--- Read the name of the output file ---*/
    
    str = config[ZONE_0]->GetMesh_Out_FileName();
    
    unsigned short lastindex = str.find_last_of(".");
    str = str.substr(0, lastindex);
    
    str += ".su2";
    
    strcpy (mesh_file, str.c_str());
    strcpy (cstr, mesh_file);
    
    output_file.precision(15);
    output_file.open(cstr, ios::out | ios::app);
    
    if (nFFDBox != 0) {
      output_file << "FFD_NBOX= " << nFFDBox << endl;
      output_file << "FFD_NLEVEL= " << nLevel << endl;
    }
    
    for (iFFDBox = 0 ; iFFDBox < nFFDBox; iFFDBox++) {
      
      output_file << "FFD_TAG= " << FFDBox[iFFDBox]->GetTag() << endl;
      output_file << "FFD_LEVEL= " << FFDBox[iFFDBox]->GetLevel() << endl;
      
      output_file << "FFD_DEGREE_I= " << FFDBox[iFFDBox]->GetlOrder()-1 << endl;
      if (polar) output_file << "FFD_DEGREE_J= " << FFDBox[iFFDBox]->GetnOrder()-1 << endl;
      else output_file << "FFD_DEGREE_J= " << FFDBox[iFFDBox]->GetmOrder()-1 << endl;
      if (nDim == 3) output_file << "FFD_DEGREE_K= " << FFDBox[iFFDBox]->GetnOrder()-1 << endl;
      if (config[ZONE_0]->GetFFD_Blending() == BSPLINE_UNIFORM) {
        output_file << "FFD_BLENDING= BSPLINE_UNIFORM" << endl;
        output_file << "BSPLINE_ORDER_I= " << FFDBox[iFFDBox]->BlendingFunction[0]->GetOrder() << endl;
        if (polar) output_file << "BSPLINE_ORDER_J= " << FFDBox[iFFDBox]->BlendingFunction[2]->GetOrder() << endl;
        else output_file << "BSPLINE_ORDER_J= " << FFDBox[iFFDBox]->BlendingFunction[1]->GetOrder() << endl;
        if (nDim == 3) output_file << "BSPLINE_ORDER_K= " << FFDBox[iFFDBox]->BlendingFunction[2]->GetOrder() << endl;
      }
      if (config[ZONE_0]->GetFFD_Blending() == BEZIER) {
        output_file << "FFD_BLENDING= BEZIER" << endl;
      }

      output_file << "FFD_PARENTS= " << FFDBox[iFFDBox]->GetnParentFFDBox() << endl;
      for (iParentFFDBox = 0; iParentFFDBox < FFDBox[iFFDBox]->GetnParentFFDBox(); iParentFFDBox++)
        output_file << FFDBox[iFFDBox]->GetParentFFDBoxTag(iParentFFDBox) << endl;
      output_file << "FFD_CHILDREN= " << FFDBox[iFFDBox]->GetnChildFFDBox() << endl;
      for (iChildFFDBox = 0; iChildFFDBox < FFDBox[iFFDBox]->GetnChildFFDBox(); iChildFFDBox++)
        output_file << FFDBox[iFFDBox]->GetChildFFDBoxTag(iChildFFDBox) << endl;
      
      if (nDim == 2) {
        output_file << "FFD_CORNER_POINTS= " << FFDBox[iFFDBox]->GetnCornerPoints()/SU2_TYPE::Int(2) << endl;
        if (polar) {
          coord = FFDBox[iFFDBox]->GetCoordCornerPoints(4);
          output_file << coord[0] << "\t" << sqrt(coord[1]*coord[1]+coord[2]*coord[2]) << endl;

          coord = FFDBox[iFFDBox]->GetCoordCornerPoints(0);
          output_file << coord[0] << "\t" << sqrt(coord[1]*coord[1]+coord[2]*coord[2]) << endl;

          coord = FFDBox[iFFDBox]->GetCoordCornerPoints(1);
          output_file << coord[0] << "\t" << sqrt(coord[1]*coord[1]+coord[2]*coord[2]) << endl;

          coord = FFDBox[iFFDBox]->GetCoordCornerPoints(5);
          output_file << coord[0] << "\t" << sqrt(coord[1]*coord[1]+coord[2]*coord[2]) << endl;
        }
        else {
          for (iCornerPoints = 0; iCornerPoints < FFDBox[iFFDBox]->GetnCornerPoints()/SU2_TYPE::Int(2); iCornerPoints++) {
            coord = FFDBox[iFFDBox]->GetCoordCornerPoints(iCornerPoints);
            output_file << coord[0] << "\t" << coord[1] << endl;
          }
        }
      }
      else {
        output_file << "FFD_CORNER_POINTS= " << FFDBox[iFFDBox]->GetnCornerPoints() << endl;
        for (iCornerPoints = 0; iCornerPoints < FFDBox[iFFDBox]->GetnCornerPoints(); iCornerPoints++) {
          coord = FFDBox[iFFDBox]->GetCoordCornerPoints(iCornerPoints);
          output_file << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
        }
      }
      
      /*--- Writing control points ---*/
      
      if (FFDBox[iFFDBox]->GetnControlPoints() == 0) {
        output_file << "FFD_CONTROL_POINTS= 0" << endl;
      }
      else {
        output_file << "FFD_CONTROL_POINTS= " << FFDBox[iFFDBox]->GetnControlPoints() << endl;
        for (iOrder = 0; iOrder < FFDBox[iFFDBox]->GetlOrder(); iOrder++)
          for (jOrder = 0; jOrder < FFDBox[iFFDBox]->GetmOrder(); jOrder++)
            for (kOrder = 0; kOrder < FFDBox[iFFDBox]->GetnOrder(); kOrder++) {
              coord = FFDBox[iFFDBox]->GetCoordControlPoints(iOrder, jOrder, kOrder);
              output_file << iOrder << "\t" << jOrder << "\t" << kOrder << "\t" << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
            }
      }
      
      /*--- Writing surface points ---*/
      
      if (FFDBox[iFFDBox]->GetnControlPoints() == 0) {
        output_file << "FFD_SURFACE_POINTS= 0" << endl;
      }
      else {
        output_file << "FFD_SURFACE_POINTS= " << GlobalTag[iFFDBox].size() << endl;
        
        for (iSurfacePoints = 0; iSurfacePoints < GlobalTag[iFFDBox].size(); iSurfacePoints++) {
          output_file << scientific << GlobalTag[iFFDBox][iSurfacePoints] << "\t" << GlobalPoint[iFFDBox][iSurfacePoints]
          << "\t" << GlobalCoordX[iFFDBox][iSurfacePoints] << "\t" << GlobalCoordY[iFFDBox][iSurfacePoints]
          << "\t" << GlobalCoordZ[iFFDBox][iSurfacePoints] << endl;
        }
        
      }
      
    }
    
    output_file.close();
    
  }

}

CFreeFormDefBox::CFreeFormDefBox(void) : CGridMovement() { }

CFreeFormDefBox::CFreeFormDefBox(unsigned short Degree[], unsigned short BSplineOrder[], unsigned short kind_blending) : CGridMovement() {
  
  unsigned short iCornerPoints, iOrder, jOrder, kOrder, iDim;
  
  /*--- FFD is always 3D (even in 2D problems) ---*/
  
  nDim = 3;
  nCornerPoints = 8;
  
  /*--- Allocate Corners points ---*/
  
  Coord_Corner_Points = new su2double* [nCornerPoints];
  for (iCornerPoints = 0; iCornerPoints < nCornerPoints; iCornerPoints++)
    Coord_Corner_Points[iCornerPoints] = new su2double [nDim];
  
  ParamCoord = new su2double[nDim]; ParamCoord_ = new su2double[nDim];
  cart_coord = new su2double[nDim]; cart_coord_ = new su2double[nDim];
  Gradient = new su2double[nDim];

  lDegree = Degree[0]; lOrder = lDegree+1;
  mDegree = Degree[1]; mOrder = mDegree+1;
  nDegree = Degree[2]; nOrder = nDegree+1;
  nControlPoints = lOrder*mOrder*nOrder;
  
  lDegree_Copy = Degree[0]; lOrder_Copy = lDegree+1;
  mDegree_Copy = Degree[1]; mOrder_Copy = mDegree+1;
  nDegree_Copy = Degree[2]; nOrder_Copy = nDegree+1;
  nControlPoints_Copy = lOrder_Copy*mOrder_Copy*nOrder_Copy;
  
  Coord_Control_Points = new su2double*** [lOrder];
  ParCoord_Control_Points = new su2double*** [lOrder];
  Coord_Control_Points_Copy = new su2double*** [lOrder];
  for (iOrder = 0; iOrder < lOrder; iOrder++) {
    Coord_Control_Points[iOrder] = new su2double** [mOrder];
    ParCoord_Control_Points[iOrder] = new su2double** [mOrder];
    Coord_Control_Points_Copy[iOrder] = new su2double** [mOrder];
    for (jOrder = 0; jOrder < mOrder; jOrder++) {
      Coord_Control_Points[iOrder][jOrder] = new su2double* [nOrder];
      ParCoord_Control_Points[iOrder][jOrder] = new su2double* [nOrder];
      Coord_Control_Points_Copy[iOrder][jOrder] = new su2double* [nOrder];
      for (kOrder = 0; kOrder < nOrder; kOrder++) {
        Coord_Control_Points[iOrder][jOrder][kOrder] = new su2double [nDim];
        ParCoord_Control_Points[iOrder][jOrder][kOrder] = new su2double [nDim];
        Coord_Control_Points_Copy[iOrder][jOrder][kOrder] = new su2double [nDim];
        for (iDim = 0; iDim < nDim; iDim++) {
          Coord_Control_Points[iOrder][jOrder][kOrder][iDim] = 0.0;
          ParCoord_Control_Points[iOrder][jOrder][kOrder][iDim] = 0.0;
          Coord_Control_Points_Copy[iOrder][jOrder][kOrder][iDim] = 0.0;
        }
      }
    }
  }

  BlendingFunction = new CFreeFormBlending*[nDim];

  if (kind_blending == BEZIER){
    BlendingFunction[0] = new CBezierBlending(lOrder, lOrder);
    BlendingFunction[1] = new CBezierBlending(mOrder, mOrder);
    BlendingFunction[2] = new CBezierBlending(nOrder, nOrder);
  }
  if (kind_blending == BSPLINE_UNIFORM){
    BlendingFunction[0] = new CBSplineBlending(BSplineOrder[0], lOrder);
    BlendingFunction[1] = new CBSplineBlending(BSplineOrder[1], mOrder);
    BlendingFunction[2] = new CBSplineBlending(BSplineOrder[2], nOrder);
  }

}

CFreeFormDefBox::~CFreeFormDefBox(void) {
  unsigned short iOrder, jOrder, kOrder, iCornerPoints, iDim;
  
  for (iOrder = 0; iOrder < lOrder; iOrder++) 
    for (jOrder = 0; jOrder < mOrder; jOrder++) 
      for (kOrder = 0; kOrder < nOrder; kOrder++) {
        delete [] Coord_Control_Points[iOrder][jOrder][kOrder];
        delete [] ParCoord_Control_Points[iOrder][jOrder][kOrder];
        delete [] Coord_Control_Points_Copy[iOrder][jOrder][kOrder];
      }
  delete [] Coord_Control_Points;
  delete [] ParCoord_Control_Points;
  delete [] Coord_Control_Points_Copy;

  delete [] ParamCoord;
  delete [] cart_coord;
  delete [] Gradient;
  
  for (iCornerPoints = 0; iCornerPoints < nCornerPoints; iCornerPoints++)
    delete [] Coord_Corner_Points[iCornerPoints];
  delete [] Coord_Corner_Points;

  for (iDim = 0; iDim < nDim; iDim++){
    delete BlendingFunction[iDim];
  }
  delete [] BlendingFunction;
}

void  CFreeFormDefBox::SetUnitCornerPoints(void) {
  
  unsigned short iDim;
  su2double *coord = new su2double [nDim];
  
  for (iDim = 0; iDim < nDim; iDim++) coord[iDim] = 0.0;
  
  coord [0] = 0.0; coord [1] = 0.0; coord [2] = 0.0; this->SetCoordCornerPoints(coord, 0);
  coord [0] = 1.0; coord [1] = 0.0; coord [2] = 0.0; this->SetCoordCornerPoints(coord, 1);
  coord [0] = 1.0; coord [1] = 1.0; coord [2] = 0.0; this->SetCoordCornerPoints(coord, 2);
  coord [0] = 0.0; coord [1] = 1.0; coord [2] = 0.0; this->SetCoordCornerPoints(coord, 3);
  coord [0] = 0.0; coord [1] = 0.0; coord [2] = 1.0; this->SetCoordCornerPoints(coord, 4);
  coord [0] = 1.0; coord [1] = 0.0; coord [2] = 1.0; this->SetCoordCornerPoints(coord, 5);
  coord [0] = 1.0; coord [1] = 1.0; coord [2] = 1.0; this->SetCoordCornerPoints(coord, 6);
  coord [0] = 0.0; coord [1] = 1.0; coord [2] = 1.0; this->SetCoordCornerPoints(coord, 7);
  
  delete [] coord;
  
}

void CFreeFormDefBox::SetControlPoints_Parallelepiped (void) {
  unsigned short iDim, iDegree, jDegree, kDegree;
  
  /*--- Set base control points according to the notation of Vtk for hexahedrons ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    Coord_Control_Points  [0]      [0]      [0]      [iDim]  = Coord_Corner_Points[0][iDim];
    Coord_Control_Points  [lOrder-1]  [0]      [0]      [iDim]  = Coord_Corner_Points[1][iDim];
    Coord_Control_Points  [lOrder-1]  [mOrder-1]  [0]      [iDim]  = Coord_Corner_Points[2][iDim];
    Coord_Control_Points  [0]      [mOrder-1]  [0]      [iDim]  = Coord_Corner_Points[3][iDim];
    Coord_Control_Points  [0]      [0]      [nOrder-1]  [iDim]  = Coord_Corner_Points[4][iDim];
    Coord_Control_Points  [lOrder-1]  [0]      [nOrder-1]  [iDim]  = Coord_Corner_Points[5][iDim];
    Coord_Control_Points  [lOrder-1]  [mOrder-1]  [nOrder-1]  [iDim]  = Coord_Corner_Points[6][iDim];
    Coord_Control_Points  [0]      [mOrder-1]  [nOrder-1]  [iDim]  = Coord_Corner_Points[7][iDim];
  }
  
  /*--- Fill the rest of the cubic matrix of control points with uniform spacing (parallelepiped) ---*/
  for (iDegree = 0; iDegree <= lDegree; iDegree++)
    for (jDegree = 0; jDegree <= mDegree; jDegree++)
      for (kDegree = 0; kDegree <= nDegree; kDegree++) {
        Coord_Control_Points[iDegree][jDegree][kDegree][0] = Coord_Corner_Points[0][0] 
        + su2double(iDegree)/su2double(lDegree)*(Coord_Corner_Points[1][0]-Coord_Corner_Points[0][0]);
        Coord_Control_Points[iDegree][jDegree][kDegree][1] = Coord_Corner_Points[0][1] 
        + su2double(jDegree)/su2double(mDegree)*(Coord_Corner_Points[3][1]-Coord_Corner_Points[0][1]);        
        Coord_Control_Points[iDegree][jDegree][kDegree][2] = Coord_Corner_Points[0][2] 
        + su2double(kDegree)/su2double(nDegree)*(Coord_Corner_Points[4][2]-Coord_Corner_Points[0][2]);
      }
}

void CFreeFormDefBox::SetSupportCP(CFreeFormDefBox *FFDBox) {
  unsigned short iDim, iOrder, jOrder, kOrder;
  unsigned short lOrder = FFDBox->GetlOrder();
  unsigned short mOrder = FFDBox->GetmOrder();
  unsigned short nOrder = FFDBox->GetnOrder();
  
  Coord_SupportCP = new su2double*** [lOrder];
  for (iOrder = 0; iOrder < lOrder; iOrder++) {
    Coord_SupportCP[iOrder] = new su2double** [mOrder];
    for (jOrder = 0; jOrder < mOrder; jOrder++) {
      Coord_SupportCP[iOrder][jOrder] = new su2double* [nOrder];
      for (kOrder = 0; kOrder < nOrder; kOrder++)
        Coord_SupportCP[iOrder][jOrder][kOrder] = new su2double [nDim];
    }
  }
  
  /*--- Set base support control points according to the notation of Vtk for hexahedrons ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    Coord_SupportCP  [0]      [0]      [0]      [iDim]  = Coord_Corner_Points[0][iDim];
    Coord_SupportCP  [lOrder-1]  [0]      [0]      [iDim]  = Coord_Corner_Points[1][iDim];
    Coord_SupportCP  [lOrder-1]  [mOrder-1]  [0]      [iDim]  = Coord_Corner_Points[2][iDim];
    Coord_SupportCP  [0]      [mOrder-1]  [0]      [iDim]  = Coord_Corner_Points[3][iDim];
    Coord_SupportCP  [0]      [0]      [nOrder-1]  [iDim]  = Coord_Corner_Points[4][iDim];
    Coord_SupportCP  [lOrder-1]  [0]      [nOrder-1]  [iDim]  = Coord_Corner_Points[5][iDim];
    Coord_SupportCP  [lOrder-1]  [mOrder-1]  [nOrder-1]  [iDim]  = Coord_Corner_Points[6][iDim];
    Coord_SupportCP  [0]      [mOrder-1]  [nOrder-1]  [iDim]  = Coord_Corner_Points[7][iDim];
  }
  
  /*--- Fill the rest of the cubic matrix of support control points with uniform spacing  ---*/
  for (iOrder = 0; iOrder < lOrder; iOrder++)
    for (jOrder = 0; jOrder < mOrder; jOrder++)
      for (kOrder = 0; kOrder < nOrder; kOrder++) {
        Coord_SupportCP[iOrder][jOrder][kOrder][0] = Coord_Corner_Points[0][0] 
        + su2double(iOrder)/su2double(lOrder-1)*(Coord_Corner_Points[1][0]-Coord_Corner_Points[0][0]);
        Coord_SupportCP[iOrder][jOrder][kOrder][1] = Coord_Corner_Points[0][1] 
        + su2double(jOrder)/su2double(mOrder-1)*(Coord_Corner_Points[3][1]-Coord_Corner_Points[0][1]);        
        Coord_SupportCP[iOrder][jOrder][kOrder][2] = Coord_Corner_Points[0][2] 
        + su2double(kOrder)/su2double(nOrder-1)*(Coord_Corner_Points[4][2]-Coord_Corner_Points[0][2]);
      }
}

void CFreeFormDefBox::SetSupportCPChange(CFreeFormDefBox *FFDBox) {
  unsigned short iDim, iOrder, jOrder, kOrder;
  su2double *CartCoordNew, *ParamCoord;
  unsigned short lOrder = FFDBox->GetlOrder();
  unsigned short mOrder = FFDBox->GetmOrder();
  unsigned short nOrder = FFDBox->GetnOrder();

  su2double ****ParamCoord_SupportCP = new su2double*** [lOrder];
  for (iOrder = 0; iOrder < lOrder; iOrder++) {
    ParamCoord_SupportCP[iOrder] = new su2double** [mOrder];
    for (jOrder = 0; jOrder < mOrder; jOrder++) {
      ParamCoord_SupportCP[iOrder][jOrder] = new su2double* [nOrder];
      for (kOrder = 0; kOrder < nOrder; kOrder++)
        ParamCoord_SupportCP[iOrder][jOrder][kOrder] = new su2double [nDim];
    }
  }
  
  for (iOrder = 0; iOrder < lOrder; iOrder++)
    for (jOrder = 0; jOrder < mOrder; jOrder++)
      for (kOrder = 0; kOrder < nOrder; kOrder++)
        for (iDim = 0; iDim < nDim; iDim++)
          ParamCoord_SupportCP[iOrder][jOrder][kOrder][iDim] = 
          Coord_SupportCP[iOrder][jOrder][kOrder][iDim];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Coord_Control_Points[0][0][0][iDim]  = FFDBox->GetCoordCornerPoints(iDim, 0);
    Coord_Control_Points[1][0][0][iDim]  = FFDBox->GetCoordCornerPoints(iDim, 1);
    Coord_Control_Points[1][1][0][iDim]  = FFDBox->GetCoordCornerPoints(iDim, 2);
    Coord_Control_Points[0][1][0][iDim]  = FFDBox->GetCoordCornerPoints(iDim, 3);
    Coord_Control_Points[0][0][1][iDim]  = FFDBox->GetCoordCornerPoints(iDim, 4);
    Coord_Control_Points[1][0][1][iDim]  = FFDBox->GetCoordCornerPoints(iDim, 5);
    Coord_Control_Points[1][1][1][iDim]  = FFDBox->GetCoordCornerPoints(iDim, 6);
    Coord_Control_Points[0][1][1][iDim]  = FFDBox->GetCoordCornerPoints(iDim, 7);
  }
  
  for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++) {
    for (jOrder = 0; jOrder < FFDBox->GetmOrder(); jOrder++) {
      for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
        ParamCoord = ParamCoord_SupportCP[iOrder][jOrder][kOrder];
        CartCoordNew = EvalCartesianCoord(ParamCoord);
        FFDBox->SetCoordControlPoints(CartCoordNew, iOrder, jOrder, kOrder);
        FFDBox->SetCoordControlPoints_Copy(CartCoordNew, iOrder, jOrder, kOrder);
      }
    }
  }

}

void CFreeFormDefBox::SetCart2Cyl_ControlPoints(CConfig *config) {
  
  unsigned short iDegree, jDegree, kDegree;
  su2double CartCoord[3];
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
  
  X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
  
  for (kDegree = 0; kDegree <= nDegree; kDegree++) {
    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
      for (iDegree = 0; iDegree <= lDegree; iDegree++) {
        
        CartCoord[0] = Coord_Control_Points[iDegree][jDegree][kDegree][0];
        CartCoord[1] = Coord_Control_Points[iDegree][jDegree][kDegree][1];
        CartCoord[2] = Coord_Control_Points[iDegree][jDegree][kDegree][2];
        
        Xbar =  CartCoord[0] - X_0; Ybar =  CartCoord[1] - Y_0; Zbar =  CartCoord[2] - Z_0;
        
        Coord_Control_Points[iDegree][jDegree][kDegree][0] = sqrt(Ybar*Ybar + Zbar*Zbar);
        Coord_Control_Points[iDegree][jDegree][kDegree][1] = atan2 ( Zbar, Ybar);
        if (Coord_Control_Points[iDegree][jDegree][kDegree][1] > PI_NUMBER/2.0) Coord_Control_Points[iDegree][jDegree][kDegree][1] -= 2.0*PI_NUMBER;
        Coord_Control_Points[iDegree][jDegree][kDegree][2] = Xbar;
        
        CartCoord[0] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][0];
        CartCoord[1] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1];
        CartCoord[2] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][2];
        
        Xbar =  CartCoord[0] - X_0; Ybar =  CartCoord[1] - Y_0; Zbar =  CartCoord[2] - Z_0;
        
        Coord_Control_Points_Copy[iDegree][jDegree][kDegree][0] = sqrt(Ybar*Ybar + Zbar*Zbar);
        Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1] = atan2 (Zbar, Ybar);
        if (Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1] > PI_NUMBER/2.0) Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1] -= 2.0*PI_NUMBER;
        Coord_Control_Points_Copy[iDegree][jDegree][kDegree][2] = Xbar;
        
      }
    }
  }
  
}

void CFreeFormDefBox::SetCyl2Cart_ControlPoints(CConfig *config) {
  
  unsigned short iDegree, jDegree, kDegree;
  su2double PolarCoord[3];
  
   su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
  X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
  
  for (kDegree = 0; kDegree <= nDegree; kDegree++) {
    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
      for (iDegree = 0; iDegree <= lDegree; iDegree++) {
        
        
         PolarCoord[0] = Coord_Control_Points[iDegree][jDegree][kDegree][0];
         PolarCoord[1] = Coord_Control_Points[iDegree][jDegree][kDegree][1];
         PolarCoord[2] = Coord_Control_Points[iDegree][jDegree][kDegree][2];
        
        
        Xbar = PolarCoord[2];
        Ybar = PolarCoord[0] * cos(PolarCoord[1]);
        Zbar = PolarCoord[0] * sin(PolarCoord[1]);
        
        PolarCoord[0] =  Xbar +X_0;  PolarCoord[1] = Ybar +Y_0; PolarCoord[2] = Zbar +Z_0;
        
        Coord_Control_Points[iDegree][jDegree][kDegree][0] = PolarCoord[0];
        Coord_Control_Points[iDegree][jDegree][kDegree][1] = PolarCoord[1];
        Coord_Control_Points[iDegree][jDegree][kDegree][2] = PolarCoord[2];
        
      }
    }
  }
  
}

void CFreeFormDefBox::SetCart2Cyl_CornerPoints(CConfig *config) {
  
  unsigned short iCornerPoint;
  su2double *CartCoord;
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
  
  X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
  
  for (iCornerPoint = 0; iCornerPoint < 8; iCornerPoint++) {
    
    CartCoord = GetCoordCornerPoints(iCornerPoint);
    Xbar =  CartCoord[0] - X_0; Ybar =  CartCoord[1] - Y_0; Zbar =  CartCoord[2] - Z_0;
    
    CartCoord[0] = sqrt(Ybar*Ybar + Zbar*Zbar);
    CartCoord[1] = atan2 ( Zbar, Ybar); if (CartCoord[1] > PI_NUMBER/2.0) CartCoord[1] -= 2.0*PI_NUMBER;
    CartCoord[2] =  Xbar;
    
  }
  
}


void CFreeFormDefBox::SetCyl2Cart_CornerPoints(CConfig *config) {
  
  unsigned short iCornerPoint;
  su2double *PolarCoord;
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
  
  X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
  
  for (iCornerPoint = 0; iCornerPoint < 8; iCornerPoint++) {
    
    PolarCoord = GetCoordCornerPoints(iCornerPoint);
    
    Xbar = PolarCoord[2];
    Ybar = PolarCoord[0] * cos(PolarCoord[1]);
    Zbar = PolarCoord[0] * sin(PolarCoord[1]);
    
    PolarCoord[0] =  Xbar + X_0;  PolarCoord[1] = Ybar + Y_0; PolarCoord[2] = Zbar + Z_0;
    
  }
  
}

void CFreeFormDefBox::SetCart2Sphe_ControlPoints(CConfig *config) {
  
  unsigned short iDegree, jDegree, kDegree;
  su2double CartCoord[3];
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
  
  X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
  
  for (kDegree = 0; kDegree <= nDegree; kDegree++) {
    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
      for (iDegree = 0; iDegree <= lDegree; iDegree++) {
        
        CartCoord[0] = Coord_Control_Points[iDegree][jDegree][kDegree][0];
        CartCoord[1] = Coord_Control_Points[iDegree][jDegree][kDegree][1];
        CartCoord[2] = Coord_Control_Points[iDegree][jDegree][kDegree][2];
        
        Xbar =  CartCoord[0] - X_0; Ybar =  CartCoord[1] - Y_0; Zbar =  CartCoord[2] - Z_0;
        
        Coord_Control_Points[iDegree][jDegree][kDegree][0] = sqrt(Xbar*Xbar + Ybar*Ybar + Zbar*Zbar);
        Coord_Control_Points[iDegree][jDegree][kDegree][1] = atan2 ( Zbar, Ybar);
        if (Coord_Control_Points[iDegree][jDegree][kDegree][1] > PI_NUMBER/2.0) Coord_Control_Points[iDegree][jDegree][kDegree][1] -= 2.0*PI_NUMBER;
        Coord_Control_Points[iDegree][jDegree][kDegree][2] = acos(Xbar/Coord_Control_Points[iDegree][jDegree][kDegree][0] );
        
        CartCoord[0] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][0];
        CartCoord[1] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1];
        CartCoord[2] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][2];
        
        Xbar =  CartCoord[0] - X_0; Ybar =  CartCoord[1] - Y_0; Zbar =  CartCoord[2] - Z_0;
        
        Coord_Control_Points_Copy[iDegree][jDegree][kDegree][0] = sqrt(Xbar*Xbar + Ybar*Ybar + Zbar*Zbar);
        Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1] = atan2 ( Zbar, Ybar);
        if (Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1] > PI_NUMBER/2.0)
          Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1]  -= 2.0*PI_NUMBER;
        Coord_Control_Points_Copy[iDegree][jDegree][kDegree][2] = acos(Xbar/Coord_Control_Points_Copy[iDegree][jDegree][kDegree][0]);
        
      }
    }
  }
  
}

void CFreeFormDefBox::SetSphe2Cart_ControlPoints(CConfig *config) {
  
  unsigned short iDegree, jDegree, kDegree;
  su2double PolarCoord[3];
  
   su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
  X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
  
  for (kDegree = 0; kDegree <= nDegree; kDegree++) {
    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
      for (iDegree = 0; iDegree <= lDegree; iDegree++) {
        
        PolarCoord[0] = Coord_Control_Points[iDegree][jDegree][kDegree][0];
        PolarCoord[1] = Coord_Control_Points[iDegree][jDegree][kDegree][1];
        PolarCoord[2] = Coord_Control_Points[iDegree][jDegree][kDegree][2];
        
        Xbar = PolarCoord[0] * cos(PolarCoord[2]);
        Ybar = PolarCoord[0] * cos(PolarCoord[1]) * sin(PolarCoord[2]);
        Zbar = PolarCoord[0] * sin(PolarCoord[1]) * sin(PolarCoord[2]);
        
        PolarCoord[0] = Xbar + X_0;  PolarCoord[1] = Ybar + Y_0; PolarCoord[2] = Zbar + Z_0;
        
        Coord_Control_Points[iDegree][jDegree][kDegree][0] = PolarCoord[0];
        Coord_Control_Points[iDegree][jDegree][kDegree][1] = PolarCoord[1];
        Coord_Control_Points[iDegree][jDegree][kDegree][2] = PolarCoord[2];
        
      }
    }
  }
  
}

void CFreeFormDefBox::SetCart2Sphe_CornerPoints(CConfig *config) {
  
  unsigned short iCornerPoint;
  su2double *CartCoord;
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
  
  X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
  
  for (iCornerPoint = 0; iCornerPoint < 8; iCornerPoint++) {
    
    CartCoord = GetCoordCornerPoints(iCornerPoint);
    Xbar =  CartCoord[0] - X_0; Ybar =  CartCoord[1] - Y_0; Zbar =  CartCoord[2] - Z_0;
    
    CartCoord[0] = sqrt(Xbar*Xbar + Ybar*Ybar + Zbar*Zbar);
    CartCoord[1] = atan2(Zbar, Ybar);  if (CartCoord[1] > PI_NUMBER/2.0) CartCoord[1] -= 2.0*PI_NUMBER;
    CartCoord[2] = acos(Xbar/CartCoord[0]);
    
  }
  
}


void CFreeFormDefBox::SetSphe2Cart_CornerPoints(CConfig *config) {
  
  unsigned short iCornerPoint;
  su2double *PolarCoord;
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
  
  X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
  
  for (iCornerPoint = 0; iCornerPoint < 8; iCornerPoint++) {
    
    PolarCoord = GetCoordCornerPoints(iCornerPoint);
    
    Xbar = PolarCoord[0] * cos(PolarCoord[2]);
    Ybar = PolarCoord[0] * cos(PolarCoord[1]) * sin(PolarCoord[2]);
    Zbar = PolarCoord[0] * sin(PolarCoord[1]) * sin(PolarCoord[2]);
    
    PolarCoord[0] =  Xbar + X_0;  PolarCoord[1] = Ybar + Y_0; PolarCoord[2] = Zbar + Z_0;
    
  }
  
}


void CFreeFormDefBox::SetCGNS(CGeometry *geometry, unsigned short iFFDBox, bool original) {
#ifdef HAVE_CGNS

  char FFDBox_filename[MAX_STRING_SIZE];
  bool new_file;
  unsigned short iDim, iDegree, jDegree, kDegree;
  unsigned short outDim;
  unsigned int pos;
  char zonename[33];
  int FFDBox_cgns_file;
  int cell_dim, phys_dim;
  int cgns_base=0, cgns_family, cgns_zone, cgns_err, dummy;
  const char * basename;

  /*--- FFD output is always 3D (even in 2D problems),
   this is important for debuging ---*/
  nDim = geometry->GetnDim();
  cell_dim = nDim,
  phys_dim = 3;

  SPRINTF (FFDBox_filename, "ffd_boxes.cgns");

  if ((original) && (iFFDBox == 0)) new_file = true;
  else new_file = false;
  
  if (new_file) {
    cgns_err = cg_open(FFDBox_filename, CG_MODE_WRITE, &FFDBox_cgns_file);
    if (cgns_err) cg_error_print();
    cgns_err = cg_descriptor_write("Title", "Visualization of the FFD boxes generated by SU2_DEF." );
    if (cgns_err) cg_error_print();
  }
  else {
    cgns_err = cg_open(FFDBox_filename, CG_MODE_MODIFY, &FFDBox_cgns_file);
    if (cgns_err) cg_error_print();
  }

  if (original) {
    basename = "Original_FFD";
  }
  else {
    basename = "Deformed_FFD";
  }

  if (iFFDBox == 0){
    cgns_err = cg_base_write(FFDBox_cgns_file, basename, cell_dim, phys_dim, &cgns_base);
    if (cgns_err) cg_error_print();
  }
  cgns_err = cg_family_write(FFDBox_cgns_file, cgns_base, Tag.c_str(), &cgns_family);
  if (cgns_err) cg_error_print();
  
  cgsize_t dims[9];
  dims[0] = lDegree+1;
  dims[1] = mDegree+1;
  if (cell_dim == 3){
     dims[2] = nDegree+1;
  }
  cgsize_t pointlen = 1;
  for(int ii=0; ii<cell_dim; ii++){
    dims[ii+cell_dim] = dims[ii] -1;
    dims[ii+2*cell_dim] = 0;
    pointlen *= dims[ii];
  }

  passivedouble *buffer = new passivedouble[pointlen];
  SPRINTF (zonename, "SU2_Zone_%d", SU2_TYPE::Int(iFFDBox));
  
  cgns_err = cg_zone_write(FFDBox_cgns_file, cgns_base, zonename, dims, CGNS_ENUMV(Structured), &cgns_zone);
  if (cgns_err) cg_error_print();
  cgns_err = cg_goto(FFDBox_cgns_file, cgns_base, zonename, 0, NULL);
  if (cgns_err) cg_error_print();
  cgns_err = cg_famname_write(Tag.c_str());
  if (cgns_err) cg_error_print();


  const char* coord_names[3] = { "CoordinateX", "CoordinateY", "CoordinateZ" };
  for (iDim=0; iDim<nDim; iDim++)
  {
    outDim = nDim == 2 ? 0: nDegree;
    for (kDegree = 0; kDegree <= outDim; kDegree++) {
      for (jDegree = 0; jDegree <= mDegree; jDegree++) {
        for (iDegree = 0; iDegree <= lDegree; iDegree++) {
          pos = iDegree + jDegree*(lDegree+1)+ kDegree*(lDegree+1)*(mDegree+1);
          buffer[pos] = SU2_TYPE::GetValue(Coord_Control_Points[iDegree][jDegree][kDegree][iDim]);
        }
      }
    }
    cgns_err = cg_coord_write(FFDBox_cgns_file, cgns_base, cgns_zone, CGNS_ENUMV(RealDouble), coord_names[iDim], buffer, &dummy);
    if (cgns_err) cg_error_print();

  }
  if (nDim==2){
    std::fill_n(buffer, pointlen, 0.0);
    cgns_err = cg_coord_write(FFDBox_cgns_file, cgns_base, cgns_zone, CGNS_ENUMV(RealDouble), coord_names[2], buffer, &dummy);
    if (cgns_err) cg_error_print();
  }

  delete [] buffer;
  cgns_err = cg_close(FFDBox_cgns_file);
  if (cgns_err) cg_error_print();
#else // Not built with CGNS support
  cout << "CGNS file requested for FFD but SU2 was built without CGNS support. No file written" << "\n"; 
#endif
}


void CFreeFormDefBox::SetTecplot(CGeometry *geometry, unsigned short iFFDBox, bool original) {
  
  ofstream FFDBox_file;
  char FFDBox_filename[MAX_STRING_SIZE];
  bool new_file;
  unsigned short iDim, iDegree, jDegree, kDegree;
  

  /*--- FFD output is always 3D (even in 2D problems),
   this is important for debuging ---*/

  nDim = 3;
  
  SPRINTF (FFDBox_filename, "ffd_boxes.dat");
  
  if ((original) && (iFFDBox == 0)) new_file = true;
  else new_file = false;
  
  if (new_file) {
    FFDBox_file.open(FFDBox_filename, ios::out);
    FFDBox_file << "TITLE = \"Visualization of the FFD boxes generated by SU2_DEF.\"" << endl;
    if (nDim == 2) FFDBox_file << "VARIABLES = \"x\", \"y\"" << endl;
    else FFDBox_file << "VARIABLES = \"x\", \"y\", \"z\"" << endl;
  }
  else FFDBox_file.open(FFDBox_filename, ios::out | ios::app);

  FFDBox_file << "ZONE T= \"" << Tag;
  if (original) FFDBox_file << " (Original FFD)\"";
  else FFDBox_file << " (Deformed FFD)\"";
  if (nDim == 2) FFDBox_file << ", I="<<lDegree+1<<", J="<<mDegree+1<<", DATAPACKING=POINT" << endl;
  else FFDBox_file << ", I="<<lDegree+1<<", J="<<mDegree+1<<", K="<<nDegree+1<<", DATAPACKING=POINT" << endl;

  FFDBox_file.precision(15);
  
  if (nDim == 2) {
    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
      for (iDegree = 0; iDegree <= lDegree; iDegree++) {
        for (iDim = 0; iDim < nDim; iDim++)
          FFDBox_file << scientific << Coord_Control_Points[iDegree][jDegree][0][iDim] << "\t";
        FFDBox_file << "\n";
      }
    }
  }
  else {
    for (kDegree = 0; kDegree <= nDegree; kDegree++) {
      for (jDegree = 0; jDegree <= mDegree; jDegree++) {
        for (iDegree = 0; iDegree <= lDegree; iDegree++) {
          for (iDim = 0; iDim < nDim; iDim++)
            FFDBox_file << scientific << Coord_Control_Points[iDegree][jDegree][kDegree][iDim] << "\t";
          FFDBox_file << "\n";
        }
      }
    }
  }
    
  FFDBox_file.close();
}

void CFreeFormDefBox::SetParaview(CGeometry *geometry, unsigned short iFFDBox, bool original) {
  
  ofstream FFDBox_file;
  char FFDBox_filename[MAX_STRING_SIZE];
  bool new_file;
  unsigned short iDim, iDegree, jDegree, kDegree;
  
  nDim = geometry->GetnDim();
  
  if (original) new_file = true;
  else new_file = false;

  if (new_file) SPRINTF (FFDBox_filename, "ffd_boxes_%d.vtk", SU2_TYPE::Int(iFFDBox));
  else SPRINTF (FFDBox_filename, "ffd_boxes_def_%d.vtk", SU2_TYPE::Int(iFFDBox));
  
  FFDBox_file.open(FFDBox_filename, ios::out);
  FFDBox_file << "# vtk DataFile Version 3.0" << endl;
  FFDBox_file << "vtk output" << endl;
  FFDBox_file << "ASCII" << endl;
  FFDBox_file << "DATASET STRUCTURED_GRID" << endl;
  
  if (nDim == 2) FFDBox_file << "DIMENSIONS "<<lDegree+1<<" "<<mDegree+1<<" "<<1<< endl;
  else FFDBox_file << "DIMENSIONS "<<lDegree+1<<" "<<mDegree+1<<" "<<nDegree+1<< endl;
  if (nDim == 2) FFDBox_file << "POINTS "<<(lDegree+1)*(mDegree+1)<<" float"<< endl;
  else FFDBox_file << "POINTS "<<(lDegree+1)*(mDegree+1)*(nDegree+1)<<" float"<< endl;

  FFDBox_file.precision(15);
  
  if (nDim == 2) {
    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
      for (iDegree = 0; iDegree <= lDegree; iDegree++) {
        for (iDim = 0; iDim < nDim; iDim++)
        FFDBox_file << scientific << Coord_Control_Points[iDegree][jDegree][0][iDim] << "\t";
        FFDBox_file << " 0.0 \n";
      }
    }
  }
  else {
    for (kDegree = 0; kDegree <= nDegree; kDegree++) {
      for (jDegree = 0; jDegree <= mDegree; jDegree++) {
        for (iDegree = 0; iDegree <= lDegree; iDegree++) {
          for (iDim = 0; iDim < nDim; iDim++)
          FFDBox_file << scientific << Coord_Control_Points[iDegree][jDegree][kDegree][iDim] << "\t";
          FFDBox_file << "\n";
        }
      }
    }
  }
    
  FFDBox_file.close();
  
}

su2double *CFreeFormDefBox::GetParametricCoord_Analytical(su2double *cart_coord) {
  unsigned short iDim;
  su2double *e1, *e2, *e3, *e12, *e23, *e13, *p;
  
  /*--- Auxiliary Basis Vectors of the deformed FFDBox ---*/
  e1 = new su2double[3]; e2 = new su2double[3]; e3 = new su2double[3];
  for (iDim = 0; iDim < nDim; iDim++) {
    e1[iDim] = Coord_Corner_Points[1][iDim]-Coord_Corner_Points[0][iDim];
    e2[iDim] = Coord_Corner_Points[3][iDim]-Coord_Corner_Points[0][iDim];
    e3[iDim] = Coord_Corner_Points[4][iDim]-Coord_Corner_Points[0][iDim];
  }
  
  /*--- Respective Cross-Products ---*/
  e12 = new su2double[3]; e23 = new su2double[3]; e13 = new su2double[3];
  CrossProduct(e1, e2, e12);
  CrossProduct(e1, e3, e13);
  CrossProduct(e2, e3, e23);
  
  /*--- p is Tranlated vector from the origin ---*/
  p = new su2double[3];
  for (iDim = 0; iDim < nDim; iDim++)
    p[iDim] = cart_coord[iDim] - Coord_Corner_Points[0][iDim];
  
  ParamCoord[0] = DotProduct(e23, p)/DotProduct(e23, e1);
  ParamCoord[1] = DotProduct(e13, p)/DotProduct(e13, e2);
  ParamCoord[2] = DotProduct(e12, p)/DotProduct(e12, e3);
  
  delete [] e1;
  delete [] e2;
  delete [] e3;
  delete [] e12;
  delete [] e23;
  delete [] e13;
  delete [] p;
  
  return ParamCoord;
}

su2double *CFreeFormDefBox::EvalCartesianCoord(su2double *ParamCoord) {
  unsigned short iDim, iDegree, jDegree, kDegree;
  
  for (iDim = 0; iDim < nDim; iDim++)
    cart_coord[iDim] = 0.0;
  
  for (iDegree = 0; iDegree <= lDegree; iDegree++)
    for (jDegree = 0; jDegree <= mDegree; jDegree++)
      for (kDegree = 0; kDegree <= nDegree; kDegree++)
        for (iDim = 0; iDim < nDim; iDim++) {
          cart_coord[iDim] += Coord_Control_Points[iDegree][jDegree][kDegree][iDim]
          * BlendingFunction[0]->GetBasis(iDegree, ParamCoord[0])
          * BlendingFunction[1]->GetBasis(jDegree, ParamCoord[1])
          * BlendingFunction[2]->GetBasis(kDegree, ParamCoord[2]);
        }
  
  return cart_coord;
}


su2double *CFreeFormDefBox::GetFFDGradient(su2double *val_coord, su2double *xyz) {
  
  unsigned short iDim, jDim, lmn[3];
  
  /*--- Set the Degree of the spline ---*/
  
  lmn[0] = lDegree; lmn[1] = mDegree; lmn[2] = nDegree;
  
  for (iDim = 0; iDim < nDim; iDim++) Gradient[iDim] = 0.0;
  
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      Gradient[jDim] += GetDerivative2(val_coord, iDim, xyz,  lmn) *
      GetDerivative3(val_coord, iDim, jDim, lmn);
  
  return Gradient;
  
}

void CFreeFormDefBox::GetFFDHessian(su2double *uvw, su2double *xyz, su2double **val_Hessian) {
  
  unsigned short iDim, jDim, lmn[3];
  
  /*--- Set the Degree of the spline ---*/
  
  lmn[0] = lDegree; lmn[1] = mDegree; lmn[2] = nDegree;
  
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      val_Hessian[iDim][jDim] = 0.0;
  
  /*--- Note that being all the functions linear combinations of polynomials, they are C^\infty,
   and the Hessian will be symmetric; no need to compute the under-diagonal part, for example ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    val_Hessian[0][0] += 2.0 * GetDerivative3(uvw, iDim,0, lmn) * GetDerivative3(uvw, iDim,0, lmn) +
    GetDerivative2(uvw, iDim,xyz, lmn) * GetDerivative5(uvw, iDim,0,0, lmn);
    
    val_Hessian[1][1] += 2.0 * GetDerivative3(uvw, iDim,1, lmn) * GetDerivative3(uvw, iDim,1, lmn) +
    GetDerivative2(uvw, iDim,xyz, lmn) * GetDerivative5(uvw, iDim,1,1, lmn);
    
    val_Hessian[2][2] += 2.0 * GetDerivative3(uvw, iDim,2, lmn) * GetDerivative3(uvw, iDim,2, lmn) +
    GetDerivative2(uvw, iDim,xyz, lmn) * GetDerivative5(uvw, iDim,2,2, lmn);
    
    val_Hessian[0][1] += 2.0 * GetDerivative3(uvw, iDim,0, lmn) * GetDerivative3(uvw, iDim,1, lmn) +
    GetDerivative2(uvw, iDim,xyz, lmn) * GetDerivative5(uvw, iDim,0,1, lmn);
    
    val_Hessian[0][2] += 2.0 * GetDerivative3(uvw, iDim,0, lmn) * GetDerivative3(uvw, iDim,2, lmn) +
    GetDerivative2(uvw, iDim,xyz, lmn) * GetDerivative5(uvw, iDim,0,2, lmn);
    
    val_Hessian[1][2] += 2.0 * GetDerivative3(uvw, iDim,1, lmn) * GetDerivative3(uvw, iDim,2, lmn) +
    GetDerivative2(uvw, iDim,xyz, lmn) * GetDerivative5(uvw, iDim,1,2, lmn);
  }
  
  val_Hessian[1][0] = val_Hessian[0][1];
  val_Hessian[2][0] = val_Hessian[0][2];
  val_Hessian[2][1] = val_Hessian[1][2];
  
}

su2double *CFreeFormDefBox::GetParametricCoord_Iterative(unsigned long iPoint, su2double *xyz, su2double *ParamCoordGuess, CConfig *config) {
  
  su2double *IndepTerm, SOR_Factor = 1.0, MinNormError, NormError, Determinant, AdjHessian[3][3], Temp[3] = {0.0,0.0,0.0};
  unsigned short iDim, jDim, RandonCounter;
  unsigned long iter;
  
  su2double tol = config->GetFFD_Tol()*1E-3;
  unsigned short it_max = config->GetnFFD_Iter();
  unsigned short Random_Trials = 500;

  /*--- Allocate the Hessian ---*/
  
  Hessian = new su2double* [nDim];
  IndepTerm = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Hessian[iDim] = new su2double[nDim];
    ParamCoord[iDim] = ParamCoordGuess[iDim];
    IndepTerm [iDim] = 0.0;
  }
  
  RandonCounter = 0; MinNormError = 1E6;
  
  /*--- External iteration ---*/

  for (iter = 0; iter < (unsigned long)it_max*Random_Trials; iter++) {
      
    /*--- The independent term of the solution of our system is -Gradient(sol_old) ---*/

    Gradient = GetFFDGradient(ParamCoord, xyz);
    
    for (iDim = 0; iDim < nDim; iDim++) IndepTerm[iDim] = - Gradient[iDim];

    /*--- Hessian = The Matrix of our system, getHessian(sol_old,xyz,...) ---*/
    
    GetFFDHessian(ParamCoord, xyz, Hessian);
    
    /*--- Adjoint to Hessian ---*/

    AdjHessian[0][0] = Hessian[1][1]*Hessian[2][2]-Hessian[1][2]*Hessian[2][1];
    AdjHessian[0][1] = Hessian[0][2]*Hessian[2][1]-Hessian[0][1]*Hessian[2][2];
    AdjHessian[0][2] = Hessian[0][1]*Hessian[1][2]-Hessian[0][2]*Hessian[1][1];
    AdjHessian[1][0] = Hessian[1][2]*Hessian[2][0]-Hessian[1][0]*Hessian[2][2];
    AdjHessian[1][1] = Hessian[0][0]*Hessian[2][2]-Hessian[0][2]*Hessian[2][0];
    AdjHessian[1][2] = Hessian[0][2]*Hessian[1][0]-Hessian[0][0]*Hessian[1][2];
    AdjHessian[2][0] = Hessian[1][0]*Hessian[2][1]-Hessian[1][1]*Hessian[2][0];
    AdjHessian[2][1] = Hessian[0][1]*Hessian[2][0]-Hessian[0][0]*Hessian[2][1];
    AdjHessian[2][2] = Hessian[0][0]*Hessian[1][1]-Hessian[0][1]*Hessian[1][0];
    
    /*--- Determinant of Hessian ---*/
    
    Determinant = Hessian[0][0]*AdjHessian[0][0]+Hessian[0][1]*AdjHessian[1][0]+Hessian[0][2]*AdjHessian[2][0];
    
    /*--- Hessian inverse ---*/
    
    if (Determinant != 0) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Temp[iDim] = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          Temp[iDim] += AdjHessian[iDim][jDim]*IndepTerm[jDim]/Determinant;
        }
      }
      for (iDim = 0; iDim < nDim; iDim++) {
        IndepTerm[iDim] = Temp[iDim];
      }
    }
    
    /*--- Update with Successive over-relaxation ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      ParamCoord[iDim] = (1.0-SOR_Factor)*ParamCoord[iDim] + SOR_Factor*(ParamCoord[iDim] + IndepTerm[iDim]);
    }

    /*--- If the gradient is small, we have converged ---*/
    
    if ((fabs(IndepTerm[0]) < tol) && (fabs(IndepTerm[1]) < tol) && (fabs(IndepTerm[2]) < tol))  break;

    /*--- Compute the norm of the error ---*/
    
    NormError = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      NormError += IndepTerm[iDim]*IndepTerm[iDim];
    NormError = sqrt(NormError);

    MinNormError = min(NormError, MinNormError);
      
    /*--- If we have no convergence with Random_Trials iterations probably we are in a local minima. ---*/
    
    if (((iter % it_max) == 0) && (iter != 0)) {
      
      RandonCounter++;
      if (RandonCounter == Random_Trials) {
        cout << endl << "Unknown point: "<< iPoint <<" (" << xyz[0] <<", "<< xyz[1] <<", "<< xyz[2] <<"). Min Error: "<< MinNormError <<". Iter: "<< iter <<"."<< endl;
      }
      else {
        SOR_Factor = 0.1;
        for (iDim = 0; iDim < nDim; iDim++)
          ParamCoord[iDim] = su2double(rand())/su2double(RAND_MAX);
      }

    }

    /* --- Splines are not defined outside of [0,1]. So if the parametric coords are outside of
     *  [0,1] the step was too big and we have to use a smaller relaxation factor. ---*/

    if ((config->GetFFD_Blending() == BSPLINE_UNIFORM)     &&
        (((ParamCoord[0] < 0.0) || (ParamCoord[0] > 1.0))  ||
         ((ParamCoord[1] < 0.0) || (ParamCoord[1] > 1.0))  ||
         ((ParamCoord[2] < 0.0) || (ParamCoord[2] > 1.0)))) {

      for (iDim = 0; iDim < nDim; iDim++){
        ParamCoord[iDim] = ParamCoordGuess[iDim];
      }
      SOR_Factor = 0.9*SOR_Factor;
    }

  }

  for (iDim = 0; iDim < nDim; iDim++)
    delete [] Hessian[iDim];
  delete [] Hessian;
  delete [] IndepTerm;

  /*--- The code has hit the max number of iterations ---*/

  if (iter == (unsigned long)it_max*Random_Trials) {
    cout << "Unknown point: (" << xyz[0] <<", "<< xyz[1] <<", "<< xyz[2] <<"). Increase the value of FFD_ITERATIONS." << endl;
  }
  
  /*--- Real Solution is now ParamCoord; Return it ---*/

  return ParamCoord;
  
}

bool CFreeFormDefBox::GetPointFFD(CGeometry *geometry, CConfig *config, unsigned long iPoint) {
  su2double Coord[3] = {0.0, 0.0, 0.0};
  unsigned short iVar, jVar, iDim;
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;

  bool Inside = false;
  bool cylindrical = (config->GetFFD_CoordSystem() == CYLINDRICAL);
  bool spherical = (config->GetFFD_CoordSystem() == SPHERICAL);
  bool polar = (config->GetFFD_CoordSystem() == POLAR);

  unsigned short Index[5][7] = {
    {0, 1, 2, 5, 0, 1, 2},
    {0, 2, 7, 5, 0, 2, 7},
    {0, 2, 3, 7, 0, 2, 3},
    {0, 5, 7, 4, 0, 5, 7},
    {2, 7, 5, 6, 2, 7, 5}};
  unsigned short nDim = geometry->GetnDim();
  
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = geometry->node[iPoint]->GetCoord(iDim);
  
  if (cylindrical) {
    
    X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
    
    Xbar =  Coord[0] - X_0; Ybar =  Coord[1] - Y_0; Zbar =  Coord[2] - Z_0;
    
    Coord[0] = sqrt(Ybar*Ybar + Zbar*Zbar);
    Coord[1] = atan2(Zbar, Ybar); if (Coord[1] > PI_NUMBER/2.0) Coord[1] -= 2.0*PI_NUMBER;
    Coord[2] = Xbar;
    
  }
  
  else if (spherical || polar) {
    
    X_0 = config->GetFFD_Axis(0); Y_0 = config->GetFFD_Axis(1);  Z_0 = config->GetFFD_Axis(2);
    
    Xbar =  Coord[0] - X_0; Ybar =  Coord[1] - Y_0; Zbar =  Coord[2] - Z_0;
    
    Coord[0] = sqrt(Xbar*Xbar + Ybar*Ybar + Zbar*Zbar);
    Coord[1] = atan2(Zbar, Ybar);  if (Coord[1] > PI_NUMBER/2.0) Coord[1] -= 2.0*PI_NUMBER;
    Coord[2] = acos(Xbar/Coord[0]);
    
  }
  
  /*--- 1st tetrahedron {V0, V1, V2, V5}
   2nd tetrahedron {V0, V2, V7, V5}
   3th tetrahedron {V0, V2, V3, V7}
   4th tetrahedron {V0, V5, V7, V4}
   5th tetrahedron {V2, V7, V5, V6} ---*/
  
  for (iVar = 0; iVar < 5; iVar++) {
    Inside = true;
    for (jVar = 0; jVar < 4; jVar++) {
      su2double Distance_Point = geometry->Point2Plane_Distance(Coord,
                                                                Coord_Corner_Points[Index[iVar][jVar+1]],
                                                                Coord_Corner_Points[Index[iVar][jVar+2]],
                                                                Coord_Corner_Points[Index[iVar][jVar+3]]);
      
      su2double Distance_Vertex = geometry->Point2Plane_Distance(Coord_Corner_Points[Index[iVar][jVar]],
                                                                 Coord_Corner_Points[Index[iVar][jVar+1]],
                                                                 Coord_Corner_Points[Index[iVar][jVar+2]],
                                                                 Coord_Corner_Points[Index[iVar][jVar+3]]);
      if (Distance_Point*Distance_Vertex < 0.0) Inside = false;          
    }
    if (Inside) break;
  }
  
  return Inside;
  
}

void CFreeFormDefBox::SetDeformationZone(CGeometry *geometry, CConfig *config, unsigned short iFFDBox) {
  su2double *Coord;
  unsigned short iMarker, iVar, jVar;
  unsigned long iVertex, iPoint;
  bool Inside = false;
  
  unsigned short Index[5][7] = {
    {0, 1, 2, 5, 0, 1, 2},
    {0, 2, 7, 5, 0, 2, 7},
    {0, 2, 3, 7, 0, 2, 3},
    {0, 5, 7, 4, 0, 5, 7},
    {2, 7, 5, 6, 2, 7, 5}};
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_DV(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {  
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        geometry->node[iPoint]->SetMove(false);
        
        Coord = geometry->node[iPoint]->GetCoord();
        
        /*--- 1st tetrahedron {V0, V1, V2, V5}
         2nd tetrahedron {V0, V2, V7, V5}
         3th tetrahedron {V0, V2, V3, V7}
         4th tetrahedron {V0, V5, V7, V4}
         5th tetrahedron {V2, V7, V5, V6} ---*/
        
        for (iVar = 0; iVar < 5; iVar++) {
          Inside = true;
          for (jVar = 0; jVar < 4; jVar++) {
            su2double Distance_Point = geometry->Point2Plane_Distance(Coord, 
                                                                   Coord_Corner_Points[Index[iVar][jVar+1]], 
                                                                   Coord_Corner_Points[Index[iVar][jVar+2]], 
                                                                   Coord_Corner_Points[Index[iVar][jVar+3]]);        
            su2double Distance_Vertex = geometry->Point2Plane_Distance(Coord_Corner_Points[Index[iVar][jVar]], 
                                                                    Coord_Corner_Points[Index[iVar][jVar+1]], 
                                                                    Coord_Corner_Points[Index[iVar][jVar+2]], 
                                                                    Coord_Corner_Points[Index[iVar][jVar+3]]);
            if (Distance_Point*Distance_Vertex < 0.0) Inside = false;          
          }
          if (Inside) break;
        }
        
        if (Inside) {
          geometry->node[iPoint]->SetMove(true);
        }
        
      }
}

su2double CFreeFormDefBox::GetDerivative1(su2double *uvw, unsigned short val_diff, unsigned short *ijk, unsigned short *lmn) {
  
  unsigned short iDim;
  su2double value = 0.0;
  
  value = BlendingFunction[val_diff]->GetDerivative(ijk[val_diff], uvw[val_diff], 1);
  for (iDim = 0; iDim < nDim; iDim++)
    if (iDim != val_diff)
      value *= BlendingFunction[iDim]->GetBasis(ijk[iDim], uvw[iDim]);
  
  return value;
  
}

su2double CFreeFormDefBox::GetDerivative2 (su2double *uvw, unsigned short dim, su2double *xyz, unsigned short *lmn) {
  
  unsigned short iDegree, jDegree, kDegree;
  su2double value = 0.0;
  
  for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
    for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
      for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
        value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] 
        * BlendingFunction[0]->GetBasis(iDegree, uvw[0])
        * BlendingFunction[1]->GetBasis(jDegree, uvw[1])
        * BlendingFunction[2]->GetBasis(kDegree, uvw[2]);
      }
   
  return 2.0*(value - xyz[dim]);
}

su2double CFreeFormDefBox::GetDerivative3(su2double *uvw, unsigned short dim, unsigned short diff_this, unsigned short *lmn) {
  
  unsigned short iDegree, jDegree, kDegree, iDim;
  su2double value = 0;
  
  unsigned short *ijk = new unsigned short[nDim];
  
  for (iDim = 0; iDim < nDim; iDim++) ijk[iDim] = 0;

  for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
    for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
      for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
        ijk[0] = iDegree; ijk[1] = jDegree; ijk[2] = kDegree;
        value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] * 
        GetDerivative1(uvw, diff_this, ijk, lmn);
      }
  
  delete [] ijk;

  return value;
}

su2double CFreeFormDefBox::GetDerivative4(su2double *uvw, unsigned short val_diff, unsigned short val_diff2,
                                       unsigned short *ijk, unsigned short *lmn) {
  unsigned short iDim;
  su2double value = 0.0;
  
  if (val_diff == val_diff2) {
    value = BlendingFunction[val_diff]->GetDerivative(ijk[val_diff], uvw[val_diff], 2);
    for (iDim = 0; iDim < nDim; iDim++)
      if (iDim != val_diff)
        value *= BlendingFunction[iDim]->GetBasis(ijk[iDim], uvw[iDim]);
  }
  else {
    value = BlendingFunction[val_diff]->GetDerivative(ijk[val_diff],  uvw[val_diff],1) *
    BlendingFunction[val_diff2]->GetDerivative(ijk[val_diff2], uvw[val_diff2], 1);
    for (iDim = 0; iDim < nDim; iDim++)
      if ((iDim != val_diff) && (iDim != val_diff2))
        value *= BlendingFunction[iDim]->GetBasis(ijk[iDim], uvw[iDim]);
  }
  
  return value;
}

su2double CFreeFormDefBox::GetDerivative5(su2double *uvw, unsigned short dim, unsigned short diff_this, unsigned short diff_this_also, 
                                      unsigned short *lmn) {
  
  unsigned short iDegree, jDegree, kDegree, iDim;
  su2double value = 0.0;
  
  unsigned short *ijk = new unsigned short[nDim];
  
  for (iDim = 0; iDim < nDim; iDim++) ijk[iDim] = 0;
  
  for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
    for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
      for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
        ijk[0] = iDegree; ijk[1] = jDegree; ijk[2] = kDegree;
        value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] *
        GetDerivative4(uvw, diff_this, diff_this_also, ijk, lmn);
      }
  
  delete [] ijk;
  
  return value;
}



CElasticityMovement::CElasticityMovement(CGeometry *geometry, CConfig *config) : CVolumetricMovement(), System(true) {
  
    size = SU2_MPI::GetSize();
    rank = SU2_MPI::GetRank();

    /*--- Initialize the number of spatial dimensions, length of the state
     vector (same as spatial dimensions for grid deformation), and grid nodes. ---*/

    unsigned short iDim, jDim;

    nDim         = geometry->GetnDim();
    nVar         = geometry->GetnDim();
    nPoint       = geometry->GetnPoint();
    nPointDomain = geometry->GetnPointDomain();

    nIterMesh   = 0;
    valResidual = 0.0;

    MinVolume = 0.0;
    MaxVolume = 0.0;

    Residual = new su2double[nDim];   for (iDim = 0; iDim < nDim; iDim++) Residual[iDim] = 0.0;
    Solution = new su2double[nDim];   for (iDim = 0; iDim < nDim; iDim++) Solution[iDim] = 0.0;

    /*--- Initialize matrix, solution, and r.h.s. structures for the linear solver. ---*/

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    StiffMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);

    /*--- Matrices to impose boundary conditions ---*/

    matrixZeros = new su2double *[nDim];
    matrixId    = new su2double *[nDim];
    for(iDim = 0; iDim < nDim; iDim++){
      matrixZeros[iDim] = new su2double[nDim];
      matrixId[iDim]    = new su2double[nDim];
    }

    for(iDim = 0; iDim < nDim; iDim++){
      for (jDim = 0; jDim < nDim; jDim++){
        matrixZeros[iDim][jDim] = 0.0;
        matrixId[iDim][jDim]    = 0.0;
      }
      matrixId[iDim][iDim] = 1.0;
    }

    /*--- Structural parameters ---*/

    E      = config->GetDeform_ElasticityMod();
    Nu     = config->GetDeform_PoissonRatio();

    Mu     = E / (2.0*(1.0 + Nu));
    Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));

    /*--- Element container structure ---*/

    element_container = new CElement* [MAX_FE_KINDS];
    for (unsigned short iKind = 0; iKind < MAX_FE_KINDS; iKind++) {
      element_container[iKind] = NULL;
    }
    if (nDim == 2){
      element_container[EL_TRIA] = new CTRIA1();
      element_container[EL_QUAD] = new CQUAD4();
    }
    else {
      element_container[EL_TETRA] = new CTETRA1();
      element_container[EL_HEXA]  = new CHEXA8();
      element_container[EL_PYRAM] = new CPYRAM5();
      element_container[EL_PRISM] = new CPRISM6();
    }

    /*--- Term ij of the Jacobian ---*/

    Jacobian_ij = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Jacobian_ij[iDim] = new su2double [nDim];
      for (jDim = 0; jDim < nDim; jDim++) {
        Jacobian_ij[iDim][jDim] = 0.0;
      }
    }

    KAux_ab = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      KAux_ab[iDim] = new su2double[nDim];
    }

    unsigned short iVar;

    if (nDim == 2){
      Ba_Mat  = new su2double* [3];
      Bb_Mat  = new su2double* [3];
      D_Mat   = new su2double* [3];
      GradNi_Ref_Mat = new su2double* [4];
      for (iVar = 0; iVar < 3; iVar++) {
        Ba_Mat[iVar]    = new su2double[nDim];
        Bb_Mat[iVar]    = new su2double[nDim];
        D_Mat[iVar]     = new su2double[3];
      }
      for (iVar = 0; iVar < 4; iVar++) {
        GradNi_Ref_Mat[iVar]  = new su2double[nDim];
      }
    }
    else if (nDim == 3){
      Ba_Mat  = new su2double* [6];
      Bb_Mat  = new su2double* [6];
      D_Mat   = new su2double* [6];
      GradNi_Ref_Mat = new su2double* [8];
      for (iVar = 0; iVar < 6; iVar++) {
        Ba_Mat[iVar]      = new su2double[nDim];
        Bb_Mat[iVar]      = new su2double[nDim];
        D_Mat[iVar]       = new su2double[6];
      }
      for (iVar = 0; iVar < 8; iVar++) {
        GradNi_Ref_Mat[iVar]  = new su2double[nDim];
      }
    }



}

CElasticityMovement::~CElasticityMovement(void) {

  unsigned short iDim, iVar;

  delete [] Residual;
  delete [] Solution;

  for (iDim = 0; iDim < nDim; iDim++) {
    delete [] matrixZeros[iDim];
    delete [] matrixId[iDim];
    delete [] Jacobian_ij[iDim];
    delete [] KAux_ab[iDim];
  }
  delete [] matrixZeros;
  delete [] matrixId;
  delete [] Jacobian_ij;
  delete [] KAux_ab;

  if (nDim == 2){
    for (iVar = 0; iVar < 3; iVar++){
      delete [] Ba_Mat[iVar];
      delete [] Bb_Mat[iVar];
      delete [] D_Mat[iVar];
    }
    for (iVar = 0; iVar < 4; iVar++){
      delete [] GradNi_Ref_Mat[iVar];
    }
  }
  else if (nDim == 3){
    for (iVar = 0; iVar < 6; iVar++){
      delete [] Ba_Mat[iVar];
      delete [] Bb_Mat[iVar];
      delete [] D_Mat[iVar];
    }
    for (iVar = 0; iVar < 8; iVar++){
      delete [] GradNi_Ref_Mat[iVar];
    }
  }

  delete [] Ba_Mat;
  delete [] Bb_Mat;
  delete [] D_Mat;
  delete [] GradNi_Ref_Mat;

  if (element_container != NULL) {
    for (iVar = 0; iVar < MAX_FE_KINDS; iVar++){
      if (element_container[iVar] != NULL) delete element_container[iVar];
    }
    delete [] element_container;
  }

}


void CElasticityMovement::SetVolume_Deformation_Elas(CGeometry *geometry, CConfig *config, bool UpdateGeo, bool screen_output, bool Derivative){

  unsigned long iNonlinear_Iter, Nonlinear_Iter = 0;

  bool discrete_adjoint = config->GetDiscrete_Adjoint();

  /*--- Retrieve number or internal iterations from config ---*/

  Nonlinear_Iter = config->GetGridDef_Nonlinear_Iter();

  /*--- Loop over the total number of grid deformation iterations. The surface
   deformation can be divided into increments to help with stability. ---*/

  for (iNonlinear_Iter = 0; iNonlinear_Iter < Nonlinear_Iter; iNonlinear_Iter++) {

    /*--- Initialize vector and sparse matrix ---*/
    LinSysSol.SetValZero();
    LinSysRes.SetValZero();
    StiffMatrix.SetValZero();
    
    if ((rank == MASTER_NODE) && (!discrete_adjoint) && screen_output)
      cout << "Computing volumes of the grid elements." << endl;

    /*--- Compute the minimum and maximum area/volume for the mesh. ---*/
    SetMinMaxVolume(geometry, config);
    if ((rank == MASTER_NODE) && (!discrete_adjoint) && screen_output) {
      if (nDim == 2) cout << scientific << "Min. area: "<< MinVolume <<", max. area: " << MaxVolume <<"." << endl;
      else           cout << scientific << "Min. volume: "<< MinVolume <<", max. volume: " << MaxVolume <<"." << endl;
    }

    /*--- Compute the stiffness matrix. ---*/
    SetStiffnessMatrix(geometry, config);

    /*--- Impose boundary conditions (all of them are ESSENTIAL BC's - displacements). ---*/
    SetBoundaryDisplacements(geometry, config);

    /*--- Solve the linear system. ---*/

#ifdef CODI_REVERSE_TYPE
    /*--- We need to guard the SendReceive_Solution otherwise the FSI adjoint breaks. ---*/
    bool TapeActive = NO;
    if (config->GetDiscrete_Adjoint()) {
      TapeActive = AD::globalTape.isActive();
      AD::StopRecording();
    }
#endif
    StiffMatrix.InitiateComms(LinSysSol, geometry, config, SOLUTION_MATRIX);
    StiffMatrix.CompleteComms(LinSysSol, geometry, config, SOLUTION_MATRIX);
    
    StiffMatrix.InitiateComms(LinSysRes, geometry, config, SOLUTION_MATRIX);
    StiffMatrix.CompleteComms(LinSysRes, geometry, config, SOLUTION_MATRIX);
#ifdef CODI_REVERSE_TYPE
    if (TapeActive) AD::StartRecording();
#endif
    nIterMesh = System.Solve(StiffMatrix, LinSysRes, LinSysSol, geometry, config);
    valResidual = System.GetResidual();

    /*--- Update the grid coordinates and cell volumes using the solution
     of the linear system (usol contains the x, y, z displacements). ---*/
    UpdateGridCoord(geometry, config);

    if (UpdateGeo)
      UpdateDualGrid(geometry, config);

    /*--- Check for failed deformation (negative volumes). ---*/
    /*--- In order to do this, we recompute the minimum and maximum area/volume for the mesh. ---*/
    SetMinMaxVolume(geometry, config);

    if ((rank == MASTER_NODE) && (!discrete_adjoint) && screen_output) {
      cout << scientific << "Non-linear iter.: " << iNonlinear_Iter+1 << "/" << Nonlinear_Iter  << ". Linear iter.: " << nIterMesh << ". ";
      if (nDim == 2) cout << "Min. area: " << MinVolume << ". Error: " << valResidual << "." << endl;
      else cout << "Min. volume: " << MinVolume << ". Error: " << valResidual << "." << endl;
    }

  }

}


void CElasticityMovement::UpdateGridCoord(CGeometry *geometry, CConfig *config){

  unsigned short iDim;
  unsigned long iPoint, total_index;
  su2double new_coord;

  /*--- Update the grid coordinates using the solution of the linear system
   after grid deformation (LinSysSol contains the x, y, z displacements). ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iDim = 0; iDim < nDim; iDim++) {
      total_index = iPoint*nDim + iDim;
      new_coord = geometry->node[iPoint]->GetCoord(iDim)+LinSysSol[total_index];
      if (fabs(new_coord) < EPS*EPS) new_coord = 0.0;
      geometry->node[iPoint]->SetCoord(iDim, new_coord);
    }

  /* --- LinSysSol contains the non-transformed displacements in the periodic halo cells.
   * Hence we still need a communication of the transformed coordinates, otherwise periodicity
   * is not maintained. ---*/

  geometry->InitiateComms(geometry, config, COORDINATES);
  geometry->CompleteComms(geometry, config, COORDINATES);
  
}


void CElasticityMovement::UpdateDualGrid(CGeometry *geometry, CConfig *config){

  /*--- After moving all nodes, update the dual mesh. Recompute the edges and
   dual mesh control volumes in the domain and on the boundaries. ---*/

  geometry->SetCoord_CG();
  geometry->SetControlVolume(config, UPDATE);
  geometry->SetBoundControlVolume(config, UPDATE);
  geometry->SetMaxLength(config);

}


void CElasticityMovement::UpdateMultiGrid(CGeometry **geometry, CConfig *config){

  unsigned short iMGfine, iMGlevel, nMGlevel = config->GetnMGLevels();

  /*--- Update the multigrid structure after moving the finest grid,
   including computing the grid velocities on the coarser levels. ---*/

  for (iMGlevel = 1; iMGlevel <= nMGlevel; iMGlevel++) {
    iMGfine = iMGlevel-1;
    geometry[iMGlevel]->SetControlVolume(config, geometry[iMGfine], UPDATE);
    geometry[iMGlevel]->SetBoundControlVolume(config, geometry[iMGfine],UPDATE);
    geometry[iMGlevel]->SetCoord(geometry[iMGfine]);
    if (config->GetGrid_Movement())
      geometry[iMGlevel]->SetRestricted_GridVelocity(geometry[iMGfine], config);
  }

}

void CElasticityMovement::SetBoundaryDisplacements(CGeometry *geometry, CConfig *config){

  /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
   deforming meshes (MARKER_FSI_INTERFACE). ---*/

  unsigned short Kind_SU2 = config->GetKind_SU2();

  /*--- Share halo vertex displacements across ranks using the solution vector.
   The transfer routines do not do this, i.e. halo vertices have wrong values. ---*/

  unsigned short iMarker, iDim;
  unsigned long iNode, iVertex;

  su2double VarIncrement = 1.0/su2double(config->GetGridDef_Nonlinear_Iter());

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_ZoneInterface(iMarker) != 0 || 
         config->GetMarker_All_Moving(iMarker)) 
        && (Kind_SU2 == SU2_CFD)) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        /*--- Get node index ---*/
        iNode = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->node[iNode]->GetDomain()) {

          /*--- Get the displacement on the vertex ---*/
          for (iDim = 0; iDim < nDim; iDim++)
             Solution[iDim] = geometry->vertex[iMarker][iVertex]->GetVarCoord()[iDim] * VarIncrement;

          /*--- Initialize the solution vector ---*/
          LinSysSol.SetBlock(iNode, Solution);
        }
      }
    }
  }
  StiffMatrix.InitiateComms(LinSysSol, geometry, config, SOLUTION_MATRIX);
  StiffMatrix.CompleteComms(LinSysSol, geometry, config, SOLUTION_MATRIX);

  /*--- Apply displacement boundary conditions to the FSI interfaces. ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_ZoneInterface(iMarker) != 0 || 
         config->GetMarker_All_Moving(iMarker)) 
        && (Kind_SU2 == SU2_CFD)) {
      SetMoving_Boundary(geometry, config, iMarker);
    }
  }

  /*--- Now, set to zero displacements of all the other boundary conditions, except the symmetry
   plane, the receive boundaries and periodic boundaries. ---*/


  /*--- As initialization, set to zero displacements of all the surfaces except the symmetry
   plane, the receive boundaries and periodic boundaries. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (((config->GetMarker_All_KindBC(iMarker) != SYMMETRY_PLANE) &&
         (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) &&
         (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY))
        && !config->GetMarker_All_Moving(iMarker)) {

      /*--- We must note that the FSI surfaces are not clamped ---*/
      if (config->GetMarker_All_ZoneInterface(iMarker) == 0){
        SetClamped_Boundary(geometry, config, iMarker);
      }
    }
  }

  /*--- All others are pending. ---*/

}


void CElasticityMovement::SetClamped_Boundary(CGeometry *geometry, CConfig *config, unsigned short val_marker){

  unsigned long iNode, iVertex;

  Solution[0] = 0.0;
  Solution[1] = 0.0;
  if (nDim==3) Solution[2] = 0.0;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/
    iNode = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Set and enforce solution ---*/
    LinSysSol.SetBlock(iNode, Solution);
    StiffMatrix.EnforceSolutionAtNode(iNode, Solution, LinSysRes);

  }

}

void CElasticityMovement::SetMoving_Boundary(CGeometry *geometry, CConfig *config, unsigned short val_marker){

  unsigned long iNode, iVertex;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/
    iNode = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Enforce solution ---*/
    StiffMatrix.EnforceSolutionAtNode(iNode, LinSysSol.GetBlock(iNode), LinSysRes);

  }

}

void CElasticityMovement::SetMinMaxVolume(CGeometry *geometry, CConfig *config) {

  unsigned long iElem, ElemCounter = 0;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;

  bool RightVol = true;

  su2double ElemVolume;

  MaxVolume = -1E22; MinVolume = 1E22;

  /*--- Loops over all the elements ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       {nNodes = 5; EL_KIND = EL_PYRAM;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         {nNodes = 6; EL_KIND = EL_PRISM;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    {nNodes = 8; EL_KIND = EL_HEXA;}

    /*--- For the number of nodes, we get the coordinates from the connectivity matrix and the geometry structure ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {

      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        element_container[EL_KIND]->SetRef_Coord(iNode, iDim, val_Coord);
      }

    }

    /*--- Compute the volume of the element (or the area in 2D cases ) ---*/

    if (nDim == 2)  ElemVolume = element_container[EL_KIND]->ComputeArea();
    else            ElemVolume = element_container[EL_KIND]->ComputeVolume();

    RightVol = true;
    if (ElemVolume < 0.0) RightVol = false;

    MaxVolume = max(MaxVolume, ElemVolume);
    MinVolume = min(MinVolume, ElemVolume);
    geometry->elem[iElem]->SetVolume(ElemVolume);

    if (!RightVol) ElemCounter++;

  }

#ifdef HAVE_MPI
  unsigned long ElemCounter_Local = ElemCounter; ElemCounter = 0;
  su2double MaxVolume_Local = MaxVolume; MaxVolume = 0.0;
  su2double MinVolume_Local = MinVolume; MinVolume = 0.0;
  SU2_MPI::Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MaxVolume_Local, &MaxVolume, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MinVolume_Local, &MinVolume, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

  /*--- Volume from  0 to 1 ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    ElemVolume = geometry->elem[iElem]->GetVolume()/MaxVolume;
    geometry->elem[iElem]->SetVolume(ElemVolume);
  }

  if ((ElemCounter != 0) && (rank == MASTER_NODE))
    cout <<"There are " << ElemCounter << " elements with negative volume.\n" << endl;

}


void CElasticityMovement::SetStiffnessMatrix(CGeometry *geometry, CConfig *config){

  unsigned long iElem;
  unsigned short iNode, iDim, jDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;

  const su2double *Kab = NULL;
  unsigned short NelNodes, jNode;

  su2double ElemVolume;

  /*--- Loops over all the elements ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       {nNodes = 5; EL_KIND = EL_PYRAM;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         {nNodes = 6; EL_KIND = EL_PRISM;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    {nNodes = 8; EL_KIND = EL_HEXA;}

    /*--- For the number of nodes, we get the coordinates from the connectivity matrix and the geometry structure ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {

      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        element_container[EL_KIND]->SetRef_Coord(iNode, iDim, val_Coord);
      }

    }

    /*--- Retrieve the volume of the element (previously computed in SetMinMaxVolume()) ---*/

    ElemVolume = geometry->elem[iElem]->GetVolume();

    /*--- Compute the stiffness of the element ---*/
    Set_Element_Stiffness(ElemVolume, config);

    /*--- Compute the element contribution to the stiffness matrix ---*/

    Compute_Element_Contribution(element_container[EL_KIND], config);

    /*--- Retrieve number of nodes ---*/

    NelNodes = element_container[EL_KIND]->GetnNodes();

    /*--- Assemble the stiffness matrix ---*/

    for (iNode = 0; iNode < NelNodes; iNode++){

      for (jNode = 0; jNode < NelNodes; jNode++){

        Kab = element_container[EL_KIND]->Get_Kab(iNode, jNode);

        for (iDim = 0; iDim < nDim; iDim++){
          for (jDim = 0; jDim < nDim; jDim++){
            Jacobian_ij[iDim][jDim] = Kab[iDim*nDim+jDim];
          }
        }

        StiffMatrix.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);

      }

    }

  }

}


void CElasticityMovement::Set_Element_Stiffness(su2double ElemVolume, CConfig *config) {

  switch (config->GetDeform_Stiffness_Type()) {
    case INVERSE_VOLUME:
      E = 1.0 / ElemVolume;           // Stiffness inverse of the volume of the element
      Nu = config->GetDeform_Coeff(); // Nu is normally a very large number, for rigid-body rotations, see Dwight (2009)
      break;
    case SOLID_WALL_DISTANCE:
      SU2_MPI::Error("SOLID_WALL DISTANCE METHOD NOT YET IMPLEMENTED FOR THIS APPROACH!!", CURRENT_FUNCTION);
      break;
    case CONSTANT_STIFFNESS:
      E      = config->GetDeform_ElasticityMod();
      Nu     = config->GetDeform_PoissonRatio();
      break;
  }

  /*--- Lamé parameters ---*/

  Mu     = E / (2.0*(1.0 + Nu));
  Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));

}


void CElasticityMovement::Compute_Element_Contribution(CElement *element, CConfig *config){

  unsigned short iVar, jVar, kVar;
  unsigned short iGauss, nGauss;
  unsigned short iNode, jNode, nNode;
  unsigned short iDim;
  unsigned short bDim;

  su2double Weight, Jac_X;

  su2double AuxMatrix[3][6];

  /*--- Initialize auxiliary matrices ---*/

  if (nDim == 2) bDim = 3;
  else bDim = 6;

  for (iVar = 0; iVar < bDim; iVar++){
    for (jVar = 0; jVar < nDim; jVar++){
      Ba_Mat[iVar][jVar] = 0.0;
      Bb_Mat[iVar][jVar] = 0.0;
    }
  }

  for (iVar = 0; iVar < 3; iVar++){
    for (jVar = 0; jVar < 6; jVar++){
      AuxMatrix[iVar][jVar] = 0.0;
    }
  }

  element->ClearElement();      /*--- Restarts the element: avoids adding over previous results in other elements --*/
  element->ComputeGrad_Linear();
  nNode = element->GetnNodes();
  nGauss = element->GetnGaussPoints();

  /*--- Compute the constitutive matrix (D_Mat) for the element - it only depends on lambda and mu, constant for the element ---*/
  Compute_Constitutive_Matrix();

  for (iGauss = 0; iGauss < nGauss; iGauss++){

    Weight = element->GetWeight(iGauss);
    Jac_X = element->GetJ_X(iGauss);

    /*--- Retrieve the values of the gradients of the shape functions for each node ---*/
    /*--- This avoids repeated operations ---*/
    for (iNode = 0; iNode < nNode; iNode++){
      for (iDim = 0; iDim < nDim; iDim++){
        GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
      }
    }

    for (iNode = 0; iNode < nNode; iNode++){

      if (nDim == 2){
        Ba_Mat[0][0] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[1][1] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[2][0] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[2][1] = GradNi_Ref_Mat[iNode][0];
      }
      else if (nDim == 3){
        Ba_Mat[0][0] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[1][1] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[2][2] = GradNi_Ref_Mat[iNode][2];
        Ba_Mat[3][0] = GradNi_Ref_Mat[iNode][1];
        Ba_Mat[3][1] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[4][0] = GradNi_Ref_Mat[iNode][2];
        Ba_Mat[4][2] = GradNi_Ref_Mat[iNode][0];
        Ba_Mat[5][1] = GradNi_Ref_Mat[iNode][2];
        Ba_Mat[5][2] = GradNi_Ref_Mat[iNode][1];
      }

        /*--- Compute the BT.D Matrix ---*/

      for (iVar = 0; iVar < nDim; iVar++){
        for (jVar = 0; jVar < bDim; jVar++){
          AuxMatrix[iVar][jVar] = 0.0;
          for (kVar = 0; kVar < bDim; kVar++){
            AuxMatrix[iVar][jVar] += Ba_Mat[kVar][iVar]*D_Mat[kVar][jVar];
          }
        }
      }

      /*--- Assumming symmetry ---*/
      for (jNode = iNode; jNode < nNode; jNode++){
        if (nDim == 2){
          Bb_Mat[0][0] = GradNi_Ref_Mat[jNode][0];
          Bb_Mat[1][1] = GradNi_Ref_Mat[jNode][1];
          Bb_Mat[2][0] = GradNi_Ref_Mat[jNode][1];
          Bb_Mat[2][1] = GradNi_Ref_Mat[jNode][0];
        }
        else if (nDim ==3){
          Bb_Mat[0][0] = GradNi_Ref_Mat[jNode][0];
          Bb_Mat[1][1] = GradNi_Ref_Mat[jNode][1];
          Bb_Mat[2][2] = GradNi_Ref_Mat[jNode][2];
          Bb_Mat[3][0] = GradNi_Ref_Mat[jNode][1];
          Bb_Mat[3][1] = GradNi_Ref_Mat[jNode][0];
          Bb_Mat[4][0] = GradNi_Ref_Mat[jNode][2];
          Bb_Mat[4][2] = GradNi_Ref_Mat[jNode][0];
          Bb_Mat[5][1] = GradNi_Ref_Mat[jNode][2];
          Bb_Mat[5][2] = GradNi_Ref_Mat[jNode][1];
        }

        for (iVar = 0; iVar < nDim; iVar++){
          for (jVar = 0; jVar < nDim; jVar++){
            KAux_ab[iVar][jVar] = 0.0;
            for (kVar = 0; kVar < bDim; kVar++){
              KAux_ab[iVar][jVar] += Weight * AuxMatrix[iVar][kVar] * Bb_Mat[kVar][jVar] * Jac_X;
            }
          }
        }

        element->Add_Kab(iNode, jNode, KAux_ab);
        /*--- Symmetric terms --*/
        if (iNode != jNode){
          element->Add_Kab_T(jNode, iNode, KAux_ab);
        }

      }

    }

  }

}

void CElasticityMovement::Compute_Constitutive_Matrix(void){

  /*--- Compute the D Matrix (for plane strain and 3-D)---*/

  if (nDim == 2){

    /*--- Assuming plane strain ---*/
    D_Mat[0][0] = Lambda + 2.0*Mu;  D_Mat[0][1] = Lambda;            D_Mat[0][2] = 0.0;
    D_Mat[1][0] = Lambda;           D_Mat[1][1] = Lambda + 2.0*Mu;   D_Mat[1][2] = 0.0;
    D_Mat[2][0] = 0.0;              D_Mat[2][1] = 0.0;               D_Mat[2][2] = Mu;

  }
  else if (nDim == 3){

    D_Mat[0][0] = Lambda + 2.0*Mu;  D_Mat[0][1] = Lambda;           D_Mat[0][2] = Lambda;           D_Mat[0][3] = 0.0;  D_Mat[0][4] = 0.0;  D_Mat[0][5] = 0.0;
    D_Mat[1][0] = Lambda;           D_Mat[1][1] = Lambda + 2.0*Mu;  D_Mat[1][2] = Lambda;           D_Mat[1][3] = 0.0;  D_Mat[1][4] = 0.0;  D_Mat[1][5] = 0.0;
    D_Mat[2][0] = Lambda;           D_Mat[2][1] = Lambda;           D_Mat[2][2] = Lambda + 2.0*Mu;  D_Mat[2][3] = 0.0;  D_Mat[2][4] = 0.0;  D_Mat[2][5] = 0.0;
    D_Mat[3][0] = 0.0;              D_Mat[3][1] = 0.0;              D_Mat[3][2] = 0.0;              D_Mat[3][3] = Mu;   D_Mat[3][4] = 0.0;  D_Mat[3][5] = 0.0;
    D_Mat[4][0] = 0.0;              D_Mat[4][1] = 0.0;              D_Mat[4][2] = 0.0;              D_Mat[4][3] = 0.0;  D_Mat[4][4] = Mu;   D_Mat[4][5] = 0.0;
    D_Mat[5][0] = 0.0;              D_Mat[5][1] = 0.0;              D_Mat[5][2] = 0.0;              D_Mat[5][3] = 0.0;  D_Mat[5][4] = 0.0;  D_Mat[5][5] = Mu;

  }

}

void CElasticityMovement::Transfer_Boundary_Displacements(CGeometry *geometry, CConfig *config, unsigned short val_marker){

  unsigned short iDim;
  unsigned long iNode, iVertex;

  su2double *VarCoord;
  su2double new_coord;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/

    iNode = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Get the displacement on the vertex ---*/

    VarCoord = geometry->vertex[val_marker][iVertex]->GetVarCoord();

    if (geometry->node[iNode]->GetDomain()) {

      /*--- Update the grid coordinates using the solution of the structural problem
       *--- recorded in VarCoord. */

      for (iDim = 0; iDim < nDim; iDim++) {
        new_coord = geometry->node[iNode]->GetCoord(iDim)+VarCoord[iDim];
        if (fabs(new_coord) < EPS*EPS) new_coord = 0.0;
        geometry->node[iNode]->SetCoord(iDim, new_coord);
      }
    }

  }

}

void CElasticityMovement::Boundary_Dependencies(CGeometry **geometry, CConfig *config){


  unsigned short iMarker;

  /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
   deforming meshes (MARKER_FSI_INTERFACE). ---*/

  unsigned short Kind_SU2 = config->GetKind_SU2();

  /*--- Set the dependencies on the FSI interfaces. ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_ZoneInterface(iMarker) != 0) && (Kind_SU2 == SU2_CFD)) {
      Transfer_Boundary_Displacements(geometry[MESH_0], config, iMarker);
    }
  }

  UpdateDualGrid(geometry[MESH_0], config);

}


CFreeFormBlending::CFreeFormBlending(){}

CFreeFormBlending::~CFreeFormBlending(){}

CBSplineBlending::CBSplineBlending(short val_order, short n_controlpoints): CFreeFormBlending(){
  SetOrder(val_order, n_controlpoints);
}

CBSplineBlending::~CBSplineBlending(){}

void CBSplineBlending::SetOrder(short val_order, short n_controlpoints){

  unsigned short iKnot;

  Order    = val_order;
  Degree   = Order - 1;
  nControl = n_controlpoints;

  KnotSize = Order + nControl;

  U.resize(KnotSize, 0.0);

  /*--- Set up the knot vectors for open uniform B-Splines ---*/

  /*--- Note: the first knots are zero now.---*/

  /*--- The next knots are equidistantly distributed in [0,1] ---*/

  for (iKnot = 0; iKnot < nControl - Order; iKnot++){
    U[Order + iKnot] = su2double(iKnot+1)/su2double(nControl - Order + 1);
  }
  for (iKnot = nControl - Order; iKnot < nControl; iKnot++){
    U[Order + iKnot]  = 1.0;
  }

  /*--- Allocate the temporary vectors for the basis evaluation ---*/

  N.resize(Order, vector<su2double>(Order, 0.0));
}

su2double CBSplineBlending::GetBasis(short val_i, su2double val_t){

  /*--- Evaluation is based on the algorithm from "The NURBS Book (Les Piegl and Wayne Tiller)" ---*/

  /*--- Special cases ---*/

  if ((val_i == 0 && val_t == U[0]) || (val_i == (short)U.size()-1 && val_t == U.back())) {return 1.0;}

  /*--- Local property of BSplines ---*/

  if ((val_t < U[val_i]) || (val_t >= U[val_i+Order])){ return 0.0;}

  unsigned short j,k;
  su2double saved, temp;

  for (j = 0; j < Order; j++){
    if ((val_t >= U[val_i+j]) && (val_t < U[val_i+j+1])) N[j][0] = 1.0;
    else N[j][0] = 0;
  }

  for (k = 1; k < Order; k++){
    if (N[0][k-1] == 0.0) saved = 0.0;
    else saved = ((val_t - U[val_i])*N[0][k-1])/(U[val_i+k] - U[val_i]);
    for (j = 0; j < Order-k; j++){
      if (N[j+1][k-1] == 0.0){
        N[j][k] = saved; saved = 0.0;
      } else {
        temp     = N[j+1][k-1]/(U[val_i+j+k+1] - U[val_i+j+1]);
        N[j][k]  = saved+(U[val_i+j+k+1] - val_t)*temp;
        saved    = (val_t - U[val_i+j+1])*temp;
      }
    }
  }
  return N[0][Order-1];
}

su2double CBSplineBlending::GetDerivative(short val_i, su2double val_t, short val_order_der){

  if ((val_t < U[val_i]) || (val_t >= U[val_i+Order])){ return 0.0;}

  /*--- Evaluate the i+p basis functions up to the order p (stored in the matrix N). ---*/

  GetBasis(val_i, val_t);

  /*--- Use the recursive definition for the derivative (hardcoded for 1st and 2nd derivative). ---*/

  if (val_order_der == 0){ return N[0][Order-1];}

  if (val_order_der == 1){
    return (Order-1.0)/(1e-10 + U[val_i+Order-1] - U[val_i]  )*N[0][Order-2]
         - (Order-1.0)/(1e-10 + U[val_i+Order]   - U[val_i+1])*N[1][Order-2];
  }

  if (val_order_der == 2 && Order > 2){
    const su2double left = (Order-2.0)/(1e-10 + U[val_i+Order-2] - U[val_i])  *N[0][Order-3]
                         - (Order-2.0)/(1e-10 + U[val_i+Order-1] - U[val_i+1])*N[1][Order-3];

    const su2double right = (Order-2.0)/(1e-10 + U[val_i+Order-1] - U[val_i+1])*N[1][Order-3]
                          - (Order-2.0)/(1e-10 + U[val_i+Order]   - U[val_i+2])*N[2][Order-3];

    return (Order-1.0)/(1e-10 + U[val_i+Order-1] - U[val_i]  )*left
         - (Order-1.0)/(1e-10 + U[val_i+Order]   - U[val_i+1])*right;
  }

  /*--- Higher order derivatives are not implemented, so we exit if they are requested. ---*/

  if (val_order_der > 2){
    SU2_MPI::Error("Higher order derivatives for BSplines are not implemented.", CURRENT_FUNCTION);
  }
  return 0.0;
}

CBezierBlending::CBezierBlending(short val_order, short n_controlpoints){
  SetOrder(val_order, n_controlpoints);
}

CBezierBlending::~CBezierBlending(){}

void CBezierBlending::SetOrder(short val_order, short n_controlpoints){
  Order  = val_order;
  Degree = Order - 1;
  binomial.resize(Order+1, 0.0);
}

su2double CBezierBlending::GetBasis(short val_i, su2double val_t){
  return GetBernstein(Degree, val_i, val_t);
}

su2double CBezierBlending::GetBernstein(short val_n, short val_i, su2double val_t){

  su2double value = 0.0;

  if (val_i > val_n) { value = 0.0; return value; }

  if (val_i == 0) {
    if (val_t == 0) value = 1.0;
    else if (val_t == 1) value = 0.0;
    else value = Binomial(val_n, val_i) * pow(val_t, val_i) * pow(1.0 - val_t, val_n - val_i);
  }
  else if (val_i == val_n) {
    if (val_t == 0) value = 0.0;
    else if (val_t == 1) value = 1.0;
    else value = pow(val_t, val_n);
  }
  else {
    if ((val_t == 0) || (val_t == 1)) value = 0.0;
    value = Binomial(val_n, val_i) * pow(val_t, val_i) * pow(1.0-val_t, val_n - val_i);
  }

  return value;
}

su2double CBezierBlending::GetDerivative(short val_i, su2double val_t, short val_order_der){
  return GetBernsteinDerivative(Degree, val_i, val_t, val_order_der);
}

su2double CBezierBlending::GetBernsteinDerivative(short val_n, short val_i, su2double val_t, short val_order_der){

  su2double value = 0.0;

  /*--- Verify this subroutine, it provides negative val_n,
   which is a wrong value for GetBernstein ---*/

  if (val_order_der == 0) {
    value = GetBernstein(val_n, val_i, val_t); return value;
  }

  if (val_i == 0) {
    value = val_n*(-GetBernsteinDerivative(val_n-1, val_i, val_t, val_order_der-1)); return value;
  }
  else {
    if (val_n == 0) {
      value = val_t; return value;
    }
    else {
      value = val_n*(GetBernsteinDerivative(val_n-1, val_i-1, val_t, val_order_der-1) - GetBernsteinDerivative(val_n-1, val_i, val_t, val_order_der-1));
      return value;
    }
  }

  return value;
}

su2double CBezierBlending::Binomial(unsigned short n, unsigned short m){

  unsigned short i, j;
  su2double result;

  binomial[0] = 1.0;
  for (i = 1; i <= n; ++i) {
    binomial[i] = 1.0;
    for (j = i-1U; j > 0; --j) {
      binomial[j] += binomial[j-1U];
    }
  }

  result = binomial[m];
  if (fabs(result) < EPS*EPS) { result = 0.0; }

  return result;

}
