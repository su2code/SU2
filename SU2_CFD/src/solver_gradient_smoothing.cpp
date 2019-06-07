/*!
 * \file solver_gradient_smoothing.cpp
 * \brief Main subroutines for the gradient smoothing problem.
 * \author T. Dick
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../include/solver_structure.hpp"
#include "../../Common/include/adt_structure.hpp"
#include <algorithm>


CGradientSmoothingSolver::CGradientSmoothingSolver(CGeometry *geometry, CConfig *config) : CSolver(), System(false, true) {

  nDim         = geometry->GetnDim();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  su2double **matrixId;

  LinSysSol.Initialize(nPoint, nPointDomain, nDim, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nDim, 0.0);
  Jacobian.Initialize(nPoint, nPointDomain, nDim, nDim, false, geometry, config);

}


CGradientSmoothingSolver::~CGradientSmoothingSolver(void) {

  delete [] Residual;
  delete [] matrixId;
}


void CGradientSmoothingSolver::ApplyGradientSmoothing(CGeometry *geometry, CSolver *solver, CConfig *config) {

  /*--- Initialize vector and sparse matrix ---*/
  LinSysSol.SetValZero();
  LinSysRes.SetValZero();
  Jacobian.SetValZero();

  Compute_Residual(geometry, solver, config);

  Compute_StiffMatrix(geometry, NULL, config);

  Solve_Linear_System(geometry, config);

  Set_Sensitivities(geometry, solver, config);

}


void CGradientSmoothingSolver::Compute_Residual(CGeometry *geometry, CSolver *solver, CConfig *config){

  unsigned long iPoint;
  unsigned short iDim;
  Residual = new su2double [nDim];

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Residual[iDim] = solver->node[iPoint]->GetSensitivity(iDim);
    }
    LinSysRes.SetBlock(iPoint, Residual);
  }
}


void CGradientSmoothingSolver::Compute_StiffMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config){

  unsigned long iPoint;
  unsigned short iDim, jDim;

  matrixId    = new su2double *[nDim];
  for(iDim = 0; iDim < nDim; iDim++){
    matrixId[iDim]    = new su2double[nDim];
  }
  for(iDim = 0; iDim < nDim; iDim++){
    for (jDim = 0; jDim < nDim; jDim++){
      matrixId[iDim][jDim]    = 0.0;
    }
    matrixId[iDim][iDim] = 2.0;
  }

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    Jacobian.SetBlock(iPoint,iPoint,matrixId);
  }

}


void CGradientSmoothingSolver::Solve_Linear_System(CGeometry *geometry, CConfig *config){

  unsigned long IterLinSol = 0;

  IterLinSol = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  SetIterLinSolver(IterLinSol);

  // valResidual = System.GetResidual();

}


void CGradientSmoothingSolver::Set_Sensitivities(CGeometry *geometry, CSolver *solver, CConfig *config){

  unsigned long iPoint, total_index;
  unsigned short iDim;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      total_index = iPoint*nDim + iDim;
      solver->node[iPoint]->SetSensitivity(iDim,LinSysSol[total_index]);
    }
  }

}

#ifdef SOMEBIZAREFLAG
// new functions from mesh solver

void CGradientSmoothingSolver::SetMinMaxVolume(CGeometry *geometry, CConfig *config, bool updated) {

  unsigned long iElem, ElemCounter = 0;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  su2double MaxVolume, MinVolume;
  int EL_KIND = 0;

  bool discrete_adjoint = config->GetDiscrete_Adjoint();

  bool RightVol = true;

  su2double ElemVolume;

  if ((rank == MASTER_NODE) && (!discrete_adjoint))
    cout << "Computing volumes of the grid elements." << endl;

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

      /*--- Compute the volume with the reference or with the current coordinates ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        if (updated) val_Coord = node[indexNode[iNode]]->GetMesh_Coord(iDim)
                               + node[indexNode[iNode]]->GetDisplacement(iDim);
        else val_Coord = node[indexNode[iNode]]->GetMesh_Coord(iDim);
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }
    }

    /*--- Compute the volume of the element (or the area in 2D cases ) ---*/

    if (nDim == 2)  ElemVolume = element_container[FEA_TERM][EL_KIND]->ComputeArea();
    else            ElemVolume = element_container[FEA_TERM][EL_KIND]->ComputeVolume();

    RightVol = true;
    if (ElemVolume < 0.0) RightVol = false;

    MaxVolume = max(MaxVolume, ElemVolume);
    MinVolume = min(MinVolume, ElemVolume);
    if (updated) element[iElem].SetCurr_Volume(ElemVolume);
    else element[iElem].SetRef_Volume(ElemVolume);

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
    if (updated){
      ElemVolume = element[iElem].GetCurr_Volume()/MaxVolume;
      element[iElem].SetCurr_Volume(ElemVolume);
    }
    else{
      ElemVolume = element[iElem].GetRef_Volume()/MaxVolume;
      element[iElem].SetRef_Volume(ElemVolume);
    }
  }

  /*--- Store the maximum and minimum volume ---*/
  if (updated){
    MaxVolume_Curr = MaxVolume;
    MinVolume_Curr = MinVolume;
  }
  else{
    MaxVolume_Ref = MaxVolume;
    MinVolume_Ref = MinVolume;
  }

  if ((ElemCounter != 0) && (rank == MASTER_NODE))
    cout <<"There are " << ElemCounter << " elements with negative volume.\n" << endl;

}

void CGradientSmoothingSolver::Compute_StiffMatrix(CGeometry *geometry, CNumerics **numerics, CConfig *config){

  unsigned long iElem, iVar, jVar;
  unsigned short iNode, iDim, jDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol;
  int EL_KIND = 0;

  su2double *Kab = NULL, *Ta = NULL;
  unsigned short NelNodes, jNode;

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
        val_Coord = Get_ValCoord(geometry, indexNode[iNode], iDim); //geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = Get_ValSol(indexNode[iNode], iDim) + val_Coord; //node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
      }

    }

    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);

    numerics[FEA_TERM]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);

    /*--- Retrieve number of nodes ---*/

    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();

    for (iNode = 0; iNode < NelNodes; iNode++) {

      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];

      LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);

      for (jNode = 0; jNode < NelNodes; jNode++) {

        Kab = element_container[FEA_TERM][EL_KIND]->Get_Kab(iNode, jNode);
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            Jacobian_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
          }
        }
        Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);
      }

    }

  }

}

void CGradientSmoothingSolver::SetMesh_Stiffness(CGeometry **geometry, CNumerics **numerics, CConfig *config){

  unsigned long iElem;

  if(!stiffness_set){
    for (iElem = 0; iElem < nElement; iElem++) {

      switch (config->GetDeform_Stiffness_Type()) {
      /*--- Stiffness inverse of the volume of the element ---*/
      case INVERSE_VOLUME: E = 1.0 / element[iElem].GetRef_Volume();  break;
      /*--- Stiffness inverse of the distance of the element to the closest wall ---*/
      case SOLID_WALL_DISTANCE: E = 1.0 / element[iElem].GetWallDistance(); break;
      }

      /*--- Set the element elastic properties in the numerics container ---*/
      numerics[FEA_TERM]->SetMeshElasticProperties(iElem, E);

    }

    stiffness_set = true;
  }

}

void CGradientSmoothingSolver::DeformMesh(CGeometry **geometry, CNumerics **numerics, CConfig *config){

  unsigned long iNonlinear_Iter, Nonlinear_Iter = 0;

  bool discrete_adjoint = config->GetDiscrete_Adjoint();

  if (multizone) SetDisplacement_Old();

  /*--- Retrieve number or internal iterations from config ---*/

  Nonlinear_Iter = config->GetGridDef_Nonlinear_Iter();

  /*--- Loop over the total number of grid deformation iterations. The surface
   deformation can be divided into increments to help with stability. ---*/

  for (iNonlinear_Iter = 0; iNonlinear_Iter < Nonlinear_Iter; iNonlinear_Iter++) {

    /*--- Initialize vector and sparse matrix ---*/
    LinSysSol.SetValZero();
    LinSysRes.SetValZero();
    Jacobian.SetValZero();

    /*--- Compute the minimum and maximum area/volume for the mesh. ---*/
    if ((rank == MASTER_NODE) && (!discrete_adjoint)) {
      if (nDim == 2) cout << scientific << "Min. area in the undeformed mesh: "<< MinVolume_Ref <<", max. area: " << MaxVolume_Ref <<"." << endl;
      else           cout << scientific << "Min. volume in the undeformed mesh: "<< MinVolume_Ref <<", max. volume: " << MaxVolume_Ref <<"." << endl;
    }

    /*--- Compute the stiffness matrix. ---*/
    Compute_StiffMatrix(geometry[MESH_0], numerics, config);

    /*--- Impose boundary conditions (all of them are ESSENTIAL BC's - displacements). ---*/
    SetBoundaryDisplacements(geometry[MESH_0], numerics[FEA_TERM], config);

    /*--- Solve the linear system. ---*/
    Solve_System_Mesh(geometry[MESH_0], config);

    /*--- Update the grid coordinates and cell volumes using the solution
     of the linear system (usol contains the x, y, z displacements). ---*/
    UpdateGridCoord(geometry[MESH_0], config);

    /*--- Update the dual grid. ---*/
    UpdateDualGrid(geometry[MESH_0], config);

    /*--- Check for failed deformation (negative volumes). ---*/
    /*--- In order to do this, we recompute the minimum and maximum area/volume for the mesh using the current coordinates. ---*/
    SetMinMaxVolume(geometry[MESH_0], config, true);

    if ((rank == MASTER_NODE) && (!discrete_adjoint)) {
      cout << scientific << "Non-linear iter.: " << iNonlinear_Iter+1 << "/" << Nonlinear_Iter  << ". Linear iter.: " << nIterMesh << ". ";
      if (nDim == 2) cout << "Min. area in the deformed mesh: " << MinVolume_Curr << ". Error: " << valResidual << "." << endl;
      else cout << "Min. volume in the deformed mesh: " << MinVolume_Curr << ". Error: " << valResidual << "." << endl;
    }

  }

  /*--- The Grid Velocity is only computed if the problem is time domain ---*/
  if (time_domain) ComputeGridVelocity(geometry[MESH_0], config);

  /*--- Update the dual grid. ---*/
  UpdateMultiGrid(geometry, config);

}

void CGradientSmoothingSolver::UpdateGridCoord(CGeometry *geometry, CConfig *config){

  unsigned short iDim;
  unsigned long iPoint, total_index;
  su2double val_disp, val_coord;

  /*--- Update the grid coordinates using the solution of the linear system
     after grid deformation (LinSysSol contains the x, y, z displacements). ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iDim = 0; iDim < nDim; iDim++) {
      total_index = iPoint*nDim + iDim;
      /*--- Retrieve the displacement from the solution of the linear system ---*/
      val_disp = LinSysSol[total_index];
      /*--- Store the displacement of the mesh node ---*/
      node[iPoint]->SetDisplacement(iDim, val_disp);
      /*--- Compute the current coordinate as Mesh_Coord + Displacement ---*/
      val_coord = node[iPoint]->GetMesh_Coord(iDim) + val_disp;
      /*--- Update the geometry container ---*/
      geometry->node[iPoint]->SetCoord(iDim, val_coord);
    }
  }

  /*--- LinSysSol contains the non-transformed displacements in the periodic halo cells.
   Hence we still need a communication of the transformed coordinates, otherwise periodicity
   is not maintained. ---*/
  geometry->Set_MPI_Coord(config);

  /*--- In the same way, communicate the displacements in the solver to make sure the halo
   nodes receive the correct value of the displacement. ---*/
  Set_MPI_Displacement(geometry,config);

}

void CGradientSmoothingSolver::UpdateDualGrid(CGeometry *geometry, CConfig *config){

  /*--- After moving all nodes, update the dual mesh. Recompute the edges and
   dual mesh control volumes in the domain and on the boundaries. ---*/

  geometry->SetCoord_CG();
  geometry->SetControlVolume(config, UPDATE);
  geometry->SetBoundControlVolume(config, UPDATE);
  geometry->SetMaxLength(config);

}

void CGradientSmoothingSolver::ComputeGridVelocity(CGeometry *geometry, CConfig *config){

  /*--- Local variables ---*/

  su2double *Disp_nP1 = NULL, *Disp_n = NULL, *Disp_nM1 = NULL;
  su2double TimeStep, GridVel = 0.0;
  unsigned long iPoint;
  unsigned short iDim;

  /*--- Compute the velocity of each node in the domain of the current rank
   (halo nodes are not computed as the grid velocity is later communicated) ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Coordinates of the current point at n+1, n, & n-1 time levels ---*/

    Disp_nM1 = node[iPoint]->GetDisplacement_n1();
    Disp_n   = node[iPoint]->GetDisplacement_n();
    Disp_nP1 = node[iPoint]->GetDisplacement();

    /*--- Unsteady time step ---*/

    TimeStep = config->GetDelta_UnstTimeND();

    /*--- Compute mesh velocity with 1st or 2nd-order approximation ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        GridVel = ( Disp_nP1[iDim] - Disp_n[iDim] ) / TimeStep;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
        GridVel = ( 3.0*Disp_nP1[iDim] - 4.0*Disp_n[iDim]
                    +  1.0*Disp_nM1[iDim] ) / (2.0*TimeStep);

      /*--- Store grid velocity for this point ---*/

      geometry->node[iPoint]->SetGridVel(iDim, GridVel);

    }
  }

  /*--- The velocity was computed for nPointDomain, now we communicate it ---*/
  geometry->Set_MPI_GridVel(config);

}

void CGradientSmoothingSolver::UpdateMultiGrid(CGeometry **geometry, CConfig *config){

  unsigned short iMGfine, iMGlevel, nMGlevel = config->GetnMGLevels();

  /*--- Update the multigrid structure after moving the finest grid,
   including computing the grid velocities on the coarser levels
   when the problem is solved in unsteady conditions. ---*/

  for (iMGlevel = 1; iMGlevel <= nMGlevel; iMGlevel++) {
    iMGfine = iMGlevel-1;
    geometry[iMGlevel]->SetControlVolume(config, geometry[iMGfine], UPDATE);
    geometry[iMGlevel]->SetBoundControlVolume(config, geometry[iMGfine],UPDATE);
    geometry[iMGlevel]->SetCoord(geometry[iMGfine]);
    if (time_domain)
      geometry[iMGlevel]->SetRestricted_GridVelocity(geometry[iMGfine], config);
  }

}

void CGradientSmoothingSolver::SetBoundaryDisplacements(CGeometry *geometry, CNumerics *numerics, CConfig *config){

  unsigned short iMarker;

  /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
   deforming meshes (MARKER_FSI_INTERFACE). ---*/

  unsigned short Kind_SU2 = config->GetKind_SU2();

  /*--- First of all, move the FSI interfaces. ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_ZoneInterface(iMarker) != 0) && (Kind_SU2 == SU2_CFD)) {
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
         (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY))) {

      /*--- We must note that the FSI surfaces are not clamped ---*/
      if (config->GetMarker_All_ZoneInterface(iMarker) == 0){
        BC_Clamped(geometry, numerics, config, iMarker);
      }
    }
  }

  /*--- All others are pending. ---*/

}

void CGradientSmoothingSolver::BC_Clamped(CGeometry *geometry, CNumerics *numerics, CConfig *config, unsigned short val_marker){

  unsigned long iNode, iVertex;
  unsigned long iPoint, jPoint;

  su2double valJacobian_ij_00 = 0.0;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/

    iNode = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iNode]->GetDomain()) {

      if (nDim == 2) {
        Solution[0] = 0.0;  Solution[1] = 0.0;
        Residual[0] = 0.0;  Residual[1] = 0.0;
      }
      else {
        Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
        Residual[0] = 0.0;  Residual[1] = 0.0;  Residual[2] = 0.0;
      }

      /*--- Initialize the reaction vector ---*/

      LinSysRes.SetBlock(iNode, Residual);
      LinSysSol.SetBlock(iNode, Solution);

      /*--- STRONG ENFORCEMENT OF THE CLAMPED BOUNDARY CONDITION ---*/

      /*--- Delete the full row for node iNode ---*/
      for (jPoint = 0; jPoint < nPoint; jPoint++){

        /*--- Check whether the block is non-zero ---*/
        valJacobian_ij_00 = Jacobian.GetBlock(iNode, jPoint,0,0);

        if (valJacobian_ij_00 != 0.0 ){
          /*--- Set the rest of the row to 0 ---*/
          if (iNode != jPoint) {
            Jacobian.SetBlock(iNode,jPoint,matrixZeros);
          }
          /*--- And the diagonal to 1.0 ---*/
          else{
            Jacobian.SetBlock(iNode,jPoint,matrixId);
          }
        }
      }

      /*--- Delete the full column for node iNode ---*/
      for (iPoint = 0; iPoint < nPoint; iPoint++){

        /*--- Check whether the block is non-zero ---*/
        valJacobian_ij_00 = Jacobian.GetBlock(iPoint, iNode,0,0);

        if (valJacobian_ij_00 != 0.0 ){
          /*--- Set the rest of the row to 0 ---*/
          if (iNode != iPoint) {
            Jacobian.SetBlock(iPoint,iNode,matrixZeros);
          }
        }
      }
    }
  }

}

void CGradientSmoothingSolver::SetMoving_Boundary(CGeometry *geometry, CConfig *config, unsigned short val_marker){

  unsigned short iDim, jDim;

  su2double *VarDisp = NULL;

  unsigned long iNode, iVertex;
  unsigned long iPoint, jPoint;

  su2double VarIncrement = 1.0/((su2double)config->GetGridDef_Nonlinear_Iter());

  su2double valJacobian_ij_00 = 0.0;
  su2double auxJacobian_ij[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  su2double VarCoord[3] = {0.0, 0.0, 0.0};

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/

    iNode = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Get the displacement on the vertex ---*/
    VarDisp = geometry->vertex[val_marker][iVertex]->GetVarCoord();

    /*--- Add it to the current displacement (this will be replaced in the transfer routines
     so that the displacements at the interface are directly transferred instead of incrementally) ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      VarCoord[iDim] = node[iNode]->GetDisplacement(iDim) + VarDisp[iDim];

    if (geometry->node[iNode]->GetDomain()) {

      if (nDim == 2) {
        Solution[0] = VarCoord[0] * VarIncrement;  Solution[1] = VarCoord[1] * VarIncrement;
        Residual[0] = VarCoord[0] * VarIncrement;  Residual[1] = VarCoord[1] * VarIncrement;
      }
      else {
        Solution[0] = VarCoord[0] * VarIncrement;  Solution[1] = VarCoord[1] * VarIncrement;  Solution[2] = VarCoord[2] * VarIncrement;
        Residual[0] = VarCoord[0] * VarIncrement;  Residual[1] = VarCoord[1] * VarIncrement;  Residual[2] = VarCoord[2] * VarIncrement;
      }

      /*--- Initialize the reaction vector ---*/

      LinSysRes.SetBlock(iNode, Residual);
      LinSysSol.SetBlock(iNode, Solution);

      /*--- STRONG ENFORCEMENT OF THE DISPLACEMENT BOUNDARY CONDITION ---*/

      /*--- Delete the full row for node iNode ---*/
      for (jPoint = 0; jPoint < nPoint; jPoint++){

        /*--- Check whether the block is non-zero ---*/
        valJacobian_ij_00 = Jacobian.GetBlock(iNode, jPoint,0,0);

        if (valJacobian_ij_00 != 0.0 ){
          /*--- Set the rest of the row to 0 ---*/
          if (iNode != jPoint) {
            Jacobian.SetBlock(iNode,jPoint,matrixZeros);
          }
          /*--- And the diagonal to 1.0 ---*/
          else{
            Jacobian.SetBlock(iNode,jPoint,matrixId);
          }
        }
      }

      /*--- Delete the columns for a particular node ---*/

      for (iPoint = 0; iPoint < nPoint; iPoint++){

        /*--- Check if the term K(iPoint, iNode) is 0 ---*/
        valJacobian_ij_00 = Jacobian.GetBlock(iPoint,iNode,0,0);

        /*--- If the node iNode has a crossed dependency with the point iPoint ---*/
        if (valJacobian_ij_00 != 0.0 ){

          /*--- Retrieve the Jacobian term ---*/
          for (iDim = 0; iDim < nDim; iDim++){
            for (jDim = 0; jDim < nDim; jDim++){
              auxJacobian_ij[iDim][jDim] = Jacobian.GetBlock(iPoint,iNode,iDim,jDim);
            }
          }

          /*--- Multiply by the imposed displacement ---*/
          for (iDim = 0; iDim < nDim; iDim++){
            Residual[iDim] = 0.0;
            for (jDim = 0; jDim < nDim; jDim++){
              Residual[iDim] += auxJacobian_ij[iDim][jDim] * VarCoord[jDim];
            }
          }

          /*--- For the whole column, except the diagonal term ---*/
          if (iNode != iPoint) {
            /*--- The term is substracted from the residual (right hand side) ---*/
            LinSysRes.SubtractBlock(iPoint, Residual);
            /*--- The Jacobian term is now set to 0 ---*/
            Jacobian.SetBlock(iPoint,iNode,matrixZeros);
          }
        }
      }
    }
  }

}

void CGradientSmoothingSolver::Solve_System_Mesh(CGeometry *geometry, CConfig *config){


  unsigned long IterLinSol = 0, iPoint, total_index;
  unsigned short iVar;

  /*--- Initialize residual and solution at the ghost points ---*/

  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {

    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }

  }

  IterLinSol = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  /*--- The the number of iterations of the linear solver ---*/

  SetIterLinSolver(IterLinSol);

  /*--- Store the value of the residual. ---*/

  valResidual = System.GetResidual();

}

void CGradientSmoothingSolver::Transfer_Boundary_Displacements(CGeometry *geometry, CConfig *config, unsigned short val_marker){

  unsigned short iDim;
  unsigned long iPoint, iVertex;

  su2double *VarCoord = NULL, *VarDisp = NULL;
  su2double new_coord;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Get node index ---*/
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Get the displacement on the vertex ---*/
    VarDisp = geometry->vertex[val_marker][iVertex]->GetVarCoord();

    /*--- Add it to the current displacement (this will be replaced in the transfer routines
     so that the displacements at the interface are directly transferred instead of incrementally) ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      VarCoord[iDim] = node[iPoint]->GetDisplacement(iDim) + VarDisp[iDim];

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Update the grid coordinates using the solution of the structural problem recorded in VarCoord. */
      for (iDim = 0; iDim < nDim; iDim++) {
        new_coord = node[iPoint]->GetMesh_Coord(iDim)+VarCoord[iDim];
        if (fabs(new_coord) < EPS*EPS) new_coord = 0.0;
        geometry->node[iPoint]->SetCoord(iDim, new_coord);
      }
    }

  }

}

void CGradientSmoothingSolver::Boundary_Dependencies(CGeometry **geometry, CConfig *config){

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

void CGradientSmoothingSolver::ComputeResidual_Multizone(CGeometry *geometry, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++){
      SetRes_BGS(iVar,0.0);
      SetRes_Max_BGS(iVar,0.0,0);
  }

  /*--- Set the residuals ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
          residual = node[iPoint]->GetDisplacement(iVar) - node[iPoint]->GetDisplacement_Old(iVar);
          AddRes_BGS(iVar,residual*residual);
          AddRes_Max_BGS(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
  }

  SetResidual_BGS(geometry, config);

}

void CGradientSmoothingSolver::Set_MPI_Displacement(CGeometry *geometry, CConfig *config) {

  unsigned short iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_Displ = NULL, *Buffer_Send_Displ = NULL, *Displ = NULL, *newDispl = NULL;
  su2double *translation;
  newDispl = new su2double[nDim];

#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nDim;        nBufferR_Vector = nVertexR*nDim;

      /*--- Allocate Receive and send buffers  ---*/

      Buffer_Receive_Displ = new su2double [nBufferR_Vector];
      Buffer_Send_Displ = new su2double[nBufferS_Vector];

      /*--- Copy the coordinates that should be sent ---*/

      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Displ = node[iPoint]->GetDisplacement();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Displ[iDim*nVertexS+iVertex] = Displ[iDim];
      }

#ifdef HAVE_MPI
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Displ, nBufferS_Vector, MPI_DOUBLE, send_to,0,
                   Buffer_Receive_Displ, nBufferR_Vector, MPI_DOUBLE, receive_from,0, MPI_COMM_WORLD, &status);
#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Receive_Displ[iDim*nVertexR+iVertex] = Buffer_Send_Displ[iDim*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/

      delete [] Buffer_Send_Displ;

      /*--- Do the coordinate transformation ---*/

      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/

        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

        /*--- Retrieve the supplied periodic information. ---*/

        angles = config->GetPeriodicRotation(iPeriodic_Index);
        translation = config->GetPeriodicTranslate(iPeriodic_Index);

        /*--- Store angles separately for clarity. ---*/

        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/

        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy coordinates before performing transformation. ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          newDispl[iDim] = Buffer_Receive_Displ[iDim*nVertexR+iVertex];

        /*--- Rotate the coordinates. ---*/

        if (nDim == 2) {
          newDispl[0] = (rotMatrix[0][0]*Buffer_Receive_Displ[0*nVertexR+iVertex] +
                         rotMatrix[0][1]*Buffer_Receive_Displ[1*nVertexR+iVertex]) - translation[0];
          newDispl[1] = (rotMatrix[1][0]*Buffer_Receive_Displ[0*nVertexR+iVertex] +
                         rotMatrix[1][1]*Buffer_Receive_Displ[1*nVertexR+iVertex]) - translation[1];
        }
        else {
          newDispl[0] = (rotMatrix[0][0]*Buffer_Receive_Displ[0*nVertexR+iVertex] +
                         rotMatrix[0][1]*Buffer_Receive_Displ[1*nVertexR+iVertex] +
                         rotMatrix[0][2]*Buffer_Receive_Displ[2*nVertexR+iVertex]);
          newDispl[1] = (rotMatrix[1][0]*Buffer_Receive_Displ[0*nVertexR+iVertex] +
                         rotMatrix[1][1]*Buffer_Receive_Displ[1*nVertexR+iVertex] +
                         rotMatrix[1][2]*Buffer_Receive_Displ[2*nVertexR+iVertex]);
          newDispl[2] = (rotMatrix[2][0]*Buffer_Receive_Displ[0*nVertexR+iVertex] +
                         rotMatrix[2][1]*Buffer_Receive_Displ[1*nVertexR+iVertex] +
                         rotMatrix[2][2]*Buffer_Receive_Displ[2*nVertexR+iVertex]);
        }

        /*--- Copy transformed coordinates back into buffer. ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          node[iPoint]->SetDisplacement(iDim, newDispl[iDim]);

      }

      /*--- Deallocate receive buffer. ---*/

      delete [] Buffer_Receive_Displ;

    }

  }

  delete [] newDispl;

}


// old functions

void CGradientSmoothingSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {


  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.

  unsigned short nSolVar;

  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;

#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];

      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
        if (dynamic) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Vel(iVar);
            Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Accel(iVar);
          }
        }
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
        if (dynamic) {
          for (iVar = nVar; iVar < 3*nVar; iVar++)
            Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
        }
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nSolVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, SolRest[iVar]);

        if (dynamic) {

          for (iVar = 0; iVar < nVar; iVar++) {
            node[iPoint]->SetSolution_Vel(iVar, SolRest[iVar+nVar]);
            node[iPoint]->SetSolution_Accel(iVar, SolRest[iVar+2*nVar]);
          }

        }

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;

    }

  }

}

void CGradientSmoothingSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {


  unsigned long iPoint;

  /*--- Set vector entries to zero ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    LinSysRes.SetBlock_Zero(iPoint);
    LinSysSol.SetBlock_Zero(iPoint);
  }

  /*--- Set matrix entries to zero ---*/

  StiffnessMatrix.SetValZero();

}

void CGradientSmoothingSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {

  unsigned long iPoint, nPoint;
  bool incremental_load = config->GetIncrementalLoad();              // If an incremental load is applied

  nPoint = geometry[MESH_0]->GetnPoint();

  /*--- We store the current solution as "Solution Old", for the case that we need to retrieve it ---*/

  if (incremental_load) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Set_OldSolution();
  }


}

void CGradientSmoothingSolver::Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {

  unsigned long iElem, iVar, jVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol;
  int EL_KIND = 0;

  su2double *Kab = NULL, *Ta  = NULL;
  unsigned short NelNodes, jNode;

  /*--- Loops over all the elements ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    {nNodes = 8; EL_KIND = EL_HEXA;}

    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {

      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
      }
    }

    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);

    /*--- Compute the components of the jacobian and the stress term ---*/
    numerics[FEA_TERM]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);


    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();

    for (iNode = 0; iNode < NelNodes; iNode++) {

      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];

      LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);

      for (jNode = 0; jNode < NelNodes; jNode++) {

        Kab = element_container[FEA_TERM][EL_KIND]->Get_Kab(iNode, jNode);

        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            Jacobian_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
          }
        }

        StiffnessMatrix.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);
      }

    }

  }


}

void CGradientSmoothingSolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

}

void CGradientSmoothingSolver::Compute_IntegrationConstants(CConfig *config) {

  su2double Delta_t= config->GetDelta_DynTime();

  su2double gamma = config->GetNewmark_gamma(), beta = config->GetNewmark_beta();

  switch (config->GetKind_TimeIntScheme_FEA()) {
    case (CD_EXPLICIT):
      cout << "NOT IMPLEMENTED YET" << endl;
      break;
    case (NEWMARK_IMPLICIT):

      /*--- Integration constants for Newmark scheme ---*/

      a_dt[0]= 1 / (beta*pow(Delta_t,2.0));
      a_dt[1]= gamma / (beta*Delta_t);
      a_dt[2]= 1 / (beta*Delta_t);
      a_dt[3]= 1 /(2*beta) - 1;
      a_dt[4]= gamma/beta - 1;
      a_dt[5]= (Delta_t/2) * (gamma/beta - 2);
      a_dt[6]= Delta_t * (1-gamma);
      a_dt[7]= gamma * Delta_t;
      a_dt[8]= 0.0;

      break;

    case (GENERALIZED_ALPHA):

      /*--- Integration constants for Generalized Alpha ---*/
      /*--- Needs to be updated if accounting for structural damping ---*/

      //      su2double beta = config->Get_Int_Coeffs(0);
      //      //  su2double gamma =  config->Get_Int_Coeffs(1);
      //      su2double alpha_f = config->Get_Int_Coeffs(2), alpha_m =  config->Get_Int_Coeffs(3);
      //
      //      a_dt[0]= (1 / (beta*pow(Delta_t,2.0))) * ((1 - alpha_m) / (1 - alpha_f)) ;
      //      a_dt[1]= 0.0 ;
      //      a_dt[2]= (1 - alpha_m) / (beta*Delta_t);
      //      a_dt[3]= ((1 - 2*beta)*(1-alpha_m) / (2*beta)) - alpha_m;
      //      a_dt[4]= 0.0;
      //      a_dt[5]= 0.0;
      //      a_dt[6]= Delta_t * (1-delta);
      //      a_dt[7]= delta * Delta_t;
      //      a_dt[8]= (1 - alpha_m) / (beta*pow(Delta_t,2.0));

      break;
  }


}

void CGradientSmoothingSolver::BC_Impose(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {}

void CGradientSmoothingSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
                                           unsigned short iMesh) {

  unsigned short iVar;
  unsigned long iPoint, total_index;

  bool first_iter = (config->GetIntIter() == 0);
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);    // Nonlinear analysis.
  bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);

  su2double solNorm = 0.0, solNorm_recv = 0.0;

  if (disc_adj_fem) {

    if (nonlinear_analysis) {

      /*--- For nonlinear discrete adjoint, we have 3 convergence criteria ---*/

      /*--- UTOL = norm(Delta_U(k)): ABSOLUTE, norm of the incremental displacements ------------*/
      /*--- RTOL = norm(Residual(k): ABSOLUTE, norm of the residual (T-F) -----------------------*/
      /*--- ETOL = Delta_U(k) * Residual(k): ABSOLUTE, norm of the product disp * res -----------*/

      Conv_Check[0] = LinSysSol.norm();               // Norm of the delta-solution vector
      Conv_Check[1] = LinSysRes.norm();               // Norm of the residual
      Conv_Check[2] = dotProd(LinSysSol, LinSysRes);  // Position for the energy tolerance

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);
    }
    else {
      /*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/
      /*---  Compute the residual Ax-f ---*/

      Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);

      /*--- Set maximum residual to zero ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
      }

      /*--- Compute the residual ---*/

      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
          AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
        }
      }

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);

      /*--- Compute the root mean square residual ---*/

      SetResidual_RMS(geometry, config);

    }

  }
  else {
    if (nonlinear_analysis){

      /*--- If the problem is nonlinear, we have 3 convergence criteria ---*/

      /*--- UTOL = norm(Delta_U(k)) / norm(U(k)) --------------------------*/
      /*--- RTOL = norm(Residual(k)) / norm(Residual(0)) ------------------*/
      /*--- ETOL = Delta_U(k) * Residual(k) / Delta_U(0) * Residual(0) ----*/

      if (first_iter){
        Conv_Ref[0] = 1.0;                                        // Position for the norm of the solution
        Conv_Ref[1] = max(LinSysRes.norm(), EPS);                 // Position for the norm of the residual
        Conv_Ref[2] = max(dotProd(LinSysSol, LinSysRes), EPS);    // Position for the energy tolerance

        /*--- Make sure the computation runs at least 2 iterations ---*/
        Conv_Check[0] = 1.0;
        Conv_Check[1] = 1.0;
        Conv_Check[2] = 1.0;

        /*--- If absolute, we check the norms ---*/
        switch (config->GetResidual_Criteria_FEM()) {
          case RESFEM_ABSOLUTE:
            Conv_Check[0] = LinSysSol.norm();         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm();         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes);  // Position for the energy tolerance
            break;
        }
      }
      else {
        /*--- Compute the norm of the solution vector Uk ---*/
        for (iPoint = 0; iPoint < nPointDomain; iPoint++){
          for (iVar = 0; iVar < nVar; iVar++){
            solNorm += node[iPoint]->GetSolution(iVar) * node[iPoint]->GetSolution(iVar);
          }
        }

        // We need to communicate the norm of the solution and compute the RMS throughout the different processors

#ifdef HAVE_MPI
        /*--- We sum the squares of the norms across the different processors ---*/
        SU2_MPI::Allreduce(&solNorm, &solNorm_recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        solNorm_recv         = solNorm;
#endif

        Conv_Ref[0] = max(sqrt(solNorm_recv), EPS);           // Norm of the solution vector

        switch (config->GetResidual_Criteria_FEM()) {
          case RESFEM_RELATIVE:
            Conv_Check[0] = LinSysSol.norm() / Conv_Ref[0];         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm() / Conv_Ref[1];         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes) / Conv_Ref[2];  // Position for the energy tolerance
            break;
          case RESFEM_ABSOLUTE:
            Conv_Check[0] = LinSysSol.norm();         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm();         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes);  // Position for the energy tolerance
            break;
          default:
            Conv_Check[0] = LinSysSol.norm() / Conv_Ref[0];         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm() / Conv_Ref[1];         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes) / Conv_Ref[2];  // Position for the energy tolerance
            break;
        }

      }

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);

    } else {

      /*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/

      /*---  Compute the residual Ax-f ---*/

      Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);

      /*--- Set maximum residual to zero ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
      }

      /*--- Compute the residual ---*/

      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
          AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
        }
      }

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);

      /*--- Compute the root mean square residual ---*/

      SetResidual_RMS(geometry, config);
    }

  }

}

void CGradientSmoothingSolver::Solve_System(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned long IterLinSol = 0, iPoint, total_index;
  unsigned short iVar;

  /*--- Initialize residual and solution at the ghost points ---*/

  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {

    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }

  }

  CSysSolve LinSystem;
  IterLinSol = LinSystem.Solve(StiffnessMatrix, LinSysRes, LinSysSol, geometry, config);

  /*--- The the number of iterations of the linear solver ---*/

  SetIterLinSolver(IterLinSol);

}

su2double CGradientSmoothingSolver::GetRes(unsigned short val_var) {}

void CGradientSmoothingSolver::Set_ElementProperties(CGeometry *geometry, CConfig *config) {

  unsigned long iElem;
  unsigned long index;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();

  bool topology_mode = config->GetTopology_Optimization();

  string filename;
  ifstream properties_file;

  element_properties = new CElementProperty*[nElement];

  /*--- Restart the solution from file information ---*/

  filename = config->GetFEA_FileName();

  /*--- If multizone, append zone name ---*/
  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  if (rank == MASTER_NODE) cout << "Filename: " << filename << "." << endl;

  properties_file.open(filename.data(), ios::in);

  /*--- In case there is no file, all elements get the same property (0) ---*/

  if (properties_file.fail()) {
    if (rank == MASTER_NODE){
      cout << "There is no element-based properties file." << endl;
      cout << "The structural domain has uniform properties." << endl;

      if (topology_mode)
        SU2_MPI::Error("Topology mode requires an element-based properties file.",CURRENT_FUNCTION);
    }

    for (iElem = 0; iElem < nElement; iElem++){
      element_properties[iElem] = new CElementProperty(0, 0, 0, 0);
    }

    element_based = false;

  }
  else{

    element_based = true;

    /*--- In case this is a parallel simulation, we need to perform the
       Global2Local index transformation first. ---*/

    long *Global2Local = new long[geometry->GetGlobal_nElemDomain()];

    /*--- First, set all indices to a negative value by default ---*/

    for (iElem = 0; iElem < geometry->GetGlobal_nElemDomain(); iElem++)
      Global2Local[iElem] = -1;

    /*--- Now fill array with the transform values only for the points in the rank (including halos) ---*/

    for (iElem = 0; iElem < nElement; iElem++)
      Global2Local[geometry->elem[iElem]->GetGlobalIndex()] = iElem;

    /*--- Read all lines in the restart file ---*/

    long iElem_Local;
    unsigned long iElem_Global_Local = 0, iElem_Global = 0; string text_line;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

    /*--- The first line is the header ---*/

    getline (properties_file, text_line);

    for (iElem_Global = 0; iElem_Global < geometry->GetGlobal_nElemDomain(); iElem_Global++ ) {

      getline (properties_file, text_line);

      istringstream point_line(text_line);

      /*--- Retrieve local index. If this element from the restart file lives
         only on a different processor, the value of iPoint_Local will be -1.
         Otherwise, the local index for this node on the current processor
         will be returned and used to instantiate the vars. ---*/

      iElem_Local = Global2Local[iElem_Global];

      if (iElem_Local >= 0) {

        if (config->GetDE_Effects())
          point_line >> index >> elProperties[0] >> elProperties[1] >> elProperties[2] >> elProperties[3];
        else
          point_line >> index >> elProperties[0] >> elProperties[1] >> elProperties[2] >> elProperties[3];

        element_properties[iElem_Local] = new CElementProperty(elProperties[0],
                                                         elProperties[1],
                                                         elProperties[2],
                                                         elProperties[3]);

        /*--- For backwards compatibility we only read a fifth column in topology mode ---*/
        if (topology_mode) {
          su2double elDensity;
          point_line >> elDensity;
          element_properties[iElem_Local]->SetDesignDensity(elDensity);
        }

        iElem_Global_Local++;
      }

    }

    /*--- Detect a wrong solution file ---*/

    if (iElem_Global_Local < nElement) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
#else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (rbuf_NotMatching != 0) {
      SU2_MPI::Error(string("The properties file ") + filename + string(" doesn't match with the mesh file!\n")  +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
    }

    /*--- Close the restart file ---*/

    properties_file.close();

    /*--- Free memory needed for the transformation ---*/

    delete [] Global2Local;


  }


}

void CGradientSmoothingSolver::ExtractAdjoint_Variables(CGeometry *geometry, CSolver **solver_container, CConfig *config)
{

    /* possible MPI reduce?? */

    /* get the adjoint solution from solver_container */
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++){
        node[iPoint]->SetSensitivity(iDim,solver_container[ADJFLOW_SOL]->node[iPoint]->GetSensitivity(iDim);
    }
}



#endif
