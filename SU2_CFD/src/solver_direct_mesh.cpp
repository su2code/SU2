/*!
 * \file solver_direct_mesh.cpp
 * \brief Main subroutines to solve moving meshes using a pseudo-linear elastic approach.
 * \author R. Sanchez
 * \version 6.1.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

CMeshSolver::CMeshSolver(CGeometry *geometry, CConfig *config) : CSolver() {

    /*--- Initialize some booleans that determine the kind of problem at hand. ---*/

    bool time_domain = config->GetTime_Domain();
    bool multizone = config->GetMultizone_Problem();

    /*--- Initialize the number of spatial dimensions, length of the state
     vector (same as spatial dimensions for grid deformation), and grid nodes. ---*/

    unsigned short iDim, jDim;
    unsigned long iPoint, iElem;

    nDim         = geometry->GetnDim();
    nVar         = geometry->GetnDim();
    nPoint       = geometry->GetnPoint();
    nPointDomain = geometry->GetnPointDomain();
    nElem        = geometry->GetnElem();

    nIterMesh   = 0;
    valResidual = 0.0;

    MinVolume_Ref = 0.0;
    MinVolume_Curr = 0.0;

    MaxVolume_Ref = 0.0;
    MaxVolume_Curr = 0.0;

    /*--- Initialize the node structure ---*/
    Coordinate = new su2double[nDim];
    node       = new CMeshVariable*[nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++){

      /*--- We store directly the reference coordinates ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Coordinate[iDim] = geometry->node[iPoint]->GetCoord(iDim);

      node[iPoint] = new CMeshVariable(Coordinate, nDim, config);
    }

    /*--- Initialize the element structure ---*/
    element = new CMeshElement[nElem];
    for (iElem = 0; iElem < nElem; iElem++)
        element[iElem] = CMeshElement();

    Residual = new su2double[nDim];   for (iDim = 0; iDim < nDim; iDim++) Residual[iDim] = 0.0;
    Solution = new su2double[nDim];   for (iDim = 0; iDim < nDim; iDim++) Solution[iDim] = 0.0;

    /*--- Initialize matrix, solution, and r.h.s. structures for the linear solver. ---*/

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    StiffMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);

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
      element_container[EL_TRIA] = new CTRIA1(nDim, config);
      element_container[EL_QUAD] = new CQUAD4(nDim, config);
    }
    else if (nDim == 3){
      element_container[EL_TETRA] = new CTETRA1(nDim, config);
      element_container[EL_HEXA] = new CHEXA8(nDim, config);
      element_container[EL_PYRAM] = new CPYRAM5(nDim, config);
      element_container[EL_PRISM] = new CPRISM6(nDim, config);
    }

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

    /*--- Compute the element volumes using the reference coordinates ---*/
    SetMinMaxVolume(geometry, config, false);

    /*--- Compute the wall distance using the reference coordinates ---*/
    SetWallDistance(geometry, config);

}

CMeshSolver::~CMeshSolver(void) {

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

void CMeshSolver::SetMinMaxVolume(CGeometry *geometry, CConfig *config, bool updated) {

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
        element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }
    }

    /*--- Compute the volume of the element (or the area in 2D cases ) ---*/

    if (nDim == 2)  ElemVolume = element_container[EL_KIND]->ComputeArea();
    else            ElemVolume = element_container[EL_KIND]->ComputeVolume();

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

void CMeshSolver::SetWallDistance(CGeometry *geometry, CConfig *config) {

  unsigned long nVertex_SolidWall, ii, jj, iVertex, iPoint, pointID;
  unsigned long iElem, PointCorners[8];
  unsigned short iNodes, nNodes;
  unsigned short iMarker, iDim;
  su2double dist, MaxDistance_Local, MinDistance_Local;
  su2double nodeDist, ElemDist;
  int rankID;

  /*--- Initialize min and max distance ---*/

  MaxDistance = -1E22; MinDistance = 1E22;

  /*--- Compute the total number of nodes on no-slip boundaries ---*/

  nVertex_SolidWall = 0;
  for(iMarker=0; iMarker<config->GetnMarker_All(); ++iMarker) {
    if( (config->GetMarker_All_KindBC(iMarker) == EULER_WALL ||
         config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)  ||
       (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ) {
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
       (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ) {
      for (iVertex=0; iVertex<geometry->GetnVertex(iMarker); ++iVertex) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        PointIDs[jj++] = iPoint;
        for (iDim=0; iDim<nDim; ++iDim){
          Coord_bound[ii++] = node[iPoint]->GetMesh_Coord(iDim);
        }
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

    for(iPoint=0; iPoint< nPoint; ++iPoint) {

      WallADT.DetermineNearestNode(node[iPoint]->GetMesh_Coord(), dist,
                                   pointID, rankID);
      node[iPoint]->SetWallDistance(dist);

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

  /*--- Normalize distance from 0 to 1 ---*/
  for (iPoint=0; iPoint < nPoint; ++iPoint) {
    nodeDist = node[iPoint]->GetWallDistance()/MaxDistance;
    node[iPoint]->SetWallDistance(nodeDist);
  }

  /*--- Compute the element distances ---*/

  for (iElem = 0; iElem < nElem; iElem++) {

    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      nNodes = 3;
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       nNodes = 5;
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         nNodes = 6;
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    nNodes = 8;

    for (iNodes = 0; iNodes < nNodes; iNodes++) {
      PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
    }

    /*--- Average the distance of the nodes in the element ---*/

    ElemDist = 0.0;
    for (iNodes = 0; iNodes < nNodes; iNodes++){
      ElemDist += node[PointCorners[iNodes]]->GetWallDistance();
    }
    ElemDist = ElemDist/(su2double)nNodes;

    element[iElem].SetWallDistance(ElemDist);

  }

}

void CMeshSolver::SetStiffnessMatrix(CGeometry *geometry, CConfig *config){

  unsigned long iElem;
  unsigned short iNode, iDim, jDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord;
  int EL_KIND = 0;

  su2double *Kab = NULL;
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
        val_Coord = node[indexNode[iNode]]->GetMesh_Coord(iDim);
        element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
      }

    }

    /*--- Compute the stiffness of the element ---*/
    Set_Element_Stiffness(iElem, geometry, config);

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

void CMeshSolver::Set_Element_Stiffness(unsigned long iElem, CGeometry *geometry, CConfig *config) {

  su2double ElemVolume, ElemDistance;

  ElemVolume = element[iElem].GetRef_Volume();
  ElemDistance = element[iElem].GetWallDistance();

  switch (config->GetDeform_Stiffness_Type()) {
    case INVERSE_VOLUME:
      E = 1.0 / ElemVolume;           // Stiffness inverse of the volume of the element
      Nu = config->GetDeform_Coeff(); // Nu is normally a very large number, for rigid-body rotations, see Dwight (2009)
      break;
    case SOLID_WALL_DISTANCE:
      E = 1.0 / ElemDistance;         // Stiffness inverse of the distance of the element to the closest wall
      Nu = config->GetDeform_Coeff(); // Nu is normally a very large number, for rigid-body rotations, see Dwight (2009)
      break;
    case CONSTANT_STIFFNESS:
      E  = config->GetDeform_ElasticityMod();
      Nu = config->GetDeform_PoissonRatio();
      break;
    case VOLUME_DISTANCE:
      E  = 1.0 / (ElemDistance * ElemVolume);
      Nu = config->GetDeform_PoissonRatio();
      break;
  }

  /*--- Lamé parameters ---*/

  Mu     = E / (2.0*(1.0 + Nu));
  Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));

}

void CMeshSolver::DeformMesh(CGeometry **geometry, CConfig *config){

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

    /*--- Compute the minimum and maximum area/volume for the mesh. ---*/
    if ((rank == MASTER_NODE) && (!discrete_adjoint)) {
      if (nDim == 2) cout << scientific << "Min. area in the undeformed mesh: "<< MinVolume_Ref <<", max. area: " << MaxVolume_Ref <<"." << endl;
      else           cout << scientific << "Min. volume in the undeformed mesh: "<< MinVolume_Ref <<", max. volume: " << MaxVolume_Ref <<"." << endl;
    }

    /*--- Compute the stiffness matrix. ---*/
    SetStiffnessMatrix(geometry[MESH_0], config);

    /*--- Impose boundary conditions (all of them are ESSENTIAL BC's - displacements). ---*/
    SetBoundaryDisplacements(geometry[MESH_0], config);

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

void CMeshSolver::UpdateGridCoord(CGeometry *geometry, CConfig *config){

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

  /* --- LinSysSol contains the non-transformed displacements in the periodic halo cells.
   * Hence we still need a communication of the transformed coordinates, otherwise periodicity
   * is not maintained. ---*/

  geometry->Set_MPI_Coord(config);

}

void CMeshSolver::UpdateDualGrid(CGeometry *geometry, CConfig *config){

  /*--- After moving all nodes, update the dual mesh. Recompute the edges and
   dual mesh control volumes in the domain and on the boundaries. ---*/

  geometry->SetCoord_CG();
  geometry->SetControlVolume(config, UPDATE);
  geometry->SetBoundControlVolume(config, UPDATE);
  geometry->SetMaxLength(config);

}

void CMeshSolver::ComputeGridVelocity(CGeometry *geometry, CConfig *config){

  /*--- Local variables ---*/

  su2double *Disp_nP1 = NULL, *Disp_n = NULL, *Disp_nM1 = NULL;
  su2double TimeStep, GridVel = 0.0;
  unsigned long iPoint;
  unsigned short iDim;

  /*--- Compute the velocity of each node in the volume mesh ---*/

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

  /*--- The velocity is computed for nPointDomain, now we communicate it ---*/
  geometry->Set_MPI_GridVel(config);

}

void CMeshSolver::UpdateMultiGrid(CGeometry **geometry, CConfig *config){

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

void CMeshSolver::SetBoundaryDisplacements(CGeometry *geometry, CConfig *config){

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
        SetClamped_Boundary(geometry, config, iMarker);
      }
    }
  }

  /*--- All others are pending. ---*/

}

void CMeshSolver::SetClamped_Boundary(CGeometry *geometry, CConfig *config, unsigned short val_marker){

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
        valJacobian_ij_00 = StiffMatrix.GetBlock(iNode, jPoint,0,0);

        if (valJacobian_ij_00 != 0.0 ){
          /*--- Set the rest of the row to 0 ---*/
          if (iNode != jPoint) {
            StiffMatrix.SetBlock(iNode,jPoint,matrixZeros);
          }
          /*--- And the diagonal to 1.0 ---*/
          else{
            StiffMatrix.SetBlock(iNode,jPoint,matrixId);
          }
        }
      }

      /*--- Delete the full column for node iNode ---*/
      for (iPoint = 0; iPoint < nPoint; iPoint++){

        /*--- Check whether the block is non-zero ---*/
        valJacobian_ij_00 = StiffMatrix.GetBlock(iPoint, iNode,0,0);

        if (valJacobian_ij_00 != 0.0 ){
          /*--- Set the rest of the row to 0 ---*/
          if (iNode != iPoint) {
            StiffMatrix.SetBlock(iPoint,iNode,matrixZeros);
          }
        }
      }
    }
  }

}

void CMeshSolver::SetMoving_Boundary(CGeometry *geometry, CConfig *config, unsigned short val_marker){

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

    /*--- Add it to the current displacement ---*/
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
        valJacobian_ij_00 = StiffMatrix.GetBlock(iNode, jPoint,0,0);

        if (valJacobian_ij_00 != 0.0 ){
          /*--- Set the rest of the row to 0 ---*/
          if (iNode != jPoint) {
            StiffMatrix.SetBlock(iNode,jPoint,matrixZeros);
          }
          /*--- And the diagonal to 1.0 ---*/
          else{
            StiffMatrix.SetBlock(iNode,jPoint,matrixId);
          }
        }
      }

      /*--- Delete the columns for a particular node ---*/

      for (iPoint = 0; iPoint < nPoint; iPoint++){

        /*--- Check if the term K(iPoint, iNode) is 0 ---*/
        valJacobian_ij_00 = StiffMatrix.GetBlock(iPoint,iNode,0,0);

        /*--- If the node iNode has a crossed dependency with the point iPoint ---*/
        if (valJacobian_ij_00 != 0.0 ){

          /*--- Retrieve the Jacobian term ---*/
          for (iDim = 0; iDim < nDim; iDim++){
            for (jDim = 0; jDim < nDim; jDim++){
              auxJacobian_ij[iDim][jDim] = StiffMatrix.GetBlock(iPoint,iNode,iDim,jDim);
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
            StiffMatrix.SetBlock(iPoint,iNode,matrixZeros);
          }
        }
      }
    }
  }

}

void CMeshSolver::Solve_System_Mesh(CGeometry *geometry, CConfig *config){

  /*--- Retrieve number or iterations, tol, output, etc. from config ---*/

  su2double SolverTol = config->GetDeform_Linear_Solver_Error(), System_Residual = 1.0;

  unsigned long MaxIter = config->GetDeform_Linear_Solver_Iter();
  unsigned long IterLinSol = 0;

  bool Screen_Output= config->GetDeform_Output();

  /*--- Initialize the structures to solve the system ---*/

  CMatrixVectorProduct* mat_vec = NULL;
  CSysSolve *system  = new CSysSolve();

  bool TapeActive = NO;

  if (config->GetDiscrete_Adjoint()){
#ifdef CODI_REVERSE_TYPE

    /*--- Check whether the tape is active, i.e. if it is recording and store the status ---*/

    TapeActive = AD::globalTape.isActive();


    /*--- Stop the recording for the linear solver ---*/

    AD::StopRecording();
#endif
  }

  /*--- Communicate any prescribed boundary displacements via MPI,
   so that all nodes have the same solution and r.h.s. entries
   across all partitions. ---*/

  StiffMatrix.SendReceive_Solution(LinSysSol, geometry, config);
  StiffMatrix.SendReceive_Solution(LinSysRes, geometry, config);

  /*--- Solve the linear system using a Krylov subspace method ---*/

  if (config->GetKind_Deform_Linear_Solver() == BCGSTAB ||
      config->GetKind_Deform_Linear_Solver() == FGMRES ||
      config->GetKind_Deform_Linear_Solver() == RESTARTED_FGMRES ||
      config->GetKind_Deform_Linear_Solver() == CONJUGATE_GRADIENT) {

    /*--- Independently of whether we are using or not derivatives,
     *--- as the matrix is now symmetric, the matrix-vector product
     *--- can be done in the same way in the forward and the reverse mode
     */

    /*--- Definition of the preconditioner matrix vector multiplication, and linear solver ---*/
    mat_vec = new CSysMatrixVectorProduct(StiffMatrix, geometry, config);
    CPreconditioner* precond = NULL;

    switch (config->GetKind_Deform_Linear_Solver_Prec()) {
    case JACOBI:
      StiffMatrix.BuildJacobiPreconditioner();
      precond = new CJacobiPreconditioner(StiffMatrix, geometry, config);
      break;
    case ILU:
      StiffMatrix.BuildILUPreconditioner();
      precond = new CILUPreconditioner(StiffMatrix, geometry, config);
      break;
    case LU_SGS:
      precond = new CLU_SGSPreconditioner(StiffMatrix, geometry, config);
      break;
    case LINELET:
      StiffMatrix.BuildJacobiPreconditioner();
      precond = new CLineletPreconditioner(StiffMatrix, geometry, config);
      break;
    default:
      StiffMatrix.BuildJacobiPreconditioner();
      precond = new CJacobiPreconditioner(StiffMatrix, geometry, config);
      break;
    }

    switch (config->GetKind_Deform_Linear_Solver()) {
    case BCGSTAB:
      IterLinSol = system->BCGSTAB_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &System_Residual, Screen_Output);
      break;
    case FGMRES:
      IterLinSol = system->FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &System_Residual, Screen_Output);
      break;
    case CONJUGATE_GRADIENT:
      IterLinSol = system->CG_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &System_Residual, Screen_Output);
      break;
    case RESTARTED_FGMRES:
      IterLinSol = 0;
      while (IterLinSol < config->GetLinear_Solver_Iter()) {
        if (IterLinSol + config->GetLinear_Solver_Restart_Frequency() > config->GetLinear_Solver_Iter())
          MaxIter = config->GetLinear_Solver_Iter() - IterLinSol;
        IterLinSol += system->FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &System_Residual, Screen_Output);
        if (LinSysRes.norm() < SolverTol) break;
        SolverTol = SolverTol*(1.0/LinSysRes.norm());
      }
      break;
    }

    /*--- Dealocate memory of the Krylov subspace method ---*/

    delete mat_vec;
    delete precond;

  }

  if(TapeActive){
    /*--- Start recording if it was stopped for the linear solver ---*/

    AD::StartRecording();

    /*--- Prepare the externally differentiated linear solver ---*/

    system->SetExternalSolve_Mesh(StiffMatrix, LinSysRes, LinSysSol, geometry, config);

  }

  delete system;

  /*--- Set number of iterations in the mesh update. ---*/

  nIterMesh = IterLinSol;

  /*--- Store the value of the residual. ---*/

  valResidual = System_Residual;


}

void CMeshSolver::Compute_Element_Contribution(CElement *element, CConfig *config){

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

  element->clearElement();      /*--- Restarts the element: avoids adding over previous results in other elements --*/
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

        element->Add_Kab(KAux_ab,iNode, jNode);
        /*--- Symmetric terms --*/
        if (iNode != jNode){
          element->Add_Kab_T(KAux_ab, jNode, iNode);
        }

      }

    }

  }

}

void CMeshSolver::Compute_Constitutive_Matrix(void){

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

//void CMeshSolver::Transfer_Boundary_Displacements(CGeometry *geometry, CConfig *config, unsigned short val_marker){

//  unsigned short iDim;
//  unsigned long iNode, iVertex;

//  su2double *VarCoord;
//  su2double new_coord;

//  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

//    /*--- Get node index ---*/

//    iNode = geometry->vertex[val_marker][iVertex]->GetNode();

//    /*--- Get the displacement on the vertex ---*/

//    VarCoord = geometry->vertex[val_marker][iVertex]->GetVarCoord();

//    if (geometry->node[iNode]->GetDomain()) {

//      /*--- Update the grid coordinates using the solution of the structural problem
//       *--- recorded in VarCoord. */

//      for (iDim = 0; iDim < nDim; iDim++) {
//        new_coord = geometry->node[iNode]->GetCoord(iDim)+VarCoord[iDim];
//        if (fabs(new_coord) < EPS*EPS) new_coord = 0.0;
//        geometry->node[iNode]->SetCoord(iDim, new_coord);
//      }
//    }

//  }

//}

//void CMeshSolver::Boundary_Dependencies(CGeometry **geometry, CConfig *config){


//  unsigned short iMarker;

//  /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
//   deforming meshes (MARKER_FSI_INTERFACE). ---*/

//  unsigned short Kind_SU2 = config->GetKind_SU2();

//  /*--- Set the dependencies on the FSI interfaces. ---*/

//  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
//    if ((config->GetMarker_All_ZoneInterface(iMarker) != 0) && (Kind_SU2 == SU2_CFD)) {
//      Transfer_Boundary_Displacements(geometry[MESH_0], config, iMarker);
//    }
//  }

//  UpdateDualGrid(geometry[MESH_0], config);

//}


void CMeshSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iDim;
  unsigned long index;
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;

  su2double curr_coord, displ;

  string UnstExt, text_line;
  ifstream restart_file;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  string restart_filename = config->GetSolution_FlowFileName();

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    restart_filename = config->GetMultizone_FileName(restart_filename, iZone);

  /*--- Modify file name for an unsteady restart ---*/

  if (dual_time || time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1];
      for (iDim = 0; iDim < nDim; iDim++){
        /*--- Update the coordinates of the mesh ---*/
        curr_coord = Restart_Data[index+iDim];
        geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, curr_coord);
        /*--- Store the displacements computed as the current coordinates
         minus the coordinates of the reference mesh file ---*/
        displ = curr_coord - node[iPoint_Local]->GetMesh_Coord(iDim);
        node[iPoint_Local]->SetDisplacement(iDim, displ);
      }
      iPoint_Global_Local++;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
      SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Communicate the loaded displacements. ---*/
  solver[MESH_0][MESH_SOL]->Set_MPI_Solution(geometry[MESH_0], config);

  /*--- Communicate the new coordinates and grid velocities at the halos ---*/
  geometry[MESH_0]->Set_MPI_Coord(config);

  /*--- Recompute the edges and dual mesh control volumes in the
   domain and on the boundaries. ---*/
  UpdateDualGrid(geometry[MESH_0], config);

  /*--- For time-domain problems, we need to compute the grid velocities ---*/
  if (time_domain){
    /*--- Update the old geometry (coordinates n and n-1) ---*/
    Restart_OldGeometry(geometry[MESH_0], config);
    /*--- Once Displacement_n and Displacement_n1 are filled,
     we can compute the Grid Velocity ---*/
    ComputeGridVelocity(geometry[MESH_0], config);
  }

  /*--- Update the multigrid structure after setting up the finest grid,
   including computing the grid velocities on the coarser levels
   when the problem is unsteady. ---*/
  UpdateMultiGrid(geometry, config);

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

void CMeshSolver::Restart_OldGeometry(CGeometry *geometry, CConfig *config) {

  /*--- This function is intended for dual time simulations ---*/

  unsigned long index;

  int Unst_RestartIter;
  ifstream restart_file_n;
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  unsigned short iDim;
  string filename = config->GetSolution_FlowFileName();
  string filename_n;

  /*--- Auxiliary variables for storing the coordinates ---*/
  su2double curr_coord, displ;

  /*--- Variables for reading the restart files ---*/
  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  /*-------------------------------------------------------------------------------------------*/
  /*----------------------- First, load the restart file for time n ---------------------------*/
  /*-------------------------------------------------------------------------------------------*/

  /*--- Modify file name for an unsteady restart ---*/
  Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
  filename_n = config->GetUnsteady_FileName(filename, Unst_RestartIter);

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry, config, filename_n);
  } else {
    Read_SU2_Restart_ASCII(geometry, config, filename_n);
  }

  /*--- Load data from the restart into correct containers. ---*/
  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1];
      for (iDim = 0; iDim < nDim; iDim++){
        curr_coord = Restart_Data[index+iDim];
        displ = curr_coord - node[iPoint_Local]->GetMesh_Coord(iDim);
        node[iPoint_Local]->SetDisplacement_n(iDim, displ);
      }
      iPoint_Global_Local++;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
      SU2_MPI::Error(string("The solution file ") + filename_n + string(" doesn't match with the mesh file!\n") +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

  /*-------------------------------------------------------------------------------------------*/
  /*------------ Now, load the restart file for time n-1, if the simulation is 2nd Order ------*/
  /*-------------------------------------------------------------------------------------------*/

  if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) {

    ifstream restart_file_n1;
    string filename_n1;

    /*--- Restart the variables for reading the restart files ---*/
    counter = 0; iPoint_Local = 0; iPoint_Global = 0; iPoint_Global_Local = 0;
    rbuf_NotMatching = 0; sbuf_NotMatching = 0;

    /*--- Modify file name for an unsteady restart ---*/
    Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
    filename_n1 = config->GetUnsteady_FileName(filename, Unst_RestartIter);

    /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

    if (config->GetRead_Binary_Restart()) {
      Read_SU2_Restart_Binary(geometry, config, filename_n1);
    } else {
      Read_SU2_Restart_ASCII(geometry, config, filename_n1);
    }

    /*--- Load data from the restart into correct containers. ---*/
    counter = 0;
    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {

        /*--- We need to store this point's data, so jump to the correct
         offset in the buffer of data from the restart file and load it. ---*/

        index = counter*Restart_Vars[1];
        for (iDim = 0; iDim < nDim; iDim++){
          curr_coord = Restart_Data[index+iDim];
          displ = curr_coord - node[iPoint_Local]->GetMesh_Coord(iDim);
          node[iPoint_Local]->SetDisplacement_n1(iDim, displ);
        }
        iPoint_Global_Local++;

        /*--- Increment the overall counter for how many points have been loaded. ---*/
        counter++;
      }

    }

    /*--- Detect a wrong solution file ---*/

    if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

  #ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
  #else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
  #endif
    if (rbuf_NotMatching != 0) {
        SU2_MPI::Error(string("The solution file ") + filename_n1 + string(" doesn't match with the mesh file!\n") +
                       string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
    }

    /*--- Delete the class memory that is used to load the restart. ---*/

    if (Restart_Vars != NULL) delete [] Restart_Vars;
    if (Restart_Data != NULL) delete [] Restart_Data;
    Restart_Vars = NULL; Restart_Data = NULL;

  }

  /*--- It's necessary to communicate this information ---*/
  geometry->Set_MPI_OldCoord(config);

}

