/*!
 * \file CRadialBasisFunctionInterpolation.cpp
 * \brief Subroutines for moving mesh volume elements using Radial Basis Function interpolation.
 * \author F. van Steen
 * \version 8.0.1 "Harrier"
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

#include "../../include/grid_movement/CRadialBasisFunctionInterpolation.hpp"
#include "../../include/interface_interpolation/CRadialBasisFunction.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"


CRadialBasisFunctionInterpolation::CRadialBasisFunctionInterpolation(CGeometry* geometry, CConfig* config) : CVolumetricMovement(geometry) {
  /*--- Retrieve type of RBF and if applicable its support radius ---*/
  kindRBF =  config->GetKindRadialBasisFunction();
  radius = config->GetRadialBasisFunctionParameter();

  // controlNodes = &boundaryNodes;//TODO start here, for data reduction the greedy nodes should be selected as control nodes and nr of control nodes should be specified.

  dataReduction  = config->GetRBF_DataReduction();

  if(dataReduction){
    controlNodes = &greedyNodes;
  }else{
    controlNodes = &boundaryNodes;
  }
}

CRadialBasisFunctionInterpolation::~CRadialBasisFunctionInterpolation(void) = default;

void CRadialBasisFunctionInterpolation::SetVolume_Deformation(CGeometry* geometry, CConfig* config, bool UpdateGeo, bool Derivative,
                                                bool ForwardProjectionDerivative){
  
  
  su2double MinVolume, MaxVolume;

  /*--- Retrieving number of deformation steps and screen output from config ---*/

  auto Nonlinear_Iter = config->GetGridDef_Nonlinear_Iter();

  auto Screen_Output = config->GetDeform_Output();
  
  /*--- Disable the screen output if we're running SU2_CFD ---*/

  if (config->GetKind_SU2() == SU2_COMPONENT::SU2_CFD && !Derivative) Screen_Output = false;
  if (config->GetSmoothGradient()) Screen_Output = true;

  
  su2double StartTime = SU2_MPI::Wtime();

  /*--- Assigning the node types ---*/
  SetControlNodes(geometry, config);
  su2double StopTime = SU2_MPI::Wtime();
  auto UsedTimeCompute = StopTime - StartTime;
  
  if(rank==MASTER_NODE){
    cout << "Setting control nodes time: " << UsedTimeCompute << " seconds" << endl;
  }

  StartTime = SU2_MPI::Wtime();
  SetInternalNodes(geometry, config); 
  StopTime = SU2_MPI::Wtime();
  UsedTimeCompute = StopTime - StartTime;
  
  if(rank==MASTER_NODE){
    cout << "Setting internal nodes time: " << UsedTimeCompute << " seconds" << endl;
  }

  /*--- Looping over the number of deformation iterations ---*/
  for (auto iNonlinear_Iter = 0ul; iNonlinear_Iter < Nonlinear_Iter; iNonlinear_Iter++) {
    
    /*--- Compute min volume in the entire mesh. ---*/

    ComputeDeforming_Element_Volume(geometry, MinVolume, MaxVolume, Screen_Output);
    if (rank == MASTER_NODE && Screen_Output)
      cout << "Min. volume: " << MinVolume << ", max. volume: " << MaxVolume << "." << endl;

    /*--- Obtaining the interpolation coefficients of the control nodes ---*/
    GetInterpolationCoefficients(geometry, config, iNonlinear_Iter);
    

    StartTime = SU2_MPI::Wtime();
    /*--- Updating the coordinates of the grid ---*/
    UpdateGridCoord(geometry, config);
    StopTime = SU2_MPI::Wtime();
    UsedTimeCompute = StopTime - StartTime;
    
    if(rank==MASTER_NODE){
      cout << "Updating grid coords time: " << UsedTimeCompute << " seconds" << endl;
    }

    if(UpdateGeo){
      UpdateDualGrid(geometry, config);
    }

    /*--- Check for failed deformation (negative volumes). ---*/

    ComputeDeforming_Element_Volume(geometry, MinVolume, MaxVolume, Screen_Output);

    /*--- Calculate amount of nonconvex elements ---*/

    ComputenNonconvexElements(geometry, Screen_Output);

    if (rank == MASTER_NODE && Screen_Output) {
      cout << "Non-linear iter.: " << iNonlinear_Iter + 1 << "/" << Nonlinear_Iter << ". ";
      if (nDim == 2)
        cout << "Min. area: " << MinVolume <<  "." << endl;
      else
        cout << "Min. volume: " << MinVolume <<  "." << endl;
    }
  }
}

void CRadialBasisFunctionInterpolation::GetInterpolationCoefficients(CGeometry* geometry, CConfig* config, unsigned long iNonlinear_Iter){
  

  if(dataReduction){
    
    GreedyIteration(geometry, config);
  }else{
    //TODO find more elegant way to assign this variable
    Global_nControlNodes = controlNodes->size();
    #ifdef HAVE_MPI
      MPI_Operations(geometry);
    #endif
    

    /*--- Deformation vector only has to be set once ---*/
    su2double StartTime = SU2_MPI::Wtime();
    if(iNonlinear_Iter == 0){
      SetDeformationVector(geometry, config);
    }
    
    su2double StopTime = SU2_MPI::Wtime();
    auto UsedTimeCompute = StopTime - StartTime;
    
    if(rank==MASTER_NODE){
      cout << "Setting deformation vector time: " << UsedTimeCompute << " seconds" << endl;
    }


    /*--- Computing the interpolation matrix with RBF evaluations based on Euclidean distance ---*/
    SetInterpolationMatrix(geometry, config);

    /*--- Solving the RBF system to get the interpolation coefficients ---*/
    StartTime = SU2_MPI::Wtime();
    SolveRBF_System();  
    StopTime = SU2_MPI::Wtime();
    UsedTimeCompute = StopTime - StartTime;
    
    if(rank==MASTER_NODE){
      cout << "Obtaining interpolation coefficients time: " << UsedTimeCompute << " seconds" << endl;
    }
  }

}


void CRadialBasisFunctionInterpolation::SetControlNodes(CGeometry* geometry, CConfig* config){
  
  unsigned short iMarker; 
  unsigned long iVertex, iNode; 

  /*--- Storing of the node, marker and vertex information ---*/

  /*--- Looping over the markers ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    /*--- Checking if not internal or send/receive marker ---*/
    if (!config->GetMarker_All_Deform_Mesh_Internal(iMarker) && !config->GetMarker_All_SendRecv(iMarker)) {

      /*--- Looping over the vertices of marker ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        /*--- Node in consideration ---*/
        iNode = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Check whether node is part of the subdomain and not shared with a receiving marker (for parallel computation)*/
        if (geometry->nodes->GetDomain(iNode)) {
          boundaryNodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));        
        }        
      }
    }
  }

  /*--- Sorting of the boundary nodes based on their index ---*/
  sort(boundaryNodes.begin(), boundaryNodes.end(), Compare);

  /*--- Obtaining unique set of boundary nodes ---*/
  boundaryNodes.resize(std::distance(boundaryNodes.begin(), unique(boundaryNodes.begin(), boundaryNodes.end(), Equal)));

  /*--- Updating the number of boundary nodes ---*/
  nBoundaryNodes = boundaryNodes.size();
}

void CRadialBasisFunctionInterpolation::SetInterpolationMatrix(CGeometry* geometry, CConfig* config){
  su2double StartTime = SU2_MPI::Wtime();
  unsigned long iNode, jNode;
  unsigned long interpMatSize;

  /*--- In case of parallel computation, the interpolation coefficients are computed on the master node.
          In order to do so the coordinates of all control nodes are collected on the master node ---*/


  
  if(rank == MASTER_NODE){
    /*--- Initialization of the interpolation matrix ---*/
    interpMat.Initialize(Global_nControlNodes);

    /*--- Construction of the interpolation matrix. 
      Since this matrix is symmetric only upper halve has to be considered ---*/

  
    /*--- Looping over the target nodes ---*/
    for(iNode = 0; iNode < Global_nControlNodes; iNode++ ){

      /*--- Looping over the control nodes ---*/
      for (jNode = iNode; jNode < Global_nControlNodes; jNode++){
        
        /*--- Distance between nodes ---*/
        #ifdef HAVE_MPI
          su2double dist(0);
          
          for (unsigned short i = 0; i < nDim; i++) dist += pow(GlobalCoords[iNode*nDim+i] - GlobalCoords[jNode*nDim+i], 2);
          dist = sqrt(dist);
          
        #else
          auto dist = GeometryToolbox::Distance(nDim, geometry->nodes->GetCoord((*controlNodes)[iNode]->GetIndex()), geometry->nodes->GetCoord((*controlNodes)[jNode]->GetIndex()));
        #endif
        /*--- Evaluation of RBF ---*/
        interpMat(iNode, jNode) = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(kindRBF, radius, dist));
      }
    }

    // /*--- Obtaining lower halve using symmetry ---*/
    const bool kernelIsSPD = (kindRBF == RADIAL_BASIS::WENDLAND_C2) || (kindRBF == RADIAL_BASIS::GAUSSIAN) ||
                            (kindRBF == RADIAL_BASIS::INV_MULTI_QUADRIC);


    su2double StopTime = SU2_MPI::Wtime();
    auto UsedTimeCompute = StopTime - StartTime;
    
    
    cout << "setting interp matrix time: " << UsedTimeCompute << " seconds" << endl;
    // SU2_MPI::Send(interpMat.data(), interpMat.size(), MPI_DOUBLE, rank, MASTER_NODE, SU2_MPI::GetComm());
    StartTime = SU2_MPI::Wtime();

    
    interpMat.Invert(kernelIsSPD);
    StopTime = SU2_MPI::Wtime();
    UsedTimeCompute = StopTime - StartTime;
    
    
    cout << "Inverting matrix time: " << UsedTimeCompute << " seconds" << endl;
    
  }
}

void CRadialBasisFunctionInterpolation::SetDeformationVector(CGeometry* geometry, CConfig* config){

  /* --- Initialization of the deformation vector ---*/
  deformationVector.resize(controlNodes->size()*nDim, 0.0);

  /*--- If requested (no by default) impose the surface deflections in
    increments and solve the grid deformation with
    successive small deformations. ---*/
  su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());

  /*--- Setting nonzero displacements of the moving markers ---*/
  for (auto i = 0; i < controlNodes->size(); i++) {
    
    if (config->GetMarker_All_Moving((*controlNodes)[i]->GetMarker())) {     
      for (auto iDim = 0; iDim < nDim; iDim++) {
        deformationVector[i+iDim*controlNodes->size()] = SU2_TYPE::GetValue(geometry->vertex[(*controlNodes)[i]->GetMarker()][(*controlNodes)[i]->GetVertex()]->GetVarCoord()[iDim] * VarIncrement);
      }
    }
  }


  #ifdef HAVE_MPI
    
    /*--- Gathering all deformation vectors on the master node ---*/
    if(rank==MASTER_NODE){

      /*--- resizing the global deformation vector ---*/
      deformationVector.resize(Global_nControlNodes*nDim);

      /*--- Receiving the local deformation vector from other processes ---*/
      unsigned long start_idx = 0;
      for (auto iProc = 0; iProc < size; iProc++) {
        if (iProc != MASTER_NODE) {
          SU2_MPI::Recv(&deformationVector[0] + start_idx, Local_nControlNodesVec[iProc]*nDim, MPI_DOUBLE, iProc, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE); // TODO can status ignore be used? 
        }
        start_idx += Local_nControlNodesVec[iProc]*nDim;
      }      

    }else{

      /*--- Sending the local deformation vector to the master node ---*/
      SU2_MPI::Send(deformationVector.data(), Local_nControlNodes*nDim, MPI_DOUBLE, MASTER_NODE, 0, SU2_MPI::GetComm());  
      
    }

    /*--- The global deformation vector is now ordered as d_1, d_2, ..., d_n, where n is the number of processes.
            Here the deformation vector is reordered to obtain an order x_1, ..., x_n, y_1, ..., y_n, z_1, ..., z_n,
                 where x_n is the deformation in x of deformation vector n */
    if(rank == MASTER_NODE){
      
      
      for (unsigned short iDim = nDim-1; iDim > 0; iDim--) {

        unsigned long start_idx = 0;

        for (unsigned short processor = 0; processor < SU2_MPI::GetSize(); processor++) {
          
          if ( processor == 0){
            start_idx += Local_nControlNodesVec[processor];
          }
          else{
            start_idx += (Local_nControlNodesVec[processor-1]*(iDim-1) + Local_nControlNodesVec[processor]);
          }

          /*--- inserting part of vector at end of deformationVector ---*/
          deformationVector.insert(deformationVector.end(), deformationVector.begin()+start_idx, deformationVector.begin()+start_idx+Local_nControlNodesVec[processor]);

          /*--- erasing moved part of the vector ---*/
          deformationVector.erase(deformationVector.begin()+start_idx, deformationVector.begin()+start_idx+Local_nControlNodesVec[processor]);
        }
      }
    }

  #endif   
}

void CRadialBasisFunctionInterpolation::SetInternalNodes(CGeometry* geometry, CConfig* config ){ 

  unsigned long node;

  /*--- resizing the internal nodes vector ---*/
  nInternalNodes = geometry->GetnPoint();// -  nBoundaryNodes; // vector has max size of nPoints
  internalNodes.resize(nInternalNodes);
  

  /*--- Looping over all nodes and check if part of domain and not on boundary ---*/
  unsigned long idx_cnt = 0, idx_control = 0;
  for (unsigned long iNode = 0; iNode < geometry->GetnPoint(); iNode++) {    
    if (!geometry->nodes->GetBoundary(iNode)) {
      internalNodes[idx_cnt++] = iNode; 
    }   
  }  

  

  #ifdef HAVE_MPI
    /*--- Looping over the markers ---*/
    for (unsigned short iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) { 

      /*--- If send or receive marker ---*/
      if (config->GetMarker_All_SendRecv(iMarker)) { 

        /*--- Loop over marker vertices ---*/
        for (unsigned long iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) { 
          
          /*--- Local node index ---*/
          node = geometry->vertex[iMarker][iVertex]->GetNode();

          //   /*--- if not among the boundary nodes ---*/
          if (find_if (boundaryNodes.begin(), boundaryNodes.end(), [&](CRadialBasisFunctionNode* i){return i->GetIndex() == node;}) == boundaryNodes.end()) {
            internalNodes[idx_cnt++] = node;
          }             
        }
      }
    }

    /*--- sorting of the local indices ---*/
    sort(internalNodes.begin(), internalNodes.begin() + idx_cnt);

    /*--- Obtaining unique set of internal nodes ---*/
    internalNodes.resize(std::distance(internalNodes.begin(), unique(internalNodes.begin(), internalNodes.begin() + idx_cnt)));

    /*--- Updating the number of internal nodes ---*/
    nInternalNodes = internalNodes.size();

  #else
    nInternalNodes = idx_cnt;
    internalNodes.resize(nInternalNodes);    
  #endif
}


void CRadialBasisFunctionInterpolation::SolveRBF_System(){
  auto rank = SU2_MPI::GetRank();
  
  unsigned long nControlNodes_local = controlNodes->size();
  
  #ifdef HAVE_MPI
    unsigned long size;
    SU2_MPI::Allreduce(&nControlNodes_local, &size, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  #else
    size = controlNodes->size();
  #endif
  /*--- resizing the interpolation coefficient vector ---*/
  coefficients.resize(nDim*size);

  if(rank == MASTER_NODE){    
    /*--- Looping through the dimensions in order to find the interpolation coefficients for each direction ---*/
    unsigned short iDim;
    for(iDim = 0; iDim < nDim; iDim++){
      interpMat.MatVecMult(deformationVector.begin()+iDim*size, coefficients.begin()+iDim*size);
    }    
  }

  #ifdef HAVE_MPI
    SU2_MPI::Bcast(coefficients.data(), coefficients.size(), MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
  #endif    
}

void CRadialBasisFunctionInterpolation::UpdateGridCoord(CGeometry* geometry, CConfig* config){
  if(rank == MASTER_NODE){
    cout << "updating the grid coordinates" << endl;
  }
  unsigned long iNode, cNode;
  unsigned short iDim;

  /*--- Vector for storing the coordinate variation ---*/
  su2double var_coord[nDim];
  
  /*--- Loop over the internal nodes ---*/
  for(iNode = 0; iNode < nInternalNodes; iNode++){

    /*--- Loop for contribution of each control node ---*/
    for(cNode = 0; cNode < Global_nControlNodes; cNode++){
     
      /*--- Determine distance between considered internal and control node ---*/
      su2double dist;

       #ifdef HAVE_MPI
        dist = 0;        
        for ( iDim = 0; iDim < nDim; iDim++) dist += pow(GlobalCoords[cNode * nDim + iDim] - geometry->nodes->GetCoord(internalNodes[iNode])[iDim], 2);
        dist = sqrt(dist);
      #else
        dist = GeometryToolbox::Distance(nDim, geometry->nodes->GetCoord((*controlNodes)[cNode]->GetIndex()), geometry->nodes->GetCoord(internalNodes[iNode]));
      #endif
      
      /*--- Evaluate RBF based on distance ---*/
      auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(kindRBF, radius, dist));
      
      /*--- Add contribution to total coordinate variation ---*/
      for( iDim = 0; iDim < nDim; iDim++){
        var_coord[iDim] += rbf*coefficients[cNode + iDim*Global_nControlNodes];
      }
    }

    /*--- Apply the coordinate variation and resetting the var_coord vector to zero ---*/
    for(iDim = 0; iDim < nDim; iDim++){
      geometry->nodes->AddCoord(internalNodes[iNode], iDim, var_coord[iDim]);
      var_coord[iDim] = 0;
    } 
  }  
  
  /*--- Applying the surface deformation, which are stored in the deformation vector ---*/

  unsigned long nControlNodes = deformationVector.size()/nDim; // size on master_node is different in case of mpi
  for(cNode = 0; cNode < nControlNodes; cNode++){
    if(config->GetMarker_All_Moving((*controlNodes)[cNode]->GetMarker())){
      for(iDim = 0; iDim < nDim; iDim++){
        #ifdef HAVE_MPI //TODO these are now the same statements
          geometry->nodes->AddCoord((*controlNodes)[cNode]->GetIndex(), iDim, deformationVector[cNode + iDim*nControlNodes]); 
        #else
          geometry->nodes->AddCoord((*controlNodes)[cNode]->GetIndex(), iDim, deformationVector[cNode + iDim*nControlNodes]);
        #endif
      }
    }
  }  
}

void CRadialBasisFunctionInterpolation::GreedyIteration(CGeometry* geometry, CConfig* config) {
  if (rank == MASTER_NODE) {
    cout << "Starting greedy iteration..." << endl;
  }

  GetInitMaxErrorNode(geometry, config);

  unsigned short greedyIter = 0;

  // Gathering the init greedy nodes
  unsigned long MaxErrorNodes[size];
  SU2_MPI::Gather(&MaxErrorNode, 1, MPI_UNSIGNED_LONG, MaxErrorNodes, 1, MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());
  

  //gathering the coordinates of the selected greedy nodes. This array should be available on the masternode throughout the greedy iteration.
 
  vector<su2double> selectedCoords(nDim); // for now a single node is selected; can be an array? 
  auto coord = geometry->nodes->GetCoord(boundaryNodes[MaxErrorNode]->GetIndex());

  for(auto iDim = 0; iDim < nDim; iDim++){
    selectedCoords[iDim] = coord[iDim];
  }
  vector<su2double> greedyCoords;
  if(rank==MASTER_NODE){
    greedyCoords.resize(nDim*size);
  }
  SU2_MPI::Gather(selectedCoords.data(), nDim, MPI_DOUBLE, greedyCoords.data(), nDim, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
  //TODO what if a different number of greedy nodes are selected? 

  if (rank == MASTER_NODE){
    for(auto x : greedyCoords){cout << x << endl;}
  }


  SU2_MPI::Barrier(SU2_MPI::GetComm());
  SU2_MPI::Abort(SU2_MPI::GetComm(), 0);
} 

void CRadialBasisFunctionInterpolation::GetInitMaxErrorNode(CGeometry* geometry, CConfig* config){
  unsigned short iNode;

  su2double maxDeformation = 0.0;
  su2double normSquaredDeformation;
  su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());

  for(iNode = 0; iNode < nBoundaryNodes; iNode++){
    normSquaredDeformation = GeometryToolbox::SquaredNorm(nDim, geometry->vertex[boundaryNodes[iNode]->GetMarker()][boundaryNodes[iNode]->GetVertex()]->GetVarCoord());    
    if(normSquaredDeformation > maxDeformation){
      maxDeformation = normSquaredDeformation;
      MaxErrorNode = iNode;
    }
  }
  MaxError = sqrt(maxDeformation) / ((su2double)config->GetGridDef_Nonlinear_Iter());
  // cout << "rank: " << rank << " Max error node: " << MaxErrorNode << endl;
}


void CRadialBasisFunctionInterpolation::MPI_Operations(CGeometry* geometry){
  Local_nControlNodes = controlNodes->size();

  Local_nControlNodesVec.resize(size);

  /*--- gathering local control node coordinate sizes on all processes. ---*/
  SU2_MPI::Allgather(&Local_nControlNodes, 1, MPI_UNSIGNED_LONG, Local_nControlNodesVec.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm()); 


  Global_nControlNodes = 0;
  for( auto& n : Local_nControlNodesVec) Global_nControlNodes += n;

  /*--- array containing the global control node coordinates. ---*/
  GlobalCoords.resize(Global_nControlNodes*nDim);
  
  /*--- array containing the local control node coordinates. ---*/
  vector<su2double> LocalCoords(nDim*Local_nControlNodes);

  
  /*--- storing local control node coordinates ---*/
  for(unsigned long iNode = 0; iNode < controlNodes->size(); iNode++){
    auto coord = geometry->nodes->GetCoord((*controlNodes)[iNode]->GetIndex());  
    for ( unsigned short iDim = 0 ; iDim < nDim; iDim++ ){
      LocalCoords[ iNode * nDim + iDim ] = coord[iDim];
    }
  }

  /*--- array containing size of local control node coordinates. ---*/
  int LocalCoordsSizes[SU2_MPI::GetSize()];

  int localCoordsSize = LocalCoords.size();
  /*--- gathering local control node coordinate sizes on all processes. ---*/
  SU2_MPI::Allgather(&localCoordsSize, 1, MPI_INT, LocalCoordsSizes, 1, MPI_INT, SU2_MPI::GetComm()); 
  

  /*--- array containing the starting indices for the allgatherv operation*/
  int disps[SU2_MPI::GetSize()];    

  for(auto x = 0; x < SU2_MPI::GetSize(); x++){
    if(x == 0){
      disps[x] = 0;
    }else{
      disps[x] = disps[x-1]+LocalCoordsSizes[x-1];
    }
  }
  
  /*--- making global control node coordinates available on all processes ---*/
  SU2_MPI::Allgatherv(LocalCoords.data(), localCoordsSize, MPI_DOUBLE, GlobalCoords.data(), LocalCoordsSizes, disps, MPI_DOUBLE, SU2_MPI::GetComm()); //TODO local coords can be deleted after this operation
};