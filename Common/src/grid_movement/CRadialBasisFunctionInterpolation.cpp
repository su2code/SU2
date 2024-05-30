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

  controlNodes = &boundaryNodes;

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


  /*--- Assigning the node types ---*/
  SetControlNodes(geometry, config);
  SetInternalNodes(geometry, config); 


  /*--- Looping over the number of deformation iterations ---*/
  for (auto iNonlinear_Iter = 0ul; iNonlinear_Iter < Nonlinear_Iter; iNonlinear_Iter++) {
    
    /*--- Compute min volume in the entire mesh. ---*/

    ComputeDeforming_Element_Volume(geometry, MinVolume, MaxVolume, Screen_Output);
    if (rank == MASTER_NODE && Screen_Output)
      cout << "Min. volume: " << MinVolume << ", max. volume: " << MaxVolume << "." << endl;

    /*--- Obtaining the interpolation coefficients of the control nodes ---*/
    GetInterpolationCoefficients(geometry, config, iNonlinear_Iter);
    
    /*--- Updating the coordinates of the grid ---*/
    UpdateGridCoord(geometry, config);


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

  /*--- Deformation vector only has to be set once ---*/
  if(iNonlinear_Iter == 0){
    SetDeformationVector(geometry, config);
  }
 
  /*--- Computing the interpolation matrix with RBF evaluations based on Euclidean distance ---*/
  SetInterpolationMatrix(geometry, config);
  
  /*--- Solving the RBF system to get the interpolation coefficients ---*/
  SolveRBF_System();  

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
  
  unsigned long iNode, jNode;
  unsigned long interpMatSize;

  /*--- In case of parallel computation, the interpolation coefficients are computed on the master node.
          In order to do so the coordinates of all control nodes are collected on the master node ---*/

  #ifdef HAVE_MPI

    /*--- Obtaining global number of control nodes on master node ---*/
    unsigned long size; //TODO use different variable here, as size is used as nr of processes outside of this scope
    SU2_MPI::Reduce(&nBoundaryNodes, &size, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());

    /*--- For other processes size is equal to local number of control nodes ---*/
    if( rank != MASTER_NODE ){
      size = nBoundaryNodes;
    }
    
    /*--- array containing coordinates of control nodes. Global control nodes for master node, 
            local control nodes for other processes. ---*/
    su2double coords[size*nDim];

    /*--- Storing coordinates in array ---*/
    for (iNode = 0; iNode < controlNodes->size(); iNode++) {

      auto coord = geometry->nodes->GetCoord((*controlNodes)[iNode]->GetIndex()); 

      for (unsigned short iDim = 0; iDim < nDim; iDim++) { 
        coords[iNode*nDim+iDim] = coord[iDim];
      }
    }

    /*--- Array containing the sizes of the local coordinates. ---*/
    unsigned long localSizes[SU2_MPI::GetSize()];
    
    /*--- local size of the coordinates ---*/
    auto local_size = nBoundaryNodes*nDim;

    if( rank == MASTER_NODE ){

      /*--- gathering the local sizes ---*/
      SU2_MPI::Gather(&local_size, 1, MPI_UNSIGNED_LONG, localSizes, 1, MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());

      /*--- receiving local coordinates from other processes ---*/
      unsigned long start_idx = 0;
      for(auto iProc = 0; iProc < SU2_MPI::GetSize(); iProc++){

        if(iProc != MASTER_NODE){
          SU2_MPI::Recv(&coords[0] + start_idx, localSizes[iProc], MPI_DOUBLE, iProc, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE); // TODO can status ignore be used? 
        }
        start_idx += localSizes[iProc];
      }
      
    }else{

      /*--- gathering local coordinate size ---*/
       SU2_MPI::Gather(&local_size, 1, MPI_UNSIGNED_LONG, NULL, 1, MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());

      /*--- sending local coordinates to the master node ---*/
      SU2_MPI::Send(coords, local_size, MPI_DOUBLE, MASTER_NODE, 0, SU2_MPI::GetComm());
    }
    /*--- setting size of the interpolation matrix ---*/
    interpMatSize = size;

  #else
    /*--- setting size of the interpolation matrix ---*/
    interpMatSize = nBoundaryNodes;
  #endif
  
  if(rank == MASTER_NODE){
    /*--- Initialization of the interpolation matrix ---*/
    interpMat.Initialize(interpMatSize);

    /*--- Construction of the interpolation matrix. 
      Since this matrix is symmetric only upper halve has to be considered ---*/

  
    /*--- Looping over the target nodes ---*/
    for(iNode = 0; iNode < interpMatSize; iNode++ ){

      /*--- Looping over the control nodes ---*/
      for (jNode = iNode; jNode < interpMatSize; jNode++){
        
        /*--- Distance between nodes ---*/
        #ifdef HAVE_MPI
          su2double dist(0);
          
          for (unsigned short i = 0; i < nDim; i++) dist += pow(coords[iNode*nDim+i] - coords[jNode*nDim+i], 2);
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

    // SU2_MPI::Send(interpMat.data(), interpMat.size(), MPI_DOUBLE, rank, MASTER_NODE, SU2_MPI::GetComm());

    interpMat.Invert(kernelIsSPD);
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
    
    /*--- define local deformation vector ---*/
    deformationVector_local = deformationVector;

    /*--- sizes for the global and local deformation vectors ---*/
    unsigned long deformationVectorSize, deformationVectorSize_local = deformationVector.size();

    /*--- Obtaining global deformation vector size on master node by summing the local sizes ---*/
    SU2_MPI::Reduce(&deformationVectorSize_local, &deformationVectorSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());

    /*--- Array containing local deformation vector sizes ---*/
    unsigned long defVecSizes[size];

    /*--- Gathering all deformation vectors on the master node ---*/
    if(rank==MASTER_NODE){

      /*--- resizing the global deformation vector ---*/
      deformationVector.resize(deformationVectorSize);
      
      /*--- Gathering the local deformation vector sizes in array defVecSizes ---*/
      SU2_MPI::Gather(&deformationVectorSize_local, 1, MPI_UNSIGNED_LONG, defVecSizes, 1, MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());


      /*--- Receiving the local deformation vector from other processes ---*/
      unsigned long start_idx = 0;
      for (auto iProc = 0; iProc < size; iProc++) {
        if (iProc != MASTER_NODE) {
          SU2_MPI::Recv(&deformationVector[0] + start_idx, defVecSizes[iProc], MPI_DOUBLE, iProc, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE); // TODO can status ignore be used? 
        }
        start_idx += defVecSizes[iProc];
      }      

    }else{

      /*--- Gathering the local deformation vector sizes in array defVecSizes ---*/
      SU2_MPI::Gather(&deformationVectorSize_local, 1, MPI_UNSIGNED_LONG, NULL, 1, MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());

      /*--- Sending the local deformation vector to the master node ---*/
      SU2_MPI::Send(deformationVector_local.data(), deformationVectorSize_local, MPI_DOUBLE, MASTER_NODE, 0, SU2_MPI::GetComm());  
      
    }

    /*--- The global deformation vector is now ordered as d_1, d_2, ..., d_n, where n is the number of processes.
            Here the deformation vector is reordered to obtain an order x_1, ..., x_n, y_1, ..., y_n, z_1, ..., z_n,
                 where x_n is the deformation in x of deformation vector n */
    if(rank == MASTER_NODE){
      
      
      for (unsigned short iDim = nDim-1; iDim > 0; iDim--) {

        unsigned long start_idx = 0;

        for (unsigned short processor = 0; processor < SU2_MPI::GetSize(); processor++) {
          
          if ( processor == 0){
            start_idx += defVecSizes[processor]/nDim;
          }
          else{
            start_idx += (defVecSizes[processor-1]/nDim*(iDim-1) + defVecSizes[processor]/nDim);
          }

          /*--- inserting part of vector at end of deformationVector ---*/
          deformationVector.insert(deformationVector.end(), deformationVector.begin()+start_idx, deformationVector.begin()+start_idx+defVecSizes[processor]/nDim);

          /*--- erasing moved part of the vector ---*/
          deformationVector.erase(deformationVector.begin()+start_idx, deformationVector.begin()+start_idx+defVecSizes[processor]/nDim);
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
  #ifdef HAVE_MPI
    unsigned long size;
    SU2_MPI::Allreduce(&nBoundaryNodes, &size, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
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

  
  unsigned long size; //TODO change variable name, overwrites variable size (nr of parallel processes)
  #ifdef HAVE_MPI

    /*--- Obtaining global nr of control nodes on all processes. ---*/
    SU2_MPI::Allreduce(&nBoundaryNodes, &size, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

    /*--- array containing the global control node coordinates. ---*/
    su2double coords[size*nDim];

    /*--- local size control node coordinates. ---*/
    int local_size = nBoundaryNodes*nDim;
    
    /*--- array containing the local control node coordinates. ---*/
    su2double coords_local[local_size];
    
    /*--- storing local control node coordinates ---*/
    for(iNode = 0; iNode < controlNodes->size(); iNode++){
      auto coord = geometry->nodes->GetCoord((*controlNodes)[iNode]->GetIndex());  
      for ( unsigned short iDim = 0 ; iDim < nDim; iDim++ ){
        coords_local[ iNode * nDim + iDim ] = coord[iDim];
      }
    }

    /*--- array containing size of local control node coordinates. ---*/
    int sizes[SU2_MPI::GetSize()];

    /*--- gathering local control node coordinate sizes on all processes. ---*/
    SU2_MPI::Allgather(&local_size, 1, MPI_INT, sizes, 1, MPI_INT, SU2_MPI::GetComm()); 
    
  
    /*--- array containing the starting indices for the allgatherv operation*/
    int disps[SU2_MPI::GetSize()];    

    for(auto x = 0; x < SU2_MPI::GetSize(); x++){
      if(x == 0){
        disps[x] = 0;
      }else{
        disps[x] = disps[x-1]+sizes[x-1];
      }
    }
    
    

    /*--- making global control node coordinates available on all processes ---*/
    SU2_MPI::Allgatherv(&coords_local, local_size, MPI_DOUBLE, &coords, sizes, disps, MPI_DOUBLE, SU2_MPI::GetComm());
  #else
    size = nBoundaryNodes;
  #endif

  /*--- Vector for storing the coordinate variation ---*/
  su2double var_coord[nDim];
  
  /*--- Loop over the internal nodes ---*/
  for(iNode = 0; iNode < nInternalNodes; iNode++){

    /*--- Loop for contribution of each control node ---*/
    for(cNode = 0; cNode < size; cNode++){
     
      /*--- Determine distance between considered internal and control node ---*/
      su2double dist;

       #ifdef HAVE_MPI
        dist = 0;        
        for (unsigned short i = 0; i < nDim; i++) dist += pow(coords[cNode*nDim+i] - geometry->nodes->GetCoord(internalNodes[iNode])[i], 2);
        dist = sqrt(dist);
      #else
        dist = GeometryToolbox::Distance(nDim, geometry->nodes->GetCoord((*controlNodes)[cNode]->GetIndex()), geometry->nodes->GetCoord(internalNodes[iNode]));
      #endif
      
      /*--- Evaluate RBF based on distance ---*/
      auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(kindRBF, radius, dist));
      
      /*--- Add contribution to total coordinate variation ---*/
      for( iDim = 0; iDim < nDim; iDim++){
        var_coord[iDim] += rbf*coefficients[cNode + iDim*size];
      }
    }

    /*--- Apply the coordinate variation and resetting the var_coord vector to zero ---*/
    for(iDim = 0; iDim < nDim; iDim++){
      geometry->nodes->AddCoord(internalNodes[iNode], iDim, var_coord[iDim]);
      var_coord[iDim] = 0;
    } 
  }  
  
  /*--- Applying the surface deformation, which are stored in the deformation vector ---*/
  for(cNode = 0; cNode < nBoundaryNodes; cNode++){
    if(config->GetMarker_All_Moving((*controlNodes)[cNode]->GetMarker())){
      for(iDim = 0; iDim < nDim; iDim++){
        #ifdef HAVE_MPI
          geometry->nodes->AddCoord((*controlNodes)[cNode]->GetIndex(), iDim, deformationVector_local[cNode + iDim*nBoundaryNodes]); 
        #else
          geometry->nodes->AddCoord((*controlNodes)[cNode]->GetIndex(), iDim, deformationVector[cNode + iDim*nBoundaryNodes]);
        #endif
      }
    }
  }  
}
