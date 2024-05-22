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
  SetInternalNodes(config, geometry); //TODO change order

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
    

      // #ifdef HAVE_MPI

    // SU2_MPI::Barrier(SU2_MPI::GetComm());
    // SU2_MPI::Abort(SU2_MPI::GetComm(), 0);
  // #else
    // std::exit(0);
  // #endif
  

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
  
  vector<unsigned long> vertices;
  vector<unsigned short> markers;
  unsigned short iMarker; 
  unsigned short iVertex; 


  /*--- Storing of the global, marker and vertex indices ---*/
  unsigned long node;
  for(iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    if(!config->GetMarker_All_Deform_Mesh_Internal(iMarker) && !config->GetMarker_All_SendRecv(iMarker)){
      for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
        node = geometry->vertex[iMarker][iVertex]->GetNode();
        if(geometry->nodes->GetDomain(node)){
          boundaryNodes.push_back(new CRadialBasisFunctionNode(node, iMarker, iVertex));        
        }        
      }
    }
  }

  /*--- Sorting of the boundary nodes based on global index ---*/
  sort(boundaryNodes.begin(), boundaryNodes.end(), Compare);

  /*--- Obtaining unique set of boundary nodes ---*/
  boundaryNodes.resize(std::distance(boundaryNodes.begin(), unique(boundaryNodes.begin(), boundaryNodes.end(), Equal)));

  /*--- Updating the number of boundary nodes ---*/
  nBoundaryNodes = boundaryNodes.size();

  // for(auto x : boundaryNodes){cout << rank << " " << x->GetIndex() << " " << geometry->nodes->GetGlobalIndex(x->GetIndex()) << endl;}
}

void CRadialBasisFunctionInterpolation::SetInterpolationMatrix(CGeometry* geometry, CConfig* config){
  auto rank = SU2_MPI::GetRank();
  unsigned long iNode, jNode;

  unsigned long interpMatSize;
  #ifdef HAVE_MPI

    //obtaining the global size of the coordinates
    unsigned long size;
    SU2_MPI::Reduce(&nBoundaryNodes, &size, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
    if(rank != MASTER_NODE){
      size = 0;
    }
    su2double coords[size*nDim];

    su2double coords_local[nBoundaryNodes*nDim];
    
    for(iNode = 0; iNode < controlNodes->size(); iNode++){
      // auto coord = geometry->vertex[(*controlNodes)[iNode]->GetMarker()][(*controlNodes)[iNode]->GetVertex()]->GetCoord();
      auto coord = geometry->nodes->GetCoord((*controlNodes)[iNode]->GetIndex());  
      coords_local[iNode*nDim] = coord[0];
      coords_local[iNode*nDim+1] = coord[1];  
    }


    SU2_MPI::Gather(&coords_local, nBoundaryNodes*nDim, MPI_DOUBLE, &coords, nBoundaryNodes*nDim, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());

    interpMatSize = size;
  #else
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
        // auto dist = GeometryToolbox::Distance(nDim, geometry->nodes->GetCoord((*controlNodes)[iNode]->GetIndex()), geometry->nodes->GetCoord((*controlNodes)[jNode]->GetIndex()));
        #ifdef HAVE_MPI
          su2double dist(0);
          
          for (unsigned short i = 0; i < nDim; i++) dist += pow(coords[iNode*nDim+i] - coords[jNode*nDim+i], 2);
          dist = sqrt(dist);
          
        #else
          // auto dist = GeometryToolbox::Distance(nDim, geometry->vertex[(*controlNodes)[iNode]->GetMarker()][(*controlNodes)[iNode]->GetVertex()]->GetCoord(), geometry->vertex[(*controlNodes)[jNode]->GetMarker()][(*controlNodes)[jNode]->GetVertex()]->GetCoord());
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
  auto rank = SU2_MPI::GetRank();

  /* --- Initialization of the deformation vector ---*/
  deformationVector.resize(controlNodes->size()*nDim, 0.0);

  /*--- If requested (no by default) impose the surface deflections in
    increments and solve the grid deformation with
    successive small deformations. ---*/
  su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());

  /*--- Setting nonzero displacements of the moving markers ---*/
  for(auto i = 0; i < controlNodes->size(); i++){
    
    if(config->GetMarker_All_Moving((*controlNodes)[i]->GetMarker())){
      for(auto iDim = 0; iDim < nDim; iDim++){
        deformationVector[i+iDim*controlNodes->size()] = SU2_TYPE::GetValue(geometry->vertex[(*controlNodes)[i]->GetMarker()][(*controlNodes)[i]->GetVertex()]->GetVarCoord()[iDim] * VarIncrement);
      }
    }
  }  


  

  #ifdef HAVE_MPI
    deformationVector_local = deformationVector;

    unsigned long deformationVectorSize, deformationVectorSize_local = deformationVector.size();

    SU2_MPI::Allreduce(&deformationVectorSize_local, &deformationVectorSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
    
    if(rank==MASTER_NODE){
      auto size = SU2_MPI::GetSize();
      deformationVector.resize(deformationVectorSize);
      unsigned long defVecSizes[size];
      
      SU2_MPI::Gather(&deformationVectorSize_local, 1, MPI_UNSIGNED_LONG, defVecSizes, 1, MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());
      
      int counts_recv[size];
      int displacements[size];

      for(int i = 0; i < size; i++){
        counts_recv[i] = 1;
        if(i == 0){
          displacements[i] = 0;
        }else{
          displacements[i] = displacements[i-1] + defVecSizes[i-1];
        }
      }  
  
    }else{
      SU2_MPI::Gather(&deformationVectorSize_local, 1, MPI_UNSIGNED_LONG, NULL, 1, MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());
    }
    
    SU2_MPI::Gather(deformationVector_local.data(), deformationVectorSize_local, MPI_DOUBLE, deformationVector.data(), deformationVectorSize_local, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());  

    
    if(rank == MASTER_NODE){
      // setting correct order of the deformation vector
      for(unsigned short processor = 0; processor < SU2_MPI::GetSize(); processor++){
        //todo include loop for dimensions
        for(unsigned short iDim = 0; iDim < nDim-1; iDim++){
          unsigned long start_idx = nBoundaryNodes*(processor+1);

          deformationVector.insert(deformationVector.end(), deformationVector.begin()+start_idx, deformationVector.begin()+start_idx+nBoundaryNodes);
          deformationVector.erase(deformationVector.begin()+start_idx, deformationVector.begin()+start_idx+nBoundaryNodes);

        }
      }
    }
  #endif
}

void CRadialBasisFunctionInterpolation::SetInternalNodes(CConfig* config, CGeometry* geometry){ 

  /*--- resizing the internal nodes vector ---*/
  nInternalNodes = geometry->GetnPoint() - nBoundaryNodes; // vector has max size of nPoints
  internalNodes.resize(nInternalNodes);
  

  /*--- Looping over all nodes and check if present in boundary nodes vector ---*/
  unsigned long idx_cnt = 0, idx_control = 0;
  for(unsigned long iNode = 0; iNode < geometry->GetnPoint(); iNode++){    
    if(!geometry->nodes->GetBoundary(iNode) && geometry->nodes->GetDomain(iNode)){
      internalNodes[idx_cnt++] = iNode; 
    }   
  }  

  #ifdef HAVE_MPI
    for(unsigned short iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++){ //TODO cleanup
      if(config->GetMarker_All_SendRecv(iMarker)){ // if send receive marker
        for(unsigned long iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){ 
          if(geometry->nodes->GetDomain(geometry->vertex[iMarker][iVertex]->GetNode())){ 
            auto node = geometry->vertex[iMarker][iVertex]->GetNode();
            if(find_if(boundaryNodes.begin(), boundaryNodes.end(), [&](CRadialBasisFunctionNode* i){return i->GetIndex() == node;}) == boundaryNodes.end()){
              internalNodes[idx_cnt++] = node;
            }
            

          }else{
            auto node = geometry->vertex[iMarker][iVertex]->GetNode();
            internalNodes[idx_cnt++] = node;
          }
        }
      }
    }
  #endif
  nInternalNodes = idx_cnt;
  internalNodes.resize(nInternalNodes);  
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

    cout << endl;  
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

  
  unsigned long size;
  #ifdef HAVE_MPI

    //obtaining the global size of the coordinates
    SU2_MPI::Allreduce(&nBoundaryNodes, &size, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

    su2double coords[size*nDim];

    su2double coords_local[nBoundaryNodes*nDim];
    
    for(iNode = 0; iNode < controlNodes->size(); iNode++){
      // auto coord = geometry->vertex[(*controlNodes)[iNode]->GetMarker()][(*controlNodes)[iNode]->GetVertex()]->GetCoord();
      auto coord = geometry->nodes->GetCoord((*controlNodes)[iNode]->GetIndex());  
      coords_local[iNode*nDim] = coord[0];
      coords_local[iNode*nDim+1] = coord[1];  
    }


    SU2_MPI::Allgather(&coords_local, nBoundaryNodes*nDim, MPI_DOUBLE, &coords, nBoundaryNodes*nDim, MPI_DOUBLE, SU2_MPI::GetComm());
    
  #else
    size = nBoundaryNodes;
  #endif

  // if(rank == MASTER_NODE){
  //   for(auto x = 0; x < geometry->GetnPoint(); x++){
  //     cout << x << " " << geometry->nodes->GetCoord(x)[0] << " " << geometry->nodes->GetCoord(x)[1] << endl;
  //   }
  // }

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

  // if(rank == MASTER_NODE){
  //   for(auto x = 0; x < geometry->GetnPoint(); x++){
  //     cout << x << " " << geometry->nodes->GetCoord(x)[0] << " " << geometry->nodes->GetCoord(x)[1] << endl;
  //   }
  // }
}
