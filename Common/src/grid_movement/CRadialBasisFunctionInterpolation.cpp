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
  /*--- Retrieving number of deformation steps. ---*/
  auto Nonlinear_Iter = config->GetGridDef_Nonlinear_Iter();

  /*--- Assigning the node types ---*/
  SetControlNodes(geometry, config);
  SetInternalNodes(geometry);

  /*--- Looping over the number of deformation iterations ---*/
  for (auto iNonlinear_Iter = 0ul; iNonlinear_Iter < Nonlinear_Iter; iNonlinear_Iter++) {

    /*--- Obtaining the interpolation coefficients of the control nodes ---*/
    GetInterpolationCoefficients(geometry, config, iNonlinear_Iter);
    
    /*--- Updating the coordinates of the grid ---*/
    UpdateGridCoord(geometry, config);
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
  unsigned short iVertex; 


  /*--- Total number of boundary nodes (including duplicates of shared boundaries) ---*/
  unsigned long nBoundNodes = 0;
  for(iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
      nBoundNodes += geometry->nVertex[iMarker];
  }

  /*--- Vector with boudary nodes has at most nBoundNodes ---*/
  boundaryNodes.resize(nBoundNodes);

  /*--- Storing of the global, marker and vertex indices ---*/
  unsigned long count = 0;
  for(iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
      boundaryNodes[count+iVertex] = new CRadialBasisFunctionNode(geometry, iMarker, iVertex);
      }
    count += geometry->nVertex[iMarker];
  }

  /*--- Sorting of the boundary nodes based on global index ---*/
  sort(boundaryNodes.begin(), boundaryNodes.end(), Compare);

  /*--- Obtaining unique set of boundary nodes ---*/
  boundaryNodes.resize(std::distance(boundaryNodes.begin(), unique(boundaryNodes.begin(), boundaryNodes.end(), Equal)));
  
  /*--- Updating the number of boundary nodes ---*/
  nBoundaryNodes = boundaryNodes.size();
  
}

void CRadialBasisFunctionInterpolation::SetInterpolationMatrix(CGeometry* geometry, CConfig* config){
  unsigned long iNode, jNode;

  /*--- Initialization of the interpolation matrix ---*/
  interpMat.Initialize(controlNodes->size());

  /*--- Construction of the interpolation matrix. 
    Since this matrix is symmetric only upper halve has to be considered ---*/

  /*--- Looping over the target nodes ---*/
  for(iNode = 0; iNode < controlNodes->size(); iNode++ ){
    /*--- Looping over the control nodes ---*/
    for (jNode = iNode; jNode < controlNodes->size(); jNode++){

      /*--- Distance between nodes ---*/
      auto dist = GeometryToolbox::Distance(nDim, geometry->nodes->GetCoord((*controlNodes)[iNode]->GetIndex()), geometry->nodes->GetCoord((*controlNodes)[jNode]->GetIndex()));
      
      /*--- Evaluation of RBF ---*/
      interpMat(iNode, jNode) = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(kindRBF, radius, dist));
    }
  }

  /*--- Obtaining lower halve using symmetry ---*/
  const bool kernelIsSPD = (kindRBF == RADIAL_BASIS::WENDLAND_C2) || (kindRBF == RADIAL_BASIS::GAUSSIAN) ||
                          (kindRBF == RADIAL_BASIS::INV_MULTI_QUADRIC);
  interpMat.Invert(kernelIsSPD);

}

void CRadialBasisFunctionInterpolation::SetDeformationVector(CGeometry* geometry, CConfig* config){
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
}

void CRadialBasisFunctionInterpolation::SetInternalNodes(CGeometry* geometry){

  /*--- resizing the internal nodes vector ---*/
  nInternalNodes = geometry->GetnPoint() - nBoundaryNodes;
  internalNodes.resize(nInternalNodes);

  /*--- Looping over all nodes and check if present in boundary nodes vector ---*/
  unsigned long idx_cnt = 0, idx_control = 0;
  for(unsigned long iNode = 0; iNode < geometry->GetnPoint(); iNode++){
    
    /*--- If iNode is equal to boundaryNodes[idx_control] 
      then this node is a boundary node and idx_control can be updated ---*/
    if(idx_control < nBoundaryNodes && iNode == boundaryNodes[idx_control]->GetIndex()){idx_control++;}
    
    /*--- If not equal then the node is an internal node ---*/
    else{
      internalNodes[idx_cnt] = iNode;
      idx_cnt++;
    }   
  }  
}


void CRadialBasisFunctionInterpolation::SolveRBF_System(){

  /*--- resizing the interpolation coefficient vector ---*/
  coefficients.resize(nDim*controlNodes->size());

  /*--- Looping through the dimensions in order to find the interpolation coefficients for each direction ---*/
  unsigned short iDim;
  for(iDim = 0; iDim < nDim; iDim++){
    interpMat.MatVecMult(deformationVector.begin()+iDim*controlNodes->size(), coefficients.begin()+iDim*controlNodes->size());
  }
}

void CRadialBasisFunctionInterpolation::UpdateGridCoord(CGeometry* geometry, CConfig* config){
  unsigned long iNode, cNode;
  unsigned short iDim;
  
  /*--- Vector for storing the coordinate variation ---*/
  su2double var_coord[nDim];
  
  /*--- Loop over the internal nodes ---*/
  for(iNode = 0; iNode < nInternalNodes; iNode++){

    /*--- Loop for contribution of each control node ---*/
    for(cNode = 0; cNode < controlNodes->size(); cNode++){
      
      /*--- Determine distance between considered internal and control node ---*/
      auto dist = GeometryToolbox::Distance(nDim, geometry->nodes->GetCoord((*controlNodes)[cNode]->GetIndex()), geometry->nodes->GetCoord(internalNodes[iNode]));

      /*--- Evaluate RBF based on distance ---*/
      auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(kindRBF, radius, dist));

      /*--- Add contribution to total coordinate variation ---*/
      for( iDim = 0; iDim < nDim; iDim++){
        var_coord[iDim] += rbf*coefficients[cNode + iDim*controlNodes->size()];
      }
    }

    /*--- Apply the coordinate variation and resetting the var_coord vector to zero ---*/
    for(iDim = 0; iDim < nDim; iDim++){
      geometry->nodes->AddCoord(internalNodes[iNode], iDim, var_coord[iDim]);
      var_coord[iDim] = 0;
    } 
  }  
  
  /*--- Applying the surface deformation, which are stored in the deformation vector ---*/
  for(cNode = 0; cNode < controlNodes->size(); cNode++){
    if(config->GetMarker_All_Moving((*controlNodes)[cNode]->GetMarker())){
      for(iDim = 0; iDim < nDim; iDim++){
        geometry->nodes->AddCoord((*controlNodes)[cNode]->GetIndex(), iDim, deformationVector[cNode + iDim*controlNodes->size()]);
      }
    }
  }  

}
