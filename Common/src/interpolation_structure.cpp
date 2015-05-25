/*!
 * \file interpolation_structure.cpp
 * \brief Main subroutines used by SU2_FSI
 * \author H. Kline
 * \version 3.2.9 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/interpolation_structure.hpp"

CInterpolator::CInterpolator(void){

}

CInterpolator::CInterpolator(CGeometry **geometry_container, CConfig **config, unsigned short val_nZone){
	/* Store pointers*/
	Geometry = geometry_container;
	nZone = val_nZone;
	/*Create Transfer Matrix, find nearest neighbors, initialize memory*/
	TransferMatrix->Initialize(Geometry,config);
	/*Set the values of the transfer matrix*/
	Set_TransferMatrix(ZONE_0, ZONE_1);

}


CInterpolator::Interpolate_Force(unsigned short iZone_0, unsigned short iZone_1){

	/*
	  geometry_container[ZONE_0], grid_movement[ZONE_0],
	  config_container[ZONE_0], config_container[ZONE_1],
	  geometry_container[ZONE_1], solver_container[ZONE_1])
	*/
  unsigned long iPoint, jPoint;
  double weight=0.0;
  /*Loop by i then by j to more efficiently call memory*/
  for (iPoint=0; iPoint<Geometry[ZONE_0]->GetnPoint(); iPoint++){
    for (jPoint=0; jPoint<Geometry[ZONE_1]->GetnPoint(); jPoint++){
      /*find the weight stored in the transfer matrix (returns NULL if not connected*/
      weight = TransferMatrix->GetBlock(iPoint,jPoint);
      if (weight!=NULL){
        /* TODO: implement AddForce to increment the force value, and SetForceZero, within geometry class? Add force to the Vertex type? Otherwise need solver container*/
        force_1+=force_0*weight;
      }
    }
  }

}

CInterpolator::Interpolate_Displacement(unsigned short iZone_0, unsigned short iZone_1, CConfig **config){

  unsigned long GlobalIndex, iPoint, i2Point, jPoint, j2Point, iVertex;
  unsigned short iMarker, iDim, nDim;
  double *NewCoord = {0.0,0.0,0.0};
  double weight;
  nDim = Geometry[iZone_0]->GetnDim();
  //GetMarker_All_FSIinterface(iMarker)
  /*--- Loop over points---*/
  for (iPoint=0; iPoint<Geometry[iZone_0]->GetnPoint(); iPoint++){
    /*--- Set NewCoord to 0 ---*/
    for (iDim=0; iDim<nDim; iDim++) NewCoord[iDim]=0.0;
    /*--- Calculate NewCoord Value ---*/
    for (jPoint=0; jPoint<Geometry[iZone_1]->GetnPoint(); jPoint++){
      weight = TransferMatrix->GetBlock(iPoint,jPoint);
      for (iDim=0; iDim<nDim; iDim++) NewCoord[iDim]+=Geometry[iZone_1]->node[jPoint]->GetDiscplacement(iDim)*weight;
      /*---TODO: also add rotation (cross product of rotation and distance)
       *pointlist_a[i].t[0] += w*(dz*pointlist_s[j].r[1] -dy*pointlist_s[j].r[2])
        pointlist_a[i].t[1] += w*(-dz*pointlist_s[j].r[0] +dx*pointlist_s[j].r[2])
        pointlist_a[i].t[2] += w*(dy*pointlist_s[j].r[0] -dx*pointlist_s[j].r[1])
       * ---*/
    }
    /*---Set VarCoord of the vertex type ---*/
    for (iMarker = 0; iMarker < config[iZone_0]->GetMarker_n_FSIinterface(); iMarker++) {
      if (config[iZone_0]->GetMarker_All_FSIinterface(iMarker) == YES) {
        for (iVertex = 0; iVertex < Geometry[iZone_0]->nVertex[iMarker]; iVertex++) {
          i2Point = Geometry[iZone_0]->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex = Geometry[iZone_0]->node[i2Point]->GetGlobalIndex();
          if (GlobalIndex == iPoint) {
            Geometry[iZone_0]->vertex[iMarker][iVertex]->SetVarCoord(NewCoord);
            break;
          }
        }
      }
    }
  }
  // must be called later:
  //flow_grid_movement->SetVolume_Deformation(Geometry[ZONE_0][MESH_0], config[ZONE_0], true);

}

CInterpolator::Set_TransferMatrix(unsigned short iZone_0, unsigned short iZone_1){
  cout<<"base class set transfer matrix"<<endl;

}

CInterpolator::~CInterpolator(void){

}


CNearestNeighbor::Set_TransferMatrix(unsigned short iZone_0, unsigned short iZone_1){
  /*Loop through vertices in FSIinterface, set nearest neighbor value to 1.0, all others to 0.0*/
  cout<<"Nearest neighbor set transfer matrix"<<endl;
}

