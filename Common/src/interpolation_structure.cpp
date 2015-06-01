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
  unsigned short nDim = geometry_container[ZONE_0]->GetnDim();
  /* Store pointers*/
	Geometry = geometry_container;
	nZone = val_nZone;

	/*--- Set matching between zones ---*/
	/*---Create Transfer Matrix, find nearest neighbors, initialize memory---*/
	TransferMatrix->Initialize(Geometry,config);
	/*---Set the values of the transfer matrix---*/
	Set_TransferMatrix(ZONE_0, ZONE_1,config);
	Force = new double**[val_nZone];
	Force[ZONE_0] = new double*[Geometry[ZONE_0]->GetnPoint()];
	Force[ZONE_1] = new double*[Geometry[ZONE_1]->GetnPoint()];

	for (unsigned long iPoint =0; iPoint< Geometry[ZONE_0]->GetnPoint(); iPoint++){
	  Force[ZONE_0][iPoint] = new double[nDim];
	  for (unsigned short iDim=0; iDim<nDim; iDim++){
	    Force[ZONE_0][iPoint][iDim]=0.0;
	  }
	}

  for (unsigned long iPoint =0; iPoint< Geometry[ZONE_1]->GetnPoint(); iPoint++){
    Force[ZONE_1][iPoint] = new double[nDim];
    for (unsigned short iDim=0; iDim<nDim; iDim++){
      Force[ZONE_1][iPoint][iDim]=0.0;
    }
  }




}


CInterpolator::Interpolate_Force(unsigned short iZone_0, unsigned short iZone_1){
  unsigned long iPoint, jPoint;
  unsigned short nDim = Geometry[ZONE_0]->GetnDim();
  double weight=0.0;
  /*Loop by i then by j to more efficiently call memory*/
  for (iPoint=0; iPoint<Geometry[iZone_0]->GetnPoint(); iPoint++){
    for (jPoint=0; jPoint<Geometry[iZone_1]->GetnPoint(); jPoint++){
      /*find the weight stored in the transfer matrix (returns NULL if zero entry)*/
      weight = TransferMatrix->GetBlock(iPoint,jPoint);
      if (weight!=NULL){
        for (unsigned short iDim=0; iDim<nDim; iDim++){
          Force[iZone_1][jPoint][iDim]+=Force[iZone_0][iPoint][iDim]*weight;
        }
      }
    }
  }

}

CInterpolator::Interpolate_Displacement(unsigned short iZone_0, unsigned short iZone_1, CConfig **config){

  unsigned long GlobalIndex, iPoint, i2Point, jPoint, j2Point, iVertex, jVertex;
  unsigned short iMarker, jMarker, iDim;
  double *NewVarCoord = {0.0,0.0,0.0}, *VarCoord, *VarRot, *distance={0.0,0.0,0.0};
  double weight;
  unsigned short nDim = Geometry[iZone_0]->GetnDim();
  /*--- Loop over vertices in the interface marker (zone 0) ---*/
  for (iMarker = 0; iMarker < config[iZone_0]->GetnMarker_All(); iMarker++){
      if (config[iZone_0]->GetMarker_All_FSIinterface(iMarker) == YES){
        for (iVertex = 0; iVertex<Geometry[iZone_0]->GetnVertex(iMarker); iVertex++) {
          iPoint =Geometry[iZone_0]->vertex[iMarker][iVertex]->GetNode();
          /*--- Set NewCoord to 0 ---*/
          for (iDim=0; iDim<nDim; iDim++) NewVarCoord[iDim]=0.0;
          /*--- Loop over vertices in the interface marker (zone 1) --*/
          for (jMarker = 0; jMarker < config[iZone_1]->GetnMarker_All(); jMarker++){
              if (config[iZone_1]->GetMarker_All_FSIinterface(jMarker) == YES){
                for (jVertex = 0; jVertex<Geometry[iZone_1]->GetnVertex(jMarker); jVertex++) {
                  jPoint =Geometry[iZone_1]->vertex[jMarker][jVertex]->GetNode();
                  /*--- Add to the NewCoord value ---*/
                  weight = TransferMatrix->GetBlock(iPoint,jPoint);
                  VarCoord = Geometry[iZone_1]->vertex[jMarker][jVertex]->GetVarCoord();
                  VarRot = Geometry[iZone_1]->vertex[jMarker][jVertex]->GetAuxVar(); // Use Aux var to store rotation vector.
                  distance = {0.0,0.0,0.0};
                  for (iDim=0; iDim<nDim; iDim++){
                    NewVarCoord[iDim]+=VarCoord[iDim]*weight;
                    distance[iDim] = Geometry[iZone_1]->vertex[jMarker][jVertex]->GetCoord(iDim)-Geometry[iZone_0]->vertex[iMarker][iVertex]->GetCoord(iDim);
                  }
                  /*--- Add contribution of rotation ---*/
                  if (nDim==2){
                    VarCoord[0]+=weight*(-distance[1]*VarRot[2]);
                    VarCoord[1]+=weight*(distance[0]*VarRot[2]);
                  }
                  if (nDim==3){
                    VarCoord[0]+=weight*(distance[2]*VarRot[1]-distance[1]*VarRot[2]);
                    VarCoord[1]+=weight*(distance[0]*VarRot[2]-distance[2]*VarRot[0]);
                    VarCoord[2]+=weight*(distance[1]*VarRot[0]-distance[0]*VarRot[1]);
                  }
                }
              }
          }
          Geometry[iZone_0]->vertex[iMarker][iVertex]->SetVarCoord(NewVarCoord);
        }
      }
  }

  // must be called later:
  //flow_grid_movement->SetVolume_Deformation(Geometry[ZONE_0][MESH_0], config[ZONE_0], true);

}

CInterpolator::Set_TransferMatrix(unsigned short iZone_0, unsigned short iZone_1, CConfig **config){
  cout<<"base class set transfer matrix: all zeros, no interpolation will be done."<<endl;
}

CInterpolator::~CInterpolator(void){}

/* Nearest Neighbor Interpolator */
CNearestNeighbor::CNearestNeighbor(CGeometry **geometry_container, CConfig **config, unsigned short nZone) :
    CInterpolator(geometry_container, config, nZone){}

CNearestNeighbor::~CNearestNeighbor(){}

CNearestNeighbor::Set_TransferMatrix(unsigned short iZone_0, unsigned short iZone_1, CConfig **config){
  unsigned long iPoint, jPoint, iVertex, jVertex,nn;
  unsigned short iMarker, iDim, jMarker;
  unsigned short nDim = Geometry[iZone_0]->GetnDim();
  double distance = 0.0, last_distance=-1.0;
  double *val = {1.0};

  /*Loop through vertices in FSIinterface, set nearest neighbor value to 1.0, all others to 0.0*/
    cout<<"Nearest neighbor set transfer matrix"<<endl;
    for (iMarker = 0; iMarker < config[iZone_0]->GetnMarker_All(); iMarker++){
      if (config[iZone_0]->GetMarker_All_FSIinterface(iMarker) == YES){
        for (iVertex = 0; iVertex<Geometry[iZone_0]->GetnVertex(iMarker); iVertex++) {
          iPoint =Geometry[iZone_0]->vertex[iMarker][iVertex]->GetNode();
          last_distance=-1.0;
          /*--- Loop over vertices in the interface marker (zone 1) --*/
          for (jMarker = 0; jMarker < config[iZone_1]->GetnMarker_All(); jMarker++){
            if (config[iZone_1]->GetMarker_All_FSIinterface(jMarker) == YES){
              for (jVertex = 0; jVertex<Geometry[iZone_1]->GetnVertex(jMarker); jVertex++) {
                jPoint =Geometry[iZone_1]->vertex[jMarker][jVertex]->GetNode();
                distance = 0.0;
                for (iDim=0; iDim<nDim;iDim++)
                  distance+=pow(Geometry[iZone_1]->vertex[jMarker][jVertex]->GetCoord(iDim)-Geometry[iZone_0]->vertex[iMarker][iVertex]->GetCoord(iDim),2.0);
                if ((last_distance==-1.0) or (distance<last_distance)){
                  last_distance=distance;
                  nn =jPoint;
                }
              }
            }
          }
          TransferMatrix->SetBlock(iPoint,nn,val);
        }
      }
    }

}

