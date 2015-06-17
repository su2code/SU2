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
  Data = NULL;
}

CInterpolator::~CInterpolator(void){}


CInterpolator::CInterpolator(CGeometry ***geometry_container, CConfig **config, unsigned short val_nZone){
  unsigned short nDim = geometry_container[ZONE_0][MESH_0]->GetnDim();
  /* Store pointers*/
	Geometry = geometry_container;
	nZone = val_nZone;


	/*--- Set matching between zones ---*/
	//TransferMatrix = new CSysTransferMatrix();
	/*---Create Transfer Matrix, find nearest neighbors, initialize memory---*/
	//TransferMatrix->Initialize(Geometry,config);

	/*---Set the values of the transfer matrix---*/
//	Set_TransferMatrix(ZONE_0, ZONE_1,config);

	/*--- Initialize Data vectors to 0 ---*/
	Data = new double**[val_nZone];
	Data[ZONE_0] = new double*[Geometry[ZONE_0][MESH_0]->GetnPoint()];
	Data[ZONE_1] = new double*[Geometry[ZONE_1][MESH_0]->GetnPoint()];

	for (unsigned long iPoint =0; iPoint< Geometry[ZONE_0][MESH_0]->GetnPoint(); iPoint++){
	  Data[ZONE_0][iPoint] = new double[nDim];
	  for (unsigned short iDim=0; iDim<nDim; iDim++){
	    Data[ZONE_0][iPoint][iDim]=0.0;
	  }
	}

  for (unsigned long iPoint =0; iPoint< Geometry[ZONE_1][MESH_0]->GetnPoint(); iPoint++){
    Data[ZONE_1][iPoint] = new double[nDim];
    for (unsigned short iDim=0; iDim<nDim; iDim++){
      Data[ZONE_1][iPoint][iDim]=0.0;
    }
  }




}


void CInterpolator::Interpolate_Data(unsigned short iZone_0, unsigned short iZone_1, CConfig **config){
  unsigned long iPoint, jPoint, jVertex, iMarker, iVertex;
  unsigned short nDim = Geometry[ZONE_0][MESH_0]->GetnDim(), jMarker;
  double weight=0.0;

  /*--- Loop through points, increment Data by the weight in the transfer matrix ---*/

  /*Loop by i then by j to more efficiently call memory*/
  for (iMarker = 0; iMarker < config[iZone_0]->GetnMarker_All(); iMarker++){
      if (config[iZone_0]->GetMarker_All_FSIinterface(iMarker) == YES){
        for (iVertex = 0; iVertex<Geometry[iZone_0][MESH_0]->GetnVertex(iMarker); iVertex++) {
          iPoint =Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
          for (unsigned short jDonor = 0; jDonor< Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetnDonorPoints(); jDonor++){
            /* Unpack info */
            iZone_1 = Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,0);
            jPoint = Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,1);
            jMarker = Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,2);
            jVertex = Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,3);
            weight = Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetDonorCoeff(jDonor);
            for (unsigned short iDim=0; iDim<nDim; iDim++){
              Data[iZone_1][jPoint][iDim]+=Data[iZone_0][iPoint][iDim]*weight;
            }
          }
        }
      }
  }

}

void CInterpolator::Interpolate_Deformation(unsigned short iZone_0, unsigned short iZone_1, CConfig **config){

  unsigned long GlobalIndex, iPoint, i2Point, jPoint, j2Point, iVertex, jVertex;
  unsigned short iMarker, jMarker, iDim;
  double *NewVarCoord = NULL, *VarCoord, *VarRot, *distance = NULL;
  double weight;
  unsigned short nDim = Geometry[iZone_0][MESH_0]->GetnDim();
  /*--- Loop over vertices in the interface marker (zone 0) ---*/
  for (iMarker = 0; iMarker < config[iZone_0]->GetnMarker_All(); iMarker++){
      if (config[iZone_0]->GetMarker_All_FSIinterface(iMarker) == YES){
        for (iVertex = 0; iVertex<Geometry[iZone_0][MESH_0]->GetnVertex(iMarker); iVertex++) {
          iPoint =Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
          /*--- Set NewCoord to 0 ---*/
          for (iDim=0; iDim<nDim; iDim++) NewVarCoord[iDim]=0.0;
          /*--- Loop over vertices in the interface marker (zone 1) --*/
          for (unsigned short jDonor = 0; jDonor< Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetnDonorPoints(); jDonor++){
            iZone_1 = Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,0);
            jPoint = Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,1);
            jMarker = Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,2);
            jVertex = Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,3);
            weight = Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetDonorCoeff(jDonor);
            /* Get translation and rotation from the solution */
            VarCoord = Geometry[iZone_1][MESH_0]->vertex[jMarker][jVertex]->GetVarCoord();
            // This is a fix so it compiles... But it still needs to be developed.
            VarRot[0] = Geometry[iZone_1][MESH_0]->vertex[jMarker][jVertex]->GetAuxVar(); // Use Aux var to store rotation vector.
            VarRot[1] = Geometry[iZone_1][MESH_0]->vertex[jMarker][jVertex]->GetAuxVar(); // Use Aux var to store rotation vector.
            VarRot[2] = Geometry[iZone_1][MESH_0]->vertex[jMarker][jVertex]->GetAuxVar(); // Use Aux var to store rotation vector.
            distance[0] = 0.0;
            distance[1] = 0.0;
            distance[2] = 0.0;
            for (iDim=0; iDim<nDim; iDim++){
              NewVarCoord[iDim]+=VarCoord[iDim]*weight;
              distance[iDim] = Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetCoord(iDim)-Geometry[iZone_1][MESH_0]->node[jPoint]->GetCoord(iDim);
            }
            /*--- Add contribution of rotation ---*/
            if (nDim==2){
              NewVarCoord[0]+=weight*(-distance[1]*VarRot[2]);
              NewVarCoord[1]+=weight*(distance[0]*VarRot[2]);
            }
            if (nDim==3){
              NewVarCoord[0]+=weight*(distance[2]*VarRot[1]-distance[1]*VarRot[2]);
              NewVarCoord[1]+=weight*(distance[0]*VarRot[2]-distance[2]*VarRot[0]);
              NewVarCoord[2]+=weight*(distance[1]*VarRot[0]-distance[0]*VarRot[1]);
            }
          }
          // Or introduce deformation vector that stores this.
          Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(NewVarCoord);
        }
      }
  }

  // must be called later:
  //flow_grid_movement->SetVolume_Deformation(Geometry[ZONE_0][MESH_0], config[ZONE_0], true);

}

//void CInterpolator::Interpolate_Solution( unsigned short iZone_dest, CConfig **config, CSolver **solver_container){
//  unsigned long iPoint, jPoint, jVertex, iMarker, iVertex;
//  unsigned short jMarker, iZone_source, nVar = solver_container[iZone_dest]->GetnVar();
//  double weight=0.0, dest_val=0.0, src_val=0.0;
//  /*--- Loop through the interface vertices in the destination zone ---*/
//  for (iMarker = 0; iMarker < config[iZone_dest]->GetnMarker_All(); iMarker++){
//    if (config[iZone_dest]->GetMarker_All_FSIinterface(iMarker) == YES){
//      for (iVertex = 0; iVertex<Geometry[iZone_dest]->GetnVertex(iMarker); iVertex++) {
//        /*--- Set the values at the interface point to 0 initially ---*/
//        iPoint =Geometry[iZone_dest]->vertex[iMarker][iVertex]->GetNode();
//        for (unsigned short iVar=0; iVar<nVar; iVar++)
//          solver_container[iZone_dest]->node[iPoint]->SetSolution(iVar,0.0);
//        /*--- Loop through donor points ---*/
//        for (unsigned short jDonor = 0; jDonor< Geometry[iZone_dest]->vertex[iMarker][iVertex]->GetnDonorPoints(); jDonor++){
//          /*--- unpack Donor Point info ---*/
//          iZone_source = Geometry[iZone_dest]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,0);
//          jPoint = Geometry[iZone_dest]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,1);
//          weight = Geometry[iZone_dest]->vertex[iMarker][iVertex]->GetDonorCoeff(jDonor);
//          for (unsigned short iVar=0; iVar<nVar; iVar++){
//            /*--- Increment the value of the solution ---*/
//            dest_val=solver_container[iZone_dest]->node[iPoint]->GetSolution(iVar);
//            dest_val+=solver_container[iZone_source]->node[jPoint]->GetSolution(iVar)*weight;
//            /*--- Set the value in the solution container ---*/
//            solver_container[iZone_dest]->node[iPoint]->SetSolution(iVar,dest_val);
//          }
//        }
//      }
//    }
//  }
//}


void CInterpolator::Set_TransferMatrix(unsigned short iZone_0, unsigned short iZone_1, CConfig **config){
  cout<<"base class set transfer matrix: all zeros, no interpolation will be done."<<endl;
}

double CInterpolator::GetData(unsigned short iZone, unsigned long iPoint, unsigned short iDim){
  if (Data !=NULL)
    return Data[iZone][iPoint][iDim];
  else
    return 0.0; // Check this.
}

double* CInterpolator::GetData(unsigned short iZone, unsigned long iPoint){
  if (Data !=NULL)
    return Data[iZone][iPoint];
  else
    return NULL;
}

void CInterpolator::SetData(unsigned short iZone, unsigned long iPoint, unsigned short iDim, double val){
  if (Data !=NULL)
    Data[iZone][iPoint][iDim]=val;
  else
    cout <<" CInterpolator object has not been initialized"<<endl;
}


/* Nearest Neighbor Interpolator */
CNearestNeighbor::CNearestNeighbor(CGeometry ***geometry_container, CConfig **config, unsigned short nZone) :
    CInterpolator(geometry_container, config, nZone){

	Set_TransferMatrix(ZONE_0, ZONE_1,config);

}

CNearestNeighbor::~CNearestNeighbor(){}

void CNearestNeighbor::Set_TransferMatrix(unsigned short iZone_0, unsigned short iZone_1, CConfig **config){
  unsigned long iPoint, jPoint, iVertex, jVertex,*nn;
  unsigned short iMarker, iDim, jMarker;
  unsigned short nDim = Geometry[iZone_0][MESH_0]->GetnDim(), iDonor, jDonor;
  double distance = 0.0, last_distance=-1.0;
//  double *val = 1.0;
  unsigned short int donorindex = 0;
  unsigned short nMarkerFSIint, nMarkerFEA, nMarkerFlow;
  unsigned short iMarkerFSIint, iMarkerFEA, iMarkerFlow;
  unsigned short markFEA, markFlow;

  nn = new unsigned long[4];
  /*--- Loop through the vertices in Interface of both zones
   * for Nearest Neighbor each vertex has only one donor point, but for other types of
   * interpolation the number of donor points must be determined first. ---*/

	/*--- Number of markers on the FSI interface ---*/
	nMarkerFSIint = (config[iZone_0]->GetMarker_n_FSIinterface())/2;

	/*--- For the number of markers on the interface... ---*/
	for (iMarkerFSIint=0; iMarkerFSIint < nMarkerFSIint; iMarkerFSIint++){

		nMarkerFEA  =  config[iZone_1]->GetnMarker_All();
		nMarkerFlow =  config[iZone_0]->GetnMarker_All();

		/*--- ... the marker markFEA ... ---*/
		for (iMarkerFEA=0; iMarkerFEA < nMarkerFEA; iMarkerFEA++){
			if ( config[iZone_1]->GetMarker_All_FSIinterface(iMarkerFEA) == (iMarkerFSIint+1)){
				markFEA=iMarkerFEA;
			}
		}
		/*--- ... corresponds to the marker markFlow. ---*/
		for (iMarkerFlow=0; iMarkerFlow < nMarkerFlow; iMarkerFlow++){
			if (config[iZone_0]->GetMarker_All_FSIinterface(iMarkerFlow) == (iMarkerFSIint+1)){
				markFlow=iMarkerFlow;
			}
		}

		/*--- For the markers on the fluid side ---*/
		/*--- Loop over the vertices on the marker ---*/
        for (iVertex = 0; iVertex<Geometry[iZone_0][MESH_0]->GetnVertex(markFlow); iVertex++) {
          iPoint =Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetNode();
          last_distance=-1.0;
          /*--- Allocate memory with known number of donor points (1 for nearest neighbor) ---*/
          Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->SetnDonorPoints(1);
          Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->Allocate_DonorInfo();
          /*--- Loop over the vertices in the corresponding interface marker (zone 1) --*/

		  for (jVertex = 0; jVertex<Geometry[iZone_1][MESH_0]->GetnVertex(markFEA); jVertex++) {
			jPoint =Geometry[iZone_1][MESH_0]->vertex[markFEA][jVertex]->GetNode();
			distance = 0.0;
			for (iDim=0; iDim<nDim;iDim++)
			  distance+=pow(Geometry[iZone_1][MESH_0]->vertex[markFEA][jVertex]->GetCoord(iDim)-Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetCoord(iDim),2.0);
			if ((last_distance==-1.0) or (distance<last_distance)){
			  last_distance=distance;
			  //nn ={iZone_1, jPoint,jMarker,jVertex};
			  nn[0] = iZone_1;
			  nn[1] = jPoint;
			  nn[2] = markFEA;
			  nn[3] = jVertex;
			}
		  }

          /*--- Set the information of the nearest neighbor ---*/
		  /*--- Enable this to check that we are doing it fine ---*/
//     	  cout << "The distance from the vertex " << iVertex << " in the Flow marker " << markFlow << " to the vertex " << nn[3] << " in the FEA marker " << markFEA << " is " << last_distance << endl;
          /*--- Check what donorindex has to be... Vertex? Node? ---*/
          Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->SetDonorInfo(donorindex,nn);
          Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->SetDonorCoeff(donorindex,1.0);
        }

		/*--- For the marker on the FEA side ---*/
		/*--- Loop over the vertices on the marker ---*/
        for (iVertex = 0; iVertex<Geometry[iZone_1][MESH_0]->GetnVertex(markFEA); iVertex++) {
          iPoint =Geometry[iZone_1][MESH_0]->vertex[markFEA][iVertex]->GetNode();
          last_distance=-1.0;
          /*--- Allocate memory with known number of donor points (1 for nearest neighbor) ---*/
          Geometry[iZone_1][MESH_0]->vertex[markFEA][iVertex]->SetnDonorPoints(1);
          Geometry[iZone_1][MESH_0]->vertex[markFEA][iVertex]->Allocate_DonorInfo();

          /*--- Loop over vertices in the interface marker (zone 0) --*/

		  for (jVertex = 0; jVertex<Geometry[iZone_0][MESH_0]->GetnVertex(markFlow); jVertex++) {
			jPoint =Geometry[iZone_0][MESH_0]->vertex[markFlow][jVertex]->GetNode();
			distance = 0.0;
			for (iDim=0; iDim<nDim;iDim++)
			  distance+=pow(Geometry[iZone_0][MESH_0]->vertex[markFlow][jVertex]->GetCoord(iDim)-Geometry[iZone_1][MESH_0]->vertex[markFEA][iVertex]->GetCoord(iDim),2.0);
			if ((last_distance==-1.0) or (distance<last_distance)){
			  last_distance=distance;
			  //nn ={iZone_1, jPoint,jMarker,jVertex};
			  nn[0] = iZone_1;
			  nn[1] = jPoint;
			  nn[2] = jMarker;
			  nn[3] = jVertex;
			}
		  }

          /*--- Set the information of the nearest neighbor ---*/
		  /*--- Enable this to check that we are doing it fine ---*/
//     	  cout << "The distance from the vertex " << iVertex << " in the FEA marker " << markFEA << " to the vertex " << nn[3] << " in the Flow marker " << markFlow << " is " << last_distance << endl;
          /*--- Check what donorindex has to be... Vertex? Node? ---*/
          Geometry[iZone_1][MESH_0]->vertex[markFEA][iVertex]->SetDonorInfo(donorindex,nn);
          Geometry[iZone_1][MESH_0]->vertex[markFEA][iVertex]->SetDonorCoeff(donorindex,1.0);
        }

	}

//    for (iMarker = 0; iMarker < config[iZone_0]->GetnMarker_All(); iMarker++){
//      if (config[iZone_0]->GetMarker_All_FSIinterface(iMarker) == YES){
//        for (iVertex = 0; iVertex<Geometry[iZone_0][MESH_0]->GetnVertex(iMarker); iVertex++) {
//          iPoint =Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
//          last_distance=-1.0;
//          /*--- Allocate memory with known number of donor points (1 for nearest neighbor) ---*/
//          Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->SetnDonorPoints(1);
//          Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->Allocate_DonorInfo();
//          /*--- Loop over vertices in the interface marker (zone 1) --*/
//          for (jMarker = 0; jMarker < config[iZone_1]->GetnMarker_All(); jMarker++){
//            if (config[iZone_1]->GetMarker_All_FSIinterface(jMarker) == YES){
//              for (jVertex = 0; jVertex<Geometry[iZone_1][MESH_0]->GetnVertex(jMarker); jVertex++) {
//                jPoint =Geometry[iZone_1][MESH_0]->vertex[jMarker][jVertex]->GetNode();
//                distance = 0.0;
//                for (iDim=0; iDim<nDim;iDim++)
//                  distance+=pow(Geometry[iZone_1][MESH_0]->vertex[jMarker][jVertex]->GetCoord(iDim)-Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetCoord(iDim),2.0);
//                if ((last_distance==-1.0) or (distance<last_distance)){
//                  last_distance=distance;
//                  //nn ={iZone_1, jPoint,jMarker,jVertex};
//                  nn[0] = iZone_1;
//                  nn[1] = jPoint;
//                  nn[2] = jMarker;
//                  nn[3] = jVertex;
//                }
//              }
//            }
//          }
//          /*--- Set the information of the nearest neighbor ---*/
//          /*--- Check what donorindex has to be... Vertex? Node? ---*/
//     	  cout << "Distance node " << iMarker << " to " << nn[2] << " is " << last_distance << endl;
//          Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->SetDonorInfo(donorindex,nn);
//          Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->SetDonorCoeff(donorindex,1.0);
//        }
//      }
//    }
//    /*--- Do the same for the next zone ---*/
//    for (iMarker = 0; iMarker < config[iZone_1]->GetnMarker_All(); iMarker++){
//      if (config[iZone_1]->GetMarker_All_FSIinterface(iMarker) == YES){
//        for (iVertex = 0; iVertex<Geometry[iZone_1][MESH_0]->GetnVertex(iMarker); iVertex++) {
//          iPoint =Geometry[iZone_1][MESH_0]->vertex[iMarker][iVertex]->GetNode();
//          last_distance=-1.0;
//          /*--- Allocate memory with known number of donor points (1 for nearest neighbor) ---*/
//          Geometry[iZone_1][MESH_0]->vertex[iMarker][iVertex]->SetnDonorPoints(1);
//          Geometry[iZone_1][MESH_0]->vertex[iMarker][iVertex]->Allocate_DonorInfo();
//          /*--- Loop over vertices in the interface marker (zone 1) --*/
//          for (jMarker = 0; jMarker < config[iZone_1]->GetnMarker_All(); jMarker++){
//            if (config[iZone_1]->GetMarker_All_FSIinterface(jMarker) == YES){
//              for (jVertex = 0; jVertex<Geometry[iZone_0][MESH_0]->GetnVertex(jMarker); jVertex++) {
//                jPoint =Geometry[iZone_0][MESH_0]->vertex[jMarker][jVertex]->GetNode();
//                distance = 0.0;
//                for (iDim=0; iDim<nDim;iDim++)
//                  distance+=pow(Geometry[iZone_0][MESH_0]->vertex[jMarker][jVertex]->GetCoord(iDim)-Geometry[iZone_1][MESH_0]->vertex[iMarker][iVertex]->GetCoord(iDim),2.0);
//                if ((last_distance==-1.0) or (distance<last_distance)){
//                  last_distance=distance;
//                  //nn ={iZone_1, jPoint,jMarker,jVertex};
//                  nn[0] = iZone_1;
//                  nn[1] = jPoint;
//                  nn[2] = jMarker;
//                  nn[3] = jVertex;
//                }
//              }
//            }
//          }
//          /*--- Set the information of the nearest neighbor ---*/
//          /*--- Check what donorindex has to be... Vertex? Node? ---*/
//     	  cout << "Distance node " << iMarker << " to " << nn[2] << " is " << last_distance << endl;
//          Geometry[iZone_1][MESH_0]->vertex[iMarker][iVertex]->SetDonorInfo(donorindex,nn);
//          Geometry[iZone_1][MESH_0]->vertex[iMarker][iVertex]->SetDonorCoeff(donorindex,1.0);
//        }
//      }
//    }

}

