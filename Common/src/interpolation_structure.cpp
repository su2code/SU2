/*!
 * \file interpolation_structure.cpp
 * \brief Main subroutines used by SU2_FSI
 * \author H. Kline
 * \version 4.0.0 "Cardinal"
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


CInterpolator::CInterpolator(CGeometry ***geometry_container, CConfig **config, unsigned int* Zones, unsigned int val_nZone){

  /* Store pointers*/
	Geometry = geometry_container;
	nZone = val_nZone;

  /*--- Initialize transfer coefficients between the zones ---*/
	/* Since this is a virtual function, call it in the child class constructor  */
  //Set_TransferCoeff(iZone_0,iZone_1,config);

  /*--- Initialize Data vectors to 0, by default with length = nDim ---*/
  /* This should be done in the child class ---*/
  //InitializeData(iZone_0,iZone_1,nDim);
	Data = NULL;

}

void CInterpolator::InitializeData(unsigned int* Zones, unsigned short val_nVar){
  nVar=val_nVar;
  unsigned int iZone;
  unsigned short it;
  if (nVar>0){
    /*--- Initialize Data vectors to 0 ---*/
    Data = new su2double**[nZone];
    for (it=0; it<nZone; it++){
      iZone = Zones[it];
      Data[iZone] = new su2double*[Geometry[iZone][MESH_0]->GetnPoint()];
      for (unsigned long iPoint =0; iPoint< Geometry[iZone][MESH_0]->GetnPoint(); iPoint++){
        Data[iZone][iPoint] = new su2double[nVar];
        for (unsigned short iVar=0; iVar<nVar; iVar++){
          Data[iZone][iPoint][iVar]=0.0;
        }
      }
    }
  }else{
    Data = NULL;
  }

}


void CInterpolator::Interpolate_Data(unsigned int iZone, CConfig **config){
  unsigned long iPoint, jPoint, jVertex, iMarker, iVertex;
  unsigned short jMarker;
  unsigned int iZone_1;
  su2double weight=0.0;

  /*--- Loop through points, increment Data in the input zone by the weight in the transfer matrix ---*/

  /*Loop by i then by j to more efficiently call memory*/
  for (iMarker = 0; iMarker < config[iZone]->GetnMarker_All(); iMarker++){
    if (config[iZone]->GetMarker_All_FSIinterface(iMarker) == YES){
      for (iVertex = 0; iVertex<Geometry[iZone][MESH_0]->GetnVertex(iMarker); iVertex++) {
        iPoint =Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetNode();
        /*--- Set Data to 0 before interpolation ---*/
        for (unsigned short iVar=0; iVar<nVar; iVar++){
          Data[iZone][iPoint][iVar]=0;
        }
        /*--- Interpolate ---*/
        for (unsigned short jDonor = 0; jDonor< Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetnDonorPoints(); jDonor++){
          /* Unpack info about the donor point */
          iZone_1 = Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,0);
          jPoint = Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,1);
          jMarker = Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,2);
          jVertex = Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,3);
          weight = Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetDonorCoeff(jDonor);
          /*--- Increment the value of the data ---*/
          for (unsigned short iVar=0; iVar<nVar; iVar++){
            Data[iZone][iPoint][iVar]+=Data[iZone_1][jPoint][iVar]*weight;
          }
        }
      }
    }
  }
}

void CInterpolator::Interpolate_Deformation(unsigned int iZone, CConfig **config){

  unsigned long GlobalIndex, iPoint, i2Point, jPoint, j2Point, iVertex, jVertex;
  unsigned short iMarker, jMarker, iDim;
  unsigned int iZone_1;
  su2double *NewVarCoord = NULL, *VarCoord, *VarRot, *distance = NULL;
  su2double weight;
  unsigned short nDim = Geometry[iZone][MESH_0]->GetnDim();
  /*--- Loop over vertices in the interface marker (zone 0) ---*/
  for (iMarker = 0; iMarker < config[iZone]->GetnMarker_All(); iMarker++){
    if (config[iZone]->GetMarker_All_FSIinterface(iMarker) == YES){
      for (iVertex = 0; iVertex<Geometry[iZone][MESH_0]->GetnVertex(iMarker); iVertex++) {
        iPoint =Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetNode();
        /*--- Set NewCoord to 0 ---*/
        for (iDim=0; iDim<nDim; iDim++) NewVarCoord[iDim]=0.0;
        /*--- Loop over vertices in the interface marker (zone 1) --*/
        for (unsigned short jDonor = 0; jDonor< Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetnDonorPoints(); jDonor++){
          iZone_1 = Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,0);
          jPoint = Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,1);
          jMarker = Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,2);
          jVertex = Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetDonorInfo(jDonor,3);
          weight = Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetDonorCoeff(jDonor);
          /* Get translation and rotation from the solution */
          VarCoord = Geometry[iZone_1][MESH_0]->vertex[jMarker][jVertex]->GetVarCoord();
          VarRot =   Geometry[iZone_1][MESH_0]->vertex[jMarker][jVertex]->GetVarRot();

          for (iDim=0; iDim<nDim; iDim++) distance[iDim]=0.0;

          for (iDim=0; iDim<nDim; iDim++){
            NewVarCoord[iDim]+=VarCoord[iDim]*weight;
            distance[iDim] = Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->GetCoord(iDim)-Geometry[iZone_1][MESH_0]->node[jPoint]->GetCoord(iDim);
          }
          /*--- Add contribution of rotation (cross product of donor point rotation and distance to donor point) ---*/
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
        /*--- Set the varcoord information ---*/
        Geometry[iZone][MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(NewVarCoord);
      }
    }
  }

  // must be called later:
  //flow_grid_movement->SetVolume_Deformation(Geometry[ZONE_0][MESH_0], config[ZONE_0], true);

}

su2double CInterpolator::GetData(unsigned int iZone, unsigned long iPoint, unsigned short iVar){
  if (Data !=NULL)
    return Data[iZone][iPoint][iVar];
  else
    return 0.0; // Check this.
}

su2double* CInterpolator::GetData(unsigned int iZone, unsigned long iPoint){
  if (Data !=NULL)
    return Data[iZone][iPoint];
  else
    return NULL;
}

void CInterpolator::SetData(unsigned int iZone, unsigned long iPoint, unsigned short iVar, su2double val){
  if (Data !=NULL)
    Data[iZone][iPoint][iVar]=val;
  else
    cout <<" CInterpolator object has not been initialized"<<endl;
}


/* Nearest Neighbor Interpolator */
CNearestNeighbor::CNearestNeighbor(CGeometry ***geometry_container, CConfig **config,  unsigned int* Zones,unsigned int nZone) :  CInterpolator(geometry_container, config, Zones,nZone){
  unsigned short nDim = geometry_container[Zones[0]][MESH_0]->GetnDim();
  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(Zones,config);

  /*--- For fluid-structure interaction data interpolated with have nDim dimensions ---*/
  InitializeData(Zones,nDim);

}

CNearestNeighbor::~CNearestNeighbor(){}

void CNearestNeighbor::Set_TransferCoeff(unsigned int* Zones, CConfig **config){
  unsigned long iPoint, jPoint, iVertex, jVertex,*nn;
  unsigned short iMarker, iDim, jMarker;
  unsigned short nDim = Geometry[Zones[0]][MESH_0]->GetnDim(), iDonor, jDonor;
  su2double distance = 0.0, last_distance=-1.0;

  unsigned short int donorindex = 0;
  unsigned short nMarkerFSIint, nMarkerFEA, nMarkerFlow;
  unsigned short iMarkerFSIint, iMarkerFEA, iMarkerFlow;
  unsigned short markFEA, markFlow;

  /*--- Restricted to 2-zone fluid-structure for now ---*/
  unsigned int iZone_0 = Zones[0];
  unsigned int iZone_1 = Zones[1];

  nn = new unsigned long[4];

  /*--- Loop through the vertices in Interface of both zones
   * for Nearest Neighbor each vertex has only one donor point, but for other types of
   * interpolation the number of donor points must be determined first. ---*/

	/*--- Number of markers on the FSI interface ---*/
	nMarkerFSIint = (config[iZone_0]->GetMarker_n_FSIinterface())/2;
  nMarkerFEA  =  config[iZone_1]->GetnMarker_All();
  nMarkerFlow =  config[iZone_0]->GetnMarker_All();

	/*--- For the number of markers on the interface... ---*/
	for (iMarkerFSIint=0; iMarkerFSIint < nMarkerFSIint; iMarkerFSIint++){

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
        for (iDim=0; iDim<nDim; iDim++)
          distance+=pow(Geometry[iZone_1][MESH_0]->vertex[markFEA][jVertex]->GetCoord(iDim)-Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetCoord(iDim),2.0);
        if ((last_distance==-1.0) or (distance<last_distance)){
          last_distance=distance;
          nn[0] = iZone_1; /* Zone of the donor point */
          nn[1] = jPoint; /* global index of the donor point */
          nn[2] = markFEA; /* marker of the donor point */
          nn[3] = jVertex; /* vertex index within marker of the donor point */
        }
		  }
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
      Geometry[iZone_1][MESH_0]->vertex[markFEA][iVertex]->SetDonorInfo(donorindex,nn);
	  for (jVertex = 0; jVertex<Geometry[iZone_0][MESH_0]->GetnVertex(markFlow); jVertex++) {
        jPoint =Geometry[iZone_0][MESH_0]->vertex[markFlow][jVertex]->GetNode();
        distance = 0.0;
        for (iDim=0; iDim<nDim; iDim++)
          distance+=pow(Geometry[iZone_0][MESH_0]->vertex[markFlow][jVertex]->GetCoord(iDim)-Geometry[iZone_1][MESH_0]->vertex[markFEA][iVertex]->GetCoord(iDim),2.0);
        if ((jVertex==0) or (distance<last_distance)){
          last_distance=distance;
          nn[0] = iZone_0; /* Zone of the donor point */
          nn[1] = jPoint; /* global index of the donor point */
          nn[2] = markFlow; /* marker of the donor point */
          nn[3] = jVertex; /* vertex index within marker of the donor point */
        }
      }
      /*--- Set the information of the nearest neighbor ---*/
      Geometry[iZone_1][MESH_0]->vertex[markFEA][iVertex]->SetDonorInfo(donorindex,nn);
      Geometry[iZone_1][MESH_0]->vertex[markFEA][iVertex]->SetDonorCoeff(donorindex,1.0);
    }
	}
}


CConsistConserve::CConsistConserve(CGeometry ***geometry_container, CConfig **config,  unsigned int* Zones,unsigned int nZone) :  CInterpolator(geometry_container, config, Zones,nZone){
  unsigned short nDim = geometry_container[Zones[0]][MESH_0]->GetnDim();
  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(Zones,config);

  /*--- For fluid-structure interaction data interpolated with have nDim dimensions ---*/
  InitializeData(Zones,nDim);
}

CConsistConserve::~CConsistConserve(){}

void CConsistConserve::Set_TransferCoeff(unsigned int* Zones, CConfig **config){
  unsigned long iPoint, jPoint, iVertex, jVertex,*nn, inode, jElem;
  long ivtx;
  unsigned short iMarker, iDim, jMarker, it;
  unsigned short nDim = Geometry[Zones[0]][MESH_0]->GetnDim();
  unsigned short iDonor, jDonor;
  su2double distance = 0.0, last_distance=-1.0, *Coord;
  su2double* myCoeff;
  su2double* myCoefftemp;
  su2double* donorCoord;
  su2double coeff;
  long donor_elem=0, temp_donor;
  unsigned int donor_face;

  unsigned short int donorindex = 0;
  unsigned short nMarkerFSIint, nMarkerFEA, nMarkerFlow;
  unsigned short iMarkerFSIint, iMarkerFEA, iMarkerFlow;
  unsigned short markFEA, markFlow, iFace;
  unsigned short index = 3; // index of the vertex info in the donorinfo array
  unsigned int nNodes;
  /*--- Restricted to 2-zone fluid-structure for now ---*/
  unsigned int iZone_0 = Zones[0];
  unsigned int iZone_1 = Zones[1];
  unsigned int nDonor=0;

  nn = new unsigned long[4];
  /*--- Number of markers on the FSI interface ---*/
  nMarkerFSIint = (config[iZone_0]->GetMarker_n_FSIinterface())/2;
  nMarkerFEA  =  config[iZone_1]->GetnMarker_All();
  nMarkerFlow =  config[iZone_0]->GetnMarker_All();
  /*--- For the number of markers on the interface... ---*/
  for (iMarkerFSIint=0; iMarkerFSIint < nMarkerFSIint; iMarkerFSIint++){
    /*--- Procedure:
     * -Loop through vertices of the aero grid
     * -Find nearest element and allocate enough space in the aero grid donor point info
     *    -set the transfer coefficient values
     *    -increment nDonor for each of the element vertices
     * -Loop through vertices of the structure grid
     *    -Allocate enough space for the donor info
     *    -Loop through the aero vertices and set the donor info at the structure vertices
     */

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
    /*--Same for all points: -*/
    nn[0] = iZone_1; /* Zone of the donor point */
    nn[2] = markFEA; /* marker of the donor point */
    /*--- For the markers on the fluid side ---*/
    /*--- Loop over the vertices on the marker ---*/
    for (iVertex = 0; iVertex<Geometry[iZone_0][MESH_0]->GetnVertex(markFlow); iVertex++) {
      iPoint =Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetNode();
      last_distance=-1.0;
      nDonor = 0;
      /*--- Loop over the vertices in the corresponding interface marker (zone 1), find the closest vertex --*/
      for (jVertex = 0; jVertex<Geometry[iZone_1][MESH_0]->GetnVertex(markFEA); jVertex++) {
        jPoint =Geometry[iZone_1][MESH_0]->vertex[markFEA][jVertex]->GetNode();
        distance = 0.0;
        for (iDim=0; iDim<nDim; iDim++)
          distance+=pow(Geometry[iZone_1][MESH_0]->vertex[markFEA][jVertex]->GetCoord(iDim)-Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetCoord(iDim),2.0);
        if ((last_distance==-1.0) or (distance<last_distance)){
          last_distance=distance;
          /*--- Info stored at the node: zone, point, marker, vertex ---*/
          nn[1] = jPoint; /* global index of the donor point */
        }
      }
      donor_elem=0;
      donor_face=0;
      last_distance=-1;
      /*--- Now that we know the closest vertex, the closest element must be one of the ones connected to the vertex--*/
      for (jElem=0; jElem<Geometry[iZone_1][MESH_0]->node[nn[1]]->GetnElem(); jElem++){
        temp_donor = Geometry[iZone_1][MESH_0]->node[nn[1]]->GetElem(jElem);
        unsigned int nFaces=1; //For 2D cases, we want to look at edges, not faces, as the 'interface'
        if (nDim==3)
          nFaces = Geometry[iZone_1][MESH_0]->elem[temp_donor]->GetnFaces();
        /*--- Loop over all the faces of this element to find one(s) on the interface boundary ---*/
        for (iFace=0; iFace<nFaces; iFace++){
          bool face_on_marker=true;
          /*--- If 3D loop over a face. if 2D loop over an element ---*/
          if (nDim==3){
            nNodes = Geometry[iZone_1][MESH_0]->elem[temp_donor]->GetnNodesFace(iFace);
            /*-- Check if on marker of interface---*/
            for (unsigned int iNode=0; iNode<nNodes; iNode++){
              jPoint = Geometry[iZone_1][MESH_0]->elem[jElem]->GetFaces(iFace, iNode);
              face_on_marker = (face_on_marker && (Geometry[iZone_1][MESH_0]->node[jPoint]->GetVertex(markFEA) !=-1));
            }
          }
          else
            nNodes = Geometry[iZone_1][MESH_0]->elem[temp_donor]->GetnNodes();

          /*--- face_on_marker is true iff all nodes on face iFace are in marker markFEA ---*/
          /*--- if iFace is part of markFEA, calculate the isoparametric coefficients ---*/
          if (face_on_marker){
            /*--- use Isoparametric rep. to find distance to projected point on the surface ---*/
            myCoefftemp = new su2double[nNodes];
            Isoparametric( myCoefftemp, iZone_0,  markFlow, iVertex, nDim, iZone_1, markFEA, temp_donor, iFace, nNodes);
            /*--- If closer than last closest projected point, save. ---*/
            Coord = Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetCoord();
            for (it=0; it< nNodes; it++){
              /*--- If 3D loop over a face. if 2D loop over an element ---*/
              if (nDim==3)
                jPoint = Geometry[iZone_1][MESH_0]->elem[donor_elem]->GetFaces(iFace,it);
              else
                jPoint = Geometry[iZone_1][MESH_0]->elem[donor_elem]->GetNode(it);

              donorCoord = Geometry[iZone_1][MESH_0]->node[jPoint]->GetCoord();
              for (iDim=0; iDim<nDim; iDim++)
                Coord[iDim]-=myCoefftemp[it]*donorCoord[iDim];
            }
            distance = 0;
            for (iDim=0; iDim<nDim; iDim++)
              distance+=Coord[iDim]*Coord[iDim];

            /*--- If the distance is shorter than last closest (and nonzero nodes are on the boundary), update ---*/
            if ((last_distance==-1) or (distance<last_distance )){
              /*--- update last distance ---*/
              last_distance = distance;
              /*--- Store info ---*/
              donor_elem = temp_donor;
              // mem leak?
              myCoeff = new su2double[nNodes];
              for (it=0; it< nNodes; it++){
                myCoeff[it] = myCoefftemp[it];
              }
              Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->SetDonorElem(temp_donor);
              Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->SetDonorFace(iFace);
              Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->SetnDonorPoints(nNodes);
            }
          }
        }
      }

      /*--- If nDonorPoints ==0, no match was found, set nearest neighbor ---*/
      if (Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetnDonorPoints()==0){
        nNodes=1;
        Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->SetnDonorPoints(nNodes);
        donor_elem = -1;
        myCoeff = new su2double[1];
        myCoeff[0] = 1;
      }
      /*--- print the eventual matched point ---*/
      /*
      cout << " Final case for "<< iVertex << " : " << distance <<", " <<
          Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetDonorFace()
          <<", "<<Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetDonorElem()
          <<", iso: ";
      for (it=0; it< nNodes; it++){
         cout << myCoeff[it] << " ";
       }
      cout << endl;
      */
      /*--- Set the appropriate amount of memory ---*/
      Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->Allocate_DonorInfo();
      iFace = Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetDonorFace();
      unsigned int it2=0; // offset introduced to reject points not on the marker
      /*--- Loop over vertices of the element ---*/
      for (it=0; it< Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetnDonorPoints(); it++){
        if (donor_elem!=-1){
          //Important: jPoint is set to the nearest neighbor above. Be careful not to overwrite it.
          if (nDim==3)
            jPoint = Geometry[iZone_1][MESH_0]->elem[donor_elem]->GetFaces(iFace,it);
          else
            jPoint = Geometry[iZone_1][MESH_0]->elem[donor_elem]->GetNode(it);
        }
        else{
          jPoint = nn[1];
        }
        ivtx = Geometry[iZone_1][MESH_0]->node[jPoint]->GetVertex(markFEA);
        if (ivtx!=-1){
          nn[1] = jPoint; /* global index of the donor point */
          nn[3] = ivtx; /* vertex index within marker of the donor point */
          Geometry[iZone_1][MESH_0]->vertex[markFEA][ivtx]->IncrementnDonor();
          Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->SetDonorInfo(it-it2,nn);
          Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->SetDonorCoeff(it-it2,myCoeff[it]);
        }
        else{
          /*---If this node of the element is not a vertex of markFEA GetVertex will return -1
          and we should neglect this point (TODO: make sure this is also taken into account in isoparam)
          Reduce the number of donor points --*/
          Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->SetnDonorPoints(
              Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetnDonorPoints()-1);
          if (abs(myCoeff[it])>1e-12)
            cout <<"Warning: interior pt assigned nonzero coeff."<< endl;
          it2++;
        }
      }
    }
    /*--- Now that all the transfer coefficients have been calculated, loop through the structure vertices
     * and set the same transfer coefficients at the matching points
     */
    index=3;
    for (jVertex = 0; jVertex<Geometry[iZone_1][MESH_0]->GetnVertex(markFEA); jVertex++) {
      ivtx=0;
      Geometry[iZone_1][MESH_0]->vertex[markFEA][jVertex]->Allocate_DonorInfo();
      /*--- Loop over all aero points ---*/
      for (iVertex=0; iVertex<Geometry[iZone_0][MESH_0]->GetnVertex(markFlow); iVertex++){
        /*--- Loop over the donor vertices for iVertex (flow vertex) ---*/
        for (inode=0; inode<Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetnDonorPoints(); inode++){
          /*--- If one of the donor points is the same as the FEA vertex, add information ---*/
          if (Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]-> GetDonorInfo(inode,index) == jVertex){
            iPoint =Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetNode();
            nn[0] = iZone_0; /* Zone of the donor point */
            nn[1] = iPoint; /* global index of the donor point */
            nn[2] = markFlow; /* marker of the donor point */
            nn[3] = iVertex; /* vertex index within marker of the donor point */
            coeff = Geometry[iZone_0][MESH_0]->vertex[markFlow][iVertex]->GetDonorCoeff(inode);
            Geometry[iZone_1][MESH_0]->vertex[markFEA][jVertex]->SetDonorInfo(ivtx,nn);
            Geometry[iZone_1][MESH_0]->vertex[markFEA][jVertex]->SetDonorCoeff(ivtx,coeff);
            ivtx++;
          }
          if (ivtx>=Geometry[iZone_1][MESH_0]->vertex[markFEA][jVertex]->GetnDonorPoints()-1){
            break;
          }
        }
        if (ivtx>=Geometry[iZone_1][MESH_0]->vertex[markFEA][jVertex]->GetnDonorPoints()-1){
          break;
        }

      }
      if (ivtx>=Geometry[iZone_1][MESH_0]->vertex[markFEA][jVertex]->GetnDonorPoints()-1){
        break;
      }
    }
  }
  /*-- Delete locally allocated memory ---*/

  if (myCoeff!=NULL)
    delete[] myCoeff;
  if (myCoefftemp!=NULL)
    delete[] myCoefftemp;

}

void CConsistConserve::Isoparametric(su2double* isoparams, unsigned int iZone_0, unsigned short iMarker, unsigned long iVertex, unsigned int nDim, unsigned int iZone_1, unsigned short jMarker, long donor_elem, unsigned int nDonorPoints, int* temp2){
  int i,j,k;
  int n0 = nDim+1, n;
  su2double tmp, tmp2, distance;
  unsigned long jVertex, inode;

  /*--- Number of neighbor points to interpolate between ---*/
  unsigned int m0 =  nDonorPoints, m;
  su2double x[m0], x_tmp[m0];
  /*--- Q R matrix system ---*/
  su2double Q[m0*m0], R[m0*m0], A[n0*m0];
  bool test[n0];
  bool testi[n0];

  su2double x2[n0];
  int offset;
  /*--- n0: # of dimensions + 1. n: n0 less any degenerate rows ---*/
  n=n0;
  m0=m;

  /*--- Create Matrix A: 1st row all 1's, 2nd row x coordinates, 3rd row y coordinates, etc ---*/
  /*--- Right hand side is [1, \vec{x}']'---*/
  for (i=0; i<m*n; i++)
    A[i]=0;
  for (i=0; i<m; i++)
    A[i]=1.0;

  /*j,n: dimension. i,m: # neighbor point*/
  x[0]=1;
  for (j=1; j<n0; j++)
    x[j]=Geometry[iZone_0][MESH_0]->vertex[iMarker][iVertex]->GetCoord(j-1);

  /*--- temp2 contains node #s that are on the surface (aka, we are excluding interior points) ---*/
  for (k=0; k<m0; k++){
    isoparams[k]=0;
    /*--- If 3D loop over a face. if 2D loop over an element ---*/
    if (nDim==3)
      jPoint = Geometry[iZone_1][MESH_0]->elem[donor_elem]->GetFaces(iFace,k);
    else
      jPoint = Geometry[iZone_1][MESH_0]->elem[donor_elem]->GetNode(k);
    // Neglect donor points that are not members of the matching marker.
    on_marker[k] =(Geometry[iZone_1][MESH_0]->node[jPoint]->GetVertex(jMarker)!=-1);
    if (on_marker[k]){
      // jth coordinate of the ith donor vertex
      for (j=1; j<n0; j++)
        A[j*m+i] = Geometry[iZone_1][MESH_0]->node[jPoint]->GetCoord(j-1);
      i++;
    }
    else
      m--;
  }

  /*--- Eliminate degenerate rows:
   * for example, if z constant including the z values will make the system degenerate
   * TODO: improve efficiency of this loop---*/
  test[0]=true; // always keep the 1st row
  for (i=1; i<n0; i++){
    // Test this row against all previous
    test[i]=true; // Assume that it is not degenerate
    for (k=0; k<i; k++){
      tmp=0; tmp2=0;
      for (j=0;j<m;j++){
        tmp+= A[i*m+j]*A[i*m+j];
        tmp2+=A[k*m+j]*A[k*m+j];
      }
      tmp  = pow(tmp,0.5);
      tmp2 = pow(tmp2,0.5);
      testi[k]=false;
      for (j=0; j<m; j++){
        // If at least one ratio is non-matching row i is not degenerate w/ row k
        if (A[i*m+j]/tmp != A[k*m+j]/tmp2)
          testi[k]=true;
      }
      // If any of testi (k<i) are false, row i is degenerate
      test[i]=(test[i] and testi[k]);
    }
    if (!test[i]) n--;
  }
  /*--- Initialize A2 now that we might have a smaller system --*/
  su2double A2[n*m];
  j=0;
  /*--- Copy only the rows that are non-degenerate ---*/
  for (i=0; i<n0; i++){
    if (test[i]){
      for (k=0;k<m;k++ ){
        A2[m*j+k]=A[m*i+k];
      }
      x2[j]=x[i];
      j++;
    }
  }

  /*--- Initialize Q,R to 0 --*/
  for (i=0; i<m*m; i++){
    Q[i]=0;
    R[i]=0;
  }
  /*--- TODO: make this loop more efficient ---*/
  /*--- Solve for rectangular Q1 R1 ---*/
  for (i=0; i<m; i++){
    tmp=0;
    for (j=0; j<n; j++)
      tmp += (A2[j*m+i])*(A2[j*m+i]);

    R[i*m+i]= pow(tmp,0.5);
    if (tmp>eps && i<n){
      for (j=0; j<n; j++)
        Q[j*m+i]=A2[j*m+i]/R[i*m+i];
    }
    else if (tmp!=0){
      for (j=0; j<n; j++)
        Q[j*m+i]=A2[j*m+i]/tmp;
    }
    for (j=i+1; j<m; j++){
      tmp=0;
      for (k=0; k<n; k++)
        tmp+=A2[k*m+j]*Q[k*m+i];

      R[i*m+j]=tmp;

      for (k=0; k<n; k++)
        A2[k*m+j]=A2[k*m+j]-Q[k*m+i]*R[i*m+j];
    }
  }

  /*--- x_tmp = Q^T * x2 ---*/
  for (i=0; i<m; i++)
    x_tmp[i]=0.0;
  for (i=0; i<m; i++){
    for (j=0; j<n; j++)
      x_tmp[i]+=Q[j*m+i]*x2[j];
  }
  /*--- solve x_tmp = R*isoparams for isoparams: upper triangular system ---*/
  for (i=n-1; i>=0; i--){
    if (R[i*m+i]>eps)
      isoparams[i]=x_tmp[i]/R[i*m+i];
    else
      isoparams[i]=0;
    for (j=0; j<i; j++)
      x_tmp[j]=x_tmp[j]-R[j*m+i]*isoparams[i];
  }

  /*--- Isoparametric coefficients have been calculated. Run checks ---*/

  su2double tol = 1e-13; // hardcoded tolerance

  /*--- Check 1: if close to 0, replace with 0 ---*/
  for (i=0; i<m; i++){
    if (abs(isoparams[i])< tol )
      isoparams[i]=0;
  }

  /*--- Check 2: if > 1, point is ouside face, not really represented accurately ---*/
  bool inside_face = true;
  for (i=0; i<m; i++){
    if (abs(isoparams[i])> 1.1 )
      inside_face = false;
  }
  if (!inside_face){
    /*--- Revert to nearest neighbor ---*/
    tmp=-1; tmp2=0.0; k=0;
    for (i=0; i<m; i++){
      /*--- If 3D loop over a face. if 2D loop over an element ---*/
      if (nDim==3)
        jPoint = Geometry[iZone_1][MESH_0]->elem[donor_elem]->GetFaces(iFace,i);
      else
        jPoint = Geometry[iZone_1][MESH_0]->elem[donor_elem]->GetNode(i);

      for (j=0;j<nDim;j++)
        tmp2+=pow((Geometry[iZone_1][MESH_0]->node[jPoint]->GetCoord(j)-x[j]),2.0);
      if (tmp==-1 or tmp2<tmp){
        tmp=tmp2;
        k=i;
      }
      isoparams[i]=0;
    }
    isoparams[k]=1;
  }
  /*--- Check 3: reorg for neglected points---   */
  i=m-1;
  for (k=m0-1;k>=0;k--){
    if (on_marker[k]){
      isoparams[k] = isoparams[i];
      i--;
    }
    else
      isoparams[k] = 0;
  }
}
