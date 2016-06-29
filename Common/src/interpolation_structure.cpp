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
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

	nZone = 0;
	Geometry = NULL;

	donor_geometry  = NULL;
	target_geometry = NULL;

	donorZone  = 0;
	targetZone = 0;

  Buffer_Receive_nVertex_Donor		= NULL;
  Buffer_Receive_nFace_Donor		= NULL;
  Buffer_Receive_nFaceNodes_Donor	= NULL;
  Buffer_Send_nVertex_Donor			= NULL;
  Buffer_Send_nFace_Donor			= NULL;
  Buffer_Send_nFaceNodes_Donor		= NULL;
  Buffer_Receive_GlobalPoint		= NULL;
  Buffer_Send_GlobalPoint			= NULL;
  Buffer_Send_FaceIndex				= NULL;
  Buffer_Receive_FaceIndex			= NULL;
  Buffer_Send_FaceNodes				= NULL;
  Buffer_Receive_FaceNodes			= NULL;
  Buffer_Send_FaceProc				= NULL;
  Buffer_Receive_FaceProc			= NULL;

  Buffer_Send_Coord					= NULL;
  Buffer_Send_Normal				= NULL;
  Buffer_Receive_Coord				= NULL;
  Buffer_Receive_Normal				= NULL;

}

CInterpolator::~CInterpolator(void){

  //if (Buffer_Receive_nVertex_Donor!=NULL) delete[] Buffer_Receive_nVertex_Donor;
}


CInterpolator::CInterpolator(CGeometry ***geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone){

  /* Store pointers*/
	Geometry = geometry_container;

	//nZone = 2;

	donorZone  = iZone;
	targetZone = jZone;

	donor_geometry  = geometry_container[donorZone][MESH_0];
	target_geometry = geometry_container[targetZone][MESH_0];

  /*--- Initialize transfer coefficients between the zones ---*/
	/* Since this is a virtual function, call it in the child class constructor  */
  //Set_TransferCoeff(targetZone,donorZone,config);
  /*--- Initialize transfer coefficients between the zones ---*/
  //Set_TransferCoeff(Zones,config);

  //Buffer_Receive_nVertex_Donor=NULL;

}

inline void CInterpolator::Set_TransferCoeff(CConfig **config) { }

void CInterpolator::Determine_ArraySize(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim){
	
  unsigned long nLocalVertex_Donor = 0, nLocalFaceNodes_Donor=0, nLocalFace_Donor=0;
  unsigned long iVertex, iPointDonor = 0;
  /* Only needed if face data is also collected */
  unsigned long inode;
  unsigned long donor_elem, jElem, jPoint;
  unsigned short iDonor;
  unsigned int nFaces=0, iFace, nNodes=0;
  bool face_on_marker = true;

#ifdef HAVE_MPI
  int rank = MASTER_NODE;
  int nProcessor = SINGLE_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

  for (iVertex = 0; iVertex < nVertexDonor; iVertex++) {
    iPointDonor = donor_geometry->vertex[markDonor][iVertex]->GetNode();
    if (donor_geometry->node[iPointDonor]->GetDomain()){
      nLocalVertex_Donor++;
      if (faces){
        /*--- On Donor geometry also communicate face info ---*/
        if (nDim == 3){
          for (jElem = 0; jElem < donor_geometry->node[iPointDonor]->GetnElem(); jElem++){
            donor_elem = donor_geometry->node[iPointDonor]->GetElem(jElem);
            nFaces = donor_geometry->elem[donor_elem]->GetnFaces();
            for (iFace = 0; iFace < nFaces; iFace++){
              face_on_marker = true;
              nNodes = donor_geometry->elem[donor_elem]->GetnNodesFace(iFace);
              for (iDonor = 0; iDonor < nNodes; iDonor++){
                /*--- Local index of the node on face --*/
                inode = donor_geometry->elem[donor_elem]->GetFaces(iFace, iDonor);
                jPoint = donor_geometry->elem[donor_elem]->GetNode(inode);

                face_on_marker = (face_on_marker && (donor_geometry->node[jPoint]->GetVertex(markDonor) !=-1));
              }
              if ( face_on_marker ){
                nLocalFace_Donor++;
                nLocalFaceNodes_Donor += nNodes;
              }
            }
          }
        }
        else{
          /*--- in 2D we use the edges --- */
          nNodes = 2;
          nFaces = donor_geometry->node[iPointDonor]->GetnPoint();
          for (iFace = 0; iFace < nFaces; iFace++){
            face_on_marker = true;
            for (iDonor = 0; iDonor < nNodes; iDonor++){
              inode = donor_geometry->node[iPointDonor]->GetEdge(iFace);
              jPoint = donor_geometry->edge[inode]->GetNode(iDonor);

              face_on_marker = (face_on_marker && (donor_geometry->node[jPoint]->GetVertex(markDonor) !=-1));
            }
            if (face_on_marker){
              nLocalFace_Donor++;
              nLocalFaceNodes_Donor += nNodes;
            }
          }
        }
      }
    }
  }

  Buffer_Send_nVertex_Donor[0] = nLocalVertex_Donor;
  if (faces){
    Buffer_Send_nFace_Donor[0] = nLocalFace_Donor;
    Buffer_Send_nFaceNodes_Donor[0] = nLocalFaceNodes_Donor;
  }

  /*--- Send Interface vertex information --*/
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalVertex_Donor, &MaxLocalVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_nVertex_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  if (faces){
    SU2_MPI::Allreduce(&nLocalFace_Donor,      &nGlobalFace_Donor,      1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFace_Donor,      &MaxFace_Donor,          1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &nGlobalFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &MaxFaceNodes_Donor,     1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFace_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MaxFace_Donor++;
  }
#else
  MaxLocalVertex_Donor    = nLocalVertex_Donor;
  Buffer_Receive_nVertex_Donor[0] = Buffer_Send_nVertex_Donor[0];
  if (faces){
    nGlobalFace_Donor       = nLocalFace_Donor;
    nGlobalFaceNodes_Donor  = nLocalFaceNodes_Donor;
    MaxFaceNodes_Donor      = nLocalFaceNodes_Donor;
    MaxFace_Donor           = nLocalFace_Donor+1;
    Buffer_Receive_nFace_Donor[0] = Buffer_Send_nFace_Donor[0];
    Buffer_Receive_nFaceNodes_Donor[0] = Buffer_Send_nFaceNodes_Donor[0];
  }
#endif

}

void CInterpolator::Collect_VertexInfo(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim){
  unsigned long nLocalVertex_Donor = 0;
  unsigned long iVertex, iPointDonor = 0, iVertexDonor, nBuffer_Coord, nBuffer_Point;
  unsigned short iDim;
  /* Only needed if face data is also collected */
  su2double  *Normal;

#ifdef HAVE_MPI
  int rank = MASTER_NODE;
  int nProcessor = SINGLE_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif


  for (iVertex = 0; iVertex < MaxLocalVertex_Donor; iVertex++) {
    Buffer_Send_GlobalPoint[iVertex] = 0;
    for (iDim = 0; iDim < nDim; iDim++){
      Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
      if (faces)
        Buffer_Send_Normal[iVertex*nDim+iDim] = 0.0;
    }
  }

  /*--- Copy coordinates and point to the auxiliar vector --*/
  nLocalVertex_Donor = 0;

  for (iVertexDonor = 0; iVertexDonor < nVertexDonor; iVertexDonor++) {
    iPointDonor = donor_geometry->vertex[markDonor][iVertexDonor]->GetNode();
    if (donor_geometry->node[iPointDonor]->GetDomain()) {
      Buffer_Send_GlobalPoint[nLocalVertex_Donor] = donor_geometry->node[iPointDonor]->GetGlobalIndex();
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[nLocalVertex_Donor*nDim+iDim] = donor_geometry->node[iPointDonor]->GetCoord(iDim);

      if (faces){
        Normal =  donor_geometry->vertex[markDonor][iVertexDonor]->GetNormal();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Normal[nLocalVertex_Donor*nDim+iDim] = Normal[iDim];
      }
      nLocalVertex_Donor++;
    }
  }
  nBuffer_Coord = MaxLocalVertex_Donor*nDim;
  nBuffer_Point = MaxLocalVertex_Donor;

#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_GlobalPoint, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  if (faces){
    SU2_MPI::Allgather(Buffer_Send_Normal, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Normal, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
  }
#else
  for (iVertex = 0; iVertex < nBuffer_Coord; iVertex++)
    Buffer_Receive_Coord[iVertex] = Buffer_Send_Coord[iVertex];

  for (iVertex = 0; iVertex < nBuffer_Point; iVertex++)
    Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];

  if (faces){
    for (iVertex = 0; iVertex < nBuffer_Coord; iVertex++)
      Buffer_Receive_Normal[iVertex] = Buffer_Send_Normal[iVertex];
  }
#endif
}

int CInterpolator::Find_InterfaceMarker(CConfig *config, unsigned short val_marker_interface){
	
	unsigned short nMarker = config->GetnMarker_All();
	unsigned short iMarker;
	
	for (iMarker = 0; iMarker < nMarker; iMarker++){
		
		  /*--- If the tag GetMarker_All_FSIinterface(iMarker) equals the index we are looping at ---*/
		  if (config->GetMarker_All_FSIinterface(iMarker) == val_marker_interface ){
			  /*--- We have identified the identifier for the interface marker ---*/
			  return iMarker;
		  }
	}
	
	return -1;
}


/* Nearest Neighbor Interpolator */
CNearestNeighbor::CNearestNeighbor(void):  CInterpolator(){ }

CNearestNeighbor::CNearestNeighbor(CGeometry ***geometry_container, CConfig **config,  unsigned int iZone, unsigned int jZone) :  CInterpolator(geometry_container, config, iZone, jZone){

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

}

CNearestNeighbor::~CNearestNeighbor(){}

void CNearestNeighbor::Set_TransferCoeff(CConfig **config){

  unsigned long jVertex;
  unsigned short iDim;
  unsigned short nDim = donor_geometry->GetnDim();

  unsigned short nMarkerInt, nMarkerDonor, nMarkerTarget;		// Number of markers on the interface, donor and target side
  unsigned short iMarkerInt;		// Variables for iteration over markers
  int markDonor = -1, markTarget = -1;

  unsigned long nVertexDonor = 0, nVertexTarget= 0;
  unsigned long Point_Target = 0;

  unsigned long iVertexTarget;
  unsigned long pGlobalPoint = 0;
  int iProcessor, pProcessor = 0;

  unsigned short int iDonor = 0;
  unsigned long Global_Point_Donor;

  /*--- Number of markers on the FSI interface ---*/
  nMarkerInt     = (int) (config[donorZone]->GetMarker_n_FSIinterface() /2);
  nMarkerTarget  = target_geometry->GetnMarker();
  nMarkerDonor   = donor_geometry->GetnMarker();

  su2double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist;

  int nProcessor = SINGLE_NODE;

#ifdef HAVE_MPI
  int rank = MASTER_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

  Buffer_Receive_nVertex_Donor = new unsigned long [nProcessor];

  // For the markers on the interface
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {

	nVertexDonor = 0;
	nVertexTarget= 0;

	/*--- On the donor side ---*/
	markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);
	  
	/*--- On the target side ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

	if(markTarget == -1 || markDonor == -1)
		continue;
		
	nVertexDonor  =  donor_geometry->GetnVertex(markDonor);
	nVertexTarget = target_geometry->GetnVertex(markTarget);

	Buffer_Send_nVertex_Donor = new unsigned long [1];
    Buffer_Receive_nVertex_Donor = new unsigned long [nProcessor];

	  /* Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor */
	  Determine_ArraySize(false, markDonor,markTarget,nVertexDonor,nDim);

    Buffer_Send_Coord = new su2double [MaxLocalVertex_Donor*nDim];
    Buffer_Send_GlobalPoint = new unsigned long [MaxLocalVertex_Donor];

    Buffer_Receive_Coord = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
    Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Donor];

	  /*-- Collect coordinates, global points, and normal vectors ---*/
	  Collect_VertexInfo(false, markDonor,markTarget,nVertexDonor,nDim);

	  /*--- Compute the closest point to a Near-Field boundary point ---*/
	  maxdist = 0.0;
	  for (iVertexTarget = 0; iVertexTarget < nVertexTarget; iVertexTarget++) {

		  Point_Target = target_geometry->vertex[markTarget][iVertexTarget]->GetNode();

		  if (target_geometry->node[Point_Target]->GetDomain()) {

			  target_geometry->vertex[markTarget][iVertexTarget]->SetnDonorPoints(1);
			  target_geometry->vertex[markTarget][iVertexTarget]->Allocate_DonorInfo();

			  /*--- Coordinates of the boundary point ---*/
			  Coord_i = target_geometry->node[Point_Target]->GetCoord();
			  mindist = 1E6; pProcessor = 0;

			  /*--- Loop over all the boundaries to find the pair ---*/
			  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){
				  for (jVertex = 0; jVertex < nVertexDonor; jVertex++) {
					  Global_Point_Donor = iProcessor*MaxLocalVertex_Donor+jVertex;
						
					  /*--- Compute the dist ---*/
					  dist = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
						  Coord_j[iDim] = Buffer_Receive_Coord[ Global_Point_Donor*nDim+iDim];
						  dist += pow(Coord_j[iDim]-Coord_i[iDim],2.0);
					  }

					  if (dist < mindist) {
						  mindist = dist; pProcessor = iProcessor; pGlobalPoint = Buffer_Receive_GlobalPoint[Global_Point_Donor];
					  }

					  if (dist == 0.0) break;
				  }
			  }

			  /*--- Store the value of the pair ---*/
			  maxdist = max(maxdist, mindist);
			  target_geometry->vertex[markTarget][iVertexTarget]->SetInterpDonorPoint(iDonor, pGlobalPoint);
			  target_geometry->vertex[markTarget][iVertexTarget]->SetInterpDonorProcessor(iDonor, pProcessor);
			  target_geometry->vertex[markTarget][iVertexTarget]->SetDonorCoeff(iDonor,1.0);
	      //cout <<" Nearest Neighbor for target g.i " << target_geometry->node[Point_Target]->GetGlobalIndex() <<" is "<< pGlobalPoint << "; d = " << mindist<< endl;
	     }
	  }

	  delete[] Buffer_Send_Coord;
	  delete[] Buffer_Send_GlobalPoint;

	  delete[] Buffer_Receive_Coord;
	  delete[] Buffer_Receive_GlobalPoint;

	  delete[] Buffer_Send_nVertex_Donor;
	  delete[] Buffer_Receive_nVertex_Donor;


  }

}



CIsoparametric::CIsoparametric(CGeometry ***geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone)  :  CInterpolator(geometry_container, config, iZone, jZone){

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

  /*--- For fluid-structure interaction data interpolated with have nDim dimensions ---*/
 // InitializeData(Zones,nDim);
}

CIsoparametric::~CIsoparametric(){}

void CIsoparametric::Set_TransferCoeff(CConfig **config){
  unsigned long iVertex, jVertex;
  unsigned long  dPoint, inode, jElem, nElem;
  unsigned short iDim, iDonor=0, iFace;

  unsigned short nDim = donor_geometry->GetnDim();

  unsigned short nMarkerInt, nMarkerDonor, nMarkerTarget;
  unsigned short iMarkerInt;

  int markDonor=0, markTarget=0;

  long donor_elem=0, temp_donor=0;
  unsigned int nNodes=0;
  /*--- Restricted to 2-zone for now ---*/
  unsigned int nFaces=1; //For 2D cases, we want to look at edges, not faces, as the 'interface'
  bool face_on_marker=true;

  unsigned long nVertexDonor = 0, nVertexTarget= 0;
  unsigned long Point_Target = 0;

  unsigned long iVertexDonor, iPointDonor = 0;
  unsigned long jGlobalPoint = 0;
  int iProcessor;

  unsigned long nLocalFace_Donor = 0, nLocalFaceNodes_Donor=0;

  unsigned long faceindex;

  su2double dist = 0.0, mindist=1E6, *Coord, *Coord_i;
  su2double myCoeff[10]; // Maximum # of donor points
  su2double  *Normal;
  su2double *projected_point = new su2double[nDim];
  su2double tmp, tmp2;
  su2double storeCoeff[10];
  unsigned long storeGlobal[10];
  int storeProc[10];

  int rank = MASTER_NODE;
  int nProcessor = SINGLE_NODE;
  Coord = new su2double[nDim];
  Normal = new su2double[nDim];

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

  /*--- Number of markers on the FSI interface ---*/
  nMarkerInt     = (int)(config[donorZone]->GetMarker_n_FSIinterface())/2;
  nMarkerTarget  = target_geometry->GetnMarker();
  nMarkerDonor   = donor_geometry->GetnMarker();


  /*--- For the number of markers on the interface... ---*/
  for (iMarkerInt=1; iMarkerInt <= nMarkerInt; iMarkerInt++){
    /*--- Procedure:
     * -Loop through vertices of the aero grid
     * -Find nearest element and allocate enough space in the aero grid donor point info
     *    -set the transfer coefficient values
     */
    nVertexDonor = 0;
	nVertexTarget= 0;

	/*--- On the donor side ---*/
	markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);
	  
	/*--- On the target side ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

	if(markTarget == -1 || markDonor == -1)
		continue;

	nVertexDonor  = donor_geometry->GetnVertex(markDonor);
	nVertexTarget = target_geometry->GetnVertex(markTarget);
	
    Buffer_Send_nVertex_Donor	 = new unsigned long [1];
    Buffer_Send_nFace_Donor		 = new unsigned long [1];
    Buffer_Send_nFaceNodes_Donor = new unsigned long [1];

    Buffer_Receive_nVertex_Donor 	= new unsigned long [nProcessor];
    Buffer_Receive_nFace_Donor 		= new unsigned long [nProcessor];
    Buffer_Receive_nFaceNodes_Donor = new unsigned long [nProcessor];

    /* Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor */
    Determine_ArraySize(true, markDonor, markTarget, nVertexDonor, nDim);

    Buffer_Send_Coord 		= new su2double [MaxLocalVertex_Donor*nDim];
    Buffer_Send_Normal		= new su2double [MaxLocalVertex_Donor*nDim];
    Buffer_Send_GlobalPoint = new unsigned long [MaxLocalVertex_Donor];

    Buffer_Receive_Coord 	   = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
    Buffer_Receive_Normal 	   = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
    Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Donor];

    /*-- Collect coordinates, global points, and normal vectors ---*/
    Collect_VertexInfo(true, markDonor, markTarget, nVertexDonor, nDim);

    Buffer_Send_FaceIndex   = new unsigned long[MaxFace_Donor];
    Buffer_Send_FaceNodes   = new unsigned long[MaxFaceNodes_Donor];
    Buffer_Send_FaceProc    = new unsigned long[MaxFaceNodes_Donor];

    Buffer_Receive_FaceIndex = new unsigned long[MaxFace_Donor*nProcessor];
    Buffer_Receive_FaceNodes = new unsigned long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_FaceProc  = new unsigned long[MaxFaceNodes_Donor*nProcessor];

    nLocalFace_Donor 	  = 0;
    nLocalFaceNodes_Donor = 0;

    /*--- Collect Face info ---*/
    for (iVertex=0; iVertex < MaxFace_Donor; iVertex++){
      Buffer_Send_FaceIndex[iVertex] = 0;
    }
    for (iVertex=0; iVertex < MaxFaceNodes_Donor; iVertex++){
      Buffer_Send_FaceNodes[iVertex] = 0;
      Buffer_Send_FaceProc[iVertex]  = 0;
    }

    Buffer_Send_FaceIndex[0] = rank * MaxFaceNodes_Donor;

    if (nDim==2) nNodes=2;

    for (iVertexDonor = 0; iVertexDonor < nVertexDonor; iVertexDonor++) {
      iPointDonor = donor_geometry->vertex[markDonor][iVertexDonor]->GetNode();

      if (donor_geometry->node[iPointDonor]->GetDomain()) {

        if (nDim==3)  nElem = donor_geometry->node[iPointDonor]->GetnElem();
        else          nElem =donor_geometry->node[iPointDonor]->GetnPoint();

        for (jElem=0; jElem < nElem; jElem++){
          if (nDim==3){
            temp_donor = donor_geometry->node[iPointDonor]->GetElem(jElem);
            nFaces = donor_geometry->elem[temp_donor]->GetnFaces();
            for (iFace = 0; iFace < nFaces; iFace++){
              /*-- Determine whether this face/edge is on the marker --*/
              face_on_marker = true;
              nNodes = donor_geometry->elem[temp_donor]->GetnNodesFace(iFace);
              for (iDonor = 0; iDonor < nNodes; iDonor++){
                inode = donor_geometry->elem[temp_donor]->GetFaces(iFace, iDonor);
                dPoint = donor_geometry->elem[temp_donor]->GetNode(inode);
                face_on_marker = (face_on_marker && (donor_geometry->node[dPoint]->GetVertex(markDonor) !=-1));
              }

              if (face_on_marker ){
                for (iDonor = 0; iDonor < nNodes; iDonor++){
                  inode = donor_geometry->elem[temp_donor]->GetFaces(iFace, iDonor);
                  dPoint = donor_geometry->elem[temp_donor]->GetNode(inode);
                  // Match node on the face to the correct global index
                  jGlobalPoint = donor_geometry->node[dPoint]->GetGlobalIndex();
                  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){
                    for (jVertex = 0; jVertex < Buffer_Receive_nVertex_Donor[iProcessor]; jVertex++) {
                      if (jGlobalPoint == Buffer_Receive_GlobalPoint[  MaxLocalVertex_Donor*iProcessor+jVertex ]){
                        Buffer_Send_FaceNodes[nLocalFaceNodes_Donor] = MaxLocalVertex_Donor*iProcessor+jVertex;
                        Buffer_Send_FaceProc[nLocalFaceNodes_Donor] =  iProcessor;
                      }
                    }
                  }
                  nLocalFaceNodes_Donor++; // Increment total number of face-nodes / processor
                }
                /* Store the indices */
                Buffer_Send_FaceIndex[nLocalFace_Donor+1] = Buffer_Send_FaceIndex[nLocalFace_Donor] + nNodes;
                nLocalFace_Donor++; // Increment number of faces / processor
              }
            }
          }
          else{
            /*-- Determine whether this face/edge is on the marker --*/
            face_on_marker = true;
            for (iDonor=0; iDonor < nNodes; iDonor++){
              inode = donor_geometry->node[iPointDonor]->GetEdge(jElem);
              dPoint = donor_geometry->edge[inode]->GetNode(iDonor);

              face_on_marker = (face_on_marker && (donor_geometry->node[dPoint]->GetVertex(markDonor) !=-1));
            }
            if (face_on_marker ){
              for (iDonor=0; iDonor < nNodes; iDonor++){
                inode = donor_geometry->node[iPointDonor]->GetEdge(jElem);
                dPoint = donor_geometry->edge[inode]->GetNode(iDonor);
                // Match node on the face to the correct global index
                jGlobalPoint=donor_geometry->node[dPoint]->GetGlobalIndex();
                for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){
                  for (jVertex = 0; jVertex < Buffer_Receive_nVertex_Donor[iProcessor]; jVertex++) {
                    if (jGlobalPoint ==Buffer_Receive_GlobalPoint[MaxLocalVertex_Donor*iProcessor+jVertex]){
                      Buffer_Send_FaceNodes[nLocalFaceNodes_Donor]=MaxLocalVertex_Donor*iProcessor+jVertex;
                      Buffer_Send_FaceProc[nLocalFaceNodes_Donor]=iProcessor;
                    }
                  }
                }
                nLocalFaceNodes_Donor++; // Increment total number of face-nodes / processor
              }
              /* Store the indices */
              Buffer_Send_FaceIndex[nLocalFace_Donor+1] = Buffer_Send_FaceIndex[nLocalFace_Donor]+nNodes;
              nLocalFace_Donor++; // Increment number of faces / processor
            }
          }
        }
      }
    }

    //Buffer_Send_FaceIndex[nLocalFace_Donor+1] = MaxFaceNodes_Donor*rank+nLocalFaceNodes_Donor;
#ifdef HAVE_MPI
    SU2_MPI::Allgather(Buffer_Send_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_FaceProc,  MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceProc,  MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_FaceIndex, MaxFace_Donor,      MPI_UNSIGNED_LONG, Buffer_Receive_FaceIndex, MaxFace_Donor,      MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
    for (iFace = 0; iFace < MaxFace_Donor; iFace++){
      Buffer_Receive_FaceIndex[iFace] = Buffer_Send_FaceIndex[iFace];
    }
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertex++)
      Buffer_Receive_FaceNodes[iVertex] = Buffer_Send_FaceNodes[iVertex];
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertex++)
      Buffer_Receive_FaceProc[iVertex] = Buffer_Send_FaceProc[iVertex];
#endif

    /*--- Loop over the vertices on the target Marker ---*/
    for (iVertex = 0; iVertex < nVertexTarget; iVertex++) {
      mindist=1E6;
      for (unsigned short iCoeff = 0; iCoeff < 10; iCoeff++){
        storeCoeff[iCoeff] = 0;
      }
      Point_Target = target_geometry->vertex[markTarget][iVertex]->GetNode();

      if (target_geometry->node[Point_Target]->GetDomain()) {

        Coord_i = target_geometry->node[Point_Target]->GetCoord();
        /*---Loop over the faces previously communicated/stored ---*/
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){

          nFaces = (unsigned int)Buffer_Receive_nFace_Donor[iProcessor];

          for (iFace = 0; iFace < nFaces; iFace++){
            /*--- ---*/

            nNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1] -
                    (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace];

            su2double *X = new su2double[nNodes*nDim];
            faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
            for (iDonor = 0; iDonor < nNodes; iDonor++){
              jVertex = Buffer_Receive_FaceNodes[iDonor+faceindex]; // index which points to the stored coordinates, global points
              for (iDim=0; iDim < nDim; iDim++){
                X[iDim*nNodes+iDonor]=
                    Buffer_Receive_Coord[jVertex*nDim+iDim];
              }
            }
            jVertex = Buffer_Receive_FaceNodes[faceindex];

            for (iDim=0; iDim<nDim; iDim++){
              Normal[iDim] = Buffer_Receive_Normal[jVertex*nDim+iDim];
            }

            /* Project point used for case where surfaces are not exactly coincident, where
             * the point is assumed connected by a rigid rod normal to the surface.
             */
            tmp  = 0;
            tmp2 = 0;
            for (iDim=0; iDim < nDim; iDim++){
              tmp  += Normal[iDim]*Normal[iDim];
              tmp2 += Normal[iDim]*(Coord_i[iDim] - X[iDim*nNodes]);
            }
            tmp  = 1/tmp;
            tmp2 = tmp2*sqrt(tmp);
            for (iDim = 0; iDim < nDim; iDim++){
              // projection of \vec{q} onto plane defined by \vec{n} and \vec{p}:
              // \vec{q} - \vec{n} ( (\vec{q}-\vec{p} ) \cdot \vec{n})
              // tmp2 = ( (\vec{q}-\vec{p} ) \cdot \vec{N})
              // \vec{n} = \vec{N}/(|N|), tmp = 1/|N|^2
              projected_point[iDim] = Coord_i[iDim] + Normal[iDim]*tmp2*tmp;
            }

            Isoparameters(nDim, nNodes, X, projected_point, myCoeff);

            /*--- Find distance to the interpolated point ---*/
            dist = 0.0;
            for (iDim = 0; iDim < nDim; iDim++){
              Coord[iDim] = Coord_i[iDim];
              for(iDonor=0; iDonor< nNodes; iDonor++){
                Coord[iDim] -= myCoeff[iDonor]*X[iDim*nNodes+iDonor];
              }
              dist += pow(Coord[iDim],2.0);
            }

            /*--- If the dist is shorter than last closest (and nonzero nodes are on the boundary), update ---*/
            if (dist < mindist ){
              /*--- update last dist ---*/
              mindist = dist;
              /*--- Store info ---*/
              donor_elem = temp_donor;
              target_geometry->vertex[markTarget][iVertex]->SetDonorElem(donor_elem); // in 2D is nearest neighbor
              target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nNodes);
              for (iDonor=0; iDonor<nNodes; iDonor++){
                storeCoeff[iDonor] = myCoeff[iDonor];
                jVertex = Buffer_Receive_FaceNodes[faceindex+iDonor];
                storeGlobal[iDonor] =Buffer_Receive_GlobalPoint[jVertex];
                storeProc[iDonor] = (int)Buffer_Receive_FaceProc[faceindex+iDonor];
              }
            }
          
            delete [] X;
          }
        }
        /*--- Set the appropriate amount of memory and fill ---*/
        nNodes =target_geometry->vertex[markTarget][iVertex]->GetnDonorPoints();
        target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();

        for (iDonor=0; iDonor<nNodes; iDonor++){
          target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint(iDonor,storeGlobal[iDonor]);
          //cout <<rank << " Global Point " << Global_Point<<" iDonor " << iDonor <<" coeff " << coeff <<" gp " << pGlobalPoint << endl;
          target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(iDonor,storeCoeff[iDonor]);
          target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(iDonor, storeProc[iDonor]);
        }
      }
    }

    delete[] Buffer_Send_nVertex_Donor;
    delete[] Buffer_Send_nFace_Donor;
    delete[] Buffer_Send_nFaceNodes_Donor;

    delete[] Buffer_Receive_nVertex_Donor;
    delete[] Buffer_Receive_nFace_Donor;
    delete[] Buffer_Receive_nFaceNodes_Donor;

    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_Normal;
    delete[] Buffer_Send_GlobalPoint;

    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_Normal;
    delete[] Buffer_Receive_GlobalPoint;

    delete[] Buffer_Send_FaceIndex;
    delete[] Buffer_Send_FaceNodes;
    delete[] Buffer_Send_FaceProc;

    delete[] Buffer_Receive_FaceIndex;
    delete[] Buffer_Receive_FaceNodes;
    delete[] Buffer_Receive_FaceProc;
  }
  delete [] Coord;
  delete [] Normal;
  
  delete [] projected_point;
}

void CIsoparametric::Isoparameters(unsigned short nDim, unsigned short nDonor,
    su2double *X, su2double *xj, su2double *isoparams){
  short iDonor,iDim,k; // indices
  su2double tmp, tmp2;
  
  su2double *x     = new su2double[nDim+1];
  su2double *x_tmp = new su2double[nDim+1];
  su2double *Q     = new su2double[nDonor*nDonor];
  su2double *R     = new su2double[nDonor*nDonor];
  su2double *A     = new su2double[nDim+1*nDonor];
  su2double *A2    = NULL;
  su2double *x2    = new su2double[nDim+1];
  
  bool *test  = new bool[nDim+1];
  bool *testi = new bool[nDim+1];
  
  su2double eps = 1E-10;
  
  short n = nDim+1;

  if (nDonor>2){
    /*--- Create Matrix A: 1st row all 1's, 2nd row x coordinates, 3rd row y coordinates, etc ---*/
    /*--- Right hand side is [1, \vec{x}']'---*/
    for (iDonor=0; iDonor<nDonor; iDonor++){
      isoparams[iDonor]=0;
      A[iDonor]=1.0;
      for (iDim=0; iDim<n; iDim++)
        A[(iDim+1)*nDonor+iDonor]=X[iDim*nDonor+iDonor];
    }

    x[0]=1.0;
    for (iDim=0; iDim<nDim; iDim++)
      x[iDim+1]=xj[iDim];

    /*--- Eliminate degenerate rows:
     * for example, if z constant including the z values will make the system degenerate
     * TODO: improve efficiency of this loop---*/
    test[0]=true; // always keep the 1st row
    for (iDim=1; iDim<nDim+1; iDim++){
      // Test this row against all previous
      test[iDim]=true; // Assume that it is not degenerate
      for (k=0; k<iDim; k++){
        tmp=0; tmp2=0;
        for (iDonor=0;iDonor<nDonor;iDonor++){
          tmp+= A[iDim*nDonor+iDonor]*A[iDim*nDonor+iDonor];
          tmp2+=A[k*nDonor+iDonor]*A[k*nDonor+iDonor];
        }
        tmp  = pow(tmp,0.5);
        tmp2 = pow(tmp2,0.5);
        testi[k]=false;
        for (iDonor=0; iDonor<nDonor; iDonor++){
          // If at least one ratio is non-matching row iDim is not degenerate w/ row k
          if (A[iDim*nDonor+iDonor]/tmp != A[k*nDonor+iDonor]/tmp2)
            testi[k]=true;
        }
        // If any of testi (k<iDim) are false, row iDim is degenerate
        test[iDim]=(test[iDim] && testi[k]);
      }
      if (!test[iDim]) n--;
    }

    /*--- Initialize A2 now that we might have a smaller system --*/
    A2 = new su2double[n*nDonor];
    iDim=0;
    /*--- Copy only the rows that are non-degenerate ---*/
    for (k=0; k<nDim+1; k++){
      if (test[k]){
        for (iDonor=0;iDonor<nDonor;iDonor++ ){
          A2[nDonor*iDim+iDonor]=A[nDonor*k+iDonor];
        }
        x2[iDim]=x[k];
        iDim++;
      }
    }
    /*--- Initialize Q,R to 0 --*/
    for (k=0; k<nDonor*nDonor; k++){
      Q[k]=0;
      R[k]=0;
    }
    /*--- TODO: make this loop more efficient ---*/
    /*--- Solve for rectangular Q1 R1 ---*/
    for (iDonor=0; iDonor<nDonor; iDonor++){
      tmp=0;
      for (iDim=0; iDim<n; iDim++)
        tmp += (A2[iDim*nDonor+iDonor])*(A2[iDim*nDonor+iDonor]);

      R[iDonor*nDonor+iDonor]= pow(tmp,0.5);
      if (tmp>eps && iDonor<n){
        for (iDim=0; iDim<n; iDim++)
          Q[iDim*nDonor+iDonor]=A2[iDim*nDonor+iDonor]/R[iDonor*nDonor+iDonor];
      }
      else if (tmp!=0){
        for (iDim=0; iDim<n; iDim++)
          Q[iDim*nDonor+iDonor]=A2[iDim*nDonor+iDonor]/tmp;
      }
      for (iDim=iDonor+1; iDim<nDonor; iDim++){
        tmp=0;
        for (k=0; k<n; k++)
          tmp+=A2[k*nDonor+iDim]*Q[k*nDonor+iDonor];

        R[iDonor*nDonor+iDim]=tmp;

        for (k=0; k<n; k++)
          A2[k*nDonor+iDim]=A2[k*nDonor+iDim]-Q[k*nDonor+iDonor]*R[iDonor*nDonor+iDim];
      }
    }
    /*--- x_tmp = Q^T * x2 ---*/
    for (iDonor=0; iDonor<nDonor; iDonor++)
      x_tmp[iDonor]=0.0;
    for (iDonor=0; iDonor<nDonor; iDonor++){
      for (iDim=0; iDim<n; iDim++)
        x_tmp[iDonor]+=Q[iDim*nDonor+iDonor]*x2[iDim];
    }

    /*--- solve x_tmp = R*isoparams for isoparams: upper triangular system ---*/
    for (iDonor=n-1; iDonor>=0; iDonor--){
      if (R[iDonor*nDonor+iDonor]>eps)
        isoparams[iDonor]=x_tmp[iDonor]/R[iDonor*nDonor+iDonor];
      else
        isoparams[iDonor]=0;
      for (k=0; k<iDonor; k++)
        x_tmp[k]=x_tmp[k]-R[k*nDonor+iDonor]*isoparams[iDonor];
    }
  }
  else{
    /*-- For 2-donors (lines) it is simpler: */
    tmp =  pow(X[0*nDonor+0]- X[0*nDonor+1],2.0);
    tmp += pow(X[1*nDonor+0]- X[1*nDonor+1],2.0);
    tmp = sqrt(tmp);

    tmp2 = pow(X[0*nDonor+0] - xj[0],2.0);
    tmp2 += pow(X[1*nDonor+0] - xj[1],2.0);
    tmp2 = sqrt(tmp2);
    isoparams[1] = tmp2/tmp;

    tmp2 = pow(X[0*nDonor+1] - xj[0],2.0);
    tmp2 += pow(X[1*nDonor+1] - xj[1],2.0);
    tmp2 = sqrt(tmp2);
    isoparams[0] = tmp2/tmp;
  }

  /*--- Isoparametric coefficients have been calculated. Run checks to eliminate outside-element issues ---*/
  if (nDonor==4){
    /*-- Bilinear coordinates, bounded by [-1,1] ---*/
    su2double xi, eta;
    xi = (1.0-isoparams[0]/isoparams[1])/(1.0+isoparams[0]/isoparams[1]);
    eta = 1- isoparams[2]*4/(1+xi);
    if (xi>1.0) xi=1.0;
    if (xi<-1.0) xi=-1.0;
    if (eta>1.0) eta=1.0;
    if (eta<-1.0) eta=-1.0;
    isoparams[0]=0.25*(1-xi)*(1-eta);
    isoparams[1]=0.25*(1+xi)*(1-eta);
    isoparams[2]=0.25*(1+xi)*(1+eta);
    isoparams[3]=0.25*(1-xi)*(1+eta);

  }
  if (nDonor<4){
    tmp = 0.0; // value for normalization
    tmp2=0; // check for maximum value, to be used to id nearest neighbor if necessary
    k=0; // index for maximum value
    for (iDonor=0; iDonor< nDonor; iDonor++){
      if (isoparams[iDonor]>tmp2){
        k=iDonor;
        tmp2=isoparams[iDonor];
      }
      // [0,1]
      if (isoparams[iDonor]<0) isoparams[iDonor]=0;
      if (isoparams[iDonor]>1) isoparams[iDonor]=1;
      tmp +=isoparams[iDonor];
    }
    if (tmp>0)
      for (iDonor=0; iDonor< nDonor; iDonor++)
        isoparams[iDonor]=isoparams[iDonor]/tmp;
    else{
      isoparams[k]=1.0;
    }
  }
  
  delete [] x;
  delete [] x_tmp;
  delete [] Q;
  delete [] R;
  delete [] A;
  if (A2 != NULL) delete [] A2;
  delete [] x2;
  
  delete [] test;
  delete [] testi;

}


/* Mirror Interpolator */
CMirror::CMirror(CGeometry ***geometry_container, CConfig **config,  unsigned int iZone, unsigned int jZone) :  CInterpolator(geometry_container, config, iZone, jZone){

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

}

CMirror::~CMirror(){}

void CMirror::Set_TransferCoeff(CConfig **config){
  unsigned long iVertex, jVertex;
  unsigned long iPoint;
  unsigned short iDonor=0, iFace=0, iTarget=0;

  unsigned short nMarkerInt, nMarkerDonor, nMarkerTarget;
  unsigned short iMarkerInt;

  int markDonor=0, markTarget=0;

  unsigned int nNodes=0, iNodes=0;
  unsigned long nVertexDonor = 0, nVertexTarget= 0;
  unsigned long Point_Donor = 0;
  unsigned long Global_Point = 0;
  unsigned long pGlobalPoint = 0;
  int iProcessor;

  unsigned long nLocalFace_Donor = 0, nLocalFaceNodes_Donor=0;

  unsigned long faceindex;

  int rank = MASTER_NODE;
  int nProcessor = SINGLE_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
  su2double *Buffer_Send_Coeff, *Buffer_Receive_Coeff;
  su2double coeff;

  /*--- Number of markers on the interface ---*/
  nMarkerInt = (config[targetZone]->GetMarker_n_FSIinterface())/2;
  nMarkerDonor  =  config[donorZone]->GetnMarker_All();
  nMarkerTarget =  config[targetZone]->GetnMarker_All();
  /*--- For the number of markers on the interface... ---*/
  for (iMarkerInt=1; iMarkerInt <= nMarkerInt; iMarkerInt++){
    /*--- Procedure:
     * -Loop through vertices of the aero grid
     * -Find nearest element and allocate enough space in the aero grid donor point info
     *    -set the transfer coefficient values
     */
    nVertexDonor = 0;
	nVertexTarget= 0;

	/*--- On the donor side ---*/
	markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);
	  
	/*--- On the target side ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

	if(markTarget == -1 || markDonor == -1)
		continue;
		    
    nVertexDonor  = donor_geometry->GetnVertex(markDonor);
    nVertexTarget = target_geometry->GetnVertex(markTarget);

    /*-- Collect the number of donor nodes: re-use 'Face' containers --*/
    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;
    for (jVertex = 0; jVertex<nVertexDonor; jVertex++) {
      Point_Donor =donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex

      if (donor_geometry->node[Point_Donor]->GetDomain()) {
        nNodes = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();
        nLocalFaceNodes_Donor+=nNodes;
        nLocalFace_Donor++;
      }
    }
    Buffer_Send_nFace_Donor= new unsigned long [1];
    Buffer_Send_nFaceNodes_Donor= new unsigned long [1];

    Buffer_Receive_nFace_Donor = new unsigned long [nProcessor];
    Buffer_Receive_nFaceNodes_Donor = new unsigned long [nProcessor];

    Buffer_Send_nFace_Donor[0] = nLocalFace_Donor;
    Buffer_Send_nFaceNodes_Donor[0] = nLocalFaceNodes_Donor;

    /*--- Send Interface vertex information --*/
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &MaxFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFace_Donor, &MaxFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFace_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MaxFace_Donor++;
#else
    nGlobalFace_Donor       = nLocalFace_Donor;
    nGlobalFaceNodes_Donor  = nLocalFaceNodes_Donor;
    MaxFaceNodes_Donor      = nLocalFaceNodes_Donor;
    MaxFace_Donor           = nLocalFace_Donor+1;
    Buffer_Receive_nFace_Donor[0] = Buffer_Send_nFace_Donor[0];
    Buffer_Receive_nFaceNodes_Donor[0] = Buffer_Send_nFaceNodes_Donor[0];
#endif

    /*-- Send donor info --*/
    Buffer_Send_FaceIndex   = new unsigned long[MaxFace_Donor];
    Buffer_Send_FaceNodes   = new unsigned long[MaxFaceNodes_Donor];
    //Buffer_Send_FaceProc    = new unsigned long[MaxFaceNodes_Donor];
    Buffer_Send_GlobalPoint = new unsigned long[MaxFaceNodes_Donor];
    Buffer_Send_Coeff       = new su2double[MaxFaceNodes_Donor];

    Buffer_Receive_FaceIndex= new unsigned long[MaxFace_Donor*nProcessor];
    Buffer_Receive_FaceNodes= new unsigned long[MaxFaceNodes_Donor*nProcessor];
    //Buffer_Receive_FaceProc = new unsigned long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_GlobalPoint = new unsigned long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_Coeff    = new su2double[MaxFaceNodes_Donor*nProcessor];

    for (iVertex=0; iVertex<MaxFace_Donor; iVertex++){
      Buffer_Send_FaceIndex[iVertex]=0;
    }
    for (iVertex=0; iVertex<MaxFaceNodes_Donor; iVertex++){
      Buffer_Send_FaceNodes[iVertex]=0;
      //Buffer_Send_FaceProc[iVertex]=0;
      Buffer_Send_GlobalPoint[iVertex]=0;
      Buffer_Send_Coeff[iVertex]=0.0;
    }
    for (iVertex=0; iVertex<MaxFace_Donor; iVertex++){
      Buffer_Send_FaceIndex[iVertex]=0;
    }

    Buffer_Send_FaceIndex[0]=rank*MaxFaceNodes_Donor;
    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;

    for (jVertex = 0; jVertex<nVertexDonor; jVertex++) {

      Point_Donor =donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex
      if (donor_geometry->node[Point_Donor]->GetDomain()) {
        nNodes = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();
        for (iDonor=0; iDonor<nNodes; iDonor++){
          Buffer_Send_FaceNodes[nLocalFaceNodes_Donor] = donor_geometry->node[Point_Donor]->GetGlobalIndex();
          Buffer_Send_GlobalPoint[nLocalFaceNodes_Donor] =
              donor_geometry->vertex[markDonor][jVertex]->GetInterpDonorPoint(iDonor);
          Buffer_Send_Coeff[nLocalFaceNodes_Donor] =
              donor_geometry->vertex[markDonor][jVertex]->GetDonorCoeff(iDonor);
          nLocalFaceNodes_Donor++;
        }
        Buffer_Send_FaceIndex[nLocalFace_Donor+1] =Buffer_Send_FaceIndex[nLocalFace_Donor]+nNodes;
        nLocalFace_Donor++;
      }
    }

#ifdef HAVE_MPI
    SU2_MPI::Allgather(Buffer_Send_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_GlobalPoint, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG,Buffer_Receive_GlobalPoint, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Coeff, MaxFaceNodes_Donor, MPI_DOUBLE,Buffer_Receive_Coeff, MaxFaceNodes_Donor, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
    for (iFace=0; iFace<MaxFace_Donor; iFace++){
      Buffer_Receive_FaceIndex[iFace] = Buffer_Send_FaceIndex[iFace];
    }
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertex++){
      Buffer_Receive_FaceNodes[iVertex] = Buffer_Send_FaceNodes[iVertex];
      Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];
      Buffer_Receive_Coeff[iVertex] = Buffer_Send_Coeff[iVertex];
    }
#endif
    /*--- Loop over the vertices on the target Marker ---*/
    for (iVertex = 0; iVertex<nVertexTarget; iVertex++) {

      iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();
      if (target_geometry->node[iPoint]->GetDomain()) {
        Global_Point = target_geometry->node[iPoint]->GetGlobalIndex();
        nNodes = 0;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){
          for (iFace = 0; iFace < Buffer_Receive_nFace_Donor[iProcessor]; iFace++) {
            faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
            iNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1]- (unsigned int)faceindex;
            for (iTarget=0; iTarget<iNodes; iTarget++){
              if (Global_Point == Buffer_Receive_GlobalPoint[faceindex+iTarget])
                nNodes++;
              //coeff =Buffer_Receive_Coeff[faceindex+iDonor];
            }
          }
        }

        target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nNodes);
        target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();

        iDonor = 0;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){
          for (iFace = 0; iFace < Buffer_Receive_nFace_Donor[iProcessor]; iFace++) {

            faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
            iNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1]- (unsigned int)faceindex;
            for (iTarget=0; iTarget<iNodes; iTarget++){
              if (Global_Point == Buffer_Receive_GlobalPoint[faceindex+iTarget]){
                coeff =Buffer_Receive_Coeff[faceindex+iTarget];
                pGlobalPoint = Buffer_Receive_FaceNodes[faceindex+iTarget];
                target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint(iDonor,pGlobalPoint);
                target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(iDonor,coeff);
                target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(iDonor, iProcessor);
                //cout <<rank << " Global Point " << Global_Point<<" iDonor " << iDonor <<" coeff " << coeff <<" gp " << pGlobalPoint << endl;
                iDonor++;
              }
            }
          }
        }
      }
    }
    delete[] Buffer_Send_nFace_Donor;
    delete[] Buffer_Send_nFaceNodes_Donor;

    delete[] Buffer_Receive_nFace_Donor;
    delete[] Buffer_Receive_nFaceNodes_Donor;

    delete[] Buffer_Send_FaceIndex;
    delete[] Buffer_Send_FaceNodes;
    delete[] Buffer_Send_GlobalPoint;
    delete[] Buffer_Send_Coeff;

    delete[] Buffer_Receive_FaceIndex;
    delete[] Buffer_Receive_FaceNodes;
    delete[] Buffer_Receive_GlobalPoint;
    delete[] Buffer_Receive_Coeff;

  }
}

CSlidingmesh::CSlidingmesh(CGeometry ***geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone)  :  CInterpolator(geometry_container, config, iZone, jZone){

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

  /*--- For fluid-structure interaction data interpolated with have nDim dimensions ---*/
 // InitializeData(Zones,nDim);
}

CSlidingmesh::~CSlidingmesh(){}

void CSlidingmesh::Set_TransferCoeff(CConfig **config){


	unsigned long jVertex, iVertex;

	unsigned short iDim, nDim = donor_geometry->GetnDim();

	unsigned short iMarkerInt, nMarkerInt, nMarkerDonor, nMarkerTarget;	

	int markDonor, markTarget, NearestNode;

	unsigned long nVertexDonor , nVertexTarget;



	unsigned long iVertexTarget, iVertexDonor;


	su2double *Coord_i, *Coord_j, dist, mindist, maxdist, *Normal;

	int rank = MASTER_NODE;
	int nProcessor = SINGLE_NODE;
	
	bool check;
	int target_StartIndex, donor_StartIndex;
	int target_forward_point, target_backward_point, donor_forward_point, donor_backward_point;
	int donor_iPoint, target_iPoint;
	int donor_OldiPoint, target_OldiPoint;
	
	int ii, jj;

	unsigned long nDonorPoints, iDonor;
	unsigned long *Donor_Vect, *tmp_Donor_Vect;
		
	su2double Intersection, length, dTMP;
	su2double *target_iMidEdge_point, *target_jMidEdge_point, *donor_iMidEdge_point, *donor_jMidEdge_point, *Direction;
	su2double *Coeff_Vect, *tmp_Coeff_Vect;
	
	int nAlreadyVisited, nToVisit, StartVisited, vPoint;
					
	int *alreadyVisitedDonor, *ToVisit, *tmpVect;
					
	int iEdgeVisited, EgdeIndex, nEdgeVisited, iNodeVisited;


	#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
	#endif


	/*--- Number of markers on the FSI interface ---*/
	
	nMarkerInt     = (int) (config[donorZone]->GetMarker_n_FSIinterface() /2);
	nMarkerTarget  = target_geometry->GetnMarker();
	nMarkerDonor   = donor_geometry->GetnMarker();

	/*--- Setting up auxiliary vectors ---*/
	
	Donor_Vect = NULL;
	Coeff_Vect = NULL;
	
	tmp_Donor_Vect = NULL;
	tmp_Coeff_Vect = NULL;
	
	Normal = new su2double[nDim];
	
	Coord_i = new su2double[nDim];
	Coord_j = new su2double[nDim];
	
	target_iMidEdge_point = new su2double[nDim];
	target_jMidEdge_point = new su2double[nDim];
	
	donor_iMidEdge_point = new su2double[nDim];
	donor_jMidEdge_point = new su2double[nDim];
	
	Direction = new su2double[nDim];
	
/* 1- Find boundary tag between touching grids */

	/*--- Number of markers on the FSI interface ---*/
	nMarkerInt    = (int)(config[donorZone]->GetMarker_n_FSIinterface())/2;
	nMarkerTarget = target_geometry->GetnMarker();
	nMarkerDonor  = donor_geometry->GetnMarker();


	/*--- For the number of markers on the interface... ---*/
	for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++){
		
		/*--- Procedure:
		* -Loop through vertices of the aero grid
		* -Find nearest element and allocate enough space in the aero grid donor point info
		*    -set the transfer coefficient values
		*/
		
		nVertexDonor = 0;
		nVertexTarget= 0;

		/*--- On the donor side ---*/
		markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);

		/*--- On the target side ---*/
		markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

		if(markTarget == -1 || markDonor == -1)
			continue;

		nVertexDonor  = donor_geometry->GetnVertex(markDonor);
		nVertexTarget = target_geometry->GetnVertex(markTarget);


/*
2- For each element find all the elements from the outer zone that are touching its boundary face
* - Find outer nearest neighbour, from it retrieve the boundary face of elements sharing that vertex and build a local supermesh until the area of the initial boundary face is covered
*/
		
		/*--- Setting up parallel data structures ---*/

		Buffer_Send_nVertex_Donor    = new unsigned long [1];
		Buffer_Receive_nVertex_Donor = new unsigned long [nProcessor];

		/*--- Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor ---*/
		Determine_ArraySize(false, markDonor, markTarget, nVertexDonor,nDim);

		Buffer_Send_Coord       = new su2double [MaxLocalVertex_Donor*nDim];
		Buffer_Send_GlobalPoint = new unsigned long [MaxLocalVertex_Donor];

		Buffer_Receive_Coord       = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
		Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Donor];

		/*--- Collect coordinates, global points, and normal vectors ---*/
		Collect_VertexInfo(false, markDonor, markTarget, nVertexDonor,nDim);


		if(nDim == 2){
			/*--- Starts with supermesh reconstruction ---*/

			target_StartIndex = target_geometry->vertex[markTarget][0]->GetNode();
					
			for (iVertex = 0; iVertex < nVertexTarget; iVertex++) {

				nDonorPoints = 0;
				
				/*--- Stores coordinates of the target node ---*/
				
				target_iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();
				
				for (iDim = 0; iDim < nDim; iDim++)
					Coord_i[iDim] = target_geometry->node[target_iPoint]->GetCoord(iDim);
				
				/*--- Brute force to find the closest donor_node ---*/
				
				mindist = 1E6;
				
				for (jVertex = 0; jVertex < nVertexDonor; jVertex++) {
					
					donor_iPoint = donor_geometry->vertex[markDonor][jVertex]->GetNode();

					for (iDim = 0; iDim < nDim; iDim++) 
						Coord_j[iDim] = donor_geometry->node[donor_iPoint]->GetCoord(iDim);
		
					dist = PointsDistance(Coord_i, Coord_j);

					if (dist < mindist) {
						mindist = dist;  
						donor_StartIndex = donor_iPoint;
					}

					if (dist == 0.0){
						donor_StartIndex = donor_iPoint;
						break;
					}
				}
				
				/*--- Contruct information regarding the target cell ---*/
				
				donor_iPoint    = donor_StartIndex;
				donor_OldiPoint = -1;
				
				target_forward_point  = FindNextNode_2D(target_geometry,                    -1, target_iPoint, markTarget);
				target_backward_point = FindNextNode_2D(target_geometry,  target_forward_point, target_iPoint, markTarget);
				
				dTMP = 0;
				for(iDim = 0; iDim < nDim; iDim++){
					target_iMidEdge_point[iDim] = ( target_geometry->node[target_forward_point ]->GetCoord(iDim) + target_geometry->node[ target_iPoint ]->GetCoord(iDim) ) / 2;
					target_jMidEdge_point[iDim] = ( target_geometry->node[target_backward_point]->GetCoord(iDim) + target_geometry->node[ target_iPoint ]->GetCoord(iDim) ) / 2;
					
					Direction[iDim] = target_jMidEdge_point[iDim] - target_iMidEdge_point[iDim];
					dTMP += Direction[iDim] * Direction[iDim];
				}
				
				dTMP = sqrt(dTMP);
				for(iDim = 0; iDim < nDim; iDim++)
					Direction[iDim] /= dTMP;
					
				length = PointsDistance(target_iMidEdge_point, target_jMidEdge_point);

				check = false;


				/*--- Proceeds along the forward direction (depending on which connected boundary node is found first) ---*/
				
				while( !check ){
					
					/*--- Proceeds until the value of the intersection area is null ---*/
					
					donor_forward_point  = FindNextNode_2D(donor_geometry,      donor_OldiPoint, donor_iPoint, markDonor);
					donor_backward_point = FindNextNode_2D(donor_geometry,  donor_forward_point, donor_iPoint, markDonor);
					
					for(iDim = 0; iDim < nDim; iDim++){
						donor_iMidEdge_point[iDim] = ( donor_geometry->node[donor_forward_point ]->GetCoord(iDim) + donor_geometry->node[ donor_iPoint ]->GetCoord(iDim) ) / 2;
						donor_jMidEdge_point[iDim] = ( donor_geometry->node[donor_backward_point]->GetCoord(iDim) + donor_geometry->node[ donor_iPoint ]->GetCoord(iDim) ) / 2;
					}
								
					Intersection = Compute_Intersection_2D(target_iMidEdge_point, target_jMidEdge_point, donor_iMidEdge_point, donor_jMidEdge_point, Direction);
			
					if ( Intersection == 0.0 ){
						check = true;
						continue;
					}
					
					/*--- In case the element intersect the target cell update the auxiliary communication data structure ---*/
						
					tmp_Coeff_Vect = new     su2double[ nDonorPoints ];
					tmp_Donor_Vect = new unsigned long[ nDonorPoints ];
					
					for( iDonor = 0; iDonor < nDonorPoints; iDonor++){
						tmp_Donor_Vect[iDonor] = Donor_Vect[iDonor];
						tmp_Coeff_Vect[iDonor] = Coeff_Vect[iDonor];
					}
					
					if (Donor_Vect != NULL)
						delete [] Donor_Vect;
						
					if (Coeff_Vect != NULL)
						delete [] Coeff_Vect;
					
					Coeff_Vect = new     su2double[ nDonorPoints + 1 ];
					Donor_Vect = new unsigned long[ nDonorPoints + 1 ];
					
					for( iDonor = 0; iDonor < nDonorPoints; iDonor++){
						Donor_Vect[iDonor] = tmp_Donor_Vect[iDonor];
						Coeff_Vect[iDonor] = tmp_Coeff_Vect[iDonor];
					}

					Coeff_Vect[ nDonorPoints ] = Intersection / length;					
					Donor_Vect[ nDonorPoints ] = donor_iPoint;

					if (tmp_Donor_Vect != NULL)
						delete [] tmp_Donor_Vect;
			
					if (tmp_Coeff_Vect != NULL)
						delete [] tmp_Coeff_Vect;

					
					donor_OldiPoint = donor_iPoint;
					donor_iPoint = donor_forward_point;
						
					nDonorPoints++;
				}
			
			
				donor_iPoint    = donor_StartIndex;
				donor_OldiPoint = FindNextNode_2D(donor_geometry, -1, donor_iPoint, markDonor);
				
				donor_iPoint = FindNextNode_2D(donor_geometry, donor_OldiPoint, donor_iPoint, markDonor);
				donor_OldiPoint = donor_StartIndex;
				
				check = false;
				
				/*--- Proceeds along the backward direction (depending on which connected boundary node is found first) ---*/
				
				while( !check ){
					
						/*--- Proceeds until the value of the intersection area is null ---*/
						
						donor_forward_point  = FindNextNode_2D(donor_geometry,     donor_OldiPoint, donor_iPoint, markDonor);
						donor_backward_point = FindNextNode_2D(donor_geometry, donor_forward_point, donor_iPoint, markDonor);
						
						for(iDim = 0; iDim < nDim; iDim++){
							donor_iMidEdge_point[iDim] = ( donor_geometry->node[donor_forward_point ]->GetCoord(iDim) + donor_geometry->node[ donor_iPoint ]->GetCoord(iDim) ) / 2;
							donor_jMidEdge_point[iDim] = ( donor_geometry->node[donor_backward_point]->GetCoord(iDim) + donor_geometry->node[ donor_iPoint ]->GetCoord(iDim) ) / 2;
						}		
						
						Intersection = Compute_Intersection_2D(target_iMidEdge_point, target_jMidEdge_point, donor_iMidEdge_point, donor_jMidEdge_point, Direction);

						if ( Intersection == 0.0 ){
							check = true;
							continue;
						}
						
						/*--- In case the element intersect the target cell update the auxiliary communication data structure ---*/
						
						tmp_Coeff_Vect = new     su2double[ nDonorPoints ];
						tmp_Donor_Vect = new unsigned long[ nDonorPoints ];
						
						for( iDonor = 0; iDonor < nDonorPoints; iDonor++){
							tmp_Donor_Vect[iDonor] = Donor_Vect[iDonor];
							tmp_Coeff_Vect[iDonor] = Coeff_Vect[iDonor];
						}
						
						if (Donor_Vect != NULL)
							delete [] Donor_Vect;
							
						if (Coeff_Vect != NULL)
							delete [] Coeff_Vect;
						
						Coeff_Vect = new     su2double[ nDonorPoints + 1 ];
						Donor_Vect = new unsigned long[ nDonorPoints + 1 ];
						
						for( iDonor = 0; iDonor < nDonorPoints; iDonor++){
							Donor_Vect[iDonor] = tmp_Donor_Vect[iDonor];
							Coeff_Vect[iDonor] = tmp_Coeff_Vect[iDonor];
						}
							
						Donor_Vect[ nDonorPoints ] = donor_iPoint;
						Coeff_Vect[ nDonorPoints ] = Intersection / length;
							
						if (tmp_Donor_Vect != NULL)
							delete [] tmp_Donor_Vect;
						
						if (tmp_Coeff_Vect != NULL)
							delete [] tmp_Coeff_Vect;
							
							
						donor_OldiPoint = donor_iPoint;
						donor_iPoint = donor_forward_point;
							
						nDonorPoints++;
				}
				
				/*--- Set the communication data structure and copy data from the auxiliary vectors ---*/
				
				target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nDonorPoints);

				target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();

				for ( iDonor = 0; iDonor < nDonorPoints; iDonor++ ){			  
					target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(iDonor, Coeff_Vect[iDonor]);
					target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint( iDonor, donor_geometry->node[ Donor_Vect[iDonor] ]->GetGlobalIndex() );
					target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(iDonor, 0);//storeProc[iDonor]);

					//cout <<rank << " Global Point " << Global_Point<<" iDonor " << iDonor <<" coeff " << coeff <<" gp " << pGlobalPoint << endl;
				}
			
			}

		}
		else{
			
			
			su2double **target_element, **donor_element;
			int nEdges_target, nNode_target;
			int nEdges_donor,  nNode_donor;
			
			
			
			su2double Area, Area_old, tmp_Area, target_area;
			
			
					
			for (iVertex = 0; iVertex < nVertexTarget; iVertex++) {
				//cout << "iVertex  " << iVertex << endl;
				nDonorPoints = 0;
				
				/*--- Stores coordinates of the target node ---*/
				
				target_iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();
				
				target_geometry->vertex[markTarget][iVertex]->GetNormal(Normal);
				
				Area = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) 
					Area += Normal[iDim]*Normal[iDim];
				Area = sqrt(Area);
				
				target_area = Area;
				
				//cout << "Area target " << Area << endl;

				for (iDim = 0; iDim < nDim; iDim++)
					Normal[iDim] /= Area;

				for (iDim = 0; iDim < nDim; iDim++)
					Coord_i[iDim] = target_geometry->node[target_iPoint]->GetCoord(iDim);
					
				/*--- Build local surface dual mesh for target element---*/
				
				nEdges_target = target_geometry->node[target_iPoint]->GetnPoint();
				
				nNode_target = nEdges_target*2 + 1;

				target_element = new su2double*[nNode_target];
				for (ii = 0; ii < nNode_target; ii++)
					target_element[ii] = new su2double[nDim];
					
					
				nNode_target = Build_3D_surface_element(target_geometry, target_iPoint, markTarget, target_element);

				/*--- Brute force to find the closest donor_node ---*/
			
				mindist = 1E6;
				
				for (jVertex = 0; jVertex < nVertexDonor; jVertex++) {
					
					donor_iPoint = donor_geometry->vertex[markDonor][jVertex]->GetNode();

					for (iDim = 0; iDim < nDim; iDim++) 
						Coord_j[iDim] = donor_geometry->node[donor_iPoint]->GetCoord(iDim);
		
					dist = PointsDistance(Coord_i, Coord_j);

					if (dist < mindist) {
						mindist = dist;  
						donor_StartIndex = donor_iPoint;
					}

					if (dist == 0.0){
						donor_StartIndex = donor_iPoint;
						break;
					}
				}
				
				donor_iPoint = donor_StartIndex;
				
				nEdges_donor = donor_geometry->node[donor_iPoint]->GetnPoint();
				
				donor_element = new su2double*[ nEdges_donor*2 + 1 ];
				for (ii = 0; ii < nEdges_donor*2 + 1; ii++)
					donor_element[ii] = new su2double[nDim];				
					
				nNode_donor = Build_3D_surface_element(donor_geometry, donor_iPoint, markDonor, donor_element);

				Area = 0;
				for (ii = 1; ii < nNode_target; ii++){
					for (jj = 1; jj < nNode_donor; jj++){
						Area += Compute_Triangle_Intersection(target_element[0], target_element[ii], target_element[ii+1], donor_element[0], donor_element[jj], donor_element[jj+1], Normal);
						//cout << Compute_Triangle_Intersection(target_element[0], target_element[ii], target_element[ii+1], donor_element[0], donor_element[jj], donor_element[jj+1], Normal) << endl;
					}
				}
				
				for (ii = 0; ii < nEdges_donor*2 + 1; ii++)
							delete [] donor_element[ii];
				delete [] donor_element;
				
				nDonorPoints = 1;
				
				/*--- In case the element intersect the target cell update the auxiliary communication data structure ---*/
						
				Coeff_Vect = new     su2double[ nDonorPoints ];
				Donor_Vect = new unsigned long[ nDonorPoints ];
				
				Coeff_Vect[0] = Area;
				Donor_Vect[0] = donor_iPoint;
				
				//cout << "Area zero " << scientific << Area << "  " << Area - target_area << endl;
				
				nAlreadyVisited = 0;
				
				if( Area/target_area < 1 ){
					
					alreadyVisitedDonor = new int[1];
					
					alreadyVisitedDonor[0] = donor_iPoint;
					nAlreadyVisited = 1;
					StartVisited = 0;
					
					Area_old = -1;
					
					while( Area != Area_old ){
						/*
							for( jj = 0; jj < nAlreadyVisited; jj++ )
								cout << alreadyVisitedDonor[jj] << endl;
						*/	
						
							Area_old = Area;
							
							ToVisit = NULL;
							nToVisit = 0;
							
							for( iNodeVisited = StartVisited; iNodeVisited < nAlreadyVisited; iNodeVisited++ ){
							//for( iNodeVisited = 0; iNodeVisited < nAlreadyVisited; iNodeVisited++ ){
							
								vPoint = alreadyVisitedDonor[ iNodeVisited ];
								
								//cout << "init " << vPoint << endl;
								
								nEdgeVisited = donor_geometry->node[ vPoint ]->GetnPoint();
							
								for (iEdgeVisited = 0; iEdgeVisited < nEdgeVisited; iEdgeVisited++){

									EgdeIndex = donor_geometry->node[ vPoint ]->GetEdge( iEdgeVisited );
										
									if( vPoint == donor_geometry->edge[ EgdeIndex ]->GetNode(0) )
										donor_iPoint = donor_geometry->edge[ EgdeIndex ]->GetNode(1);
									else
										donor_iPoint = donor_geometry->edge[ EgdeIndex ]->GetNode(0);

									if ( donor_geometry->node[donor_iPoint]->GetVertex(markDonor) != -1 ){
											
										check = 0;
										
										for( jj = 0; jj < nAlreadyVisited; jj++ ){
											if( donor_iPoint == alreadyVisitedDonor[jj] ){
													check = 1; break;
											}
										}
										
										if( check == 0 ){
											for( jj = 0; jj < nToVisit; jj++ )
												if( donor_iPoint == ToVisit[jj] ){
														check = 1; break;
												}		
										}
										
										if( check == 0 ){ // Add node to visit and visit it
											
											//cout << "check " << vPoint << endl;
											
											tmpVect = new int[ nToVisit ];
											
											for( jj = 0; jj < nToVisit; jj++ )
												tmpVect[jj] = ToVisit[jj];
												
											if( ToVisit != NULL )
												delete [] ToVisit;
												
											ToVisit = new int[nToVisit + 1];
											
											for( jj = 0; jj < nToVisit; jj++ )
												ToVisit[jj] = tmpVect[jj];
															
											delete [] tmpVect;				
											
											ToVisit[nToVisit] = donor_iPoint;
											
											nToVisit++;	
											
											// Find partial coefficient

											nEdges_donor = donor_geometry->node[donor_iPoint]->GetnPoint();
											
											donor_element = new su2double*[ nEdges_donor*2 + 1 ];
											for (ii = 0; ii < nEdges_donor*2 + 1; ii++)
												donor_element[ii] = new su2double[nDim];				
												
											nNode_donor = Build_3D_surface_element(donor_geometry, donor_iPoint, markDonor, donor_element);

											tmp_Area = 0;
											for (ii = 1; ii < nNode_target; ii++){
												for (jj = 1; jj < nNode_donor; jj++){
													tmp_Area += Compute_Triangle_Intersection(target_element[0], target_element[ii], target_element[ii+1], donor_element[0], donor_element[jj], donor_element[jj+1], Normal);
													//cout << Compute_Triangle_Intersection(target_element[0], target_element[ii], target_element[ii+1], donor_element[0], donor_element[jj], donor_element[jj+1], Normal) << endl;
												}
											}
											
											for (ii = 0; ii < nEdges_donor*2 + 1; ii++)
														delete [] donor_element[ii];
											delete [] donor_element;
											
											
											/*--- In case the element intersect the target cell update the auxiliary communication data structure ---*/
						
											tmp_Coeff_Vect = new     su2double[ nDonorPoints ];
											tmp_Donor_Vect = new unsigned long[ nDonorPoints ];
											
											for( iDonor = 0; iDonor < nDonorPoints; iDonor++){
												tmp_Donor_Vect[iDonor] = Donor_Vect[iDonor];
												tmp_Coeff_Vect[iDonor] = Coeff_Vect[iDonor];
											}
											
											if (Donor_Vect != NULL)
												delete [] Donor_Vect;
												
											if (Coeff_Vect != NULL)
												delete [] Coeff_Vect;
											
											Coeff_Vect = new     su2double[ nDonorPoints + 1 ];
											Donor_Vect = new unsigned long[ nDonorPoints + 1 ];
											
											for( iDonor = 0; iDonor < nDonorPoints; iDonor++){
												Donor_Vect[iDonor] = tmp_Donor_Vect[iDonor];
												Coeff_Vect[iDonor] = tmp_Coeff_Vect[iDonor];
											}

											Coeff_Vect[ nDonorPoints ] = tmp_Area;					
											Donor_Vect[ nDonorPoints ] = donor_iPoint;

											if (tmp_Donor_Vect != NULL)
												delete [] tmp_Donor_Vect;
									
											if (tmp_Coeff_Vect != NULL)
												delete [] tmp_Coeff_Vect;
												
											nDonorPoints++;
											
											Area += tmp_Area;
						
										}
										
									}
										
								}	
							}
							
							StartVisited = nAlreadyVisited;
							
							tmpVect = new int[ nAlreadyVisited ];
											
							for( jj = 0; jj < nAlreadyVisited; jj++ )
								tmpVect[jj] = alreadyVisitedDonor[jj];
								
							if( alreadyVisitedDonor != NULL )
								delete [] alreadyVisitedDonor;
								
							alreadyVisitedDonor = new int[ nAlreadyVisited + nToVisit ];
							
							for( jj = 0; jj < nAlreadyVisited; jj++ )
								alreadyVisitedDonor[jj] = tmpVect[jj];
								
							for( jj = 0; jj < nToVisit; jj++ )
								alreadyVisitedDonor[ nAlreadyVisited + jj ] = ToVisit[jj];
											
							delete [] tmpVect;				
							
							nAlreadyVisited += nToVisit;							
							/*
							if( Area/target_area > 1 ){
								cout << "Area " << scientific << Area << "  " << Area - target_area << "  " << Area/target_area << "  " << (Area - target_area)/target_area  << endl;
								getchar();
							}
							*/
							/*
							cout << "Visited  " << nAlreadyVisited << endl;
							cout << "Area " << scientific << Area << "  " << Area - target_area << endl;
							getchar();
							*/
					}
				
				
					delete [] alreadyVisitedDonor;
					delete [] ToVisit;
					
				}
				/*
				cout << "Visited  " << nAlreadyVisited << endl;
				cout << "Area " << scientific << Area << "  " << Area - target_area << endl;
				cout << endl;
				*/
			/*
				if( Area/target_area < 1 ){
					cout << "Area " << scientific << Area << "  " << Area - target_area << endl;
					getchar();
				}
			*/
				
				
				/*--- Set the communication data structure and copy data from the auxiliary vectors ---*/
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nDonorPoints);
	
				target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();

				for ( iDonor = 0; iDonor < nDonorPoints; iDonor++ ){			  
					target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(iDonor, Coeff_Vect[iDonor] / Area);
					target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint( iDonor, donor_geometry->node[ Donor_Vect[iDonor] ]->GetGlobalIndex() );
					target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(iDonor, 0);//storeProc[iDonor]);

					//cout <<rank << " Global Point " << Global_Point<<" iDonor " << iDonor <<" coeff " << coeff <<" gp " << pGlobalPoint << endl;				 
				}

				
				
				
				
				
				
				
				for (ii = 0; ii < nEdges_target*2 + 1; ii++)
							delete [] target_element[ii];
				delete [] target_element;
			}
		}
		
		
		
		delete [] Buffer_Send_Coord;
		delete [] Buffer_Send_GlobalPoint;

		delete [] Buffer_Receive_Coord;
		delete [] Buffer_Receive_GlobalPoint;

		delete [] Buffer_Send_nVertex_Donor;
		delete [] Buffer_Receive_nVertex_Donor;
	}


	//Compute_Triangle_Intersection(Coord_i, Coord_i, Coord_i, Coord_i, Coord_i, Coord_i, Coord_i);

	delete [] Normal;
	
	delete [] Coord_i;
	delete [] Coord_j;
	
	delete [] target_iMidEdge_point;
	delete [] target_jMidEdge_point;
	
	delete [] donor_iMidEdge_point;
	delete [] donor_jMidEdge_point;
	
	delete [] Direction;
}

su2double CSlidingmesh::PointsDistance(su2double *point_i, su2double *point_j){

	/*--- Compute distance between 2 points ---*/
	
	unsigned short iDim, nDim = donor_geometry->GetnDim();
	su2double m;
	
	m = 0 ;
	for(iDim = 0; iDim < nDim; iDim++)
		m += (point_j[iDim] - point_i[iDim])*(point_j[iDim] - point_i[iDim]);
	
	return sqrt(m);
}

int CSlidingmesh::Build_3D_surface_element(CGeometry *geometry, unsigned long centralNode, unsigned long markID, su2double** element){
	
	int nEdges, iEdge, jEdge, kEdge, nNode, iNode, iEdgeIndex, jEdgeIndex, kEdgeIndex, iPoint, jPoint, kPoint, i, OldiPoint, StartIndex;
	
	unsigned short nDim = 3, iDim;
	
	nEdges = geometry->node[centralNode]->GetnPoint();
	
	for (iDim = 0; iDim < nDim; iDim++)
		element[0][iDim] = geometry->node[centralNode]->GetCoord(iDim);
				
	iNode = 1;
				
	for (jEdge = 0; jEdge < nEdges; jEdge++){

		iEdgeIndex = geometry->node[centralNode]->GetEdge(jEdge);
			
		if( centralNode == geometry->edge[iEdgeIndex]->GetNode(0) )
			iPoint = geometry->edge[iEdgeIndex]->GetNode(1);
		else
			iPoint = geometry->edge[iEdgeIndex]->GetNode(0);

		if ( geometry->node[iPoint]->GetVertex(markID) != -1 )
			break;
	}
				
	for (iDim = 0; iDim < nDim; iDim++) 
			element[iNode][iDim] = ( element[0][iDim] + geometry->node[iPoint]->GetCoord(iDim) )/2;
			
	StartIndex = iPoint;
	
	OldiPoint = -1;
				
	while( iNode != nNode-1 ){	
						
		for (jEdge = 0; jEdge < nEdges; jEdge++){
			
			jEdgeIndex = geometry->node[centralNode]->GetEdge(jEdge);
			
			
			if( jEdgeIndex != iEdgeIndex ){
			
				if( centralNode == geometry->edge[jEdgeIndex]->GetNode(0))
					jPoint = geometry->edge[jEdgeIndex]->GetNode(1);
				else
					jPoint = geometry->edge[jEdgeIndex]->GetNode(0);
					
		
				if ( geometry->node[jPoint]->GetVertex(markID) != -1 && jPoint != OldiPoint ){
					
					for ( kEdge = 0; kEdge < geometry->node[iPoint]->GetnPoint(); kEdge++ ){
										
						kEdgeIndex = geometry->node[iPoint]->GetEdge(kEdge);
						
						if( iPoint == geometry->edge[kEdgeIndex]->GetNode(0) )
							kPoint = geometry->edge[kEdgeIndex]->GetNode(1);
						else
							kPoint = geometry->edge[kEdgeIndex]->GetNode(0);
						
						if( kPoint == jPoint ){
							
							iNode++;
							
							for (iDim = 0; iDim < nDim; iDim++) 
								element[iNode][iDim] = ( element[0][iDim] + geometry->node[iPoint]->GetCoord(iDim) + geometry->node[jPoint]->GetCoord(iDim))/3;
							
							iNode++;

							for (iDim = 0; iDim < nDim; iDim++) 
								element[iNode][iDim] = ( element[0][iDim] + geometry->node[jPoint]->GetCoord(iDim) )/2;

							OldiPoint = iPoint;	
							iPoint = jPoint;
							iEdgeIndex = jEdgeIndex;
						
							break;
						}
					}
				}				
				
			}
			
			if(iPoint == jPoint)
				break;
		}
		
		if(iPoint == StartIndex)
			break;
	}
	
	for (iDim = 0; iDim < nDim; iDim++)
		element[iNode+1][iDim] = element[1][iDim];
		
				
	return iNode;
}

su2double CSlidingmesh::Compute_Intersection_2D(su2double* A1, su2double* A2, su2double* B1, su2double* B2, su2double* Direction){
	
	/*--- Given 2 segments, each defined by 2 points, it projects them along a given direction and it computes the length of the segment resulting from their intersection ---*/
	/*--- The algorithm works both for 2D and 3D problems ---*/
	
	unsigned short iDim;
	unsigned short nDim = donor_geometry->GetnDim();
	
	su2double Intersection;
	
	su2double dotA2, dotB1, dotB2, MaxArea;
	
	dotA2 = 0;
	for(iDim = 0; iDim < nDim; iDim++)
		dotA2 += ( A2[iDim] - A1[iDim] ) * Direction[iDim];

	if( dotA2 >= 0 ){
		dotB1 = 0;
		dotB2 = 0;
		for(iDim = 0; iDim < nDim; iDim++){
			dotB1 += ( B1[iDim] - A1[iDim] ) * Direction[iDim];
			dotB2 += ( B2[iDim] - A1[iDim] ) * Direction[iDim];
		}
	}
	else{
		dotA2 *= -1;
		
		dotB1 = 0;
		dotB2 = 0;
		for(iDim = 0; iDim < nDim; iDim++){
			dotB1 -= ( B1[iDim] - A1[iDim] ) * Direction[iDim];
			dotB2 -= ( B2[iDim] - A1[iDim] ) * Direction[iDim];
		}
	}
	
	
	if( dotB1 >= 0 && dotB1 <= dotA2 ){
		if ( dotB2 < 0 )
			return fabs( dotB1 );
		if ( dotB2 > dotA2 )
			return fabs( dotA2 - dotB1 );
			
		return fabs( dotB1 - dotB2 );
	}
	
	if( dotB2 >= 0 && dotB2 <= dotA2 ){
		if ( dotB1 < 0 )
			return fabs(dotB2);
		if ( dotB1 > dotA2 )
			return fabs( dotA2 - dotB2 );
	}
	
	if( ( dotB1 <= 0 && dotA2 <= dotB2 ) || ( dotB2 <= 0 && dotA2 <= dotB1 ) )
		return fabs( dotA2 );
		
	return 0.0;
}

int CSlidingmesh::FindNextNode_2D(CGeometry *geometry, int PreviousNode, unsigned long NodeID, unsigned long markID){
		
		/*--- ONLY for 2D grids ---*/
		/*--- It takes as input the grid, then starts from the NodeID and searches for its neigbours on the specified markID ---*/
		/*--- Since in 2D there are 2 possible neighbours, PreviousNode is needed to move along a certain direction along the boundary ---*/

		int iPoint, jPoint, iEdge, nEdges, EdgeIndex;
		bool check;

		check = false;
			
		nEdges = geometry->node[NodeID]->GetnPoint();

		for (iEdge = 0; iEdge < nEdges; iEdge++){
			
			EdgeIndex = geometry->node[NodeID]->GetEdge(iEdge);
			
			if( NodeID == geometry->edge[EdgeIndex]->GetNode(0) ){
				jPoint =  geometry->edge[EdgeIndex]->GetNode(1);
			}
			else
				jPoint = geometry->edge[EdgeIndex]->GetNode(0);

			if ( geometry->node[jPoint]->GetVertex(markID) != -1 && jPoint != PreviousNode )
				return jPoint;
		}

		return -1;
}

su2double CSlidingmesh::Compute_Triangle_Intersection(su2double* A1, su2double* A2, su2double* A3, su2double* B1, su2double* B2, su2double* B3, su2double* Direction){
	
	/* --- This routine is ONLY for 3D grids --- */
	
	unsigned short iDim;
	unsigned short nDim = donor_geometry->GetnDim();
	
	su2double I[3], J[3], K[3];
	su2double a1[3], a2[3], a3[3];
	su2double b1[3], b2[3], b3[3];
	su2double m1, m2;
	su2double dot, r[3], s[3];
	
	/* --- Reference frame is determined by: x = A1A2 y = x ^ ( -Direction ) --- */
	
	for(iDim = 0; iDim < 3; iDim++){
		I[iDim] = 0;
		K[iDim] = 0;
		
		s[iDim] = 0;
		r[iDim] = 0;
		
		a1[iDim] = 0;
		a2[iDim] = 0;
		a3[iDim] = 0;
		
		b1[iDim] = 0;
		b2[iDim] = 0;
		b3[iDim] = 0;
	}
	
	
	m1 = 0;
	for(iDim = 0; iDim < nDim; iDim++){
		K[iDim] = Direction[iDim];

		m1 += K[iDim] * K[iDim];
	}
	
	for(iDim = 0; iDim < nDim; iDim++)
		K[iDim] /= sqrt(m1);
		
	m2 = 0;
	for(iDim = 0; iDim < nDim; iDim++)
		m2 += (A2[iDim] - A1[iDim]) * K[iDim];

	m1 = 0;
	for(iDim = 0; iDim < nDim; iDim++){
		
		I[iDim] = (A2[iDim] - A1[iDim]) - m2 * K[iDim];
		
		m1 += I[iDim] * I[iDim];
	}
	
	for(iDim = 0; iDim < nDim; iDim++)
		I[iDim] /= sqrt(m1);

	
	// Cross product to find Y, works both for 2D and 3D because versor where initially set to zero.
	J[0] =   I[1]*(-K[2]) - (-K[1])*I[2];
	J[1] = -(I[0]*(-K[2]) - (-K[0])*I[2]);
	J[2] =   I[0]*(-K[1]) - (-K[0])*I[1];
	
	/* --- Project all points on the plane specified by Direction and change their reference frame taking A1 as origin --- */
	
	for(iDim = 0; iDim < nDim; iDim++){
		a2[0] += (A2[iDim] - A1[iDim]) * I[iDim];
		a2[1] += (A2[iDim] - A1[iDim]) * J[iDim];
		a2[2] += (A2[iDim] - A1[iDim]) * K[iDim];
		
		a3[0] += (A3[iDim] - A1[iDim]) * I[iDim];
		a3[1] += (A3[iDim] - A1[iDim]) * J[iDim];
		a3[2] += (A3[iDim] - A1[iDim]) * K[iDim];
		
		b1[0] += (B1[iDim] - A1[iDim]) * I[iDim];
		b1[1] += (B1[iDim] - A1[iDim]) * J[iDim];
		b1[2] += (B1[iDim] - A1[iDim]) * K[iDim];
		
		b2[0] += (B2[iDim] - A1[iDim]) * I[iDim];
		b2[1] += (B2[iDim] - A1[iDim]) * J[iDim];
		b2[2] += (B2[iDim] - A1[iDim]) * K[iDim];
		
		b3[0] += (B3[iDim] - A1[iDim]) * I[iDim];
		b3[1] += (B3[iDim] - A1[iDim]) * J[iDim];
		b3[2] += (B3[iDim] - A1[iDim]) * K[iDim];		
	}
	/*
	for(iDim = 0; iDim < nDim; iDim++)
		cout << a1[iDim] << "  ";
	cout << endl;
	
	for(iDim = 0; iDim < nDim; iDim++)
		cout << a2[iDim] << "  ";
	cout << endl;
	
	for(iDim = 0; iDim < nDim; iDim++)
		cout << a3[iDim] << "  ";
	cout << endl;
	
	for(iDim = 0; iDim < nDim; iDim++)
		cout << b1[iDim] << "  ";
	cout << endl;
	
	for(iDim = 0; iDim < nDim; iDim++)
		cout << b2[iDim] << "  ";
	cout << endl;
	
	for(iDim = 0; iDim < nDim; iDim++)
		cout << b3[iDim] << "  ";
	cout << endl;
	*/
	/* --- Find a B point inside triangle A --- */
	/*
	a1[0] = 0;
	a1[1] = 0;
	a1[2] = 0;
	
	a2[0] = 1;
	a2[1] = 0;
	a2[2] = 0;
	
	a3[0] = 0;
	a3[1] = 1;
	a3[2] = 0;
	
	b1[0] = 0.5;
	b1[1] = 0;
	b1[2] = 0;
	
	b2[0] = 1.5;
	b2[1] = 0;
	b2[2] = 0;
	
	b3[0] = 0.5;
	b3[1] = 1;
	b3[2] = 0;
	
	cout << "Area  " << ComputeIntersectionArea( a1, a2, a3, b1, b2, b3 ) << endl;
	getchar();
	*/
	/* --- Compute intersection area --- */
	
	return ComputeIntersectionArea( a1, a2, a3, b1, b2, b3 );
}

su2double CSlidingmesh::ComputeIntersectionArea( su2double* P1, su2double* P2, su2double* P3, su2double* Q1, su2double* Q2, su2double* Q3 ){
	
	bool check;
	
	unsigned short iDim, iPoints, nPoints, i, j, k, count;
	unsigned short nDim, IntersectionCounter, min_theta_index;

	su2double points[16][2], IntersectionPoint[2], theta[6];
	su2double TriangleP[4][2], TriangleQ[4][2];
	su2double Area, det, dot1, dot2, ref_cross, com_cross, dtmp, min_theta;
	
	nDim  = 2;
	check = 1;
	nPoints = 0;
	
	for(iDim = 0; iDim < nDim; iDim++){
		TriangleP[0][iDim] = P1[iDim] - P1[iDim];
		TriangleP[1][iDim] = P2[iDim] - P1[iDim];
		TriangleP[2][iDim] = P3[iDim] - P1[iDim];
		TriangleP[3][iDim] = P1[iDim] - P1[iDim];
		
		TriangleQ[0][iDim] = Q1[iDim] - P1[iDim];
		TriangleQ[1][iDim] = Q2[iDim] - P1[iDim];
		TriangleQ[2][iDim] = Q3[iDim] - P1[iDim];
		TriangleQ[3][iDim] = Q1[iDim] - P1[iDim];
	}
	
	
	for( j = 0; j < 3; j++){
		if( CheckPointInsideTriangle(TriangleP[j], TriangleQ[0], TriangleQ[1], TriangleQ[2]) ){

			// Then P1 is also inside triangle Q, so store it
			for(iDim = 0; iDim < nDim; iDim++)
				points[nPoints][iDim] = TriangleP[j][iDim];

			nPoints++;		
		}
	}
	
	for( j = 0; j < 3; j++){	
		if( CheckPointInsideTriangle(TriangleQ[j], TriangleP[0], TriangleP[1], TriangleP[2]) ){

			// Then Q1 is also inside triangle P, so store it
			for(iDim = 0; iDim < nDim; iDim++)
				points[nPoints][iDim] = TriangleQ[j][iDim];

			nPoints++;		
		}
	}
		
		
	// Compute all edge intersections

	for( j = 0; j < 3; j++){

		for( i = 0; i < 3; i++){
	
				det = (TriangleP[j][0] - TriangleP[j+1][0]) * ( TriangleQ[i][1] - TriangleQ[i+1][1] ) - (TriangleP[j][1] - TriangleP[j+1][1]) * (TriangleQ[i][0] - TriangleQ[i+1][0]);
				
				if ( det != 0.0 ){
					ComputeLineIntersectionPoint( TriangleP[j], TriangleP[j+1], TriangleQ[i], TriangleQ[i+1], IntersectionPoint );
					
					dot1 = 0;
					dot2 = 0;
					for(iDim = 0; iDim < nDim; iDim++){
							dot1 += ( TriangleP[j][iDim] - IntersectionPoint[iDim] ) * ( TriangleP[j+1][iDim] - IntersectionPoint[iDim] );
							dot2 += ( TriangleQ[i][iDim] - IntersectionPoint[iDim] ) * ( TriangleQ[i+1][iDim] - IntersectionPoint[iDim] );
					}
					
					if( dot1 <= 0 && dot2 <= 0 ){ // It found one intersection

						// Store temporarily the intersection point

						for(iDim = 0; iDim < nDim; iDim++)
							points[nPoints][iDim] = IntersectionPoint[iDim];
							
						nPoints++;
					}	
				}
			
		}
	}
		
	// Remove double points, if any

	for( i = 0; i < nPoints; i++){
			for( j = i+1; j < nPoints; j++){
				if(points[j][0] == points[i][0] && points[j][1] == points[i][1]){
						for( k = j; k < nPoints-1; k++){
							points[k][0] = points[k+1][0];
							points[k][1] = points[k+1][1];
						}
						nPoints--;
						j--;
				}
			}
	}		
	
	// Re-order nodes	

	for( i = 1; i < nPoints; i++){ // Change again reference frame
		
		for(iDim = 0; iDim < nDim; iDim++)
			points[i][iDim] -= points[0][iDim];
		
		
		// Compute polar azimuth for each node but the first
		
		if (i != 0)
			theta[i] = atan2(points[i][1], points[i][0]);

	}
	
	for(iDim = 0; iDim < nDim; iDim++)
			points[0][iDim] = 0;
	
	for( i = 1; i < nPoints; i++){
		
		min_theta = theta[i];
		min_theta_index = 0;
		
		for( j = i + 1; j < nPoints; j++){ 
			
			if( theta[j] < min_theta ){
				min_theta = theta[j];
				min_theta_index = j;
			}
		}
		
		if( min_theta_index != 0 ){
			dtmp = theta[i];
			theta[i] = theta[min_theta_index];
			theta[min_theta_index] = dtmp;
			
			dtmp = points[i][0];
			points[i][0] = points[min_theta_index][0];
			points[min_theta_index][0] = dtmp;
			
			dtmp = points[i][1];
			points[i][1] = points[min_theta_index][1];
			points[min_theta_index][1] = dtmp;
		}
	}

	// compute area using cross product rule, points position are referred to the 2-dimensional, local, reference frame centered in points[0]
	
	Area = 0;
	
	if (nPoints > 2){
		for( i = 1; i < nPoints-1; i++ ){
			
			// AxBy
			Area += ( points[i][0] - points[0][0] ) * ( points[i+1][1] - points[0][1] );
			
			// AyBx
			Area -= ( points[i][1] - points[0][1] ) * ( points[i+1][0] - points[0][0] );		
		}
	}
	
	return fabs(Area)/2;
}

void CSlidingmesh::ComputeLineIntersectionPoint( su2double* A1, su2double* A2, su2double* B1, su2double* B2, su2double* IntersectionPoint ){
	
	/* --- Uses determinant rule --- */
	
	unsigned short iDim, iPoints, nPoints;
	unsigned short nDim = donor_geometry->GetnDim();
	
	su2double det;
	
	det = (A1[0] - A2[0]) * (B1[1] - B2[1]) - (A1[1] - A2[1]) * (B1[0] - B2[0]);
	
	if ( det != 0.0 ){ // else there is no intersection point
		IntersectionPoint[0] = ( ( A1[0]*A2[1] - A1[1]*A2[0] ) * ( B1[0] - B2[0] ) - ( B1[0]*B2[1] - B1[1]*B2[0] ) * ( A1[0] - A2[0] ) ) / det;
		IntersectionPoint[1] = ( ( A1[0]*A2[1] - A1[1]*A2[0] ) * ( B1[1] - B2[1] ) - ( B1[0]*B2[1] - B1[1]*B2[0] ) * ( A1[1] - A2[1] ) ) / det;
	}
	
	return;
}

bool CSlidingmesh::CheckPointInsideTriangle(su2double* Point, su2double* T1, su2double* T2, su2double* T3){

	unsigned short iDim;
	unsigned short nDim;
	unsigned short check;
	
	su2double vect1[2], vect2[2], r[2];
	su2double dot;
	
	check = 0;
	nDim  = 2;
	
	/* --- Check first edge --- */
	
	dot = 0;
	for(iDim = 0; iDim < nDim; iDim++){
		vect1[iDim] = T3[iDim] - T1[iDim];
		vect2[iDim] = T2[iDim] - T1[iDim];
		
		r[iDim] = Point[iDim] - T1[iDim];
		
		dot += vect2[iDim] * vect2[iDim];
	}
	dot = sqrt(dot);
	
	for(iDim = 0; iDim < nDim; iDim++)
		vect2[iDim] /= dot;
		
	dot = 0;
	for(iDim = 0; iDim < nDim; iDim++)
		dot += vect1[iDim] * vect2[iDim];
	
	for(iDim = 0; iDim < nDim; iDim++)
		vect1[iDim] = T3[iDim] - (T1[iDim] + dot * vect2[iDim]);
		
	dot = 0;
	for(iDim = 0; iDim < nDim; iDim++)
		dot += vect1[iDim] * r[iDim];
		
	if (dot >= 0)
		check++;

	/* --- Check second edge --- */
		
	dot = 0;
	for(iDim = 0; iDim < nDim; iDim++){
		vect1[iDim] = T1[iDim] - T2[iDim];
		vect2[iDim] = T3[iDim] - T2[iDim];
		
		r[iDim] = Point[iDim] - T2[iDim];
		
		dot += vect2[iDim] * vect2[iDim];
	}
	dot = sqrt(dot);
	
	for(iDim = 0; iDim < nDim; iDim++)
		vect2[iDim] /= dot;
		
	dot = 0;
	for(iDim = 0; iDim < nDim; iDim++)
		dot += vect1[iDim] * vect2[iDim];
	
	for(iDim = 0; iDim < nDim; iDim++)
		vect1[iDim] = T1[iDim] - (T2[iDim] + dot * vect2[iDim]);
		
	dot = 0;
	for(iDim = 0; iDim < nDim; iDim++)
		dot += vect1[iDim] * r[iDim];
		
	if (dot >= 0)
		check++;

	/* --- Check third edge --- */
	
	dot = 0;
	for(iDim = 0; iDim < nDim; iDim++){
		vect1[iDim] = T2[iDim] - T3[iDim];
		vect2[iDim] = T1[iDim] - T3[iDim];
		
		r[iDim] = Point[iDim] - T3[iDim];
		
		dot += vect2[iDim] * vect2[iDim];
	}
	dot = sqrt(dot);
	
	for(iDim = 0; iDim < nDim; iDim++)
		vect2[iDim] /= dot;
		
	dot = 0;
	for(iDim = 0; iDim < nDim; iDim++)
		dot += vect1[iDim] * vect2[iDim];
	
	for(iDim = 0; iDim < nDim; iDim++)
		vect1[iDim] = T2[iDim] - (T3[iDim] + dot * vect2[iDim]);
		
	dot = 0;
	for(iDim = 0; iDim < nDim; iDim++)
		dot += vect1[iDim] * r[iDim];
		
	if (dot >= 0)
		check++;

	return (check == 3);
		
}

bool CSlidingmesh::CheckPointOwner(CGeometry *geometry, su2double *QueryPoint, unsigned long VertexID, unsigned long markID){

	int jPoint, iPoint, nEdges, iEdge, iNode, EdgeIndex, checkTargetOwner, nEdgesOnMarker;

	
	unsigned short iDim;
	unsigned short nDim = geometry->GetnDim();
  
	su2double EdgeCG[3], *EdgeNormal;
	su2double dotproduct;

	nEdges = geometry->node[VertexID]->GetnPoint();

	checkTargetOwner = 0;
	nEdgesOnMarker   = 0;

	for (iEdge = 0; iEdge < nEdges; iEdge++){

		EdgeIndex  = geometry->node[VertexID]->GetEdge(iEdge);
		EdgeNormal = geometry->edge[EdgeIndex]->GetNormal();
		
		if( VertexID == geometry->edge[EdgeIndex]->GetNode(0) ){
			iPoint = geometry->edge[EdgeIndex]->GetNode(0);
			jPoint = geometry->edge[EdgeIndex]->GetNode(1);
		}
		else{
			jPoint = geometry->edge[EdgeIndex]->GetNode(0);
			iPoint = geometry->edge[EdgeIndex]->GetNode(1);
			
			for(iDim = 0; iDim < nDim; iDim++)
				EdgeNormal[iDim] *= -1;
		}
		
		
		if ( donor_geometry->node[jPoint]->GetVertex(markID) != -1 ){

			nEdgesOnMarker++;
			
			bool face_on_marker, edge_on_face, element_flag;
			int iElem, jElem, iFace, iNode, node_Idx, elem_Idx, nNodes, nFaces;
			su2double face_CG[3];
			
			dotproduct = 0;
			
			if(nDim == 3){
				
				/* Find boundary faces sharing the edge */
				element_flag = false;
				
				for (iElem = 0; iElem < geometry->node[iPoint]->GetnElem(); iElem++){
					elem_Idx = geometry->node[iPoint]->GetElem(iElem);
		
					//cout << elem_Idx << "  " ;
					
					nFaces = geometry->elem[elem_Idx]->GetnFaces();
					
					for (iFace = 0; iFace < nFaces; iFace++){
						face_on_marker = true;
						edge_on_face = false;
						
						nNodes = geometry->elem[elem_Idx]->GetnNodesFace(iFace);
						
						for (iNode = 0; iNode < nNodes; iNode++){
							
							/*--- Local index of the node on face --*/
							node_Idx = geometry->elem[elem_Idx]->GetNode( geometry->elem[elem_Idx]->GetFaces(iFace, iNode) );
							
							if(node_Idx == jPoint);
								edge_on_face = true;
								
							face_on_marker = (face_on_marker and (geometry->node[node_Idx]->GetVertex(markID) != -1));
						}
						
						if ( face_on_marker && edge_on_face ){
							
							for(iDim = 0; iDim < nDim; iDim++)
								face_CG[iDim] = 0;
							
							for (iNode = 0; iNode < nNodes; iNode++){
								node_Idx = geometry->elem[elem_Idx]->GetNode( geometry->elem[elem_Idx]->GetFaces(iFace, iNode) );
								for(iDim = 0; iDim < nDim; iDim++)
									face_CG[iDim] += geometry->node[node_Idx]->GetCoord(iDim);
							}
							
							for(iDim = 0; iDim < nDim; iDim++)
								face_CG[iDim] /= nNodes;
													
							break;
						}
						
						if ( face_on_marker && edge_on_face )
							break;
					}
					
					if ( face_on_marker && edge_on_face )
							break;

				}
				
				for(iDim = 0; iDim < nDim; iDim++)
					dotproduct += (QueryPoint[iDim] - face_CG[iDim]) * EdgeNormal[iDim];

			}
			else{
				for(iDim = 0; iDim < nDim; iDim++)
					dotproduct += (QueryPoint[iDim] - geometry->edge[EdgeIndex]->GetCG(iDim)) * EdgeNormal[iDim];
			}
			
			/* if the target node is on the inner side of edge face, then the scalar product between (X-edgeCG) and edge normal must be negative */
			/* In realta non dovrebbe essere il CG dell'edge ma il cg di una delle facce di boundary adiacenti all'edge */
			
			if( dotproduct < 0 ) // target point on the inner side of edge
				checkTargetOwner++;

			}
	}

	if( checkTargetOwner == nEdgesOnMarker )
		return true;
	else
		return false;

}
