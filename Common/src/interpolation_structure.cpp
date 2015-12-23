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

	nZone = 0;
	Geometry = NULL;

	donor_geometry  = NULL;
	target_geometry = NULL;

	donorZone  = 0;
	targetZone = 0;

}

CInterpolator::~CInterpolator(void){}


CInterpolator::CInterpolator(CGeometry ***geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone){

  /* Store pointers*/
	Geometry = geometry_container;

	nZone = 2;

	donorZone  = iZone;
	targetZone = jZone;

	donor_geometry  = geometry_container[donorZone][MESH_0];
	target_geometry = geometry_container[targetZone][MESH_0];

  /*--- Initialize transfer coefficients between the zones ---*/
	/* Since this is a virtual function, call it in the child class constructor  */
  //Set_TransferCoeff(targetZone,donorZone,config);
  /*--- Initialize transfer coefficients between the zones ---*/
  //Set_TransferCoeff(Zones,config);

}

inline void CInterpolator::Set_TransferCoeff(CConfig **config) { }


CNearestNeighbor::CNearestNeighbor(void):  CInterpolator(){ }

/* Nearest Neighbor Interpolator */
CNearestNeighbor::CNearestNeighbor(CGeometry ***geometry_container, CConfig **config,  unsigned int iZone, unsigned int jZone) :  CInterpolator(geometry_container, config, iZone, jZone){

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

}

CNearestNeighbor::~CNearestNeighbor(){}

void CNearestNeighbor::Set_TransferCoeff(CConfig **config){

  unsigned long iVertex, jVertex;
  unsigned short iDim;
  unsigned short nDim = donor_geometry->GetnDim();

  unsigned short nMarkerInt, nMarkerDonor, nMarkerTarget;		// Number of markers on the interface, donor and target side
  unsigned short iMarkerInt, iMarkerDonor, iMarkerTarget;		// Variables for iteration over markers
  int markDonor = -1, markTarget = -1;

  unsigned long nVertexDonor = 0, nVertexTarget= 0;
  unsigned long Point_Donor = 0, Point_Target = 0;
  unsigned long Vertex_Donor = 0, Vertex_Target = 0;

  unsigned long iVertexDonor, iPointDonor = 0;
  unsigned long iVertexTarget, iPointTarget = 0;
  unsigned long pPoint = 0, Global_Point = 0;
  unsigned long jGlobalPoint = 0, pGlobalPoint = 0;
  int iProcessor, pProcessor = 0;

  unsigned long nLocalVertex_Donor = 0, nLocalVertex_Target = 0;
  unsigned long nGlobalVertex_Donor = 0, MaxLocalVertex_Donor = 0;

  unsigned long Global_Point_Donor;
  int Donor_Processor;
  unsigned short int donorindex = 0;

  /*--- Number of markers on the FSI interface ---*/
  nMarkerInt     = (config[donorZone]->GetMarker_n_FSIinterface())/2;
  nMarkerTarget  = target_geometry->GetnMarker();
  nMarkerDonor   = donor_geometry->GetnMarker();

  su2double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist;

  int rank = MASTER_NODE;
  int nProcessor = SINGLE_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

  // For the markers on the interface
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {

	  markDonor = -1;
	  markTarget = -1;

	  /*--- On the donor side ---*/

	  for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++){
		  /*--- If the tag GetMarker_All_FSIinterface(iMarkerDonor) equals the index we are looping at ---*/
		  if (config[donorZone]->GetMarker_All_FSIinterface(iMarkerDonor) == iMarkerInt ){
			  /*--- We have identified the identifier for the structural marker ---*/
			  markDonor = iMarkerDonor;
			  /*--- Store the number of local points that belong to markDonor ---*/
			  nVertexDonor = donor_geometry->GetnVertex(iMarkerDonor);
			  break;
		  }
		  else {
			  /*--- If the tag hasn't matched any tag within the donor markers ---*/
			  markDonor = -1;
			  nVertexDonor = 0;
		  }
	  }

	  /*--- On the target side ---*/
	  for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++){
		  /*--- If the tag GetMarker_All_FSIinterface(iMarkerFlow) equals the index we are looping at ---*/
		  if (config[targetZone]->GetMarker_All_FSIinterface(iMarkerTarget) == iMarkerInt ){
			  /*--- We have identified the identifier for the target marker ---*/
			  markTarget = iMarkerTarget;
			  /*--- Store the number of local points that belong to markTarget ---*/
			  nVertexTarget = target_geometry->GetnVertex(iMarkerTarget);
			  break;
		  }
		  else {
			  /*--- If the tag hasn't matched any tag within the Flow markers ---*/
			  nVertexTarget = 0;
			  markTarget = -1;
		  }
	  }

	  unsigned long *Buffer_Send_nVertex_Donor = new unsigned long [1];
	  unsigned long *Buffer_Receive_nVertex_Donor = new unsigned long [nProcessor];

	  unsigned long *Buffer_Send_nVertex_Target = new unsigned long [1];
	  unsigned long *Buffer_Receive_nVertex_Target = new unsigned long [nProcessor];

	  nLocalVertex_Donor = 0;
	  for (iVertex = 0; iVertex < nVertexDonor; iVertex++) {
		iPointDonor = donor_geometry->vertex[markDonor][iVertex]->GetNode();
		if (donor_geometry->node[iPointDonor]->GetDomain()) nLocalVertex_Donor++;
	  }

	  nLocalVertex_Target = 0;
	  for (iVertex = 0; iVertex < nVertexTarget; iVertex++) {
		iPointTarget = target_geometry->vertex[markTarget][iVertex]->GetNode();
		if (target_geometry->node[iPointTarget]->GetDomain()) nLocalVertex_Target++;
	  }

	  Buffer_Send_nVertex_Donor[0] = nLocalVertex_Donor;
	  Buffer_Send_nVertex_Target[0] = nLocalVertex_Target;

	  /*--- Send Interface vertex information --*/
#ifdef HAVE_MPI
	  SU2_MPI::Allreduce(&nLocalVertex_Donor, &nGlobalVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	  SU2_MPI::Allreduce(&nLocalVertex_Donor, &MaxLocalVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
	  SU2_MPI::Allgather(Buffer_Send_nVertex_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
	  SU2_MPI::Allgather(Buffer_Send_nVertex_Target, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex_Target, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
	  nGlobalVertex_Donor 		= nLocalVertex_Donor;
	  MaxLocalVertex_Donor 		= nLocalVertex_Donor;
	  Buffer_Receive_nVertex_Donor[0] = Buffer_Send_nVertex_Donor[0];
	  Buffer_Receive_nVertex_Target[0] = Buffer_Send_nVertex_Target[0];
#endif

	  su2double *Buffer_Send_Coord = new su2double [MaxLocalVertex_Donor*nDim];
	  unsigned long *Buffer_Send_GlobalPoint = new unsigned long [MaxLocalVertex_Donor];

	  su2double *Buffer_Receive_Coord = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
	  unsigned long *Buffer_Receive_Point = new unsigned long [nProcessor*MaxLocalVertex_Donor];
	  unsigned long *Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Donor];

	  unsigned long nBuffer_Coord = MaxLocalVertex_Donor*nDim;
	  unsigned long nBuffer_Point = MaxLocalVertex_Donor;

	  for (iVertex = 0; iVertex < MaxLocalVertex_Donor; iVertex++) {
      Buffer_Send_GlobalPoint[iVertex] = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
	  }

	  /*--- Copy coordinates and point to the auxiliar vector --*/
	  nLocalVertex_Donor = 0;
		for (iVertexDonor = 0; iVertexDonor < nVertexDonor; iVertexDonor++) {
		  iPointDonor = donor_geometry->vertex[markDonor][iVertexDonor]->GetNode();
		  if (donor_geometry->node[iPointDonor]->GetDomain()) {
			Buffer_Send_GlobalPoint[nLocalVertex_Donor] = donor_geometry->node[iPointDonor]->GetGlobalIndex();
			for (iDim = 0; iDim < nDim; iDim++)
			  Buffer_Send_Coord[nLocalVertex_Donor*nDim+iDim] = donor_geometry->node[iPointDonor]->GetCoord(iDim);
			nLocalVertex_Donor++;
		  }
		}

#ifdef HAVE_MPI
	  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
	  SU2_MPI::Allgather(Buffer_Send_GlobalPoint, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
	  for (iVertex = 0; iVertex < nBuffer_Coord; iVertex++){
		  Buffer_Receive_Coord[iVertex] = Buffer_Send_Coord[iVertex];
	  }
	  for (iVertex = 0; iVertex < nBuffer_Point; iVertex++){
		  Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];
	  }
#endif

	  /*--- Compute the closest point to a Near-Field boundary point ---*/
	  maxdist = 0.0;
	  for (iVertexTarget = 0; iVertexTarget < nVertexTarget; iVertexTarget++) {

		  Point_Target = target_geometry->vertex[markTarget][iVertexTarget]->GetNode();

		  if (target_geometry->node[Point_Target]->GetDomain()) {

			  target_geometry->vertex[markTarget][iVertexTarget]->SetnDonorPoints(1);
			  target_geometry->vertex[markTarget][iVertexTarget]->Allocate_DonorInfo();

			  /*--- Coordinates of the boundary point ---*/
			  Coord_i = target_geometry->node[Point_Target]->GetCoord();
			  mindist = 1E6; pProcessor = 0; pPoint = 0;

			  /*--- Loop over all the boundaries to find the pair ---*/
			  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){
				  for (jVertex = 0; jVertex < Buffer_Receive_nVertex_Donor[iProcessor]; jVertex++) {
					  Point_Donor = Buffer_Receive_Point[iProcessor*MaxLocalVertex_Donor+jVertex];
					  Global_Point_Donor = Buffer_Receive_GlobalPoint[iProcessor*MaxLocalVertex_Donor+jVertex];

					  /*--- Compute the dist ---*/
					  dist = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
						  Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_Donor+jVertex)*nDim+iDim];
						  dist += pow(Coord_j[iDim]-Coord_i[iDim],2.0);
					  }

					  if (dist < mindist) {
						  mindist = dist; pProcessor = iProcessor; pGlobalPoint = Global_Point_Donor;
					  }

					  if (dist == 0.0) break;
				  }
			  }

			  /*--- Store the value of the pair ---*/
			  maxdist = max(maxdist, mindist);
			  target_geometry->vertex[markTarget][iVertexTarget]->SetInterpDonorPoint(donorindex, pGlobalPoint);
			  target_geometry->vertex[markTarget][iVertexTarget]->SetInterpDonorProcessor(donorindex, pProcessor);
			  target_geometry->vertex[markTarget][iVertexTarget]->SetDonorCoeff(donorindex,1.0);

//	      unsigned long gpoint = target_geometry->vertex[markTarget][iVertexTarget]->GetNode();
//	      cout <<" Nearest Neighbor for target g.i " << target_geometry->node[gpoint]->GetGlobalIndex() <<" is "<< pGlobalPoint << "; d = " << mindist<< endl;
		  }
	  }

	  delete[] Buffer_Send_Coord;
	  delete[] Buffer_Send_GlobalPoint;

	  delete[] Buffer_Receive_Coord;
	  delete[] Buffer_Receive_Point;
	  delete[] Buffer_Receive_GlobalPoint;

	  delete[] Buffer_Send_nVertex_Donor;
	  delete[] Buffer_Receive_nVertex_Donor;
	  delete[] Buffer_Send_nVertex_Target;
	  delete[] Buffer_Receive_nVertex_Target;

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
  unsigned long jPoint, inode, jElem, iNearestNode=0, iNearestVertex=0;
  unsigned short iDim, donorindex=0, iFace;

  unsigned short nDim = donor_geometry->GetnDim();

  unsigned short nMarkerInt, nMarkerDonor, nMarkerTarget;
  unsigned short iMarkerInt, iMarkerDonor, iMarkerTarget;

  int markDonor=0, markTarget=0;

  long donor_elem=0, temp_donor=0;
  unsigned int nNodes=0;
  /*--- Restricted to 2-zone for now ---*/
  unsigned int nFaces=1; //For 2D cases, we want to look at edges, not faces, as the 'interface'
  bool face_on_marker=true;

  unsigned long nVertexDonor = 0, nVertexTarget= 0;
  unsigned long Point_Donor = 0, Point_Target = 0;
  unsigned long Vertex_Donor = 0, Vertex_Target = 0;

  unsigned long iVertexDonor, iPointDonor = 0;
  unsigned long iVertexTarget, iPointTarget = 0;
  unsigned long pPoint = 0, Global_Point = 0;
  unsigned long jGlobalPoint = 0, pGlobalPoint = 0;
  int iProcessor, pProcessor = 0;

  unsigned long nLocalVertex_Donor = 0, nLocalVertex_Target = 0;
  unsigned long nLocalFace_Donor = 0, nLocalFaceNodes_Donor=0;
  unsigned long nGlobalVertex_Donor = 0, MaxLocalVertex_Donor = 0;
  unsigned long nGlobalFace_Donor = 0,nGlobalFaceNodes_Donor = 0, MaxFace_Donor=0, MaxFaceNodes_Donor=0;

  unsigned long Global_Point_Donor;
  int Donor_Processor;

  su2double dist = 0.0, mindist=1E6, *Coord, *Coord_i, *Coord_j;
  su2double myCoeff[10]; // Maximum # of donor points
  su2double  *Normal;
  su2double projected_point[nDim];
  su2double tmp, tmp2;
  su2double storeCoeff[10];
  unsigned long storeGlobal[10];
  int storeProc[10];

  su2double *Buffer_Send_Coord, *Buffer_Send_Normal,*Buffer_Receive_Coord , *Buffer_Receive_Normal;
  unsigned long *Buffer_Send_GlobalPoint, *Buffer_Receive_GlobalPoint;
  unsigned long nBuffer_Coord, nBuffer_Point;
  unsigned long *Buffer_Send_FaceIndex, *Buffer_Receive_FaceIndex,
  *Buffer_Send_FaceNodes, *Buffer_Receive_FaceNodes;
  unsigned long *Buffer_Send_nVertex_Donor, *Buffer_Receive_nVertex_Donor,
  *Buffer_Send_nFace_Donor, *Buffer_Receive_nFace_Donor,
  *Buffer_Send_nFaceNodes_Donor, *Buffer_Receive_nFaceNodes_Donor;


  int rank = MASTER_NODE;
  int nProcessor = SINGLE_NODE;
  Coord = new su2double[nDim];

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

  /*--- Number of markers on the FSI interface ---*/
  nMarkerInt     = (config[donorZone]->GetMarker_n_FSIinterface())/2;
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
    nVertexTarget = 0;
    /*--- ... the marker markDonor ... ---*/
    for (iMarkerDonor=0; iMarkerDonor < nMarkerDonor; iMarkerDonor++){
      if ( config[donorZone]->GetMarker_All_FSIinterface(iMarkerDonor) == (iMarkerInt)){
        markDonor=iMarkerDonor;
        nVertexDonor = donor_geometry->GetnVertex(iMarkerDonor);
        break;
      }
      else{
        nVertexDonor = 0;
      }
    }

    /*--- ... the marker markTarget. ---*/
    for (iMarkerTarget=0; iMarkerTarget < nMarkerTarget; iMarkerTarget++){
      if (config[targetZone]->GetMarker_All_FSIinterface(iMarkerTarget) == (iMarkerInt)){
        markTarget=iMarkerTarget;
        nVertexTarget = 0;
        break;
      }
      else{
        nVertexTarget = 0;
      }
    }

    Buffer_Send_nVertex_Donor = new unsigned long [1];
    Buffer_Receive_nVertex_Donor = new unsigned long [nProcessor];
    Buffer_Send_nFace_Donor = new unsigned long [1];
    Buffer_Receive_nFace_Donor = new unsigned long [nProcessor];
    Buffer_Send_nFaceNodes_Donor = new unsigned long [1];
    Buffer_Receive_nFaceNodes_Donor = new unsigned long [nProcessor];

    nLocalVertex_Donor = 0;
    nLocalFaceNodes_Donor=0;
    nLocalFace_Donor=0;
    for (iVertex = 0; iVertex < nVertexDonor; iVertex++) {
      iPointDonor = donor_geometry->vertex[markDonor][iVertex]->GetNode();
      if (donor_geometry->node[iPointDonor]->GetDomain()){
        nLocalVertex_Donor++;
        /*--- On Donor geometry also communicate face info ---*/
        if (nDim==3){
          for (jElem=0; jElem<donor_geometry->node[iPointDonor]->GetnElem(); jElem++){
            temp_donor = donor_geometry->node[iPointDonor]->GetElem(jElem);
            nFaces = donor_geometry->elem[temp_donor]->GetnFaces();
            for (iFace=0; iFace<nFaces; iFace++){
              face_on_marker=true;
              nNodes = donor_geometry->elem[temp_donor]->GetnNodesFace(iFace);
              for (unsigned int ifacenode=0; ifacenode<nNodes; ifacenode++){
                /*--- Local index of the node on face --*/
                inode = donor_geometry->elem[temp_donor]->GetFaces(iFace, ifacenode);
                jPoint = donor_geometry->elem[temp_donor]->GetNode(inode);
                face_on_marker = (face_on_marker and (donor_geometry->node[jPoint]->GetVertex(markDonor) !=-1));
              }
              if (face_on_marker ){
                nLocalFace_Donor++;
                nLocalFaceNodes_Donor+=nNodes;
              }
            }
          }
        }
        else{
         /*--- in 2D we use the edges --- */
         nNodes=2;
         nFaces = donor_geometry->node[iPointDonor]->GetnPoint();
         for (iFace=0; iFace<nFaces; iFace++){
           face_on_marker=true;
           for (unsigned int ifacenode=0; ifacenode<nNodes; ifacenode++){
             inode = donor_geometry->node[iPointDonor]->GetEdge(iFace);
             jPoint = donor_geometry->edge[inode]->GetNode(ifacenode);
             face_on_marker = (face_on_marker and (donor_geometry->node[jPoint]->GetVertex(markDonor) !=-1));
           }
           if (face_on_marker ){
             nLocalFace_Donor++;
             nLocalFaceNodes_Donor+=nNodes;
           }
         }
       }
     }
    }
    Buffer_Send_nVertex_Donor[0] = nLocalVertex_Donor;
    Buffer_Send_nFace_Donor[0] = nLocalFace_Donor;
    Buffer_Send_nFaceNodes_Donor[0] = nLocalFaceNodes_Donor;


    /*--- Send Interface vertex information --*/
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalVertex_Donor, &nGlobalVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFace_Donor, &nGlobalFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &nGlobalFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalVertex_Donor, &MaxLocalVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &MaxFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFace_Donor, &MaxFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nVertex_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFace_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MaxFace_Donor++;
#else
    nGlobalVertex_Donor     = nLocalVertex_Donor;
    MaxLocalVertex_Donor    = nLocalVertex_Donor;
    nGlobalFace_Donor       = nLocalFace_Donor;
    nGlobalFaceNodes_Donor  = nLocalFaceNodes_Donor;
    MaxFace_Donor           = nLocalFace_Donor+1;
    MaxFaceNodes_Donor      = nLocalFaceNodes_Donor;
    Buffer_Receive_nVertex_Donor[0] = Buffer_Send_nVertex_Donor[0];
    Buffer_Receive_nFace_Donor[0] = Buffer_Send_nFace_Donor[0];
    Buffer_Receive_nFaceNodes_Donor[0] = Buffer_Send_nFaceNodes_Donor[0];

#endif
    Buffer_Send_Coord = new su2double [MaxLocalVertex_Donor*nDim];
    Buffer_Send_Normal = new su2double [MaxLocalVertex_Donor*nDim];
    Buffer_Send_GlobalPoint = new unsigned long [MaxLocalVertex_Donor];

    Buffer_Receive_Coord = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
    Buffer_Receive_Normal = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
    Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Donor];

    nBuffer_Coord = MaxLocalVertex_Donor*nDim;
    nBuffer_Point = MaxLocalVertex_Donor;

    Buffer_Send_FaceIndex = new unsigned long [MaxFace_Donor];
    Buffer_Receive_FaceIndex = new unsigned long[MaxFace_Donor*nProcessor];
    Buffer_Send_FaceNodes = new unsigned long [MaxFaceNodes_Donor];
    Buffer_Receive_FaceNodes = new unsigned long [MaxFaceNodes_Donor*nProcessor];

    for (iVertex = 0; iVertex < MaxLocalVertex_Donor; iVertex++) {
      Buffer_Send_GlobalPoint[iVertex] = 0;
      for (iDim = 0; iDim < nDim; iDim++){
        Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
        Buffer_Send_Normal[iVertex*nDim+iDim] = 0.0;
      }
    }
    Buffer_Send_FaceIndex[0]=rank*MaxFace_Donor;
    /*--- Copy coordinates and point to the auxiliar vector --*/
    nLocalVertex_Donor = 0;
    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;
    for (iVertexDonor = 0; iVertexDonor < nVertexDonor; iVertexDonor++) {
      iPointDonor = donor_geometry->vertex[markDonor][iVertexDonor]->GetNode();
      if (donor_geometry->node[iPointDonor]->GetDomain()) {
        Buffer_Send_GlobalPoint[nLocalVertex_Donor] = donor_geometry->node[iPointDonor]->GetGlobalIndex();
        Normal =  donor_geometry->vertex[markDonor][iVertexDonor]->GetNormal();
        for (iDim = 0; iDim < nDim; iDim++){
          Buffer_Send_Coord[nLocalVertex_Donor*nDim+iDim] = donor_geometry->node[iPointDonor]->GetCoord(iDim);
          Buffer_Send_Normal[nLocalVertex_Donor*nDim+iDim] = Normal[iDim];
        }
        nLocalVertex_Donor++;

        if (nDim==3){
          for (jElem=0; jElem<donor_geometry->node[iPointDonor]->GetnElem(); jElem++){
            temp_donor = donor_geometry->node[iPointDonor]->GetElem(jElem);
            nFaces = donor_geometry->elem[temp_donor]->GetnFaces();
            /*-- Loop over all faces on this element and find out how many are on the marker --*/
            for (iFace=0; iFace<nFaces; iFace++){
              face_on_marker=true;
              nNodes = donor_geometry->elem[temp_donor]->GetnNodesFace(iFace);
              for (unsigned int ifacenode=0; ifacenode<nNodes; ifacenode++){
                /*--- Local index of the node on face --*/
                inode = donor_geometry->elem[temp_donor]->GetFaces(iFace, ifacenode);
                jPoint = donor_geometry->elem[temp_donor]->GetNode(inode);
                face_on_marker = (face_on_marker and (donor_geometry->node[jPoint]->GetVertex(markDonor) !=-1));
              }
              if (face_on_marker ){
                //nLocalFaceNodes_Donor=0;
                for (unsigned int ifacenode=0; ifacenode<nNodes; ifacenode++){
                  inode = donor_geometry->elem[temp_donor]->GetFaces(iFace, ifacenode);
                  jPoint = donor_geometry->elem[temp_donor]->GetNode(inode);
                  jVertex = donor_geometry->node[jPoint]->GetVertex(markDonor);
                  Buffer_Send_FaceNodes[nLocalFaceNodes_Donor]=rank*MaxFaceNodes_Donor+jVertex;
                  nLocalFaceNodes_Donor++;
                }
                // index pointing to where this face's nodes start in the FaceNodes array
                Buffer_Send_FaceIndex[nLocalFace_Donor+1] = MaxFaceNodes_Donor*rank+nLocalFaceNodes_Donor;
                nLocalFace_Donor++;
              }
            }
          }
        }
        else{
          /*--- in 2D we use the edges --- */
          nNodes=2;
          nFaces = donor_geometry->node[iPointDonor]->GetnPoint();
          for (iFace=0; iFace<nFaces; iFace++){
            face_on_marker=true;
            for (unsigned int ifacenode=0; ifacenode<nNodes; ifacenode++){
              inode = donor_geometry->node[iPointDonor]->GetEdge(iFace);
              jPoint = donor_geometry->edge[inode]->GetNode(ifacenode);
              face_on_marker = (face_on_marker and (donor_geometry->node[jPoint]->GetVertex(markDonor) !=-1));
            }
            if (face_on_marker ){

              //nLocalFaceNodes_Donor=0;
              for (unsigned int ifacenode=0; ifacenode<nNodes; ifacenode++){
                inode = donor_geometry->node[iPointDonor]->GetEdge(iFace);
                jPoint = donor_geometry->edge[inode]->GetNode(ifacenode);
                jVertex = donor_geometry->node[jPoint]->GetVertex(markDonor);
                Buffer_Send_FaceNodes[nLocalFaceNodes_Donor]=rank*MaxFaceNodes_Donor+jVertex;
                nLocalFaceNodes_Donor++;
              }
              Buffer_Send_FaceIndex[nLocalFace_Donor+1] = MaxFaceNodes_Donor*rank+nLocalFaceNodes_Donor;
              nLocalFace_Donor++;
            }
          }
        }
      }
    }
    //Buffer_Send_FaceIndex[nLocalFace_Donor+1] = MaxFaceNodes_Donor*rank+nLocalFaceNodes_Donor;
#ifdef HAVE_MPI
    SU2_MPI::Allgather(Buffer_Send_Normal, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Normal, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_GlobalPoint, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    SU2_MPI::Allgather(Buffer_Send_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

#else
    for (iVertex = 0; iVertex < nBuffer_Coord; iVertex++){
      Buffer_Receive_Coord[iVertex] = Buffer_Send_Coord[iVertex];
      Buffer_Receive_Normal[iVertex] = Buffer_Send_Normal[iVertex];
    }
    for (iVertex = 0; iVertex < nBuffer_Point; iVertex++){
      Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];
    }
    for (iFace=0; iFace<MaxFace_Donor; iFace++){
      Buffer_Recieve_FaceIndex[iFace] = Buffer_Send_FaceIndex[iFace];
    }
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertx++)
      Buffer_Recieve_FaceNodes[iVertex] = Buffer_Send_FaceNodes[iVertex];

#endif


    /*--- Loop over the vertices on the target Marker ---*/
    for (iVertex = 0; iVertex<target_geometry->GetnVertex(markTarget); iVertex++) {
      mindist=1E6;
      for (unsigned short iCoeff=0; iCoeff<10; iCoeff++){
        storeCoeff[iCoeff]=0;
      }
      Point_Target = target_geometry->vertex[markTarget][iVertex]->GetNode();
      if (target_geometry->node[Point_Target]->GetDomain()) {
        Coord_i = target_geometry->vertex[markTarget][iVertex]->GetCoord();
        /*---Loop over the faces previously communicated/stored ---*/
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){
          nFaces = Buffer_Receive_nFace_Donor[iProcessor];
          for (iFace = 0; iFace< nFaces; iFace++){
            /*--- ---*/
            nNodes = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1] -
                    Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace];
            su2double X[nNodes*nDim];
            unsigned long faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace];
            for (donorindex=0; donorindex<nNodes; donorindex++){
              jVertex = Buffer_Receive_FaceNodes[donorindex+faceindex];
              for (iDim=0; iDim<nDim; iDim++){
                X[iDim*nNodes+donorindex]=
                    Buffer_Receive_Coord[jVertex*nDim+iDim];
              }
            }
            jVertex = Buffer_Receive_FaceNodes[faceindex];
            for (iDim=0; iDim<nDim; iDim++){
              Normal[iDim] = Buffer_Receive_Normal[jVertex*nDim+iDim];
            }
            // Project point
            tmp = 0;
            tmp2=0;
            for (iDim=0; iDim<nDim; iDim++){
              tmp+=Normal[iDim]*Normal[iDim];
              tmp2+=Normal[iDim]*(Coord_i[iDim]-X[iDim*nNodes]);
            }
            tmp = 1/tmp;
            tmp2 = tmp2*sqrt(tmp);
            for (iDim=0; iDim<nDim; iDim++){
              // projection of \vec{q} onto plane defined by \vec{n} and \vec{p}:
              // \vec{q} - \vec{n} ( (\vec{q}-\vec{p} ) \cdot \vec{n})
              // tmp2 = ( (\vec{q}-\vec{p} ) \cdot \vec{N})
              // \vec{n} = \vec{N}/(|N|), tmp = 1/|N|^2
              projected_point[iDim]=Coord_i[iDim] + Normal[iDim]*tmp2*tmp;
            }
            Isoparameters(myCoeff,nDim, nNodes, X, projected_point);
            /*--- Recalculate the dist using the isoparametric representation ---*/
            dist = 0.0;
            for (iDim=0; iDim<nDim; iDim++){
              Coord[iDim] = target_geometry->vertex[markTarget][iVertex]->GetCoord(iDim);
              for(donorindex=0; donorindex< nNodes; donorindex++){
                Coord[iDim]-=myCoeff[donorindex]*X[iDim*nNodes+donorindex];
              }
              dist+=pow(Coord[iDim],2.0);
            }
            /*--- If the dist is shorter than last closest (and nonzero nodes are on the boundary), update ---*/
            if (dist<mindist ){
              /*--- update last dist ---*/
              /*-- Debug output ---*/
              cout << iVertex <<  " dist " << dist << endl;
              for (iDim=0; iDim<nDim; iDim++)
                cout << target_geometry->vertex[markTarget][iVertex]->GetCoord(iDim) << " ";
              cout << endl;
              for (donorindex=0; donorindex<nNodes; donorindex++){
                jVertex = Buffer_Receive_FaceNodes[donorindex+faceindex];
                cout <<  myCoeff[donorindex] <<", " << Buffer_Receive_GlobalPoint[jVertex]<< " -> " ;
                for (iDim=0; iDim<nDim; iDim++){
                  cout << X[iDim*nNodes+donorindex] <<" ";
                }
                cout << endl;
              }
              /*--  ---*/
              mindist = dist;
              /*--- Store info ---*/
              donor_elem = temp_donor;
              target_geometry->vertex[markTarget][iVertex]->SetDonorElem(donor_elem); // in 2D is nearest neighbor
              target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nNodes);
              for (donorindex=0; donorindex<nNodes; donorindex++){
                storeCoeff[donorindex] = myCoeff[donorindex];
                jVertex = Buffer_Receive_FaceNodes[faceindex+donorindex];
                storeGlobal[donorindex] =Buffer_Receive_GlobalPoint[jVertex];
                storeProc[donorindex] = iProcessor;
              }
            }
          }
        }
        /*--- Set the appropriate amount of memory and fill ---*/
        target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();
        nNodes =target_geometry->vertex[markTarget][iVertex]->GetnDonorPoints();
        for (donorindex=0; donorindex<nNodes; donorindex++){
          target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint(donorindex,storeGlobal[donorindex]);
          target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(donorindex,storeCoeff[donorindex]);
          target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(donorindex, storeProc[donorindex]);
        }
      }
    }
    delete[] Buffer_Send_nVertex_Donor;
    delete[] Buffer_Receive_nVertex_Donor;
    delete[] Buffer_Send_nFace_Donor;
    delete[] Buffer_Receive_nFace_Donor;
    delete[] Buffer_Send_nFaceNodes_Donor;
    delete[] Buffer_Receive_nFaceNodes_Donor;
    //
    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_Normal;
    delete[] Buffer_Send_GlobalPoint;

    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_Normal;
    delete[] Buffer_Receive_GlobalPoint;

    delete[] Buffer_Send_FaceIndex;
    delete[] Buffer_Receive_FaceIndex;
    delete[] Buffer_Send_FaceNodes;
    delete[] Buffer_Receive_FaceNodes;
  }
  delete [] Coord;
}

void CIsoparametric::Isoparameters(su2double *isoparams, unsigned short nDim,
    unsigned short nDonor, su2double *X, su2double *xj){
  short m0 = nDonor, n0=nDim+1,m,n;
  short i,j,k, iedge, iDim;
  unsigned long jPoint, jPoint2;
  su2double tmp, tmp2;
  su2double x[n0], x_tmp[n0];
  su2double Q[m0*m0], R[m0*m0], A[n0*m0];
  su2double x2[n0];
  su2double eps = 1E-10;
  bool test[n0], testi[n0];
  m=m0; n=n0;
  if (nDonor>2){
    /*--- Create Matrix A: 1st row all 1's, 2nd row x coordinates, 3rd row y coordinates, etc ---*/
    /*--- Right hand side is [1, \vec{x}']'---*/
    for (i=0; i<m; i++){
      isoparams[i]=0;
      A[i]=1.0;
      for (j=0; j<n; j++)
        A[(j+1)*m+i]=X[j*m+i];
    }

    /*j,n: dimension. i,m: # neighbor point*/
    x[0]=1.0;
    for (iDim=0; iDim<nDim; iDim++)
      x[iDim+1]=xj[iDim];

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
  }
  else{
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
    j=0; // index for maximum value
    for (i=0; i< nDonor; i++){
      if (isoparams[i]>tmp2){
        j=i;
        tmp2=isoparams[i];
      }
      // [0,1]
      if (isoparams[i]<0) isoparams[i]=0;
      if (isoparams[i]>1) isoparams[i]=1;
      tmp +=isoparams[i];
    }
    if (tmp>0)
      for (i=0; i< nDonor; i++)
        isoparams[i]=isoparams[i]/tmp;
    else{
      isoparams[j]=1.0;
    }
  }

}


/* Mirror Interpolator */
CMirror::CMirror(CGeometry ***geometry_container, CConfig **config,  unsigned int iZone, unsigned int jZone) :  CInterpolator(geometry_container, config, iZone, jZone){

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

}

CMirror::~CMirror(){}

void CMirror::Set_TransferCoeff(CConfig **config){
  unsigned long iVertex, jVertex;
   unsigned long jPoint, iPoint, tPoint, nTargets;
   unsigned short iDim, donorindex=0, iFace;

   unsigned short nDim = donor_geometry->GetnDim();

   unsigned short nMarkerInt, nMarkerDonor, nMarkerTarget;
   unsigned short iMarkerInt, iMarkerDonor, iMarkerTarget;

   int markDonor=0, markTarget=0;

   long donor_elem=0, temp_donor=0;
   unsigned int nNodes=0;
   /*--- Restricted to 2-zone for now ---*/
   unsigned int nFaces=1; //For 2D cases, we want to look at edges, not faces, as the 'interface'
   bool face_on_marker=true;

   unsigned long iVertexDonor, iPointDonor = 0;
   unsigned long iVertexTarget, iPointTarget = 0;
   unsigned long pPoint = 0, Global_Point = 0;
   unsigned long jGlobalPoint = 0, pGlobalPoint = 0;
   int iProcessor, pProcessor = 0;

   unsigned long nLocalVertex_Donor = 0, nLocalVertex_Target = 0;
   unsigned long nGlobalVertex_Donor = 0, MaxLocalVertex_Donor = 0;

   unsigned long Global_Point_Donor;
   int Donor_Processor;

   int rank = MASTER_NODE;
   int nProcessor = SINGLE_NODE;

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

    /*--- ... the marker markDonor ... ---*/
    for (iMarkerDonor=0; iMarkerDonor < nMarkerDonor; iMarkerDonor++){
      if ( config[donorZone]->GetMarker_All_FSIinterface(iMarkerDonor) == (iMarkerInt)){
        markDonor=iMarkerDonor;
      }
    }

    /*--- ... the marker markTarget. ---*/
    for (iMarkerTarget=0; iMarkerTarget < nMarkerTarget; iMarkerTarget++){
      if (config[targetZone]->GetMarker_All_FSIinterface(iMarkerTarget) == (iMarkerInt)){
        markTarget=iMarkerTarget;
      }
    }

    /*--- Loop over the vertices on the target Marker ---*/
    for (iVertex = 0; iVertex<target_geometry->GetnVertex(markTarget); iVertex++) {

      iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();
      Global_Point = target_geometry->node[iPoint]->GetGlobalIndex();
      nNodes = 0;
      for (jVertex = 0; jVertex<donor_geometry->GetnVertex(markDonor); jVertex++) {
        jPoint =donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex
        nTargets = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();
        for (unsigned short iTarget=0; iTarget<nTargets; iTarget++){
          tPoint = donor_geometry->vertex[markDonor][jVertex]->GetInterpDonorPoint(iTarget);
          if (tPoint==Global_Point){
            nNodes++;
          }
        }
      }

      target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nNodes);

      if (nNodes>0){
        /*--- Set the appropriate amount of memory ---*/
        target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();
        /*--- Find the coefficient info from the donor geometry --- */
        donorindex=0;
        for (jVertex = 0; jVertex<donor_geometry->GetnVertex(markDonor); jVertex++) {
          jPoint =donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex
          nTargets = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();

          for (unsigned short iTarget=0; iTarget<nTargets; iTarget++){
            tPoint = donor_geometry->vertex[markDonor][jVertex]->GetInterpDonorPoint(iTarget);
            if (tPoint==Global_Point){
              pGlobalPoint = donor_geometry->node[jPoint]->GetGlobalIndex();
              target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint(donorindex,pGlobalPoint);
              coeff = donor_geometry->vertex[markDonor][jVertex]->GetDonorCoeff(iTarget);
              target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(donorindex,coeff);
              target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(donorindex, MASTER_NODE);
              donorindex++;
            }
          }
        }
      }
      // FOR PARALELL:
      //target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProc(donorindex,proc);
      //cout <<" myCoeff  " << myCoeff[donorindex] << " ";
      //cout << endl;
    }

  }

}
