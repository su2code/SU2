/*!
 * \file interpolation_structure.cpp
 * \brief Main subroutines used by SU2_FSI
 * \author H. Kline
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

CInterpolator::CInterpolator(void) {

  nZone = 0;
  Geometry = NULL;

  donor_geometry  = NULL;
  target_geometry = NULL;

  donorZone  = 0;
  targetZone = 0;

  Buffer_Receive_nVertex_Donor     = NULL;
  Buffer_Receive_nFace_Donor       = NULL;
  Buffer_Receive_nFaceNodes_Donor  = NULL;
  Buffer_Send_nVertex_Donor        = NULL;
  Buffer_Send_nFace_Donor          = NULL;
  Buffer_Send_nFaceNodes_Donor     = NULL;
  Buffer_Receive_GlobalPoint       = NULL;
  Buffer_Send_GlobalPoint          = NULL;
  Buffer_Send_FaceIndex            = NULL;
  Buffer_Receive_FaceIndex         = NULL;
  Buffer_Send_FaceNodes            = NULL;
  Buffer_Receive_FaceNodes         = NULL;
  Buffer_Send_FaceProc             = NULL;
  Buffer_Receive_FaceProc          = NULL;

  Buffer_Send_Coord                = NULL;
  Buffer_Send_Normal               = NULL;
  Buffer_Receive_Coord             = NULL;
  Buffer_Receive_Normal            = NULL;

}

CInterpolator::~CInterpolator(void) {

  //if (Buffer_Receive_nVertex_Donor!= NULL) delete[] Buffer_Receive_nVertex_Donor;
}


CInterpolator::CInterpolator(CGeometry ***geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone) {

  /* Store pointers*/
  Geometry = geometry_container;

  donorZone  = iZone;
  targetZone = jZone;

  donor_geometry  = geometry_container[donorZone][MESH_0];
  target_geometry = geometry_container[targetZone][MESH_0];

  /*--- Initialize transfer coefficients between the zones ---*/
    /* Since this is a virtual function, call it in the child class constructor  */
  //Set_TransferCoeff(targetZone,donorZone,config);
  /*--- Initialize transfer coefficients between the zones ---*/
  //Set_TransferCoeff(Zones,config);

  //Buffer_Receive_nVertex_Donor = NULL;

}

inline void CInterpolator::Set_TransferCoeff(CConfig **config) { }

void CInterpolator::Determine_ArraySize(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim) {
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
    if (donor_geometry->node[iPointDonor]->GetDomain()) {
      nLocalVertex_Donor++;
      if (faces) {
        /*--- On Donor geometry also communicate face info ---*/
        if (nDim==3) {
          for (jElem=0; jElem<donor_geometry->node[iPointDonor]->GetnElem(); jElem++) {
            donor_elem = donor_geometry->node[iPointDonor]->GetElem(jElem);
            nFaces = donor_geometry->elem[donor_elem]->GetnFaces();
            for (iFace=0; iFace<nFaces; iFace++) {
              face_on_marker=true;
              nNodes = donor_geometry->elem[donor_elem]->GetnNodesFace(iFace);
              for (iDonor=0; iDonor<nNodes; iDonor++) {
                /*--- Local index of the node on face --*/
                inode = donor_geometry->elem[donor_elem]->GetFaces(iFace, iDonor);
                jPoint = donor_geometry->elem[donor_elem]->GetNode(inode);
                face_on_marker = (face_on_marker && (donor_geometry->node[jPoint]->GetVertex(markDonor) !=-1));
              }
              if (face_on_marker ) {
                nLocalFace_Donor++;
                nLocalFaceNodes_Donor+=nNodes;
              }
            }
          }
        }
        else {
          /*--- in 2D we use the edges ---*/
          nNodes=2;
          nFaces = donor_geometry->node[iPointDonor]->GetnPoint();
          for (iFace=0; iFace<nFaces; iFace++) {
            face_on_marker=true;
            for (iDonor=0; iDonor<nNodes; iDonor++) {
              inode = donor_geometry->node[iPointDonor]->GetEdge(iFace);
              jPoint = donor_geometry->edge[inode]->GetNode(iDonor);
              face_on_marker = (face_on_marker && (donor_geometry->node[jPoint]->GetVertex(markDonor) !=-1));
            }
            if (face_on_marker ) {
              nLocalFace_Donor++;
              nLocalFaceNodes_Donor+=nNodes;
            }
          }
        }
      }
    }
  }

  Buffer_Send_nVertex_Donor[0] = nLocalVertex_Donor;
  if (faces) {
    Buffer_Send_nFace_Donor[0] = nLocalFace_Donor;
    Buffer_Send_nFaceNodes_Donor[0] = nLocalFaceNodes_Donor;
  }

  /*--- Send Interface vertex information --*/
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalVertex_Donor, &MaxLocalVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_nVertex_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  if (faces) {
    SU2_MPI::Allreduce(&nLocalFace_Donor, &nGlobalFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFace_Donor, &MaxFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &nGlobalFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &MaxFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFace_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MaxFace_Donor++;
  }
#else
  MaxLocalVertex_Donor    = nLocalVertex_Donor;
  Buffer_Receive_nVertex_Donor[0] = Buffer_Send_nVertex_Donor[0];
  if (faces) {
    nGlobalFace_Donor       = nLocalFace_Donor;
    nGlobalFaceNodes_Donor  = nLocalFaceNodes_Donor;
    MaxFaceNodes_Donor      = nLocalFaceNodes_Donor;
    MaxFace_Donor           = nLocalFace_Donor+1;
    Buffer_Receive_nFace_Donor[0] = Buffer_Send_nFace_Donor[0];
    Buffer_Receive_nFaceNodes_Donor[0] = Buffer_Send_nFaceNodes_Donor[0];
  }
#endif

}

void CInterpolator::Collect_VertexInfo(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim) {
  unsigned long nLocalVertex_Donor = 0;
  unsigned long iVertex, iPointDonor = 0, iVertexDonor, nBuffer_Coord, nBuffer_Point;
  unsigned short iDim;
  /* Only needed if face data is also collected */
  su2double  *Normal;

#ifdef HAVE_MPI
  int rank;
  int nProcessor = SINGLE_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif


  for (iVertex = 0; iVertex < MaxLocalVertex_Donor; iVertex++) {
    Buffer_Send_GlobalPoint[iVertex] = 0;
    for (iDim = 0; iDim < nDim; iDim++) {
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

      if (faces) {
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
  if (faces) {
    SU2_MPI::Allgather(Buffer_Send_Normal, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Normal, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
  }
#else
  for (iVertex = 0; iVertex < nBuffer_Coord; iVertex++)
    Buffer_Receive_Coord[iVertex] = Buffer_Send_Coord[iVertex];

  for (iVertex = 0; iVertex < nBuffer_Point; iVertex++)
    Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];

  if (faces) {
    for (iVertex = 0; iVertex < nBuffer_Coord; iVertex++)
      Buffer_Receive_Normal[iVertex] = Buffer_Send_Normal[iVertex];
  }
#endif
}

int CInterpolator::Find_InterfaceMarker(CConfig *config, unsigned short val_marker_interface) {
    
  unsigned short nMarker = config->GetnMarker_All();
  unsigned short iMarker;

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    /*--- If the tag GetMarker_All_FSIinterface(iMarker) equals the index we are looping at ---*/
    if (config->GetMarker_All_FSIinterface(iMarker) == val_marker_interface ) {

      /*--- We have identified the identifier for the interface marker ---*/
      return iMarker;
    }
  }

  return -1;
}


/* Nearest Neighbor Interpolator */
CNearestNeighbor::CNearestNeighbor(void):  CInterpolator() { }

CNearestNeighbor::CNearestNeighbor(CGeometry ***geometry_container, CConfig **config,  unsigned int iZone, unsigned int jZone) :  CInterpolator(geometry_container, config, iZone, jZone) {

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

}

CNearestNeighbor::~CNearestNeighbor() {}

void CNearestNeighbor::Set_TransferCoeff(CConfig **config) {

  int iProcessor, pProcessor, nProcessor;
  int markDonor, markTarget, Target_check, Donor_check;

  unsigned short iDim, nDim, iMarkerInt, nMarkerInt, iDonor;    

  unsigned long nVertexDonor, nVertexTarget, Point_Target, jVertex, iVertexTarget;
  unsigned long Global_Point_Donor, pGlobalPoint=0;

  su2double *Coord_i, Coord_j[3], dist, mindist, maxdist;

#ifdef HAVE_MPI

  int rank = MASTER_NODE;
  int *Buffer_Recv_mark=NULL, iRank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  if (rank == MASTER_NODE) 
    Buffer_Recv_mark = new int[nProcessor];

#else

  nProcessor = SINGLE_NODE;

#endif

  /*--- Initialize variables --- */
  
  nMarkerInt = (int) ( config[donorZone]->GetMarker_n_FSIinterface() / 2 );
  
  nDim = donor_geometry->GetnDim();

  iDonor = 0;
  
  Buffer_Receive_nVertex_Donor = new unsigned long [nProcessor];


  /*--- Cycle over nMarkersInt interface to determine communication pattern ---*/

  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);
      
    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    #ifdef HAVE_MPI
    
    Donor_check  = -1;
    Target_check = -1;
        
    /*--- We gather a vector in MASTER_NODE to determines whether the boundary is not on the processor because of the partition or because the zone does not include it ---*/
    
    SU2_MPI::Gather(&markDonor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Donor_check = Buffer_Recv_mark[iRank];
          break;
        }
    
    SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    
    SU2_MPI::Gather(&markTarget, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Target_check = Buffer_Recv_mark[iRank];
          break;
        }

    SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
        
    #else
    Donor_check  = markDonor;
    Target_check = markTarget;  
    #endif
    
    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if(Target_check == -1 || Donor_check == -1)
      continue;

    if(markDonor != -1)
      nVertexDonor  = donor_geometry->GetnVertex( markDonor );
    else
      nVertexDonor  = 0;
    
    if(markTarget != -1)
      nVertexTarget = target_geometry->GetnVertex( markTarget );
    else
      nVertexTarget  = 0;
    
    Buffer_Send_nVertex_Donor  = new unsigned long [ 1 ];

    /* Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor */
    Determine_ArraySize(false, markDonor, markTarget, nVertexDonor, nDim);

    Buffer_Send_Coord          = new su2double     [ MaxLocalVertex_Donor * nDim ];
    Buffer_Send_GlobalPoint    = new unsigned long [ MaxLocalVertex_Donor ];
    Buffer_Receive_Coord       = new su2double     [ nProcessor * MaxLocalVertex_Donor * nDim ];
    Buffer_Receive_GlobalPoint = new unsigned long [ nProcessor * MaxLocalVertex_Donor ];

    /*-- Collect coordinates, global points, and normal vectors ---*/
    Collect_VertexInfo( false, markDonor, markTarget, nVertexDonor, nDim );

    /*--- Compute the closest point to a Near-Field boundary point ---*/
    maxdist = 0.0;

    for (iVertexTarget = 0; iVertexTarget < nVertexTarget; iVertexTarget++) {

      Point_Target = target_geometry->vertex[markTarget][iVertexTarget]->GetNode();

      if ( target_geometry->node[Point_Target]->GetDomain() ) {

        target_geometry->vertex[markTarget][iVertexTarget]->SetnDonorPoints(1);
        target_geometry->vertex[markTarget][iVertexTarget]->Allocate_DonorInfo(); // Possible meme leak?

        /*--- Coordinates of the boundary point ---*/
        Coord_i = target_geometry->node[Point_Target]->GetCoord();

        mindist    = 1E6; 
        pProcessor = 0;

        /*--- Loop over all the boundaries to find the pair ---*/

        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (jVertex = 0; jVertex < MaxLocalVertex_Donor; jVertex++) {
            Global_Point_Donor = iProcessor*MaxLocalVertex_Donor+jVertex;

            /*--- Compute the dist ---*/
            dist = 0.0; 
            for (iDim = 0; iDim < nDim; iDim++) {
              Coord_j[iDim] = Buffer_Receive_Coord[ Global_Point_Donor*nDim+iDim];
              dist += pow(Coord_j[iDim] - Coord_i[iDim], 2.0);
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
        target_geometry->vertex[markTarget][iVertexTarget]->SetDonorCoeff(iDonor, 1.0);
      }
    }

    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_GlobalPoint;
    
    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_GlobalPoint;

    delete[] Buffer_Send_nVertex_Donor;

  }

  delete[] Buffer_Receive_nVertex_Donor;

  #ifdef HAVE_MPI
  if (rank == MASTER_NODE) 
    delete [] Buffer_Recv_mark;
  #endif
}



CIsoparametric::CIsoparametric(CGeometry ***geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone)  :  CInterpolator(geometry_container, config, iZone, jZone) {

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

  /*--- For fluid-structure interaction data interpolated with have nDim dimensions ---*/
 // InitializeData(Zones,nDim);
}

CIsoparametric::~CIsoparametric() {}

void CIsoparametric::Set_TransferCoeff(CConfig **config) {
  unsigned long iVertex, jVertex;
  unsigned long  dPoint, inode, jElem, nElem;
  unsigned short iDim, iDonor=0, iFace;

  unsigned short nDim = donor_geometry->GetnDim();

  unsigned short nMarkerInt;
  unsigned short iMarkerInt;

  int markDonor=0, markTarget=0;
  int Target_check, Donor_check;

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

  int *Buffer_Recv_mark=NULL, iRank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  if (rank == MASTER_NODE) 
    Buffer_Recv_mark = new int[nProcessor];

#else

  nProcessor = SINGLE_NODE;

#endif

  /*--- Number of markers on the FSI interface ---*/
  nMarkerInt = (config[donorZone]->GetMarker_n_FSIinterface())/2;

  /*--- For the number of markers on the interface... ---*/
  for (iMarkerInt=1; iMarkerInt <= nMarkerInt; iMarkerInt++) {
    /*--- Procedure:
    * -Loop through vertices of the aero grid
    * -Find nearest element and allocate enough space in the aero grid donor point info
    *    -set the transfer coefficient values
    */

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);
      
    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    #ifdef HAVE_MPI
    
    Donor_check  = -1;
    Target_check = -1;
        
    /*--- We gather a vector in MASTER_NODE to determines whether the boundary is not on the processor because of the partition or because the zone does not include it ---*/
    
    SU2_MPI::Gather(&markDonor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Donor_check = Buffer_Recv_mark[iRank];
          break;
        }
    
    SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    SU2_MPI::Gather(&markTarget, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Target_check = Buffer_Recv_mark[iRank];
          break;
        }

    SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
        
    #else
    Donor_check  = markDonor;
    Target_check = markTarget;  
    #endif
    
    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if(Target_check == -1 || Donor_check == -1)
      continue;

    if(markDonor != -1)
      nVertexDonor  = donor_geometry->GetnVertex( markDonor );
    else
      nVertexDonor  = 0;

    if(markTarget != -1)
      nVertexTarget = target_geometry->GetnVertex( markTarget );
    else
      nVertexTarget  = 0;
    
    Buffer_Send_nVertex_Donor    = new unsigned long [1];
    Buffer_Send_nFace_Donor      = new unsigned long [1];
    Buffer_Send_nFaceNodes_Donor = new unsigned long [1];

    Buffer_Receive_nVertex_Donor    = new unsigned long [nProcessor];
    Buffer_Receive_nFace_Donor      = new unsigned long [nProcessor];
    Buffer_Receive_nFaceNodes_Donor = new unsigned long [nProcessor];

    /* Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor */
    Determine_ArraySize(true, markDonor, markTarget, nVertexDonor, nDim);

    Buffer_Send_Coord       = new su2double [MaxLocalVertex_Donor*nDim];
    Buffer_Send_Normal      = new su2double [MaxLocalVertex_Donor*nDim];
    Buffer_Send_GlobalPoint = new unsigned long [MaxLocalVertex_Donor];

    Buffer_Receive_Coord       = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
    Buffer_Receive_Normal      = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
    Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Donor];

    /*-- Collect coordinates, global points, and normal vectors ---*/
    Collect_VertexInfo(true, markDonor,markTarget,nVertexDonor,nDim);

    Buffer_Send_FaceIndex    = new unsigned long[MaxFace_Donor];
    Buffer_Send_FaceNodes    = new unsigned long[MaxFaceNodes_Donor];
    Buffer_Send_FaceProc     = new unsigned long[MaxFaceNodes_Donor];

    Buffer_Receive_FaceIndex = new unsigned long[MaxFace_Donor*nProcessor];
    Buffer_Receive_FaceNodes = new unsigned long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_FaceProc  = new unsigned long[MaxFaceNodes_Donor*nProcessor];

    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;

    /*--- Collect Face info ---*/

    for (iVertex = 0; iVertex < MaxFace_Donor; iVertex++) {
      Buffer_Send_FaceIndex[iVertex] = 0;
    }
    for (iVertex=0; iVertex<MaxFaceNodes_Donor; iVertex++) {
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

    for (jElem=0; jElem < nElem; jElem++) {
      if (nDim==3) {
        temp_donor = donor_geometry->node[iPointDonor]->GetElem(jElem);
        nFaces = donor_geometry->elem[temp_donor]->GetnFaces();
        for (iFace=0; iFace<nFaces; iFace++) {
          /*-- Determine whether this face/edge is on the marker --*/
          face_on_marker=true;
          nNodes = donor_geometry->elem[temp_donor]->GetnNodesFace(iFace);
          for (iDonor=0; iDonor<nNodes; iDonor++) {
            inode = donor_geometry->elem[temp_donor]->GetFaces(iFace, iDonor);
            dPoint = donor_geometry->elem[temp_donor]->GetNode(inode);
            face_on_marker = (face_on_marker && (donor_geometry->node[dPoint]->GetVertex(markDonor) !=-1));
          }

          if (face_on_marker ) {
            for (iDonor=0; iDonor<nNodes; iDonor++) {
              inode = donor_geometry->elem[temp_donor]->GetFaces(iFace, iDonor);
              dPoint = donor_geometry->elem[temp_donor]->GetNode(inode);
              // Match node on the face to the correct global index
              jGlobalPoint=donor_geometry->node[dPoint]->GetGlobalIndex();
              for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                for (jVertex = 0; jVertex < Buffer_Receive_nVertex_Donor[iProcessor]; jVertex++) {
                  if (jGlobalPoint ==Buffer_Receive_GlobalPoint[MaxLocalVertex_Donor*iProcessor+jVertex]) {
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
      else {
        /*-- Determine whether this face/edge is on the marker --*/
        face_on_marker=true;
        for (iDonor=0; iDonor<nNodes; iDonor++) {
          inode = donor_geometry->node[iPointDonor]->GetEdge(jElem);
          dPoint = donor_geometry->edge[inode]->GetNode(iDonor);
          face_on_marker = (face_on_marker && (donor_geometry->node[dPoint]->GetVertex(markDonor) !=-1));
        }
        if (face_on_marker ) {
          for (iDonor=0; iDonor<nNodes; iDonor++) {
            inode = donor_geometry->node[iPointDonor]->GetEdge(jElem);
            dPoint = donor_geometry->edge[inode]->GetNode(iDonor);
            // Match node on the face to the correct global index
            jGlobalPoint=donor_geometry->node[dPoint]->GetGlobalIndex();
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
              for (jVertex = 0; jVertex < Buffer_Receive_nVertex_Donor[iProcessor]; jVertex++) {
                if (jGlobalPoint ==Buffer_Receive_GlobalPoint[MaxLocalVertex_Donor*iProcessor+jVertex]) {
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
    SU2_MPI::Allgather(Buffer_Send_FaceProc, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceProc, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
    for (iFace=0; iFace<MaxFace_Donor; iFace++) {
      Buffer_Receive_FaceIndex[iFace] = Buffer_Send_FaceIndex[iFace];
    }
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertex++)
      Buffer_Receive_FaceNodes[iVertex] = Buffer_Send_FaceNodes[iVertex];
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertex++)
      Buffer_Receive_FaceProc[iVertex] = Buffer_Send_FaceProc[iVertex];
#endif

    /*--- Loop over the vertices on the target Marker ---*/
    for (iVertex = 0; iVertex<nVertexTarget; iVertex++) {
      mindist=1E6;
      for (unsigned short iCoeff=0; iCoeff<10; iCoeff++) {
    storeCoeff[iCoeff]=0;
      }
      Point_Target = target_geometry->vertex[markTarget][iVertex]->GetNode();

      if (target_geometry->node[Point_Target]->GetDomain()) {

    Coord_i = target_geometry->node[Point_Target]->GetCoord();
    /*---Loop over the faces previously communicated/stored ---*/
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {

      nFaces = (unsigned int)Buffer_Receive_nFace_Donor[iProcessor];

      for (iFace = 0; iFace< nFaces; iFace++) {
        /*--- ---*/

        nNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1] -
                (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace];

        su2double *X = new su2double[nNodes*nDim];
        faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
        for (iDonor=0; iDonor<nNodes; iDonor++) {
          jVertex = Buffer_Receive_FaceNodes[iDonor+faceindex]; // index which points to the stored coordinates, global points
          for (iDim=0; iDim<nDim; iDim++) {
            X[iDim*nNodes+iDonor]=
                Buffer_Receive_Coord[jVertex*nDim+iDim];
          }
        }
        jVertex = Buffer_Receive_FaceNodes[faceindex];

        for (iDim=0; iDim<nDim; iDim++) {
          Normal[iDim] = Buffer_Receive_Normal[jVertex*nDim+iDim];
        }

        /* Project point used for case where surfaces are not exactly coincident, where
         * the point is assumed connected by a rigid rod normal to the surface.
         */
        tmp = 0;
        tmp2=0;
        for (iDim=0; iDim<nDim; iDim++) {
          tmp+=Normal[iDim]*Normal[iDim];
          tmp2+=Normal[iDim]*(Coord_i[iDim]-X[iDim*nNodes]);
        }
        tmp = 1/tmp;
        tmp2 = tmp2*sqrt(tmp);
        for (iDim=0; iDim<nDim; iDim++) {
          // projection of \vec{q} onto plane defined by \vec{n} and \vec{p}:
          // \vec{q} - \vec{n} ( (\vec{q}-\vec{p} ) \cdot \vec{n})
          // tmp2 = ( (\vec{q}-\vec{p} ) \cdot \vec{N})
          // \vec{n} = \vec{N}/(|N|), tmp = 1/|N|^2
          projected_point[iDim]=Coord_i[iDim] + Normal[iDim]*tmp2*tmp;
        }

        Isoparameters(nDim, nNodes, X, projected_point,myCoeff);

        /*--- Find distance to the interpolated point ---*/
        dist = 0.0;
        for (iDim=0; iDim<nDim; iDim++) {
          Coord[iDim] = Coord_i[iDim];
          for(iDonor=0; iDonor< nNodes; iDonor++) {
            Coord[iDim]-=myCoeff[iDonor]*X[iDim*nNodes+iDonor];
          }
          dist+=pow(Coord[iDim],2.0);
        }

        /*--- If the dist is shorter than last closest (and nonzero nodes are on the boundary), update ---*/
        if (dist<mindist ) {
          /*--- update last dist ---*/
          mindist = dist;
          /*--- Store info ---*/
          donor_elem = temp_donor;
          target_geometry->vertex[markTarget][iVertex]->SetDonorElem(donor_elem); // in 2D is nearest neighbor
          target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nNodes);
          for (iDonor=0; iDonor<nNodes; iDonor++) {
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

    for (iDonor=0; iDonor<nNodes; iDonor++) {
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

  #ifdef HAVE_MPI
  if (rank == MASTER_NODE) 
    delete [] Buffer_Recv_mark;
  #endif
}

void CIsoparametric::Isoparameters(unsigned short nDim, unsigned short nDonor,
    su2double *X, su2double *xj, su2double *isoparams) {
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

  if (nDonor>2) {
    /*--- Create Matrix A: 1st row all 1's, 2nd row x coordinates, 3rd row y coordinates, etc ---*/
    /*--- Right hand side is [1, \vec{x}']'---*/
    for (iDonor=0; iDonor<nDonor; iDonor++) {
      isoparams[iDonor]=0;
      A[iDonor] = 1.0;
      for (iDim=0; iDim<n; iDim++)
        A[(iDim+1)*nDonor+iDonor]=X[iDim*nDonor+iDonor];
    }

    x[0] = 1.0;
    for (iDim=0; iDim<nDim; iDim++)
      x[iDim+1]=xj[iDim];

    /*--- Eliminate degenerate rows:
     * for example, if z constant including the z values will make the system degenerate
     * TODO: improve efficiency of this loop---*/
    test[0]=true; // always keep the 1st row
    for (iDim=1; iDim<nDim+1; iDim++) {
      // Test this row against all previous
      test[iDim]=true; // Assume that it is not degenerate
      for (k=0; k<iDim; k++) {
        tmp=0; tmp2=0;
        for (iDonor=0;iDonor<nDonor;iDonor++) {
          tmp+= A[iDim*nDonor+iDonor]*A[iDim*nDonor+iDonor];
          tmp2+=A[k*nDonor+iDonor]*A[k*nDonor+iDonor];
        }
        tmp  = pow(tmp,0.5);
        tmp2 = pow(tmp2,0.5);
        testi[k]=false;
        for (iDonor=0; iDonor<nDonor; iDonor++) {
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
    for (k=0; k<nDim+1; k++) {
      if (test[k]) {
        for (iDonor=0;iDonor<nDonor;iDonor++ ) {
          A2[nDonor*iDim+iDonor]=A[nDonor*k+iDonor];
        }
        x2[iDim]=x[k];
        iDim++;
      }
    }
    /*--- Initialize Q,R to 0 --*/
    for (k=0; k<nDonor*nDonor; k++) {
      Q[k]=0;
      R[k]=0;
    }
    /*--- TODO: make this loop more efficient ---*/
    /*--- Solve for rectangular Q1 R1 ---*/
    for (iDonor=0; iDonor<nDonor; iDonor++) {
      tmp=0;
      for (iDim=0; iDim<n; iDim++)
        tmp += (A2[iDim*nDonor+iDonor])*(A2[iDim*nDonor+iDonor]);

      R[iDonor*nDonor+iDonor]= pow(tmp,0.5);
      if (tmp>eps && iDonor<n) {
        for (iDim=0; iDim<n; iDim++)
          Q[iDim*nDonor+iDonor]=A2[iDim*nDonor+iDonor]/R[iDonor*nDonor+iDonor];
      }
      else if (tmp!=0) {
        for (iDim=0; iDim<n; iDim++)
          Q[iDim*nDonor+iDonor]=A2[iDim*nDonor+iDonor]/tmp;
      }
      for (iDim=iDonor+1; iDim<nDonor; iDim++) {
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
    for (iDonor=0; iDonor<nDonor; iDonor++) {
      for (iDim=0; iDim<n; iDim++)
        x_tmp[iDonor]+=Q[iDim*nDonor+iDonor]*x2[iDim];
    }

    /*--- solve x_tmp = R*isoparams for isoparams: upper triangular system ---*/
    for (iDonor = n-1; iDonor>=0; iDonor--) {
      if (R[iDonor*nDonor+iDonor]>eps)
        isoparams[iDonor]=x_tmp[iDonor]/R[iDonor*nDonor+iDonor];
      else
        isoparams[iDonor]=0;
      for (k=0; k<iDonor; k++)
        x_tmp[k]=x_tmp[k]-R[k*nDonor+iDonor]*isoparams[iDonor];
    }
  }
  else {
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
  if (nDonor==4) {
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
  if (nDonor<4) {
    tmp = 0.0; // value for normalization
    tmp2=0; // check for maximum value, to be used to id nearest neighbor if necessary
    k=0; // index for maximum value
    for (iDonor=0; iDonor< nDonor; iDonor++) {
      if (isoparams[iDonor]>tmp2) {
        k=iDonor;
        tmp2=isoparams[iDonor];
      }
      // [0,1]
      if (isoparams[iDonor]<0) isoparams[iDonor]=0;
      if (isoparams[iDonor]>1) isoparams[iDonor] = 1;
      tmp +=isoparams[iDonor];
    }
    if (tmp>0)
      for (iDonor=0; iDonor< nDonor; iDonor++)
        isoparams[iDonor]=isoparams[iDonor]/tmp;
    else {
      isoparams[k] = 1.0;
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
CMirror::CMirror(CGeometry ***geometry_container, CConfig **config,  unsigned int iZone, unsigned int jZone) :  CInterpolator(geometry_container, config, iZone, jZone) {

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

}

CMirror::~CMirror() {}

void CMirror::Set_TransferCoeff(CConfig **config) {
  unsigned long iVertex, jVertex;
  unsigned long iPoint;
  unsigned short iDonor=0, iFace=0, iTarget=0;

  unsigned short nMarkerInt;
  unsigned short iMarkerInt;

  int markDonor=0, markTarget=0;
  int Target_check, Donor_check;

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
  int *Buffer_Recv_mark=NULL, iRank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  if (rank == MASTER_NODE) 
    Buffer_Recv_mark = new int[nProcessor];

#else

  nProcessor = SINGLE_NODE;

#endif

  su2double *Buffer_Send_Coeff, *Buffer_Receive_Coeff;
  su2double coeff;

  /*--- Number of markers on the interface ---*/
  nMarkerInt = (config[targetZone]->GetMarker_n_FSIinterface())/2;

  /*--- For the number of markers on the interface... ---*/
  for (iMarkerInt=1; iMarkerInt <= nMarkerInt; iMarkerInt++) {
   /*--- Procedure:
    * -Loop through vertices of the aero grid
    * -Find nearest element and allocate enough space in the aero grid donor point info
    *    -set the transfer coefficient values
    */

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    #ifdef HAVE_MPI

    Donor_check  = -1;
    Target_check = -1;

    /*--- We gather a vector in MASTER_NODE to determines whether the boundary is not on the processor because of the partition or because the zone does not include it ---*/

    SU2_MPI::Gather(&markDonor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Donor_check = Buffer_Recv_mark[iRank];
          break;
        }

    SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    SU2_MPI::Gather(&markTarget, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Target_check = Buffer_Recv_mark[iRank];
          break;
        }

    SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    #else
    Donor_check  = markDonor;
    Target_check = markTarget;  
    #endif

    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if(Target_check == -1 || Donor_check == -1)
      continue;

    if(markDonor != -1)
      nVertexDonor  = donor_geometry->GetnVertex( markDonor );
    else
      nVertexDonor  = 0;

    if(markTarget != -1)
      nVertexTarget = target_geometry->GetnVertex( markTarget );
    else
      nVertexTarget  = 0;

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

    for (iVertex=0; iVertex<MaxFace_Donor; iVertex++) {
      Buffer_Send_FaceIndex[iVertex]=0;
    }
    for (iVertex=0; iVertex<MaxFaceNodes_Donor; iVertex++) {
      Buffer_Send_FaceNodes[iVertex]=0;
      //Buffer_Send_FaceProc[iVertex]=0;
      Buffer_Send_GlobalPoint[iVertex]=0;
      Buffer_Send_Coeff[iVertex]=0.0;
    }
    for (iVertex=0; iVertex<MaxFace_Donor; iVertex++) {
      Buffer_Send_FaceIndex[iVertex]=0;
    }

    Buffer_Send_FaceIndex[0]=rank*MaxFaceNodes_Donor;
    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;

    for (jVertex = 0; jVertex<nVertexDonor; jVertex++) {

      Point_Donor =donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex
      if (donor_geometry->node[Point_Donor]->GetDomain()) {
        nNodes = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();
        for (iDonor=0; iDonor<nNodes; iDonor++) {
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
    for (iFace=0; iFace<MaxFace_Donor; iFace++) {
      Buffer_Receive_FaceIndex[iFace] = Buffer_Send_FaceIndex[iFace];
    }
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertex++) {
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
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (iFace = 0; iFace < Buffer_Receive_nFace_Donor[iProcessor]; iFace++) {
            faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
            iNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1]- (unsigned int)faceindex;
            for (iTarget=0; iTarget<iNodes; iTarget++) {
              if (Global_Point == Buffer_Receive_GlobalPoint[faceindex+iTarget])
                nNodes++;
              //coeff =Buffer_Receive_Coeff[faceindex+iDonor];
            }
          }
        }

        target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nNodes);
        target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();

        iDonor = 0;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (iFace = 0; iFace < Buffer_Receive_nFace_Donor[iProcessor]; iFace++) {

            faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
            iNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1]- (unsigned int)faceindex;
            for (iTarget=0; iTarget<iNodes; iTarget++) {
              if (Global_Point == Buffer_Receive_GlobalPoint[faceindex+iTarget]) {
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

  #ifdef HAVE_MPI
  if (rank == MASTER_NODE) 
    delete [] Buffer_Recv_mark;
  #endif
}
