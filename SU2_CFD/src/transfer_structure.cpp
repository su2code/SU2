/*!
 * \file transfer_structure.cpp
 * \brief Main subroutines for MPI transfer of information between zones
 * \author R. Sanchez
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

#include "../include/transfer_structure.hpp"

CTransfer::CTransfer(void) {
  
  Physical_Constants = NULL;
  Donor_Variable     = NULL;
  Target_Variable    = NULL;
  
  nVar = 0;
  
}

CTransfer::CTransfer(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) {
  
  unsigned short iVar;
  
  Physical_Constants = new su2double[val_nConst];
  Donor_Variable     = new su2double[val_nVar];
  Target_Variable    = new su2double[val_nVar];
  
  nVar = val_nVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Donor_Variable[iVar]  = 0.0;
    Target_Variable[iVar] = 0.0;
  }
  
  for (iVar = 0; iVar < val_nConst; iVar++) {
    Physical_Constants[iVar] = 0.0;
  }
  
}

CTransfer::~CTransfer(void) {
  
  if (Physical_Constants   != NULL) delete [] Physical_Constants;
  if (Donor_Variable       != NULL) delete [] Donor_Variable;
  if (Target_Variable      != NULL) delete [] Target_Variable;
  
}

void CTransfer::Scatter_InterfaceData(CSolver *donor_solution, CSolver *target_solution,
                                      CGeometry *donor_geometry, CGeometry *target_geometry,
                                      CConfig *donor_config, CConfig *target_config) {
  
  unsigned short nMarkerInt, nMarkerDonor, nMarkerTarget;       // Number of markers on the interface, donor and target side
  unsigned short iMarkerInt, iMarkerDonor, iMarkerTarget;       // Variables for iteration over markers
  int Marker_Donor = -1, Marker_Target = -1;
  int Target_check, Donor_check;
  
  unsigned long iVertex;                            // Variables for iteration over vertices and nodes
  
  unsigned short iVar;
  
  GetPhysical_Constants(donor_solution, target_solution, donor_geometry, target_geometry,
                        donor_config, target_config);
  
  unsigned long Point_Donor, Point_Target;
  
  bool fsi = donor_config->GetFSI_Simulation();
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int *Buffer_Recv_mark=NULL, iRank;
  
  if (rank == MASTER_NODE) 
    Buffer_Recv_mark = new int[size];
#endif
  
  unsigned long nLocalVertexDonor   = 0, nLocalVertexTarget   = 0;
  unsigned long MaxLocalVertexDonor = 0, MaxLocalVertexTarget = 0;
  
  unsigned long nBuffer_DonorVariables = 0, nBuffer_TargetVariables = 0;
  unsigned long nBuffer_DonorIndices   = 0, nBuffer_TargetIndices   = 0;

  unsigned long Processor_Target;
  
  int iProcessor, nProcessor = 0;
  
  /*--- Number of markers on the FSI interface ---*/
  
  nMarkerInt     = (donor_config->GetMarker_n_FSIinterface())/2;
  nMarkerTarget  = target_geometry->GetnMarker();
  nMarkerDonor   = donor_geometry->GetnMarker();
  
  nProcessor = size;
  
  /*--- Outer loop over the markers on the FSI interface: compute one by one ---*/
  /*--- The tags are always an integer greater than 1: loop from 1 to nMarkerFSI ---*/
  
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {
    
    Marker_Donor  = -1;
    Marker_Target = -1;
    
    /*--- Initialize pointer buffers inside the loop, so we can delete for each marker. ---*/
    unsigned long Buffer_Send_nVertexDonor[1], *Buffer_Recv_nVertexDonor = NULL;
    unsigned long Buffer_Send_nVertexTarget[1], *Buffer_Recv_nVertexTarget = NULL;
    
    /*--- The donor and target markers are tagged with the same index.
     *--- This is independent of the MPI domain decomposition.
     *--- We need to loop over all markers on both sides and get the number of nodes
     *--- that belong to each FSI marker for each processor ---*/
    
    /*--- On the donor side ---*/
    
    for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++) {
      /*--- If the tag GetMarker_All_FSIinterface(iMarkerDonor) equals the index we are looping at ---*/
      if ( donor_config->GetMarker_All_FSIinterface(iMarkerDonor) == iMarkerInt ) {
        Marker_Donor = iMarkerDonor;
        /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
        break;
      }
    }
    
    /*--- On the target side ---*/
    
    for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++) {
      /*--- If the tag GetMarker_All_FSIinterface(iMarkerFlow) equals the index we are looping at ---*/
      if ( target_config->GetMarker_All_FSIinterface(iMarkerTarget) == iMarkerInt ) {
        Marker_Target = iMarkerTarget;
        /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
        break;
      }
    }

    #ifdef HAVE_MPI

    Donor_check  = -1;
    Target_check = -1;

    /*--- We gather a vector in MASTER_NODE that determines if the boundary is not on the processor because of the partition or because the zone does not include it  ---*/

    SU2_MPI::Gather(&Marker_Donor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE) {
      for (iRank = 0; iRank < nProcessor; iRank++) {
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Donor_check = Buffer_Recv_mark[iRank];
          break;
        }       
      }
    }

    SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    SU2_MPI::Gather(&Marker_Target, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE) {
      for (iRank = 0; iRank < nProcessor; iRank++) {
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Target_check = Buffer_Recv_mark[iRank];
          break;
        }   
      }
    }

    SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    #else
    Donor_check  = Marker_Donor;
    Target_check = Marker_Target;   
    #endif

    if(Target_check == -1 || Donor_check == -1) {
      continue;
    }

    nLocalVertexDonor  = 0;
    nLocalVertexTarget = 0;

    if( Marker_Donor != -1 )
      nLocalVertexDonor = donor_geometry->GetnVertex(Marker_Donor);

    if( Marker_Target != -1 )
      nLocalVertexTarget = target_geometry->GetnVertex(Marker_Target);

    Buffer_Send_nVertexDonor[0] = nLocalVertexDonor;                               // Retrieve total number of vertices on Donor marker
    Buffer_Send_nVertexTarget[0] = nLocalVertexTarget;                             // Retrieve total number of vertices on Target marker
    if (rank == MASTER_NODE) Buffer_Recv_nVertexDonor = new unsigned long[size];   // Allocate memory to receive how many vertices are on each rank on the structural side
    if (rank == MASTER_NODE) Buffer_Recv_nVertexTarget = new unsigned long[size];  // Allocate memory to receive how many vertices are on each rank on the fluid side
#ifdef HAVE_MPI
    /*--- We receive MaxLocalVertexFEA as the maximum number of vertices in one single processor on the structural side---*/
    SU2_MPI::Allreduce(&nLocalVertexDonor, &MaxLocalVertexDonor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    /*--- We receive MaxLocalVertexFlow as the maximum number of vertices in one single processor on the fluid side ---*/
    SU2_MPI::Allreduce(&nLocalVertexTarget, &MaxLocalVertexTarget, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    
    /*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the structural side ---*/
    SU2_MPI::Gather(&Buffer_Send_nVertexDonor, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexDonor, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    /*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the fluid side ---*/
    SU2_MPI::Gather(&Buffer_Send_nVertexTarget, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexTarget, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
    MaxLocalVertexDonor  = nLocalVertexDonor;
    MaxLocalVertexTarget = nLocalVertexTarget;
    
    Buffer_Recv_nVertexDonor[0] = Buffer_Send_nVertexDonor[0];
    Buffer_Recv_nVertexTarget[0] = Buffer_Send_nVertexTarget[0];
    
#endif
    
    /*--- We will be gathering the structural coordinates into the master node ---*/
    /*--- Then we will distribute them using a scatter operation into the appropriate fluid processor ---*/
    nBuffer_DonorVariables = MaxLocalVertexDonor * nVar;
    nBuffer_TargetVariables = MaxLocalVertexTarget * nVar;
    
    /*--- We will be gathering donor index and donor processor (for flow -> donor = structure) ---*/
    /*--- Then we will pass on to the structural side the index (fea point) to the appropriate processor ---*/
    nBuffer_DonorIndices = 2 * MaxLocalVertexDonor;
    nBuffer_TargetIndices = MaxLocalVertexTarget;
    
    /*--- Send and Recv buffers ---*/
    
    /*--- Buffers to send and receive the variables in the donor mesh ---*/
    su2double *Buffer_Send_DonorVariables = new su2double[nBuffer_DonorVariables];
    su2double *Buffer_Recv_DonorVariables = NULL;
    
    /*--- Buffers to send and receive the indices in the donor mesh ---*/
    long *Buffer_Send_DonorIndices = new long[nBuffer_DonorIndices];
    long *Buffer_Recv_DonorIndices = NULL;
    
    /*--- Buffers to send and receive the variables in the target mesh---*/
    su2double *Buffer_Send_TargetVariables = NULL;
    su2double *Buffer_Recv_TargetVariables = new su2double[nBuffer_TargetVariables];
    
    /*--- Buffers to send and receive the target indices ---*/
    long *Buffer_Send_TargetIndices = NULL;
    long *Buffer_Recv_TargetIndices = new long[nBuffer_TargetIndices];
    
    /*--- Prepare the receive buffers (1st step) and send buffers (2nd step) on the master node only. ---*/
    
    if (rank == MASTER_NODE) {
      Buffer_Recv_DonorVariables  = new su2double[size*nBuffer_DonorVariables];
      Buffer_Recv_DonorIndices    = new long[size*nBuffer_DonorIndices];
      Buffer_Send_TargetVariables = new su2double[size*nBuffer_TargetVariables];
      Buffer_Send_TargetIndices   = new long[size*nBuffer_TargetIndices];
    }
    
    /*--- On the fluid side ---*/
    
    /*--- If this processor owns the marker we are looping at on the structural side ---*/
    
    /*--- First we initialize all of the indices and processors to -1 ---*/
    /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
    for (iVertex = 0; iVertex < nBuffer_DonorIndices; iVertex++)
      Buffer_Send_DonorIndices[iVertex] = -1;
     
    /*--- We have identified the local index of the FEA marker ---*/
    /*--- We loop over all the vertices in that marker and in that particular processor ---*/

    for (iVertex = 0; iVertex < nLocalVertexDonor; iVertex++) {

      Point_Donor = donor_geometry->vertex[Marker_Donor][iVertex]->GetNode();

      /*--- If this processor owns the node ---*/
      if (donor_geometry->node[Point_Donor]->GetDomain()) {
        Point_Target = donor_geometry->vertex[Marker_Donor][iVertex]->GetDonorPoint();

        Processor_Target = donor_geometry->vertex[Marker_Donor][iVertex]->GetDonorProcessor();

        GetDonor_Variable(donor_solution, donor_geometry, donor_config, Marker_Donor, iVertex, Point_Donor);

        for (iVar = 0; iVar < nVar; iVar++) 
          Buffer_Send_DonorVariables[iVertex*nVar+iVar] = Donor_Variable[iVar];


        Buffer_Send_DonorIndices[2*iVertex]     = Point_Target;
        Buffer_Send_DonorIndices[2*iVertex + 1] = Processor_Target;
      }

    }

    
#ifdef HAVE_MPI
    /*--- Once all the messages have been sent, we gather them all into the MASTER_NODE ---*/
    SU2_MPI::Gather(Buffer_Send_DonorVariables, nBuffer_DonorVariables, MPI_DOUBLE, Buffer_Recv_DonorVariables, nBuffer_DonorVariables, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_DonorIndices, nBuffer_DonorIndices, MPI_LONG, Buffer_Recv_DonorIndices, nBuffer_DonorIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);
    
#else
    for (unsigned long iVariable = 0; iVariable < nBuffer_DonorVariables; iVariable++)
      Buffer_Recv_DonorVariables[iVariable] = Buffer_Send_DonorVariables[iVariable];
    for (unsigned long iVariable = 0; iVariable < nBuffer_DonorIndices; iVariable++)
      Buffer_Recv_DonorIndices[iVariable] = Buffer_Send_DonorIndices[iVariable];
#endif
    
    /*--- Counter to determine where in the array we have to set the information ---*/
    long *Counter_Processor_Target = NULL;
    long iProcessor_Donor = 0, iIndex_Donor = 0;
    long iProcessor_Target = 0, iPoint_Target = 0, iIndex_Target = 0;
    long Point_Target_Send = 0, Processor_Target_Send = 0;
    
    /*--- Now we pack the information to send it over to the different processors ---*/
    
    if (rank == MASTER_NODE) {
      
      /*--- We set the counter to 0 ---*/
      Counter_Processor_Target = new long[nProcessor];
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        Counter_Processor_Target[iProcessor] = 0;
      }
      
      /*--- First we initialize the index vector to -1 ---*/
      /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
      for (iVertex = 0; iVertex < nProcessor*nBuffer_TargetIndices; iVertex++)
        Buffer_Send_TargetIndices[iVertex] = -2;
      
      /*--- As of now we do the loop over the flow points ---*/
      /*--- The number of points for flow and structure does not necessarily have to match ---*/
      /*--- In fact, it's possible that a processor asks for nStruct nodes and there are only ---*/
      /*--- nFlow < nStruct available; this is due to halo nodes ---*/
      
      /*--- For every processor from which we have received information ---*/
      /*--- (This is, for every processor on the structural side) ---*/
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        
        /*--- This is the initial index on the coordinates buffer for that particular processor on the structural side ---*/
        iProcessor_Donor = iProcessor*nBuffer_DonorVariables;
        /*--- This is the initial index on the donor index/processor buffer for that particular processor on the structural side ---*/
        iIndex_Donor = iProcessor*nBuffer_DonorIndices;
        
        /*--- For every vertex in the information retreived from iProcessor ---*/
        for (iVertex = 0; iVertex < Buffer_Recv_nVertexDonor[iProcessor]; iVertex++) {
          
          /*--- The processor and index for the flow are: ---*/
          Processor_Target_Send = Buffer_Recv_DonorIndices[iIndex_Donor+iVertex*2+1];
          Point_Target_Send     = Buffer_Recv_DonorIndices[iIndex_Donor+iVertex*2];
          
          /*--- Load the buffer at the appropriate position ---*/
          /*--- This is determined on the fluid side by:
           *--- Processor_Target*nBuffer_StructTraction -> Initial position of the processor array (fluid side)
           *--- +
           *--- Counter_Processor_Struct*nVar -> Initial position of the nVar array for the particular point on the fluid side
           *--- +
           *--- iVar -> Position within the nVar array that corresponds to a point
           *---
           *--- While on the structural side is:
           *--- iProcessor*nBuffer_FlowTraction -> Initial position on the processor array (structural side)
           *--- +
           *--- iVertex*nVar -> Initial position of the nVar array for the particular point on the structural side
           */
          
          /*--- We check that we are not setting the value for a halo node ---*/
          if (Point_Target_Send != -1) {
            iProcessor_Target = Processor_Target_Send*nBuffer_TargetVariables;
            iIndex_Target = Processor_Target_Send*nBuffer_TargetIndices;
            iPoint_Target = Counter_Processor_Target[Processor_Target_Send]*nVar;
            
            for (iVar = 0; iVar < nVar; iVar++)
              Buffer_Send_TargetVariables[iProcessor_Target + iPoint_Target + iVar] = Buffer_Recv_DonorVariables[iProcessor_Donor + iVertex*nVar + iVar];
            
            /*--- We set the fluid index at an appropriate position matching the coordinates ---*/
            Buffer_Send_TargetIndices[iIndex_Target + Counter_Processor_Target[Processor_Target_Send]] = Point_Target_Send;
            
            Counter_Processor_Target[Processor_Target_Send]++;
          }
          
        }
        
      }
      
    }
    
#ifdef HAVE_MPI
    /*--- Once all the messages have been prepared, we scatter them all from the MASTER_NODE ---*/
    SU2_MPI::Scatter(Buffer_Send_TargetVariables, nBuffer_TargetVariables, MPI_DOUBLE, Buffer_Recv_TargetVariables, nBuffer_TargetVariables, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Scatter(Buffer_Send_TargetIndices, nBuffer_TargetIndices, MPI_LONG, Buffer_Recv_TargetIndices, nBuffer_TargetIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
    for (unsigned long iVariable = 0; iVariable < nBuffer_TargetVariables; iVariable++)
      Buffer_Recv_TargetVariables[iVariable] = Buffer_Send_TargetVariables[iVariable];
    for (unsigned long iVariable = 0; iVariable < nBuffer_TargetIndices; iVariable++)
      Buffer_Recv_TargetIndices[iVariable] = Buffer_Send_TargetIndices[iVariable];
#endif
    
    long indexPoint_iVertex, Point_Target_Check =0;
    
    /*--- For the target marker we are studying ---*/
    if (Marker_Target >= 0) {
      
      /*--- We have identified the local index of the Structural marker ---*/
      /*--- We loop over all the vertices in that marker and in that particular processor ---*/
      
      for (iVertex = 0; iVertex < nLocalVertexTarget; iVertex++) {
        
        Point_Target = target_geometry->vertex[Marker_Target][iVertex]->GetNode();
        
        if (target_geometry->node[Point_Target]->GetDomain()) {
          /*--- Find the index of the point Point_Struct in the buffer Buffer_Recv_SetIndex ---*/
          indexPoint_iVertex = std::distance(Buffer_Recv_TargetIndices, std::find(Buffer_Recv_TargetIndices, Buffer_Recv_TargetIndices + MaxLocalVertexTarget, Point_Target));
          
          Point_Target_Check = Buffer_Recv_TargetIndices[indexPoint_iVertex];
          
          if (Point_Target_Check < 0 && fsi) {
            cout << "WARNING: A nonphysical point is being considered for traction transfer." << endl;
            exit(EXIT_FAILURE);
          }
          
          for (iVar = 0; iVar < nVar; iVar++)
            Target_Variable[iVar] = Buffer_Recv_TargetVariables[indexPoint_iVertex*nVar+iVar];
          
          SetTarget_Variable(target_solution, target_geometry, target_config, Marker_Target, iVertex, Point_Target); 
    
        }
        
      }
      
    }
    
    delete [] Buffer_Send_DonorVariables;
    delete [] Buffer_Send_DonorIndices;
    delete [] Buffer_Recv_TargetVariables;
    delete [] Buffer_Recv_TargetIndices;
    
    if (rank == MASTER_NODE) {
      delete [] Buffer_Recv_nVertexDonor;
      delete [] Buffer_Recv_nVertexTarget;
      delete [] Buffer_Recv_DonorVariables;
      delete [] Buffer_Recv_DonorIndices;
      delete [] Buffer_Send_TargetVariables;
      delete [] Buffer_Send_TargetIndices;
      delete [] Counter_Processor_Target;
    }
    
  }
  
  #ifdef HAVE_MPI
  if (rank == MASTER_NODE) 
    delete [] Buffer_Recv_mark;
  #endif
  
}

void CTransfer::Broadcast_InterfaceData_Matching(CSolver *donor_solution, CSolver *target_solution,
                                                 CGeometry *donor_geometry, CGeometry *target_geometry,
                                                 CConfig *donor_config, CConfig *target_config) {
  
  unsigned short nMarkerInt, nMarkerDonor, nMarkerTarget;       // Number of markers on the interface, donor and target side
  unsigned short iMarkerInt, iMarkerDonor, iMarkerTarget;       // Variables for iteration over markers
  int Marker_Donor = -1, Marker_Target = -1;
  int Target_check, Donor_check;
  
  unsigned long iVertex;                                // Variables for iteration over vertices and nodes
  
  unsigned short iVar;
  
  GetPhysical_Constants(donor_solution, target_solution, donor_geometry, target_geometry,
                        donor_config, target_config);
  
  unsigned long Point_Donor_Global, Donor_Global_Index;
  unsigned long Point_Donor, Point_Target;
  
  bool fsi = donor_config->GetFSI_Simulation();
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int *Buffer_Recv_mark=NULL, iRank;
  
  if (rank == MASTER_NODE) 
    Buffer_Recv_mark = new int[size];
#endif
  
  unsigned long iLocalVertex      = 0;
  unsigned long nLocalVertexDonor = 0, nLocalVertexDonorOwned = 0;
  
  unsigned long MaxLocalVertexDonor = 0;
  unsigned long TotalVertexDonor    = 0;
  
  unsigned long nBuffer_DonorVariables = 0;
  unsigned long nBuffer_DonorIndices   = 0;
  
  unsigned long nBuffer_BcastVariables = 0, nBuffer_BcastIndices = 0;
  
  int nProcessor = 0;
  
  /*--- Number of markers on the FSI interface ---*/
  
  nMarkerInt     = ( donor_config->GetMarker_n_FSIinterface() ) / 2;
  nMarkerTarget  = target_geometry->GetnMarker();
  nMarkerDonor   = donor_geometry->GetnMarker();
  
  nProcessor = size;
  
  /*--- Outer loop over the markers on the FSI interface: compute one by one ---*/
  /*--- The tags are always an integer greater than 1: loop from 1 to nMarkerFSI ---*/
  
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {
    
    Marker_Donor  = -1;
    Marker_Target = -1;
    
    /*--- Initialize pointer buffers inside the loop, so we can delete for each marker. ---*/
    unsigned long Buffer_Send_nVertexDonor[1], *Buffer_Recv_nVertexDonor = NULL;
    
    for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++) {
      /*--- If the tag GetMarker_All_FSIinterface(iMarkerDonor) equals the index we are looping at ---*/
      if ( donor_config->GetMarker_All_FSIinterface(iMarkerDonor) == iMarkerInt ) {
        Marker_Donor = iMarkerDonor;
        /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
        break;
      }
    }
    
    /*--- On the target side we only have to identify the marker; then we'll loop over it and retrieve from the fluid points ---*/
    
    for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++) {
      /*--- If the tag GetMarker_All_FSIinterface(iMarkerFlow) equals the index we are looping at ---*/
      if ( target_config->GetMarker_All_FSIinterface(iMarkerTarget) == iMarkerInt ) {
        /*--- Store the identifier for the fluid marker ---*/
        Marker_Target = iMarkerTarget;
        /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
        break;
      }
    }
    
    #ifdef HAVE_MPI

    Donor_check  = -1;
    Target_check = -1;

    /*--- We gather a vector in MASTER_NODE that determines if the boundary is not on the processor because of the partition or because the zone does not include it  ---*/

    SU2_MPI::Gather(&Marker_Donor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE) {
      for (iRank = 0; iRank < nProcessor; iRank++) {
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Donor_check = Buffer_Recv_mark[iRank];
          break;
        }       
      }
    }

    SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    SU2_MPI::Gather(&Marker_Target, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE) {
      for (iRank = 0; iRank < nProcessor; iRank++) {
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Target_check = Buffer_Recv_mark[iRank];
          break;
        }   
      }
    }

    SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    #else
    Donor_check  = Marker_Donor;
    Target_check = Marker_Target;   
    #endif

    if(Target_check == -1 || Donor_check == -1) {
      continue;
    }

    nLocalVertexDonorOwned = 0;
    nLocalVertexDonor      = 0;

    if( Marker_Donor != -1 ) {
      nLocalVertexDonor = donor_geometry->GetnVertex(Marker_Donor);

      for (iVertex = 0; iVertex < nLocalVertexDonor; iVertex++) {
        Point_Donor = donor_geometry->vertex[Marker_Donor][iVertex]->GetNode();
        if (donor_geometry->node[Point_Donor]->GetDomain())
          nLocalVertexDonorOwned++;
        }
    }

    Buffer_Send_nVertexDonor[0] = nLocalVertexDonor;                               // Retrieve total number of vertices on Donor marker
    if (rank == MASTER_NODE) Buffer_Recv_nVertexDonor = new unsigned long[size];   // Allocate memory to receive how many vertices are on each rank on the structural side

#ifdef HAVE_MPI
    /*--- We receive MaxLocalVertexDonor as the maximum number of vertices in one single processor on the donor side---*/
    SU2_MPI::Allreduce(&nLocalVertexDonor, &MaxLocalVertexDonor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    /*--- We receive TotalVertexDonorOwned as the total (real) number of vertices in one single interface marker on the donor side ---*/
    SU2_MPI::Allreduce(&nLocalVertexDonorOwned, &TotalVertexDonor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    /*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the structural side ---*/
    SU2_MPI::Gather(&Buffer_Send_nVertexDonor, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexDonor, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
    MaxLocalVertexDonor         = nLocalVertexDonor;
    TotalVertexDonor            = nLocalVertexDonorOwned;
    Buffer_Recv_nVertexDonor[0] = Buffer_Send_nVertexDonor[0];
#endif
    
    /*--- We will be gathering the donor information into the master node ---*/
    nBuffer_DonorVariables = MaxLocalVertexDonor * nVar;
    nBuffer_DonorIndices = MaxLocalVertexDonor;
    
    /*--- Then we will broadcasting it to all the processors so they can retrieve the info they need ---*/
    /*--- We only broadcast those nodes that we need ---*/
    nBuffer_BcastVariables = TotalVertexDonor * nVar;
    nBuffer_BcastIndices = TotalVertexDonor;
    
    /*--- Send and Recv buffers ---*/
    
    /*--- Buffers to send and receive the variables in the donor mesh ---*/
    su2double *Buffer_Send_DonorVariables = new su2double[nBuffer_DonorVariables];
    su2double *Buffer_Recv_DonorVariables = NULL;
    
    /*--- Buffers to send and receive the indices in the donor mesh ---*/
    long *Buffer_Send_DonorIndices = new long[nBuffer_DonorIndices];
    long *Buffer_Recv_DonorIndices = NULL;
    
    /*--- Buffers to broadcast the variables and the indices ---*/
    su2double *Buffer_Bcast_Variables = new su2double[nBuffer_BcastVariables];
    long *Buffer_Bcast_Indices        = new long[nBuffer_BcastIndices];
    
    /*--- Prepare the receive buffers (1st step) and send buffers (2nd step) on the master node only. ---*/
    
    if (rank == MASTER_NODE) {
      Buffer_Recv_DonorVariables  = new su2double[size*nBuffer_DonorVariables];
      Buffer_Recv_DonorIndices    = new long[size*nBuffer_DonorIndices];
    }
    
    /*--- On the donor side ---*/
    /*--- First we initialize all of the indices and processors to -1 ---*/
    /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
    for (iVertex = 0; iVertex < nBuffer_DonorIndices; iVertex++)
      Buffer_Send_DonorIndices[iVertex] = -1;
    
    for (iVertex = 0; iVertex < nLocalVertexDonor; iVertex++) {

      Point_Donor = donor_geometry->vertex[Marker_Donor][iVertex]->GetNode();

      /*--- If this processor owns the node ---*/
      if (donor_geometry->node[Point_Donor]->GetDomain()) {
        
        GetDonor_Variable(donor_solution, donor_geometry, donor_config, Marker_Donor, iVertex, Point_Donor);

        for (iVar = 0; iVar < nVar; iVar++) 
          Buffer_Send_DonorVariables[iVertex*nVar+iVar] = Donor_Variable[iVar];
          
        Point_Donor_Global = donor_geometry->node[Point_Donor]->GetGlobalIndex();

        Buffer_Send_DonorIndices[iVertex] = Point_Donor_Global;
      }

    }
    
#ifdef HAVE_MPI
    /*--- Once all the messages have been prepared, we gather them all into the MASTER_NODE ---*/
    SU2_MPI::Gather(Buffer_Send_DonorVariables, nBuffer_DonorVariables, MPI_DOUBLE, Buffer_Recv_DonorVariables, nBuffer_DonorVariables, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_DonorIndices, nBuffer_DonorIndices, MPI_LONG, Buffer_Recv_DonorIndices, nBuffer_DonorIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);
    
#else
    for (unsigned long iVariable = 0; iVariable < nBuffer_DonorVariables; iVariable++)
      Buffer_Recv_DonorVariables[iVariable] = Buffer_Send_DonorVariables[iVariable];
    for (unsigned long iVariable = 0; iVariable < nBuffer_DonorIndices; iVariable++)
      Buffer_Recv_DonorIndices[iVariable] = Buffer_Send_DonorIndices[iVariable];
#endif
    
    /*--- Now we pack the information to send it over to the different processors ---*/
    
    if (rank == MASTER_NODE) {
      
      /*--- For all the data we have received ---*/
      /*--- We initialize a counter to determine the position in the broadcast vector ---*/
      iLocalVertex = 0;
      
      for (iVertex = 0; iVertex < nProcessor*nBuffer_DonorIndices; iVertex++) {
        
        /*--- If the donor index is not -1 (this is, if the node is not originally a halo node) ---*/
        if (Buffer_Recv_DonorIndices[iVertex] != -1) {
          
          /*--- We set the donor index ---*/
          Buffer_Bcast_Indices[iLocalVertex] = Buffer_Recv_DonorIndices[iVertex];
          
          for (iVar = 0; iVar < nVar; iVar++) {
            Buffer_Bcast_Variables[iLocalVertex*nVar+iVar] = Buffer_Recv_DonorVariables[iVertex*nVar + iVar];
          }
          
          iLocalVertex++;
          
        }
        
        if (iLocalVertex == TotalVertexDonor) break;
        
      }
      
    }
    
#ifdef HAVE_MPI
    SU2_MPI::Bcast(Buffer_Bcast_Variables, nBuffer_BcastVariables, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(Buffer_Bcast_Indices, nBuffer_BcastIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);
#endif
    
    long indexPoint_iVertex, Point_Target_Check=0;

    /*--- For the target marker we are studying ---*/
    if (Marker_Target >= 0) {
      
      /*--- We have identified the local index of the Structural marker ---*/
      /*--- We loop over all the vertices in that marker and in that particular processor ---*/
      
      for (iVertex = 0; iVertex < target_geometry->GetnVertex(Marker_Target); iVertex++) {
        
        Point_Target = target_geometry->vertex[Marker_Target][iVertex]->GetNode();
        
        /*--- If this processor owns the node ---*/
        if (target_geometry->node[Point_Target]->GetDomain()) {
          
          /*--- Find the global index of the donor point for Point_Target ---*/
          Donor_Global_Index = target_geometry->vertex[Marker_Target][iVertex]->GetGlobalDonorPoint();
          
          /*--- Find the index of the global donor point in the buffer Buffer_Bcast_Indices ---*/
          indexPoint_iVertex = std::distance(Buffer_Bcast_Indices, std::find(Buffer_Bcast_Indices, Buffer_Bcast_Indices + nBuffer_BcastIndices, Donor_Global_Index));
          
          Point_Target_Check = Buffer_Bcast_Indices[indexPoint_iVertex];
          
          if (Point_Target_Check < 0 && fsi) {
            cout << "WARNING: A nonphysical point is being considered for traction transfer." << endl;
            exit(EXIT_FAILURE);
          }
          
          for (iVar = 0; iVar < nVar; iVar++)
            Target_Variable[iVar] = Buffer_Bcast_Variables[indexPoint_iVertex*nVar+iVar];
          
          if (Point_Target_Check >= 0)
            SetTarget_Variable(target_solution, target_geometry, target_config, Marker_Target, iVertex, Point_Target);
          
        }
        
      }
      
    }
    
    delete [] Buffer_Send_DonorVariables;
    delete [] Buffer_Send_DonorIndices;
    delete [] Buffer_Bcast_Variables;
    delete [] Buffer_Bcast_Indices;
    
    if (rank == MASTER_NODE) {
      delete [] Buffer_Recv_nVertexDonor;
      delete [] Buffer_Recv_DonorVariables;
      delete [] Buffer_Recv_DonorIndices;
    }  
  }
  
  #ifdef HAVE_MPI
  if (rank == MASTER_NODE) 
    delete [] Buffer_Recv_mark;
  #endif
  
}

void CTransfer::Broadcast_InterfaceData_Interpolate(CSolver *donor_solution, CSolver *target_solution,
                                                    CGeometry *donor_geometry, CGeometry *target_geometry,
                                                    CConfig *donor_config, CConfig *target_config) {
  
  
  unsigned short nMarkerInt, nMarkerDonor, nMarkerTarget;       // Number of markers on the interface, donor and target side
  unsigned short iMarkerInt, iMarkerDonor, iMarkerTarget;       // Variables for iteration over markers
  int Marker_Donor, Marker_Target;
  int Target_check, Donor_check;
  
  unsigned long iVertex;                            // Variables for iteration over vertices and nodes
  
  unsigned short iVar;
  
  GetPhysical_Constants(donor_solution, target_solution, donor_geometry, target_geometry,
                        donor_config, target_config);
  
  unsigned long Point_Donor_Global, Donor_Global_Index;
  unsigned long Point_Donor, Point_Target;
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
  bool fsi = donor_config->GetFSI_Simulation();
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int *Buffer_Recv_mark=NULL, iRank;
  
  if (rank == MASTER_NODE) 
    Buffer_Recv_mark = new int[size];
#endif
  
  unsigned long Buffer_Send_nVertexDonor[1], *Buffer_Recv_nVertexDonor;
  unsigned long iLocalVertex = 0;
  unsigned long nLocalVertexDonor = 0, nLocalVertexDonorOwned = 0;
  
  unsigned long MaxLocalVertexDonor = 0, TotalVertexDonor = 0;
  
  unsigned long nBuffer_DonorVariables = 0;
  unsigned long nBuffer_DonorIndices = 0;
  
  unsigned long nBuffer_BcastVariables = 0, nBuffer_BcastIndices = 0;
  
  int nProcessor = 0;
  
  /*--- Number of markers on the FSI interface ---*/
  
  nMarkerInt     = (donor_config->GetMarker_n_FSIinterface())/2;
  nMarkerTarget  = target_config->GetnMarker_All();
  nMarkerDonor   = donor_config->GetnMarker_All();
  
  nProcessor = size;
  
  /*--- Outer loop over the markers on the FSI interface: compute one by one ---*/
  /*--- The tags are always an integer greater than 1: loop from 1 to nMarkerFSI ---*/
  
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {
    
    Buffer_Recv_nVertexDonor = NULL;
    
    Marker_Donor = -1;
    Marker_Target = -1;
    
    /*--- The donor and target markers are tagged with the same index.
     *--- This is independent of the MPI domain decomposition.
     *--- We need to loop over all markers on both sides and get the number of nodes
     *--- that belong to each FSI marker for each processor ---*/
    
    /*--- On the donor side ---*/
    
    for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++) {
      /*--- If the tag GetMarker_All_FSIinterface(iMarkerDonor) equals the index we are looping at ---*/
      if ( donor_config->GetMarker_All_FSIinterface(iMarkerDonor) == iMarkerInt ) {
        /*--- Store the identifier for the structural marker ---*/
        Marker_Donor = iMarkerDonor;
        /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
        break;
      }
    }
    
    /*--- On the target side we only have to identify the marker; then we'll loop over it and retrieve from the donor points ---*/
    
    for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++) {
      /*--- If the tag GetMarker_All_FSIinterface(iMarkerFlow) equals the index we are looping at ---*/
      if ( target_config->GetMarker_All_FSIinterface(iMarkerTarget) == iMarkerInt ) {
        /*--- Store the identifier for the fluid marker ---*/
        Marker_Target = iMarkerTarget;
        /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
        break;
      }
    }
    
    #ifdef HAVE_MPI

    Donor_check  = -1;
    Target_check = -1;

    /*--- We gather a vector in MASTER_NODE that determines if the boundary is not on the processor because of the partition or because the zone does not include it  ---*/

    SU2_MPI::Gather(&Marker_Donor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE) {
      for (iRank = 0; iRank < nProcessor; iRank++) {
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Donor_check = Buffer_Recv_mark[iRank];
          break;
        }       
      }
    }

    SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    SU2_MPI::Gather(&Marker_Target, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE) {
      for (iRank = 0; iRank < nProcessor; iRank++) {
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Target_check = Buffer_Recv_mark[iRank];
          break;
        }   
      }
    }

    SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    #else
    Donor_check  = Marker_Donor;
    Target_check = Marker_Target;   
    #endif

    if(Target_check == -1 || Donor_check == -1) {
      continue;
    }

    nLocalVertexDonorOwned = 0;
    nLocalVertexDonor      = 0;
    
    if( Marker_Donor != -1 ) {
        nLocalVertexDonor = donor_geometry->GetnVertex(Marker_Donor);
    
        for (iVertex = 0; iVertex < nLocalVertexDonor; iVertex++) {
              Point_Donor = donor_geometry->vertex[Marker_Donor][iVertex]->GetNode();
              if (donor_geometry->node[Point_Donor]->GetDomain())
                nLocalVertexDonorOwned++;
            }
    }
    
    Buffer_Send_nVertexDonor[0] = nLocalVertexDonor;                   // Retrieve total number of vertices on Donor marker
    if (rank == MASTER_NODE) Buffer_Recv_nVertexDonor = new unsigned long[size];   // Allocate memory to receive how many vertices are on each rank on the structural side
    
#ifdef HAVE_MPI
    /*--- We receive MaxLocalVertexDonor as the maximum number of vertices in one single processor on the donor side---*/
    SU2_MPI::Allreduce(&nLocalVertexDonor, &MaxLocalVertexDonor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    /*--- We receive TotalVertexDonorOwned as the total (real) number of vertices in one single interface marker on the donor side ---*/
    SU2_MPI::Allreduce(&nLocalVertexDonorOwned, &TotalVertexDonor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    /*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the structural side ---*/
    SU2_MPI::Gather(&Buffer_Send_nVertexDonor, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexDonor, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
    MaxLocalVertexDonor         = nLocalVertexDonor;
    TotalVertexDonor            = nLocalVertexDonorOwned;
    Buffer_Recv_nVertexDonor[0] = Buffer_Send_nVertexDonor[0];
#endif
    
    /*--- We will be gathering the donor information into the master node ---*/
    nBuffer_DonorVariables = MaxLocalVertexDonor * nVar;
    nBuffer_DonorIndices   = MaxLocalVertexDonor;
    
    /*--- Then we will broadcasting it to all the processors so they can retrieve the info they need ---*/
    /*--- We only broadcast those nodes that we need ---*/
    nBuffer_BcastVariables = TotalVertexDonor * nVar;
    nBuffer_BcastIndices   = TotalVertexDonor;
    
    /*--- Send and Recv buffers ---*/
    
    /*--- Buffers to send and receive the variables in the donor mesh ---*/
    su2double *Buffer_Send_DonorVariables = new su2double[nBuffer_DonorVariables];
    su2double *Buffer_Recv_DonorVariables = NULL;
    
    /*--- Buffers to send and receive the indices in the donor mesh ---*/
    long *Buffer_Send_DonorIndices = new long[nBuffer_DonorIndices];
    long *Buffer_Recv_DonorIndices = NULL;
    
    /*--- Buffers to broadcast the variables and the indices ---*/
    su2double *Buffer_Bcast_Variables = new su2double[nBuffer_BcastVariables];
    long *Buffer_Bcast_Indices        = new long[nBuffer_BcastIndices];
    
    /*--- Prepare the receive buffers (1st step) and send buffers (2nd step) on the master node only. ---*/
    
    if (rank == MASTER_NODE) {
      Buffer_Recv_DonorVariables  = new su2double[size*nBuffer_DonorVariables];
      Buffer_Recv_DonorIndices    = new long[size*nBuffer_DonorIndices];
    }
    
    /*--- On the donor side ---*/
    /*--- First we initialize all of the indices and processors to -1 ---*/
    /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
    for (iVertex = 0; iVertex < nBuffer_DonorIndices; iVertex++)
      Buffer_Send_DonorIndices[iVertex] = -1;
    

    for (iVertex = 0; iVertex < nLocalVertexDonor; iVertex++) {
    Point_Donor = donor_geometry->vertex[Marker_Donor][iVertex]->GetNode();
    
    /*--- If this processor owns the node ---*/
        
    if (donor_geometry->node[Point_Donor]->GetDomain()) {
    
        GetDonor_Variable(donor_solution, donor_geometry, donor_config, Marker_Donor, iVertex, Point_Donor);
    
        for (iVar = 0; iVar < nVar; iVar++) 
            Buffer_Send_DonorVariables[iVertex*nVar+iVar] = Donor_Variable[iVar];
        
        Point_Donor_Global = donor_geometry->node[Point_Donor]->GetGlobalIndex();
        Buffer_Send_DonorIndices[iVertex]     = Point_Donor_Global;
    }

    }
      
#ifdef HAVE_MPI
    /*--- Once all the messages have been prepared, we gather them all into the MASTER_NODE ---*/
    SU2_MPI::Gather(Buffer_Send_DonorVariables, nBuffer_DonorVariables, MPI_DOUBLE, Buffer_Recv_DonorVariables, nBuffer_DonorVariables, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_DonorIndices, nBuffer_DonorIndices, MPI_LONG, Buffer_Recv_DonorIndices, nBuffer_DonorIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);
    
#else
    for (unsigned long iVariable = 0; iVariable < nBuffer_DonorVariables; iVariable++)
      Buffer_Recv_DonorVariables[iVariable] = Buffer_Send_DonorVariables[iVariable];
    for (unsigned long iVariable = 0; iVariable < nBuffer_DonorIndices; iVariable++)
      Buffer_Recv_DonorIndices[iVariable] = Buffer_Send_DonorIndices[iVariable];
#endif
    
    /*--- Now we pack the information to send it over to the different processors ---*/
    
    if (rank == MASTER_NODE) {

        /*--- For all the data we have received ---*/
        /*--- We initialize a counter to determine the position in the broadcast vector ---*/
        iLocalVertex = 0;
      
        for (iVertex = 0; iVertex < nProcessor*nBuffer_DonorIndices; iVertex++) {
        
            /*--- If the donor index is not -1 (this is, if the node is not originally a halo node) ---*/
            if (Buffer_Recv_DonorIndices[iVertex] != -1) {
          
                /*--- We set the donor index ---*/
                Buffer_Bcast_Indices[iLocalVertex] = Buffer_Recv_DonorIndices[iVertex];
          
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Bcast_Variables[iLocalVertex*nVar+iVar] = Buffer_Recv_DonorVariables[iVertex*nVar + iVar];
                
          
                iLocalVertex++;
          
            }
        
            if (iLocalVertex == TotalVertexDonor) break;
        
        }
      
    }
    
#ifdef HAVE_MPI
    SU2_MPI::Bcast(Buffer_Bcast_Variables, nBuffer_BcastVariables, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(Buffer_Bcast_Indices, nBuffer_BcastIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);
#endif

    long indexPoint_iVertex, Point_Target_Check=0;
    unsigned short iDonorPoint, nDonorPoints;
    su2double donorCoeff;
    
    /*--- For the target marker we are studying ---*/
    if (Marker_Target >= 0) {
      
      /*--- We have identified the local index of the Structural marker ---*/
      /*--- We loop over all the vertices in that marker and in that particular processor ---*/
      
      for (iVertex = 0; iVertex < target_geometry->GetnVertex(Marker_Target); iVertex++) {
        
        Point_Target = target_geometry->vertex[Marker_Target][iVertex]->GetNode();

        /*--- If this processor owns the node ---*/
        if (target_geometry->node[Point_Target]->GetDomain()) {
          TotalVertexDonor++;
          nDonorPoints = target_geometry->vertex[Marker_Target][iVertex]->GetnDonorPoints();
          
          /*--- As we will be adding data, we need to set the variable to 0 ---*/
          for (iVar = 0; iVar < nVar; iVar++) Target_Variable[iVar] = 0.0;
          
          /*--- For the number of donor points ---*/
          for (iDonorPoint = 0; iDonorPoint < nDonorPoints; iDonorPoint++) {
            
            /*--- Find the global index of the donor points for Point_Target ---*/
            Donor_Global_Index = target_geometry->vertex[Marker_Target][iVertex]->GetInterpDonorPoint(iDonorPoint);
            
            /*--- We need to get the donor coefficient in a way like this: ---*/
            donorCoeff = target_geometry->vertex[Marker_Target][iVertex]->GetDonorCoeff(iDonorPoint);
            
            /*--- Find the index of the global donor point in the buffer Buffer_Bcast_Indices ---*/
           
            indexPoint_iVertex = std::distance(Buffer_Bcast_Indices, std::find(Buffer_Bcast_Indices, Buffer_Bcast_Indices + nBuffer_BcastIndices, Donor_Global_Index));
           
            Point_Target_Check = Buffer_Bcast_Indices[indexPoint_iVertex];

            if (Point_Target_Check < 0 && fsi) {
              cout << "WARNING: A nonphysical point is being considered for traction transfer." << endl;
              exit(EXIT_FAILURE);
            }

            for (iVar = 0; iVar < nVar; iVar++)
              Target_Variable[iVar] += donorCoeff * Buffer_Bcast_Variables[indexPoint_iVertex*nVar+iVar];
          }

          if (Point_Target_Check >= 0)
            SetTarget_Variable(target_solution, target_geometry, target_config, Marker_Target, iVertex, Point_Target);
        }
        
      }
      
    }

    delete [] Buffer_Send_DonorVariables;
    delete [] Buffer_Send_DonorIndices;
    delete [] Buffer_Bcast_Variables;
    delete [] Buffer_Bcast_Indices;
    
    if (rank == MASTER_NODE) {
      delete [] Buffer_Recv_nVertexDonor;
      delete [] Buffer_Recv_DonorVariables;
      delete [] Buffer_Recv_DonorIndices;
    }
      
  }
  
  #ifdef HAVE_MPI
  if (rank == MASTER_NODE) 
    delete [] Buffer_Recv_mark;
  #endif
}

void CTransfer::Allgather_InterfaceData(CSolver *donor_solution, CSolver *target_solution,
                                        CGeometry *donor_geometry, CGeometry *target_geometry,
                                        CConfig *donor_config, CConfig *target_config) {
  
  
  unsigned short nMarkerInt, nMarkerDonor, nMarkerTarget;       // Number of markers on the interface, donor and target side
  unsigned short iMarkerInt, iMarkerDonor, iMarkerTarget;       // Variables for iteration over markers
  int Marker_Donor = -1, Marker_Target = -1;
  int Target_check, Donor_check;
  
  unsigned long iVertex;                                // Variables for iteration over vertices and nodes
  unsigned short iVar;
  
  GetPhysical_Constants(donor_solution, target_solution, donor_geometry, target_geometry,
                        donor_config, target_config);
  
  unsigned long Point_Donor_Global, Donor_Global_Index;
  unsigned long Point_Donor, Point_Target;
  
  bool fsi = donor_config->GetFSI_Simulation();
  
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  int rank = MASTER_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int *Buffer_Recv_mark=NULL, iRank;
  
  if (rank == MASTER_NODE) 
    Buffer_Recv_mark = new int[size];
#endif
  
  unsigned long iLocalVertex = 0;
  unsigned long nLocalVertexDonor = 0;
  
  unsigned long MaxLocalVertexDonor = 0;
  
  unsigned long nBuffer_DonorVariables = 0;
  unsigned long nBuffer_DonorIndices = 0;
  
  int nProcessor = 0;
  
  /*--- Number of markers on the FSI interface ---*/
  
  nMarkerInt     = (donor_config->GetMarker_n_FSIinterface())/2;
  nMarkerTarget  = target_geometry->GetnMarker();
  nMarkerDonor   = donor_geometry->GetnMarker();
  
  nProcessor = size;
  
  /*--- Outer loop over the markers on the FSI interface: compute one by one ---*/
  /*--- The tags are always an integer greater than 1: loop from 1 to nMarkerFSI ---*/
  
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {
    
    Marker_Donor = -1;
    Marker_Target = -1;
    
    /*--- Initialize pointer buffers inside the loop, so we can delete for each marker. ---*/
    /*--- We are only sending the values the processor owns ---*/
    unsigned long Buffer_Send_nVertexDonor[1], *Buffer_Recv_nVertexDonor = NULL;
    
    /*--- The donor and target markers are tagged with the same index.
     *--- This is independent of the MPI domain decomposition.
     *--- We need to loop over all markers on both sides and get the number of nodes
     *--- that belong to each FSI marker for each processor ---*/
    
    /*--- On the donor side ---*/
    
    for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++) {
      /*--- If the tag GetMarker_All_FSIinterface(iMarkerDonor) equals the index we are looping at ---*/
      if ( donor_config->GetMarker_All_FSIinterface(iMarkerDonor) == iMarkerInt ) {
        /*--- Store the identifier for the structural marker ---*/
        Marker_Donor = iMarkerDonor;
        /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
        break;
      }
    }
    
    /*--- On the target side we only have to identify the marker; then we'll loop over it and retrieve from the donor points ---*/
    
    for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++) {
      /*--- If the tag GetMarker_All_FSIinterface(iMarkerFlow) equals the index we are looping at ---*/
      if ( target_config->GetMarker_All_FSIinterface(iMarkerTarget) == iMarkerInt ) {
        /*--- Store the identifier for the fluid marker ---*/
        Marker_Target = iMarkerTarget;
        /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
        break;
      }
    }
    
    #ifdef HAVE_MPI

    Donor_check  = -1;
    Target_check = -1;

    /*--- We gather a vector in MASTER_NODE that determines if the boundary is not on the processor because of the partition or because the zone does not include it  ---*/

    SU2_MPI::Gather(&Marker_Donor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE) {
      for (iRank = 0; iRank < nProcessor; iRank++) {
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Donor_check = Buffer_Recv_mark[iRank];
          break;
        }       
      }
    }

    SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    SU2_MPI::Gather(&Marker_Target, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE) {
      for (iRank = 0; iRank < nProcessor; iRank++) {
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Target_check = Buffer_Recv_mark[iRank];
          break;
        }   
      }
    }

    SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    #else
    Donor_check  = Marker_Donor;
    Target_check = Marker_Target;   
    #endif

    if(Target_check == -1 || Donor_check == -1) {
      continue;
    }

    nLocalVertexDonor  = 0;

    if( Marker_Donor != -1 )
        nLocalVertexDonor = donor_geometry->GetnVertex(Marker_Donor);
    
    Buffer_Send_nVertexDonor[0] = nLocalVertexDonor;      // Retrieve total number of vertices on Donor marker
    Buffer_Recv_nVertexDonor = new unsigned long[size];   // Allocate memory to receive how many vertices are on each rank on the structural side
    
#ifdef HAVE_MPI
    /*--- We receive MaxLocalVertexDonor as the maximum number of vertices in one single processor on the donor side---*/
    SU2_MPI::Allreduce(&nLocalVertexDonor, &MaxLocalVertexDonor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    /*--- We gather a vector in all processors that determines how many elements are there on each processor on the structural side ---*/
    SU2_MPI::Allgather(&Buffer_Send_nVertexDonor, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexDonor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
    MaxLocalVertexDonor         = nLocalVertexDonor;
    Buffer_Recv_nVertexDonor[0] = Buffer_Send_nVertexDonor[0];
#endif
    
    /*--- We will be gathering the donor information into the master node ---*/
    nBuffer_DonorVariables = MaxLocalVertexDonor * nVar;
    nBuffer_DonorIndices = MaxLocalVertexDonor;
    
    /*--- Send and Recv buffers ---*/
    
    /*--- Buffers to send and receive the variables in the donor mesh ---*/
    su2double *Buffer_Send_DonorVariables = new su2double[nBuffer_DonorVariables];
    su2double *Buffer_Recv_DonorVariables = new su2double[size*nBuffer_DonorVariables];
    
    /*--- Buffers to send and receive the indices in the donor mesh ---*/
    long *Buffer_Send_DonorIndices = new long[nBuffer_DonorIndices];
    long *Buffer_Recv_DonorIndices = new long[size*nBuffer_DonorIndices];
    
    /*--- On the donor side ---*/
    /*--- First we initialize all of the indices and processors to -1 ---*/
    /*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
    for (iVertex = 0; iVertex < nBuffer_DonorIndices; iVertex++)
      Buffer_Send_DonorIndices[iVertex] = -1;
    
    /*--- Also to avoid having random values in the variables vector ---*/
    for (iVertex = 0; iVertex < nBuffer_DonorIndices; iVertex++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Buffer_Send_DonorVariables[iVertex*nVar + iVar] = 0.0;
      }
    }
    
    if (Marker_Donor >= 0) {
      
      iLocalVertex = 0;
      
      for (iVertex = 0; iVertex < donor_geometry->GetnVertex(Marker_Donor); iVertex++) {
        
        Point_Donor = donor_geometry->vertex[Marker_Donor][iVertex]->GetNode();
        
        GetDonor_Variable(donor_solution, donor_geometry, donor_config, Marker_Donor, iVertex, Point_Donor);
        
        /*--- If this processor owns the node ---*/
        if (donor_geometry->node[Point_Donor]->GetDomain()) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Buffer_Send_DonorVariables[iLocalVertex*nVar+iVar] = Donor_Variable[iVar];
          }
          
          Point_Donor_Global = donor_geometry->node[Point_Donor]->GetGlobalIndex();
          Buffer_Send_DonorIndices[iLocalVertex] = Point_Donor_Global;
          
          iLocalVertex++;
        }
        
      }
      
    }
    
#ifdef HAVE_MPI
    /*--- Once all the messages have been prepared, we gather them all into all the processors ---*/
    SU2_MPI::Allgather(Buffer_Send_DonorVariables, nBuffer_DonorVariables, MPI_DOUBLE, Buffer_Recv_DonorVariables, nBuffer_DonorVariables, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_DonorIndices, nBuffer_DonorIndices, MPI_LONG, Buffer_Recv_DonorIndices, nBuffer_DonorIndices, MPI_LONG, MPI_COMM_WORLD);
#else
    for (unsigned long iVariable = 0; iVariable < nBuffer_DonorVariables; iVariable++)
      Buffer_Recv_DonorVariables[iVariable] = Buffer_Send_DonorVariables[iVariable];
    for (unsigned long iVariable = 0; iVariable < nBuffer_DonorIndices; iVariable++)
      Buffer_Recv_DonorIndices[iVariable] = Buffer_Send_DonorIndices[iVariable];
#endif
    
    long indexPoint_iVertex, Point_Target_Check = 0;
    unsigned short iDonorPoint, nDonorPoints;
    su2double donorCoeff;
    
    /*--- For the target marker we are studying ---*/
    if (Marker_Target >= 0) {
      
      /*--- We have identified the local index of the Structural marker ---*/
      /*--- We loop over all the vertices in that marker and in that particular processor ---*/
      
      for (iVertex = 0; iVertex < target_geometry->GetnVertex(Marker_Target); iVertex++) {
        
        Point_Target = target_geometry->vertex[Marker_Target][iVertex]->GetNode();
        
        /*--- If this processor owns the node ---*/
        if (target_geometry->node[Point_Target]->GetDomain()) {
          
          nDonorPoints = target_geometry->vertex[Marker_Target][iVertex]->GetnDonorPoints();
          
          /*--- As we will be adding data, we need to set the variable to 0 ---*/
          for (iVar = 0; iVar < nVar; iVar++) Target_Variable[iVar] = 0.0;
          
          /*--- For the number of donor points ---*/
          for (iDonorPoint = 0; iDonorPoint < nDonorPoints; iDonorPoint++) {
            
            /*--- Find the global index of the donor points for Point_Target ---*/
            Donor_Global_Index = target_geometry->vertex[Marker_Target][iVertex]->GetInterpDonorPoint(iDonorPoint);
            
            /*--- We need to get the donor coefficient in a way like this: ---*/
            donorCoeff = target_geometry->vertex[Marker_Target][iVertex]->GetDonorCoeff(iDonorPoint);
            
            /*--- Find the index of the global donor point in the buffer Buffer_Bcast_Indices ---*/
            indexPoint_iVertex = std::distance(Buffer_Recv_DonorIndices, std::find(Buffer_Recv_DonorIndices, Buffer_Recv_DonorIndices + nProcessor*nBuffer_DonorIndices, Donor_Global_Index));
            
            Point_Target_Check = Buffer_Recv_DonorIndices[indexPoint_iVertex];
            
            if (Point_Target_Check < 0 && fsi) {
              cout << "WARNING: A nonphysical point is being considered for traction transfer." << endl;
              exit(EXIT_FAILURE);
            }
            
            for (iVar = 0; iVar < nVar; iVar++)
              Target_Variable[iVar] += donorCoeff * Buffer_Recv_DonorVariables[indexPoint_iVertex*nVar+iVar];
          }
          
          if (Point_Target_Check >= 0)
            SetTarget_Variable(target_solution, target_geometry, target_config, Marker_Target, iVertex, Point_Target);          
          
        }
        
      }
      
    }
    
    delete [] Buffer_Send_DonorVariables;
    delete [] Buffer_Send_DonorIndices;
    
    delete [] Buffer_Recv_DonorVariables;
    delete [] Buffer_Recv_DonorIndices;
    
    delete [] Buffer_Recv_nVertexDonor;
    
  }

  #ifdef HAVE_MPI
  if (rank == MASTER_NODE) 
    delete [] Buffer_Recv_mark;
  #endif
  
}
