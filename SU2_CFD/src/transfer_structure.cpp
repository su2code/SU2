/*!
 * \file transfer_structure.cpp
 * \brief Main subroutines for MPI transfer of information between zones
 * \author R. Sanchez
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  Physical_Constants = NULL;
  Donor_Variable     = NULL;
  Target_Variable    = NULL;
  SpanLevelDonor     = NULL;
  SpanValueCoeffTarget = NULL;
  
  nVar = 0;
  
}

CTransfer::CTransfer(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) {
  
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  Physical_Constants = NULL;
  Donor_Variable     = NULL;
  Target_Variable    = NULL;

  unsigned short iVar;
  
  Physical_Constants = new su2double[val_nConst];
  Donor_Variable     = new su2double[val_nVar];
  Target_Variable    = new su2double[val_nVar];

  /*--- By default, the value is aggregated in the transfer routine ---*/
  valAggregated      = true;
    
  nVar = val_nVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Donor_Variable[iVar]  = 0.0;
    Target_Variable[iVar] = 0.0;
  }
  
  for (iVar = 0; iVar < val_nConst; iVar++) {
    Physical_Constants[iVar] = 0.0;
  }

  SpanLevelDonor       = NULL;
  SpanValueCoeffTarget = NULL;

}

CTransfer::~CTransfer(void) {
  
  if (Physical_Constants   != NULL) delete [] Physical_Constants;
  if (Donor_Variable       != NULL) delete [] Donor_Variable;
  if (Target_Variable      != NULL) delete [] Target_Variable;

  if (SpanValueCoeffTarget != NULL) delete[] SpanValueCoeffTarget;
  if (SpanLevelDonor       != NULL) delete[] SpanLevelDonor;

  
}

void CTransfer::Broadcast_InterfaceData(CSolver *donor_solution, CSolver *target_solution,
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
  
#ifdef HAVE_MPI
  int *Buffer_Recv_mark = NULL, iRank;

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
  
  nMarkerInt     = (donor_config->GetMarker_n_ZoneInterface())/2;
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
      /*--- If the tag GetMarker_All_ZoneInterface(iMarkerDonor) equals the index we are looping at ---*/
      if ( donor_config->GetMarker_All_ZoneInterface(iMarkerDonor) == iMarkerInt ) {
        /*--- Store the identifier for the structural marker ---*/
        Marker_Donor = iMarkerDonor;
        /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
        break;
      }
    }
    
    /*--- On the target side we only have to identify the marker; then we'll loop over it and retrieve from the donor points ---*/
    
    for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++) {
      /*--- If the tag GetMarker_All_ZoneInterface(iMarkerFlow) equals the index we are looping at ---*/
      if ( target_config->GetMarker_All_ZoneInterface(iMarkerTarget) == iMarkerInt ) {
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

    long indexPoint_iVertex;
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
          
          InitializeTarget_Variable(target_solution, Marker_Target, iVertex, nDonorPoints);
          
          /*--- For the number of donor points ---*/
          for (iDonorPoint = 0; iDonorPoint < nDonorPoints; iDonorPoint++) {
            
            /*--- Find the global index of the donor points for Point_Target ---*/
            Donor_Global_Index = target_geometry->vertex[Marker_Target][iVertex]->GetInterpDonorPoint(iDonorPoint);
            
            /*--- We need to get the donor coefficient in a way like this: ---*/
            donorCoeff = target_geometry->vertex[Marker_Target][iVertex]->GetDonorCoeff(iDonorPoint);
            
            /*--- Find the index of the global donor point in the buffer Buffer_Bcast_Indices ---*/
           
            indexPoint_iVertex = std::distance(Buffer_Bcast_Indices, std::find(Buffer_Bcast_Indices, Buffer_Bcast_Indices + nBuffer_BcastIndices, Donor_Global_Index));
           
            /*--- Recover the Target_Variable from the buffer of variables ---*/
            RecoverTarget_Variable(indexPoint_iVertex, Buffer_Bcast_Variables, donorCoeff);

            /*--- If the value is not directly aggregated in the previous function ---*/
            if (!valAggregated) SetTarget_Variable(target_solution, target_geometry, target_config, Marker_Target, iVertex, Point_Target);

          }

          /*--- If we have aggregated the values in the function RecoverTarget_Variable, the set is outside the loop ---*/
          if (valAggregated) SetTarget_Variable(target_solution, target_geometry, target_config, Marker_Target, iVertex, Point_Target);
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
  if (rank == MASTER_NODE && Buffer_Recv_mark != NULL) 
    delete [] Buffer_Recv_mark;
  #endif
}

void CTransfer::Preprocessing_InterfaceAverage(CGeometry *donor_geometry, CGeometry *target_geometry,
    CConfig *donor_config, CConfig *target_config, unsigned short iMarkerInt){

  unsigned short  nMarkerDonor, nMarkerTarget;		// Number of markers on the interface, donor and target side
  unsigned short  iMarkerDonor, iMarkerTarget;		// Variables for iteration over markers
  unsigned short iSpan,jSpan, tSpan = 0, kSpan = 0, nSpanDonor, nSpanTarget, Donor_Flag = 0, Target_Flag = 0;
  int Marker_Donor = -1, Marker_Target = -1;

  su2double *SpanValuesDonor, *SpanValuesTarget, dist, test, dist2, test2;

#ifdef HAVE_MPI
  int iSize;
  int *BuffMarkerDonor, *BuffDonorFlag;
#endif


  nMarkerDonor   = donor_geometry->GetnMarker();
  nMarkerTarget  = target_geometry->GetnMarker();
  //TODO turbo this approach only works if all the turboamchinery marker of all zones have the same amount of span wise sections.
  //TODO turbo initialization needed for the MPI routine should be place somewhere else.
  nSpanDonor     = donor_config->GetnSpanWiseSections();
  nSpanTarget    = target_config->GetnSpanWiseSections();

  /*--- On the donor side ---*/
  for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++){
    /*--- If the tag GetMarker_All_MixingPlaneInterface equals the index we are looping at ---*/
    if ( donor_config->GetMarker_All_MixingPlaneInterface(iMarkerDonor) == iMarkerInt ){
      /*--- We have identified the local index of the Donor marker ---*/
      /*--- Now we are going to store the average values that belong to Marker_Donor on each processor ---*/
      /*--- Store the identifier for the structural marker ---*/
      Marker_Donor = iMarkerDonor;
      Donor_Flag = donor_config->GetMarker_All_TurbomachineryFlag(iMarkerDonor);
      //							cout << " donor is "<< donor_config->GetMarker_All_TagBound(Marker_Donor)<<" in imarker interface "<< iMarkerInt <<endl;
      /*--- Exit the for loop: we have found the local index for Mixing-Plane interface ---*/
      break;
    }
    else {
      /*--- If the tag hasn't matched any tag within the donor markers ---*/
      Marker_Donor = -1;
      Donor_Flag   = -1;
    }
  }

#ifdef HAVE_MPI
  BuffMarkerDonor          = new int[size];
  BuffDonorFlag            = new int[size];
  for (iSize=0; iSize<size;iSize++){
    BuffMarkerDonor[iSize]            = -1;
    BuffDonorFlag[iSize]              = -1;
  }

  SU2_MPI::Allgather(&Marker_Donor, 1 , MPI_INT, BuffMarkerDonor, 1, MPI_INT, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Donor_Flag, 1 , MPI_INT, BuffDonorFlag, 1, MPI_INT, MPI_COMM_WORLD);


  Marker_Donor= -1;
  Donor_Flag= -1;


  for (iSize=0; iSize<size;iSize++){
    if(BuffMarkerDonor[iSize] > 0.0){
      Marker_Donor = BuffMarkerDonor[iSize];
      Donor_Flag   = BuffDonorFlag[iSize];
      break;
    }
  }
  delete [] BuffMarkerDonor;
  delete [] BuffDonorFlag;
#endif

  /*--- On the target side we have to identify the marker as well ---*/

  for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++){
    /*--- If the tag GetMarker_All_MixingPlaneInterface(iMarkerTarget) equals the index we are looping at ---*/
    if ( target_config->GetMarker_All_MixingPlaneInterface(iMarkerTarget) == iMarkerInt ){
      /*--- Store the identifier for the fluid marker ---*/

      // here i should then store it in the target zone

      Marker_Target = iMarkerTarget;
      Target_Flag = target_config->GetMarker_All_TurbomachineryFlag(iMarkerTarget);
      //					cout << " target is "<< target_config->GetMarker_All_TagBound(Marker_Target) <<" in imarker interface "<< iMarkerInt <<endl;
      //				/*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
      break;
    }
    else {
      /*--- If the tag hasn't matched any tag within the Flow markers ---*/
      Marker_Target = -1;
    }
  }

  if (Marker_Target != -1 && Marker_Donor != -1){

    SpanValuesDonor  = donor_geometry->GetSpanWiseValue(Donor_Flag);
    SpanValuesTarget = target_geometry->GetSpanWiseValue(Target_Flag);


    for(iSpan = 1; iSpan <nSpanTarget-1; iSpan++){
      dist  = 10E+06;
      dist2 = 10E+06;
      for(jSpan = 0; jSpan < nSpanDonor;jSpan++){
        test = abs(SpanValuesTarget[iSpan] - SpanValuesDonor[jSpan]);
        test2 = abs(SpanValuesTarget[iSpan] - SpanValuesDonor[jSpan]);
        if(test < dist && SpanValuesTarget[iSpan] > SpanValuesDonor[jSpan]){
          dist = test;
          kSpan = jSpan;
        }
        if(test2 < dist2){
          dist2 = test2;
          tSpan =jSpan;
        }

      }
      switch(donor_config->GetKind_MixingPlaneInterface()){
      case MATCHING:
        SpanLevelDonor[iSpan]        = iSpan;
        SpanValueCoeffTarget[iSpan]  = 0.0;
        break;
      case NEAREST_SPAN:
        SpanLevelDonor[iSpan]        = tSpan;
        SpanValueCoeffTarget[iSpan]  = 0.0;
        break;
      case LINEAR_INTERPOLATION:
        SpanLevelDonor[iSpan]        = kSpan;
        SpanValueCoeffTarget[iSpan]  = (SpanValuesTarget[iSpan] - SpanValuesDonor[kSpan])/(SpanValuesDonor[kSpan + 1] - SpanValuesDonor[kSpan]);
        break;
      default:
        SU2_MPI::Error("MixingPlane interface option not implemented yet", CURRENT_FUNCTION);
        break;

      }
    }
  }

}


void CTransfer::Allgather_InterfaceAverage(CSolver *donor_solution, CSolver *target_solution,
    CGeometry *donor_geometry, CGeometry *target_geometry,
    CConfig *donor_config, CConfig *target_config, unsigned short iMarkerInt){
  unsigned short  nMarkerDonor, nMarkerTarget;		// Number of markers on the interface, donor and target side
  unsigned short  iMarkerDonor, iMarkerTarget;		// Variables for iteration over markers
  unsigned short iSpan, nSpanDonor, nSpanTarget;
  int Marker_Donor = -1, Marker_Target = -1;
  su2double *avgPressureDonor = NULL, *avgDensityDonor = NULL, *avgNormalVelDonor = NULL,
      *avgTangVelDonor = NULL, *avg3DVelDonor = NULL, *avgNuDonor = NULL, *avgOmegaDonor = NULL, *avgKineDonor = NULL;
  su2double *avgPressureTarget = NULL, *avgDensityTarget = NULL, *avgNormalVelTarget = NULL,
      *avg3DVelTarget = NULL, *avgTangVelTarget = NULL, *avgNuTarget = NULL, *avgOmegaTarget = NULL, *avgKineTarget = NULL;

#ifdef HAVE_MPI
  int iSize;
  su2double *BuffAvgPressureDonor = NULL, *BuffAvgDensityDonor = NULL, *BuffAvgNormalVelDonor = NULL, *BuffAvg3DVelDonor = NULL,
      *BuffAvgTangVelDonor = NULL, *BuffAvgNuDonor = NULL, *BuffAvgKineDonor = NULL, *BuffAvgOmegaDonor = NULL;
  int nSpanSize, *BuffMarkerDonor;
#endif


  nMarkerTarget  = target_geometry->GetnMarker();
  nMarkerDonor   = donor_geometry->GetnMarker();
  nSpanDonor     = donor_config->GetnSpanWiseSections() +1;
  nSpanTarget    = target_config->GetnSpanWiseSections() +1;


  avgDensityDonor                  = new su2double[nSpanDonor];
  avgPressureDonor                 = new su2double[nSpanDonor];
  avgNormalVelDonor                = new su2double[nSpanDonor];
  avgTangVelDonor                  = new su2double[nSpanDonor];
  avg3DVelDonor                    = new su2double[nSpanDonor];
  avgNuDonor                       = new su2double[nSpanDonor];
  avgKineDonor                     = new su2double[nSpanDonor];
  avgOmegaDonor                    = new su2double[nSpanDonor];

  for (iSpan = 0; iSpan < nSpanDonor; iSpan++){
    avgDensityDonor[iSpan]         = -1.0;
    avgPressureDonor[iSpan]        = -1.0;
    avgNormalVelDonor[iSpan]       = -1.0;
    avgTangVelDonor[iSpan]         = -1.0;
    avg3DVelDonor[iSpan]           = -1.0;
    avgNuDonor[iSpan]              = -1.0;
    avgKineDonor[iSpan]            = -1.0;
    avgOmegaDonor[iSpan]           = -1.0;
  }

  avgDensityTarget                 = new su2double[nSpanTarget];
  avgPressureTarget                = new su2double[nSpanTarget];
  avgNormalVelTarget               = new su2double[nSpanTarget];
  avgTangVelTarget                 = new su2double[nSpanTarget];
  avg3DVelTarget                   = new su2double[nSpanTarget];
  avgNuTarget                      = new su2double[nSpanTarget];
  avgKineTarget                    = new su2double[nSpanTarget];
  avgOmegaTarget                   = new su2double[nSpanTarget];


  for (iSpan = 0; iSpan < nSpanTarget; iSpan++){
    avgDensityTarget[iSpan]        = -1.0;
    avgPressureTarget[iSpan]       = -1.0;
    avgNormalVelTarget[iSpan]      = -1.0;
    avgTangVelTarget[iSpan]        = -1.0;
    avg3DVelTarget[iSpan]          = -1.0;
    avgNuTarget[iSpan]             = -1.0;
    avgKineTarget[iSpan]           = -1.0;
    avgOmegaTarget[iSpan]          = -1.0;
  }

  /*--- Outer loop over the markers on the Mixing-Plane interface: compute one by one ---*/
  /*--- The tags are always an integer greater than 1: loop from 1 to nMarkerMixingPlane ---*/
  Marker_Donor = -1;
  Marker_Target = -1;

  /*--- The donor and target markers are tagged with the same index.
   *--- This is independent of the MPI domain decomposition.
   *--- We need to loop over all markers on both sides  ---*/

  /*--- On the donor side ---*/

  for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++){
    /*--- If the tag GetMarker_All_MixingPlaneInterface equals the index we are looping at ---*/
    if ( donor_config->GetMarker_All_MixingPlaneInterface(iMarkerDonor) == iMarkerInt ){
      /*--- We have identified the local index of the Donor marker ---*/
      /*--- Now we are going to store the average values that belong to Marker_Donor on each processor ---*/
      /*--- Store the identifier for the structural marker ---*/
      Marker_Donor = iMarkerDonor;
      /*--- Exit the for loop: we have found the local index for Mixing-Plane interface ---*/
      break;
    }
    else {
      /*--- If the tag hasn't matched any tag within the donor markers ---*/
      Marker_Donor = -1;
    }
  }
  /*--- Here we want to make available the quantities for all the processors and collect them in a buffer
   * for each span of the donor the span-wise height vector also so that then we can interpolate on the target side  ---*/
  if (Marker_Donor != -1){
    for(iSpan = 0; iSpan < nSpanDonor; iSpan++){
      GetDonor_Variable(donor_solution, donor_geometry, donor_config, Marker_Donor, iSpan, rank);
      avgDensityDonor[iSpan]          = Donor_Variable[0];
      avgPressureDonor[iSpan]         = Donor_Variable[1];
      avgNormalVelDonor[iSpan]        = Donor_Variable[2];
      avgTangVelDonor[iSpan]          = Donor_Variable[3];
      avg3DVelDonor[iSpan]            = Donor_Variable[4];
      avgNuDonor[iSpan]               = Donor_Variable[5];
      avgKineDonor[iSpan]             = Donor_Variable[6];
      avgOmegaDonor[iSpan]            = Donor_Variable[7];
    }
  }

#ifdef HAVE_MPI
  nSpanSize = size*nSpanDonor;
  BuffAvgDensityDonor                 = new su2double[nSpanSize];
  BuffAvgPressureDonor                = new su2double[nSpanSize];
  BuffAvgNormalVelDonor               = new su2double[nSpanSize];
  BuffAvgTangVelDonor                 = new su2double[nSpanSize];
  BuffAvg3DVelDonor                   = new su2double[nSpanSize];
  BuffAvgNuDonor                      = new su2double[nSpanSize];
  BuffAvgKineDonor                    = new su2double[nSpanSize];
  BuffAvgOmegaDonor                   = new su2double[nSpanSize];
  BuffMarkerDonor                     = new int[size];

  for (iSpan=0;iSpan<nSpanSize;iSpan++){
    BuffAvgDensityDonor[iSpan]        = -1.0;
    BuffAvgPressureDonor[iSpan]       = -1.0;
    BuffAvgNormalVelDonor[iSpan]      = -1.0;
    BuffAvgTangVelDonor[iSpan]        = -1.0;
    BuffAvg3DVelDonor[iSpan]          = -1.0;
    BuffAvgNuDonor[iSpan]             = -1.0;
    BuffAvgKineDonor[iSpan]           = -1.0;
    BuffAvgOmegaDonor[iSpan]          = -1.0;
  }

  for (iSize=0; iSize<size;iSize++){
    BuffMarkerDonor[iSize]            = -1;
  }

  SU2_MPI::Allgather(avgDensityDonor, nSpanDonor , MPI_DOUBLE, BuffAvgDensityDonor, nSpanDonor, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(avgPressureDonor, nSpanDonor , MPI_DOUBLE, BuffAvgPressureDonor, nSpanDonor, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(avgNormalVelDonor, nSpanDonor , MPI_DOUBLE, BuffAvgNormalVelDonor, nSpanDonor, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(avgTangVelDonor, nSpanDonor , MPI_DOUBLE, BuffAvgTangVelDonor, nSpanDonor, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(avg3DVelDonor, nSpanDonor , MPI_DOUBLE, BuffAvg3DVelDonor, nSpanDonor, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(avgNuDonor, nSpanDonor , MPI_DOUBLE, BuffAvgNuDonor, nSpanDonor, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(avgKineDonor, nSpanDonor , MPI_DOUBLE, BuffAvgKineDonor, nSpanDonor, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(avgOmegaDonor, nSpanDonor , MPI_DOUBLE, BuffAvgOmegaDonor, nSpanDonor, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Marker_Donor, 1 , MPI_INT, BuffMarkerDonor, 1, MPI_INT, MPI_COMM_WORLD);

  for (iSpan = 0; iSpan < nSpanDonor; iSpan++){
    avgDensityDonor[iSpan]            = -1.0;
    avgPressureDonor[iSpan]           = -1.0;
    avgNormalVelDonor[iSpan]          = -1.0;
    avgTangVelDonor[iSpan]            = -1.0;
    avg3DVelDonor[iSpan]              = -1.0;
    avgNuDonor[iSpan]                 = -1.0;
    avgKineDonor[iSpan]               = -1.0;
    avgOmegaDonor[iSpan]              = -1.0;
  }

  Marker_Donor= -1;

  for (iSize=0; iSize<size;iSize++){
    if(BuffAvgDensityDonor[nSpanDonor*iSize] > 0.0){
      for (iSpan = 0; iSpan < nSpanDonor; iSpan++){
        avgDensityDonor[iSpan]        = BuffAvgDensityDonor[nSpanDonor*iSize + iSpan];
        avgPressureDonor[iSpan]       = BuffAvgPressureDonor[nSpanDonor*iSize + iSpan];
        avgNormalVelDonor[iSpan]      = BuffAvgNormalVelDonor[nSpanDonor*iSize + iSpan];
        avgTangVelDonor[iSpan]        = BuffAvgTangVelDonor[nSpanDonor*iSize + iSpan];
        avg3DVelDonor[iSpan]          = BuffAvg3DVelDonor[nSpanDonor*iSize + iSpan];
        avgNuDonor[iSpan]             = BuffAvgNuDonor[nSpanDonor*iSize + iSpan];
        avgKineDonor[iSpan]           = BuffAvgKineDonor[nSpanDonor*iSize + iSpan];
        avgOmegaDonor[iSpan]          = BuffAvgOmegaDonor[nSpanDonor*iSize + iSpan];
      }
      Marker_Donor                    = BuffMarkerDonor[iSize];
      break;
    }
  }
  delete [] BuffAvgDensityDonor;
  delete [] BuffAvgPressureDonor;
  delete [] BuffAvgNormalVelDonor;
  delete [] BuffAvgTangVelDonor;
  delete [] BuffAvg3DVelDonor;
  delete [] BuffAvgNuDonor;
  delete [] BuffAvgKineDonor;
  delete [] BuffAvgOmegaDonor;
  delete [] BuffMarkerDonor;

#endif

  /*--- On the target side we have to identify the marker as well ---*/
  for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++){
    /*--- If the tag GetMarker_All_MixingPlaneInterface(iMarkerTarget) equals the index we are looping at ---*/
    if ( target_config->GetMarker_All_MixingPlaneInterface(iMarkerTarget) == iMarkerInt ){
      /*--- Store the identifier for the fluid marker ---*/
      Marker_Target = iMarkerTarget;
      /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
      break;
    }
    else {
      /*--- If the tag hasn't matched any tag within the Flow markers ---*/
      Marker_Target = -1;
    }
  }


  if (Marker_Target != -1 && Marker_Donor != -1){

    /*--- linear interpolation of the average value of for the internal span-wise levels ---*/
    for(iSpan = 1; iSpan < nSpanTarget -2 ; iSpan++){
      avgDensityTarget[iSpan]                = SpanValueCoeffTarget[iSpan]*(avgDensityDonor[SpanLevelDonor[iSpan] + 1] - avgDensityDonor[SpanLevelDonor[iSpan]]);
      avgDensityTarget[iSpan]               += avgDensityDonor[SpanLevelDonor[iSpan]];
      avgPressureTarget[iSpan]               = SpanValueCoeffTarget[iSpan]*(avgPressureDonor[SpanLevelDonor[iSpan] + 1] - avgPressureDonor[SpanLevelDonor[iSpan]]);
      avgPressureTarget[iSpan]              += avgPressureDonor[SpanLevelDonor[iSpan]];
      avgNormalVelTarget[iSpan]              = SpanValueCoeffTarget[iSpan]*(avgNormalVelDonor[SpanLevelDonor[iSpan] + 1] - avgNormalVelDonor[SpanLevelDonor[iSpan]]);
      avgNormalVelTarget[iSpan]             += avgNormalVelDonor[SpanLevelDonor[iSpan]];
      avgTangVelTarget[iSpan]                = SpanValueCoeffTarget[iSpan]*(avgTangVelDonor[SpanLevelDonor[iSpan] + 1] - avgTangVelDonor[SpanLevelDonor[iSpan]]);
      avgTangVelTarget[iSpan]               += avgTangVelDonor[SpanLevelDonor[iSpan]];
      avg3DVelTarget[iSpan]                  = SpanValueCoeffTarget[iSpan]*(avg3DVelDonor[SpanLevelDonor[iSpan] + 1] - avg3DVelDonor[SpanLevelDonor[iSpan]]);
      avg3DVelTarget[iSpan]                 += avg3DVelDonor[SpanLevelDonor[iSpan]];
      avgNuTarget[iSpan]                     = SpanValueCoeffTarget[iSpan]*(avgNuDonor[SpanLevelDonor[iSpan] + 1] - avgNuDonor[SpanLevelDonor[iSpan]]);
      avgNuTarget[iSpan]                    += avgNuDonor[SpanLevelDonor[iSpan]];
      avgKineTarget[iSpan]                   = SpanValueCoeffTarget[iSpan]*(avgKineDonor[SpanLevelDonor[iSpan] + 1] - avgKineDonor[SpanLevelDonor[iSpan]]);
      avgKineTarget[iSpan]                  += avgKineDonor[SpanLevelDonor[iSpan]];
      avgOmegaTarget[iSpan]                  = SpanValueCoeffTarget[iSpan]*(avgOmegaDonor[SpanLevelDonor[iSpan] + 1] - avgOmegaDonor[SpanLevelDonor[iSpan] ]);
      avgOmegaTarget[iSpan]                 += avgOmegaDonor[SpanLevelDonor[iSpan]];
    }


    /*--- transfer values at the hub ---*/
    avgDensityTarget[0]                      = avgDensityDonor[0];
    avgPressureTarget[0]                     = avgPressureDonor[0];
    avgNormalVelTarget[0]                    = avgNormalVelDonor[0];
    avgTangVelTarget[0]                      = avgTangVelDonor[0];
    avg3DVelTarget[0]                        = avg3DVelDonor[0];
    avgNuTarget[0]                           = avgNuDonor[0];
    avgKineTarget[0]                         = avgKineDonor[0];
    avgOmegaTarget[0]                        = avgOmegaDonor[0];

    /*--- transfer values at the shroud ---*/
    avgDensityTarget[nSpanTarget - 2]        = avgDensityDonor[nSpanDonor - 2];
    avgPressureTarget[nSpanTarget - 2]       = avgPressureDonor[nSpanDonor - 2];
    avgNormalVelTarget[nSpanTarget - 2]      = avgNormalVelDonor[nSpanDonor - 2];
    avgTangVelTarget[nSpanTarget - 2]        = avgTangVelDonor[nSpanDonor - 2];
    avg3DVelTarget[nSpanTarget - 2]          = avg3DVelDonor[nSpanDonor - 2];
    avgNuTarget[nSpanTarget - 2]             = avgNuDonor[nSpanDonor - 2];
    avgKineTarget[nSpanTarget - 2]           = avgKineDonor[nSpanDonor - 2];
    avgOmegaTarget[nSpanTarget - 2]          = avgOmegaDonor[nSpanDonor - 2];

    /*--- transfer 1D values ---*/
    avgDensityTarget[nSpanTarget - 1]        = avgDensityDonor[nSpanDonor - 1];
    avgPressureTarget[nSpanTarget - 1]       = avgPressureDonor[nSpanDonor - 1];
    avgNormalVelTarget[nSpanTarget - 1]      = avgNormalVelDonor[nSpanDonor - 1];
    avgTangVelTarget[nSpanTarget - 1]        = avgTangVelDonor[nSpanDonor - 1];
    avg3DVelTarget[nSpanTarget - 1]          = avg3DVelDonor[nSpanDonor - 1];
    avgNuTarget[nSpanTarget - 1]             = avgNuDonor[nSpanDonor - 1];
    avgKineTarget[nSpanTarget - 1]           = avgKineDonor[nSpanDonor - 1];
    avgOmegaTarget[nSpanTarget - 1]          = avgOmegaDonor[nSpanDonor - 1];


    /*---finally, the interpolated value is sent  to the target zone ---*/
    for(iSpan = 0; iSpan < nSpanTarget ; iSpan++){
      Target_Variable[0]                     = avgDensityTarget[iSpan];
      Target_Variable[1]                     = avgPressureTarget[iSpan];
      Target_Variable[2]                     = avgNormalVelTarget[iSpan];
      Target_Variable[3]                     = avgTangVelTarget[iSpan];
      Target_Variable[4]                     = avg3DVelTarget[iSpan];
      Target_Variable[5]                     = avgNuTarget[iSpan];
      Target_Variable[6]                     = avgKineTarget[iSpan];
      Target_Variable[7]                     = avgOmegaTarget[iSpan];


      SetTarget_Variable(target_solution, target_geometry, target_config, Marker_Target, iSpan, rank);
    }
  }

  delete [] avgDensityDonor;
  delete [] avgPressureDonor;
  delete [] avgNormalVelDonor;
  delete [] avgTangVelDonor;
  delete [] avg3DVelDonor;
  delete [] avgNuDonor;
  delete [] avgKineDonor;
  delete [] avgOmegaDonor;


  delete [] avgDensityTarget;
  delete [] avgPressureTarget;
  delete [] avgNormalVelTarget;
  delete [] avgTangVelTarget;
  delete [] avg3DVelTarget;
  delete [] avgNuTarget;
  delete [] avgKineTarget;
  delete [] avgOmegaTarget;


}

void CTransfer::GatherAverageValues(CSolver *donor_solution, CSolver *target_solution, unsigned short donorZone){


  /*--- here we made the strong assumption that the mesh zone order follow the same order of the turbomachinery markers ---*/
  SetAverageValues(donor_solution, target_solution, donorZone);

}

void CTransfer::GatherAverageTurboGeoValues(CGeometry *donor_geometry, CGeometry *target_geometry, unsigned short donorZone){


  /*--- here we made the strong assumption that the mesh zone order follow the same order of the turbomachinery markers ---*/
  SetAverageTurboGeoValues(donor_geometry, target_geometry, donorZone);

}

