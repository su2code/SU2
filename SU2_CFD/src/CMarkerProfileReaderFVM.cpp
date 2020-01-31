/*!
 * \file CMarkerProfileReaderFVM.cpp
 * \brief Class that handles the reading of marker profile files.
 * \author T. Economon
 * \version 7.0.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../include/CMarkerProfileReaderFVM.hpp"

CMarkerProfileReaderFVM::CMarkerProfileReaderFVM(CGeometry      *val_geometry,
                                                 CConfig        *val_config,
                                                 string         val_filename,
                                                 unsigned short val_kind_marker,
                                                 unsigned short val_number_vars) {
  
  /*--- Store input values and pointers to class data. ---*/
  
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  this->config    = val_config;
  this->geometry  = val_geometry;
  dimension       = geometry->GetnDim();
  
  filename     = val_filename;
  markerType   = val_kind_marker;
  numberOfVars = val_number_vars;
  
  /* Attempt to open the specified file. */
  ifstream profile_file;
  profile_file.open(filename.data(), ios::in);
  
  /* If the file is not found, then we merge the information necessary
   and write a template marker profile file. Otherwise, we read and
   store the information in the marker profile file. */
  
  if (profile_file.fail()) {
    MergeProfileMarkers();
    WriteMarkerProfileTemplate();
  } else {
    ReadMarkerProfile();
  }
  
}

CMarkerProfileReaderFVM::~CMarkerProfileReaderFVM(void) { }

void CMarkerProfileReaderFVM::ReadMarkerProfile() {
  
  /*--- Open the profile file (we have already error checked) ---*/
  
  ifstream profile_file;
  profile_file.open(filename.data(), ios::in);
  
  /*--- Identify the markers and data set in the profile file ---*/
  
  string text_line;
  while (getline (profile_file, text_line)) {
    
    string::size_type position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      text_line.erase (0,6); numberOfProfiles = atoi(text_line.c_str());
      
      numberOfRowsInProfile.resize(numberOfProfiles);
      numberOfColumnsInProfile.resize(numberOfProfiles);
      
      for (unsigned short iMarker = 0 ; iMarker < numberOfProfiles; iMarker++) {
        
        getline (profile_file, text_line);
        text_line.erase (0,11);
        for (unsigned short iChar = 0; iChar < 20; iChar++) {
          position = text_line.find( " ", 0 );  if (position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\r", 0 ); if (position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\n", 0 ); if (position != string::npos) text_line.erase (position,1);
        }
        profileTags.push_back(text_line.c_str());
        
        getline (profile_file, text_line);
        text_line.erase (0,5); numberOfRowsInProfile[iMarker] = atoi(text_line.c_str());
        
        getline (profile_file, text_line);
        text_line.erase (0,5); numberOfColumnsInProfile[iMarker] = atoi(text_line.c_str());
        
        /*--- Skip the data. This is read in the next loop. ---*/
        
        for (unsigned long iRow = 0; iRow < numberOfRowsInProfile[iMarker]; iRow++) getline (profile_file, text_line);
        
      }
    } else {
      SU2_MPI::Error("While opening profile file, no \"NMARK=\" specification was found", CURRENT_FUNCTION);
    }
  }
  
  profile_file.close();
  
  /*--- Compute array bounds and offsets. Allocate data structure. ---*/
  
  profileData.resize(numberOfProfiles);
  for (unsigned short iMarker = 0; iMarker < numberOfProfiles; iMarker++) {
    profileData[iMarker].resize(numberOfRowsInProfile[iMarker]*numberOfColumnsInProfile[iMarker], 0.0);
  }
  
  /*--- Read all lines in the profile file and extract data. ---*/
  
  profile_file.open(filename.data(), ios::in);
  
  int counter = 0;
  while (getline (profile_file, text_line)) {
    
    string::size_type position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      
      for (unsigned short iMarker = 0; iMarker < numberOfProfiles; iMarker++) {
        
        /*--- Skip the tag, nRow, and nCol lines. ---*/
        
        getline (profile_file, text_line);
        getline (profile_file, text_line);
        getline (profile_file, text_line);
        
        /*--- Now read the data for each row and store. ---*/
        
        for (unsigned long iRow = 0; iRow < numberOfRowsInProfile[iMarker]; iRow++) {
          
          getline (profile_file, text_line);
          
          istringstream point_line(text_line);
          
          /*--- Store the values (starting with node coordinates) --*/
          
          for (unsigned short iVar = 0; iVar < numberOfColumnsInProfile[iMarker]; iVar++)
            point_line >> profileData[iMarker][iRow*numberOfColumnsInProfile[iMarker] + iVar];
          
          /*--- Increment our local row counter. ---*/
          
          counter++;
          
        }
      }
    }
  }
  
  profile_file.close();
  
}

void CMarkerProfileReaderFVM::MergeProfileMarkers() {
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned long iPoint, jPoint, kPoint;
  
  int iProcessor, nProcessor = size;
  
  unsigned long iVertex, iMarker;
  unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  
  unsigned long index, iChar;
  
  char str_buf[MAX_STRING_SIZE];
  vector<string> Marker_Tags;
  vector<string>::iterator it;
  
  vector<unsigned long> nRowCum_Counter;
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];
  
  /*--- Search all boundaries on the present rank to count the number
   of nodes found on profile markers. ---*/
  
  nLocalPoint = 0; numberOfProfiles = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == markerType) {
      numberOfProfiles++;
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- Only communicate owned nodes to avoid duplicates. ---*/
        
        if (geometry->node[iPoint]->GetDomain())
          nLocalPoint++;
        
      }
    }
  }
  Buffer_Send_nPoin[0] = nLocalPoint;
  
  /*--- Communicate the total number of nodes on this domain. ---*/
  
  SU2_MPI::Gather(&Buffer_Send_nPoin, 1, MPI_UNSIGNED_LONG,
                  Buffer_Recv_nPoin, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  
  /*--- Send and Recv buffers. ---*/
  
  su2double *Buffer_Send_X = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_X = NULL;
  
  su2double *Buffer_Send_Y = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Y = NULL;
  
  su2double *Buffer_Send_Z = NULL, *Buffer_Recv_Z = NULL;
  if (dimension == 3) Buffer_Send_Z = new su2double[MaxLocalPoint];
  
  char *Buffer_Send_Str = new char[MaxLocalPoint*MAX_STRING_SIZE];
  char *Buffer_Recv_Str = NULL;
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_X = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Y = new su2double[nProcessor*MaxLocalPoint];
    if (dimension == 3) Buffer_Recv_Z = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Str = new char[nProcessor*MaxLocalPoint*MAX_STRING_SIZE];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    
    unsigned long nGlobal_InletPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_InletPoint += Buffer_Recv_nPoin[iProcessor];
    }
    
    profileCoords.resize(numberOfProfiles);
    for (iMarker = 0; iMarker < numberOfProfiles; iMarker++) {
      profileCoords[iMarker].resize(dimension);
    }
    
  }
  
  /*--- Main communication routine. Loop over each coordinate and perform
   the MPI comm. Temporary 1-D buffers are used to send the coordinates at
   all nodes on each partition to the master node. These are then unpacked
   by the master and sorted by marker tag in one large n-dim. array. ---*/
  
  su2double *Coords_Local; jPoint = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == markerType) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- Only communicate owned nodes to avoid duplicates. ---*/
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Retrieve local coordinates at this node. ---*/
          
          Coords_Local = geometry->node[iPoint]->GetCoord();
          
          /*--- Load local coords into the temporary send buffer. ---*/
          
          Buffer_Send_X[jPoint] = Coords_Local[0];
          Buffer_Send_Y[jPoint] = Coords_Local[1];
          if (dimension == 3) Buffer_Send_Z[jPoint] = Coords_Local[2];
          
          /*--- If US system, the output should be in inches ---*/
          
          if (config->GetSystemMeasurements() == US) {
            Buffer_Send_X[jPoint] *= 12.0;
            Buffer_Send_Y[jPoint] *= 12.0;
            if (dimension == 3) Buffer_Send_Z[jPoint] *= 12.0;
          }
          
          /*--- Store the marker tag for this particular node. ---*/
          
          SPRINTF(&Buffer_Send_Str[jPoint*MAX_STRING_SIZE], "%s",
                  config->GetMarker_All_TagBound(iMarker).c_str());
          
          /*--- Increment jPoint as the counter. We need this because iPoint
           may include halo nodes that we skip over during this loop. ---*/
          
          jPoint++;
          
        }
      }
    }
  }
  
  /*--- Gather the coordinate data on the master node using MPI. ---*/
  
  SU2_MPI::Gather(Buffer_Send_X, (int)MaxLocalPoint, MPI_DOUBLE,
                  Buffer_Recv_X, (int)MaxLocalPoint, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Y, (int)MaxLocalPoint, MPI_DOUBLE,
                  Buffer_Recv_Y, (int)MaxLocalPoint, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (dimension == 3) {
    SU2_MPI::Gather(Buffer_Send_Z, (int)MaxLocalPoint, MPI_DOUBLE,
                    Buffer_Recv_Z, (int)MaxLocalPoint, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  }
  SU2_MPI::Gather(Buffer_Send_Str, (int)MaxLocalPoint*MAX_STRING_SIZE, MPI_CHAR,
                  Buffer_Recv_Str, (int)MaxLocalPoint*MAX_STRING_SIZE, MPI_CHAR, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- The master node unpacks and sorts this variable by marker tag. ---*/
  
  if (rank == MASTER_NODE) {
    
    profileTags.clear();
    
    /*--- First, parse the marker tags to count how many total profile markers
     we have now on the master. ---*/
    
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        index = (iProcessor*MaxLocalPoint + iPoint)*MAX_STRING_SIZE;
        for (iChar = 0; iChar < MAX_STRING_SIZE; iChar++) {
          str_buf[iChar] = Buffer_Recv_Str[index + iChar];
        }
        Marker_Tags.push_back(str_buf);
        profileTags.push_back(str_buf);
      }
    }
    
    /*--- Sort and remove the duplicate profile marker strings. ---*/
    
    sort(profileTags.begin(), profileTags.end());
    profileTags.erase(unique(profileTags.begin(),
                             profileTags.end()),
                      profileTags.end());
    
    /*--- Store the unique number of markers for writing later. ---*/
    
    numberOfProfiles = profileTags.size();
    
    /*--- Count the number of rows (nodes) per marker. ---*/
    
    numberOfRowsInProfile.resize(numberOfProfiles,0.0);
    jPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        for (iMarker = 0; iMarker < numberOfProfiles; iMarker++) {
          if (profileTags[iMarker] == Marker_Tags[jPoint]) {
            numberOfRowsInProfile[iMarker]++;
          }
        }
        jPoint++;
      }
    }

    /*--- Load up the coordinates, sorted into chunks per marker. ---*/
    
    jPoint = 0; kPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        for (iMarker = 0; iMarker < numberOfProfiles; iMarker++) {
          
          if (profileTags[iMarker] == Marker_Tags[kPoint]) {
            
            /*--- Find our current index for this marker and store coords. ---*/
            
            profileCoords[iMarker][0].push_back(Buffer_Recv_X[jPoint]);
            profileCoords[iMarker][1].push_back(Buffer_Recv_Y[jPoint]);
            if (dimension == 3)
              profileCoords[iMarker][2].push_back(Buffer_Recv_Z[jPoint]);
            
          }
        }
        
        /*--- Increment point counter for marker tags and data. ---*/
        
        kPoint++;
        jPoint++;
        
      }
      
      /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
      
      jPoint = (iProcessor+1)*MaxLocalPoint;
      
    }
  }
  
  /*--- Immediately release the temporary data buffers. ---*/
  
  delete [] Buffer_Send_X;
  delete [] Buffer_Send_Y;
  if (Buffer_Send_Z != NULL) delete [] Buffer_Send_Z;
  delete [] Buffer_Send_Str;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_X;
    delete [] Buffer_Recv_Y;
    if (Buffer_Recv_Z != NULL)  delete [] Buffer_Recv_Z;
    delete [] Buffer_Recv_nPoin;
    delete [] Buffer_Recv_Str;
  }
  
}

void CMarkerProfileReaderFVM::WriteMarkerProfileTemplate() {
  
  /*--- Count the number of columns that we have for this case.
   Here, we have dimension entries for node coordinates and then the
   total number of columns of data specified in the constructor. ---*/
  
  const unsigned short nColumns = dimension + numberOfVars;
  
  /*--- Write the profile file. Note that we have already merged
   all of the information for the markers and coordinates previously
   in the MergeProfileMarkers() routine and only the master writes. ---*/
  
  if (rank == MASTER_NODE) {
    
    ofstream node_file("profile_example.dat");
    
    node_file << "NMARK= " << numberOfProfiles << endl;
    
    for (unsigned short iMarker = 0; iMarker < numberOfProfiles; iMarker++) {
      
      /*--- Access the default data for this marker. ---*/
      
      string Marker_Tag = profileTags[iMarker];
      
      /*--- Header information for this marker. ---*/
      
      node_file << "MARKER_TAG= " << Marker_Tag              << endl;
      node_file << "NROW="        << numberOfRowsInProfile[iMarker] << endl;
      node_file << "NCOL="        << nColumns          << endl;
      
      node_file << setprecision(15);
      node_file << std::scientific;
      
      /*--- Loop over the data structure and write the coords and vars. ---*/
      
      for (unsigned long iPoint = 0; iPoint < numberOfRowsInProfile[iMarker]; iPoint++) {
        
        for (unsigned short iDim = 0; iDim < dimension; iDim++) {
          node_file << profileCoords[iMarker][iDim][iPoint] << "\t";
        }
        for (unsigned short iDim = 0; iDim < numberOfVars; iDim++) {
          node_file << 0.0 << "\t";
        }
        node_file << endl;
      }
      
    }
    node_file.close();
    
    /*--- Print a message to inform the user about the template file. ---*/
    
    stringstream err;
    err << endl;
    err << "  Could not find the input file for the marker profile." << endl;
    err << "  Looked for: " << filename << "." << endl;
    err << "  Created a template profile file with node coordinates" << endl;
    err << "  and correct number of columns at `profile_example.dat`." << endl;
    err << "  You can use this file as a guide for making your own profile" << endl;
    err << "  specification." << endl << endl;
    SU2_MPI::Error(err.str(), CURRENT_FUNCTION);
    
  }
  
}
