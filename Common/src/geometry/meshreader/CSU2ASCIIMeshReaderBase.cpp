/*!
 * \file CSU2ASCIIMeshReaderBase.cpp
 * \brief Helper class for the reading of a native SU2 ASCII grid file.
 * \author T. Economon
 * \version 7.0.5 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/geometry/meshreader/CSU2ASCIIMeshReaderBase.hpp"

CSU2ASCIIMeshReaderBase::CSU2ASCIIMeshReaderBase(CConfig        *val_config,
                                                 unsigned short val_iZone,
                                                 unsigned short val_nZone)
: CMeshReader(val_config, val_iZone, val_nZone) {
  
  /* Store the current zone to be read and the total number of zones. */
  myZone = val_iZone;
  nZones = val_nZone;
  
  /* Store the mesh filename since we will open/close multiple times. */
  meshFilename = config->GetMesh_FileName();
  
}

CSU2ASCIIMeshReaderBase::~CSU2ASCIIMeshReaderBase(void) { }

void CSU2ASCIIMeshReaderBase::ReadMetadata() {
  
  bool harmonic_balance = config->GetTime_Marching() == HARMONIC_BALANCE;
  bool multizone_file = config->GetMultizone_Mesh();
  
  /*--- Open grid file ---*/
  
  mesh_file.open(meshFilename.c_str(), ios::in);
  if (mesh_file.fail()) {
    SU2_MPI::Error(string("Error opening SU2 ASCII grid.") +
                   string(" \n Check that the file exists."), CURRENT_FUNCTION);
  }
  
  /*--- If more than one, find the curent zone in the mesh file. ---*/
  
  string text_line;
  string::size_type position;
  if ((nZones > 1 && multizone_file) || harmonic_balance) {
    if (harmonic_balance) {
      if (rank == MASTER_NODE)
        cout << "Reading time instance " << config->GetiInst()+1 << "." << endl;
    } else {
      bool foundZone = false;
      while (getline (mesh_file,text_line)) {
        /*--- Search for the current domain ---*/
        position = text_line.find ("IZONE=",0);
        if (position != string::npos) {
          text_line.erase (0,6);
          unsigned short jZone = atoi(text_line.c_str());
          if (jZone == myZone+1) {
            if (rank == MASTER_NODE)
              cout << "Reading zone " << myZone << " from native SU2 ASCII mesh." << endl;
            foundZone = true;
            break;
          }
        }
      }
      if (!foundZone) {
        SU2_MPI::Error(string("Could not find the IZONE= keyword or the zone contents.") +
                       string(" \n Check the SU2 ASCII file format."),
                       CURRENT_FUNCTION);
      }
    }
  }
  
  /*--- Read the metadata: problem dimension, offsets for angle
   of attack and angle of sideslip, global points, global elements,
   and number of markers. Perform error checks as we go. ---*/
  
  bool foundNDIME = false, foundNPOIN = false;
  bool foundNELEM = false, foundNMARK = false;
  
  while (getline (mesh_file, text_line)) {
    
    /*--- Read the dimension of the problem ---*/
    
    position = text_line.find ("NDIME=",0);
    if (position != string::npos) {
      text_line.erase (0,6);
      dimension = atoi(text_line.c_str());
      foundNDIME = true;
    }
    
    /*--- The AoA and AoS offset values are optional. ---*/
    
    position = text_line.find ("AOA_OFFSET=",0);
    if (position != string::npos) {
      su2double AoA_Offset = 0.0;
      text_line.erase (0,11);
      AoA_Offset = atof(text_line.c_str());
      
      /*--- The offset is in deg ---*/
      
      su2double AoA_Current = config->GetAoA() + AoA_Offset;
      
      if (config->GetDiscard_InFiles() == false) {
        if ((rank == MASTER_NODE) && (AoA_Offset != 0.0))  {
          cout.precision(6);
          cout << fixed <<"WARNING: AoA in the config file (" << config->GetAoA() << " deg.) +" << endl;
          cout << "         AoA offset in mesh file (" << AoA_Offset << " deg.) = " << AoA_Current << " deg." << endl;
        }
        config->SetAoA_Offset(AoA_Offset);
        config->SetAoA(AoA_Current);
      }
      else {
        if ((rank == MASTER_NODE) && (AoA_Offset != 0.0))
          cout <<"WARNING: Discarding the AoA offset in the geometry file." << endl;
      }
      
    }
    
    position = text_line.find ("AOS_OFFSET=",0);
    if (position != string::npos) {
      su2double AoS_Offset = 0.0;
      text_line.erase (0,11);
      AoS_Offset = atof(text_line.c_str());
      
      /*--- The offset is in deg ---*/
      
      su2double AoS_Current = config->GetAoS() + AoS_Offset;
      
      if (config->GetDiscard_InFiles() == false) {
        if ((rank == MASTER_NODE) && (AoS_Offset != 0.0))  {
          cout.precision(6);
          cout << fixed <<"WARNING: AoS in the config file (" << config->GetAoS() << " deg.) +" << endl;
          cout << "         AoS offset in mesh file (" << AoS_Offset << " deg.) = " << AoS_Current << " deg." << endl;
        }
        config->SetAoS_Offset(AoS_Offset);
        config->SetAoS(AoS_Current);
      }
      else {
        if ((rank == MASTER_NODE) && (AoS_Offset != 0.0))
          cout <<"WARNING: Discarding the AoS offset in the geometry file." << endl;
      }
      
    }
    
    position = text_line.find ("NPOIN=",0);
    if (position != string::npos) {
      text_line.erase (0,6);
      numberOfGlobalPoints = atoi(text_line.c_str());
      for (unsigned long iPoint = 0; iPoint < numberOfGlobalPoints; iPoint++)
        getline (mesh_file, text_line);
      foundNPOIN = true;
    }
    
    position = text_line.find ("NELEM=",0);
    if (position != string::npos) {
      text_line.erase (0,6);
      numberOfGlobalElements = atoi(text_line.c_str());
      for (unsigned long iElem = 0; iElem < numberOfGlobalElements; iElem++)
        getline (mesh_file, text_line);
      foundNELEM = true;
    }
    
    position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      text_line.erase (0,6);
      numberOfMarkers = atoi(text_line.c_str());
      foundNMARK = true;
    }
    
    /* Stop before we reach the next zone then check for errors below. */
    position = text_line.find ("IZONE=",0);
    if (position != string::npos) {
      break;
    }
  }
  
  /* Close the mesh file. */
  mesh_file.close();
  
  /* Throw an error if any of the keywords was not found. */
  if (!foundNDIME) {
    SU2_MPI::Error(string("Could not find NDIME= keyword.") +
                   string(" \n Check the SU2 ASCII file format."),
                   CURRENT_FUNCTION);
  }
  if (!foundNPOIN) {
    SU2_MPI::Error(string("Could not find NPOIN= keyword.") +
                   string(" \n Check the SU2 ASCII file format."),
                   CURRENT_FUNCTION);
  }
  if (!foundNELEM) {
    SU2_MPI::Error(string("Could not find NELEM= keyword.") +
                   string(" \n Check the SU2 ASCII file format."),
                   CURRENT_FUNCTION);
  }
  if (!foundNMARK) {
    SU2_MPI::Error(string("Could not find NMARK= keyword.") +
                   string(" \n Check the SU2 ASCII file format."),
                   CURRENT_FUNCTION);
  }
  
}

void CSU2ASCIIMeshReaderBase::FastForwardToMyZone() {
  
  /*--- If more than one, fast-forward to my zone in the mesh file.  ---*/
  
  if (nZones > 1 && (config->GetMultizone_Mesh())) {
    
    string text_line;
    string::size_type position;
    while (getline (mesh_file,text_line)) {
      
      /*--- Search for the current domain ---*/
      position = text_line.find ("IZONE=",0);
      if (position != string::npos) {
        text_line.erase (0,6);
        unsigned short jZone = atoi(text_line.c_str());
        if (jZone == myZone+1) {
          return;
        }
      }
    }
  }
  
}
