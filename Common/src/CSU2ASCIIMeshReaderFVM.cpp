/*!
 * \file CSU2ASCIIMeshReaderFVM.cpp
 * \brief Reads a native SU2 ASCII grid into linear partitions for the
 *        finite volume solver (FVM).
 * \author T. Economon
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

#include "../include/toolboxes/CLinearPartitioner.hpp"
#include "../include/CSU2ASCIIMeshReaderFVM.hpp"

CSU2ASCIIMeshReaderFVM::CSU2ASCIIMeshReaderFVM(CConfig        *val_config,
                                               unsigned short val_iZone,
                                               unsigned short val_nZone)
: CMeshReaderFVM(val_config, val_iZone, val_nZone) {
  
  actuator_disk  = (((config->GetnMarker_ActDiskInlet() != 0) ||
                     (config->GetnMarker_ActDiskOutlet() != 0)) &&
                    ((config->GetKind_SU2() == SU2_CFD) ||
                     ((config->GetKind_SU2() == SU2_DEF) &&
                      (config->GetActDisk_SU2_DEF()))));
  if (config->GetActDisk_DoubleSurface()) actuator_disk = false;
  ActDiskNewPoints = 0;
  Xloc = 0.0; Yloc = 0.0; Zloc = 0.0;
  
  /* Store the current zone to be read and the total number of zones. */
  myZone = val_iZone;
  nZones = val_nZone;
  
  /* Store the mesh filename since we will open/close multiple times. */
  meshFilename = config->GetMesh_FileName();
  
  /* Read the basic metadata and perform some basic error checks. */
  ReadMetadata();
  
  /* If the mesh contains an actuator disk as a single surface,
   we need to first split the surface into repeated points and update
   the connectivity for each element touching the surface. */
  if (actuator_disk)
    SplitActuatorDiskSurface();
  
  /* Read and store the points, interior elements, and surface elements.
   We store only the points and interior elements on our rank's linear
   partition, but the master stores the entire set of surface connectivity. */
  ReadPointCoordinates();
  ReadVolumeElementConnectivity();
  ReadSurfaceElementConnectivity();
  
}

CSU2ASCIIMeshReaderFVM::~CSU2ASCIIMeshReaderFVM(void) { }

void CSU2ASCIIMeshReaderFVM::ReadMetadata() {
  
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

void CSU2ASCIIMeshReaderFVM::SplitActuatorDiskSurface() {
  
  /*--- Actuator disk preprocesing ---*/
  
  string Marker_Tag_Duplicate;
  bool InElem, Perimeter;
  unsigned long Counter = 0;
  Xloc = 0.0; Yloc = 0.0; Zloc = 0.0;
  unsigned long nElem_Bound_;
  
  vector<unsigned long long> EdgeBegin, EdgeEnd;
  
  unsigned long AuxEdge, iEdge, jEdge, nEdges, nPointVolume;
  unsigned long long FirstEdgeIndex, SecondEdgeIndex;
  
  vector<unsigned long> connectivity(N_POINTS_HEXAHEDRON,0);
  vector<unsigned long> ActDiskPoint_Front_Inv(numberOfGlobalPoints);
  vector<unsigned long> ActDiskPoint_Front;
  vector<unsigned long> VolumePoint;
  vector<unsigned long> PerimeterPoint;
  
  ActDiskNewPoints = 0;
  
  /* Note this routine is hard-coded to split for only one actuator disk
   boundary. Throw an error otherwise. */
  if (config->GetnMarker_ActDiskInlet() > 1) {
    SU2_MPI::Error(string("Current implementation can only split a single actuator disk.") +
                   string(" \n Remove disks or re-export your mesh with double surfaces (repeated points)."),
                   CURRENT_FUNCTION);
  }
  
  /*--- Open grid file ---*/
  
  mesh_file.open(meshFilename.c_str(), ios::in);
  FastForwardToMyZone();
  
  /*--- Read grid file with format SU2 ---*/
  
  string text_line;
  string::size_type position;
  while (getline (mesh_file, text_line)) {
    
    position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      
      for (unsigned short iMarker = 0 ; iMarker < numberOfMarkers; iMarker++) {
        
        getline (mesh_file, text_line);
        text_line.erase (0,11); string::size_type position;
        for (unsigned short iChar = 0; iChar < 20; iChar++) {
          position = text_line.find( " ", 0 );
          if (position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\r", 0 );
          if (position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\n", 0 );
          if (position != string::npos) text_line.erase (position,1);
        }
        string Marker_Tag = text_line.c_str();
        
        getline (mesh_file, text_line);
        text_line.erase (0,13);
        nElem_Bound_ = atoi(text_line.c_str());
        
        if (Marker_Tag != config->GetMarker_ActDiskInlet_TagBound(0)) {
          for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) { getline (mesh_file, text_line); }
        }
        else {
          
          if (rank == MASTER_NODE)
            cout << "Splitting the surface " << Marker_Tag << "( " << nElem_Bound_  << " boundary elements )." << endl;
          
          /*--- Create a list of edges ---*/
          
          for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {
            
            getline(mesh_file, text_line);
            
            unsigned short VTK_Type;
            istringstream bound_line(text_line);
            bound_line >> VTK_Type;
            
            switch(VTK_Type) {
              case LINE:
                bound_line >> connectivity[0];
                bound_line >> connectivity[1];
                
                EdgeBegin.push_back(connectivity[0]); EdgeEnd.push_back(connectivity[1]);
                break;
              case TRIANGLE:
                bound_line >> connectivity[0];
                bound_line >> connectivity[1];
                bound_line >> connectivity[2];
                EdgeBegin.push_back(connectivity[0]); EdgeEnd.push_back(connectivity[1]);
                EdgeBegin.push_back(connectivity[1]); EdgeEnd.push_back(connectivity[2]);
                EdgeBegin.push_back(connectivity[2]); EdgeEnd.push_back(connectivity[0]);
                break;
              case QUADRILATERAL:
                bound_line >> connectivity[0];
                bound_line >> connectivity[1];
                bound_line >> connectivity[2];
                bound_line >> connectivity[3];
                EdgeBegin.push_back(connectivity[0]); EdgeEnd.push_back(connectivity[1]);
                EdgeBegin.push_back(connectivity[1]); EdgeEnd.push_back(connectivity[2]);
                EdgeBegin.push_back(connectivity[2]); EdgeEnd.push_back(connectivity[3]);
                EdgeBegin.push_back(connectivity[3]); EdgeEnd.push_back(connectivity[0]);
                break;
            }
            
          }
          
          /*--- Set the total number of edges ---*/
          
          nEdges = EdgeBegin.size();
          
          /*--- Sort edges based on local point index, first index is always the largest ---*/
          
          for (iEdge = 0; iEdge <  nEdges; iEdge++) {
            if (EdgeEnd[iEdge] < EdgeBegin[iEdge]) {
              AuxEdge = EdgeEnd[iEdge]; EdgeEnd[iEdge] = EdgeBegin[iEdge]; EdgeBegin[iEdge] = AuxEdge;
            }
          }
          
          /*--- Bubble sort of the points based on the first index   ---*/
          
          for (iEdge = 0; iEdge < nEdges; iEdge++) {
            for (jEdge = iEdge+1; jEdge < nEdges; jEdge++) {
              
              FirstEdgeIndex = EdgeBegin[jEdge] << 31;
              FirstEdgeIndex += EdgeEnd[jEdge];
              
              SecondEdgeIndex = EdgeBegin[iEdge] << 31;
              SecondEdgeIndex += EdgeEnd[iEdge];
              
              if (FirstEdgeIndex <= SecondEdgeIndex) {
                AuxEdge = EdgeBegin[iEdge]; EdgeBegin[iEdge] = EdgeBegin[jEdge]; EdgeBegin[jEdge] = AuxEdge;
                AuxEdge = EdgeEnd[iEdge];  EdgeEnd[iEdge] = EdgeEnd[jEdge]; EdgeEnd[jEdge] = AuxEdge;
              }
            }
          }
          
          if (dimension == 3) {
            
            /*--- Check the begning of the list ---*/
            
            if (!((EdgeBegin[0] == EdgeBegin[1]) && (EdgeEnd[0] == EdgeEnd[1]))) {
              PerimeterPoint.push_back(EdgeBegin[0]);
              PerimeterPoint.push_back(EdgeEnd[0]);
            }
            
            for (iEdge = 1; iEdge < nEdges-1; iEdge++) {
              bool Check_1 = !((EdgeBegin[iEdge] == EdgeBegin[iEdge-1]) && (EdgeEnd[iEdge] == EdgeEnd[iEdge-1]));
              bool Check_2 = !((EdgeBegin[iEdge] == EdgeBegin[iEdge+1]) && (EdgeEnd[iEdge] == EdgeEnd[iEdge+1]));
              if ((Check_1 && Check_2)) {
                PerimeterPoint.push_back(EdgeBegin[iEdge]);
                PerimeterPoint.push_back(EdgeEnd[iEdge]);
              }
            }
            
            /*--- Check the  end of the list ---*/
            
            if (!((EdgeBegin[nEdges-1] == EdgeBegin[nEdges-2]) && (EdgeEnd[nEdges-1] == EdgeEnd[nEdges-2]))) {
              PerimeterPoint.push_back(EdgeBegin[nEdges-1]);
              PerimeterPoint.push_back(EdgeEnd[nEdges-1]);
            }
            
          } else {
            
            
            /*--- Create a list with all the points ---*/
            
            for (iEdge = 0; iEdge < nEdges; iEdge++) {
              ActDiskPoint_Front.push_back(EdgeBegin[iEdge]);
              ActDiskPoint_Front.push_back(EdgeEnd[iEdge]);
            }
            
            sort(ActDiskPoint_Front.begin(), ActDiskPoint_Front.end());
            vector<unsigned long>::iterator it = unique(ActDiskPoint_Front.begin(), ActDiskPoint_Front.end());
            ActDiskPoint_Front.resize(it - ActDiskPoint_Front.begin());
            
            /*--- Check the begning of the list ---*/
            
            if (!(ActDiskPoint_Front[0] == ActDiskPoint_Front[1]) ) { PerimeterPoint.push_back(ActDiskPoint_Front[0]); }
            
            for (unsigned long iPoint = 1; iPoint < ActDiskPoint_Front.size()-1; iPoint++) {
              bool Check_1 = !((ActDiskPoint_Front[iPoint] == ActDiskPoint_Front[iPoint-1]) );
              bool Check_2 = !((ActDiskPoint_Front[iPoint] == ActDiskPoint_Front[iPoint+1]) );
              if ((Check_1 && Check_2)) { PerimeterPoint.push_back(ActDiskPoint_Front[iEdge]); }
            }
            
            /*--- Check the  end of the list ---*/
            
            if (!((EdgeBegin[ActDiskPoint_Front.size()-1] == EdgeBegin[ActDiskPoint_Front.size()-2]) )) {
              PerimeterPoint.push_back(ActDiskPoint_Front[ActDiskPoint_Front.size()-1]);
            }
            
            ActDiskPoint_Front.clear();
            
          }
          
          vector<unsigned long>::iterator it;
          sort(PerimeterPoint.begin(), PerimeterPoint.end());
          it = unique(PerimeterPoint.begin(), PerimeterPoint.end());
          PerimeterPoint.resize(it - PerimeterPoint.begin());
          
          for (iEdge = 0; iEdge < nEdges; iEdge++) {
            
            Perimeter = false;
            for (unsigned long iPoint = 0; iPoint < PerimeterPoint.size(); iPoint++) {
              if (EdgeBegin[iEdge] == PerimeterPoint[iPoint]) {
                Perimeter = true; break;
              }
            }
            
            if (!Perimeter) ActDiskPoint_Front.push_back(EdgeBegin[iEdge]);
            
            Perimeter = false;
            for (unsigned long iPoint = 0; iPoint < PerimeterPoint.size(); iPoint++) {
              if (EdgeEnd[iEdge] == PerimeterPoint[iPoint]) {
                Perimeter = true; break;
              }
            }
            
            if (!Perimeter) ActDiskPoint_Front.push_back(EdgeEnd[iEdge]);
            
          }
          
          /*--- Sort, and remove repeated points from the disk list of points ---*/
          
          sort(ActDiskPoint_Front.begin(), ActDiskPoint_Front.end());
          it = unique(ActDiskPoint_Front.begin(), ActDiskPoint_Front.end());
          ActDiskPoint_Front.resize(it - ActDiskPoint_Front.begin());
          ActDiskNewPoints = ActDiskPoint_Front.size();
          
          if (rank == MASTER_NODE)
            cout << "Splitting the surface " << Marker_Tag << "( " << ActDiskPoint_Front.size()  << " internal points )." << endl;
          
          /*--- Create a map from original point to the new ones (back plane) ---*/
          
          ActDiskPoint_Back.resize(numberOfGlobalPoints);
          ActDisk_Bool.resize(numberOfGlobalPoints);
          
          for (unsigned long iPoint = 0; iPoint < numberOfGlobalPoints; iPoint++) {
            ActDisk_Bool[iPoint] = false;
            ActDiskPoint_Back[iPoint] = 0;
          }
          
          unsigned long kPoint = numberOfGlobalPoints;
          for (unsigned long iPoint = 0; iPoint < ActDiskPoint_Front.size(); iPoint++) {
            ActDiskPoint_Front_Inv[ActDiskPoint_Front[iPoint]] = iPoint;
            ActDisk_Bool[ActDiskPoint_Front[iPoint]] = true;
            ActDiskPoint_Back[ActDiskPoint_Front[iPoint]] = kPoint;
            kPoint++;
          }
          
        }
      }
    }
  }
  
  mesh_file.close();
  
  /*--- Store the coordinates of the new points ---*/
  
  CoordXActDisk.resize(ActDiskNewPoints);
  CoordYActDisk.resize(ActDiskNewPoints);
  CoordZActDisk.resize(ActDiskNewPoints);
  
  /* Open the mesh file again to read the coordinates of the new points. */
  
  mesh_file.open(meshFilename, ios::in);
  
  while (getline (mesh_file, text_line)) {
    
    position = text_line.find ("NPOIN=",0);
    if (position != string::npos) {
      for (unsigned long iPoint = 0; iPoint < numberOfGlobalPoints; iPoint++) {
        getline (mesh_file, text_line);
        istringstream point_line(text_line);
        
        su2double Coords[3] = {0.0,0.0,0.0};
        if (dimension == 2) {
          point_line >> Coords[0];
          point_line >> Coords[1];
        } else {
          point_line >> Coords[0];
          point_line >> Coords[1];
          point_line >> Coords[2];
        }
        
        /*--- Compute the CG of the actuator disk surface ---*/
        
        if (ActDisk_Bool[iPoint]) {
          CoordXActDisk[ActDiskPoint_Front_Inv[iPoint]] = Coords[0];
          CoordYActDisk[ActDiskPoint_Front_Inv[iPoint]] = Coords[1];
          Xloc += Coords[0]; Yloc += Coords[1];
          if (dimension == 3) {
            CoordZActDisk[ActDiskPoint_Front_Inv[iPoint]] = Coords[2];
            Zloc += Coords[2];
          }
          Counter++;
        }
        
      }
    }
    
    /*--- Locate and tag points that touch the actuator disk surface. ---*/
    
    position = text_line.find ("NELEM=",0);
    if (position != string::npos) {
      for (unsigned long iElem = 0; iElem < numberOfGlobalElements; iElem++) {
        
        getline(mesh_file, text_line);
        istringstream elem_line(text_line);
        
        unsigned short VTK_Type;
        elem_line >> VTK_Type;
        
        switch(VTK_Type) {
          case TRIANGLE:
            elem_line >> connectivity[0];
            elem_line >> connectivity[1];
            elem_line >> connectivity[2];
            InElem = false;
            for (unsigned long i = 0; i < (unsigned long)N_POINTS_TRIANGLE; i++) {
              if (ActDisk_Bool[connectivity[i]]) { InElem = true; break; } }
            if (InElem) {
              for (unsigned long i = 0; i < (unsigned long)N_POINTS_TRIANGLE; i++) {
                VolumePoint.push_back(connectivity[i]); } }
            break;
          case QUADRILATERAL:
            elem_line >> connectivity[0];
            elem_line >> connectivity[1];
            elem_line >> connectivity[2];
            elem_line >> connectivity[3];
            InElem = false;
            for (unsigned long i = 0; i < (unsigned long)N_POINTS_QUADRILATERAL; i++) {
              if (ActDisk_Bool[connectivity[i]]) { InElem = true; break; } }
            if (InElem) {
              for (unsigned long i = 0; i < (unsigned long)N_POINTS_QUADRILATERAL; i++) {
                VolumePoint.push_back(connectivity[i]); } }
            break;
          case TETRAHEDRON:
            elem_line >> connectivity[0];
            elem_line >> connectivity[1];
            elem_line >> connectivity[2];
            elem_line >> connectivity[3];
            InElem = false;
            for (unsigned long i = 0; i < (unsigned long)N_POINTS_TETRAHEDRON; i++) {
              if (ActDisk_Bool[connectivity[i]]) { InElem = true; break; }              }
            if (InElem) {
              for (unsigned long i = 0; i < (unsigned long)N_POINTS_TETRAHEDRON; i++) {
                VolumePoint.push_back(connectivity[i]); } }
            break;
          case HEXAHEDRON:
            elem_line >> connectivity[0];
            elem_line >> connectivity[1];
            elem_line >> connectivity[2];
            elem_line >> connectivity[3];
            elem_line >> connectivity[4];
            elem_line >> connectivity[5];
            elem_line >> connectivity[6];
            elem_line >> connectivity[7];
            InElem = false;
            for (unsigned long i = 0; i < (unsigned long)N_POINTS_HEXAHEDRON; i++) {
              if (ActDisk_Bool[connectivity[i]]) { InElem = true; break; } }
            if (InElem) {
              for (unsigned long i = 0; i < (unsigned long)N_POINTS_HEXAHEDRON; i++) {
                VolumePoint.push_back(connectivity[i]); } }
            break;
          case PRISM:
            elem_line >> connectivity[0];
            elem_line >> connectivity[1];
            elem_line >> connectivity[2];
            elem_line >> connectivity[3];
            elem_line >> connectivity[4];
            elem_line >> connectivity[5];
            InElem = false;
            for (unsigned long i = 0; i < (unsigned long)N_POINTS_PRISM; i++) {
              if (ActDisk_Bool[connectivity[i]]) { InElem = true; break; } }
            if (InElem) {
              for (unsigned long i = 0; i < (unsigned long)N_POINTS_PRISM; i++) {
                VolumePoint.push_back(connectivity[i]); } }
            break;
          case PYRAMID:
            elem_line >> connectivity[0];
            elem_line >> connectivity[1];
            elem_line >> connectivity[2];
            elem_line >> connectivity[3];
            elem_line >> connectivity[4];
            InElem = false;
            for (unsigned long i = 0; i < (unsigned long)N_POINTS_PYRAMID; i++) {
              if (ActDisk_Bool[connectivity[i]]) { InElem = true; break; } }
            if (InElem) {
              for (unsigned long i = 0; i < (unsigned long)N_POINTS_PYRAMID; i++) {
                VolumePoint.push_back(connectivity[i]); } }
            break;
        }
      }
    }
    
  }
  
  mesh_file.close();
  
  /*--- Compute the CG of the surface ---*/
  
  Xloc /= su2double(Counter);
  Yloc /= su2double(Counter);
  Zloc /= su2double(Counter);
  
  /*--- Sort and remove repeated points from the disk list of points. ---*/
  
  vector<unsigned long>::iterator it;
  sort(VolumePoint.begin(), VolumePoint.end());
  it = unique(VolumePoint.begin(), VolumePoint.end());
  VolumePoint.resize(it - VolumePoint.begin());
  nPointVolume = VolumePoint.size();
  
  /*--- Prepare some class data vectors for storage. ---*/
  
  CoordXVolumePoint.resize(nPointVolume);
  CoordYVolumePoint.resize(nPointVolume);
  CoordZVolumePoint.resize(nPointVolume);
  VolumePoint_Inv.resize(numberOfGlobalPoints);
  
  vector<bool> MapVolumePointBool(numberOfGlobalPoints);
  for (unsigned long iPoint = 0; iPoint < numberOfGlobalPoints; iPoint++) {
    MapVolumePointBool[iPoint] = false;
  }
  
  for (unsigned long iPoint = 0; iPoint < nPointVolume; iPoint++) {
    MapVolumePointBool[VolumePoint[iPoint]] = true;
    VolumePoint_Inv[VolumePoint[iPoint]]    = iPoint;
  }
  
  /*--- Store the coordinates of all the surface and volume
   points that touch the actuator disk ---*/
  
  mesh_file.open(meshFilename, ios::in);
  while (getline (mesh_file, text_line)) {
    position = text_line.find ("NPOIN=",0);
    if (position != string::npos) {
      for (unsigned long iPoint = 0; iPoint < numberOfGlobalPoints; iPoint++) {
        getline (mesh_file, text_line);
        istringstream point_line(text_line);
        su2double Coords[3] = {0.0,0.0,0.0};
        if (dimension == 2) {
          point_line >> Coords[0];
          point_line >> Coords[1];
        } else {
          point_line >> Coords[0];
          point_line >> Coords[1];
          point_line >> Coords[2];
        }
        
        if (MapVolumePointBool[iPoint]) {
          CoordXVolumePoint[VolumePoint_Inv[iPoint]] = Coords[0];
          CoordYVolumePoint[VolumePoint_Inv[iPoint]] = Coords[1];
          if (dimension == 3) {
            CoordZVolumePoint[VolumePoint_Inv[iPoint]] = Coords[2];
          }
        }
      }
    }
  }
  
  /* Lastly, increment the total number of points in order to add the
   new repeated points on the actuator disk. We also increment the
   number of markers by one. */
  numberOfGlobalPoints += ActDiskNewPoints;
  numberOfMarkers++;
  
  /* Close the mesh file. */
  mesh_file.close();
  
}

void CSU2ASCIIMeshReaderFVM::ReadPointCoordinates() {
  
  /* Get a partitioner to help with linear partitioning. */
  CLinearPartitioner pointPartitioner(numberOfGlobalPoints,0);
  
  /* Determine number of local points */
  for (unsigned long globalIndex=0; globalIndex < numberOfGlobalPoints; globalIndex++) {
    if ((int)pointPartitioner.GetRankContainingIndex(globalIndex) == rank) {
      numberOfLocalPoints++;
    }
  }
  
  /* Prepare our data structure for the point coordinates. */
  localPointCoordinates.resize(dimension);
  for (int k = 0; k < dimension; k++)
    localPointCoordinates[k].reserve(numberOfLocalPoints);
  
  /*--- Open the mesh file and jump to our zone. ---*/
  
  mesh_file.open(meshFilename, ios::in);
  FastForwardToMyZone();
  
  /*--- Read the point coordinates into our data structure. ---*/
  
  string text_line;
  string::size_type position;
  while (getline (mesh_file, text_line)) {
    
    position = text_line.find ("NPOIN=",0);
    if (position != string::npos) {
      
      unsigned long GlobalIndex = 0;
      while (GlobalIndex < numberOfGlobalPoints) {
        
        if (!actuator_disk) { getline(mesh_file, text_line); }
        else {
          if (GlobalIndex < numberOfGlobalPoints-ActDiskNewPoints) {
            getline(mesh_file, text_line);
          }
          else {
            /* This is a new actuator disk point, so we must construct a
             string with the new point's coordinates. */
            ostringstream strsX, strsY, strsZ;
            unsigned long BackActDisk_Index = GlobalIndex;
            unsigned long LocalIndex = BackActDisk_Index - (numberOfGlobalPoints-ActDiskNewPoints);
            strsX.precision(20); strsY.precision(20); strsZ.precision(20);
            su2double CoordX = CoordXActDisk[LocalIndex]; strsX << scientific << CoordX;
            su2double CoordY = CoordYActDisk[LocalIndex]; strsY << scientific << CoordY;
            su2double CoordZ = CoordZActDisk[LocalIndex]; strsZ << scientific << CoordZ;
            text_line = strsX.str() + "\t" + strsY.str() + "\t" + strsZ.str();
          }
        }
        
        /*--- We only read information for this node if it is owned by this
         rank based upon our initial linear partitioning. ---*/
        
        passivedouble Coords[3] = {0.0,0.0,0.0};
        if ((int)pointPartitioner.GetRankContainingIndex(GlobalIndex) == rank) {
          
          istringstream point_line(text_line);
          
          /* Store the coordinates more clearly. */
          if (dimension == 2) {
            point_line >> Coords[0];
            point_line >> Coords[1];
          } else {
            point_line >> Coords[0];
            point_line >> Coords[1];
            point_line >> Coords[2];
          }
          
          /* Load into the coordinate class data structure. */
          for (unsigned short iDim = 0; iDim < dimension; iDim++) {
            localPointCoordinates[iDim].push_back(Coords[iDim]);
            
          }
        }
        GlobalIndex++;
      }
    }
  }
  
  mesh_file.close();
  
}

void CSU2ASCIIMeshReaderFVM::ReadVolumeElementConnectivity() {
  
  /* Get a partitioner to help with linear partitioning. */
  CLinearPartitioner pointPartitioner(numberOfGlobalPoints,0);
  
  /* Loop over our analytically defined of elements and store only those
   that contain a node within our linear partition of points. */
  numberOfLocalElements  = 0;
  vector<unsigned long> connectivity(N_POINTS_HEXAHEDRON,0);
  
  /*--- Open the mesh file and jump to our zone. ---*/
  
  mesh_file.open(meshFilename, ios::in);
  FastForwardToMyZone();
  
  string text_line;
  string::size_type position;
  while (getline (mesh_file, text_line)) {
    
    /*--- Find the section containing the interior elements. ---*/
    
    position = text_line.find ("NELEM=",0);
    if (position != string::npos) {
      
      /*--- Loop over all the volumetric elements and store any element that
       contains at least one of an owned node for this rank (i.e., there will
       be element redundancy, since multiple ranks will store the same elems
       on the boundaries of the initial linear partitioning. ---*/
      
      numberOfLocalElements = 0;
      unsigned long GlobalIndex = 0;
      bool isOwned = false;
      while (GlobalIndex < numberOfGlobalElements) {
        getline(mesh_file, text_line);
        istringstream elem_line(text_line);
        
        /*--- Decide whether this rank needs each element. ---*/
        
        unsigned short VTK_Type;
        elem_line >> VTK_Type;
        
        switch(VTK_Type) {
            
          case TRIANGLE:
            
            /*--- Store the connectivity for this element more clearly. ---*/
            
            elem_line >> connectivity[0];
            elem_line >> connectivity[1];
            elem_line >> connectivity[2];
            
            /*--- Adjust for actuator disk splitting if necessary. ---*/
            
            if (actuator_disk) {
              for (unsigned short i = 0; i<N_POINTS_TRIANGLE; i++) {
                if (ActDisk_Bool[connectivity[i]]) {
                  
                  su2double Xcg = 0.0; unsigned long Counter = 0;
                  for (unsigned short j = 0; j<N_POINTS_TRIANGLE; j++) {
                    if (connectivity[j] < numberOfGlobalPoints-ActDiskNewPoints) {
                      Xcg += CoordXVolumePoint[VolumePoint_Inv[connectivity[j]]];
                      Counter++;
                    }
                  }
                  Xcg = Xcg / su2double(Counter);
                  
                  if (Counter != 0)  {
                    if (Xcg > Xloc) {
                      connectivity[i] = ActDiskPoint_Back[connectivity[i]];
                    }
                    else { connectivity[i] = connectivity[i]; }
                  }
                  
                }
              }
            }
            
            /* Check whether any of the points reside in our linear partition. */
            isOwned = false;
            for (unsigned short i = 0; i < N_POINTS_TRIANGLE; i++) {
              if ((int)pointPartitioner.GetRankContainingIndex(connectivity[i]) == rank) {
                isOwned = true;
              }
            }
            
            /* If so, we need to store the element locally. */
            if (isOwned) {
              localVolumeElementConnectivity.push_back(GlobalIndex);
              localVolumeElementConnectivity.push_back(VTK_Type);
              for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
                localVolumeElementConnectivity.push_back(connectivity[i]);
              }
              numberOfLocalElements++;
            }
            GlobalIndex++;
            break;
            
          case QUADRILATERAL:
            
            /*--- Load the connectivity for this element. ---*/
            
            elem_line >> connectivity[0];
            elem_line >> connectivity[1];
            elem_line >> connectivity[2];
            elem_line >> connectivity[3];
            
            if (actuator_disk) {
              for (unsigned short  i = 0; i<N_POINTS_QUADRILATERAL; i++) {
                if (ActDisk_Bool[connectivity[i]]) {
                  
                  su2double Xcg = 0.0; unsigned long Counter = 0;
                  for (unsigned short j = 0; j<N_POINTS_QUADRILATERAL; j++) {
                    if (connectivity[j] < numberOfGlobalPoints-ActDiskNewPoints) {
                      Xcg += CoordXVolumePoint[VolumePoint_Inv[connectivity[j]]];
                      Counter++;
                    }
                  }
                  Xcg = Xcg / su2double(Counter);
                  
                  
                  if (Counter != 0)  {
                    if (Xcg > Xloc) {
                      connectivity[i] = ActDiskPoint_Back[connectivity[i]];
                    }
                    else { connectivity[i] = connectivity[i]; }
                  }
                  
                }
              }
            }
            
            /* Check whether any of the points reside in our linear partition. */
            isOwned = false;
            for (unsigned short i = 0; i < N_POINTS_QUADRILATERAL; i++) {
              if ((int)pointPartitioner.GetRankContainingIndex(connectivity[i]) == rank) {
                isOwned = true;
              }
            }
            
            /* If so, we need to store the element locally. */
            if (isOwned) {
              localVolumeElementConnectivity.push_back(GlobalIndex);
              localVolumeElementConnectivity.push_back(VTK_Type);
              for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
                localVolumeElementConnectivity.push_back(connectivity[i]);
              }
              numberOfLocalElements++;
            }
            GlobalIndex++;
            break;
            
          case TETRAHEDRON:
            
            /*--- Load the connectivity for this element. ---*/
            
            elem_line >> connectivity[0];
            elem_line >> connectivity[1];
            elem_line >> connectivity[2];
            elem_line >> connectivity[3];
            
            if (actuator_disk) {
              for (unsigned short  i = 0; i<N_POINTS_TETRAHEDRON; i++) {
                if (ActDisk_Bool[connectivity[i]]) {
                  
                  su2double Xcg = 0.0; unsigned long Counter = 0;
                  for (unsigned short j = 0; j<N_POINTS_TETRAHEDRON; j++) {
                    if (connectivity[j] < numberOfGlobalPoints-ActDiskNewPoints) {
                      Xcg += CoordXVolumePoint[VolumePoint_Inv[connectivity[j]]];
                      Counter++;
                    }
                  }
                  Xcg = Xcg / su2double(Counter);
                  
                  
                  if (Counter != 0)  {
                    if (Xcg > Xloc) {
                      connectivity[i] = ActDiskPoint_Back[connectivity[i]];
                    }
                    else { connectivity[i] = connectivity[i]; }
                  }
                  
                }
              }
            }
            
            /* Check whether any of the points reside in our linear partition. */
            isOwned = false;
            for (unsigned short i = 0; i < N_POINTS_TETRAHEDRON; i++) {
              if ((int)pointPartitioner.GetRankContainingIndex(connectivity[i]) == rank) {
                isOwned = true;
              }
            }
            
            /* If so, we need to store the element locally. */
            if (isOwned) {
              localVolumeElementConnectivity.push_back(GlobalIndex);
              localVolumeElementConnectivity.push_back(VTK_Type);
              for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
                localVolumeElementConnectivity.push_back(connectivity[i]);
              }
              numberOfLocalElements++;
            }
            GlobalIndex++;
            break;
            
          case HEXAHEDRON:
            
            /*--- Load the connectivity for this element. ---*/
            
            elem_line >> connectivity[0];
            elem_line >> connectivity[1];
            elem_line >> connectivity[2];
            elem_line >> connectivity[3];
            elem_line >> connectivity[4];
            elem_line >> connectivity[5];
            elem_line >> connectivity[6];
            elem_line >> connectivity[7];
            
            if (actuator_disk) {
              for (unsigned short  i = 0; i<N_POINTS_HEXAHEDRON; i++) {
                if (ActDisk_Bool[connectivity[i]]) {
                  
                  su2double Xcg = 0.0; unsigned long Counter = 0;
                  for (unsigned short j = 0; j<N_POINTS_HEXAHEDRON; j++) {
                    if (connectivity[j] < numberOfGlobalPoints-ActDiskNewPoints) {
                      Xcg += CoordXVolumePoint[VolumePoint_Inv[connectivity[j]]];
                      Counter++;
                    }
                  }
                  Xcg = Xcg / su2double(Counter);
                  
                  if (Counter != 0)  {
                    if (Xcg > Xloc) { connectivity[i] = ActDiskPoint_Back[connectivity[i]]; }
                    else { connectivity[i] = connectivity[i]; }
                  }
                }
              }
            }
            
            /* Check whether any of the points reside in our linear partition. */
            isOwned = false;
            for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
              if ((int)pointPartitioner.GetRankContainingIndex(connectivity[i]) == rank) {
                isOwned = true;
              }
            }
            
            /* If so, we need to store the element locally. */
            if (isOwned) {
              localVolumeElementConnectivity.push_back(GlobalIndex);
              localVolumeElementConnectivity.push_back(VTK_Type);
              for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
                localVolumeElementConnectivity.push_back(connectivity[i]);
              }
              numberOfLocalElements++;
            }
            GlobalIndex++;
            break;
            
          case PRISM:
            
            /*--- Load the connectivity for this element. ---*/
            
            elem_line >> connectivity[0];
            elem_line >> connectivity[1];
            elem_line >> connectivity[2];
            elem_line >> connectivity[3];
            elem_line >> connectivity[4];
            elem_line >> connectivity[5];
            
            if (actuator_disk) {
              for (unsigned short i = 0; i<N_POINTS_PRISM; i++) {
                if (ActDisk_Bool[connectivity[i]]) {
                  
                  su2double Xcg = 0.0; unsigned long Counter = 0;
                  for (unsigned short j = 0; j<N_POINTS_PRISM; j++) {
                    if (connectivity[j] < numberOfGlobalPoints-ActDiskNewPoints) {
                      Xcg += CoordXVolumePoint[VolumePoint_Inv[connectivity[j]]];
                      Counter++;
                    }
                  }
                  Xcg = Xcg / su2double(Counter);
                  
                  if (Counter != 0)  {
                    if (Xcg > Xloc) { connectivity[i] = ActDiskPoint_Back[connectivity[i]]; }
                    else { connectivity[i] = connectivity[i]; }
                  }
                }
              }
            }
            
            /* Check whether any of the points reside in our linear partition. */
            isOwned = false;
            for (unsigned short i = 0; i < N_POINTS_PRISM; i++) {
              if ((int)pointPartitioner.GetRankContainingIndex(connectivity[i]) == rank) {
                isOwned = true;
              }
            }
            
            /* If so, we need to store the element locally. */
            if (isOwned) {
              localVolumeElementConnectivity.push_back(GlobalIndex);
              localVolumeElementConnectivity.push_back(VTK_Type);
              for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
                localVolumeElementConnectivity.push_back(connectivity[i]);
              }
              numberOfLocalElements++;
            }
            GlobalIndex++;
            break;
            
          case PYRAMID:
            
            /*--- Load the connectivity for this element. ---*/
            
            elem_line >> connectivity[0];
            elem_line >> connectivity[1];
            elem_line >> connectivity[2];
            elem_line >> connectivity[3];
            elem_line >> connectivity[4];
            
            if (actuator_disk) {
              for (unsigned short i = 0; i<N_POINTS_PYRAMID; i++) {
                if (ActDisk_Bool[connectivity[i]]) {
                  
                  su2double Xcg = 0.0; unsigned long Counter = 0;
                  for (unsigned short j = 0; j<N_POINTS_PYRAMID; j++) {
                    if (connectivity[j] < numberOfGlobalPoints-ActDiskNewPoints) {
                      Xcg += CoordXVolumePoint[VolumePoint_Inv[connectivity[j]]];
                      Counter++;
                    }
                  }
                  Xcg = Xcg / su2double(Counter);
                  
                  if (Counter != 0)  {
                    if (Xcg > Xloc) { connectivity[i] = ActDiskPoint_Back[connectivity[i]]; }
                    else { connectivity[i] = connectivity[i]; }
                  }
                }
              }
            }
            
            /* Check whether any of the points reside in our linear partition. */
            isOwned = false;
            for (unsigned short i = 0; i < N_POINTS_PYRAMID; i++) {
              if ((int)pointPartitioner.GetRankContainingIndex(connectivity[i]) == rank) {
                isOwned = true;
              }
            }
            
            /* If so, we need to store the element locally. */
            if (isOwned) {
              localVolumeElementConnectivity.push_back(GlobalIndex);
              localVolumeElementConnectivity.push_back(VTK_Type);
              for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
                localVolumeElementConnectivity.push_back(connectivity[i]);
              }
              numberOfLocalElements++;
            }
            GlobalIndex++;
            break;
        }
      }
      if (GlobalIndex == numberOfGlobalElements) break;
    }
  }
  
  mesh_file.close();
  
}

void CSU2ASCIIMeshReaderFVM::ReadSurfaceElementConnectivity() {
  
  /* We already read in the number of markers with the metadata. */
  surfaceElementConnectivity.resize(numberOfMarkers);
  markerNames.resize(numberOfMarkers);
  
  vector<unsigned long> connectivity(N_POINTS_HEXAHEDRON,0);
  
  /*--- In this routine, the boundary info is read by all ranks,
   however, the surface connectivity is still handled by the
   master node (and eventually distributed by the master as well). ---*/
  
  mesh_file.open(meshFilename, ios::in);
  FastForwardToMyZone();
  
  string text_line;
  string::size_type position;
  bool foundMarkers = false;
  while (getline (mesh_file, text_line)) {
    
    /*--- Jump to the section containing the markers. ---*/
    
    position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      
      foundMarkers = true;
      unsigned short iMarker = 0;
      while (iMarker < numberOfMarkers) {
        getline (mesh_file, text_line);
        text_line.erase (0,11);
        string::size_type position;
        
        for (unsigned short iChar = 0; iChar < 20; iChar++) {
          position = text_line.find( " ", 0 );
          if (position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\r", 0 );
          if (position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\n", 0 );
          if (position != string::npos) text_line.erase (position,1);
        }
        markerNames[iMarker] = text_line.c_str();

        bool duplicate = false;
        string Marker_Tag_Duplicate;
        if ((actuator_disk) &&
            (markerNames[iMarker] == config->GetMarker_ActDiskInlet_TagBound(0))) {
          duplicate = true;
          markerNames[iMarker+1] = config->GetMarker_ActDiskOutlet_TagBound(0);
        }
        
        /*--- Physical boundaries definition ---*/
        
        if (markerNames[iMarker]  != "SEND_RECEIVE") {
          
          getline (mesh_file, text_line);
          text_line.erase (0,13);
          unsigned long nElem_Bound = atoi(text_line.c_str());
          
          /*--- Allocate space for elements ---*/
          
          for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound; iElem_Bound++) {
            getline(mesh_file, text_line);
            istringstream bound_line(text_line);
            
            unsigned short VTK_Type;
            bound_line >> VTK_Type;
            switch(VTK_Type) {
              case LINE:
                
                if (dimension == 3) {
                  SU2_MPI::Error(string("Line boundary conditions are not possible for 3D calculations.") +
                                 string("Please check the SU2 ASCII mesh file."), CURRENT_FUNCTION);
                }
                
                bound_line >> connectivity[0];
                bound_line >> connectivity[1];
                
                surfaceElementConnectivity[iMarker].push_back(0);
                surfaceElementConnectivity[iMarker].push_back(VTK_Type);
                for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++)
                  surfaceElementConnectivity[iMarker].push_back(connectivity[i]);
                
                /* Duplicate the boundary element on an actuator disk if necessary. */
                
                if (duplicate) {
                  if (ActDisk_Bool[connectivity[0]]) {
                    connectivity[0] = ActDiskPoint_Back[connectivity[0]];
                  }
                  if (ActDisk_Bool[connectivity[1]]) {
                    connectivity[1] = ActDiskPoint_Back[connectivity[1]];
                  }
                  surfaceElementConnectivity[iMarker+1].push_back(0);
                  surfaceElementConnectivity[iMarker+1].push_back(VTK_Type);
                  for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++)
                    surfaceElementConnectivity[iMarker+1].push_back(connectivity[i]);
                }
                
                break;
                
              case TRIANGLE:
                
                bound_line >> connectivity[0];
                bound_line >> connectivity[1];
                bound_line >> connectivity[2];
                surfaceElementConnectivity[iMarker].push_back(0);
                surfaceElementConnectivity[iMarker].push_back(VTK_Type);
                for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++)
                  surfaceElementConnectivity[iMarker].push_back(connectivity[i]);
                
                if (duplicate) {
                  if (ActDisk_Bool[connectivity[0]]) {
                    connectivity[0] = ActDiskPoint_Back[connectivity[0]];
                  }
                  if (ActDisk_Bool[connectivity[1]]) {
                    connectivity[1] = ActDiskPoint_Back[connectivity[1]];
                  }
                  if (ActDisk_Bool[connectivity[2]]) {
                    connectivity[2] = ActDiskPoint_Back[connectivity[2]];
                  }
                  surfaceElementConnectivity[iMarker+1].push_back(0);
                  surfaceElementConnectivity[iMarker+1].push_back(VTK_Type);
                  for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++)
                    surfaceElementConnectivity[iMarker+1].push_back(connectivity[i]);
                  
                }
                
                break;
                
              case QUADRILATERAL:
                
                bound_line >> connectivity[0];
                bound_line >> connectivity[1];
                bound_line >> connectivity[2];
                bound_line >> connectivity[3];
                surfaceElementConnectivity[iMarker].push_back(0);
                surfaceElementConnectivity[iMarker].push_back(VTK_Type);
                for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++)
                  surfaceElementConnectivity[iMarker].push_back(connectivity[i]);
                
                if (duplicate) {
                  if (ActDisk_Bool[connectivity[0]]) {
                    connectivity[0] = ActDiskPoint_Back[connectivity[0]];
                    
                  }
                  if (ActDisk_Bool[connectivity[1]]) {
                    connectivity[1] = ActDiskPoint_Back[connectivity[1]];
                    
                  }
                  if (ActDisk_Bool[connectivity[2]]) {
                    connectivity[2] = ActDiskPoint_Back[connectivity[2]];
                    
                  }
                  if (ActDisk_Bool[connectivity[3]]) {
                    connectivity[3] = ActDiskPoint_Back[connectivity[3]];
                    
                  }
                  surfaceElementConnectivity[iMarker+1].push_back(0);
                  surfaceElementConnectivity[iMarker+1].push_back(VTK_Type);
                  for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++)
                    surfaceElementConnectivity[iMarker+1].push_back(connectivity[i]);
                  
                }
                
                break;
                
            }
          }
          
          /*--- Increment our counter if we stored a duplicate. ---*/
          
          iMarker++;
          if (duplicate) {
            iMarker++;
          }
          
        } else {
          /*--- Throw an error if we find deprecated references to SEND_RECEIVE
           boundaries in the mesh. ---*/
          SU2_MPI::Error(string("Mesh file contains deprecated SEND_RECEIVE marker!\n\n") +
                         string("Please remove any SEND_RECEIVE markers from the SU2 ASCII mesh."),
                         CURRENT_FUNCTION);
        }
        
        if (iMarker == numberOfMarkers) break;
      }
    }
    if (foundMarkers) break;
  }
  
  /*--- Final error check for deprecated periodic BC format. ---*/
  
  while (getline (mesh_file, text_line)) {
    
    /*--- Find any periodic transformation information. ---*/
    
    position = text_line.find ("NPERIODIC=",0);
    if (position != string::npos) {
      unsigned short nPeriodic;
      
      /*--- Read and store the number of transformations. ---*/
      text_line.erase (0,10); nPeriodic = atoi(text_line.c_str());
      if (rank == MASTER_NODE) {
        if (nPeriodic - 1 != 0)
          SU2_MPI::Error(string("Mesh file contains deprecated periodic format!\n\n") +
                         string("For SU2 v7.0.0 and later, preprocessing of periodic grids by SU2_MSH\n") +
                         string("is no longer necessary. Please use the original mesh file (prior to SU2_MSH)\n") +
                         string("with the same MARKER_PERIODIC definition in the configuration file.") , CURRENT_FUNCTION);
      }
    }
  }
  
  /*--- Close the input file ---*/
  
  mesh_file.close();
  
}

void CSU2ASCIIMeshReaderFVM::FastForwardToMyZone() {
  
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
