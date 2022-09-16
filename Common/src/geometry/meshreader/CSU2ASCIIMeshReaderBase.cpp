/*!
 * \file CSU2ASCIIMeshReaderBase.cpp
 * \brief Helper class for the reading of a native SU2 ASCII grid file.
 * \author T. Economon
 * \version 7.1.1 "Blackbird"
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

#include "../../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../../include/geometry/meshreader/CSU2ASCIIMeshReaderBase.hpp"

CSU2ASCIIMeshReaderBase::CSU2ASCIIMeshReaderBase(CConfig        *val_config,
                                                 unsigned short val_iZone,
                                                 unsigned short val_nZone)
: CMeshReaderFVM(val_config, val_iZone, val_nZone) {
  
  /* Store the current zone to be read and the total number of zones. */
  myZone = val_iZone;
  nZones = val_nZone;
  
  /* Store the mesh filename since we will open/close multiple times. */
  meshFilename = config->GetMesh_FileName();
  
}

CSU2ASCIIMeshReaderBase::~CSU2ASCIIMeshReaderBase(void) { }

bool CSU2ASCIIMeshReaderBase::ReadMetadata(const bool single_pass, CConfig *config) {

  const bool harmonic_balance = config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE;
  const bool multizone_file = config->GetMultizone_Mesh();

  /*--- Open grid file ---*/

  mesh_file.open(meshFilename);
  if (mesh_file.fail()) {
    SU2_MPI::Error("Error opening SU2 ASCII grid.\n"
                   "Check that the file exists.", CURRENT_FUNCTION);
  }

  /*--- If more than one, find the curent zone in the mesh file. ---*/

  string text_line;
  if ((nZones > 1 && multizone_file) || harmonic_balance) {
    if (harmonic_balance) {
      if (rank == MASTER_NODE) cout << "Reading time instance " << config->GetiInst()+1 << "." << endl;
    }
    else {
      bool foundZone = false;
      while (getline (mesh_file,text_line)) {
        /*--- Search for the current domain ---*/
        if (text_line.find ("IZONE=",0) != string::npos) {
          text_line.erase (0,6);
          unsigned short jZone = atoi(text_line.c_str());
          if (jZone == myZone+1) {
            if (rank == MASTER_NODE) cout << "Reading zone " << myZone << " from native SU2 ASCII mesh." << endl;
            foundZone = true;
            break;
          }
        }
      }
      if (!foundZone) {
        SU2_MPI::Error("Could not find the IZONE= keyword or the zone contents.\n"
                       "Check the SU2 ASCII file format.", CURRENT_FUNCTION);
      }
    }
  }

  /*--- Read the metadata: problem dimension, offsets for angle
   of attack and angle of sideslip, global points, global elements,
   and number of markers. Perform error checks as we go. ---*/

  bool foundNDIME = false, foundNPOIN = false;
  bool foundNELEM = false, foundNMARK = false;

  int current_section_idx = 0;
  bool single_pass_active = false;

  while (getline (mesh_file, text_line)) {

    /*--- Read the dimension of the problem ---*/

    if (!foundNDIME && text_line.find ("NDIME=",0) != string::npos) {
      text_line.erase (0,6);
      dimension = atoi(text_line.c_str());
      foundNDIME = true;
      continue;
    }

    /*--- The AoA and AoS offset values are optional. ---*/

    if (text_line.find ("AOA_OFFSET=",0) != string::npos) {
      text_line.erase (0,11);
      su2double AoA_Offset = atof(text_line.c_str());

      /*--- The offset is in deg ---*/
      const su2double AoA_Current = config->GetAoA() + AoA_Offset;
      config->SetAoA_Offset(AoA_Offset);
      config->SetAoA(AoA_Current);

      if (AoA_Offset != 0.0) {
        if (!config->GetDiscard_InFiles()) {
          cout << "WARNING: AoA in the config file (" << config->GetAoA() << " deg.) +\n";
          cout << "         AoA offset in mesh file (" << AoA_Offset << " deg.) = " << AoA_Current << " deg." << endl;
        }
        else {
          cout << "WARNING: Discarding the AoA offset in the mesh file." << endl;
        }
      }
      continue;
    }

    if (text_line.find ("AOS_OFFSET=",0) != string::npos) {
      text_line.erase (0,11);
      su2double AoS_Offset = atof(text_line.c_str());

      /*--- The offset is in deg ---*/
      const su2double AoS_Current = config->GetAoS() + AoS_Offset;
      config->SetAoS_Offset(AoS_Offset);
      config->SetAoS(AoS_Current);

      if (AoS_Offset != 0.0) {
        if (!config->GetDiscard_InFiles()) {
          cout << "WARNING: AoS in the config file (" << config->GetAoS() << " deg.) +\n";
          cout << "         AoS offset in mesh file (" << AoS_Offset << " deg.) = " << AoS_Current << " deg." << endl;
        }
        else {
          cout << "WARNING: Discarding the AoS offset in the mesh file." << endl;
        }
      }
      continue;
    }

    if (!foundNPOIN && text_line.find ("NPOIN=",0) != string::npos) {
      text_line.erase (0,6);
      numberOfGlobalPoints = atoi(text_line.c_str());

      /* If the points were found first, read them, otherwise just consume the lines. */
      if (single_pass && foundNDIME && current_section_idx == 0) {
        single_pass_active = true;
        ReadPointCoordinates(true);
      }
      else {
        for (auto iPoint = 0ul; iPoint < numberOfGlobalPoints; iPoint++)
          getline (mesh_file, text_line);
      }
      SectionOrder[current_section_idx++] = FileSection::POINTS;
      foundNPOIN = true;
      continue;
    }

    if (!foundNELEM && text_line.find ("NELEM=",0) != string::npos) {
      text_line.erase (0,6);
      numberOfGlobalElements = atoi(text_line.c_str());

      if (single_pass_active) {
        ReadVolumeElementConnectivity(true);
      }
      else {
        for (auto iElem = 0ul; iElem < numberOfGlobalElements; iElem++)
          getline (mesh_file, text_line);
      }
      SectionOrder[current_section_idx++] = FileSection::ELEMENTS;
      foundNELEM = true;
      continue;
    }

    if (!foundNMARK && text_line.find ("NMARK=",0) != string::npos) {
      text_line.erase (0,6);
      numberOfMarkers = atoi(text_line.c_str());

      if (current_section_idx != 2) {
        SU2_MPI::Error("Markers must be listed after points and elements in the SU2 mesh file.", CURRENT_FUNCTION);
      }

      if (single_pass_active) ReadSurfaceElementConnectivity(true);

      SectionOrder[current_section_idx++] = FileSection::MARKERS;
      foundNMARK = true;
      continue;
    }

    /* Stop before we reach the next zone then check for errors below. */
    if (text_line.find ("IZONE=",0) != string::npos) {
      break;
    }
  }

  mesh_file.close();

  /* Throw an error if any of the keywords was not found. */
  if (!foundNDIME) {
    SU2_MPI::Error("Could not find the keyword \"NDIME=\".\n"
                   "Check the SU2 ASCII file format.", CURRENT_FUNCTION);
  }
  if (!foundNPOIN) {
    SU2_MPI::Error("Could not find the keyword \"NPOIN=\".\n"
                   "Check the SU2 ASCII file format.", CURRENT_FUNCTION);
  }
  if (!foundNELEM) {
    SU2_MPI::Error("Could not find the keyword \"NELEM=\".\n"
                   "Check the SU2 ASCII file format.", CURRENT_FUNCTION);
  }
  if (!foundNMARK) {
    SU2_MPI::Error("Could not find the keyword \"NMARK=\".\n"
                   "Check the SU2 ASCII file format.", CURRENT_FUNCTION);
  }

  return single_pass_active;
}

void CSU2ASCIIMeshReaderBase::ReadPointCoordinates(const bool single_pass) {

  /* Get a partitioner to help with linear partitioning. */
  CLinearPartitioner pointPartitioner(numberOfGlobalPoints,0);

  /* Determine number of local points */
  numberOfLocalPoints = pointPartitioner.GetSizeOnRank(rank);

  /* Prepare our data structure for the point coordinates. */
  localPointCoordinates.resize(dimension);
  for (int k = 0; k < dimension; k++)
    localPointCoordinates[k].reserve(numberOfLocalPoints);

  /*--- Read the point coordinates into our data structure. ---*/

  while (true) {
    string text_line;
    if (!single_pass) {
      getline(mesh_file, text_line);
      if (text_line.find("NPOIN=",0) == string::npos) continue;
    }

    for (unsigned long GlobalIndex = 0; GlobalIndex < numberOfGlobalPoints; ++GlobalIndex) {

      if (!actuator_disk) {
        getline(mesh_file, text_line);
      }
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
      if (pointPartitioner.IndexBelongsToRank(GlobalIndex, rank)) {

        istringstream point_line(text_line);

        /* Store the coordinates more clearly. */
        point_line >> Coords[0];
        point_line >> Coords[1];
        if (dimension == 3) {
          point_line >> Coords[2];
        }

        /* Load into the coordinate class data structure. */
        for (unsigned short iDim = 0; iDim < dimension; iDim++) {
          localPointCoordinates[iDim].push_back(Coords[iDim]);
        }
      }
    }
    break;
  }
}

void CSU2ASCIIMeshReaderBase::ReadVolumeElementConnectivity(const bool single_pass) {

  /* Get a partitioner to help with linear partitioning. */
  CLinearPartitioner pointPartitioner(numberOfGlobalPoints,0);

  /* Loop over our analytically defined of elements and store only those
   that contain a node within our linear partition of points. */
  numberOfLocalElements  = 0;
  array<unsigned long, N_POINTS_HEXAHEDRON> connectivity{};

  while (true) {
    string text_line;
    if (!single_pass) {
      if (!getline(mesh_file, text_line)) break;
      if (text_line.find("NELEM=",0) == string::npos) continue;
    }

    /*--- Loop over all the volumetric elements and store any element that
     contains at least one of an owned node for this rank (i.e., there will
     be element redundancy, since multiple ranks will store the same elems
     on the boundaries of the initial linear partitioning. ---*/

    numberOfLocalElements = 0;

    for (unsigned long GlobalIndex = 0; GlobalIndex < numberOfGlobalElements; ++GlobalIndex) {
      getline(mesh_file, text_line);
      istringstream elem_line(text_line);

      /*--- Decide whether this rank needs each element. ---*/

      unsigned short VTK_Type;
      elem_line >> VTK_Type;

      const auto nPointsElem = nPointsOfElementType(VTK_Type);

      for (unsigned short  i = 0; i < nPointsElem; i++) {
        elem_line >> connectivity[i];
      }

      if (actuator_disk) {
        for (unsigned short  i = 0; i<nPointsElem; i++) {
          if (ActDisk_Bool[connectivity[i]]) {

            su2double Xcg = 0.0; unsigned long Counter = 0;
            for (unsigned short j = 0; j<nPointsElem; j++) {
              if (connectivity[j] < numberOfGlobalPoints-ActDiskNewPoints) {
                Xcg += CoordXVolumePoint[VolumePoint_Inv[connectivity[j]]];
                Counter++;
              }
            }
            Xcg = Xcg / su2double(Counter);

            if (Counter != 0 && Xcg > Xloc)  {
              connectivity[i] = ActDiskPoint_Back[connectivity[i]];
            }
          }
        }
      }

      /* Check whether any of the points reside in our linear partition. */
      bool isOwned = false;
      for (unsigned short i = 0; i < nPointsElem; i++) {
        if (pointPartitioner.IndexBelongsToRank(connectivity[i], rank)) {
          isOwned = true;
          break;
        }
      }

      /* If element is owned, we need to store it locally. */
      if (isOwned) {
        localVolumeElementConnectivity.push_back(GlobalIndex);
        localVolumeElementConnectivity.push_back(VTK_Type);
        /// TODO: Use a compressed format.
        for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
          localVolumeElementConnectivity.push_back(connectivity[i]);
        }
        numberOfLocalElements++;
      }
    }
    break;
  }
}

void CSU2ASCIIMeshReaderBase::ReadSurfaceElementConnectivity(const bool single_pass) {

  /* We already read in the number of markers with the metadata. */
  surfaceElementConnectivity.resize(numberOfMarkers);
  markerNames.resize(numberOfMarkers);

  array<unsigned long, N_POINTS_HEXAHEDRON> connectivity{};

  /*--- In this routine, the boundary info is read by all ranks,
   however, the surface connectivity is still handled by the
   master node (and eventually distributed by the master as well). ---*/

  while (true) {
    string text_line;
    if (!single_pass) {
      if (!getline(mesh_file, text_line)) break;
      if (text_line.find("NMARK=",0) == string::npos) continue;
    }

    for (unsigned short iMarker = 0; iMarker < numberOfMarkers; ++iMarker) {
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
      markerNames[iMarker] = text_line;

      bool duplicate = false;
      if ((actuator_disk) &&
          (markerNames[iMarker] == config->GetMarker_ActDiskInlet_TagBound(0))) {
        duplicate = true;
        markerNames[iMarker+1] = config->GetMarker_ActDiskOutlet_TagBound(0);
      }

      /*--- Physical boundaries definition ---*/

      if (markerNames[iMarker] == "SEND_RECEIVE") {
        /*--- Throw an error if we find deprecated references to SEND_RECEIVE
         boundaries in the mesh. ---*/
        SU2_MPI::Error("Mesh file contains deprecated SEND_RECEIVE marker!\n"
                       "Please remove any SEND_RECEIVE markers from the SU2 ASCII mesh.",
                       CURRENT_FUNCTION);
      }

      getline (mesh_file, text_line);
      text_line.erase (0,13);
      unsigned long nElem_Bound = atoi(text_line.c_str());

      /*--- Allocate space for elements ---*/

      for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound; iElem_Bound++) {
        getline(mesh_file, text_line);
        istringstream bound_line(text_line);

        unsigned short VTK_Type;
        bound_line >> VTK_Type;

        const auto nPointsElem = nPointsOfElementType(VTK_Type);

        if (dimension == 3 && VTK_Type == LINE) {
          SU2_MPI::Error("Line boundary conditions are not possible for 3D calculations.\n"
                         "Please check the SU2 ASCII mesh file.", CURRENT_FUNCTION);
        }

        for (unsigned short i = 0; i < nPointsElem; i++) {
          bound_line >> connectivity[i];
        }

        surfaceElementConnectivity[iMarker].push_back(0);
        surfaceElementConnectivity[iMarker].push_back(VTK_Type);
        for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
          surfaceElementConnectivity[iMarker].push_back(connectivity[i]);
        }

        if (duplicate) {
          for (unsigned short i = 0; i < nPointsElem; i++) {
            if (ActDisk_Bool[connectivity[i]]) {
              connectivity[i] = ActDiskPoint_Back[connectivity[i]];
            }
          }
          surfaceElementConnectivity[iMarker + 1].push_back(0);
          surfaceElementConnectivity[iMarker + 1].push_back(VTK_Type);
          for (unsigned short i = 0; i < N_POINTS_HEXAHEDRON; i++) {
            surfaceElementConnectivity[iMarker + 1].push_back(connectivity[i]);
          }
        }
      }
      /*--- Increment the counter an extra time if we stored a duplicate. ---*/
      iMarker += duplicate;
    }
    break;
  }

  if (rank != MASTER_NODE) return;

  /*--- Final error check for deprecated periodic BC format. ---*/

  string text_line;
  while (getline (mesh_file, text_line)) {

    /*--- Find any periodic transformation information. ---*/

    if (text_line.find ("NPERIODIC=",0) != string::npos) {

      /*--- Read and store the number of transformations. ---*/
      text_line.erase(0,10);
      unsigned short nPeriodic = atoi(text_line.c_str());
      if (nPeriodic - 1 != 0) {
        SU2_MPI::Error("Mesh file contains deprecated periodic format!\n\n"
                       "For SU2 v7.0.0 and later, preprocessing of periodic grids by SU2_MSH\n"
                       "is no longer necessary. Please use the original mesh file (prior to SU2_MSH)\n"
                       "with the same MARKER_PERIODIC definition in the configuration file.", CURRENT_FUNCTION);
      }
    }

    /*--- Stop before we reach the next zone. ---*/
    if (text_line.find ("IZONE=",0) != string::npos) break;
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
