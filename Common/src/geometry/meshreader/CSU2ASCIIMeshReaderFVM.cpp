/*!
 * \file CSU2ASCIIMeshReaderFVM.cpp
 * \brief Reads a native SU2 ASCII grid into linear partitions for the
 *        finite volume solver (FVM).
 * \author T. Economon
 * \version 7.5.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../../include/geometry/meshreader/CSU2ASCIIMeshReaderFVM.hpp"

CSU2ASCIIMeshReaderFVM::CSU2ASCIIMeshReaderFVM(CConfig *val_config,
                                               unsigned short val_iZone,
                                               unsigned short val_nZone)
: CSU2ASCIIMeshReaderBase(val_config, val_iZone, val_nZone) {
  
  actuator_disk  = (((config->GetnMarker_ActDiskInlet() != 0) ||
                     (config->GetnMarker_ActDiskOutlet() != 0)) &&
                    ((config->GetKind_SU2() == SU2_COMPONENT::SU2_CFD) ||
                     ((config->GetKind_SU2() == SU2_COMPONENT::SU2_DEF) &&
                      (config->GetActDisk_SU2_DEF()))));
  if (config->GetActDisk_DoubleSurface()) actuator_disk = false;
  ActDiskNewPoints = 0;
  Xloc = 0.0; Yloc = 0.0; Zloc = 0.0;

  /* Read the basic metadata and perform some basic error checks. */
  const auto try_single_pass = !actuator_disk;

  if (ReadMetadata(try_single_pass, val_config)) {
    /* The file contents were read together with the metadata. */
    return;
  }

  /* If the mesh contains an actuator disk as a single surface,
   we need to first split the surface into repeated points and update
   the connectivity for each element touching the surface. */
  if (actuator_disk)
    SplitActuatorDiskSurface();

  /* Read and store the points, interior elements, and surface elements.
   We store only the points and interior elements on our rank's linear
   partition, but the master stores the entire set of surface connectivity. */

  mesh_file.open(meshFilename);
  FastForwardToMyZone();

  for (auto section : SectionOrder) {
    switch (section) {
    case FileSection::ELEMENTS:
      ReadVolumeElementConnectivity();
      break;
    case FileSection::POINTS:
      ReadPointCoordinates();
      break;
    case FileSection::MARKERS:
      ReadSurfaceElementConnectivity();
      break;
    }
  }
  mesh_file.close();
}

void CSU2ASCIIMeshReaderFVM::SplitActuatorDiskSurface() {

  /*--- Actuator disk preprocesing ---*/

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

  mesh_file.open(meshFilename);
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
      break;
    }
  }

  mesh_file.close();

  /*--- Store the coordinates of the new points ---*/

  CoordXActDisk.resize(ActDiskNewPoints);
  CoordYActDisk.resize(ActDiskNewPoints);
  CoordZActDisk.resize(ActDiskNewPoints);

  /* Open the mesh file again to read the coordinates of the new points. */

  mesh_file.open(meshFilename);
  FastForwardToMyZone();

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

  mesh_file.open(meshFilename);
  FastForwardToMyZone();

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
      break;
    }
  }

  /* Lastly, increment the total number of points in order to add the
   new repeated points on the actuator disk. We also increment the
   number of markers by one. */
  numberOfGlobalPoints += ActDiskNewPoints;
  numberOfMarkers++;

  mesh_file.close();

}