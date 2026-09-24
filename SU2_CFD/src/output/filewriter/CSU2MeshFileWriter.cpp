/*!
 * \file CSU2MeshFileWriter.cpp
 * \brief Filewriter class SU2 native mesh format.
 * \author T. Albring
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/output/filewriter/CSU2MeshFileWriter.hpp"

#include <numeric>
#include "../../../../Common/include/toolboxes/printing_toolbox.hpp"

const string CSU2MeshFileWriter::fileExt = ".su2";

CSU2MeshFileWriter::CSU2MeshFileWriter(CParallelDataSorter *valDataSorter,
                                       unsigned short valiZone, unsigned short valnZone) :
   CFileWriter(valDataSorter, fileExt), iZone(valiZone), nZone(valnZone) {}

void CSU2MeshFileWriter::WriteData(string val_filename) {

  /*--- For multizone cases all zones are written into one file, the zones after the first are appended. ---*/

  OpenMPIFile(val_filename, iZone != 0);

  /*--- Write the header. ---*/

  ostringstream header;

  if (iZone == 0 && nZone > 1) header << "NZONE= " << nZone << endl;
  if (nZone > 1) header << "IZONE= " << iZone + 1 << endl;

  header << "NDIME= " << dataSorter->GetnDim() << endl;
  header << "NELEM= " << dataSorter->GetnElemGlobal() << endl;

  WriteMPIString(header.str(), MASTER_NODE);

  /*--- Each rank formats the data of its own elements and points into a string, and all ranks then write their
   strings to the file at the same time, one after the other in rank order. The global index of an element or
   point is its local index plus the number of elements or points of the ranks before this one. ---*/

  auto offsetOfRank = [&](unsigned long localCount) {
    vector<unsigned long> counts(size, localCount);
    SU2_MPI::Allgather(&localCount, 1, MPI_UNSIGNED_LONG, counts.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
    return std::accumulate(counts.begin(), counts.begin() + rank, 0ul);
  };

  ostringstream data;

  /*--- Write the connectivity, the type of each element is written before its nodes. ---*/

  unsigned long nElem = 0;
  for (auto type : {TRIANGLE, QUADRILATERAL, TETRAHEDRON, HEXAHEDRON, PRISM, PYRAMID})
    nElem += dataSorter->GetnElem(type);

  unsigned long offset = offsetOfRank(nElem);

  nElem = 0;
  for (auto type : {TRIANGLE, QUADRILATERAL, TETRAHEDRON, HEXAHEDRON, PRISM, PYRAMID}) {
    const auto nPoints = nPointsOfElementType(type);
    for (auto iElem = 0ul; iElem < dataSorter->GetnElem(type); iElem++) {
      data << type << "\t";
      for (auto iNode = 0u; iNode < nPoints; ++iNode)
        data << dataSorter->GetElemConnectivity(type, iElem, iNode) - 1 << "\t";
      data << nElem + offset << "\n"; nElem++;
    }
  }

  WriteMPIStringAll(data.str());

  /*--- Write the node coordinates. ---*/

  WriteMPIString("NPOIN= " + to_string(dataSorter->GetnPointsGlobal()) + "\n", MASTER_NODE);

  offset = offsetOfRank(dataSorter->GetnPoints());

  data.str("");
  data.clear();
  data.precision(15);
  data << scientific;

  for (auto iPoint = 0ul; iPoint < dataSorter->GetnPoints(); iPoint++) {

    /*--- Loop over the coordinates and write the values to file. ---*/

    for (auto iDim = 0u; iDim < dataSorter->GetnDim(); iDim++) {
      data << dataSorter->GetData(iDim, iPoint) << "\t";
    }

    /*--- Write global index. ---*/

    data << iPoint + offset << "\n";
  }

  WriteMPIStringAll(data.str());

  /*--- The boundaries are copied from the file written by the mesh deformation, only the master node has them.
   This is the last thing written to the file. ---*/

  ostringstream boundaries;

  if (rank == MASTER_NODE) {

    /*--- Read the boundary information ---*/

    string str = "boundary";
    if (nZone > 1) str += "_" + PrintingToolbox::to_string(iZone);
    str += ".dat";

    ifstream input_file;
    input_file.open(str);

    if (!input_file.is_open()) {
      SU2_MPI::Error(string("Cannot find ") + str, CURRENT_FUNCTION);
    }

    /*--- Read grid file with format SU2 ---*/

    string text_line;
    while (getline(input_file, text_line)) {

      /*--- Write the physical boundaries ---*/

      auto position = text_line.find("NMARK=",0);

      if (position == string::npos) continue;

      text_line.erase(0,6);
      const auto nMarker_ = atoi(text_line.c_str());
      boundaries << "NMARK= " << nMarker_ << endl;

      for (auto iMarker = 0; iMarker < nMarker_; iMarker++) {

        getline(input_file, text_line);
        text_line.erase(0,11);
        for (int iChar = 0; iChar < 20; iChar++) {
          position = text_line.find(' ', 0);
          if (position != string::npos) text_line.erase(position,1);
          position = text_line.find('\r', 0);
          if (position != string::npos) text_line.erase(position,1);
          position = text_line.find('\n', 0);
          if (position != string::npos) text_line.erase(position,1);
        }
        string Marker_Tag = text_line;

        /*--- Standart physical boundary ---*/

        getline (input_file, text_line);

        text_line.erase(0,13);
        const auto nElem_Bound_ = atoi(text_line.c_str());
        boundaries << "MARKER_TAG= " << Marker_Tag << endl;
        boundaries << "MARKER_ELEMS= " << nElem_Bound_<< endl;
        getline (input_file, text_line);

        text_line.erase(0,8);
        const auto SendTo = atoi(text_line.c_str());

        if (Marker_Tag == "SEND_RECEIVE") {
          boundaries << "SEND_TO= " << SendTo << endl;
        }

        for (auto iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {

          getline(input_file, text_line);
          istringstream bound_line(text_line);

          unsigned short VTK_Type;
          bound_line >> VTK_Type;
          boundaries << VTK_Type;
          unsigned long vnodes[4] = {0};

          switch (VTK_Type) {
          case LINE:
          case VERTEX:
            bound_line >> vnodes[0]; bound_line >> vnodes[1];
            boundaries << "\t" << vnodes[0] << "\t" << vnodes[1] << "\n";
            break;
          case TRIANGLE:
            bound_line >> vnodes[0]; bound_line >> vnodes[1]; bound_line >> vnodes[2];
            boundaries << "\t" << vnodes[0] << "\t" << vnodes[1] << "\t" << vnodes[2] << "\n";
            break;
          case QUADRILATERAL:
            bound_line >> vnodes[0]; bound_line >> vnodes[1]; bound_line >> vnodes[2]; bound_line >> vnodes[3];
            boundaries << "\t" << vnodes[0] << "\t" << vnodes[1] << "\t" << vnodes[2] << "\t" << vnodes[3] << "\n";
            break;
          }
        }
      }
    }

  }

  /*--- Only the master node has this text, the other ranks write nothing. ---*/

  WriteMPIString(boundaries.str(), MASTER_NODE);

  CloseMPIFile();
}
