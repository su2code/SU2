/*!
 * \file CSU2MeshFileWriter.cpp
 * \brief Filewriter class SU2 native mesh format.
 * \author T. Albring
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../../../Common/include/toolboxes/printing_toolbox.hpp"

const string CSU2MeshFileWriter::fileExt = ".su2";

CSU2MeshFileWriter::CSU2MeshFileWriter(CParallelDataSorter *valDataSorter,
                                       unsigned short valiZone, unsigned short valnZone) :
   CFileWriter(valDataSorter, fileExt), iZone(valiZone), nZone(valnZone) {}

void CSU2MeshFileWriter::WriteData(string val_filename) {

  ofstream output_file;

  /*--- We append the pre-defined suffix (extension) to the filename (prefix) ---*/
  val_filename.append(fileExt);

  /*--- Only the FIRST node writes the header (it does not matter if that is the master). ---*/

  if (rank == 0) {
    /*--- For multizone-cases this only works if the all zonal meshes are in one file.
          If the meshes are separate for each zone another solution has to be found. ---*/
    if (iZone==0) {
      output_file.open(val_filename);
    } else {
      output_file.open(val_filename, ios::app);
    }

    if (iZone==0 && nZone>1) {
      output_file << "NZONE= " << nZone << endl;
    }

    if (nZone > 1){
      output_file << "IZONE= " << iZone+1 << endl;
    }

    /*--- Write dimensions data. ---*/

    output_file << "NDIME= " << dataSorter->GetnDim() << endl;

    output_file << "NELEM= " << dataSorter->GetnElemGlobal() << endl;

    output_file.close();
  }

  unsigned long nElem = 0, offset = 0;

  for (int iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      output_file.open(val_filename, ios::app);

      for (auto iElem = 0ul; iElem < dataSorter->GetnElem(TRIANGLE); iElem++) {
        output_file << "5\t";
        for (auto iNode = 0u; iNode < N_POINTS_TRIANGLE; ++iNode)
          output_file << dataSorter->GetElemConnectivity(TRIANGLE, iElem, iNode) - 1 << "\t";
        output_file << nElem + offset << "\n"; nElem++;
      }
      for (auto iElem = 0ul; iElem < dataSorter->GetnElem(QUADRILATERAL); iElem++) {
        output_file << "9\t";
        for (auto iNode = 0u; iNode < N_POINTS_QUADRILATERAL; ++iNode)
          output_file << dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, iNode) - 1 << "\t";
        output_file << nElem + offset << "\n"; nElem++;
      }
      for (auto iElem = 0ul; iElem < dataSorter->GetnElem(TETRAHEDRON); iElem++) {
        output_file << "10\t";
        for (auto iNode = 0u; iNode < N_POINTS_TETRAHEDRON; ++iNode)
          output_file << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, iNode) - 1 << "\t";
        output_file << nElem + offset << "\n"; nElem++;
      }
      for (auto iElem = 0ul; iElem < dataSorter->GetnElem(HEXAHEDRON); iElem++) {
        output_file << "12\t";
        for (auto iNode = 0u; iNode < N_POINTS_HEXAHEDRON; ++iNode)
          output_file << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, iNode) - 1 << "\t";
        output_file << nElem + offset << "\n"; nElem++;
      }
      for (auto iElem = 0ul; iElem < dataSorter->GetnElem(PRISM); iElem++) {
        output_file << "13\t";
        for (auto iNode = 0u; iNode < N_POINTS_PRISM; ++iNode)
          output_file << dataSorter->GetElemConnectivity(PRISM, iElem, iNode) - 1 << "\t";
        output_file << nElem + offset << "\n"; nElem++;
      }

      for (auto iElem = 0ul; iElem < dataSorter->GetnElem(PYRAMID); iElem++) {
        output_file << "14\t";
        for (auto iNode = 0u; iNode < N_POINTS_PYRAMID; ++iNode)
          output_file << dataSorter->GetElemConnectivity(PYRAMID, iElem, iNode) - 1 << "\t";
        output_file << nElem + offset << "\n"; nElem++;
      }

      output_file.close();
    }

    /*--- Communicate offset, implies a barrier. ---*/
    SU2_MPI::Allreduce(&nElem, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  }

  /*--- Write the node coordinates. ---*/

  if (rank == 0) {
    output_file.open(val_filename, ios::app);
    output_file << "NPOIN= " << dataSorter->GetnPointsGlobal() << "\n";
    output_file.close();
  }

  unsigned long myPoint = 0; offset = 0;

  for (int iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      output_file.open(val_filename, ios::app);
      output_file.precision(15);

      for (auto iPoint = 0ul; iPoint < dataSorter->GetnPoints(); iPoint++) {

        /*--- Loop over the coordinates and write the values to file. ---*/

        for (auto iDim = 0u; iDim < dataSorter->GetnDim(); iDim++) {
          output_file << scientific << dataSorter->GetData(iDim, iPoint) << "\t";
        }

        /*--- Write global index. ---*/

        output_file << iPoint + offset << "\n";
        myPoint++;
      }

      output_file.close();
    }

    /*--- Communicate offset, implies a barrier. ---*/
    SU2_MPI::Allreduce(&myPoint, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  }

  if (rank == MASTER_NODE) {

    output_file.open(val_filename, ios::app);

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
      output_file << "NMARK= " << nMarker_ << endl;

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
        output_file << "MARKER_TAG= " << Marker_Tag << endl;
        output_file << "MARKER_ELEMS= " << nElem_Bound_<< endl;
        getline (input_file, text_line);

        text_line.erase(0,8);
        const auto SendTo = atoi(text_line.c_str());

        if (Marker_Tag == "SEND_RECEIVE") {
          output_file << "SEND_TO= " << SendTo << endl;
        }

        for (auto iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {

          getline(input_file, text_line);
          istringstream bound_line(text_line);

          unsigned short VTK_Type;
          bound_line >> VTK_Type;
          output_file << VTK_Type;
          unsigned long vnodes[4] = {0};

          switch (VTK_Type) {
          case LINE:
          case VERTEX:
            bound_line >> vnodes[0]; bound_line >> vnodes[1];
            output_file << "\t" << vnodes[0] << "\t" << vnodes[1] << "\n";
            break;
          case TRIANGLE:
            bound_line >> vnodes[0]; bound_line >> vnodes[1]; bound_line >> vnodes[2];
            output_file << "\t" << vnodes[0] << "\t" << vnodes[1] << "\t" << vnodes[2] << "\n";
            break;
          case QUADRILATERAL:
            bound_line >> vnodes[0]; bound_line >> vnodes[1]; bound_line >> vnodes[2]; bound_line >> vnodes[3];
            output_file << "\t" << vnodes[0] << "\t" << vnodes[1] << "\t" << vnodes[2] << "\t" << vnodes[3] << "\n";
            break;
          }
        }
      }
    }

    output_file.close();
  }

  SU2_MPI::Barrier(SU2_MPI::GetComm());
}
