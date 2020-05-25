/*!
 * \file CSU2MeshFileWriter.cpp
 * \brief Filewriter class SU2 native mesh format.
 * \author T. Albring
 * \version 7.0.4 "Blackbird"
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

#include "../../../include/output/filewriter/CSU2MeshFileWriter.hpp"
#include "../../../../Common/include/toolboxes/printing_toolbox.hpp"

const string CSU2MeshFileWriter::fileExt = ".su2";

CSU2MeshFileWriter::CSU2MeshFileWriter(string valFileName, CParallelDataSorter *valDataSorter,
                                       unsigned short valiZone, unsigned short valnZone) :
   CFileWriter(std::move(valFileName), valDataSorter, fileExt), iZone(valiZone), nZone(valnZone) {}


CSU2MeshFileWriter::~CSU2MeshFileWriter(){

}


void CSU2MeshFileWriter::Write_Data(){

  unsigned long iElem, iPoint, iElem_Bound, nElem_Bound_, vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], offset, nElem;
  unsigned short iMarker, iDim, iChar, VTK_Type, nMarker_;
  short SendTo;
  ifstream input_file;
  string text_line, Marker_Tag, str;
  string::size_type position;
  int iProcessor;

  ofstream output_file;
  char cstr[MAX_STRING_SIZE], out_file[MAX_STRING_SIZE];

  strcpy (out_file, fileName.c_str());
  strcpy (cstr, out_file);

  if (rank == MASTER_NODE){

    /*--- For multizone-cases this only works if the all zonal meshes are in one file.
          If the meshes are separate for each zone another solution has to be found. ---*/
    if (iZone==0) {
      output_file.open(cstr, ios::out);
    } else {
      output_file.open(cstr, ios::out | ios::app);
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

  output_file.open(cstr, ios::out | ios::app);
  output_file.precision(15);
  nElem = 0;
  offset = 0;

  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iElem = 0; iElem < dataSorter->GetnElem(TRIANGLE); iElem++) {
        output_file << "5\t";
        output_file << dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 0) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 1) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(TRIANGLE, iElem, 2) - 1 << "\t";
        output_file << nElem + offset << "\n"; nElem++;
      }
      for (iElem = 0; iElem < dataSorter->GetnElem(QUADRILATERAL); iElem++) {
        output_file << "9\t";
        output_file << dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 0) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 1) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 2) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(QUADRILATERAL, iElem, 3) - 1 << "\t";
        output_file << nElem + offset << "\n"; nElem++;
      }
      for (iElem = 0; iElem < dataSorter->GetnElem(TETRAHEDRON); iElem++) {
        output_file << "10\t";
        output_file << dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 0) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 1) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 2) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(TETRAHEDRON, iElem, 3) - 1 << "\t";
        output_file << nElem + offset << "\n"; nElem++;
      }
      for (iElem = 0; iElem < dataSorter->GetnElem(HEXAHEDRON); iElem++) {
        output_file << "12\t";
        output_file << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 0) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 1) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 2) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 3) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 4) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 5) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 6) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(HEXAHEDRON, iElem, 7) - 1 << "\t";
        output_file << nElem + offset << "\n"; nElem++;
      }
      for (iElem = 0; iElem < dataSorter->GetnElem(PRISM); iElem++) {
        output_file << "13\t";
        output_file << dataSorter->GetElem_Connectivity(PRISM, iElem, 0) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(PRISM, iElem, 1) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(PRISM, iElem, 2) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(PRISM, iElem, 3) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(PRISM, iElem, 4) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(PRISM, iElem, 5) - 1 << "\t";
        output_file << nElem + offset << "\n"; nElem++;
      }

      for (iElem = 0; iElem < dataSorter->GetnElem(PYRAMID); iElem++) {
        output_file << "14\t";
        output_file << dataSorter->GetElem_Connectivity(PYRAMID, iElem, 0) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(PYRAMID, iElem, 1) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(PYRAMID, iElem, 2) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(PYRAMID, iElem, 3) - 1 << "\t";
        output_file << dataSorter->GetElem_Connectivity(PYRAMID, iElem, 4) - 1 << "\t";
        output_file << nElem + offset << "\n"; nElem++;
      }
    }
    output_file.flush();
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nElem, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }


  /*--- Write the node coordinates ---*/
  if (rank == MASTER_NODE){
    output_file << "NPOIN= " << dataSorter->GetnPointsGlobal();
    output_file << endl;
    output_file.flush();
  }


  unsigned long Global_Index, myPoint = 0;
  offset = 0;

  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iPoint = 0; iPoint <  dataSorter->GetnPoints(); iPoint++) {

        /*--- Global Index of the current point. (note outer loop over procs) ---*/

        Global_Index = iPoint + offset;

        /*--- Only write original domain points, i.e., exclude any periodic
         or halo nodes, even if they are output in the viz. files. ---*/

        /*--- Loop over the variables and write the values to file ---*/

        for (iDim = 0; iDim < dataSorter->GetnDim(); iDim++) {
          output_file << scientific << dataSorter->GetData(iDim, iPoint) << "\t";
        }

        /*--- Write global index. (note outer loop over procs) ---*/

        output_file << Global_Index << "\t";
        myPoint++;

        output_file << "\n";

      }
    }
    /*--- Flush the file and wait for all processors to arrive. ---*/
    output_file.flush();
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&myPoint, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  }

  output_file.close();

  if (rank == MASTER_NODE){

    output_file.open(cstr, ios::out | ios::app);

    /*--- Read the boundary information ---*/

    if (nZone == 1){
      str = "boundary";
    } else {
      str = "boundary_" + PrintingToolbox::to_string(iZone);
    }

    str += ".dat";

    input_file.open(str.c_str(), ios::out);

    if (!input_file.is_open()){
      SU2_MPI::Error(string("Cannot find ") + str, CURRENT_FUNCTION);
    }

    /*--- Read grid file with format SU2 ---*/

    while (getline (input_file, text_line)) {

      /*--- Write the physical boundaries ---*/

      position = text_line.find ("NMARK=",0);
      if (position != string::npos) {

        text_line.erase (0,6); nMarker_ = atoi(text_line.c_str());
        output_file << "NMARK= " << nMarker_ << endl;

        for (iMarker = 0 ; iMarker < nMarker_; iMarker++) {

          getline (input_file, text_line);
          text_line.erase (0,11);
          string::size_type position;
          for (iChar = 0; iChar < 20; iChar++) {
            position = text_line.find( " ", 0 );
            if (position != string::npos) text_line.erase (position,1);
            position = text_line.find( "\r", 0 );
            if (position != string::npos) text_line.erase (position,1);
            position = text_line.find( "\n", 0 );
            if (position != string::npos) text_line.erase (position,1);
          }
          Marker_Tag = text_line.c_str();

          /*--- Standart physical boundary ---*/

          getline (input_file, text_line);

          text_line.erase (0,13); nElem_Bound_ = atoi(text_line.c_str());
          output_file << "MARKER_TAG= " << Marker_Tag << endl;
          output_file << "MARKER_ELEMS= " << nElem_Bound_<< endl;
          getline (input_file, text_line);

          text_line.erase (0,8); SendTo = atoi(text_line.c_str());

          if (Marker_Tag == "SEND_RECEIVE"){
            output_file << "SEND_TO= " << SendTo << endl;
          }
          for (iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {

            getline(input_file, text_line);
            istringstream bound_line(text_line);

            bound_line >> VTK_Type;
            output_file << VTK_Type;

            switch(VTK_Type) {
            case LINE:
              bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
              output_file << "\t" << vnodes_edge[0] << "\t" << vnodes_edge[1] << endl;
              break;
            case TRIANGLE:
              bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
              output_file << "\t" << vnodes_triangle[0] << "\t" << vnodes_triangle[1] << "\t" << vnodes_triangle[2] << endl;
              break;
            case QUADRILATERAL:
              bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
              output_file << "\t" << vnodes_quad[0] << "\t" << vnodes_quad[1] << "\t" << vnodes_quad[2] << "\t" << vnodes_quad[3] << endl;
              break;
            case VERTEX:
              bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
              output_file << "\t" << vnodes_edge[0] <<  "\t" << vnodes_edge[1] <<endl;
              break;
            }
          }
        }
      }

    }

    input_file.close();

//    remove(str.c_str());

    output_file.close();
  }
}
