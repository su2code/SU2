/*!
 * \file CParaviewFileWriter.cpp
 * \brief Filewriter class for Paraview ASCII format.
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

#include "../../../include/output/filewriter/CParaviewFileWriter.hpp"

const string CParaviewFileWriter::fileExt = ".vtk";

CParaviewFileWriter::CParaviewFileWriter(CParallelDataSorter *valDataSorter) :
  CFileWriter(valDataSorter, fileExt){}


CParaviewFileWriter::~CParaviewFileWriter()= default;

void CParaviewFileWriter::WriteData(string val_filename){

  /*--- We append the pre-defined suffix (extension) to the filename (prefix) ---*/
  val_filename.append(fileExt);

  if (!dataSorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  unsigned short iDim = 0, nDim = dataSorter->GetnDim();

  unsigned long iPoint, iElem;

  unsigned long nGlobal_Elem_Storage;

  ofstream Paraview_File;

  int iProcessor;

  const vector<string> fieldNames = dataSorter->GetFieldNames();

  /*--- Set a timer for the file writing. ---*/

  startTime = SU2_MPI::Wtime();

  /*--- Open Paraview ASCII file and write the header. ---*/

  if (rank == MASTER_NODE) {

    Paraview_File.open(val_filename.c_str(), ios::out);
    Paraview_File.precision(6);
    Paraview_File << "# vtk DataFile Version 3.0\n";
    Paraview_File << "vtk output\n";
    Paraview_File << "ASCII\n";
    Paraview_File << "DATASET UNSTRUCTURED_GRID\n";

    /*--- Write the header ---*/
    Paraview_File << "POINTS "<< dataSorter->GetnPointsGlobal() <<" double\n";

  }

  Paraview_File.close();

#ifdef HAVE_MPI
  SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif

  /*--- Each processor opens the file. ---*/

  Paraview_File.open(val_filename.c_str(), ios::out | ios::app);

  /*--- Write surface and volumetric point coordinates. ---*/

  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {

      /*--- Write the node data from this proc ---*/


      for (iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++) {
        for (iDim = 0; iDim < nDim; iDim++)
          Paraview_File << scientific << dataSorter->GetData(iDim, iPoint) << "\t";
        if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
      }
    }

    Paraview_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif
  }

  /*--- Reduce the total number of each element. ---*/

  unsigned long nParallel_Line = dataSorter->GetnElem(LINE),
                nParallel_Tria = dataSorter->GetnElem(TRIANGLE),
                nParallel_Quad = dataSorter->GetnElem(QUADRILATERAL),
                nParallel_Tetr = dataSorter->GetnElem(TETRAHEDRON),
                nParallel_Hexa = dataSorter->GetnElem(HEXAHEDRON),
                nParallel_Pris = dataSorter->GetnElem(PRISM),
                nParallel_Pyra = dataSorter->GetnElem(PYRAMID);

  if (rank == MASTER_NODE) {

    /*--- Write the header ---*/
    nGlobal_Elem_Storage = dataSorter->GetnElemGlobal() + dataSorter->GetnConnGlobal();

    Paraview_File << "\nCELLS " << dataSorter->GetnElemGlobal() << "\t" << nGlobal_Elem_Storage << "\n";

  }

  Paraview_File.flush();
#ifdef HAVE_MPI
  SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif

  /*--- Write connectivity data. ---*/

  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {


      for (iElem = 0; iElem < nParallel_Line; iElem++) {
        Paraview_File << N_POINTS_LINE << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(LINE, iElem, 0)-1 << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(LINE, iElem, 1)-1 << "\t";
      }

      for (iElem = 0; iElem < nParallel_Tria; iElem++) {
        Paraview_File << N_POINTS_TRIANGLE << "\t";
        Paraview_File <<  dataSorter->GetElemConnectivity(TRIANGLE, iElem, 0)-1 << "\t";
        Paraview_File <<  dataSorter->GetElemConnectivity(TRIANGLE, iElem, 1)-1 << "\t";
        Paraview_File <<  dataSorter->GetElemConnectivity(TRIANGLE, iElem, 2)-1 << "\t";
      }

      for (iElem = 0; iElem < nParallel_Quad; iElem++) {
        Paraview_File << N_POINTS_QUADRILATERAL << "\t";
        Paraview_File <<  dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, 0)-1 << "\t";
        Paraview_File <<  dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, 1)-1 << "\t";
        Paraview_File <<  dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, 2)-1 << "\t";
        Paraview_File <<  dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, 3)-1 << "\t";
      }


      for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
        Paraview_File << N_POINTS_TETRAHEDRON << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 0)-1 << "\t"
                      << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 1)-1 << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 2)-1 << "\t"
                      << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 3)-1 << "\t";
      }

      for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
        Paraview_File << N_POINTS_HEXAHEDRON << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 0)-1 << "\t"
                      << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 1)-1 << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 2)-1 << "\t"
                      << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 3)-1 << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 4)-1 << "\t"
                      << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 5)-1 << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 6)-1 << "\t"
                      << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 7)-1 << "\t";
      }

      for (iElem = 0; iElem < nParallel_Pris; iElem++) {
        Paraview_File << N_POINTS_PRISM << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(PRISM, iElem, 0)-1 << "\t"
                      << dataSorter->GetElemConnectivity(PRISM, iElem, 1)-1 << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(PRISM, iElem, 2)-1 << "\t"
                      << dataSorter->GetElemConnectivity(PRISM, iElem, 3)-1 << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(PRISM, iElem, 4)-1 << "\t"
                      << dataSorter->GetElemConnectivity(PRISM, iElem, 5)-1 << "\t";
      }

      for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
        Paraview_File << N_POINTS_PYRAMID << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(PYRAMID, iElem, 0)-1 << "\t"
                      << dataSorter->GetElemConnectivity(PYRAMID, iElem, 1)-1 << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(PYRAMID, iElem, 2)-1 << "\t"
                      << dataSorter->GetElemConnectivity(PYRAMID, iElem, 3)-1 << "\t";
        Paraview_File << dataSorter->GetElemConnectivity(PYRAMID, iElem, 4)-1 << "\t";
      }

    }    Paraview_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif
  }

  if (rank == MASTER_NODE) {

    /*--- Write the header ---*/
    Paraview_File << "\nCELL_TYPES " << dataSorter->GetnElemGlobal() << "\n";

  }

  Paraview_File.flush();
#ifdef HAVE_MPI
  SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif

  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iElem = 0; iElem < nParallel_Line; iElem++) Paraview_File << "3\t";
      for (iElem = 0; iElem < nParallel_Tria; iElem++) Paraview_File << "5\t";
      for (iElem = 0; iElem < nParallel_Quad; iElem++) Paraview_File << "9\t";
      for (iElem = 0; iElem < nParallel_Tetr; iElem++) Paraview_File << "10\t";
      for (iElem = 0; iElem < nParallel_Hexa; iElem++) Paraview_File << "12\t";
      for (iElem = 0; iElem < nParallel_Pris; iElem++) Paraview_File << "13\t";
      for (iElem = 0; iElem < nParallel_Pyra; iElem++) Paraview_File << "14\t";
    }
    Paraview_File.flush();
#ifdef HAVE_MPI
    SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif
  }

  if (rank == MASTER_NODE) {
    /*--- Write the header ---*/
    Paraview_File << "\nPOINT_DATA "<< dataSorter->GetnPointsGlobal() <<"\n";

  }

  Paraview_File.flush();
#ifdef HAVE_MPI
  SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif

  unsigned short varStart = 2;
  if (nDim == 3) varStart++;

  /*--- Need to adjust container location to avoid PointID tag and coords. ---*/
  unsigned short VarCounter = varStart;

  for (unsigned short iField = varStart; iField < fieldNames.size(); iField++) {

    string fieldname = fieldNames[iField];

    fieldname.erase(remove(fieldname.begin(), fieldname.end(), '"'), fieldname.end());

    bool output_variable = true, isVector = false;
    size_t found = fieldNames[iField].find("_x");
    if (found!=string::npos) {
      output_variable = true;
      isVector = true;
    }
    found = fieldNames[iField].find("_y");
    if (found!=string::npos) {
      output_variable = false;
      //skip
      Paraview_File.flush();
#ifdef HAVE_MPI
      SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif
      VarCounter++;
    }
    found = fieldNames[iField].find("_z");
    if (found!=string::npos) {
      output_variable = false;
      //skip
      Paraview_File.flush();
#ifdef HAVE_MPI
      SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif
      VarCounter++;
    }

    if (output_variable && isVector) {

      fieldname.erase(fieldname.end()-2,fieldname.end());

      if (rank == MASTER_NODE) {
        Paraview_File << "\nVECTORS " << fieldname << " double\n";
      }

      Paraview_File.flush();
#ifdef HAVE_MPI
      SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif

      /*--- Write surface and volumetric point coordinates. ---*/

      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        if (rank == iProcessor) {

          /*--- Write the node data from this proc ---*/

          for (iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++) {
            Paraview_File << scientific << dataSorter->GetData(VarCounter+0, iPoint) << "\t" << dataSorter->GetData(VarCounter+1, iPoint) << "\t";
            if (nDim == 3) Paraview_File << scientific << dataSorter->GetData(VarCounter+2, iPoint) << "\t";
            if (nDim == 2) Paraview_File << scientific << "0.0" << "\t";
          }
        }

        Paraview_File.flush();
#ifdef HAVE_MPI
        SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif
      }

      VarCounter++;

    } else if (output_variable) {

      if (rank == MASTER_NODE) {

        Paraview_File << "\nSCALARS " << fieldname << " double 1\n";
        Paraview_File << "LOOKUP_TABLE default\n";
      }

      Paraview_File.flush();
#ifdef HAVE_MPI
      SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif

      /*--- Write surface and volumetric point coordinates. ---*/

      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        if (rank == iProcessor) {

          /*--- Write the node data from this proc ---*/

          for (iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++) {
            Paraview_File << scientific << dataSorter->GetData(VarCounter, iPoint) << "\t";
          }

        }
        Paraview_File.flush();
#ifdef HAVE_MPI
        SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif
      }

      VarCounter++;
    }

  }

  Paraview_File.close();


  /*--- Compute and store the write time. ---*/

  stopTime = SU2_MPI::Wtime();

  usedTime = stopTime-startTime;

  fileSize = DetermineFilesize(val_filename);

  /*--- Compute and store the bandwidth ---*/

  bandwidth = fileSize/(1.0e6)/usedTime;
}

