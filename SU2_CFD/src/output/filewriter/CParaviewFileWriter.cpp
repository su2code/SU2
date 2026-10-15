/*!
 * \file CParaviewFileWriter.cpp
 * \brief Filewriter class for Paraview ASCII format.
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

#include "../../../include/output/filewriter/CParaviewFileWriter.hpp"

const string CParaviewFileWriter::fileExt = ".vtk";

CParaviewFileWriter::CParaviewFileWriter(CParallelDataSorter *valDataSorter) :
  CFileWriter(valDataSorter, fileExt){}


CParaviewFileWriter::~CParaviewFileWriter()= default;

void CParaviewFileWriter::WriteData(string val_filename){

  if (!dataSorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  const unsigned short nDim = dataSorter->GetnDim();
  const vector<string> fieldNames = dataSorter->GetFieldNames();

  /*--- Each rank formats the data of its own points and elements into a string, and all ranks then write
   their strings to the file at the same time, one after the other in rank order. ---*/

  ostringstream data;
  data << scientific;

  auto resetData = [&]() {
    data.str("");
    data.clear();
    data << scientific;
  };

  OpenMPIFile(val_filename);

  /*--- Write the header. ---*/

  WriteMPIString("# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n", MASTER_NODE);
  WriteMPIString("POINTS " + to_string(dataSorter->GetnPointsGlobal()) + " double\n", MASTER_NODE);

  /*--- Write surface and volumetric point coordinates. ---*/

  for (unsigned long iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) data << dataSorter->GetData(iDim, iPoint) << "\t";
    if (nDim == 2) data << "0.0" << "\t";
  }
  WriteMPIStringAll(data.str());

  /*--- Write the connectivity, the number of nodes of each element is written before its nodes. ---*/

  const unsigned long nGlobal_Elem_Storage = dataSorter->GetnElemGlobal() + dataSorter->GetnConnGlobal();

  WriteMPIString("\nCELLS " + to_string(dataSorter->GetnElemGlobal()) + "\t" + to_string(nGlobal_Elem_Storage) + "\n",
                 MASTER_NODE);

  resetData();

  for (auto type : {LINE, TRIANGLE, QUADRILATERAL, TETRAHEDRON, HEXAHEDRON, PRISM, PYRAMID}) {
    const auto nPoints = nPointsOfElementType(type);
    for (unsigned long iElem = 0; iElem < dataSorter->GetnElem(type); iElem++) {
      data << nPoints << "\t";
      for (unsigned short iNode = 0; iNode < nPoints; iNode++)
        data << dataSorter->GetElemConnectivity(type, iElem, iNode) - 1 << "\t";
    }
  }
  WriteMPIStringAll(data.str());

  /*--- Write the type of each element. ---*/

  WriteMPIString("\nCELL_TYPES " + to_string(dataSorter->GetnElemGlobal()) + "\n", MASTER_NODE);

  resetData();

  for (auto type : {LINE, TRIANGLE, QUADRILATERAL, TETRAHEDRON, HEXAHEDRON, PRISM, PYRAMID}) {
    for (unsigned long iElem = 0; iElem < dataSorter->GetnElem(type); iElem++) data << type << "\t";
  }
  WriteMPIStringAll(data.str());

  /*--- Write the fields. ---*/

  WriteMPIString("\nPOINT_DATA " + to_string(dataSorter->GetnPointsGlobal()) + "\n", MASTER_NODE);

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
      /*--- We have found a vector, so skip the Y component. ---*/
      output_variable = false;
      VarCounter++;
    }
    found = fieldNames[iField].find("_z");
    if (found!=string::npos) {
      /*--- We have found a vector, so skip the Z component. ---*/
      output_variable = false;
      VarCounter++;
    }

    if (!output_variable) continue;

    resetData();

    if (isVector) {

      fieldname.erase(fieldname.end()-2,fieldname.end());

      WriteMPIString("\nVECTORS " + fieldname + " double\n", MASTER_NODE);

      for (unsigned long iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++) {
        data << dataSorter->GetData(VarCounter+0, iPoint) << "\t" << dataSorter->GetData(VarCounter+1, iPoint) << "\t";
        if (nDim == 3) data << dataSorter->GetData(VarCounter+2, iPoint) << "\t";
        if (nDim == 2) data << "0.0" << "\t";
      }

    } else {

      WriteMPIString("\nSCALARS " + fieldname + " double 1\nLOOKUP_TABLE default\n", MASTER_NODE);

      for (unsigned long iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++) {
        data << dataSorter->GetData(VarCounter, iPoint) << "\t";
      }
    }

    WriteMPIStringAll(data.str());
    VarCounter++;
  }

  CloseMPIFile();
}
