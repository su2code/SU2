/*!
 * \file CTecplotFileWriter.cpp
 * \brief Filewriter class for Tecplot ASCII format.
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

#include "../../../include/output/filewriter/CTecplotFileWriter.hpp"

const string CTecplotFileWriter::fileExt = ".dat";

CTecplotFileWriter::CTecplotFileWriter(CParallelDataSorter *valDataSorter,
                                       unsigned long valTimeIter, su2double valTimeStep) :
  CFileWriter(valDataSorter, fileExt), timeIter(valTimeIter), timeStep(valTimeStep){}

CTecplotFileWriter::~CTecplotFileWriter()= default;

void CTecplotFileWriter::WriteData(string val_filename){

  if (!dataSorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

  const vector<string> fieldNames = dataSorter->GetFieldNames();

  unsigned short iVar;

  unsigned long iPoint, iElem;

  /*--- Reduce the total number of each element. ---*/

  unsigned long nParallel_Line = dataSorter->GetnElem(LINE),
                nParallel_Tria = dataSorter->GetnElem(TRIANGLE),
                nParallel_Quad = dataSorter->GetnElem(QUADRILATERAL),
                nParallel_Tetr = dataSorter->GetnElem(TETRAHEDRON),
                nParallel_Hexa = dataSorter->GetnElem(HEXAHEDRON),
                nParallel_Pris = dataSorter->GetnElem(PRISM),
                nParallel_Pyra = dataSorter->GetnElem(PYRAMID);

  unsigned long nTot_Line = dataSorter->GetnElemGlobal(LINE),
                nTot_Tria = dataSorter->GetnElemGlobal(TRIANGLE),
                nTot_Quad = dataSorter->GetnElemGlobal(QUADRILATERAL),
                nTot_Tetr = dataSorter->GetnElemGlobal(TETRAHEDRON),
                nTot_Hexa = dataSorter->GetnElemGlobal(HEXAHEDRON),
                nTot_Pris = dataSorter->GetnElemGlobal(PRISM),
                nTot_Pyra = dataSorter->GetnElemGlobal(PYRAMID);

  /*--- Open Tecplot ASCII file and write the header. ---*/

  OpenMPIFile(val_filename);

  ostringstream header;
  header.precision(6);
  header << "TITLE = \"Visualization of the solution\"" << endl;

  header << "VARIABLES = ";
  for (iVar = 0; iVar < fieldNames.size()-1; iVar++) {
    header << "\"" << fieldNames[iVar] << "\",";
  }
  header << "\"" << fieldNames[fieldNames.size()-1] << "\"" << endl;

  header << "ZONE ";

  if (timeStep > 0.0){
    header << "STRANDID="<<SU2_TYPE::Int(timeIter+1)<<", SOLUTIONTIME="<< timeIter*timeStep <<", ";
  }

  header << "NODES= "<< dataSorter->GetnPointsGlobal() <<", ELEMENTS= "<< dataSorter->GetnElemGlobal();

  if (dataSorter->GetnDim() == 3){
    if ((nTot_Quad > 0 || nTot_Tria > 0) && (nTot_Hexa + nTot_Pris + nTot_Pyra + nTot_Tetr == 0)){
      header << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;
    }
    else {
      header <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
    }
  }
  else {
    if (nTot_Line > 0 && (nTot_Tria + nTot_Quad == 0)){
      header << ", DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
    }
    else{
      header << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
    }
  }

  WriteMPIString(header.str(), MASTER_NODE);

  /*--- Each rank formats the data of its own points and elements into a string, and all ranks then write
   their strings to the file at the same time, one after the other in rank order. ---*/

  ostringstream data;
  data.precision(6);
  data << scientific;

  /*--- Write surface and volumetric solution data. ---*/

  for (iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++) {
    for (iVar = 0; iVar < fieldNames.size(); iVar++)
      data << dataSorter->GetData(iVar, iPoint) << "\t";
    data << endl;
  }

  WriteMPIStringAll(data.str());

  data.str("");
  data.clear();


  /*--- Write connectivity data. ---*/

  {
    {

      for (iElem = 0; iElem < nParallel_Line; iElem++) {
        data << dataSorter->GetElemConnectivity(LINE, iElem, 0) << "\t";
        data << dataSorter->GetElemConnectivity(LINE, iElem, 1)<< "\n";
      }


      for (iElem = 0; iElem < nParallel_Tria; iElem++) {
        data << dataSorter->GetElemConnectivity(TRIANGLE, iElem, 0) << "\t";
        data << dataSorter->GetElemConnectivity(TRIANGLE, iElem, 1) << "\t";
        data << dataSorter->GetElemConnectivity(TRIANGLE, iElem, 2) << "\t";
        data << dataSorter->GetElemConnectivity(TRIANGLE, iElem, 2) << "\n";
      }

      for (iElem = 0; iElem < nParallel_Quad; iElem++) {
        data << dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, 0) << "\t";
        data << dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, 1) << "\t";
        data << dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, 2) << "\t";
        data << dataSorter->GetElemConnectivity(QUADRILATERAL, iElem, 3) << "\n";
      }

      for (iElem = 0; iElem < nParallel_Tetr; iElem++) {
        data << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 0) << "\t" << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 1) << "\t";
        data << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 2) << "\t" << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 2) << "\t";
        data << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 3) << "\t" << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 3) << "\t";
        data << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 3) << "\t" << dataSorter->GetElemConnectivity(TETRAHEDRON, iElem, 3) << "\n";
      }

      for (iElem = 0; iElem < nParallel_Hexa; iElem++) {
        data << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 0) << "\t" << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 1) << "\t";
        data << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 2) << "\t" << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 3) << "\t";
        data << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 4) << "\t" << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 5) << "\t";
        data << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 6) << "\t" << dataSorter->GetElemConnectivity(HEXAHEDRON, iElem, 7) << "\n";
      }

      for (iElem = 0; iElem < nParallel_Pris; iElem++) {
        data << dataSorter->GetElemConnectivity(PRISM, iElem, 0) << "\t" << dataSorter->GetElemConnectivity(PRISM, iElem, 1) << "\t";
        data << dataSorter->GetElemConnectivity(PRISM, iElem, 1) << "\t" << dataSorter->GetElemConnectivity(PRISM, iElem, 2) << "\t";
        data << dataSorter->GetElemConnectivity(PRISM, iElem, 3) << "\t" << dataSorter->GetElemConnectivity(PRISM, iElem, 4) << "\t";
        data << dataSorter->GetElemConnectivity(PRISM, iElem, 4) << "\t" << dataSorter->GetElemConnectivity(PRISM, iElem, 5) << "\n";
      }

      for (iElem = 0; iElem < nParallel_Pyra; iElem++) {
        data << dataSorter->GetElemConnectivity(PYRAMID, iElem, 0) << "\t" << dataSorter->GetElemConnectivity(PYRAMID, iElem, 1) << "\t";
        data << dataSorter->GetElemConnectivity(PYRAMID, iElem, 2) << "\t" << dataSorter->GetElemConnectivity(PYRAMID, iElem, 3) << "\t";
        data << dataSorter->GetElemConnectivity(PYRAMID, iElem, 4) << "\t" << dataSorter->GetElemConnectivity(PYRAMID, iElem, 4) << "\t";
        data << dataSorter->GetElemConnectivity(PYRAMID, iElem, 4) << "\t" << dataSorter->GetElemConnectivity(PYRAMID, iElem, 4) << "\n";
      }


    }
  }

  WriteMPIStringAll(data.str());

  CloseMPIFile();
}


