/*!
 * \file CParaviewVTMFileWriter.cpp
 * \brief Filewriter class for Paraview binary format.
 * \author T. Albring
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/output/filewriter/CParaviewVTMFileWriter.hpp"
#include "../../../include/output/filewriter/CParaviewXMLFileWriter.hpp"
#include "../../../../Common/include/toolboxes/printing_toolbox.hpp"

const string CParaviewVTMFileWriter::fileExt = ".vtm";

CParaviewVTMFileWriter::CParaviewVTMFileWriter(string fileName, string folderName, su2double time, unsigned short iZone, unsigned short nZone)
  : CFileWriter(std::move(fileName), fileExt),
    folderName(std::move(folderName)){
  
  if (rank == MASTER_NODE)
#if defined(_WIN32)
    _mkdir(this->folderName.c_str());
#else 
    mkdir(this->folderName.c_str(), 0777); // notice that 777 is different than 0777
#endif
  
  this->iZone = iZone;
  this->nZone = nZone;
  
  curTime = time;
}


CParaviewVTMFileWriter::~CParaviewVTMFileWriter(){

}

void CParaviewVTMFileWriter::Write_Data(){
  
  /*--- If we are in the first zone, create new file and write the file header, 
   * otherwise append to already existing file ---*/
  
  if (rank == MASTER_NODE){
    ofstream multiBlockFile;
    if (iZone == 0)
      multiBlockFile.open (fileName.c_str());
    else 
      multiBlockFile.open(fileName.c_str(), ios::app);
    
    if (iZone == 0){
      multiBlockFile << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">" << endl;
      multiBlockFile << "<vtkMultiBlockDataSet>" << endl;
    }
    
    /*--- Write all blocks that have been added ---*/
    
    multiBlockFile << output.str();
    
    /*--- If we are in the last zone, write the additional data and close all blocks ---*/
    
    if (iZone == nZone-1){
      multiBlockFile << "</vtkMultiBlockDataSet>" << endl;
      multiBlockFile << "<FieldData>" << endl;
      multiBlockFile << "<DataArray type='Float32' Name='TimeValue'>" << endl;
      multiBlockFile << curTime << endl;
      multiBlockFile << "</DataArray>" << endl;
      multiBlockFile << "</FieldData>" << endl;
      multiBlockFile << "</VTKFile>" << endl;
    }
    multiBlockFile.close();
  }
  
}

void CParaviewVTMFileWriter::AddDataset(string name, string file, vector<string> fieldNames,
                                        unsigned short nDim, CParallelDataSorter* dataSorter){
  
  /*--- Construct the full file name incl. folder ---*/
  
  string fullFilename = folderName + "/" + file;
  
  /*--- Create an XML writer and dump data into file ---*/
  
  CParaviewXMLFileWriter* XMLWriter = new CParaviewXMLFileWriter(fieldNames, nDim, fullFilename, dataSorter);
  XMLWriter->Write_Data();
  delete XMLWriter;
  
  /*--- Add the dataset to the vtm file ---*/
  
  AddDataset(name, fullFilename + CParaviewXMLFileWriter::fileExt);
  
}
