/*!
 * \file CParaviewVTMFileWriter.cpp
 * \brief Filewriter class for Paraview binary format.
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

#include <utility>

#include "../../../include/output/filewriter/CParaviewVTMFileWriter.hpp"
#include "../../../include/output/filewriter/CParaviewXMLFileWriter.hpp"
#include "../../../../Common/include/toolboxes/printing_toolbox.hpp"
#if defined(_WIN32) || defined(_WIN64) || defined (__WINDOWS__)
#include <direct.h>
#endif

const string CParaviewVTMFileWriter::fileExt = ".vtm";

CParaviewVTMFileWriter::CParaviewVTMFileWriter(su2double valTime,
                                               unsigned short valiZone, unsigned short valnZone)
  : CFileWriter(fileExt), iZone(valiZone), nZone(valnZone), curTime(valTime){

  nWrittenDatasets = 0;
  accumulatedBandwidth = 0;
}


CParaviewVTMFileWriter::~CParaviewVTMFileWriter()= default;

void CParaviewVTMFileWriter::WriteData(string val_filename){

  /*--- We append the pre-defined suffix (extension) to the filename (prefix) ---*/
  val_filename.append(fileExt);

  /*--- If we are in the first zone, create new file and write the file header,
   * otherwise append to already existing file ---*/

  if (rank == MASTER_NODE){
   ofstream multiBlockFile;
    if (iZone == 0)
      multiBlockFile.open (val_filename.c_str());
    else
      multiBlockFile.open(val_filename.c_str(), ios::app);

    if (iZone == 0){
      multiBlockFile << R"(<VTKFile type="vtkMultiBlockDataSet" version="1.0">)" << endl;
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

void CParaviewVTMFileWriter::AddDataset(const string& foldername, string name, const string& file, CParallelDataSorter* dataSorter){

  /*--- Construct the full file name incl. folder ---*/
  /*--- Note that the folder name is simply the filename ---*/

  string fullFilename = foldername + "/zone_" + to_string(iZone) + "/" + file;

  /*--- Create an XML writer and dump data into file ---*/

  CParaviewXMLFileWriter XMLWriter(dataSorter);
  XMLWriter.WriteData(fullFilename);

  /*--- Add the dataset to the vtm file ---*/

  AddDataset(std::move(name), fullFilename + CParaviewXMLFileWriter::fileExt);

  /*--- Update the bandwidth ---*/

  nWrittenDatasets++;

  accumulatedBandwidth += XMLWriter.GetBandwidth();

  bandwidth = accumulatedBandwidth/nWrittenDatasets;
}



void CParaviewVTMFileWriter::WriteFolderData(const string& foldername, CConfig *config,
                                             string multiZoneHeaderString,
                                             CParallelDataSorter* volumeDataSorter,
                                             CParallelDataSorter* surfaceDataSorter,
                                             CGeometry *geometry){

  if (rank == MASTER_NODE){
#if defined(_WIN32) || defined(_WIN64) || defined (__WINDOWS__)
    _mkdir(foldername.c_str());
    _mkdir((foldername + "/zone_" + to_string(iZone)).c_str());
#else
    mkdir(foldername.c_str(), 0777);
    //mkdir((this->folderName + "/zone_" + to_string(valiZone)).c_str(), 0777);
    mkdir((foldername + "/zone_" + to_string(iZone)).c_str(), 0777);
#endif
  }

   /*--- Open a block for the zone ---*/

  StartBlock(std::move(multiZoneHeaderString));

  StartBlock("Internal");
  AddDataset(foldername,"Internal", "Internal", volumeDataSorter);
  EndBlock();

  /*--- Open a block for the boundary ---*/

  StartBlock("Boundary");

  /*--- Loop over all markers used in the config file ---*/

  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_CfgFile(); iMarker++){

    /*--- Get the name of the marker ---*/

    string markerTag = config->GetMarker_CfgFile_TagBound(iMarker);

    /*--- If the current marker can be found on this partition store its name.
     * Note that we have to provide a vector of markers to the sorter routine, although we only do
     * one marker at a time, i.e. ::marker always contains one item. ---*/

    vector<string> marker;
    for (unsigned short jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++){

      /*--- We want to write all markers except send-receive markers ---*/

      if (config->GetMarker_All_TagBound(jMarker) == markerTag &&
        config->GetMarker_All_KindBC(jMarker) != SEND_RECEIVE){
          marker.push_back(markerTag);
        }
    }

    /*--- Only sort if there is at least one processor that has this marker ---*/

    int globalMarkerSize = 0, localMarkerSize = marker.size();
    SU2_MPI::Allreduce(&localMarkerSize, &globalMarkerSize, 1, MPI_INT, MPI_SUM, SU2_MPI::GetComm());

    if (globalMarkerSize > 0){

      /*--- Sort connectivity of the current marker ---*/

      surfaceDataSorter->SortConnectivity(config, geometry, marker);
      surfaceDataSorter->SortOutputData();

      /*--- Add the dataset ---*/

      AddDataset(foldername, markerTag, markerTag, surfaceDataSorter);

    }
  }

  /*--- End "Boundary" block ---*/
  EndBlock();
  /*--- End "Zone" block ---*/
  EndBlock();


}
