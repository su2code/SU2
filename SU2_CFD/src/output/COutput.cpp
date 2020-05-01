/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information
 * \author F. Palacios, T. Economon
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

#include "../../include/output/COutput.hpp"
#include "../../include/output/filewriter/CFVMDataSorter.hpp"
#include "../../include/output/filewriter/CFEMDataSorter.hpp"
#include "../../include/output/filewriter/CSurfaceFVMDataSorter.hpp"
#include "../../include/output/filewriter/CSurfaceFEMDataSorter.hpp"
#include "../../include/output/filewriter/CParaviewFileWriter.hpp"
#include "../../include/output/filewriter/CSTLFileWriter.hpp"
#include "../../include/output/filewriter/CParaviewBinaryFileWriter.hpp"
#include "../../include/output/filewriter/CParaviewXMLFileWriter.hpp"
#include "../../include/output/filewriter/CParaviewVTMFileWriter.hpp"
#include "../../include/output/filewriter/CTecplotFileWriter.hpp"
#include "../../include/output/filewriter/CTecplotBinaryFileWriter.hpp"
#include "../../include/output/filewriter/CCSVFileWriter.hpp"
#include "../../include/output/filewriter/CSU2FileWriter.hpp"
#include "../../include/output/filewriter/CSU2BinaryFileWriter.hpp"
#include "../../include/output/filewriter/CSU2MeshFileWriter.hpp"


#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

COutput::COutput(CConfig *config, unsigned short nDim, bool fem_output): femOutput(fem_output) {

  this->nDim = nDim;

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  fieldWidth = 12;

  convergenceTable = new PrintingToolbox::CTablePrinter(&std::cout);
  multiZoneHeaderTable = new PrintingToolbox::CTablePrinter(&std::cout);
  fileWritingTable = new PrintingToolbox::CTablePrinter(&std::cout);
  historyFileTable = new PrintingToolbox::CTablePrinter(&histFile, "");

  /*--- Set default filenames ---*/

  surfaceFilename = "surface";
  volumeFilename  = "volume";
  restartFilename = "restart";

  /*--- Retrieve the history filename ---*/

  historyFilename = config->GetConv_FileName();

  /*--- Add the correct file extension depending on the file format ---*/

  string hist_ext = ".csv";
  if (config->GetTabular_FileFormat() == TAB_TECPLOT) hist_ext = ".dat";

  /*--- Append the zone ID ---*/

  historyFilename = config->GetMultizone_HistoryFileName(historyFilename, config->GetiZone(), hist_ext);

  /*--- Append the restart iteration ---*/

  if (config->GetTime_Domain() && config->GetRestart()) {
    historyFilename = config->GetUnsteady_FileName(historyFilename, config->GetRestart_Iter(), hist_ext);
  }

  historySep = ",";

  /*--- Initialize residual ---*/

  rhoResNew = EPS;
  rhoResOld = EPS;

  nRequestedHistoryFields = config->GetnHistoryOutput();
  for (unsigned short iField = 0; iField < nRequestedHistoryFields; iField++){
    requestedHistoryFields.push_back(config->GetHistoryOutput_Field(iField));
  }

  nRequestedScreenFields = config->GetnScreenOutput();
  for (unsigned short iField = 0; iField < nRequestedScreenFields; iField++){
    requestedScreenFields.push_back(config->GetScreenOutput_Field(iField));
  }

  nRequestedVolumeFields = config->GetnVolumeOutput();
  for (unsigned short iField = 0; iField < nRequestedVolumeFields; iField++){
    requestedVolumeFields.push_back(config->GetVolumeOutput_Field(iField));
  }

  gridMovement = config->GetGrid_Movement();

  multiZone     = config->GetMultizone_Problem();

  /*--- Default is to write history to file and screen --- */

  noWriting = false;

  /*--- Initialize convergence monitoring structure ---*/

  nCauchy_Elems = config->GetCauchy_Elems();
  cauchyEps = config->GetCauchy_Eps();
  minLogResidual = config->GetMinLogResidual();

  for (unsigned short iField = 0; iField < config->GetnConv_Field(); iField++){
    convFields.emplace_back(config->GetConv_Field(iField));
  }

  newFunc = vector<su2double>(convFields.size());
  oldFunc = vector<su2double>(convFields.size());
  cauchySerie = vector<vector<su2double>>(convFields.size(), vector<su2double>(nCauchy_Elems, 0.0));
  cauchyValue = 0.0;
  convergence = false;

  /*--- Initialize time convergence monitoring structure ---*/

  nWndCauchy_Elems = config->GetWnd_Cauchy_Elems();
  wndCauchyEps     = config->GetWnd_Cauchy_Eps();

  wndConvFields.reserve(config->GetnWndConv_Field());
  for (unsigned short iField = 0; iField < config->GetnWndConv_Field(); iField++){
    wndConvFields.emplace_back(config->GetWndConv_Field(iField));
  }

  WndOld_Func = vector<su2double>(wndConvFields.size());
  WndNew_Func = vector<su2double>(wndConvFields.size());
  WndCauchy_Serie = vector<vector<su2double>>(wndConvFields.size(), vector<su2double>(nWndCauchy_Elems, 0.0));
  WndCauchy_Value = 0.0;
  TimeConvergence = false;

  /*--- Check that the number of cauchy elems is not too large ---*/

  if (nCauchy_Elems > 1000){
    SU2_MPI::Error("Number of Cauchy Elems must be smaller than 1000", CURRENT_FUNCTION);
  }
  if (nWndCauchy_Elems > 1000){
    SU2_MPI::Error("Number of Time Cauchy Elems must be smaller than 1000", CURRENT_FUNCTION);
  }

  /*--- Initialize all convergence flags to false. ---*/

  convergence        = false;

  buildFieldIndexCache = false;

  curInnerIter = 0;
  curOuterIter = 0;
  curTimeIter  = 0;

  volumeDataSorter = nullptr;
  surfaceDataSorter = nullptr;

  headerNeeded = false;

}

COutput::~COutput(void) {

  delete convergenceTable;
  delete multiZoneHeaderTable;
  delete fileWritingTable;
  delete historyFileTable;

  if (volumeDataSorter != nullptr)
    delete volumeDataSorter;

  volumeDataSorter = nullptr;

  if (surfaceDataSorter != nullptr)
    delete surfaceDataSorter;

  surfaceDataSorter = nullptr;


}



void COutput::SetHistory_Output(CGeometry *geometry,
                                  CSolver **solver_container,
                                  CConfig *config,
                                  unsigned long TimeIter,
                                  unsigned long OuterIter,
                                  unsigned long InnerIter) {

  curTimeIter  = TimeIter;
  curAbsTimeIter = TimeIter - config->GetRestart_Iter();
  curOuterIter = OuterIter;
  curInnerIter = InnerIter;

  bool write_header, write_history, write_screen;

  /*--- Retrieve residual and extra data -----------------------------------------------------------------*/

  LoadCommonHistoryData(config);

  LoadHistoryData(config, geometry, solver_container);

  Convergence_Monitoring(config, curInnerIter);

  Postprocess_HistoryData(config);

  MonitorTimeConvergence(config, curTimeIter);

  /*--- Output using only the master node ---*/

  if (rank == MASTER_NODE && !noWriting) {

    /*--- Write the history file ---------------------------------------------------------------------------*/
    write_history = WriteHistoryFile_Output(config);
    if (write_history) SetHistoryFile_Output(config);

    /*--- Write the screen header---------------------------------------------------------------------------*/
    write_header = WriteScreen_Header(config);
    if (write_header) SetScreen_Header(config);

    /*--- Write the screen output---------------------------------------------------------------------------*/
    write_screen = WriteScreen_Output(config);
    if (write_screen) SetScreen_Output(config);

  }

}

void COutput::SetHistory_Output(CGeometry *geometry,
                                CSolver **solver_container,
                                CConfig *config) {

  /*--- Retrieve residual and extra data -----------------------------------------------------------------*/

  LoadCommonHistoryData(config);

  LoadHistoryData(config, geometry, solver_container);

  Convergence_Monitoring(config, curInnerIter);

  Postprocess_HistoryData(config);
}

void COutput::SetMultizoneHistory_Output(COutput **output, CConfig **config, CConfig *driver_config, unsigned long TimeIter, unsigned long OuterIter){

  curTimeIter  = TimeIter;
  curAbsTimeIter = TimeIter - driver_config->GetRestart_Iter();
  curOuterIter = OuterIter;

  bool write_header, write_screen, write_history;

  /*--- Retrieve residual and extra data -----------------------------------------------------------------*/

  LoadCommonHistoryData(driver_config);

  LoadMultizoneHistoryData(output, config);

  Convergence_Monitoring(driver_config, curOuterIter);

  Postprocess_HistoryData(driver_config);

  MonitorTimeConvergence(driver_config, curTimeIter);

  /*--- Output using only the master node ---*/

  if (rank == MASTER_NODE && !noWriting) {

    /*--- Write the history file ---------------------------------------------------------------------------*/
    write_history = WriteHistoryFile_Output(driver_config);
    if (write_history) SetHistoryFile_Output(driver_config);

    /*--- Write the screen header---------------------------------------------------------------------------*/
    write_header = WriteScreen_Header(driver_config);
    if (write_header) SetScreen_Header(driver_config);

    /*--- Write the screen output---------------------------------------------------------------------------*/
    write_screen = WriteScreen_Output(driver_config);
    if (write_screen) SetScreen_Output(driver_config);

  }

}

void COutput::AllocateDataSorters(CConfig *config, CGeometry *geometry){

  /*---- Construct a data sorter object to partition and distribute
   *  the local data into linear chunks across the processors ---*/

  if (femOutput){

    if (volumeDataSorter == nullptr)
      volumeDataSorter = new CFEMDataSorter(config, geometry, volumeFieldNames);

    if (surfaceDataSorter == nullptr)
      surfaceDataSorter = new CSurfaceFEMDataSorter(config, geometry,
                                                  dynamic_cast<CFEMDataSorter*>(volumeDataSorter));

  }  else {

    if (volumeDataSorter == nullptr)
      volumeDataSorter = new CFVMDataSorter(config, geometry, volumeFieldNames);

    if (surfaceDataSorter == nullptr)
      surfaceDataSorter = new CSurfaceFVMDataSorter(config, geometry,
                                                  dynamic_cast<CFVMDataSorter*>(volumeDataSorter));

  }

}

void COutput::Load_Data(CGeometry *geometry, CConfig *config, CSolver** solver_container){

  /*--- Check if the data sorters are allocated, if not, allocate them. --- */

  AllocateDataSorters(config, geometry);

  /*--- Loop over all points and store the requested volume output data into the data sorter objects ---*/

  LoadDataIntoSorter(config, geometry, solver_container);

  /*--- Partition and sort the volume output data -- */

  volumeDataSorter->SortOutputData();

}

void COutput::WriteToFile(CConfig *config, CGeometry *geometry, unsigned short format, string fileName){

  CFileWriter *fileWriter = NULL;

  unsigned short lastindex = fileName.find_last_of(".");
  fileName = fileName.substr(0, lastindex);

  /*--- Write files depending on the format --- */

  switch (format) {

    case SURFACE_CSV:

      if (fileName.empty())
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);

      surfaceDataSorter->SortConnectivity(config, geometry);
      surfaceDataSorter->SortOutputData();

      if (rank == MASTER_NODE) {
        (*fileWritingTable) << "CSV file" << fileName + CSU2FileWriter::fileExt;
      }

      fileWriter = new CSU2FileWriter(fileName, surfaceDataSorter);

      break;

    case RESTART_ASCII: case CSV:

      if (fileName.empty())
        fileName = config->GetFilename(restartFilename, "", curTimeIter);

      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "SU2 ASCII restart" << fileName + CSU2FileWriter::fileExt;
      }

      fileWriter = new CSU2FileWriter(fileName, volumeDataSorter);

      break;

    case RESTART_BINARY:

      if (fileName.empty())
        fileName = config->GetFilename(restartFilename, "", curTimeIter);

      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "SU2 restart" << fileName + CSU2BinaryFileWriter::fileExt;
      }

      fileWriter = new CSU2BinaryFileWriter(fileName, volumeDataSorter);

      break;

    case MESH:

      if (fileName.empty())
        fileName = volumeFilename;

      /*--- Load and sort the output data and connectivity. ---*/

      volumeDataSorter->SortConnectivity(config, geometry, true);

      /*--- Set the mesh ASCII format ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "SU2 mesh" << fileName + CSU2MeshFileWriter::fileExt;
      }

      fileWriter = new CSU2MeshFileWriter(fileName, volumeDataSorter,
                                          config->GetiZone(), config->GetnZone());


      break;

    case TECPLOT_BINARY:

      if (fileName.empty())
        fileName = config->GetFilename(volumeFilename, "", curTimeIter);

      /*--- Load and sort the output data and connectivity. ---*/

      volumeDataSorter->SortConnectivity(config, geometry, false);

      /*--- Write tecplot binary ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Tecplot binary" << fileName + CTecplotBinaryFileWriter::fileExt;
      }

      fileWriter = new CTecplotBinaryFileWriter(fileName, volumeDataSorter,
                                                curTimeIter, GetHistoryFieldValue("TIME_STEP"));

      break;

    case TECPLOT:

      if (fileName.empty())
        fileName = config->GetFilename(volumeFilename, "", curTimeIter);

      /*--- Load and sort the output data and connectivity. ---*/

      volumeDataSorter->SortConnectivity(config, geometry, true);

      /*--- Write tecplot ascii ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Tecplot ASCII" << fileName + CTecplotFileWriter::fileExt;
      }

      fileWriter = new CTecplotFileWriter(fileName, volumeDataSorter,
                                          curTimeIter, GetHistoryFieldValue("TIME_STEP"));

      break;

    case PARAVIEW_XML:

      if (fileName.empty())
        fileName = config->GetFilename(volumeFilename, "", curTimeIter);

      /*--- Load and sort the output data and connectivity. ---*/

      volumeDataSorter->SortConnectivity(config, geometry, true);

      /*--- Write paraview binary ---*/
      if (rank == MASTER_NODE) {
        (*fileWritingTable) << "Paraview" << fileName + CParaviewXMLFileWriter::fileExt;
      }

      fileWriter = new CParaviewXMLFileWriter(fileName, volumeDataSorter);

      break;

    case PARAVIEW_BINARY:

      if (fileName.empty())
        fileName = config->GetFilename(volumeFilename, "", curTimeIter);

      /*--- Load and sort the output data and connectivity. ---*/

      volumeDataSorter->SortConnectivity(config, geometry, true);

      /*--- Write paraview binary ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Paraview binary" << fileName + CParaviewBinaryFileWriter::fileExt;
      }

      fileWriter = new CParaviewBinaryFileWriter(fileName, volumeDataSorter);

      break;

    case PARAVIEW_MULTIBLOCK:
      {

        if (fileName.empty())
          fileName = config->GetFilename(volumeFilename, "", curTimeIter);

        /*--- Sort volume connectivity ---*/

        volumeDataSorter->SortConnectivity(config, geometry, true);

        /*--- The file name of the multiblock file is the case name (i.e. the config file name w/o ext.) ---*/

        fileName = config->GetUnsteady_FileName(config->GetCaseName(), curTimeIter, "");

        /*--- Allocate the vtm file writer ---*/

        fileWriter = new CParaviewVTMFileWriter(fileName, fileName, GetHistoryFieldValue("CUR_TIME"),
                                                config->GetiZone(), config->GetnZone());

        /*--- We cast the pointer to its true type, to avoid virtual functions ---*/

        CParaviewVTMFileWriter* vtmWriter = dynamic_cast<CParaviewVTMFileWriter*>(fileWriter);

        if (rank == MASTER_NODE) {
            (*fileWritingTable) << "Paraview Multiblock"
                                << fileName + CParaviewVTMFileWriter::fileExt;
        }

        /*--- Open a block for the zone ---*/

        vtmWriter->StartBlock(multiZoneHeaderString);

        fileName = "Internal";

        /*--- Open a block for the internal (volume) data and add the dataset ---*/

        vtmWriter->StartBlock(fileName);
        vtmWriter->AddDataset(fileName, fileName, volumeDataSorter);
        vtmWriter->EndBlock();

        /*--- Open a block for the boundary ---*/

        vtmWriter->StartBlock("Boundary");

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
          SU2_MPI::Allreduce(&localMarkerSize, &globalMarkerSize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

          if (globalMarkerSize > 0){

            /*--- Sort connectivity of the current marker ---*/

            surfaceDataSorter->SortConnectivity(config, geometry, marker);
            surfaceDataSorter->SortOutputData();

            /*--- Add the dataset ---*/

            vtmWriter->AddDataset(markerTag, markerTag, surfaceDataSorter);

          }
        }
        /*--- End "Boundary" block ---*/
        vtmWriter->EndBlock();
        /*--- End "Zone" block ---*/
        vtmWriter->EndBlock();
      }


      break;

    case PARAVIEW:

      if (fileName.empty())
        fileName = config->GetFilename(volumeFilename, "", curTimeIter);

      /*--- Load and sort the output data and connectivity. ---*/

      volumeDataSorter->SortConnectivity(config, geometry, true);

      /*--- Write paraview ascii ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Paraview ASCII" << fileName + CParaviewFileWriter::fileExt;
      }

      fileWriter = new CParaviewFileWriter(fileName, volumeDataSorter);

      break;

    case SURFACE_PARAVIEW:

      if (fileName.empty())
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);

      /*--- Load and sort the output data and connectivity. ---*/

      surfaceDataSorter->SortConnectivity(config, geometry);
      surfaceDataSorter->SortOutputData();

      /*--- Write surface paraview ascii ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Paraview ASCII surface" << fileName + CParaviewFileWriter::fileExt;
      }

      fileWriter = new CParaviewFileWriter(fileName, surfaceDataSorter);

      break;

    case SURFACE_PARAVIEW_BINARY:

      if (fileName.empty())
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);

      /*--- Load and sort the output data and connectivity. ---*/

      surfaceDataSorter->SortConnectivity(config, geometry);
      surfaceDataSorter->SortOutputData();

      /*--- Write surface paraview binary ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Paraview binary surface" << fileName + CParaviewBinaryFileWriter::fileExt;
      }

      fileWriter = new CParaviewBinaryFileWriter(fileName, surfaceDataSorter);

      break;

    case SURFACE_PARAVIEW_XML:

      if (fileName.empty())
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);

      /*--- Load and sort the output data and connectivity. ---*/

      surfaceDataSorter->SortConnectivity(config, geometry);
      surfaceDataSorter->SortOutputData();

      /*--- Write paraview binary ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Paraview surface" << fileName + CParaviewXMLFileWriter::fileExt;
      }

      fileWriter = new CParaviewXMLFileWriter(fileName, surfaceDataSorter);

      break;

    case SURFACE_TECPLOT:

      if (fileName.empty())
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);

      /*--- Load and sort the output data and connectivity. ---*/

      surfaceDataSorter->SortConnectivity(config, geometry);
      surfaceDataSorter->SortOutputData();

      /*--- Write surface tecplot ascii ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Tecplot ASCII surface" << fileName + CTecplotFileWriter::fileExt;
      }

      fileWriter = new CTecplotFileWriter(fileName, surfaceDataSorter,
                                          curTimeIter, GetHistoryFieldValue("TIME_STEP"));

      break;

    case SURFACE_TECPLOT_BINARY:

      if (fileName.empty())
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);

      /*--- Load and sort the output data and connectivity. ---*/

      surfaceDataSorter->SortConnectivity(config, geometry);
      surfaceDataSorter->SortOutputData();

      /*--- Write surface tecplot binary ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "Tecplot binary surface" << fileName + CTecplotBinaryFileWriter::fileExt;
      }

      fileWriter = new CTecplotBinaryFileWriter(fileName, surfaceDataSorter,
                                                curTimeIter, GetHistoryFieldValue("TIME_STEP"));

      break;

    case STL:

      if (fileName.empty())
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);

      /*--- Load and sort the output data and connectivity. ---*/

      surfaceDataSorter->SortConnectivity(config, geometry);
      surfaceDataSorter->SortOutputData();

      /*--- Write ASCII STL ---*/
      if (rank == MASTER_NODE) {
          (*fileWritingTable) << "STL ASCII" << fileName + CSTLFileWriter::fileExt;
      }

      fileWriter = new CSTLFileWriter(fileName, surfaceDataSorter);

      break;

    default:
      fileWriter = NULL;
      break;
  }

  if (fileWriter != NULL){

    /*--- Write data to file ---*/

    fileWriter->Write_Data();

    su2double BandWidth = fileWriter->Get_Bandwidth();

    /*--- Compute and store the bandwidth ---*/

    if (format == RESTART_BINARY){
      config->SetRestart_Bandwidth_Agg(config->GetRestart_Bandwidth_Agg()+BandWidth);
    }

    if (config->GetWrt_Performance() && (rank == MASTER_NODE)){
      fileWritingTable->SetAlign(PrintingToolbox::CTablePrinter::RIGHT);
      (*fileWritingTable) << " " << "(" + PrintingToolbox::to_string(BandWidth) + " MB/s)";
      fileWritingTable->SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    }

    delete fileWriter;

  }
}



bool COutput::SetResult_Files(CGeometry *geometry, CConfig *config, CSolver** solver_container,
                              unsigned long iter, bool force_writing){

  bool writeFiles = WriteVolume_Output(config, iter, force_writing);

  /*--- Check if the data sorters are allocated, if not, allocate them. --- */

  AllocateDataSorters(config, geometry);

  /*--- Collect the volume data from the solvers.
   *  If time-domain is enabled, we also load the data although we don't output it,
   *  since we might want to do time-averaging. ---*/

  if (writeFiles || config->GetTime_Domain())
    LoadDataIntoSorter(config, geometry, solver_container);

  if (writeFiles){

    /*--- Partition and sort the data --- */

    volumeDataSorter->SortOutputData();

    unsigned short nVolumeFiles = config->GetnVolumeOutputFiles();
    unsigned short *VolumeFiles = config->GetVolumeOutputFiles();

    if (rank == MASTER_NODE && nVolumeFiles != 0){
      fileWritingTable->SetAlign(PrintingToolbox::CTablePrinter::CENTER);
      fileWritingTable->PrintHeader();
      fileWritingTable->SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    }

    /*--- Loop through all requested output files and write
     * the partitioned and sorted data stored in the data sorters. ---*/

    for (unsigned short iFile = 0; iFile < nVolumeFiles; iFile++){

      WriteToFile(config, geometry, VolumeFiles[iFile]);

    }

    if (rank == MASTER_NODE && nVolumeFiles != 0){
      fileWritingTable->PrintFooter();
      headerNeeded = true;
    }

    /*--- Write any additonal files defined in the child class ----*/

    WriteAdditionalFiles(config, geometry, solver_container);

    return true;
  }

  return false;
}

void COutput::PrintConvergenceSummary(){

  PrintingToolbox::CTablePrinter  ConvSummary(&cout);

  ConvSummary.AddColumn("Convergence Field", 28);
  ConvSummary.AddColumn("Value", 14);
  ConvSummary.AddColumn("Criterion", 14);
  ConvSummary.AddColumn("Converged",12);
  ConvSummary.SetAlign(PrintingToolbox::CTablePrinter::CENTER);
  ConvSummary.PrintHeader();
  for (unsigned short iField_Conv = 0; iField_Conv < convFields.size(); iField_Conv++){
    const string &convField = convFields[iField_Conv];
    if (historyOutput_Map.at(convField).fieldType == HistoryFieldType::COEFFICIENT) {
      string convMark = "No";
      if ( historyOutput_Map.at("CAUCHY_" + convField).value < cauchyEps) convMark = "Yes";
      ConvSummary << historyOutput_Map.at("CAUCHY_" + convField).fieldName
          <<  historyOutput_Map.at("CAUCHY_" + convField).value
          << " < " + PrintingToolbox::to_string(cauchyEps) << convMark;
    }
    else if (historyOutput_Map.at(convField).fieldType == HistoryFieldType::RESIDUAL ||
        historyOutput_Map.at(convField).fieldType == HistoryFieldType::AUTO_RESIDUAL)  {
      string convMark = "No";
      if (historyOutput_Map.at(convField).value < minLogResidual) convMark = "Yes";
      ConvSummary << historyOutput_Map.at(convField).fieldName
          << historyOutput_Map.at(convField).value
          << " < " + PrintingToolbox::to_string(minLogResidual) << convMark;
    }
  }
  ConvSummary.PrintFooter();
}

bool COutput::Convergence_Monitoring(CConfig *config, unsigned long Iteration) {

  unsigned short iCounter;

  convergence = true;

  for (unsigned short iField_Conv = 0; iField_Conv < convFields.size(); iField_Conv++){

    bool fieldConverged = false;

    const string &convField = convFields[iField_Conv];
    if (historyOutput_Map.count(convField) > 0){
      su2double monitor = historyOutput_Map.at(convField).value;

      /*--- Cauchy based convergence criteria ---*/

      if (historyOutput_Map.at(convField).fieldType == HistoryFieldType::COEFFICIENT) {

        if (Iteration == 0){
          for (iCounter = 0; iCounter < nCauchy_Elems; iCounter++){
            cauchySerie[iField_Conv][iCounter] = 0.0;
          }
          newFunc[iField_Conv] = monitor;
        }

        oldFunc[iField_Conv] = newFunc[iField_Conv];
        newFunc[iField_Conv] = monitor;
        cauchyFunc = fabs(newFunc[iField_Conv] - oldFunc[iField_Conv])/fabs(monitor);

        cauchySerie[iField_Conv][Iteration % nCauchy_Elems] = cauchyFunc;
        cauchyValue = 0.0;
        for (iCounter = 0; iCounter < nCauchy_Elems; iCounter++)
          cauchyValue += cauchySerie[iField_Conv][iCounter];

        cauchyValue /= nCauchy_Elems;

        if (cauchyValue >= cauchyEps) { fieldConverged = false;}
        else { fieldConverged = true;}

        /*--- Start monitoring only if the current iteration
         *  is larger than the number of cauchy elements and
         * the number of start-up iterations --- */

        if (Iteration < max(config->GetStartConv_Iter(), nCauchy_Elems)){
          fieldConverged = false;
        }

        SetHistoryOutputValue("CAUCHY_" + convField, cauchyValue);

        if(Iteration == 0){
          SetHistoryOutputValue("CAUCHY_" + convField, 1.0);
        }
      }


      /*--- Residual based convergence criteria ---*/

      if (historyOutput_Map.at(convField).fieldType == HistoryFieldType::RESIDUAL ||
          historyOutput_Map.at(convField).fieldType == HistoryFieldType::AUTO_RESIDUAL) {

        /*--- Check the convergence ---*/

        if (Iteration != 0 && (monitor <= minLogResidual)) { fieldConverged = true;  }
        else { fieldConverged = false; }

      }

      /*--- Do not apply any convergence criteria of the number
     of iterations is less than a particular value ---*/

      if (Iteration < config->GetStartConv_Iter()) {
        fieldConverged = false;
      }

      convergence = fieldConverged && convergence;
    }
  }

  if (convFields.empty()) convergence = false;

  /*--- Apply the same convergence criteria to all the processors ---*/

#ifdef HAVE_MPI

  unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
  sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
  rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;

  /*--- Convergence criteria ---*/

  sbuf_conv[0] = convergence;
  SU2_MPI::Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);

  /*-- Compute global convergence criteria in the master node --*/

  sbuf_conv[0] = 0;
  if (rank == MASTER_NODE) {
    if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
    else sbuf_conv[0] = 0;
  }

  SU2_MPI::Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);

  if (sbuf_conv[0] == 1) { convergence = true; }
  else { convergence = false;  }

  delete [] sbuf_conv;
  delete [] rbuf_conv;

#endif

  return convergence;
}

bool COutput::MonitorTimeConvergence(CConfig *config, unsigned long TimeIteration) {

  bool Inner_IterConv = GetConvergence() || config->GetnInner_Iter()-1 <= curInnerIter; //Check, if Inner_Iter is converged

  if(TimeIteration == 0){
    for (unsigned short iField_Conv = 0; iField_Conv < wndConvFields.size(); iField_Conv++){
      const string WndConv_Field= wndConvFields[iField_Conv];
      if (historyOutput_Map.count(WndConv_Field) > 0){
        SetHistoryOutputValue("CAUCHY_"+ WndConv_Field, 1.0);
      }
    }
  }
  if(Inner_IterConv && TimeIteration >= config->GetStartWindowIteration()){
    TimeConvergence = true;
    unsigned short iCounter;

    for (unsigned short iField_Conv = 0; iField_Conv < wndConvFields.size(); iField_Conv++){
      bool fieldConverged = false;
      const string WndConv_Field= wndConvFields[iField_Conv];

      if (historyOutput_Map.count(WndConv_Field) > 0){
        su2double monitor = historyOutput_Map[WndConv_Field].value;

        /*--- Cauchy based convergence criteria ---*/

        if (historyOutput_Map[WndConv_Field].fieldType == HistoryFieldType::AUTO_COEFFICIENT) { //TAVG values are AUTO_COEFF
          if (TimeIteration == config->GetStartWindowIteration()){
            for (iCounter = 0; iCounter < nWndCauchy_Elems; iCounter++){
              WndCauchy_Serie[iField_Conv][iCounter] = 0.0;
            }
            WndNew_Func[iField_Conv] = monitor;
          }
          WndOld_Func[iField_Conv] = WndNew_Func[iField_Conv];
          WndNew_Func[iField_Conv] = monitor;
          WndCauchy_Func = fabs(WndNew_Func[iField_Conv] - WndOld_Func[iField_Conv]);
          WndCauchy_Serie[iField_Conv][TimeIteration % nWndCauchy_Elems] = WndCauchy_Func;
          WndCauchy_Value = 1.0;

          if (TimeIteration >= nWndCauchy_Elems+config->GetStartWindowIteration()){
            WndCauchy_Value = 0.0;
            for (iCounter = 0; iCounter < nWndCauchy_Elems; iCounter++){
              WndCauchy_Value += WndCauchy_Serie[iField_Conv][iCounter];
            }
            WndCauchy_Value /= nWndCauchy_Elems;
          }
          if (WndCauchy_Value >= wndCauchyEps){fieldConverged = false;}
          else{fieldConverged = true;}

          /*--- Start monitoring only if the current iteration is larger than the
           *  number of cauchy elements and the number of start-up iterations ---*/

          if (TimeIteration <  config->GetStartWindowIteration() + max(config->GetWnd_StartConv_Iter(), nWndCauchy_Elems)){
            fieldConverged = false;
          }
          SetHistoryOutputValue("CAUCHY_" + WndConv_Field, WndCauchy_Value);
        }
        TimeConvergence = fieldConverged && TimeConvergence;

        /*--- Stop the simulation in case a nan appears, do not save the solution ---*/

        if (monitor != monitor){
          SU2_MPI::Error("SU2 has diverged (NaN detected).", CURRENT_FUNCTION);}
      }
    }

    /*--- Do not apply any convergence criterion if the option is disabled. */
    if(!config->GetWnd_Cauchy_Crit()){TimeConvergence = false;}
    if(wndConvFields.empty()){TimeConvergence = false;}
  }
  return TimeConvergence;
}

void COutput::SetHistoryFile_Header(CConfig *config) {

  unsigned short iField_Output = 0,
      iReqField = 0,
      iMarker = 0;
  stringstream out;
  int width = 20;

  for (iField_Output = 0; iField_Output < historyOutput_List.size(); iField_Output++){
    const string &fieldIdentifier = historyOutput_List[iField_Output];
    const HistoryOutputField &field = historyOutput_Map.at(fieldIdentifier);
    for (iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      const string & requestedField = requestedHistoryFields[iReqField];
      if (requestedField == field.outputGroup || (requestedField == fieldIdentifier)){
        if (field.screenFormat == ScreenOutputFormat::INTEGER) width = std::max((int)field.fieldName.size()+2, 10);
        else{ width = std::max((int)field.fieldName.size()+2, 18);}
        historyFileTable->AddColumn("\"" + field.fieldName + "\"", width);
      }
    }
  }

  for (iField_Output = 0; iField_Output < historyOutputPerSurface_List.size(); iField_Output++){
    const string &fieldIdentifier = historyOutputPerSurface_List[iField_Output];
    for (iMarker = 0; iMarker < historyOutputPerSurface_Map[fieldIdentifier].size(); iMarker++){
      const HistoryOutputField &field = historyOutputPerSurface_Map.at(fieldIdentifier)[iMarker];
      for (iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
        const string &requestedField = requestedHistoryFields[iReqField];
        if (requestedField == field.outputGroup || (requestedField == fieldIdentifier)){
          if (field.screenFormat == ScreenOutputFormat::INTEGER) width = std::max((int)field.fieldName.size()+2, 10);
          else{ width = std::max((int)field.fieldName.size()+2, 18);}
          historyFileTable->AddColumn("\"" + field.fieldName + "\"", width);
        }
      }
    }
  }

  if (config->GetTabular_FileFormat() == TAB_TECPLOT) {
    histFile << "VARIABLES = \\" << endl;
  }
  historyFileTable->PrintHeader();
  histFile.flush();
}


void COutput::SetHistoryFile_Output(CConfig *config) {

  unsigned short iField_Output = 0,
      iReqField = 0,
      iMarker = 0;
  stringstream out;

  for (iField_Output = 0; iField_Output < historyOutput_List.size(); iField_Output++){
    const string &fieldIdentifier = historyOutput_List[iField_Output];
    const HistoryOutputField &field = historyOutput_Map.at(fieldIdentifier);
    for (iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      const string &RequestedField = requestedHistoryFields[iReqField];
      if (RequestedField == field.outputGroup){
        (*historyFileTable) << field.value;
      }
    }
  }

  for (iField_Output = 0; iField_Output < historyOutputPerSurface_List.size(); iField_Output++){
    const string &fieldIdentifier = historyOutputPerSurface_List[iField_Output];
    for (iMarker = 0; iMarker < historyOutputPerSurface_Map[fieldIdentifier].size(); iMarker++){
      const HistoryOutputField &field = historyOutputPerSurface_Map.at(fieldIdentifier)[iMarker];
      for (iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
        const string &RequestedField = requestedHistoryFields[iReqField];
        if (RequestedField == field.outputGroup){
          (*historyFileTable) << field.value;
        }
      }
    }
  }

  /*--- Print the string to file and remove the last two characters (a separator and a space) ---*/

  histFile.flush();
}

void COutput::SetScreen_Header(CConfig *config) {
  if (config->GetMultizone_Problem())
    multiZoneHeaderTable->PrintHeader();
  convergenceTable->PrintHeader();
}


void COutput::SetScreen_Output(CConfig *config) {

  string RequestedField;

  for (unsigned short iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
    stringstream out;
    RequestedField = requestedScreenFields[iReqField];
    if (historyOutput_Map.count(RequestedField) > 0){
      switch (historyOutput_Map.at(RequestedField).screenFormat) {
        case ScreenOutputFormat::INTEGER:
          PrintingToolbox::PrintScreenInteger(out, SU2_TYPE::Int(historyOutput_Map.at(RequestedField).value), fieldWidth);
          break;
        case ScreenOutputFormat::FIXED:
          PrintingToolbox::PrintScreenFixed(out, historyOutput_Map.at(RequestedField).value, fieldWidth);
          break;
        case ScreenOutputFormat::SCIENTIFIC:
          PrintingToolbox::PrintScreenScientific(out, historyOutput_Map.at(RequestedField).value, fieldWidth);
          break;
        case ScreenOutputFormat::PERCENT:
          PrintingToolbox::PrintScreenPercent(out, historyOutput_Map[RequestedField].value, fieldWidth);
          break;
      }
    }
    if (historyOutputPerSurface_Map.count(RequestedField) > 0){
      switch (historyOutputPerSurface_Map.at(RequestedField)[0].screenFormat) {
        case ScreenOutputFormat::INTEGER:
          PrintingToolbox::PrintScreenInteger(out, SU2_TYPE::Int(historyOutputPerSurface_Map.at(RequestedField)[0].value), fieldWidth);
          break;
        case ScreenOutputFormat::FIXED:
          PrintingToolbox::PrintScreenFixed(out, historyOutputPerSurface_Map.at(RequestedField)[0].value, fieldWidth);
          break;
        case ScreenOutputFormat::SCIENTIFIC:
          PrintingToolbox::PrintScreenScientific(out, historyOutputPerSurface_Map.at(RequestedField)[0].value, fieldWidth);
          break;
        case ScreenOutputFormat::PERCENT:
          PrintingToolbox::PrintScreenPercent(out, historyOutputPerSurface_Map[RequestedField][0].value, fieldWidth);
          break;
      }
    }
    (*convergenceTable) << out.str();
  }
  SetAdditionalScreenOutput(config);
}

void COutput::PreprocessHistoryOutput(CConfig *config, bool wrt){

    noWriting = !wrt;

    /*--- Set the common output fields ---*/

    SetCommonHistoryFields(config);

    /*--- Set the History output fields using a virtual function call to the child implementation ---*/

    SetHistoryOutputFields(config);

    /*--- Postprocess the history fields. Creates new fields based on the ones set in the child classes ---*/

    Postprocess_HistoryFields(config);

    /*--- We use a fixed size of the file output summary table ---*/

    int total_width = 72;
    fileWritingTable->AddColumn("File Writing Summary", (total_width)/2-1);
    fileWritingTable->AddColumn("Filename", total_width/2-1);
    fileWritingTable->SetAlign(PrintingToolbox::CTablePrinter::LEFT);

    /*--- Check for consistency and remove fields that are requested but not available --- */

    CheckHistoryOutput();

    if (rank == MASTER_NODE && !noWriting){

      /*--- Open history file and print the header ---*/
      if (!config->GetMultizone_Problem() || config->GetWrt_ZoneHist())
        PrepareHistoryFile(config);

      total_width = nRequestedScreenFields*fieldWidth + (nRequestedScreenFields-1);

      /*--- Set the multizone screen header ---*/

      if (config->GetMultizone_Problem()){
        multiZoneHeaderTable->AddColumn(multiZoneHeaderString, total_width);
        multiZoneHeaderTable->SetAlign(PrintingToolbox::CTablePrinter::CENTER);
        multiZoneHeaderTable->SetPrintHeaderBottomLine(false);
      }

    }

}

void COutput::PreprocessMultizoneHistoryOutput(COutput **output, CConfig **config, CConfig* driver_config, bool wrt){

  noWriting = !wrt;

  /*--- Set the common history fields for all solvers ---*/

  SetCommonHistoryFields(driver_config);

  /*--- Set the History output fields using a virtual function call to the child implementation ---*/

  SetMultizoneHistoryOutputFields(output, config);

  /*--- Postprocess the history fields. Creates new fields based on the ones set in the child classes ---*/

  Postprocess_HistoryFields(driver_config);

  /*--- We use a fixed size of the file output summary table ---*/

  int total_width = 72;
  fileWritingTable->AddColumn("File Writing Summary", (total_width-1)/2);
  fileWritingTable->AddColumn("Filename", total_width/2);
  fileWritingTable->SetAlign(PrintingToolbox::CTablePrinter::LEFT);

  /*--- Check for consistency and remove fields that are requested but not available --- */

  CheckHistoryOutput();

  if (rank == MASTER_NODE && !noWriting){

    /*--- Open history file and print the header ---*/

    PrepareHistoryFile(driver_config);

    total_width = nRequestedScreenFields*fieldWidth + (nRequestedScreenFields-1);

    /*--- Set the multizone screen header ---*/

    if (config[ZONE_0]->GetMultizone_Problem()){
      multiZoneHeaderTable->AddColumn(multiZoneHeaderString, nRequestedScreenFields*fieldWidth + (nRequestedScreenFields-1));
      multiZoneHeaderTable->SetAlign(PrintingToolbox::CTablePrinter::CENTER);
      multiZoneHeaderTable->SetPrintHeaderBottomLine(false);
    }

  }

}

void COutput::PrepareHistoryFile(CConfig *config){

  /*--- Open the history file ---*/

  histFile.open(historyFilename.c_str(), ios::out);

  /*--- Create and format the history file table ---*/

  historyFileTable->SetInnerSeparator(historySep);
  historyFileTable->SetAlign(PrintingToolbox::CTablePrinter::CENTER);
  historyFileTable->SetPrintHeaderTopLine(false);
  historyFileTable->SetPrintHeaderBottomLine(false);
  historyFileTable->SetPrecision(10);

  /*--- Add the header to the history file. ---*/

  SetHistoryFile_Header(config);

}

void COutput::CheckHistoryOutput(){


  /*--- Set screen convergence output header and remove unavailable fields ---*/

  string requestedField;
  vector<string> FieldsToRemove;
  vector<bool> FoundField(nRequestedHistoryFields, false);

  for (unsigned short iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
    requestedField = requestedScreenFields[iReqField];
    if (historyOutput_Map.count(requestedField) > 0){
      convergenceTable->AddColumn(historyOutput_Map.at(requestedField).fieldName, fieldWidth);
    }
    else if (historyOutputPerSurface_Map.count(requestedField) > 0){
      convergenceTable->AddColumn(historyOutputPerSurface_Map.at(requestedField)[0].fieldName, fieldWidth);
    }else {
      FieldsToRemove.push_back(requestedField);
    }
  }

  /*--- Remove fields which are not defined --- */

  for (unsigned short iReqField = 0; iReqField < FieldsToRemove.size(); iReqField++){
    if (rank == MASTER_NODE) {
      if (iReqField == 0){
        cout << "  Info: Ignoring the following screen output fields:" << endl;
        cout << "  ";
      }        cout << FieldsToRemove[iReqField];
      if (iReqField != FieldsToRemove.size()-1){
        cout << ", ";
      } else {
        cout << endl;
      }
    }
    requestedScreenFields.erase(std::find(requestedScreenFields.begin(),
                                          requestedScreenFields.end(), FieldsToRemove[iReqField]));
  }

  nRequestedScreenFields = requestedScreenFields.size();

  if (rank == MASTER_NODE){
    cout <<"Screen output fields: ";
    for (unsigned short iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
      requestedField = requestedScreenFields[iReqField];
      cout << requestedScreenFields[iReqField];
      if (iReqField != nRequestedScreenFields - 1) cout << ", ";
    }
    cout << endl;
  }

  /*--- Remove unavailable fields from the history file output ---*/

  FieldsToRemove.clear();
  FoundField = vector<bool>(nRequestedHistoryFields, false);

  for (unsigned short iField_Output = 0; iField_Output < historyOutput_List.size(); iField_Output++){
    const string &fieldReference = historyOutput_List[iField_Output];
    if (historyOutput_Map.count(fieldReference) > 0){
      const HistoryOutputField &field = historyOutput_Map.at(fieldReference);
      for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
        requestedField = requestedHistoryFields[iReqField];
        if (requestedField == field.outputGroup){
          FoundField[iReqField] = true;
        }
      }
    }
  }

  for (unsigned short iField_Output = 0; iField_Output < historyOutputPerSurface_List.size(); iField_Output++){
    const string &fieldReference = historyOutputPerSurface_List[iField_Output];
    if (historyOutputPerSurface_Map.count(fieldReference) > 0){
      for (unsigned short iMarker = 0; iMarker < historyOutputPerSurface_Map.at(fieldReference).size(); iMarker++){
        const HistoryOutputField &Field = historyOutputPerSurface_Map.at(fieldReference)[iMarker];
        for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
          requestedField = requestedHistoryFields[iReqField];
          if (requestedField == Field.outputGroup){
            FoundField[iReqField] = true;
          }
        }
      }
    }
  }

  for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
    if (!FoundField[iReqField]){
      FieldsToRemove.push_back(requestedHistoryFields[iReqField]);
    }
  }

  /*--- Remove fields which are not defined --- */

  for (unsigned short iReqField = 0; iReqField < FieldsToRemove.size(); iReqField++){
    if (rank == MASTER_NODE) {
      if (iReqField == 0){
        cout << "  Info: Ignoring the following history output groups:" << endl;
        cout << "  ";
      }        cout << FieldsToRemove[iReqField];
      if (iReqField != FieldsToRemove.size()-1){
        cout << ", ";
      } else {
        cout << endl;
      }
    }
    requestedHistoryFields.erase(std::find(requestedHistoryFields.begin(),
                                           requestedHistoryFields.end(), FieldsToRemove[iReqField]));
  }

  nRequestedHistoryFields = requestedHistoryFields.size();

  if (rank == MASTER_NODE){
    cout <<"History output group(s): ";
    for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      requestedField = requestedHistoryFields[iReqField];
      cout << requestedHistoryFields[iReqField];
      if (iReqField != nRequestedHistoryFields - 1) cout << ", ";
    }
    cout << endl;
  }

  /*--- Check that the requested convergence monitoring field is available ---*/
  bool removedField = false;
  FieldsToRemove.clear();
  for (unsigned short iField_Conv = 0; iField_Conv < convFields.size(); iField_Conv++){
    if (historyOutput_Map.count(convFields[iField_Conv]) == 0){
      if (!removedField) {
        if(rank == MASTER_NODE) cout << "Ignoring Convergence Field(s): ";
        removedField = true;
      }
      if(rank == MASTER_NODE) cout << convFields[iField_Conv] << " ";
      FieldsToRemove.push_back(convFields[iField_Conv]);
    }
  }
  if (removedField && (rank == MASTER_NODE)) cout << endl;
  for (unsigned short iField_Conv = 0; iField_Conv < FieldsToRemove.size(); iField_Conv++){
    convFields.erase(std::find(convFields.begin(),
                               convFields.end(), FieldsToRemove[iField_Conv]));
  }

  if (rank == MASTER_NODE){
    if(convFields.empty()){
      cout << "Warning: No (valid) fields chosen for convergence monitoring. Convergence monitoring inactive."<<  endl;
    }
    else{
      cout <<"Convergence field(s): ";
      for (unsigned short iField_Conv = 0; iField_Conv < convFields.size(); iField_Conv++){
        cout << convFields[iField_Conv];
        if (iField_Conv != convFields.size() - 1) cout << ", ";
      }
      cout << endl;
    }
  }

  /*--- Check that the requested time convergene monitoring field(s) are available*/
  removedField = false;
  FieldsToRemove.clear();
  for (unsigned short iField_Conv = 0; iField_Conv < wndConvFields.size(); iField_Conv++){
    if (historyOutput_Map.count(wndConvFields[iField_Conv]) == 0){
      if (!removedField) {
        if(rank == MASTER_NODE) cout << "Ignoring Time Convergence Field(s): ";
        removedField = true;
      }
      if(rank == MASTER_NODE)cout << wndConvFields[iField_Conv] << " ";
      FieldsToRemove.push_back(wndConvFields[iField_Conv]);
    }
  }
  if (removedField && rank ==MASTER_NODE) cout << endl;

  for (unsigned short iField_Conv = 0; iField_Conv < FieldsToRemove.size(); iField_Conv++){
    wndConvFields.erase(std::find(wndConvFields.begin(), wndConvFields.end(), FieldsToRemove[iField_Conv]));
  }
  if (rank == MASTER_NODE){
    if(wndConvFields.empty()){
      cout << "Warning: No (valid) fields chosen for time convergence monitoring. Time convergence monitoring inactive."<<  endl;
    }
    else{
      cout <<"Time Convergence field(s): ";
      for (unsigned short iField_Conv = 0; iField_Conv < wndConvFields.size(); iField_Conv++){
        cout << wndConvFields[iField_Conv];
        if (iField_Conv != wndConvFields.size() - 1) cout << ", ";
      }
      cout << endl;
    }
  }
}

void COutput::PreprocessVolumeOutput(CConfig *config){

  /*--- Set the volume output fields using a virtual function call to the child implementation ---*/

  SetVolumeOutputFields(config);

  /*---Coordinates and solution groups must be always in the output.
   * If they are not requested, add them here. ---*/

  vector<string>::iterator itCoord = std::find(requestedVolumeFields.begin(),
                                          requestedVolumeFields.end(), "COORDINATES");
  if (itCoord == requestedVolumeFields.end()){
    requestedVolumeFields.emplace_back("COORDINATES");
    nRequestedVolumeFields++;
  }
  vector<string>::iterator itSol = std::find(requestedVolumeFields.begin(),
                                          requestedVolumeFields.end(), "SOLUTION");
  if (itSol == requestedVolumeFields.end()){
    requestedVolumeFields.emplace_back("SOLUTION");
    nRequestedVolumeFields++;
  }

  nVolumeFields = 0;

  string RequestedField;
  std::vector<bool> FoundField(nRequestedVolumeFields, false);
  vector<string> FieldsToRemove;


  /*--- Loop through all fields defined in the corresponding SetVolumeOutputFields().
 * If it is also defined in the config (either as part of a group or a single field), the field
 * object gets an offset so that we know where to find the data in the Local_Data() array.
 *  Note that the default offset is -1. An index !=-1 defines this field as part of the output. ---*/

  for (unsigned short iField_Output = 0; iField_Output < volumeOutput_List.size(); iField_Output++){

    const string &fieldReference = volumeOutput_List[iField_Output];
    if (volumeOutput_Map.count(fieldReference) > 0){
      VolumeOutputField &Field = volumeOutput_Map.at(fieldReference);

      /*--- Loop through all fields specified in the config ---*/

      for (unsigned short iReqField = 0; iReqField < nRequestedVolumeFields; iReqField++){

        RequestedField = requestedVolumeFields[iReqField];

        if (((RequestedField == Field.outputGroup) || (RequestedField == fieldReference)) && (Field.offset == -1)){
          Field.offset = nVolumeFields;
          volumeFieldNames.push_back(Field.fieldName);
          nVolumeFields++;

          FoundField[iReqField] = true;
        }
      }
    }
  }

  for (unsigned short iReqField = 0; iReqField < nRequestedVolumeFields; iReqField++){
    if (!FoundField[iReqField]){
      FieldsToRemove.push_back(requestedVolumeFields[iReqField]);
    }
  }

  /*--- Remove fields which are not defined --- */

  for (unsigned short iReqField = 0; iReqField < FieldsToRemove.size(); iReqField++){
    if (rank == MASTER_NODE) {
      if (iReqField == 0){
        cout << "  Info: Ignoring the following volume output fields/groups:" << endl;
        cout << "  ";
      }
      cout << FieldsToRemove[iReqField];
      if (iReqField != FieldsToRemove.size()-1){
        cout << ", ";
      } else {
        cout << endl;
      }
    }
    requestedVolumeFields.erase(std::find(requestedVolumeFields.begin(),
                                          requestedVolumeFields.end(), FieldsToRemove[iReqField]));
  }

  nRequestedVolumeFields = requestedVolumeFields.size();

  if (rank == MASTER_NODE){
    cout <<"Volume output fields: ";
    for (unsigned short iReqField = 0; iReqField < nRequestedVolumeFields; iReqField++){
      RequestedField = requestedVolumeFields[iReqField];
      cout << requestedVolumeFields[iReqField];
      if (iReqField != nRequestedVolumeFields - 1) cout << ", ";
    }
    cout << endl;
  }
}

void COutput::LoadDataIntoSorter(CConfig* config, CGeometry* geometry, CSolver** solver){

  unsigned short iMarker = 0;
  unsigned long iPoint = 0, jPoint = 0;
  unsigned long iVertex = 0;

  /*--- Reset the offset cache and index --- */
  cachePosition = 0;
  fieldIndexCache.clear();
  curGetFieldIndex = 0;
  fieldGetIndexCache.clear();

  if (femOutput){

    /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
     geometrical information for the FEM DG solver. ---*/

    CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

    unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();

    CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

    /*--- Access the solution by looping over the owned volume elements. ---*/

    for(unsigned long l=0; l<nVolElemOwned; ++l) {

      for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {

        buildFieldIndexCache = fieldIndexCache.empty();

        LoadVolumeDataFEM(config, geometry, solver, l, jPoint, j);

        jPoint++;

      }
    }

  } else {

    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

      /*--- Load the volume data into the data sorter. --- */

      buildFieldIndexCache = fieldIndexCache.empty();

      LoadVolumeData(config, geometry, solver, iPoint);

    }

    /*--- Reset the offset cache and index --- */
    cachePosition = 0;
    fieldIndexCache.clear();
    curGetFieldIndex = 0;
    fieldGetIndexCache.clear();

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

      /*--- We only want to have surface values on solid walls ---*/

      if (config->GetSolid_Wall(iMarker)){
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          /*--- Load the surface data into the data sorter. --- */

          if(geometry->node[iPoint]->GetDomain()){

            buildFieldIndexCache = fieldIndexCache.empty();

            LoadSurfaceData(config, geometry, solver, iPoint, iMarker, iVertex);

          }
        }
      }
    }
  }
}

void COutput::SetVolumeOutputValue(string name, unsigned long iPoint, su2double value){

  if (buildFieldIndexCache){

    /*--- Build up the offset cache to speed up subsequent
     * calls of this routine since the order of calls is
     * the same for every value of iPoint --- */

    if (volumeOutput_Map.count(name) > 0){
      const short Offset = volumeOutput_Map.at(name).offset;
      fieldIndexCache.push_back(Offset);
      if (Offset != -1){
        volumeDataSorter->SetUnsorted_Data(iPoint, Offset, value);
      }
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);
    }
  } else {

    /*--- Use the offset cache for the access ---*/

    const short Offset = fieldIndexCache[cachePosition++];
    if (Offset != -1){
      volumeDataSorter->SetUnsorted_Data(iPoint, Offset, value);
    }
    if (cachePosition == fieldIndexCache.size()){
      cachePosition = 0;
    }
  }

}

su2double COutput::GetVolumeOutputValue(string name, unsigned long iPoint){

  if (buildFieldIndexCache){

    /*--- Build up the offset cache to speed up subsequent
     * calls of this routine since the order of calls is
     * the same for every value of iPoint --- */

    if (volumeOutput_Map.count(name) > 0){
      const short Offset = volumeOutput_Map.at(name).offset;
      fieldGetIndexCache.push_back(Offset);
      if (Offset != -1){
        return volumeDataSorter->GetUnsorted_Data(iPoint, Offset);
      }
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);
    }
  } else {

    /*--- Use the offset cache for the access ---*/

    const short Offset = fieldGetIndexCache[curGetFieldIndex++];

    if (curGetFieldIndex == fieldGetIndexCache.size()){
      curGetFieldIndex = 0;
    }
    if (Offset != -1){
      return volumeDataSorter->GetUnsorted_Data(iPoint, Offset);
    }
  }

  return 0.0;
}

void COutput::SetAvgVolumeOutputValue(string name, unsigned long iPoint, su2double value){

  const su2double scaling = 1.0 / su2double(curAbsTimeIter + 1);

  if (buildFieldIndexCache){

    /*--- Build up the offset cache to speed up subsequent
     * calls of this routine since the order of calls is
     * the same for every value of iPoint --- */

    if (volumeOutput_Map.count(name) > 0){
      const short Offset = volumeOutput_Map.at(name).offset;
      fieldIndexCache.push_back(Offset);
      if (Offset != -1){

        const su2double old_value = volumeDataSorter->GetUnsorted_Data(iPoint, Offset);
        const su2double new_value = value * scaling + old_value *( 1.0 - scaling);

        volumeDataSorter->SetUnsorted_Data(iPoint, Offset, new_value);
      }
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);
    }
  } else {

    /*--- Use the offset cache for the access ---*/

    const short Offset = fieldIndexCache[cachePosition++];
    if (Offset != -1){

      const su2double old_value = volumeDataSorter->GetUnsorted_Data(iPoint, Offset);
      const su2double new_value = value * scaling + old_value *( 1.0 - scaling);

      volumeDataSorter->SetUnsorted_Data(iPoint, Offset, new_value);
    }
    if (cachePosition == fieldIndexCache.size()){
      cachePosition = 0;
    }
  }

}





void COutput::Postprocess_HistoryData(CConfig *config){

  map<string, pair<su2double, int> > Average;
  map<string, int> Count;

  for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
    const string &fieldIdentifier = historyOutput_List[iField];
    const HistoryOutputField &currentField = historyOutput_Map.at(fieldIdentifier);
    if (currentField.fieldType == HistoryFieldType::RESIDUAL){
      if ( SetInit_Residuals(config) || (currentField.value > initialResiduals[fieldIdentifier]) ) {
        initialResiduals[fieldIdentifier] = currentField.value;
      }
      SetHistoryOutputValue("REL_" + fieldIdentifier,
                            currentField.value - initialResiduals[fieldIdentifier]);

      Average[currentField.outputGroup].first += currentField.value;
      Average[currentField.outputGroup].second++;

    }

    if (currentField.fieldType == HistoryFieldType::COEFFICIENT){
      if(SetUpdate_Averages(config)){
        if (config->GetTime_Domain()){
          windowedTimeAverages[historyOutput_List[iField]].addValue(currentField.value,config->GetTimeIter(), config->GetStartWindowIteration()); //Collecting Values for Windowing
          SetHistoryOutputValue("TAVG_" + fieldIdentifier, windowedTimeAverages[fieldIdentifier].WindowedUpdate(config->GetKindWindow()));
          if (config->GetDirectDiff() != NO_DERIVATIVE) {
            SetHistoryOutputValue("D_TAVG_" + fieldIdentifier, SU2_TYPE::GetDerivative(windowedTimeAverages[fieldIdentifier].GetVal()));
          }
        }
      }
      if (config->GetDirectDiff() != NO_DERIVATIVE){
        SetHistoryOutputValue("D_" + fieldIdentifier, SU2_TYPE::GetDerivative(currentField.value));
      }
    }
  }

  map<string, pair<su2double, int> >::iterator it = Average.begin();
  for (it = Average.begin(); it != Average.end(); it++){
    const su2double& value = it->second.first;
    const int& count = it->second.second;
    const su2double average = value/count;
    if (historyOutput_Map.count("AVG_" + it->first) > 0 )
      SetHistoryOutputValue("AVG_" + it->first, average);
  }
}

void COutput::Postprocess_HistoryFields(CConfig *config){

  map<string, bool> Average;
  map<string, string> AverageGroupName = {{"BGS_RES", "bgs"},{"RMS_RES","rms"},{"MAX_RES", "max"}};

  for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
    const string &fieldIdentifier = historyOutput_List[iField];
    const HistoryOutputField &currentField = historyOutput_Map.at(fieldIdentifier);
    if (currentField.fieldType == HistoryFieldType::RESIDUAL){
      AddHistoryOutput("REL_" + fieldIdentifier, "rel" + currentField.fieldName, currentField.screenFormat,
                       "REL_" + currentField.outputGroup,  "Relative residual.", HistoryFieldType::AUTO_RESIDUAL);
      Average[currentField.outputGroup] = true;
    }
  }

  map<string, bool>::iterator it = Average.begin();
  for (it = Average.begin(); it != Average.end(); it++){
    if (AverageGroupName.count(it->first) > 0) {
      AddHistoryOutput("AVG_" + it->first, "avg[" + AverageGroupName[it->first] + "]", ScreenOutputFormat::FIXED,
          "AVG_" + it->first , "Average residual over all solution variables.", HistoryFieldType::AUTO_RESIDUAL);
    }
  }

  if (config->GetTime_Domain()){
    for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
      const string &fieldIdentifier = historyOutput_List[iField];
      const HistoryOutputField &currentField = historyOutput_Map.at(fieldIdentifier);
      if (currentField.fieldType == HistoryFieldType::COEFFICIENT){
        AddHistoryOutput("TAVG_"   + fieldIdentifier, "tavg["  + currentField.fieldName + "]",
                         currentField.screenFormat, "TAVG_"   + currentField.outputGroup, "Time averaged values.",
                         HistoryFieldType::AUTO_COEFFICIENT);
      }
    }
  }

  if (config->GetDirectDiff()){
    for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
      const string &fieldIdentifier = historyOutput_List[iField];
      const HistoryOutputField &currentField = historyOutput_Map.at(fieldIdentifier);
      if (currentField.fieldType == HistoryFieldType::COEFFICIENT){
        AddHistoryOutput("D_"      + fieldIdentifier, "d["     + currentField.fieldName + "]",
                         currentField.screenFormat, "D_"      + currentField.outputGroup,
                         "Derivative value (DIRECT_DIFF=YES)", HistoryFieldType::AUTO_COEFFICIENT);
      }
    }
  }

  if (config->GetTime_Domain() && config->GetDirectDiff()){
    for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){
      const string &fieldIdentifier = historyOutput_List[iField];
      const HistoryOutputField &currentField = historyOutput_Map.at(fieldIdentifier);
      if (currentField.fieldType == HistoryFieldType::COEFFICIENT){
        AddHistoryOutput("D_TAVG_" + fieldIdentifier, "dtavg[" + currentField.fieldName + "]",
                         currentField.screenFormat, "D_TAVG_" + currentField.outputGroup,
                         "Derivative of the time averaged value (DIRECT_DIFF=YES)", HistoryFieldType::AUTO_COEFFICIENT);
      }
    }
  }

  for (unsigned short iFieldConv = 0; iFieldConv < convFields.size(); iFieldConv++){
    const string &convField = convFields[iFieldConv];
    if (historyOutput_Map.count(convField) > 0){
      if (historyOutput_Map[convField].fieldType == HistoryFieldType::COEFFICIENT){
        AddHistoryOutput("CAUCHY_" + convField, "Cauchy["  + historyOutput_Map.at(convField).fieldName + "]",
                         ScreenOutputFormat::SCIENTIFIC, "CAUCHY", "Cauchy residual value of field set with CONV_FIELD.",
                         HistoryFieldType::AUTO_COEFFICIENT);
      }
    }
  }
  for (unsigned short iFieldConv = 0; iFieldConv < wndConvFields.size(); iFieldConv++){
    const string &wndConvField = wndConvFields[iFieldConv];
    if (historyOutput_Map.count(wndConvField) > 0){
      AddHistoryOutput("CAUCHY_" + wndConvField, "Cauchy[" + historyOutput_Map[wndConvField].fieldName  + "]", ScreenOutputFormat::SCIENTIFIC, "CAUCHY", "Cauchy residual value of field set with WND_CONV_FIELD.", HistoryFieldType::AUTO_COEFFICIENT);
    }
  }
}

bool COutput::WriteScreen_Header(CConfig *config) {

  unsigned long RestartIter = 0;

  if (config->GetRestart() && config->GetTime_Domain()){
    RestartIter = config->GetRestart_Iter();
  }

  unsigned long ScreenWrt_Freq_Inner = config->GetScreen_Wrt_Freq(2);
  unsigned long ScreenWrt_Freq_Outer = config->GetScreen_Wrt_Freq(1);
  unsigned long ScreenWrt_Freq_Time  = config->GetScreen_Wrt_Freq(0);

  /*--- Header is always disabled for multizone problems unless explicitely requested --- */

  if (config->GetMultizone_Problem() && !config->GetWrt_ZoneConv()){
    return false;
  }

  /*--- Always print header if it is forced ---*/

  if (headerNeeded){
    headerNeeded = false;
    return true;
  }

  /* --- Always print header in the first iteration --- */

  if ((curInnerIter == 0) &&
      (curOuterIter == 0) &&
      (curTimeIter == RestartIter)){
    return true;
  }

  if (!PrintOutput(curTimeIter, ScreenWrt_Freq_Time)&&
      !(curTimeIter == config->GetnTime_Iter() - 1)){
    return false;
  }

  /*--- If there is no inner or outer iteration, don't print header ---*/
  if (ScreenWrt_Freq_Outer == 0 && ScreenWrt_Freq_Inner == 0){
    return false;
  }

  /*--- Print header if we are at the first inner iteration ---*/

  if (curInnerIter == 0){
    return true;
  }

  return false;
}

bool COutput::WriteScreen_Output(CConfig *config) {

  unsigned long ScreenWrt_Freq_Inner = config->GetScreen_Wrt_Freq(2);
  unsigned long ScreenWrt_Freq_Outer = config->GetScreen_Wrt_Freq(1);
  unsigned long ScreenWrt_Freq_Time  = config->GetScreen_Wrt_Freq(0);

  if (config->GetMultizone_Problem() && !config->GetWrt_ZoneConv()){

    return false;

  }

  /*--- Check if screen output should be written --- */

  if (!PrintOutput(curTimeIter, ScreenWrt_Freq_Time)&&
      !(curTimeIter == config->GetnTime_Iter() - 1)){

    return false;

  }

  if (convergence) {return true;}

  if (!PrintOutput(curOuterIter, ScreenWrt_Freq_Outer) &&
      !(curOuterIter == config->GetnOuter_Iter() - 1)){

    return false;

  }

  if (!PrintOutput(curInnerIter, ScreenWrt_Freq_Inner) &&
      !(curInnerIter == config->GetnInner_Iter() - 1)){

    return false;

  }

  return true;

}

bool COutput::WriteHistoryFile_Output(CConfig *config) {

  unsigned long HistoryWrt_Freq_Inner = config->GetHistory_Wrt_Freq(2);
  unsigned long HistoryWrt_Freq_Outer = config->GetHistory_Wrt_Freq(1);
  unsigned long HistoryWrt_Freq_Time  = config->GetHistory_Wrt_Freq(0);

  if (config->GetMultizone_Problem() && !config->GetWrt_ZoneHist()){

    return false;

  }

  /*--- Check if screen output should be written --- */

  if (!PrintOutput(curTimeIter, HistoryWrt_Freq_Time)&&
      !(curTimeIter == config->GetnTime_Iter() - 1)){

    return false;

  }

  if (convergence) {return true;}

  if (!PrintOutput(curOuterIter,HistoryWrt_Freq_Outer) &&
      !(curOuterIter == config->GetnOuter_Iter() - 1)){

    return false;

  }

  if (!PrintOutput(curInnerIter, HistoryWrt_Freq_Inner) &&
      !(curInnerIter == config->GetnInner_Iter() - 1)){

    return false;

  }

  return true;

}

bool COutput::WriteVolume_Output(CConfig *config, unsigned long Iter, bool force_writing){
  if (config->GetTime_Domain()) return ((Iter % config->GetVolume_Wrt_Freq() == 0)) || force_writing;
  else {
    return ((Iter > 0) && (Iter % config->GetVolume_Wrt_Freq() == 0)) || force_writing;
  }
}

void COutput::SetCommonHistoryFields(CConfig *config){

  /// BEGIN_GROUP: ITERATION, DESCRIPTION: Iteration identifier.
  /// DESCRIPTION: The time iteration index.
  AddHistoryOutput("TIME_ITER",     "Time_Iter",  ScreenOutputFormat::INTEGER, "ITER", "Time iteration index");
  /// DESCRIPTION: The outer iteration index.
  AddHistoryOutput("OUTER_ITER",   "Outer_Iter",  ScreenOutputFormat::INTEGER, "ITER", "Outer iteration index");
  /// DESCRIPTION: The inner iteration index.
  AddHistoryOutput("INNER_ITER",   "Inner_Iter", ScreenOutputFormat::INTEGER,  "ITER", "Inner iteration index");
  /// END_GROUP

  /// BEGIN_GROUP: TIME_DOMAIN, DESCRIPTION: Time integration information
  /// Description: The current time
  AddHistoryOutput("CUR_TIME", "Cur_Time", ScreenOutputFormat::SCIENTIFIC, "TIME_DOMAIN", "Current physical time (s)");
  /// Description: The current time step
  AddHistoryOutput("TIME_STEP", "Time_Step", ScreenOutputFormat::SCIENTIFIC, "TIME_DOMAIN", "Current time step (s)");

  /// DESCRIPTION: Currently used wall-clock time.
  AddHistoryOutput("WALL_TIME",   "Time(sec)", ScreenOutputFormat::SCIENTIFIC, "WALL_TIME", "Average wall-clock time");

  AddHistoryOutput("NONPHYSICAL_POINTS", "Nonphysical_Points", ScreenOutputFormat::INTEGER, "NONPHYSICAL_POINTS", "The number of non-physical points in the solution");
}

void COutput::LoadCommonHistoryData(CConfig *config){

  SetHistoryOutputValue("TIME_STEP", config->GetDelta_UnstTimeND()*config->GetTime_Ref());

  /*--- Update the current time only if the time iteration has changed ---*/

  if (SU2_TYPE::Int(GetHistoryFieldValue("TIME_ITER")) != static_cast<int>(curTimeIter)) {
    SetHistoryOutputValue("CUR_TIME",  GetHistoryFieldValue("CUR_TIME") + GetHistoryFieldValue("TIME_STEP"));
  }

  SetHistoryOutputValue("TIME_ITER",  curTimeIter);
  SetHistoryOutputValue("INNER_ITER", curInnerIter);
  SetHistoryOutputValue("OUTER_ITER", curOuterIter);

  su2double StopTime, UsedTime;

  StopTime = SU2_MPI::Wtime();

  UsedTime = (StopTime - config->Get_StartTime())/((curOuterIter + 1) * (curInnerIter+1));

  SetHistoryOutputValue("WALL_TIME", UsedTime);

  SetHistoryOutputValue("NONPHYSICAL_POINTS", config->GetNonphysical_Points());
}


void COutput::PrintHistoryFields(){

  if (rank == MASTER_NODE){

    PrintingToolbox::CTablePrinter HistoryFieldTable(&std::cout);

    unsigned short NameSize = 0, GroupSize = 0, DescrSize = 0;

    for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){

      HistoryOutputField &Field = historyOutput_Map.at(historyOutput_List[iField]);

      if (Field.description != ""){
        if (historyOutput_List[iField].size() > NameSize){
          NameSize = historyOutput_List[iField].size();
        }
        if (Field.outputGroup.size() > GroupSize){
          GroupSize = Field.outputGroup.size();
        }
        if (Field.description.size() > DescrSize){
          DescrSize = Field.description.size();
        }
      }
    }

    cout << "Available screen/history output fields for the current configuration in " << multiZoneHeaderString << ":" << endl;

    HistoryFieldTable.AddColumn("Name", NameSize);
    HistoryFieldTable.AddColumn("Group Name", GroupSize);
    HistoryFieldTable.AddColumn("Type",5);
    HistoryFieldTable.AddColumn("Description", DescrSize);
    HistoryFieldTable.SetAlign(PrintingToolbox::CTablePrinter::LEFT);

    HistoryFieldTable.PrintHeader();

    for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){

      HistoryOutputField &Field = historyOutput_Map.at(historyOutput_List[iField]);

      if (Field.fieldType == HistoryFieldType::DEFAULT
          || Field.fieldType == HistoryFieldType::COEFFICIENT
          || Field.fieldType == HistoryFieldType::RESIDUAL){
        string type;
        switch (Field.fieldType) {
          case HistoryFieldType::COEFFICIENT:
            type = "C";
            break;
          case HistoryFieldType::RESIDUAL:
            type = "R";
            break;
          default:
            type = "D";
            break;
        }

        if (Field.description != "")
          HistoryFieldTable << historyOutput_List[iField] << Field.outputGroup << type << Field.description;

      }
    }

    HistoryFieldTable.PrintFooter();

    cout << "Type legend: Default (D), Residual (R), Coefficient (C)" << endl;

    cout << "Generated screen/history fields (only first field of every group is shown):" << endl;

    PrintingToolbox::CTablePrinter ModifierTable(&std::cout);

    ModifierTable.AddColumn("Name", NameSize);
    ModifierTable.AddColumn("Group Name", GroupSize);
    ModifierTable.AddColumn("Type",5);
    ModifierTable.AddColumn("Description", DescrSize);
    ModifierTable.SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    ModifierTable.PrintHeader();

    std::map<string, bool> GroupVisited;

    for (unsigned short iField = 0; iField < historyOutput_List.size(); iField++){

      HistoryOutputField &Field = historyOutput_Map.at(historyOutput_List[iField]);

      if ((Field.fieldType == HistoryFieldType::AUTO_COEFFICIENT ||
           Field.fieldType == HistoryFieldType::AUTO_RESIDUAL) && (GroupVisited.count(Field.outputGroup) == 0)){
        string type;
        switch (Field.fieldType) {
          case HistoryFieldType::AUTO_COEFFICIENT:
            type = "AC";
            break;
          case HistoryFieldType::AUTO_RESIDUAL:
            type = "AR";
            break;
          default:
            type = "AD";
            break;
        }

        if (Field.description != "")
          ModifierTable << historyOutput_List[iField] << Field.outputGroup << type << Field.description;

        GroupVisited[Field.outputGroup] = true;
      }
    }
    ModifierTable.PrintFooter();

  }
}

void COutput::PrintVolumeFields(){

  if (rank == MASTER_NODE){

    PrintingToolbox::CTablePrinter VolumeFieldTable(&std::cout);

    unsigned short NameSize = 0, GroupSize = 0, DescrSize = 0;

    for (unsigned short iField = 0; iField < volumeOutput_List.size(); iField++){

      VolumeOutputField &Field = volumeOutput_Map.at(volumeOutput_List[iField]);

      if (Field.description != ""){
        if (volumeOutput_List[iField].size() > NameSize){
          NameSize = volumeOutput_List[iField].size();
        }
        if (Field.outputGroup.size() > GroupSize){
          GroupSize = Field.outputGroup.size();
        }
        if (Field.description.size() > DescrSize){
          DescrSize = Field.description.size();
        }
      }
    }

    cout << "Available volume output fields for the current configuration in " << multiZoneHeaderString << ":" << endl;
    cout << "Note: COORDINATES and SOLUTION groups are always in the volume output." << endl;
    VolumeFieldTable.AddColumn("Name", NameSize);
    VolumeFieldTable.AddColumn("Group Name", GroupSize);
    VolumeFieldTable.AddColumn("Description", DescrSize);
    VolumeFieldTable.SetAlign(PrintingToolbox::CTablePrinter::LEFT);

    VolumeFieldTable.PrintHeader();

    for (unsigned short iField = 0; iField < volumeOutput_List.size(); iField++){

      VolumeOutputField &Field = volumeOutput_Map.at(volumeOutput_List[iField]);

      if (Field.description != "")
        VolumeFieldTable << volumeOutput_List[iField] << Field.outputGroup << Field.description;

    }

    VolumeFieldTable.PrintFooter();
  }
}
