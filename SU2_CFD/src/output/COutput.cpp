/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information
 * \author F. Palacios, T. Economon
 * \version 7.0.3 "Blackbird"
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

#include <regex>

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


  historyFieldsAll.SetScope(config->GetGlobalScope());
  volumeFieldsAll.SetScope(config->GetGlobalScope());

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
                                                curTimeIter, historyFieldsAll["TIME_STEP"].value);

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
                                          curTimeIter, historyFieldsAll["TIME_STEP"].value);

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

        fileWriter = new CParaviewVTMFileWriter(fileName, fileName, historyFieldsAll["CUR_TIME"].value,
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
                                          curTimeIter, historyFieldsAll["TIME_STEP"].value);

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
                                                curTimeIter, historyFieldsAll["TIME_STEP"].value);

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

  const auto& convFieldRef = historyFieldsAll.GetFieldsByKey(convFields);

  for (const auto& field : HistoryOutFieldCollection::GetFieldsByType({FieldType::COEFFICIENT}, convFieldRef)){
    ConvSummary << field->second.fieldName
                << field->second.value
                << " < " + PrintingToolbox::to_string(cauchyEps)
                << (field->second.value < cauchyEps ? "Yes" : "No");
  }

  for (const auto& field : HistoryOutFieldCollection::GetFieldsByType({FieldType::RESIDUAL,
                                                                      FieldType::AUTO_RESIDUAL}, convFieldRef)){
    ConvSummary << field->second.fieldName
                << field->second.value
                << " < " + PrintingToolbox::to_string(minLogResidual)
                << (field->second.value < minLogResidual ? "Yes" : "No");
  }

  ConvSummary.PrintFooter();
}

bool COutput::Convergence_Monitoring(CConfig *config, unsigned long Iteration) {

  unsigned short iCounter;

  convergence = true;

  for (unsigned short iField_Conv = 0; iField_Conv < convFields.size(); iField_Conv++){

    bool fieldConverged = false;

    const string &convField = convFields[iField_Conv];

    if (GetHistoryFieldsAll().CheckKey(convField)){

      const HistoryOutputField& outField = historyFieldsAll[convField];

      su2double monitor = outField.value;

      /*--- Cauchy based convergence criteria ---*/

      if (outField.fieldType == FieldType::COEFFICIENT) {

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

      if (outField.fieldType == FieldType::RESIDUAL ||
          outField.fieldType == FieldType::AUTO_RESIDUAL) {

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
      if (GetHistoryFieldsAll().CheckKey(WndConv_Field)){
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

      if (GetHistoryFieldsAll().CheckKey(WndConv_Field)){
        su2double monitor = historyFieldsAll[WndConv_Field].value;

        /*--- Cauchy based convergence criteria ---*/

        if (historyFieldsAll[WndConv_Field].fieldType == FieldType::AUTO_COEFFICIENT) { //TAVG values are AUTO_COEFF
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

  stringstream out;
  int width = 20;

  const auto& histFieldsWithName  = GetHistoryFieldsAll().GetFieldsByKey(requestedHistoryFields);
  const auto& histFieldsWithGroup = GetHistoryFieldsAll().GetFieldsByGroup(requestedHistoryFields);
  const auto& histFields          = HistoryOutFieldCollection::Combine(histFieldsWithGroup, histFieldsWithName);

  for (const auto& field : histFields){
    if (field->second.screenFormat == ScreenOutputFormat::INTEGER){
      width = std::max((int)field->second.fieldName.size()+2, 10);
    } else {
      width = std::max((int)field->second.fieldName.size()+2, 18);
    }
    historyFileTable->AddColumn("\"" + field->second.fieldName + "\"", width);
  }

  if (config->GetTabular_FileFormat() == TAB_TECPLOT) {
    histFile << "VARIABLES = \\" << endl;
  }
  historyFileTable->PrintHeader();
  histFile.flush();
}


void COutput::SetHistoryFile_Output(CConfig *config) {

  stringstream out;

  const auto& histFieldsWithName  = GetHistoryFieldsAll().GetFieldsByKey(requestedHistoryFields);
  const auto& histFieldsWithGroup = GetHistoryFieldsAll().GetFieldsByGroup(requestedHistoryFields);
  const auto& histFields          = HistoryOutFieldCollection::Combine(histFieldsWithGroup, histFieldsWithName);

  for (const auto& field : histFields){
    (*historyFileTable) << field->second.value;
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

  const auto& histFieldsWithName  = GetHistoryFieldsAll().GetFieldsByKey(requestedScreenFields);

  for (const auto& field : histFieldsWithName) {
    stringstream out;
    switch (field->second.screenFormat) {
      case ScreenOutputFormat::INTEGER:
        PrintingToolbox::PrintScreenInteger(out, SU2_TYPE::Int(field->second.value), fieldWidth);
        break;
      case ScreenOutputFormat::FIXED:
        PrintingToolbox::PrintScreenFixed(out, field->second.value, fieldWidth);
        break;
      case ScreenOutputFormat::SCIENTIFIC:
        PrintingToolbox::PrintScreenScientific(out, field->second.value, fieldWidth);
        break;
      case ScreenOutputFormat::PERCENT:
        PrintingToolbox::PrintScreenPercent(out, field->second.value, fieldWidth);
        break;
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

  auto printInfo = [&](const std::vector<string>& fields, const std::string& info){
    if (rank == MASTER_NODE){
      for (unsigned int i = 0; i < fields.size(); i++){
        if (i == 0) {
          cout << info << endl;
        }
        cout << fields[i];
        if (i != fields.size()-1){
          cout << ", ";
        } else {
          cout << endl;
        }
      }
    }
  };

  /*--- Set screen convergence output header and remove unavailable fields ---*/

  string requestedField;
  vector<string> FieldsToRemove;
  vector<bool> FoundField(nRequestedHistoryFields, false);
  regex exp("\\{\\S*\\}");
  regex integrateExp("integrate\\{\\S*\\}");

  auto addCustomFields = [&](std::vector<string>& fields){
    /*--- Check if any of the fields is a expression ---*/
    for (const auto& field : fields){
      if (regex_match(field, exp)){
        /*--- Make sure that the expression is not already there ---*/
        if (!historyFieldsAll.CheckKey(field))
          AddCustomHistoryOutput(field);
      }
      if (regex_match(field, integrateExp)){
        if (!historyFieldsAll.CheckKey(field)){
          AddCustomHistoryOutput(field, FieldType::CUSTOM_INTEGRATE);
        }
      }
    }
  };

  addCustomFields(requestedScreenFields);
  addCustomFields(requestedHistoryFields);

  const auto& screenFields = GetHistoryFieldsAll().GetFieldsByKey(requestedScreenFields, FieldsToRemove);

  /*--- Add fields to the screen output table--- */

  for (const auto& field : screenFields){
    convergenceTable->AddColumn(field->second.fieldName, fieldWidth);
  }

  requestedScreenFields = HistoryOutFieldCollection::GetKeys(screenFields);
  nRequestedScreenFields = requestedScreenFields.size();

  printInfo(requestedScreenFields, "Monitoring fields on screen: ");

  /*--- Print info that fields have been removed to screen ---*/

  printInfo(FieldsToRemove, "Fields ignored: ");

  const auto& histFieldsWithName = GetHistoryFieldsAll().GetFieldsByKey(requestedHistoryFields);
  const auto& histFieldsWithGroup = GetHistoryFieldsAll().GetFieldsByGroup(requestedHistoryFields);

  requestedHistoryFields =
      HistoryOutFieldCollection::GetKeys(HistoryOutFieldCollection::Combine(histFieldsWithGroup, histFieldsWithName));
  nRequestedHistoryFields = requestedHistoryFields.size();

  const auto& convergenceFields = GetHistoryFieldsAll().GetFieldsByKey(convFields, FieldsToRemove);

  convFields = HistoryOutFieldCollection::GetKeys(convergenceFields);

  printInfo(convFields, "Convergence monitoring fields: ");
  printInfo(FieldsToRemove, "Fields ignored: ");

  const auto& timeConvergenceFields = GetHistoryFieldsAll().GetFieldsByKey(wndConvFields, FieldsToRemove);

  wndConvFields = HistoryOutFieldCollection::GetKeys(timeConvergenceFields);
  printInfo(wndConvFields, "Time Convergence monitoring fields: ");
  printInfo(FieldsToRemove, "Fields ignored: ");

  historyFieldsAll.UpdateTokens();
  historyFieldsAll.EvalCustomFields(historyFieldsAll.GetFieldsByType({FieldType::CUSTOM_EVAL}));

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

  regex exp("\\{\\S*\\}");
  regex integrateExp("integrate\\{\\S*\\}");

  /*--- Check if any of the fields is a expression ---*/
  for (auto& field : requestedVolumeFields){
    if (regex_match(field, exp)){
      AddCustomVolumeOutput(field);
    }
  }

  for (const auto& field : requestedHistoryFields){
    if (regex_match(field, integrateExp)){
      string customField = field;
      customField.erase(0, 9);
      cout << "Adding " + customField << endl;
      if (!volumeFieldsAll.CheckKey(customField))
        AddCustomVolumeOutput(customField, FieldType::CUSTOM_INTEGRATE);
    }
  }

  for (const auto& field : requestedScreenFields){
    if (regex_match(field, integrateExp)){
      string customField = field;
      customField.erase(0, 9);
      cout << "Adding " + customField << endl;
      AddCustomVolumeOutput(customField, FieldType::CUSTOM_INTEGRATE);
    }
  }

  /*--- Loop through all fields defined in the corresponding SetVolumeOutputFields().
 * If it is also defined in the config (either as part of a group or a single field), the field
 * object gets an offset so that we know where to find the data in the Local_Data() array.
 *  Note that the default offset is -1. An index !=-1 defines this field as part of the output. ---*/

  for (const auto& field : volumeFieldsAll.GetReferencesAll()){
    for (unsigned int iReqField = 0; iReqField < nRequestedVolumeFields; iReqField++){
      const string& RequestedField = requestedVolumeFields[iReqField];
      if (((RequestedField == field->second.outputGroup) || (RequestedField == field->first))
          && (field->second.offset == -1)){
        field->second.offset = nVolumeFields;
        volumeFieldNames.push_back(field->second.fieldName);
        nVolumeFields++;
        FoundField[iReqField] = true;
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

  volumeFieldsAll.UpdateTokens();
  volumeFieldsAll.EvalCustomFields(volumeFieldsAll.GetFieldsByType({FieldType::CUSTOM_EVAL}));
  volumeFieldsAll.IntegrateCustomFields(volumeFieldsAll.GetFieldsByType({FieldType::CUSTOM_INTEGRATE}), 0.0);

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

  const auto& customFieldRef = volumeFieldsAll.GetFieldsByType({FieldType::CUSTOM_EVAL});
  const auto& customIntegration = volumeFieldsAll.GetFieldsByType({FieldType::CUSTOM_INTEGRATE});
  for (const auto& field : customIntegration){
    field->second.value = 0.0;
  }

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

        if (!customFieldRef.empty()){

          volumeFieldsAll.UpdateTokens();

          volumeFieldsAll.EvalCustomFields(customFieldRef);

        }

        WriteToDataSorter(jPoint);

        jPoint++;

      }
    }

  } else {

    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

      /*--- Load the volume data into the data sorter. --- */

      buildFieldIndexCache = fieldIndexCache.empty();

      LoadVolumeData(config, geometry, solver, iPoint);

      if (!customFieldRef.empty() || !customIntegration.empty()){

        volumeFieldsAll.UpdateTokens();

        volumeFieldsAll.EvalCustomFields(customFieldRef);

        volumeFieldsAll.IntegrateCustomFields(customIntegration, geometry->node[iPoint]->GetVolume());

      }

      WriteToDataSorter(iPoint);

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

            LoadVolumeData(config, geometry, solver, iPoint);

            LoadSurfaceData(config, geometry, solver, iPoint, iMarker, iVertex);

            if (!customFieldRef.empty()){

              volumeFieldsAll.UpdateTokens();

              volumeFieldsAll.EvalCustomFields(customFieldRef);

            }
            WriteToDataSorter(iPoint);
          }
        }
      }
    }
  }

  for (const auto& field : customIntegration){
    historyFieldsAll.SetValueByKey("integrate" + field->first, field->second.value);
  }
}

 void COutput::WriteToDataSorter(unsigned long iPoint){
  for (const auto& field : volumeFieldsAll.GetReferencesAll()){
    if (field->second.offset != -1){
      volumeDataSorter->SetUnsorted_Data(iPoint, field->second.offset, field->second.value);
    }
  }
}

void COutput::SetVolumeOutputValue(string name, unsigned long iPoint, su2double value){

  if (buildFieldIndexCache){

    /*--- Build up the offset cache to speed up subsequent
     * calls of this routine since the order of calls is
     * the same for every value of iPoint --- */

    if (volumeFieldsAll.CheckKey(name)){
      fieldIndexCache.push_back(volumeFieldsAll.GetIndex(name));
      volumeFieldsAll.SetValueByIndex(fieldIndexCache.back(), value);
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);
    }
  } else {

    /*--- Use the offset cache for the access ---*/

    const int index = fieldIndexCache[cachePosition++];
    volumeFieldsAll.SetValueByIndex(index, value);
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

    if (volumeFieldsAll.CheckKey(name)){
      fieldGetIndexCache.push_back(volumeFieldsAll.GetIndex(name));
      return volumeFieldsAll.GetItemByIndex(fieldGetIndexCache.back()).value;
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);
    }
  } else {

    /*--- Use the offset cache for the access ---*/

    const short index = fieldGetIndexCache[curGetFieldIndex++];

    if (curGetFieldIndex == fieldGetIndexCache.size()){
      curGetFieldIndex = 0;
    }
    return volumeFieldsAll.GetItemByIndex(index).value;
  }

  return 0.0;
}

void COutput::SetAvgVolumeOutputValue(string name, unsigned long iPoint, su2double value){

  const su2double scaling = 1.0 / su2double(curAbsTimeIter + 1);

  if (buildFieldIndexCache){

    /*--- Build up the offset cache to speed up subsequent
     * calls of this routine since the order of calls is
     * the same for every value of iPoint --- */

    if (volumeFieldsAll.CheckKey(name)){

      const int index = volumeFieldsAll.GetIndex(name);
      fieldIndexCache.push_back(index);
      const su2double old_value = volumeFieldsAll.GetItemByIndex(index).value;
      const su2double new_value = value * scaling + old_value *( 1.0 - scaling);
      volumeFieldsAll.SetValueByIndex(fieldIndexCache.back(), new_value);

    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);
    }
  } else {

    /*--- Use the offset cache for the access ---*/

    const int index = fieldIndexCache[cachePosition++];

    const su2double old_value = volumeFieldsAll.GetItemByIndex(index).value;
    const su2double new_value = value * scaling + old_value *( 1.0 - scaling);
    volumeFieldsAll.SetValueByIndex(index, new_value);

    if (cachePosition == fieldIndexCache.size()){
      cachePosition = 0;
    }
  }

}





void COutput::Postprocess_HistoryData(CConfig *config){

  map<string, pair<su2double, int> > Average;
  map<string, int> Count;

  for (const auto& field : historyFieldsAll.GetFieldsByType({FieldType::RESIDUAL})){
    if (SetInit_Residuals(config) || (field->second.value > initialResiduals[field->first])) {
      initialResiduals[field->first] = field->second.value;
    }
    SetHistoryOutputValue("REL_" + field->first, field->second.value - initialResiduals[field->first]);
    Average[field->second.outputGroup].first += field->second.value;
    Average[field->second.outputGroup].second++;
  }

  for (const auto& field : historyFieldsAll.GetFieldsByType({FieldType::COEFFICIENT})){
    if (SetUpdate_Averages(config)){
      if (config->GetTime_Domain()){
        windowedTimeAverages[field->first].addValue(field->second.value, config->GetTimeIter(), config->GetStartWindowIteration());
        SetHistoryOutputValue("TAVG_" + field->first, windowedTimeAverages[field->first].WindowedUpdate(config->GetKindWindow()));
        if (config->GetDirectDiff()){
          SetHistoryOutputValue("D_TAVG_" + field->first, SU2_TYPE::GetDerivative(windowedTimeAverages[field->first].GetVal()));
        }
      }
    }
    if (config->GetDirectDiff()){
      SetHistoryOutputValue("D_" + field->first, SU2_TYPE::GetDerivative(field->second.value));
    }
  }

  map<string, pair<su2double, int> >::iterator it = Average.begin();
  for (it = Average.begin(); it != Average.end(); it++){
    const su2double& value = it->second.first;
    const int& count = it->second.second;
    const su2double average = value/count;
    if (historyFieldsAll.CheckKey("AVG_" + it->first))
      SetHistoryOutputValue("AVG_" + it->first, average);
  }

  historyFieldsAll.UpdateTokens();
  historyFieldsAll.EvalCustomFields(historyFieldsAll.GetFieldsByType({FieldType::CUSTOM_EVAL}));
}

void COutput::Postprocess_HistoryFields(CConfig *config){

  map<string, bool> Average;
  map<string, string> AverageGroupName = {{"BGS_RES", "bgs"},{"RMS_RES","rms"},{"MAX_RES", "max"}};

  for (const auto& field : historyFieldsAll.GetFieldsByType({FieldType::RESIDUAL})){
    AddHistoryOutput("REL_" + field->first, "rel" + field->second.fieldName, field->second.screenFormat,
                     "REL_" + field->second.outputGroup,  "Relative residual.", FieldType::AUTO_RESIDUAL);
    Average[field->second.outputGroup] = true;
  }

  map<string, bool>::iterator it = Average.begin();
  for (it = Average.begin(); it != Average.end(); it++){
    if (AverageGroupName.count(it->first) > 0) {
      AddHistoryOutput("AVG_" + it->first, "avg[" + AverageGroupName[it->first] + "]", ScreenOutputFormat::FIXED,
          "AVG_" + it->first , "Average residual over all solution variables.", FieldType::AUTO_RESIDUAL);
    }
  }

  const auto& coefficentFields = historyFieldsAll.GetFieldsByType({FieldType::COEFFICIENT});

  if (config->GetTime_Domain()){
    for (auto field : coefficentFields){
      AddHistoryOutput("TAVG_"   + field->first, "tavg["  + field->second.fieldName + "]",
                       field->second.screenFormat, "TAVG_"   + field->second.outputGroup, "Time averaged values.",
                       FieldType::AUTO_COEFFICIENT);
    }
  }

  if (config->GetDirectDiff()){
    for (const auto& field : coefficentFields){
      AddHistoryOutput("D_" + field->first, "d[" + field->second.fieldName + "]",
                       field->second.screenFormat, "D_" + field->second.outputGroup,
                       "Derivative value (DIRECT_DIFF=YES)", FieldType::AUTO_COEFFICIENT);
    }
  }

  if (config->GetTime_Domain() && config->GetDirectDiff()){
    for (const auto& field : coefficentFields){
      AddHistoryOutput("D_TAVG_" + field->first, "dtavg[" + field->second.fieldName + "]",
                       field->second.screenFormat, "D_TAVG_" + field->second.outputGroup,
                       "Derivative of the time averaged value (DIRECT_DIFF=YES)", FieldType::AUTO_COEFFICIENT);
    }
  }

  /*--- Filter convergence fields which are coefficients ---*/

  const auto& convergenceFields = HistoryOutFieldCollection::GetFieldsByType({FieldType::COEFFICIENT},
                                                                       historyFieldsAll.GetFieldsByKey(convFields));
  for (const auto& field : convergenceFields){
    AddHistoryOutput("CAUCHY_" + field->first, "Cauchy["  + field->second.fieldName + "]",
                     ScreenOutputFormat::SCIENTIFIC, "CAUCHY", "Cauchy residual value of field set with CONV_FIELD.",
                     FieldType::AUTO_COEFFICIENT);
  }

  const auto& wndConvergenceFields = historyFieldsAll.GetFieldsByKey(wndConvFields);

  for (const auto& field : wndConvergenceFields){
    AddHistoryOutput("CAUCHY_" + field->first, "Cauchy["  + field->second.fieldName + "]",
                     ScreenOutputFormat::SCIENTIFIC, "CAUCHY", "Cauchy residual value of field set with WND_CONV_FIELD.",
                     FieldType::AUTO_COEFFICIENT);
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

  if (SU2_TYPE::Int(historyFieldsAll["TIME_ITER"].value) != static_cast<int>(curTimeIter)) {
    SetHistoryOutputValue("CUR_TIME",  historyFieldsAll["CUR_TIME"].value + historyFieldsAll["TIME_STEP"].value);
  }

  SetHistoryOutputValue("TIME_ITER",  curTimeIter);
  SetHistoryOutputValue("INNER_ITER", curInnerIter);
  SetHistoryOutputValue("OUTER_ITER", curOuterIter);

  su2double StopTime, UsedTime;
#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif

  UsedTime = (StopTime - config->Get_StartTime())/((curOuterIter + 1) * (curInnerIter+1));

  SetHistoryOutputValue("WALL_TIME", UsedTime);

  SetHistoryOutputValue("NONPHYSICAL_POINTS", config->GetNonphysical_Points());
}


void COutput::PrintHistoryFields(){

  if (rank == MASTER_NODE){

    PrintingToolbox::CTablePrinter HistoryFieldTable(&std::cout);

    unsigned short NameSize = 0, GroupSize = 0, DescrSize = 0;

    for (const auto& field : historyFieldsAll.GetReferencesAll()){
      if (field->second.description != ""){
        if (field->first.size() > NameSize){
          NameSize =field->first.size();
        }
        if (field->second.outputGroup.size() > GroupSize){
          GroupSize = field->second.outputGroup.size();
        }
        if (field->second.description.size() > DescrSize){
          DescrSize = field->second.description.size();
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

    for (const auto& field : historyFieldsAll.GetReferencesAll()){

      if (field->second.fieldType == FieldType::DEFAULT
          || field->second.fieldType == FieldType::COEFFICIENT
          || field->second.fieldType == FieldType::RESIDUAL){
        string type;
        switch (field->second.fieldType) {
          case FieldType::COEFFICIENT:
            type = "C";
            break;
          case FieldType::RESIDUAL:
            type = "R";
            break;
          default:
            type = "D";
            break;
        }

        if (field->second.description != "")
          HistoryFieldTable << field->first << field->second.outputGroup << type << field->second.description;

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

    for (const auto& field : historyFieldsAll.GetReferencesAll()){

      if ((field->second.fieldType == FieldType::AUTO_COEFFICIENT ||
           field->second.fieldType == FieldType::AUTO_RESIDUAL) && (GroupVisited.count(field->second.outputGroup) == 0)){
        string type;
        switch (field->second.fieldType) {
          case FieldType::AUTO_COEFFICIENT:
            type = "AC";
            break;
          case FieldType::AUTO_RESIDUAL:
            type = "AR";
            break;
          default:
            type = "AD";
            break;
        }

        if (field->second.description != "")
          ModifierTable << field->first << field->second.outputGroup << type << field->second.description;

        GroupVisited[field->second.outputGroup] = true;
      }
    }

    ModifierTable.PrintFooter();

  }
}

void COutput::PrintVolumeFields(){

  if (rank == MASTER_NODE){

    PrintingToolbox::CTablePrinter VolumeFieldTable(&std::cout);

    unsigned short NameSize = 0, GroupSize = 0, DescrSize = 0;

    for (const auto& field : volumeFieldsAll.GetReferencesAll()){
      if (field->second.description != ""){
        if (field->first.size() > NameSize){
          NameSize = field->first.size();
        }
        if (field->second.outputGroup.size() > GroupSize){
          GroupSize = field->second.outputGroup.size();
        }
        if (field->second.description.size() > DescrSize){
          DescrSize = field->second.description.size();
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

    for (const auto& field : volumeFieldsAll.GetReferencesAll()){
      if (field->second.description != ""){
        VolumeFieldTable << field->first << field->second.outputGroup << field->second.description;
      }
    }
    VolumeFieldTable.PrintFooter();
  }
}


