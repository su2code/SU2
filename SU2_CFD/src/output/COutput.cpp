/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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
#include "../../include/output/filewriter/CParaviewBinaryFileWriter.hpp"
#include "../../include/output/filewriter/CTecplotFileWriter.hpp"
#include "../../include/output/filewriter/CTecplotBinaryFileWriter.hpp"
#include "../../include/output/filewriter/CCSVFileWriter.hpp"
#include "../../include/output/filewriter/CSU2FileWriter.hpp"
#include "../../include/output/filewriter/CSU2BinaryFileWriter.hpp"
#include "../../include/output/filewriter/CSU2MeshFileWriter.hpp"



COutput::COutput(CConfig *config) {
  
  
  if((!config->GetMultizone_Problem() && !config->GetSinglezone_Driver()) 
     || config->GetBoolTurbomachinery() || config->GetUnsteady_Simulation() == HARMONIC_BALANCE){
    output_legacy = new COutputLegacy(config);
  }
  
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  field_width = 12;
  
  ConvergenceTable = new PrintingToolbox::CTablePrinter(&std::cout);
  MultiZoneHeaderTable = new PrintingToolbox::CTablePrinter(&std::cout);
  
  /*--- Set default filenames ---*/
  
  SurfaceFilename = "surface";
  VolumeFilename  = "volume";
  RestartFilename = "restart";
  
  /*--- Retrieve the history filename ---*/
 
  HistoryFilename = config->GetConv_FileName();
  
  /*--- Append the zone ID ---*/
 
 if(config->GetnZone() > 1){
   HistoryFilename = config->GetMultizone_HistoryFileName(HistoryFilename, config->GetiZone());
 }
 
  HistorySep = ","; 
  
  unsigned short iZone;
 
  
  /*--- Initialize residual ---*/
  
  RhoRes_New = EPS;
  RhoRes_Old = new su2double[config->GetnZone()];
  for (iZone = 0; iZone < config->GetnZone(); iZone++) RhoRes_Old[iZone] = EPS;
  
  nRequestedHistoryFields = config->GetnHistoryOutput();
  for (unsigned short iField = 0; iField < nRequestedHistoryFields; iField++){
    RequestedHistoryFields.push_back(config->GetHistoryOutput_Field(iField));
  }
  
  nRequestedScreenFields = config->GetnScreenOutput();
  for (unsigned short iField = 0; iField < nRequestedScreenFields; iField++){
    RequestedScreenFields.push_back(config->GetScreenOutput_Field(iField));
  }
  
  nRequestedVolumeFields = config->GetnVolumeOutput();
  for (unsigned short iField = 0; iField < nRequestedVolumeFields; iField++){
    RequestedVolumeFields.push_back(config->GetVolumeOutput_Field(iField));
  }
  
  grid_movement = config->GetGrid_Movement(); 
  
  multizone     = config->GetMultizone_Problem();
  
  /*--- Default is not to use the FEM output merging --- */
  
  fem_output = false;

  /*--- Default is to write history to file and screen --- */

  no_writing = false;

  Cauchy_Serie = new su2double[config->GetCauchy_Elems()];

  Conv_Field = config->GetConv_Field();

  Cauchy_Value = 0.0;
  for (unsigned short iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
    Cauchy_Serie[iCounter] = 0.0;

  Convergence = false;
  
  /*--- Initialize all convergence flags to false. ---*/
  
  Convergence        = false;
  Convergence_FSI    = false;
  Convergence_FullMG = false;
  
}

COutput::~COutput(void) {

  if (RhoRes_Old != NULL) delete [] RhoRes_Old;
  
  delete ConvergenceTable;
  delete MultiZoneHeaderTable;

  delete [] Cauchy_Serie;
  
}



void COutput::SetHistory_Output(CGeometry *geometry,
                                  CSolver **solver_container,
                                  CConfig *config,
                                  unsigned long TimeIter,
                                  unsigned long OuterIter,
                                  unsigned long InnerIter) {

  curr_TimeIter  = TimeIter;
  curr_OuterIter = OuterIter;
  curr_InnerIter = InnerIter;
  
  bool write_header, write_history, write_screen;
  
  /*--- Retrieve residual and extra data -----------------------------------------------------------------*/
  
  LoadCommonHistoryData(config);
  
  LoadHistoryData(config, geometry, solver_container);
  
  Convergence_Monitoring(config, curr_InnerIter);
  
  Postprocess_HistoryData(config);
  
  /*--- Output using only the master node ---*/
  
  if (rank == MASTER_NODE && !no_writing) {
    
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
  
  Convergence_Monitoring(config, curr_InnerIter);  
  
  Postprocess_HistoryData(config);

}
void COutput::SetMultizoneHistory_Output(COutput **output, CConfig **config, unsigned long TimeIter, unsigned long OuterIter){
  
  curr_TimeIter  = TimeIter;
  curr_OuterIter = OuterIter;
  

  bool write_header, write_screen, write_history;
  
  /*--- Retrieve residual and extra data -----------------------------------------------------------------*/
  
  LoadMultizoneHistoryData(output, config);
  
  /*--- Output using only the master node ---*/
  
  if (rank == MASTER_NODE && !no_writing) {
    
    /*--- Write the history file ---------------------------------------------------------------------------*/
    write_history = WriteHistoryFile_Output(config[ZONE_0]);
    if (write_history) SetHistoryFile_Output(config[ZONE_0]);
    
    /*--- Write the screen header---------------------------------------------------------------------------*/
    write_header = WriteScreen_Header(config[ZONE_0]);
    if (write_header) SetScreen_Header(config[ZONE_0]);
    
    /*--- Write the screen output---------------------------------------------------------------------------*/
    write_screen = WriteScreen_Output(config[ZONE_0]);
    if (write_screen) SetScreen_Output(config[ZONE_0]);
    
  }
  
}
void COutput::SetCFL_Number(CSolver *****solver_container, CConfig **config, unsigned short val_iZone) {
  
  su2double CFLFactor = 1.0, power = 1.0, CFL = 0.0, CFLMin = 0.0, CFLMax = 0.0, Div = 1.0, Diff = 0.0, MGFactor[100];
  unsigned short iMesh;
  
  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  unsigned long ExtIter = config[val_iZone]->GetExtIter();
  unsigned short nVar = 1;

  bool energy = config[val_iZone]->GetEnergy_Equation();
  bool weakly_coupled_heat = config[val_iZone]->GetWeakly_Coupled_Heat();

  switch( config[val_iZone]->GetKind_Solver()) {
    case EULER : case NAVIER_STOKES : case RANS:
      if (energy) {
        nVar = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetnVar();
        RhoRes_New = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetRes_RMS(nVar-1);
      }
      else if (weakly_coupled_heat) {
        RhoRes_New = solver_container[val_iZone][INST_0][FinestMesh][HEAT_SOL]->GetRes_RMS(0);
      }
      else {
        RhoRes_New = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetRes_RMS(0);
      }
      break;
    case ADJ_EULER : case ADJ_NAVIER_STOKES: case ADJ_RANS:
      RhoRes_New = solver_container[val_iZone][INST_0][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(0);
      break;
    case HEAT_EQUATION_FVM:
      RhoRes_New = solver_container[val_iZone][INST_0][FinestMesh][HEAT_SOL]->GetRes_RMS(0);
      break;
  }
  
  if (RhoRes_New < EPS) RhoRes_New = EPS;
  if (RhoRes_Old[val_iZone] < EPS) RhoRes_Old[val_iZone] = RhoRes_New;

  Div = RhoRes_Old[val_iZone]/RhoRes_New;
  Diff = RhoRes_New-RhoRes_Old[val_iZone];

  /*--- Compute MG factor ---*/

  for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
    if (iMesh == MESH_0) MGFactor[iMesh] = 1.0;
    else MGFactor[iMesh] = MGFactor[iMesh-1] * config[val_iZone]->GetCFL(iMesh)/config[val_iZone]->GetCFL(iMesh-1);
  }

  if (Div < 1.0) power = config[val_iZone]->GetCFL_AdaptParam(0);
  else power = config[val_iZone]->GetCFL_AdaptParam(1);

  /*--- Detect a stall in the residual ---*/

  if ((fabs(Diff) <= RhoRes_New*1E-8) && (ExtIter != 0)) { Div = 0.1; power = config[val_iZone]->GetCFL_AdaptParam(1); }

  CFLMin = config[val_iZone]->GetCFL_AdaptParam(2);
  CFLMax = config[val_iZone]->GetCFL_AdaptParam(3);

  CFLFactor = pow(Div, power);

  for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
    CFL = config[val_iZone]->GetCFL(iMesh);
    CFL *= CFLFactor;

    if ((iMesh == MESH_0) && (CFL <= CFLMin)) {
      for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
        config[val_iZone]->SetCFL(iMesh, 1.001*CFLMin*MGFactor[iMesh]);
      }
      break;
    }
    if ((iMesh == MESH_0) && (CFL >= CFLMax)) {
      for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++)
        config[val_iZone]->SetCFL(iMesh, 0.999*CFLMax*MGFactor[iMesh]);
      break;
    }

    config[val_iZone]->SetCFL(iMesh, CFL);
  }

  switch( config[val_iZone]->GetKind_Solver()) {
  case EULER : case NAVIER_STOKES : case RANS:
    nVar = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetnVar();
    if (energy) RhoRes_Old[val_iZone] = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetRes_RMS(nVar-1);
    else if (weakly_coupled_heat) RhoRes_Old[val_iZone] = solver_container[val_iZone][INST_0][FinestMesh][HEAT_SOL]->GetRes_RMS(0);
    else RhoRes_Old[val_iZone] = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetRes_RMS(0);
    break;
  case ADJ_EULER : case ADJ_NAVIER_STOKES: case ADJ_RANS:
    RhoRes_Old[val_iZone] = solver_container[val_iZone][INST_0][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(0);
    break;
  case HEAT_EQUATION_FVM:
    RhoRes_Old[val_iZone] = solver_container[val_iZone][INST_0][FinestMesh][HEAT_SOL]->GetRes_RMS(0);
    break;
  }
  
}

void COutput::Load_Data(CGeometry *geometry, CConfig *config, CSolver** solver_container){
  
  /*--- Collect that data defined in the subclasses from the different processors ---*/
  
  if (rank == MASTER_NODE)
    cout << endl << "Loading solution output data locally on each rank." << endl;
  
  CollectVolumeData(config, geometry, solver_container);
 
  /*---- Construct a data sorter object to partition and distribute
   *  the local data into linear chunks across the processors ---*/
  
  if (fem_output){
  
    data_sorter = new CFEMDataSorter(config, geometry, GlobalField_Counter, Local_Data);

  }  else {
    
    data_sorter = new CFVMDataSorter(config, geometry, GlobalField_Counter, Local_Data);
    
  }
  
  /*--- Sort the data, needed for volume and surface output ---*/
  
  if (rank == MASTER_NODE) 
    cout << "Sorting output data across all ranks." << endl;
  
  data_sorter->SortOutputData(config, geometry);
  
}


void COutput::DeallocateData_Parallel(){
  
  if (data_sorter != NULL)
    delete data_sorter;
  
  data_sorter = NULL;
  
}


void COutput::SetFileWriter(CConfig *config, CGeometry *geometry, CParallelDataSorter *sorter, CFileWriter *&filewriter, unsigned short format){
  
  /*--- Write files depending on the format --- */
  
  switch (format) {
    
    case CSV:

      sorter->SortConnectivity(config, geometry, true);
      
      if (rank == MASTER_NODE) {
          cout << "Writing CSV file." << endl;     
      }
         
      filewriter = new CCSVFileWriter(Variable_Names, nDim);
          
      break;
    
  case SU2_RESTART_ASCII:
    
    if (rank == MASTER_NODE) {
        cout << "Writing SU2 ASCII restart file." << endl;     
    }
       
    filewriter = new CSU2FileWriter(Variable_Names, nDim);
        
    break;
    
  case SU2_RESTART_BINARY:
    
    if (rank == MASTER_NODE) {
      cout << "Writing SU2 binary restart file." << endl;   
    }

    filewriter = new CSU2BinaryFileWriter(Variable_Names, nDim);
        
    break;
  
  case SU2_MESH:
    
    /*--- Load and sort the output data and connectivity.
     *  Note that the solver container is not need in this case. ---*/
    
    sorter->SortConnectivity(config, geometry, true);
    
    /*--- Set the mesh ASCII format ---*/
    
    if (rank == MASTER_NODE) {
      cout << "Writing SU2 mesh file." << endl;
    }
    
    filewriter = new CSU2MeshFileWriter(Variable_Names, nDim, config->GetiZone(), config->GetnZone());
    
    
    break;    
    
  case TECPLOT_BINARY:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
     sorter->SortConnectivity(config, geometry, false);
    
    /*--- Write tecplot binary ---*/
    
    if (rank == MASTER_NODE) {
      cout << "Writing Tecplot binary file solution file." << endl;
    }
    
    filewriter = new CTecplotBinaryFileWriter(Variable_Names, nDim, curr_TimeIter, config->GetTime_Step());
      
    break;
    
  case TECPLOT:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
    sorter->SortConnectivity(config, geometry, true);
    
    /*--- Write tecplot binary ---*/
    
    if (rank == MASTER_NODE) {
      cout << "Writing Tecplot ASCII file solution file." << endl;
    }
    
    filewriter = new CTecplotFileWriter(Variable_Names, nDim, curr_TimeIter, config->GetTime_Step());

    break;
    
  case PARAVIEW_BINARY:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
    sorter->SortConnectivity(config, geometry, true);
    
    /*--- Write paraview binary ---*/
    if (rank == MASTER_NODE) {
      cout << "Writing Paraview binary file solution file." << endl;
    }
    
    filewriter = new CParaviewBinaryFileWriter(Variable_Names, nDim);
    
    break;
    
  case PARAVIEW:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
    sorter->SortConnectivity(config, geometry, true);
    
    /*--- Write paraview binary ---*/
    if (rank == MASTER_NODE) {
        cout << "Writing Paraview ASCII file volume solution file." << endl;      
    }
    
    filewriter = new CParaviewFileWriter(Variable_Names, nDim);

    break;
    
  default:
    SU2_MPI::Error("Requested volume output format not available.", CURRENT_FUNCTION);
    break;
  } 
}

void COutput::SetSurface_Output(CGeometry *geometry, CConfig *config, unsigned short format){
  
  CParallelDataSorter* surface_sort = NULL;
  CFileWriter*         file_writer = NULL;
  
  if (fem_output){
    surface_sort = new CSurfaceFEMDataSorter(config, GlobalField_Counter, dynamic_cast<CFEMDataSorter*>(data_sorter));
  } else {
    surface_sort = new CSurfaceFVMDataSorter(config, GlobalField_Counter, dynamic_cast<CFVMDataSorter*>(data_sorter));    
  }

  /*--- Set the file writer --- */
  
  SetFileWriter(config, geometry, surface_sort, file_writer, format);
    
  /*--- Sort the surface output data --- */
  
  surface_sort->SortOutputData(config, geometry);
  
  /*--- Write data to file --- */
  
  file_writer->Write_Data(config->GetFilename(SurfaceFilename, ""), surface_sort);
  
  delete file_writer;
  delete surface_sort;
  
}

void COutput::SetVolume_Output(CGeometry *geometry, CConfig *config, unsigned short format){
  
  string FileName = VolumeFilename;
  
  if(format == SU2_RESTART_ASCII || format == SU2_RESTART_BINARY){
    FileName = RestartFilename;
  }
  
  CFileWriter* file_writer = NULL;

  /*--- Set the file writer --- */
  
  SetFileWriter(config, geometry, data_sorter, file_writer, format);
  
  /*--- Write data to file --- */
  
  file_writer->Write_Data(config->GetFilename(FileName, ""), data_sorter);

  delete file_writer;
  
}


bool COutput::Convergence_Monitoring(CConfig *config, unsigned long Iteration) {

  unsigned short iCounter;
  
  bool Already_Converged = Convergence;
  
  su2double monitor = HistoryOutput_Map[Conv_Field].Value;
  
  /*--- Cauchy based convergence criteria ---*/
  
  if (HistoryOutput_Map[Conv_Field].FieldType == TYPE_COEFFICIENT) {
    
    if (Iteration == 0){
      for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++){
        Cauchy_Serie[iCounter] = 0.0;
      }
      New_Func = monitor;
    }
    
    Old_Func = New_Func;
    New_Func = monitor;
    Cauchy_Func = fabs(New_Func - Old_Func);
    
    Cauchy_Serie[Iteration % config->GetCauchy_Elems()] = Cauchy_Func;
    Cauchy_Value = 1.0;
    if (Iteration >= config->GetCauchy_Elems()){     
      Cauchy_Value = 0.0;
      for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
        Cauchy_Value += Cauchy_Serie[iCounter];
    }
    
    if (Cauchy_Value >= config->GetCauchy_Eps()) { Convergence = false; Convergence_FullMG = false; }
    else { Convergence = true; Convergence_FullMG = true; New_Func = 0.0;}
    
    SetHistoryOutputValue("CAUCHY", Cauchy_Value);
 
  }
  
  /*--- Residual based convergence criteria ---*/
  
  if (HistoryOutput_Map[Conv_Field].FieldType == TYPE_RESIDUAL || HistoryOutput_Map[Conv_Field].FieldType == TYPE_REL_RESIDUAL) {
    
    /*--- Check the convergence ---*/
    
    if ((monitor <= config->GetMinLogResidual())) { Convergence = true; Convergence_FullMG = true; }
    else { Convergence = false; Convergence_FullMG = false; }
    
  }
  
  /*--- Do not apply any convergence criteria of the number
     of iterations is less than a particular value ---*/
  
  if (Iteration < config->GetStartConv_Iter()) {
    Convergence = false;
    Convergence_FullMG = false;
  }
  
  if (Already_Converged) { Convergence = true; Convergence_FullMG = true; }
  
  
  /*--- Apply the same convergence criteria to all the processors ---*/
  
#ifdef HAVE_MPI
  
  unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
  sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
  rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;
  
  /*--- Convergence criteria ---*/
  
  sbuf_conv[0] = Convergence;
  SU2_MPI::Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  
  /*-- Compute global convergence criteria in the master node --*/
  
  sbuf_conv[0] = 0;
  if (rank == MASTER_NODE) {
    if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
    else sbuf_conv[0] = 0;
  }
  
  SU2_MPI::Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
  
  if (sbuf_conv[0] == 1) { Convergence = true; Convergence_FullMG = true; }
  else { Convergence = false; Convergence_FullMG = false; }
  
  delete [] sbuf_conv;
  delete [] rbuf_conv;
  
#endif
  
  /*--- Stop the simulation in case a nan appears, do not save the solution ---*/
  
  if (monitor != monitor) {
    SU2_MPI::Error("SU2 has diverged (NaN detected).", CURRENT_FUNCTION);
  }
  
  
  return Convergence;
}

void COutput::SetHistoryFile_Header(CConfig *config) { 

  
  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
      (config->GetOutput_FileFormat() == FIELDVIEW) ||
      (config->GetOutput_FileFormat() == FIELDVIEW_BINARY)) {
    HistFile << "TITLE = \"SU2 Simulation\"" << endl;
    HistFile << "VARIABLES = ";
  }
  
  stringstream out;
  string RequestedField;
  std::vector<bool> found_field(nRequestedHistoryFields, false);
    
  for (unsigned short iField_Output = 0; iField_Output < HistoryOutput_List.size(); iField_Output++){
    HistoryOutputField &Field = HistoryOutput_Map[HistoryOutput_List[iField_Output]];
    for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      RequestedField = RequestedHistoryFields[iReqField];   
      if (RequestedField == Field.OutputGroup || (RequestedField == HistoryOutput_List[iField_Output])){
        found_field[iReqField] = true;        
        out << "\"" << Field.FieldName << "\"" << HistorySep;
      }
    }  
  }
  
  for (unsigned short iField_Output = 0; iField_Output < HistoryOutputPerSurface_List.size(); iField_Output++){
    for (unsigned short iMarker = 0; iMarker < HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField_Output]].size(); iMarker++){  
      HistoryOutputField &Field = HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField_Output]][iMarker];
      for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
        RequestedField = RequestedHistoryFields[iReqField];   
        if (RequestedField == Field.OutputGroup || (RequestedField == HistoryOutputPerSurface_List[iField_Output])){
          found_field[iReqField] = true;          
          out << "\"" << Field.FieldName << "\"" << HistorySep;        
        }
      }
    }
  }
  
  /*--- Print the string to file and remove the last character (a separator) ---*/
  HistFile << out.str().substr(0, out.str().size() - 1);
  HistFile << endl;
  if (config->GetOutput_FileFormat() == TECPLOT ||
      config->GetOutput_FileFormat() == TECPLOT_BINARY ||
      config->GetOutput_FileFormat() == FIELDVIEW ||
      config->GetOutput_FileFormat() == FIELDVIEW_BINARY) {
    HistFile << "ZONE T= \"Convergence history\"" << endl;
  }
  HistFile.flush();
  
}


void COutput::SetHistoryFile_Output(CConfig *config) { 
  
  stringstream out;
  string RequestedField;
  
  
  for (unsigned short iField_Output = 0; iField_Output < HistoryOutput_List.size(); iField_Output++){
    HistoryOutputField &Field = HistoryOutput_Map[HistoryOutput_List[iField_Output]];
    for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      RequestedField = RequestedHistoryFields[iReqField];   
      if (RequestedField == Field.OutputGroup){
        out << std::setprecision(10) << Field.Value << HistorySep << " ";
      }
    }
  }
  
  for (unsigned short iField_Output = 0; iField_Output < HistoryOutputPerSurface_List.size(); iField_Output++){
    for (unsigned short iMarker = 0; iMarker < HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField_Output]].size(); iMarker++){
      HistoryOutputField &Field = HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField_Output]][iMarker];
      for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
        RequestedField = RequestedHistoryFields[iReqField];   
        if (RequestedField == Field.OutputGroup){
          out << std::setprecision(10) << Field.Value << HistorySep << " ";
        }
      }
    }
  }
  
  /*--- Print the string to file and remove the last two characters (a separator and a space) ---*/
  
  HistFile << out.str().substr(0, out.str().size()-2);
  HistFile << endl;
  HistFile.flush();
}

void COutput::SetScreen_Header(CConfig *config) {
  if (config->GetMultizone_Problem()) 
    MultiZoneHeaderTable->PrintHeader();
  ConvergenceTable->PrintHeader();
}


void COutput::SetScreen_Output(CConfig *config) {
  
  string RequestedField;
  
  for (unsigned short iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
    stringstream out;
    RequestedField = RequestedScreenFields[iReqField]; 
    if (HistoryOutput_Map.count(RequestedField) > 0){  
      switch (HistoryOutput_Map[RequestedField].ScreenFormat) {
        case FORMAT_INTEGER:
          PrintingToolbox::PrintScreenInteger(out, SU2_TYPE::Int(HistoryOutput_Map[RequestedField].Value), field_width);
          break;
        case FORMAT_FIXED:
          PrintingToolbox::PrintScreenFixed(out, HistoryOutput_Map[RequestedField].Value, field_width);
          break;
        case FORMAT_SCIENTIFIC:
          PrintingToolbox::PrintScreenScientific(out, HistoryOutput_Map[RequestedField].Value, field_width);
          break;      
      }
    }
    if (HistoryOutputPerSurface_Map.count(RequestedField) > 0){
      switch (HistoryOutputPerSurface_Map[RequestedField][0].ScreenFormat) {
        case FORMAT_INTEGER:
          PrintingToolbox::PrintScreenInteger(out, SU2_TYPE::Int(HistoryOutputPerSurface_Map[RequestedField][0].Value), field_width);
          break;
        case FORMAT_FIXED:
          PrintingToolbox::PrintScreenFixed(out, HistoryOutputPerSurface_Map[RequestedField][0].Value, field_width);
          break;
        case FORMAT_SCIENTIFIC:
          PrintingToolbox::PrintScreenScientific(out, HistoryOutputPerSurface_Map[RequestedField][0].Value, field_width);
          break;   
      }
    }      
    (*ConvergenceTable) << out.str();
  }
}

void COutput::PreprocessHistoryOutput(CConfig *config, bool wrt){
  
    no_writing = !wrt;

    /*--- Set the common output fields ---*/
    
    SetCommonHistoryFields(config);
    
    /*--- Set the History output fields using a virtual function call to the child implementation ---*/
    
    SetHistoryOutputFields(config);
    
    /*--- Postprocess the history fields. Creates new fields based on the ones set in the child classes ---*/
    
    Postprocess_HistoryFields(config);
    
    if (rank == MASTER_NODE && !no_writing){
      
      /*--- Check for consistency and remove fields that are requested but not available --- */
      
      CheckHistoryOutput();
      
      /*--- Open history file and print the header ---*/
      
      PrepareHistoryFile(config);
      
      /*--- Set the multizone screen header ---*/
      
      if (config->GetMultizone_Problem()){
        MultiZoneHeaderTable->AddColumn(MultiZoneHeaderString, nRequestedScreenFields*field_width + (nRequestedScreenFields-1));      
        MultiZoneHeaderTable->SetAlign(PrintingToolbox::CTablePrinter::CENTER);
        MultiZoneHeaderTable->SetPrintHeaderBottomLine(false);
      }
      
    }
    
}

void COutput::PreprocessMultizoneHistoryOutput(COutput **output, CConfig **config, bool wrt){
  
  no_writing = !wrt;

  /*--- Set the History output fields using a virtual function call to the child implementation ---*/
  
  SetMultizoneHistoryOutputFields(output, config);
  
  if (rank == MASTER_NODE && !no_writing){
    
    /*--- Postprocess the history fields. Creates new fields based on the ones set in the child classes ---*/
   
    //Postprocess_HistoryFields(config[ZONE_0]);
    
    /*--- Check for consistency and remove fields that are requested but not available --- */
    
    CheckHistoryOutput();
    
    /*--- Open history file and print the header ---*/
    
    PrepareHistoryFile(config[ZONE_0]);
    
    /*--- Set the multizone screen header ---*/

    if (config[ZONE_0]->GetMultizone_Problem()){
      MultiZoneHeaderTable->AddColumn(MultiZoneHeaderString, nRequestedScreenFields*field_width + (nRequestedScreenFields-1));      
      MultiZoneHeaderTable->SetAlign(PrintingToolbox::CTablePrinter::CENTER);
      MultiZoneHeaderTable->SetPrintHeaderBottomLine(false);
    }
    
  }
  
}

void COutput::PrepareHistoryFile(CConfig *config){

  char buffer[50];
  
  string history_filename;
  
  strcpy (char_histfile, HistoryFilename.data());
  
   /*--- Append the restart iteration: if dynamic problem and restart ---*/
  
  if (config->GetTime_Domain() && config->GetRestart()) {
    long iExtIter = config->GetRestart_Iter();
    if (SU2_TYPE::Int(iExtIter) < 10) SPRINTF (buffer, "_0000%d", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10) && (SU2_TYPE::Int(iExtIter) < 100)) SPRINTF (buffer, "_000%d", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100) && (SU2_TYPE::Int(iExtIter) < 1000)) SPRINTF (buffer, "_00%d", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d", SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d", SU2_TYPE::Int(iExtIter));
    strcat(char_histfile, buffer);
  }
  
  /*--- Add the correct file extension depending on the file format ---*/
  
  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == FIELDVIEW)) SPRINTF (buffer, ".dat");
  else if ((config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
           (config->GetOutput_FileFormat() == FIELDVIEW_BINARY))  SPRINTF (buffer, ".plt");
  else if (config->GetOutput_FileFormat() == PARAVIEW || config->GetOutput_FileFormat() == PARAVIEW_BINARY)  SPRINTF (buffer, ".csv");
  strcat(char_histfile, buffer);
  
  /*--- Open the history file ---*/
  
  HistFile.open(char_histfile, ios::out);
  HistFile.precision(15);
  
  /*--- Add the header to the history file. ---*/
  
  SetHistoryFile_Header(config);    
  
}

void COutput::CheckHistoryOutput(){
  
  
  /*--- Set screen convergence output header and remove unavailable fields ---*/
  
  string RequestedField;
  vector<string> FieldsToRemove;
  vector<bool> FoundField(nRequestedHistoryFields, false);
  
  for (unsigned short iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
    RequestedField = RequestedScreenFields[iReqField];  
    if (HistoryOutput_Map.count(RequestedField) > 0){ 
      ConvergenceTable->AddColumn(HistoryOutput_Map[RequestedField].FieldName, field_width);
    }
    else if (HistoryOutputPerSurface_Map.count(RequestedField) > 0){
      ConvergenceTable->AddColumn(HistoryOutputPerSurface_Map[RequestedField][0].FieldName, field_width);
    }else {
      FieldsToRemove.push_back(RequestedField);
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
    RequestedScreenFields.erase(std::find(RequestedScreenFields.begin(), RequestedScreenFields.end(), FieldsToRemove[iReqField]));
  }
  
  nRequestedScreenFields = RequestedScreenFields.size();
  
  if (rank == MASTER_NODE){
    cout <<"Screen output fields: ";
    for (unsigned short iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
      RequestedField = RequestedScreenFields[iReqField];            
      cout << RequestedScreenFields[iReqField];
      if (iReqField != nRequestedScreenFields - 1) cout << ", ";
    }
    cout << endl;
  }
  
  /*--- Remove unavailable fields from the history file output ---*/
  
  FieldsToRemove.clear();
  FoundField = vector<bool>(nRequestedHistoryFields, false);
  
  for (unsigned short iField_Output = 0; iField_Output < HistoryOutput_List.size(); iField_Output++){
    HistoryOutputField &Field = HistoryOutput_Map[HistoryOutput_List[iField_Output]];
    for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      RequestedField = RequestedHistoryFields[iReqField];   
      if (RequestedField == Field.OutputGroup){
        FoundField[iReqField] = true;
      }
    }
  }
  
  for (unsigned short iField_Output = 0; iField_Output < HistoryOutputPerSurface_List.size(); iField_Output++){
    for (unsigned short iMarker = 0; iMarker < HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField_Output]].size(); iMarker++){
      HistoryOutputField &Field = HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField_Output]][iMarker];
      for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
        RequestedField = RequestedHistoryFields[iReqField];   
        if (RequestedField == Field.OutputGroup){
          FoundField[iReqField] = true;
        }
      }
    }
  }
  
  for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
    if (!FoundField[iReqField]){
      FieldsToRemove.push_back(RequestedHistoryFields[iReqField]);
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
    RequestedHistoryFields.erase(std::find(RequestedHistoryFields.begin(), RequestedHistoryFields.end(), FieldsToRemove[iReqField]));
  }
  
  nRequestedHistoryFields = RequestedHistoryFields.size();
  
  if (rank == MASTER_NODE){
    cout <<"History output groups: ";
    for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      RequestedField = RequestedHistoryFields[iReqField];            
      cout << RequestedHistoryFields[iReqField];
      if (iReqField != nRequestedHistoryFields - 1) cout << ", ";
    }
    cout << endl;
  }
  
  /*--- Check that the requested convergence monitoring field is available ---*/

  if (HistoryOutput_Map.count(Conv_Field) == 0){
    SU2_MPI::Error(string("Convergence monitoring field ") + Conv_Field + string(" not available"), CURRENT_FUNCTION);
  }
}

void COutput::PreprocessVolumeOutput(CConfig *config, CGeometry *geometry){

//  /*--- Make sure that coordinates are always in the volume output --- */
  
//  if(!(std::find(RequestedVolumeFields.begin(), RequestedVolumeFields.end(), "COORDINATES") != RequestedVolumeFields.end())) {
//    RequestedVolumeFields.push_back("COORDINATES");
//    nRequestedVolumeFields++;
//  }
  
  
  /*--- Set the volume output fields using a virtual function call to the child implementation ---*/  
  
  SetVolumeOutputFields(config);
   
  GlobalField_Counter = 0;
  
  string RequestedField;
  std::vector<bool> FoundField(nRequestedVolumeFields, false);
  vector<string> FieldsToRemove;
  
  
  /*--- Loop through all fields defined in the corresponding SetVolumeOutputFields(). 
 * If it is also defined in the config (either as part of a group or a single field), the field 
 * object gets an offset so that we know where to find the data in the Local_Data() array.
 *  Note that the default offset is -1. An index !=-1 defines this field as part of the output. ---*/

  for (unsigned short iField_Output = 0; iField_Output < VolumeOutput_List.size(); iField_Output++){
    
    VolumeOutputField &Field = VolumeOutput_Map[VolumeOutput_List[iField_Output]];
    
    /*--- Loop through all fields specified in the config ---*/
    
    for (unsigned short iReqField = 0; iReqField < nRequestedVolumeFields; iReqField++){
      
      RequestedField = RequestedVolumeFields[iReqField];  
            
      if (((RequestedField == Field.OutputGroup) || (RequestedField == VolumeOutput_List[iField_Output])) && (Field.Offset == -1)){
        Field.Offset = GlobalField_Counter;
        Variable_Names.push_back(Field.FieldName);
        GlobalField_Counter++;
        
        FoundField[iReqField] = true;
      }
    }    
  }
  
  for (unsigned short iReqField = 0; iReqField < nRequestedVolumeFields; iReqField++){
    if (!FoundField[iReqField]){
      FieldsToRemove.push_back(RequestedVolumeFields[iReqField]);
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
    RequestedVolumeFields.erase(std::find(RequestedVolumeFields.begin(), RequestedVolumeFields.end(), FieldsToRemove[iReqField]));
  }
  
  if (rank == MASTER_NODE){
    cout <<"Volume output fields: ";
    for (unsigned short iReqField = 0; iReqField < nRequestedVolumeFields; iReqField++){
      RequestedField = RequestedVolumeFields[iReqField];            
      cout << RequestedVolumeFields[iReqField];
      if (iReqField != nRequestedVolumeFields - 1) cout << ", ";
    }
    cout << endl;
  }
  
  unsigned long nPoint = 0;
  
  if (fem_output){
    
    /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
    
    CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);
    
    unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
    
    CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();
    
    /*--- Access the solution by looping over the owned volume elements. ---*/
    
    for(unsigned long l=0; l<nVolElemOwned; ++l) {
      
      for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {        
       
        nPoint++;
        
      }
    }
  } else {
    nPoint = geometry->GetnPoint();
  }
  
  /*--- Now that we know the number of fields, create the local data array to temporarily store the volume output 
   * before writing it to file ---*/
   
  Local_Data.resize(nPoint, std::vector<su2double>(GlobalField_Counter, 0.0));
  
}

void COutput::CollectVolumeData(CConfig* config, CGeometry* geometry, CSolver** solver){
  
  bool Wrt_Halo = config->GetWrt_Halo();
  unsigned short iMarker = 0;
  unsigned long iPoint = 0, jPoint = 0;
  long iVertex = 0;
  
  if (fem_output){
    
    /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
     geometrical information for the FEM DG solver. ---*/
  
    CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);
  
    unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
    
    CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();
    
    /*--- Access the solution by looping over the owned volume elements. ---*/
  
    for(unsigned long l=0; l<nVolElemOwned; ++l) {

      for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
        
        LoadVolumeDataFEM(config, geometry, solver, l, jPoint, j);
        
        jPoint++;
        
      }
    }
    
  } else {
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      /*--- Load the volume data into the Local_Data() array. --- */
      
      LoadVolumeData(config, geometry, solver, jPoint);
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_Plotting(iMarker) == YES) {
          iVertex = geometry->node[iPoint]->GetVertex(iMarker);
          if (iVertex != -1){ 
            
            /*--- Load the surface data into the Local_Data() array. --- */
            
            LoadSurfaceData(config, geometry, solver, jPoint, iMarker, iVertex);
          }
        }
      }
      jPoint++;
    }
  }
}

void COutput::Postprocess_HistoryData(CConfig *config){
   
  map<string, su2double> Average;
  map<string, int> Count;
  
  for (unsigned short iField = 0; iField < HistoryOutput_List.size(); iField++){
    HistoryOutputField &currentField = HistoryOutput_Map[HistoryOutput_List[iField]];
    if (currentField.FieldType == TYPE_RESIDUAL){
      if ( SetInit_Residuals(config) || (currentField.Value > Init_Residuals[HistoryOutput_List[iField]]) ) {
        Init_Residuals[HistoryOutput_List[iField]] = currentField.Value;
      }
      SetHistoryOutputValue("REL_" + HistoryOutput_List[iField], currentField.Value - Init_Residuals[HistoryOutput_List[iField]]);
      
      Average[currentField.OutputGroup] += currentField.Value;
      Count[currentField.OutputGroup]++;
           
    }
    if (currentField.FieldType == TYPE_COEFFICIENT){
      if(SetUpdate_Averages(config)){
        SetHistoryOutputValue("TAVG_" + HistoryOutput_List[iField], RunningAverages[HistoryOutput_List[iField]].Update(currentField.Value));
      }
      if (config->GetDirectDiff() != NO_DERIVATIVE){
        SetHistoryOutputValue("D_" + HistoryOutput_List[iField], SU2_TYPE::GetDerivative(currentField.Value));      
        SetHistoryOutputValue("D_TAVG_" + HistoryOutput_List[iField], SU2_TYPE::GetDerivative(RunningAverages[HistoryOutput_List[iField]].Get()));
        
      }
    }
  }

  map<string, su2double>::iterator it = Average.begin();
  for (it = Average.begin(); it != Average.end(); it++){
    SetHistoryOutputValue("AVG_" + it->first, it->second/Count[it->first]);
  }
  
}

void COutput::Postprocess_HistoryFields(CConfig *config){
  
  map<string, bool> Average;
  map<string, string> AverageGroupName =  CCreateMap<string, string>("BGS_RES", "bgs")("RMS_RES","rms")("MAX_RES", "max");
  
  for (unsigned short iField = 0; iField < HistoryOutput_List.size(); iField++){
    HistoryOutputField &currentField = HistoryOutput_Map[HistoryOutput_List[iField]];
    if (currentField.FieldType == TYPE_RESIDUAL){
      AddHistoryOutput("REL_" + HistoryOutput_List[iField], "rel" + currentField.FieldName, currentField.ScreenFormat, "REL_" + currentField.OutputGroup, TYPE_REL_RESIDUAL);
      Average[currentField.OutputGroup] = true;
    }
    if (currentField.FieldType == TYPE_COEFFICIENT){
      AddHistoryOutput("TAVG_"   + HistoryOutput_List[iField], "tavg["  + currentField.FieldName + "]", currentField.ScreenFormat, "TAVG_"   + currentField.OutputGroup);
      AddHistoryOutput("D_"      + HistoryOutput_List[iField], "d["     + currentField.FieldName + "]", currentField.ScreenFormat, "D_"      + currentField.OutputGroup);  
      AddHistoryOutput("D_TAVG_" + HistoryOutput_List[iField], "dtavg[" + currentField.FieldName + "]", currentField.ScreenFormat, "D_TAVG_" + currentField.OutputGroup);  
    }
  }
  
  if (HistoryOutput_Map[Conv_Field].FieldType == TYPE_COEFFICIENT){
    AddHistoryOutput("CAUCHY", "C["  + HistoryOutput_Map[Conv_Field].FieldName + "]", FORMAT_SCIENTIFIC, "RESIDUAL");
  }
  
   map<string, bool>::iterator it = Average.begin();
   for (it = Average.begin(); it != Average.end(); it++){
     AddHistoryOutput("AVG_" + it->first, "avg[" + AverageGroupName[it->first] + "]", FORMAT_FIXED, "AVG_RES");
   }
  
}

bool COutput::WriteScreen_Header(CConfig *config) {  
  bool write_header = false;
  if (config->GetUnsteady_Simulation() == STEADY) {
    write_header = ((curr_InnerIter % (config->GetWrt_Con_Freq()*40)) == 0) || (config->GetMultizone_Problem() && curr_InnerIter == 0);
  } else if  (config->GetUnsteady_Simulation() == TIME_STEPPING) {
    if (!config->GetRestart())
      write_header = ((curr_TimeIter % (config->GetWrt_Con_Freq()*40)) == 0) || (config->GetMultizone_Problem() && curr_InnerIter == 0);    
    else {
      write_header = (((curr_TimeIter - config->GetRestart_Iter()) % (config->GetWrt_Con_Freq()*40)) == 0) || (config->GetMultizone_Problem() && curr_InnerIter == 0);    
    }
  } else {
    write_header = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND) && config->GetInnerIter() == 0;
  }

  /*--- For multizone problems, print the header only if requested explicitly (default of GetWrt_ZoneConv is false) ---*/
  if(config->GetMultizone_Problem()) write_header = (write_header && config->GetWrt_ZoneConv());

  return write_header;
}

bool COutput::WriteScreen_Output(CConfig *config) {
  bool write_output = false;
  
  write_output = config->GetnInner_Iter() - 1 == curr_InnerIter;
  
  if (((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) )){
    write_output = write_output || PrintOutput(config->GetInnerIter(), config->GetWrt_Con_Freq_DualTime());
  }
  else if (((config->GetUnsteady_Simulation() == STEADY) || (config->GetUnsteady_Simulation() == TIME_STEPPING) )){
    write_output = write_output ||  PrintOutput(config->GetInnerIter(), config->GetWrt_Con_Freq()) ;    
  } 

  /*--- For multizone problems, print the body only if requested explicitly (default of GetWrt_ZoneConv is false) ---*/
  if(config->GetMultizone_Problem()) write_output = (write_output && config->GetWrt_ZoneConv());

  return write_output;
}

bool COutput::WriteHistoryFile_Output(CConfig *config) { 
// if (!write_dualtime){
//   return true;
// }
// else {
//   return false;
// }
  return true;
}

void COutput::SetCommonHistoryFields(CConfig *config){
  
  /// BEGIN_GROUP: ITERATION, DESCRIPTION: Iteration identifier.
  /// DESCRIPTION: The time iteration index.
  AddHistoryOutput("TIME_ITER",     "Time_Iter",  FORMAT_INTEGER, "ITER"); 
  /// DESCRIPTION: The outer iteration index.
  AddHistoryOutput("OUTER_ITER",   "Outer_Iter",  FORMAT_INTEGER, "ITER"); 
  /// DESCRIPTION: The inner iteration index.
  AddHistoryOutput("INNER_ITER",   "Inner_Iter", FORMAT_INTEGER,  "ITER"); 
  /// END_GROUP
  
  /// BEGIN_GROUP: TIME_DOMAIN, DESCRIPTION: Time integration information
  /// Description: The current time
  AddHistoryOutput("CUR_TIME", "Cur_Time", FORMAT_FIXED, "TIME_DOMAIN");
  /// Description: The current time step
  AddHistoryOutput("TIME_STEP", "Time_Step", FORMAT_FIXED, "TIME_DOMAIN");
 
  /// DESCRIPTION: Currently used wall-clock time.
  AddHistoryOutput("PHYS_TIME",   "Time(sec)", FORMAT_SCIENTIFIC, "PHYS_TIME"); 
  
}

void COutput::LoadCommonHistoryData(CConfig *config){
  
  SetHistoryOutputValue("TIME_ITER",  curr_TimeIter);  
  SetHistoryOutputValue("INNER_ITER", curr_InnerIter);
  SetHistoryOutputValue("OUTER_ITER", curr_OuterIter); 
  
  SetHistoryOutputValue("CUR_TIME",  curr_TimeIter*config->GetTime_Step());
  SetHistoryOutputValue("TIME_STEP", config->GetTime_Step());
  
  su2double StopTime, UsedTime;
#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  
  UsedTime = (StopTime - config->Get_StartTime())/((curr_OuterIter + 1) * (curr_InnerIter+1));
  
  SetHistoryOutputValue("PHYS_TIME", UsedTime);
  
}

