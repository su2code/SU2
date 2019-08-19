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


#include "../../../Common/include/geometry_structure.hpp"
#include "../../include/solver_structure.hpp"

COutput::COutput(CConfig *config, unsigned short nDim) {
  
  this->nDim = nDim;

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
  
  Build_Offset_Cache = false;
  
  curr_InnerIter = 0;
  curr_OuterIter = 0;
  curr_TimeIter  = 0;
  
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
void COutput::SetMultizoneHistory_Output(COutput **output, CConfig **config, CConfig *driver_config, unsigned long TimeIter, unsigned long OuterIter){
  
  curr_TimeIter  = TimeIter;
  curr_OuterIter = OuterIter;
  
  bool write_header, write_screen, write_history;
  
  /*--- Retrieve residual and extra data -----------------------------------------------------------------*/
  
  LoadMultizoneHistoryData(output, config);
  
  Convergence_Monitoring(driver_config, curr_OuterIter);  
  
  /*--- Output using only the master node ---*/
  
  if (rank == MASTER_NODE && !no_writing) {
    
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
void COutput::SetCFL_Number(CSolver *****solver_container, CConfig **config, unsigned short val_iZone) {
  
  su2double CFLFactor = 1.0, power = 1.0, CFL = 0.0, CFLMin = 0.0, CFLMax = 0.0, Div = 1.0, Diff = 0.0, MGFactor[100];
  unsigned short iMesh;
  
  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  unsigned short nVar = 1;

  bool energy = config[val_iZone]->GetEnergy_Equation();
  bool weakly_coupled_heat = config[val_iZone]->GetWeakly_Coupled_Heat();

  switch( config[val_iZone]->GetKind_Solver()) {
    case EULER : case NAVIER_STOKES : case RANS:
    case INC_EULER : case INC_NAVIER_STOKES : case INC_RANS:
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

  if ((fabs(Diff) <= RhoRes_New*1E-8) && (curr_InnerIter != 0)) { Div = 0.1; power = config[val_iZone]->GetCFL_AdaptParam(1); }

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
  case INC_EULER : case INC_NAVIER_STOKES : case INC_RANS:      
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
  
  /*---- Construct a data sorter object to partition and distribute
   *  the local data into linear chunks across the processors ---*/
  
  if (fem_output){
  
    data_sorter = new CFEMDataSorter(config, geometry, GlobalField_Counter, Local_Data);

  }  else {
    
    data_sorter = new CFVMDataSorter(config, geometry, GlobalField_Counter, Local_Data);
    
  }  

  /*--- Now that we know the number of fields, create the local data array to temporarily store the volume output 
   * before writing it to file ---*/
   
  Local_Data.resize(data_sorter->GetnLocalPointSort(), std::vector<su2double>(GlobalField_Counter, 0.0));
  
  /*--- Collect that data defined in the subclasses from the different processors ---*/
  
  if (rank == MASTER_NODE)
    cout << endl << "Loading solution output data locally on each rank." << endl;
  
  CollectVolumeData(config, geometry, solver_container);
  
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
    
    filewriter = new CTecplotBinaryFileWriter(Variable_Names, nDim, curr_TimeIter, GetHistoryFieldValue("TIME_STEP"));
      
    break;
    
  case TECPLOT:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
    sorter->SortConnectivity(config, geometry, true);
    
    /*--- Write tecplot binary ---*/
    
    if (rank == MASTER_NODE) {
      cout << "Writing Tecplot ASCII file solution file." << endl;
    }
    
    filewriter = new CTecplotFileWriter(Variable_Names, nDim, curr_TimeIter, GetHistoryFieldValue("TIME_STEP"));

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

void COutput::SetSurface_Output(CGeometry *geometry, CConfig *config, unsigned short format, bool time_dep){
  
  CParallelDataSorter* surface_sort = NULL;
  CFileWriter*         file_writer = NULL;
  
  if (fem_output){
    surface_sort = new CSurfaceFEMDataSorter(config, geometry, GlobalField_Counter, dynamic_cast<CFEMDataSorter*>(data_sorter));
  } else {
    surface_sort = new CSurfaceFVMDataSorter(config, geometry,GlobalField_Counter, dynamic_cast<CFVMDataSorter*>(data_sorter));    
  }

  /*--- Set the file writer --- */
  
  SetFileWriter(config, geometry, surface_sort, file_writer, format);
    
  /*--- Sort the surface output data --- */
  
  surface_sort->SortOutputData(config, geometry);
  
  string FileName = SurfaceFilename;
  
  /*--- Remove extension --- */
  
  unsigned short lastindex = FileName.find_last_of(".");
  FileName = FileName.substr(0, lastindex);
  
  /*--- Add time iteration if requested --- */
  
  if (time_dep){
    FileName = config->GetFilename(FileName, "", curr_TimeIter);
  }
  
  if (surface_sort->GetnElem() > 0){
  
    /*--- Write data to file --- */
  
    file_writer->Write_Data(FileName, surface_sort);
  
  }
  
  delete file_writer;
  delete surface_sort;
  
}

void COutput::SetVolume_Output(CGeometry *geometry, CConfig *config, unsigned short format, bool time_dep){
  
  string FileName = VolumeFilename;
  
  if(format == SU2_RESTART_ASCII || format == SU2_RESTART_BINARY){
    FileName = RestartFilename;
  }
  
  /*--- Remove extension --- */
  
  unsigned short lastindex = FileName.find_last_of(".");
  FileName = FileName.substr(0, lastindex);
  
  /*--- Add time iteration if requested --- */
  
  if (time_dep){
    FileName = config->GetFilename(FileName, "", curr_TimeIter);
  }
  
  CFileWriter* file_writer = NULL;

  /*--- Set the file writer --- */
  
  SetFileWriter(config, geometry, data_sorter, file_writer, format);
  
  /*--- Write data to file --- */
  
  file_writer->Write_Data(FileName, data_sorter);
  
  if ((rank == MASTER_NODE) && config->GetWrt_Performance()) {
    cout << "Wrote " << file_writer->Get_Filesize()/(1.0e6) << " MB to disk in ";
    cout << file_writer->Get_UsedTime() << " s. (" << file_writer->Get_Bandwidth() << " MB/s)." << endl;
  }
  
  if(format == SU2_RESTART_ASCII || format == SU2_RESTART_BINARY){
    config->SetRestart_Bandwidth_Agg(config->GetRestart_Bandwidth_Agg() + file_writer->Get_Bandwidth());
  }

  delete file_writer;
  
}


bool COutput::Convergence_Monitoring(CConfig *config, unsigned long Iteration) {

  unsigned short iCounter;
    
  Convergence = false;
  
  if( HistoryOutput_Map.count(Conv_Field) > 0 ){
    
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
    
    if (HistoryOutput_Map[Conv_Field].FieldType == TYPE_RESIDUAL || HistoryOutput_Map[Conv_Field].FieldType == TYPE_AUTO_RESIDUAL) {
      
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
    long iIter = config->GetRestart_Iter();
    if (SU2_TYPE::Int(iIter) < 10) SPRINTF (buffer, "_0000%d", SU2_TYPE::Int(iIter));
    if ((SU2_TYPE::Int(iIter) >= 10) && (SU2_TYPE::Int(iIter) < 100)) SPRINTF (buffer, "_000%d", SU2_TYPE::Int(iIter));
    if ((SU2_TYPE::Int(iIter) >= 100) && (SU2_TYPE::Int(iIter) < 1000)) SPRINTF (buffer, "_00%d", SU2_TYPE::Int(iIter));
    if ((SU2_TYPE::Int(iIter) >= 1000) && (SU2_TYPE::Int(iIter) < 10000)) SPRINTF (buffer, "_0%d", SU2_TYPE::Int(iIter));
    if (SU2_TYPE::Int(iIter) >= 10000) SPRINTF (buffer, "_%d", SU2_TYPE::Int(iIter));
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

void COutput::PreprocessVolumeOutput(CConfig *config){

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
  
}

void COutput::CollectVolumeData(CConfig* config, CGeometry* geometry, CSolver** solver){
  
  unsigned short iMarker = 0;
  unsigned long iPoint = 0, jPoint = 0;
  unsigned long iVertex = 0;
  
  /*--- Reset the offset cache and index --- */
  Offset_Cache_Index = 0;
  Offset_Cache.clear();
  
  if (fem_output){
    
    /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
     geometrical information for the FEM DG solver. ---*/
  
    CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);
  
    unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
    
    CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();
    
    /*--- Access the solution by looping over the owned volume elements. ---*/
  
    for(unsigned long l=0; l<nVolElemOwned; ++l) {

      for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
        
        Build_Offset_Cache = !Offset_Cache.size() ? true : false;
        
        LoadVolumeDataFEM(config, geometry, solver, l, jPoint, j);
        
        jPoint++;
        
        CheckOffsetCache();        

      }
    }
    
  } else {
    
    for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      /*--- Load the volume data into the Local_Data() array. --- */
      
      Build_Offset_Cache = !Offset_Cache.size() ? true : false;

      LoadVolumeData(config, geometry, solver, iPoint);

    }
    
    /*--- Reset the offset cache and index --- */
    Offset_Cache_Index = 0;
    Offset_Cache.clear(); 
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
  
        if(geometry->node[iPoint]->GetDomain()){
          
          Build_Offset_Cache = !Offset_Cache.size() ? true : false;
   
          LoadSurfaceData(config, geometry, solver, iPoint, iMarker, iVertex);
          
        }
      }   
    } 
  }
}

void COutput::SetVolumeOutputValue(string name, unsigned long iPoint, su2double value){
  
  if (Build_Offset_Cache){ 
    
    /*--- Build up the offset cache to speed up subsequent 
     * calls of this routine since the order of calls is 
     * the same for every value of iPoint --- */
    
    if (VolumeOutput_Map.count(name) > 0){
      const short Offset = VolumeOutput_Map[name].Offset;
      Offset_Cache.push_back(Offset);        
      if (Offset != -1){
        Local_Data[iPoint][Offset] = value;
      }
    } else {
      SU2_MPI::Error(string("Cannot find output field with name ") + name, CURRENT_FUNCTION);    
    }
  } else {
    
    /*--- Use the offset cache for the access ---*/
    
    const short Offset = Offset_Cache[Offset_Cache_Index++];
    if (Offset != -1){
      Local_Data[iPoint][Offset] = value;
    }   
    if (Offset_Cache_Index == Offset_Cache.size()){
      Offset_Cache_Index = 0;
    }
  }
  
}

void COutput::CheckOffsetCache(){
  
  if (!Offset_Cache_Checked){
    vector<short> Offset_Cache_Copy = Offset_Cache;
    
    /*--- Remove the -1 offset --- */
    
    Offset_Cache_Copy.erase(std::remove(Offset_Cache_Copy.begin(), Offset_Cache_Copy.end(), -1), 
                            Offset_Cache_Copy.end());
    
    /*--- Check if all offsets are unique. If thats not the case, then SetVolumeOutputValue() was called
       * more than once for the same output field. --- */
    
    vector<short>::iterator it = std::unique( Offset_Cache_Copy.begin(), Offset_Cache_Copy.end() );
    if (it != Offset_Cache_Copy.end() ){
      SU2_MPI::Error("Offset cache contains duplicate entries.", CURRENT_FUNCTION);
    }
  }
  Offset_Cache_Checked = true;
  
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

    for (unsigned short iField = 0; iField < HistoryOutputPerSurface_List.size(); iField++){
      for (unsigned short iMarker = 0; iMarker < HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField]].size(); iMarker++){
        HistoryOutputField &Field = HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField]][iMarker];
        if (Field.FieldType == TYPE_COEFFICIENT){
          if (config->GetDirectDiff() != NO_DERIVATIVE){
            SetHistoryOutputValue("D_" + HistoryOutputPerSurface_List[iField][iMarker], SU2_TYPE::GetDerivative(Field.Value));
          }
        }
      }
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
      AddHistoryOutput("REL_" + HistoryOutput_List[iField], "rel" + currentField.FieldName, currentField.ScreenFormat, "REL_" + currentField.OutputGroup, "Relative residual.", TYPE_AUTO_RESIDUAL);
      Average[currentField.OutputGroup] = true;
    }
  }
  
  map<string, bool>::iterator it = Average.begin();
  for (it = Average.begin(); it != Average.end(); it++){
    AddHistoryOutput("AVG_" + it->first, "avg[" + AverageGroupName[it->first] + "]", FORMAT_FIXED, "AVG_" + it->first , "Average residual over all solution variables.", TYPE_AUTO_RESIDUAL);
  }  
  
  for (unsigned short iField = 0; iField < HistoryOutput_List.size(); iField++){
    HistoryOutputField &currentField = HistoryOutput_Map[HistoryOutput_List[iField]];
    if (currentField.FieldType == TYPE_COEFFICIENT){
      AddHistoryOutput("TAVG_"   + HistoryOutput_List[iField], "tavg["  + currentField.FieldName + "]", currentField.ScreenFormat, "TAVG_"   + currentField.OutputGroup, "Time averaged values.", TYPE_AUTO_COEFFICIENT);
    }
  }
  for (unsigned short iField = 0; iField < HistoryOutput_List.size(); iField++){
    HistoryOutputField &currentField = HistoryOutput_Map[HistoryOutput_List[iField]];
    if (currentField.FieldType == TYPE_COEFFICIENT){
      AddHistoryOutput("D_"      + HistoryOutput_List[iField], "d["     + currentField.FieldName + "]", currentField.ScreenFormat, "D_"      + currentField.OutputGroup, "Derivative value (DIRECT_DIFF=YES)", TYPE_AUTO_COEFFICIENT);  
    }
  }
  
  for (unsigned short iField = 0; iField < HistoryOutput_List.size(); iField++){
    HistoryOutputField &currentField = HistoryOutput_Map[HistoryOutput_List[iField]];
    if (currentField.FieldType == TYPE_COEFFICIENT){
      AddHistoryOutput("D_TAVG_" + HistoryOutput_List[iField], "dtavg[" + currentField.FieldName + "]", currentField.ScreenFormat, "D_TAVG_" + currentField.OutputGroup, "Derivative of the time averaged value (DIRECT_DIFF=YES)", TYPE_AUTO_COEFFICIENT);  
    }
  }

  for (unsigned short iField = 0; iField < HistoryOutputPerSurface_List.size(); iField++){
    for (unsigned short iMarker = 0; iMarker < HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField]].size(); iMarker++){
      HistoryOutputField &Field = HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField]][iMarker];
      if (Field.FieldType == TYPE_COEFFICIENT){
        AddHistoryOutput("D_"      + HistoryOutputPerSurface_List[iField][iMarker], "d["     + Field.FieldName + "]", Field.ScreenFormat, "D_"      + Field.OutputGroup, "Derivative values for per-surface output.", TYPE_AUTO_COEFFICIENT);
      }
    }
  }
  
  if (HistoryOutput_Map[Conv_Field].FieldType == TYPE_COEFFICIENT){
    AddHistoryOutput("CAUCHY", "C["  + HistoryOutput_Map[Conv_Field].FieldName + "]", FORMAT_SCIENTIFIC, "CAUCHY","Cauchy residual value of field set with CONV_FIELD." ,TYPE_AUTO_COEFFICIENT);
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

  /* --- Always print header in the first iteration --- */
  
  if ((curr_InnerIter == 0) && 
      (curr_OuterIter == 0) && 
      (curr_TimeIter == RestartIter)){
    return true;
  }
  
  if (!PrintOutput(curr_TimeIter, ScreenWrt_Freq_Time)&& 
      !(curr_TimeIter == config->GetnTime_Iter() - 1)){
    return false;
  }
   
  /*--- If there is no inner or outer iteration, don't print header ---*/
  if (ScreenWrt_Freq_Outer == 0 && ScreenWrt_Freq_Inner == 0){
    return false;
  }
  
  /*--- Print header if we are at the first inner iteration ---*/
  
  if (curr_InnerIter == 0){
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
  
  if (!PrintOutput(curr_TimeIter, ScreenWrt_Freq_Time)&& 
      !(curr_TimeIter == config->GetnTime_Iter() - 1)){
    
    return false;
    
  }
  
  if (Convergence) {return true;}
  
  if (!PrintOutput(curr_OuterIter, ScreenWrt_Freq_Outer) && 
      !(curr_OuterIter == config->GetnOuter_Iter() - 1)){
    
    return false;
    
  }
  
  if (!PrintOutput(curr_InnerIter, ScreenWrt_Freq_Inner) &&
      !(curr_InnerIter == config->GetnInner_Iter() - 1)){
    
    return false;
    
  }
 
  return true;
  
}

bool COutput::WriteHistoryFile_Output(CConfig *config) { 

  unsigned long HistoryWrt_Freq_Inner = config->GetHistory_Wrt_Freq(2);
  unsigned long HistoryWrt_Freq_Outer = config->GetHistory_Wrt_Freq(1);
  unsigned long HistoryWrt_Freq_Time  = config->GetHistory_Wrt_Freq(0);    
    
  /*--- Check if screen output should be written --- */
  
  if (!PrintOutput(curr_TimeIter, HistoryWrt_Freq_Time)&& 
      !(curr_TimeIter == config->GetnTime_Iter() - 1)){
    
    return false;
    
  }
  
  if (Convergence) {return true;}
  
  if (!PrintOutput(curr_OuterIter,HistoryWrt_Freq_Outer) && 
      !(curr_OuterIter == config->GetnOuter_Iter() - 1)){
    
    return false;
    
  }
  
  if (!PrintOutput(curr_InnerIter, HistoryWrt_Freq_Inner) &&
      !(curr_InnerIter == config->GetnInner_Iter() - 1)){
    
    return false;
    
  }
 
  return true;

}

void COutput::SetCommonHistoryFields(CConfig *config){
  
  /// BEGIN_GROUP: ITERATION, DESCRIPTION: Iteration identifier.
  /// DESCRIPTION: The time iteration index.
  AddHistoryOutput("TIME_ITER",     "Time_Iter",  FORMAT_INTEGER, "ITER", "Time iteration index"); 
  /// DESCRIPTION: The outer iteration index.
  AddHistoryOutput("OUTER_ITER",   "Outer_Iter",  FORMAT_INTEGER, "ITER", "Outer iteration index"); 
  /// DESCRIPTION: The inner iteration index.
  AddHistoryOutput("INNER_ITER",   "Inner_Iter", FORMAT_INTEGER,  "ITER", "Inner iteration index"); 
  /// END_GROUP
  
  /// BEGIN_GROUP: TIME_DOMAIN, DESCRIPTION: Time integration information
  /// Description: The current time
  AddHistoryOutput("CUR_TIME", "Cur_Time", FORMAT_SCIENTIFIC, "TIME_DOMAIN", "Current physical time (s)");
  /// Description: The current time step
  AddHistoryOutput("TIME_STEP", "Time_Step", FORMAT_SCIENTIFIC, "TIME_DOMAIN", "Current time step (s)");
 
  /// DESCRIPTION: Currently used wall-clock time.
  AddHistoryOutput("PHYS_TIME",   "Time(sec)", FORMAT_SCIENTIFIC, "PHYS_TIME", "Average wall-clock time"); 
  
}

void COutput::LoadCommonHistoryData(CConfig *config){
  
  SetHistoryOutputValue("TIME_ITER",  curr_TimeIter);  
  SetHistoryOutputValue("INNER_ITER", curr_InnerIter);
  SetHistoryOutputValue("OUTER_ITER", curr_OuterIter); 
  
  if (config->GetTime_Domain()){
    SetHistoryOutputValue("TIME_STEP", config->GetDelta_UnstTimeND()*config->GetTime_Ref());           
    if (curr_InnerIter == 0){
      SetHistoryOutputValue("CUR_TIME",  GetHistoryFieldValue("CUR_TIME") + GetHistoryFieldValue("TIME_STEP"));      
    }
  }
  
  su2double StopTime, UsedTime;
#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  
  UsedTime = (StopTime - config->Get_StartTime())/((curr_OuterIter + 1) * (curr_InnerIter+1));
  
  SetHistoryOutputValue("PHYS_TIME", UsedTime);
  
}


void COutput::PrintHistoryFields(){ 
  
  if (rank == MASTER_NODE){
    
    PrintingToolbox::CTablePrinter HistoryFieldTable(&std::cout);
    
    unsigned short NameSize = 0, GroupSize = 0, DescrSize = 0;
    
    for (unsigned short iField = 0; iField < HistoryOutput_List.size(); iField++){
      
      HistoryOutputField &Field = HistoryOutput_Map[HistoryOutput_List[iField]];
      
      if (Field.Description != ""){
        if (HistoryOutput_List[iField].size() > NameSize){
          NameSize = HistoryOutput_List[iField].size();
        }
        if (Field.OutputGroup.size() > GroupSize){
          GroupSize = Field.OutputGroup.size();
        }
        if (Field.Description.size() > DescrSize){
          DescrSize = Field.Description.size();
        }
      }
    }
    
    cout << "Available output fields for the current configuration in " << MultiZoneHeaderString << ":" << endl;
    
    HistoryFieldTable.AddColumn("Name", NameSize);
    HistoryFieldTable.AddColumn("Group Name", GroupSize);
    HistoryFieldTable.AddColumn("Type",5);
    HistoryFieldTable.AddColumn("Description", DescrSize);
    HistoryFieldTable.SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    
    HistoryFieldTable.PrintHeader();
    
    for (unsigned short iField = 0; iField < HistoryOutput_List.size(); iField++){
      
      HistoryOutputField &Field = HistoryOutput_Map[HistoryOutput_List[iField]];
      
      if (Field.FieldType == TYPE_DEFAULT || Field.FieldType == TYPE_COEFFICIENT || Field.FieldType == TYPE_RESIDUAL){
        string type;
        switch (Field.FieldType) {
          case TYPE_COEFFICIENT:
            type = "C";
            break;
          case TYPE_RESIDUAL:
            type = "R";
            break;
          default:
            type = "D";
            break;
        }
        
        if (Field.Description != "")
          HistoryFieldTable << HistoryOutput_List[iField] << Field.OutputGroup << type << Field.Description;
        
      }
    }
    
    HistoryFieldTable.PrintFooter();
    
    cout << "Type legend: Default (D), Residual (R), Coefficient (C)" << endl;
    
    cout << "Generated output fields (only first field of every group is shown):" << endl;
    
    PrintingToolbox::CTablePrinter ModifierTable(&std::cout);
    
    ModifierTable.AddColumn("Name", NameSize);
    ModifierTable.AddColumn("Group Name", GroupSize);
    ModifierTable.AddColumn("Type",5);
    ModifierTable.AddColumn("Description", DescrSize);
    ModifierTable.SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    ModifierTable.PrintHeader();
    
    std::map<string, bool> GroupVisited;
    
    for (unsigned short iField = 0; iField < HistoryOutput_List.size(); iField++){
      
      HistoryOutputField &Field = HistoryOutput_Map[HistoryOutput_List[iField]];
      
      if ((Field.FieldType == TYPE_AUTO_COEFFICIENT || Field.FieldType == TYPE_AUTO_RESIDUAL) && (GroupVisited.count(Field.OutputGroup) == 0)){
        string type;
        switch (Field.FieldType) {
          case TYPE_AUTO_COEFFICIENT:
            type = "AC";
            break;
          case TYPE_AUTO_RESIDUAL:
            type = "AR";
            break;
          default:
            type = "AD";
            break;
        }
        
        if (Field.Description != "")
          ModifierTable << HistoryOutput_List[iField] << Field.OutputGroup << type << Field.Description;
        
        GroupVisited[Field.OutputGroup] = true;
      }
    }   
    ModifierTable.PrintFooter();

  }
}

void COutput::PrintVolumeFields(){
  
  if (rank == MASTER_NODE){
    
    PrintingToolbox::CTablePrinter VolumeFieldTable(&std::cout);
    
    unsigned short NameSize = 0, GroupSize = 0, DescrSize = 0;
    
    for (unsigned short iField = 0; iField < VolumeOutput_List.size(); iField++){
      
      VolumeOutputField &Field = VolumeOutput_Map[VolumeOutput_List[iField]];
      
      if (Field.Description != ""){
        if (VolumeOutput_List[iField].size() > NameSize){
          NameSize = VolumeOutput_List[iField].size();
        }
        if (Field.OutputGroup.size() > GroupSize){
          GroupSize = Field.OutputGroup.size();
        }
        if (Field.Description.size() > DescrSize){
          DescrSize = Field.Description.size();
        }
      }
    }
    
    cout << "Available output fields for the current configuration in " << MultiZoneHeaderString << ":" << endl;
    
    VolumeFieldTable.AddColumn("Name", NameSize);
    VolumeFieldTable.AddColumn("Group Name", GroupSize);
    VolumeFieldTable.AddColumn("Description", DescrSize);
    VolumeFieldTable.SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    
    VolumeFieldTable.PrintHeader();
    
    for (unsigned short iField = 0; iField < VolumeOutput_List.size(); iField++){
      
      VolumeOutputField &Field = VolumeOutput_Map[VolumeOutput_List[iField]];

      if (Field.Description != "")
        VolumeFieldTable << VolumeOutput_List[iField] << Field.OutputGroup << Field.Description;
      
    }
    
    VolumeFieldTable.PrintFooter();
  }
}
