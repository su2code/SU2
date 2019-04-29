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

#include "../../include/output/output.hpp"

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
  
  /*--- Initialize point and connectivity counters to zero. ---*/
  
  nGlobal_Poin      = 0;
  nSurf_Poin        = 0;
  nGlobal_Elem      = 0;
  nSurf_Elem        = 0;
  nGlobal_Tria      = 0;
  nGlobal_Quad      = 0;
  nGlobal_Tetr      = 0;
  nGlobal_Hexa      = 0;
  nGlobal_Pris      = 0;
  nGlobal_Pyra      = 0;
  nGlobal_Line      = 0;
  nGlobal_BoundTria = 0;
  nGlobal_BoundQuad = 0;

  /*--- Initialize pointers to NULL ---*/
  
  Coords = NULL;
  Conn_Line = NULL;     Conn_BoundTria = NULL;  Conn_BoundQuad = NULL;
  Conn_Tria = NULL;     Conn_Quad = NULL;       Conn_Tetr = NULL;
  Conn_Hexa = NULL;     Conn_Pris = NULL;       Conn_Pyra = NULL;
  Data = NULL;
  
  /*--- Initialize parallel pointers to NULL ---*/
  
  nGlobal_Poin_Par    = 0;
  nGlobal_Elem_Par    = 0;
  nGlobal_Surf_Poin   = 0;
  nParallel_Poin      = 0;
  nSurf_Poin_Par      = 0;
  nSurf_Elem_Par      = 0;
  nParallel_Tria      = 0;
  nParallel_Quad      = 0;
  nParallel_Tetr      = 0;
  nParallel_Hexa      = 0;
  nParallel_Pris      = 0;
  nParallel_Pyra      = 0;
  nParallel_Line      = 0;
  nParallel_BoundTria = 0;
  nParallel_BoundQuad = 0;
  
  /*--- Initialize pointers to NULL ---*/
  
  Conn_BoundLine_Par = NULL;  Conn_BoundTria_Par = NULL;  Conn_BoundQuad_Par = NULL;
  Conn_Tria_Par = NULL;  Conn_Quad_Par = NULL;       Conn_Tetr_Par = NULL;
  Conn_Hexa_Par = NULL;  Conn_Pris_Par = NULL;       Conn_Pyra_Par = NULL;
  
  Local_Data_Copy    = NULL;
  Parallel_Data      = NULL;
  Parallel_Surf_Data = NULL;

  /*--- Initialize structures for storing linear partitioning offsets ---*/

  nGlobalPoint_Sort = 0;
  nLocalPoint_Sort  = 0;
  nPoint_Restart    = 0;

  Local_Halo_Sort = NULL;

  beg_node = NULL;
  end_node = NULL;

  nPoint_Lin = NULL;
  nPoint_Cum = NULL;
  
  /*--- Inlet profile data structures. ---*/

  nRow_InletFile    = NULL;
  nRowCum_InletFile = NULL;
  InletCoords       = NULL;

  Marker_Tags_InletFile.clear();

  /*--- Initialize CGNS write flag ---*/
  
  wrote_base_file = false;
  
  /*--- Initialize CGNS write flag ---*/
  
  wrote_CGNS_base = false;
  
  /*--- Initialize Tecplot surface flag ---*/
  
  wrote_surf_file = false;
  
  /*--- Initialize Paraview write flag ---*/
  
  wrote_Paraview_base = false;
  
  /*--- Initialize residual ---*/
  
  RhoRes_New = EPS;
  RhoRes_Old = new su2double[config->GetnZone()];
  for (iZone = 0; iZone < config->GetnZone(); iZone++) RhoRes_Old[iZone] = EPS;
  
  wrote_Paraview_base = false;

  
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
  
}

COutput::~COutput(void) {
  /* delete pointers initialized at construction*/
  /* Coords and Conn_*(Connectivity) have their own dealloc functions */
  /* Data is taken care of in DeallocateSolution function */

  if (RhoRes_Old != NULL) delete [] RhoRes_Old;

  
  /*--- Deallocate the structures holding the linear partitioning ---*/

  if (Local_Halo_Sort != NULL) delete [] Local_Halo_Sort;

  if (beg_node != NULL) delete [] beg_node;
  if (end_node != NULL) delete [] end_node;

  if (nPoint_Lin != NULL) delete [] nPoint_Lin;
  if (nPoint_Cum != NULL) delete [] nPoint_Cum;
  
  delete ConvergenceTable;
  delete MultiZoneHeaderTable;
  
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
  
  LoadHistoryData(config, geometry, solver_container);
  
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
  
  LoadHistoryData(config, geometry, solver_container);
  
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
  
  /*--- Sort the data, needed for volume and surface output ---*/
  
  if (rank == MASTER_NODE) 
    cout << "Sorting output data across all ranks." << endl;
  
  if (fem_output){
    SortOutputData_FEM(config, geometry);
  }
  else {
    SortOutputData(config, geometry);
  }
}

void COutput::SetSurface_Output(CGeometry *geometry, CConfig *config, unsigned short format){
  
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  /*--- Write files depending on the format --- */
  
  switch (format) {
  
  case CSV:
    
    SortConnectivity(config, geometry, true, true);
    
    if (rank == MASTER_NODE) {
        cout << "Writing surface CSV file." << endl;
    }
    
    WriteSurface_CSV(config, geometry);
  
    break;  
    
  case TECPLOT_BINARY:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
    SortConnectivity(config, geometry, true, false);
    
    /*--- Write tecplot binary ---*/
    
    if (rank == MASTER_NODE) {
      cout << "Writing Tecplot binary file surface solution file." << endl;  
    }
    
    WriteTecplotBinary_Parallel(config, geometry, iZone, config->GetnZone(), true);
      
    break;
    
  case TECPLOT:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
    SortConnectivity(config, geometry, true, true);
    
    /*--- Write tecplot binary ---*/
    
    if (rank == MASTER_NODE) {
      cout << "Writing Tecplot ASCII file surface solution file." << endl;
    }
    
    WriteTecplotASCII_Parallel(config, geometry, iZone, config->GetnZone(), true);

    break;
    
  case PARAVIEW_BINARY:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
    SortConnectivity(config, geometry, true, true);
    
    /*--- Write paraview binary ---*/
    
    if (rank == MASTER_NODE) {
      cout << "Writing Paraview binary file surface solution file." << endl;
    }
    
    WriteParaViewBinary_Parallel(config, geometry, iZone, nZone, true);
    
    break;
    
  case PARAVIEW:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
    SortConnectivity(config, geometry, true, true);
    
    /*--- Write paraview binary ---*/
    
    if (rank == MASTER_NODE) {
      cout << "Writing Paraview ASCII file surface solution file." << endl;
    }

    WriteParaViewASCII_Parallel(config, geometry, iZone, nZone, true);

    break;
  default:
    SU2_MPI::Error("Requested surface output format not available.", CURRENT_FUNCTION);
    break;
  } 

  /*--- Clean up the surface data that was only needed for output. ---*/
 
  DeallocateConnectivity_Parallel(true);
  
  DeallocateSurfaceData_Parallel();
  
}

void COutput::SetVolume_Output(CGeometry *geometry, CConfig *config, unsigned short format){
   
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  /*--- Write files depending on the format --- */
  
  switch (format) {
    
  case SU2_RESTART_ASCII:
    
    if (rank == MASTER_NODE) {
        cout << "Writing SU2 ASCII restart file." << endl;     
    }
   
    WriteRestart_Parallel_ASCII(config, geometry);
    
    break;
    
  case SU2_RESTART_BINARY:
    
    if (rank == MASTER_NODE) {
      cout << "Writing SU2 binary restart file." << endl;   
    }

    WriteRestart_Parallel_Binary(config, geometry);
    
    break;
  
  case SU2_MESH:
    
    /*--- Load and sort the output data and connectivity.
     *  Note that the solver container is not need in this case. ---*/
    
    SortConnectivity(config, geometry, false, true);
    
    /*--- Set the mesh ASCII format ---*/
    
    if (rank == MASTER_NODE) {
      cout << "Writing SU2 mesh file." << endl;
    }
    
    SetSU2_MeshASCII(config, geometry);
    
    break;    
    
  case TECPLOT_BINARY:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
    SortConnectivity(config, geometry, false, false);
    
    /*--- Write tecplot binary ---*/
    
    if (rank == MASTER_NODE) {
      cout << "Writing Tecplot binary file volume solution file." << endl;
    }
    
    WriteTecplotBinary_Parallel(config, geometry, iZone, config->GetnZone(), false);
      
    break;
    
  case TECPLOT:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
    SortConnectivity(config, geometry, false, true);
    
    /*--- Write tecplot binary ---*/
    
    if (rank == MASTER_NODE) {
      cout << "Writing Tecplot ASCII file volume solution file." << endl;
    }
    
    WriteTecplotASCII_Parallel(config, geometry, iZone, config->GetnZone(), false);

    break;
    
  case PARAVIEW_BINARY:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
    SortConnectivity(config, geometry, false, true);
    
    /*--- Write paraview binary ---*/
    if (rank == MASTER_NODE) {
      cout << "Writing Paraview binary file volume solution file." << endl;
    }
    
    WriteParaViewBinary_Parallel(config, geometry, iZone, nZone, false);
    
    break;
    
  case PARAVIEW:
    
    /*--- Load and sort the output data and connectivity. ---*/
    
    SortConnectivity(config, geometry, false, true);
    
    /*--- Write paraview binary ---*/
    if (rank == MASTER_NODE) {
        cout << "Writing Paraview ASCII file volume solution file." << endl;      
    }
    
    WriteParaViewASCII_Parallel(config, geometry, iZone, nZone, false);

    break;
    
  default:
    SU2_MPI::Error("Requested volume output format not available.", CURRENT_FUNCTION);
    break;
  } 

  /*--- Clean up the surface data that was only needed for output. ---*/
  
  if (format != SU2_RESTART_ASCII && format != SU2_RESTART_BINARY)
    DeallocateConnectivity_Parallel(false);
  
}



void COutput::SortConnectivity(CConfig *config, CGeometry *geometry, bool surf, bool val_sort) {

  /*--- Sort connectivity for each type of element (excluding halos). Note
   In these routines, we sort the connectivity into a linear partitioning
   across all processors based on the global index of the grid nodes. ---*/
  
  /*--- Sort volumetric grid connectivity. ---*/
  
  if (!surf) {
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
      cout <<"Sorting volumetric grid connectivity." << endl;
    
    if (!fem_output){
      SortVolumetricConnectivity(config, geometry, TRIANGLE,      val_sort);
      SortVolumetricConnectivity(config, geometry, QUADRILATERAL, val_sort);
      SortVolumetricConnectivity(config, geometry, TETRAHEDRON,   val_sort);
      SortVolumetricConnectivity(config, geometry, HEXAHEDRON,    val_sort);
      SortVolumetricConnectivity(config, geometry, PRISM,         val_sort);
      SortVolumetricConnectivity(config, geometry, PYRAMID,       val_sort);  
    } else {
      SortVolumetricConnectivity_FEM(config, geometry, TRIANGLE     );
      SortVolumetricConnectivity_FEM(config, geometry, QUADRILATERAL);
      SortVolumetricConnectivity_FEM(config, geometry, TETRAHEDRON  );
      SortVolumetricConnectivity_FEM(config, geometry, HEXAHEDRON   );
      SortVolumetricConnectivity_FEM(config, geometry, PRISM        );
      SortVolumetricConnectivity_FEM(config, geometry, PYRAMID      );
    }
    
    /*--- Reduce the total number of cells we will be writing in the output files. ---*/
    
    unsigned long nTotal_Elem = nParallel_Tria + nParallel_Quad + nParallel_Tetr + nParallel_Hexa + nParallel_Pris + nParallel_Pyra;
#ifndef HAVE_MPI
  nGlobal_Elem_Par = nTotal_Elem;
#else
  SU2_MPI::Allreduce(&nTotal_Elem, &nGlobal_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  }
  
  /*--- Sort surface grid connectivity. ---*/
  
  else {
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
      cout <<"Sorting surface grid connectivity." << endl;
    
    if (!fem_output){
      SortSurfaceConnectivity(config, geometry, LINE         );
      SortSurfaceConnectivity(config, geometry, TRIANGLE     );
      SortSurfaceConnectivity(config, geometry, QUADRILATERAL);
    } else {
      SortSurfaceConnectivity_FEM(config, geometry, LINE         );
      SortSurfaceConnectivity_FEM(config, geometry, TRIANGLE     );
      SortSurfaceConnectivity_FEM(config, geometry, QUADRILATERAL);   
    }    
    
    unsigned long nTotal_Surf_Elem = nParallel_Line + nParallel_BoundTria + nParallel_BoundQuad;
#ifndef HAVE_MPI
    nSurf_Elem_Par   = nTotal_Surf_Elem;
#else
    SU2_MPI::Allreduce(&nTotal_Surf_Elem, &nSurf_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    
    /*--- Sort the surface output data --- */
    if (fem_output){
      SortOutputData_Surface_FEM(config, geometry);
    } else {
      SortOutputData_Surface(config, geometry);
    }
  }
}


void COutput::SortVolumetricConnectivity(CConfig *config,
                                         CGeometry *geometry,
                                         unsigned short Elem_Type,
                                         bool val_sort) {
  
  unsigned long iProcessor;
  unsigned short NODES_PER_ELEMENT = 0;
  unsigned long iPoint, jPoint, kPoint, nLocalPoint, nTotalPoint;
  unsigned long nElem_Total = 0, Global_Index;
  
  unsigned long iVertex, iMarker;
  int SendRecv, RecvFrom;
  
  bool notPeriodic, notHalo, addedPeriodic, isPeriodic;
  
  int *Local_Halo = NULL;
  int *Conn_Elem  = NULL;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  switch (Elem_Type) {
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case PRISM:
      NODES_PER_ELEMENT = N_POINTS_PRISM;
      break;
    case PYRAMID:
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
  }
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that were part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  
  if (config->GetKind_SU2() != SU2_DEF) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
              (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
              (SendRecv < 0)) {
            Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
          }
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found. ---*/
  
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];
  
  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++)
    Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain.
         For SU2_DEF we want to keep all periodic nodes. ---*/
        
        if (config->GetKind_SU2() == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }
        
        notPeriodic = (isPeriodic && (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < (unsigned long)size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (Global_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          
          kPoint = (iProcessor+1)*maxAddedPeriodic;
          
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
        
      }
    }
  }
  
  /*--- Now that we've done the gymnastics to find any periodic points,
   compute the total number of local and global points for the output. ---*/
  
  nLocalPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif
  
  /*--- Compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  unsigned long *nPoint_Linear = new unsigned long[size+1];
  
  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }
  
  /*--- Get the number of remainder points after the even division. ---*/
  
  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }
  
  /*--- Store the local number of nodes and the beginning/end index ---*/
  
  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[size] = nTotalPoint;
  
  /*--- We start with the connectivity distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many elements we must send to each other rank in order to
   have all elements sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first nPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/
  
  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];
  
  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;
  
  for (int ii = 0; ii < (int)geometry->GetnElem(); ii++ ) {
    if (geometry->elem[ii]->GetVTK_Type() == Elem_Type) {
      for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {
        
        /*--- Get the index of the current point. ---*/
        
        iPoint = geometry->elem[ii]->GetNode(jj);
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- Search for the lowest global index in this element. We
         send the element to the processor owning the range that includes
         the lowest global index value. ---*/
        
        for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
          jPoint = geometry->elem[ii]->GetNode(kk);
          unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
          if (newID < Global_Index) Global_Index = newID;
        }
        
        /*--- Search for the processor that owns this point. If we are
         sorting the elements, we use the linear partitioning to find
         the rank, otherwise, we simply have the current rank load its
         own elements into the connectivity data structure. ---*/
        
        if (val_sort) {
          iProcessor = Global_Index/npoint_procs[0];
          if (iProcessor >= (unsigned long)size)
            iProcessor = (unsigned long)size-1;
          if (Global_Index >= nPoint_Linear[iProcessor])
            while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
          else
            while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
        } else {
          iProcessor = rank;
        }
        
        
        /*--- If we have not visited this element yet, increment our
         number of elements that must be sent to a particular proc. ---*/
        
        if ((nElem_Flag[iProcessor] != ii)) {
          nElem_Flag[iProcessor] = ii;
          nElem_Send[iProcessor+1]++;
        }
        
      }
    }
  }
  
  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif
  
  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
  
  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;
    
    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }
  
  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/
  
  unsigned long *connSend = NULL;
  connSend = new unsigned long[NODES_PER_ELEMENT*nElem_Send[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Send[size]; ii++)
    connSend[ii] = 0;
  
  /*--- Allocate arrays for storing halo flags. ---*/
  
  unsigned short *haloSend = new unsigned short[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    haloSend[ii] = false;
  
  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/
  
  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = NODES_PER_ELEMENT*nElem_Send[ii];
  
  unsigned long *haloIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) haloIndex[ii] = nElem_Send[ii];
  
  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/
  
  for (int ii = 0; ii < (int)geometry->GetnElem(); ii++) {
    if (geometry->elem[ii]->GetVTK_Type() == Elem_Type) {
      for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {
        
        /*--- Get the index of the current point. ---*/
        
        iPoint = geometry->elem[ii]->GetNode(jj);
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- Search for the lowest global index in this element. We
         send the element to the processor owning the range that includes
         the lowest global index value. ---*/
        
        for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
          jPoint = geometry->elem[ii]->GetNode(kk);
          unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
          if (newID < Global_Index) Global_Index = newID;
        }
       
        /*--- Search for the processor that owns this point. If we are
         sorting the elements, we use the linear partitioning to find
         the rank, otherwise, we simply have the current rank load its
         own elements into the connectivity data structure. ---*/
        
        if (val_sort) {
          iProcessor = Global_Index/npoint_procs[0];
          if (iProcessor >= (unsigned long)size)
            iProcessor = (unsigned long)size-1;
          if (Global_Index >= nPoint_Linear[iProcessor])
            while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
          else
            while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
        } else {
          iProcessor = rank;
        }
        
        
        /*--- Load connectivity into the buffer for sending ---*/
        
        if (nElem_Flag[iProcessor] != ii) {
          
          nElem_Flag[iProcessor] = ii;
          unsigned long nn = index[iProcessor];
          unsigned long mm = haloIndex[iProcessor];
          
          /*--- Load the connectivity values. ---*/
          
          for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
            iPoint = geometry->elem[ii]->GetNode(kk);
            connSend[nn] = geometry->node[iPoint]->GetGlobalIndex(); nn++;
            
            /*--- Check if this is a halo node. If so, flag this element
             as a halo cell. We will use this later to sort and remove
             any duplicates from the connectivity list. ---*/
            
            if (Local_Halo[iPoint]) haloSend[mm] = true;
            
          }
          
          /*--- Increment the index by the message length ---*/
          
          index[iProcessor]    += NODES_PER_ELEMENT;
          haloIndex[iProcessor]++;
          
        }
      }
    }
  }
  
  /*--- Free memory after loading up the send buffer. ---*/
  
  delete [] index;
  delete [] haloIndex;
  
  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  unsigned long *connRecv = NULL;
  connRecv = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Recv[size]; ii++)
    connRecv[ii] = 0;
  
  unsigned short *haloRecv = new unsigned short[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    haloRecv[ii] = false;
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/
  
  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];
  
  /*--- Launch the non-blocking recv's for the connectivity. ---*/
  
  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = NODES_PER_ELEMENT*nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = NODES_PER_ELEMENT*nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Repeat the process to communicate the halo flags. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(haloRecv[ll]), count, MPI_UNSIGNED_SHORT, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the halo flags. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(haloSend[ll]), count, MPI_UNSIGNED_SHORT, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif
  
  /*--- Copy my own rank's data into the recv buffer directly. ---*/
  
  int mm = NODES_PER_ELEMENT*nElem_Recv[rank];
  int ll = NODES_PER_ELEMENT*nElem_Send[rank];
  int kk = NODES_PER_ELEMENT*nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];
  
  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) haloRecv[mm] = haloSend[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. Note that we add 1 here
   to the connectivity for vizualization packages. First, allocate
   appropriate amount of memory for this section. ---*/
  
  if (nElem_Recv[size] > 0) Conn_Elem = new int[NODES_PER_ELEMENT*nElem_Recv[size]];
  int count = 0; nElem_Total = 0;
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    if (!haloRecv[ii]) {
      nElem_Total++;
      for (int jj = 0; jj < NODES_PER_ELEMENT; jj++) {
        Conn_Elem[count] = (int)connRecv[ii*NODES_PER_ELEMENT+jj] + 1;
        count++;
      }
    }
  }
  
  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
  switch (Elem_Type) {
    case TRIANGLE:
      nParallel_Tria = nElem_Total;
      if (nParallel_Tria > 0) Conn_Tria_Par = Conn_Elem;
      break;
    case QUADRILATERAL:
      nParallel_Quad = nElem_Total;
      if (nParallel_Quad > 0) Conn_Quad_Par = Conn_Elem;
      break;
    case TETRAHEDRON:
      nParallel_Tetr = nElem_Total;
      if (nParallel_Tetr > 0) Conn_Tetr_Par = Conn_Elem;
      break;
    case HEXAHEDRON:
      nParallel_Hexa = nElem_Total;
      if (nParallel_Hexa > 0) Conn_Hexa_Par = Conn_Elem;
      break;
    case PRISM:
      nParallel_Pris = nElem_Total;
      if (nParallel_Pris > 0) Conn_Pris_Par = Conn_Elem;
      break;
    case PYRAMID:
      nParallel_Pyra = nElem_Total;
      if (nParallel_Pyra > 0) Conn_Pyra_Par = Conn_Elem;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      break;
  }
  
  /*--- Free temporary memory from communications ---*/
  
  delete [] connSend;
  delete [] connRecv;
  delete [] haloSend;
  delete [] haloRecv;
  delete [] Local_Halo;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic; 
  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear;

}

void COutput::SortSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  unsigned long iProcessor;
  unsigned short NODES_PER_ELEMENT;
  unsigned long iPoint, jPoint, kPoint, nLocalPoint, nTotalPoint;
  unsigned long nElem_Total = 0, Global_Index;
  
  unsigned long iVertex, iMarker;
  int SendRecv, RecvFrom;
  
  bool notPeriodic, notHalo, addedPeriodic, isPeriodic;
  
  int *Local_Halo = NULL;
  int *Conn_Elem  = NULL;
  
#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  switch (Elem_Type) {
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      NODES_PER_ELEMENT = 0;
      break;
  }
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that were part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  
  if (config->GetKind_SU2() != SU2_DEF) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
              (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
              (SendRecv < 0)) {
            Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
          }
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found. ---*/
  
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];
  
  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++)
    Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain.
         For SU2_DEF we want to keep all periodic nodes. ---*/
        
        if (config->GetKind_SU2() == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }
        
        notPeriodic = (isPeriodic && (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < (unsigned long)size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (Global_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          
          kPoint = (iProcessor+1)*maxAddedPeriodic;
          
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
        
      }
    }
  }
  
  /*--- Now that we've done the gymnastics to find any periodic points,
   compute the total number of local and global points for the output. ---*/
  
  nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif
  
  /*--- Compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  unsigned long *nPoint_Linear = new unsigned long[size+1];
  
  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }
  
  /*--- Get the number of remainder points after the even division. ---*/
  
  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }
  
  /*--- Store the local number of nodes and the beginning/end index ---*/
  
  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[size] = nTotalPoint;
  
  /*--- We start with the connectivity distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many elements we must send to each other rank in order to
   have all elements sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first nPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/
  
  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];
  
  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      
      for (int ii = 0; ii < (int)geometry->GetnElem_Bound(iMarker); ii++) {
        
        if (geometry->bound[iMarker][ii]->GetVTK_Type() == Elem_Type) {
          for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {
            
            /*--- Get the index of the current point. ---*/
            
            iPoint = geometry->bound[iMarker][ii]->GetNode(jj);
            Global_Index = geometry->node[iPoint]->GetGlobalIndex();
            
            /*--- Search for the lowest global index in this element. We
             send the element to the processor owning the range that includes
             the lowest global index value. ---*/
            
            for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
              jPoint = geometry->bound[iMarker][ii]->GetNode(kk);
              unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
              if (newID < Global_Index) Global_Index = newID;
            }
            
            /*--- Search for the processor that owns this point ---*/
            
            iProcessor = Global_Index/npoint_procs[0];
            if (iProcessor >= (unsigned long)size)
              iProcessor = (unsigned long)size-1;
            if (Global_Index >= nPoint_Linear[iProcessor])
              while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
            else
              while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
            
            /*--- If we have not visited this element yet, increment our
             number of elements that must be sent to a particular proc. ---*/
            
            if ((nElem_Flag[iProcessor] != ii)) {
              nElem_Flag[iProcessor] = ii;
              nElem_Send[iProcessor+1]++;
            }
            
          }
        }
      }
    }
  }
  
  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif
  
  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
  
  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;
    
    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }
  
  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/
  
  unsigned long *connSend = NULL;
  connSend = new unsigned long[NODES_PER_ELEMENT*nElem_Send[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Send[size]; ii++)
    connSend[ii] = 0;
  
  /*--- Allocate arrays for storing halo flags. ---*/
  
  unsigned short *haloSend = new unsigned short[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    haloSend[ii] = false;
  
  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/
  
  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = NODES_PER_ELEMENT*nElem_Send[ii];
  
  unsigned long *haloIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) haloIndex[ii] = nElem_Send[ii];
  
  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      
      for (int ii = 0; ii < (int)geometry->GetnElem_Bound(iMarker); ii++) {
        
        if (geometry->bound[iMarker][ii]->GetVTK_Type() == Elem_Type) {
          for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {
            
            /*--- Get the index of the current point. ---*/
            
            iPoint = geometry->bound[iMarker][ii]->GetNode(jj);
            Global_Index = geometry->node[iPoint]->GetGlobalIndex();
            
            /*--- Search for the lowest global index in this element. We
             send the element to the processor owning the range that includes
             the lowest global index value. ---*/
            
            for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
              jPoint = geometry->bound[iMarker][ii]->GetNode(kk);
              unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
              if (newID < Global_Index) Global_Index = newID;
            }
            
            /*--- Search for the processor that owns this point ---*/
            
            iProcessor = Global_Index/npoint_procs[0];
            if (iProcessor >= (unsigned long)size)
              iProcessor = (unsigned long)size-1;
            if (Global_Index >= nPoint_Linear[iProcessor])
              while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
            else
              while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
            
            /*--- Load connectivity into the buffer for sending ---*/
            
            if (nElem_Flag[iProcessor] != ii) {
              
              nElem_Flag[iProcessor] = ii;
              unsigned long nn = index[iProcessor];
              unsigned long mm = haloIndex[iProcessor];
              
              /*--- Load the connectivity values. ---*/
              
              for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
                iPoint = geometry->bound[iMarker][ii]->GetNode(kk);
                connSend[nn] = geometry->node[iPoint]->GetGlobalIndex(); nn++;
                
                /*--- Check if this is a halo node. If so, flag this element
                 as a halo cell. We will use this later to sort and remove
                 any duplicates from the connectivity list. ---*/
                
                if (Local_Halo[iPoint]) haloSend[mm] = true;
                
              }
              
              /*--- Increment the index by the message length ---*/
              
              index[iProcessor]    += NODES_PER_ELEMENT;
              haloIndex[iProcessor]++;
              
            }
          }
        }
      }
    }
  }
  
  /*--- Free memory after loading up the send buffer. ---*/
  
  delete [] index;
  delete [] haloIndex;
  
  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  unsigned long *connRecv = NULL;
  connRecv = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Recv[size]; ii++)
    connRecv[ii] = 0;
  
  unsigned short *haloRecv = new unsigned short[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    haloRecv[ii] = false;
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/
  
  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];
  
  /*--- Launch the non-blocking recv's for the connectivity. ---*/
  
  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = NODES_PER_ELEMENT*nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = NODES_PER_ELEMENT*nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Repeat the process to communicate the halo flags. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(haloRecv[ll]), count, MPI_UNSIGNED_SHORT, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the halo flags. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(haloSend[ll]), count, MPI_UNSIGNED_SHORT, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif
  
  /*--- Copy my own rank's data into the recv buffer directly. ---*/
  
  int mm = NODES_PER_ELEMENT*nElem_Recv[rank];
  int ll = NODES_PER_ELEMENT*nElem_Send[rank];
  int kk = NODES_PER_ELEMENT*nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];
  
  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) haloRecv[mm] = haloSend[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. Note that we add 1 here
   to the connectivity for vizualization packages. First, allocate
   appropriate amount of memory for this section. ---*/
  
  if (nElem_Recv[size] > 0) Conn_Elem = new int[NODES_PER_ELEMENT*nElem_Recv[size]];
  int count = 0; nElem_Total = 0;
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    if (!haloRecv[ii]) {
      nElem_Total++;
      for (int jj = 0; jj < NODES_PER_ELEMENT; jj++) {
        Conn_Elem[count] = (int)connRecv[ii*NODES_PER_ELEMENT+jj] + 1;
        count++;
      }
    }
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
  switch (Elem_Type) {
    case LINE:
      nParallel_Line = nElem_Total;
      if (nParallel_Line > 0) Conn_BoundLine_Par = Conn_Elem;
      break;
    case TRIANGLE:
      nParallel_BoundTria = nElem_Total;
      if (nParallel_BoundTria > 0) Conn_BoundTria_Par = Conn_Elem;
      break;
    case QUADRILATERAL:
      nParallel_BoundQuad = nElem_Total;
      if (nParallel_BoundQuad > 0) Conn_BoundQuad_Par = Conn_Elem;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      break;
  }
  
  /*--- Free temporary memory from communications ---*/
  
  delete [] connSend;
  delete [] connRecv;
  delete [] haloSend;
  delete [] haloRecv;
  delete [] Local_Halo;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear;
  
}

void COutput::SortOutputData(CConfig *config, CGeometry *geometry) {
  
  unsigned long iProcessor;
  unsigned long iPoint, Global_Index, nTotalPoint;
  
  int VARS_PER_POINT = GlobalField_Counter;
  
#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint_Sort, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint_Sort;
#endif
  
  /*--- Now that we know the actual number of points we need to output,
   compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  unsigned long *nPoint_Linear = new unsigned long[size+1];
  
  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }
  
  /*--- Get the number of remainder points after the even division. ---*/
  
  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }
  
  /*--- Store the local number of nodes and the beginning/end index ---*/
  
  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[size] = nTotalPoint;
  
  /*--- We start with the grid nodes distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many nodes we must send to each other rank in order to
   have all nodes sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first ~ nGlobalPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/
  
  int *nPoint_Send = new int[size+1]; nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1]; nPoint_Recv[0] = 0;
  int *nPoint_Flag = new int[size];
  
  for (int ii=0; ii < size; ii++) {
    nPoint_Send[ii] = 0;
    nPoint_Recv[ii] = 0;
    nPoint_Flag[ii]= -1;
  }
  nPoint_Send[size] = 0; nPoint_Recv[size] = 0;
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++ ) {
    
    /*--- We only write interior points and recovered periodic points. ---*/
    
    if (!Local_Halo_Sort[iPoint]) {
      
      /*--- Get the global index of the current point. ---*/
      
      Global_Index = geometry->node[iPoint]->GetGlobalIndex();
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear[iProcessor])
        while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
      
      /*--- If we have not visited this node yet, increment our
       number of elements that must be sent to a particular proc. ---*/
      
      if (nPoint_Flag[iProcessor] != (int)iPoint) {
        nPoint_Flag[iProcessor] = (int)iPoint;
        nPoint_Send[iProcessor+1]++;
      }
      
    }
  }
  
  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nPoint_Recv[1] = nPoint_Send[1];
#endif
  
  /*--- Prepare to send coordinates. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
  
  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nPoint_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nPoint_Recv[ii+1] > 0)) nRecvs++;
    
    nPoint_Send[ii+1] += nPoint_Send[ii];
    nPoint_Recv[ii+1] += nPoint_Recv[ii];
  }
  
  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/
  
  su2double *connSend = NULL;
  connSend = new su2double[VARS_PER_POINT*nPoint_Send[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Send[size]; ii++)
    connSend[ii] = 0;
  
  /*--- Allocate arrays for sending the global ID. ---*/
  
  unsigned long *idSend = new unsigned long[nPoint_Send[size]];
  for (int ii = 0; ii < nPoint_Send[size]; ii++)
    idSend[ii] = 0;
  
  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/
  
  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = VARS_PER_POINT*nPoint_Send[ii];
  
  unsigned long *idIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) idIndex[ii] = nPoint_Send[ii];
  
  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- We only write interior points and recovered periodic points. ---*/
    
    if (!Local_Halo_Sort[iPoint]) {
      
      /*--- Get the index of the current point. ---*/
      
      Global_Index = geometry->node[iPoint]->GetGlobalIndex();
      
      /*--- Search for the processor that owns this point. ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear[iProcessor])
        while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
      
      /*--- Load node coordinates into the buffer for sending. ---*/
      
      if (nPoint_Flag[iProcessor] != (int)iPoint) {
        
        nPoint_Flag[iProcessor] = (int)iPoint;
        unsigned long nn = index[iProcessor];
        
        /*--- Load the data values. ---*/
        
        for (unsigned short kk = 0; kk < VARS_PER_POINT; kk++) {
          connSend[nn] = Local_Data[iPoint][kk]; nn++;
        }
        
        /*--- Load the global ID (minus offset) for sorting the
         points once they all reach the correct processor. ---*/
        
        nn = idIndex[iProcessor];
        idSend[nn] = Global_Index - starting_node[iProcessor];
        
        /*--- Increment the index by the message length ---*/
        
        index[iProcessor]  += VARS_PER_POINT;
        idIndex[iProcessor]++;
        
      }
    }
  }
  
  /*--- Free memory after loading up the send buffer. ---*/
  
  delete [] index;
  delete [] idIndex;
  
  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  su2double *connRecv = NULL;
  connRecv = new su2double[VARS_PER_POINT*nPoint_Recv[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Recv[size]; ii++)
    connRecv[ii] = 0;
  
  unsigned long *idRecv = new unsigned long[nPoint_Recv[size]];
  for (int ii = 0; ii < nPoint_Recv[size]; ii++)
    idRecv[ii] = 0;
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the global IDs. ---*/
  
  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];
  
  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = VARS_PER_POINT*nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = VARS_PER_POINT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_DOUBLE, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = VARS_PER_POINT*nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = VARS_PER_POINT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_DOUBLE, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Repeat the process to communicate the global IDs. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the global IDs. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif
  
  /*--- Copy my own rank's data into the recv buffer directly. ---*/
  
  int mm = VARS_PER_POINT*nPoint_Recv[rank];
  int ll = VARS_PER_POINT*nPoint_Send[rank];
  int kk = VARS_PER_POINT*nPoint_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];
  
  mm = nPoint_Recv[rank];
  ll = nPoint_Send[rank];
  kk = nPoint_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. First, allocate the
   appropriate amount of memory for this section. ---*/
  
  Parallel_Data = new su2double*[VARS_PER_POINT];
  for (int jj = 0; jj < VARS_PER_POINT; jj++) {
    Parallel_Data[jj] = new su2double[nPoint_Recv[size]];
    for (int ii = 0; ii < nPoint_Recv[size]; ii++) {
      Parallel_Data[jj][idRecv[ii]] = connRecv[ii*VARS_PER_POINT+jj];
    }
  }
  
  /*--- Store the total number of local points my rank has for
   the current section after completing the communications. ---*/
  
  nParallel_Poin = nPoint_Recv[size];
  
  /*--- Reduce the total number of points we will write in the output files. ---*/

#ifndef HAVE_MPI
  nGlobal_Poin_Par = nParallel_Poin;
#else
  SU2_MPI::Allreduce(&nParallel_Poin, &nGlobal_Poin_Par, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Free temporary memory from communications ---*/
  
  delete [] connSend;
  delete [] connRecv;
  delete [] idSend;
  delete [] idRecv;
  delete [] nPoint_Recv;
  delete [] nPoint_Send;
  delete [] nPoint_Flag;

  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear;
  
}

void COutput::SortOutputData_Surface(CConfig *config, CGeometry *geometry) {
  
  unsigned short iMarker;
  unsigned long iProcessor;
  unsigned long iPoint, jPoint, kPoint, iElem;
  unsigned long Global_Index, nLocalPoint, nTotalPoint, iVertex;
  
  int VARS_PER_POINT = GlobalField_Counter;
  int *Local_Halo = NULL;
  int iNode, count;
  int SendRecv, RecvFrom;
  
  bool notPeriodic, notHalo, addedPeriodic, isPeriodic;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 1: We already have the surface connectivity spread out in     ---*/
  /*---         linear partitions across all procs and the output data     ---*/
  /*---         for the entire field is similarly linearly partitioned.    ---*/
  /*---         We need to identify which nodes in the volume data are     ---*/
  /*---         also surface points. Our first step is to loop over all    ---*/
  /*---         of the sorted surface connectivity and create a data       ---*/
  /*---         structure on each proc that can identify the local surf    ---*/
  /*---         points. Note that the linear partitioning is slightly      ---*/
  /*---         different between the nodes and elements, so we will       ---*/
  /*---         have to move between the two systems in this routine.      ---*/
  /*--------------------------------------------------------------------------*/
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. This is the linear partitioning for nodes. ---*/
  
  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      
      /*--- Checking for less than or equal to the rank, because there may
       be some periodic halo nodes that send info to the same rank. ---*/
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                      (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        if (isPeriodic) Local_Halo[iPoint] = false;
      }
    }
  }
  
  /*--- Sum total number of nodes that belong to the domain ---*/
  
  nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif
  
  /*--- Now that we know the actual number of points we need to output,
   compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  
  unsigned long *nPoint_Linear_Nodes = new unsigned long[size+1];
  unsigned long *nPoint_Linear_Elems = new unsigned long[size+1];
  
  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }
  
  /*--- Get the number of remainder points after the even division. ---*/
  
  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }
  
  /*--- Store the local number of nodes and the beginning/end index ---*/
  
  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear_Nodes[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear_Nodes[ii] = nPoint_Linear_Nodes[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear_Nodes[size] = nTotalPoint;
  
  /*--- Prepare to check and communicate the nodes that each proc has
   locally from the surface connectivity. ---*/
  
  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];
  
  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;
  
  /*--- Loop through our local line elements and check where each
   of the grid nodes resides based on global index. ---*/
  
  for (int ii = 0; ii < (int)nParallel_Line; ii++) {
    for ( int jj = 0; jj < N_POINTS_LINE; jj++ ) {
      
      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/
      
      iNode = ii*N_POINTS_LINE+jj;
      Global_Index = Conn_BoundLine_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
      
      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/
      
      if ((nElem_Flag[iProcessor] != iNode)) {
        nElem_Flag[iProcessor] = iNode;
        nElem_Send[iProcessor+1]++;
      }
      
    }
  }
  
  /*--- Reset out flags and then loop through our local triangle surface
   elements performing the same check for where each grid node resides. ---*/
  
  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;
  
  for (int ii = 0; ii < (int)nParallel_BoundTria; ii++) {
    for ( int jj = 0; jj < N_POINTS_TRIANGLE; jj++ ) {
      
      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/
      
      iNode = ii*N_POINTS_TRIANGLE + jj;
      Global_Index = Conn_BoundTria_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
      
      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/
      
      if ((nElem_Flag[iProcessor] != iNode)) {
        nElem_Flag[iProcessor] = iNode;
        nElem_Send[iProcessor+1]++;
      }
      
    }
  }
  
  /*--- Reset out flags and then loop through our local quad surface
   elements performing the same check for where each grid node resides. ---*/
  
  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;
  
  for (int ii = 0; ii < (int)nParallel_BoundQuad; ii++) {
    for ( int jj = 0; jj < N_POINTS_QUADRILATERAL; jj++ ) {
      
      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/
      
      iNode = ii*N_POINTS_QUADRILATERAL+jj;
      Global_Index = Conn_BoundQuad_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
      
      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/
      
      if ((nElem_Flag[iProcessor] != iNode)) {
        nElem_Flag[iProcessor] = iNode;
        nElem_Send[iProcessor+1]++;
      }
      
    }
  }
  
  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many nodes it will receive from each other processor. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif
  
  /*--- Prepare to send. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
  
  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;
    
    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }
  
  /*--- Allocate arrays for sending the global ID. ---*/
  
  unsigned long *idSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++) idSend[ii] = 0;
  
  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/
  
  unsigned long *idIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) idIndex[ii] = nElem_Send[ii];
  
  /*--- Now loop back through the local connectivities for the surface
   elements and load up the global IDs for sending to their home proc. ---*/
  
  for (int ii = 0; ii < (int)nParallel_Line; ii++) {
    for ( int jj = 0; jj < N_POINTS_LINE; jj++ ) {
      
      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/
      
      iNode = ii*N_POINTS_LINE+jj;
      Global_Index = Conn_BoundLine_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
      
      /*--- Load global ID into the buffer for sending ---*/
      
      if (nElem_Flag[iProcessor] != iNode) {
        
        nElem_Flag[iProcessor] = iNode;
        unsigned long nn = idIndex[iProcessor];
        
        /*--- Load the connectivity values. ---*/
        
        idSend[nn] = Global_Index; nn++;
        
        /*--- Increment the index by the message length ---*/
        
        idIndex[iProcessor]++;
        
      }
      
    }
  }
  
  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;
  
  for (int ii = 0; ii < (int)nParallel_BoundTria; ii++) {
    for ( int jj = 0; jj < N_POINTS_TRIANGLE; jj++ ) {
      
      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/
      
      iNode = ii*N_POINTS_TRIANGLE + jj;
      Global_Index = Conn_BoundTria_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
      
      /*--- Load global ID into the buffer for sending ---*/
      
      if (nElem_Flag[iProcessor] != iNode) {
        
        nElem_Flag[iProcessor] = iNode;
        unsigned long nn = idIndex[iProcessor];
        
        /*--- Load the connectivity values. ---*/
        
        idSend[nn] = Global_Index; nn++;
        
        /*--- Increment the index by the message length ---*/
        
        idIndex[iProcessor]++;
        
      }
      
    }
  }
  
  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;
  
  for (int ii = 0; ii < (int)nParallel_BoundQuad; ii++) {
    for ( int jj = 0; jj < N_POINTS_QUADRILATERAL; jj++ ) {
      
      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/
      
      iNode = ii*N_POINTS_QUADRILATERAL+jj;
      Global_Index = Conn_BoundQuad_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
      
      /*--- Load global ID into the buffer for sending ---*/
      
      if (nElem_Flag[iProcessor] != iNode) {
        
        nElem_Flag[iProcessor] = iNode;
        unsigned long nn = idIndex[iProcessor];
        
        /*--- Load the connectivity values. ---*/
        
        idSend[nn] = Global_Index; nn++;
        
        /*--- Increment the index by the message length ---*/
        
        idIndex[iProcessor]++;
        
      }
      
    }
  }
  
  /*--- Allocate the memory that we need for receiving the global IDs
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  unsigned long *idRecv = NULL;
  idRecv = new unsigned long[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    idRecv[ii] = 0;
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/
  
  send_req = new SU2_MPI::Request[nSends];
  recv_req = new SU2_MPI::Request[nRecvs];
  
  /*--- Launch the non-blocking recv's for the global IDs. ---*/
  
  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the global IDs. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
#endif
  
  /*--- Copy my own rank's data into the recv buffer directly. ---*/
  
  int mm = nElem_Recv[rank];
  int ll = nElem_Send[rank];
  int kk = nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  int number = nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Each proc now knows which is its local grid nodes from     ---*/
  /*---         the entire volume solution are part of the surface. We     ---*/
  /*---         now apply a mask to extract just those points on the       ---*/
  /*---         surface. We also need to perform a renumbering so that     ---*/
  /*---         the surface data (nodes and connectivity) have their       ---*/
  /*---         own global numbering. This is important for writing        ---*/
  /*---         output files in a later routine.                           ---*/
  /*--------------------------------------------------------------------------*/
  
  /*--- Create a local data structure that acts as a mask to extract the
   set of points within the local set that are on the surface. ---*/
  
  int *surfPoint = new int[nParallel_Poin];
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) surfPoint[iPoint] = -1;
  
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    surfPoint[(int)idRecv[ii]- starting_node[rank]] = (int)idRecv[ii];
  }
  
  /*--- First, add up the number of surface points I have on my rank. ---*/
  
  nSurf_Poin_Par = 0;
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
    if (surfPoint[iPoint] != -1) {
      nSurf_Poin_Par++;
    }
  }
  
  /*--- Communicate this number of local surface points to all other
   processors so that it can be used to create offsets for the new
   global numbering for the surface points. ---*/
  
  int *nPoint_Send = new int[size+1]; nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1]; nPoint_Recv[0] = 0;
  
  for (int ii=1; ii < size+1; ii++) nPoint_Send[ii]= (int)nSurf_Poin_Par;
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nPoint_Recv[1] = nPoint_Send[1];
#endif
  
  /*--- Go to cumulative storage format to compute the offsets. ---*/
  
  for (int ii = 0; ii < size; ii++) {
    nPoint_Send[ii+1] += nPoint_Send[ii];
    nPoint_Recv[ii+1] += nPoint_Recv[ii];
  }
  
  /*--- Now that we know the number of local surface points that we have,
   we can allocate the new data structure to hold these points alone. Here,
   we also copy the data for those points from our volume data structure. ---*/
  
  Parallel_Surf_Data = new su2double*[VARS_PER_POINT];
  for (int jj = 0; jj < VARS_PER_POINT; jj++) {
    Parallel_Surf_Data[jj] = new su2double[nSurf_Poin_Par];
    count = 0;
    for (int ii = 0; ii < (int)nParallel_Poin; ii++) {
      if (surfPoint[ii] !=-1) {
        Parallel_Surf_Data[jj][count] = Parallel_Data[jj][ii];
        count++;
      }
    }
  }
  
  /*--- Reduce the total number of surf points we have. This will be
   needed for writing the surface solution files later. ---*/
  
#ifndef HAVE_MPI
  nGlobal_Surf_Poin = nSurf_Poin_Par;
#else
  SU2_MPI::Allreduce(&nSurf_Poin_Par, &nGlobal_Surf_Poin, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Now that we know every proc's global offset for the number of
   surface points, we can create the new global numbering. Here, we
   create a new mapping using two arrays, which will need to be
   communicated. We use our mask again here.  ---*/
  
  unsigned long *globalP = new unsigned long[nSurf_Poin_Par];
  unsigned long *renumbP = new unsigned long[nSurf_Poin_Par];
  
  count = 0;
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
    if (surfPoint[iPoint] != -1) {
      globalP[count] = surfPoint[iPoint];
      renumbP[count] = count + nPoint_Recv[rank];
      count++;
    }
  }
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Communicate the arrays with the new global surface point   ---*/
  /*---         numbering to the procs that hold the connectivity for      ---*/
  /*---         each element. This will be done in two phases. First,      ---*/
  /*---         we send the arrays around to the other procs based on      ---*/
  /*---         the linear partitioning for the elems. This gets us        ---*/
  /*---         most of the way there, however, due to the type of         ---*/
  /*---         linear partitioning for the elements, there may exist      ---*/
  /*---         elements that have nodes outside of the linear part.       ---*/
  /*---         bounds. This is because the elems are distributed based    ---*/
  /*---         on the node with the smallest global ID.                   ---*/
  /*--------------------------------------------------------------------------*/
  
  /*--- First, we perform the linear partitioning again as it is done
   for elements, which is slightly different than for nodes (above). ---*/
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that were part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  
  if (config->GetKind_SU2() != SU2_DEF) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
              (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
              (SendRecv < 0)) {
            Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
          }
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found. ---*/
  
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];
  
  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++)
    Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain.
         For SU2_DEF we want to keep all periodic nodes. ---*/
        
        if (config->GetKind_SU2() == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }
        
        notPeriodic = (isPeriodic && (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < (unsigned long)size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (Global_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          
          kPoint = (iProcessor+1)*maxAddedPeriodic;
          
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
        
      }
    }
  }
  
  /*--- Now that we've done the gymnastics to find any periodic points,
   compute the total number of local and global points for the output. ---*/
  
  nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif
  
  /*--- Compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }
  
  /*--- Get the number of remainder points after the even division. ---*/
  
  rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }
  
  /*--- Store the local number of nodes and the beginning/end index ---*/
  
  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear_Elems[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear_Elems[ii] = nPoint_Linear_Elems[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear_Elems[size] = nTotalPoint;
  
  /*--- Reset our flags and counters ---*/
  
  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;
  
  /*--- Loop through my local surface nodes, find which proc the global
   value lives on, then communicate the global ID and remumbered value. ---*/
  
  for (int ii = 0; ii < (int)nSurf_Poin_Par; ii++) {
    
    Global_Index = globalP[ii];
    
    /*--- Search for the processor that owns this point ---*/
    
    iProcessor = Global_Index/npoint_procs[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Linear_Elems[iProcessor])
      while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;
    
    /*--- If we have not visited this element yet, increment our
     number of elements that must be sent to a particular proc. ---*/
    
    if ((nElem_Flag[iProcessor] != ii)) {
      nElem_Flag[iProcessor] = ii;
      nElem_Send[iProcessor+1]++;
    }
    
  }
  
  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif
  
  /*--- Prepare to send. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
  
  nSends = 0; nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;
    
    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }
  
  /*--- Allocate memory to hold the globals that we are
   sending. ---*/
  
  unsigned long *globalSend = NULL;
  globalSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    globalSend[ii] = 0;
  
  /*--- Allocate memory to hold the renumbering that we are
   sending. ---*/
  
  unsigned long *renumbSend = NULL;
  renumbSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    renumbSend[ii] = 0;
  
  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/
  
  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = nElem_Send[ii];
  
  /*--- Loop back through and load up the buffers for the global IDs
   and their new renumbering values. ---*/
  
  for (int ii = 0; ii < (int)nSurf_Poin_Par; ii++) {
    
    Global_Index = globalP[ii];
    
    /*--- Search for the processor that owns this point ---*/
    
    iProcessor = Global_Index/npoint_procs[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Linear_Elems[iProcessor])
      while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;
    
    
    if (nElem_Flag[iProcessor] != ii) {
      
      nElem_Flag[iProcessor] = ii;
      unsigned long nn = index[iProcessor];
      
      globalSend[nn] = Global_Index;
      renumbSend[nn] = renumbP[ii];
      
      /*--- Increment the index by the message length ---*/
      
      index[iProcessor]++;
      
    }
  }
  
  /*--- Free memory after loading up the send buffer. ---*/
  
  delete [] index;
  
  /*--- Allocate the memory that we need for receiving the
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  unsigned long *globalRecv = NULL;
  globalRecv = new unsigned long[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    globalRecv[ii] = 0;
  
  unsigned long *renumbRecv = NULL;
  renumbRecv = new unsigned long[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    renumbRecv[ii] = 0;
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/
  
  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];
  
  /*--- Launch the non-blocking recv's for the global ID. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(globalRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the global ID. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(globalSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking recv's for the renumbered ID. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(renumbRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the renumbered ID. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(renumbSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
  
#endif
  
  /*--- Load our own procs data into the buffers directly. ---*/
  
  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) globalRecv[mm] = globalSend[nn];
  
  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) renumbRecv[mm] = renumbSend[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*-- Now update my local connectivitiy for the surface with the new
   numbering. Create a new mapping for global -> renumber for nodes. Note
   the adding of 1 back in here for the eventual viz. purposes. ---*/
  
  map<unsigned long,unsigned long> Global2Renumber;
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    Global2Renumber[globalRecv[ii]] = renumbRecv[ii] + 1;
    Renumber2Global[renumbRecv[ii] + 1] = globalRecv[ii];
  }
  
  
  /*--- The final step is one last pass over all elements to check
   for points outside of the linear partitions of the elements. Again,
   note that elems were distributed based on their smallest global ID,
   so some nodes of the elem may have global IDs lying outside of the
   linear partitioning. We need to recover the mapping for these
   outliers. We loop over all local surface elements to find these. ---*/
  
  vector<unsigned long>::iterator it;
  vector<unsigned long> outliers;
  
  for (int ii = 0; ii < (int)nParallel_Line; ii++) {
    for ( int jj = 0; jj < N_POINTS_LINE; jj++ ) {
      
      iNode = ii*N_POINTS_LINE+jj;
      Global_Index = Conn_BoundLine_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Elems[iProcessor])
        while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;
      
      /*--- Store the global ID if it is outside our own linear partition. ---*/
      
      if ((iProcessor != (unsigned long)rank)) {
        outliers.push_back(Global_Index);
      }
      
    }
  }
  
  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;
  
  for (int ii = 0; ii < (int)nParallel_BoundTria; ii++) {
    for ( int jj = 0; jj < N_POINTS_TRIANGLE; jj++ ) {
      
      iNode = ii*N_POINTS_TRIANGLE + jj;
      Global_Index = Conn_BoundTria_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Elems[iProcessor])
        while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;
      
      /*--- Store the global ID if it is outside our own linear partition. ---*/
      
      if ((iProcessor != (unsigned long)rank)) {
        outliers.push_back(Global_Index);
      }
      
    }
  }
  
  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;
  
  for (int ii = 0; ii < (int)nParallel_BoundQuad; ii++) {
    for ( int jj = 0; jj < N_POINTS_QUADRILATERAL; jj++ ) {
      
      iNode = ii*N_POINTS_QUADRILATERAL+jj;
      Global_Index = Conn_BoundQuad_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Elems[iProcessor])
        while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;
      
      /*--- Store the global ID if it is outside our own linear partition. ---*/
      
      if ((iProcessor != (unsigned long)rank)) {
        outliers.push_back(Global_Index);
      }
      
    }
  }
  
  /*--- Create a unique list of global IDs that fall outside of our procs
   linear partition. ---*/
  
  sort(outliers.begin(), outliers.end());
  it = unique(outliers.begin(), outliers.end());
  outliers.resize(it - outliers.begin());
  
  /*--- Now loop over the outliers and communicate to those procs that
   hold the new numbering for our outlier points. We need to ask for the
   new numbering from these procs. ---*/
  
  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;
  
  for (int ii = 0; ii < (int)outliers.size(); ii++) {
    
    Global_Index = outliers[ii];
    
    /*--- Search for the processor that owns this point ---*/
    
    iProcessor = Global_Index/npoint_procs[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
      while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
    
    /*--- If we have not visited this element yet, increment our
     number of elements that must be sent to a particular proc. ---*/
    
    if ((nElem_Flag[iProcessor] != ii)) {
      nElem_Flag[iProcessor] = ii;
      nElem_Send[iProcessor+1]++;
    }
    
  }
  
  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif
  
  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
  
  nSends = 0; nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;
    
    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }
  
  delete [] idSend;
  idSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    idSend[ii] = 0;
  
  /*--- Reset our index variable for reuse. ---*/
  
  for (int ii=0; ii < size; ii++) idIndex[ii] = nElem_Send[ii];
  
  /*--- Loop over the outliers again and load up the global IDs. ---*/
  
  for (int ii = 0; ii < (int)outliers.size(); ii++) {
    
    Global_Index = outliers[ii];
    
    /*--- Search for the processor that owns this point ---*/
    
    iProcessor = Global_Index/npoint_procs[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
      while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
    
    /*--- If we have not visited this element yet, increment our
     number of elements that must be sent to a particular proc. ---*/
    
    if ((nElem_Flag[iProcessor] != ii)) {
      
      nElem_Flag[iProcessor] = ii;
      unsigned long nn = idIndex[iProcessor];
      
      /*--- Load the global ID values. ---*/
      
      idSend[nn] = Global_Index; nn++;
      
      /*--- Increment the index by the message length ---*/
      
      idIndex[iProcessor]++;
      
    }
  }
  
  /*--- Allocate the memory that we need for receiving the
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  delete [] idRecv;
  idRecv = new unsigned long[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    idRecv[ii] = 0;
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/
  
  send_req = new SU2_MPI::Request[nSends];
  recv_req = new SU2_MPI::Request[nRecvs];
  
  /*--- Launch the non-blocking recv's for the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
#endif
  
  /*--- Copy my own rank's data into the recv buffer directly. ---*/
  
  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  number = nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*--- The procs holding the outlier grid nodes now have the global IDs
   that they need to have their renumbering shared. ---*/
  
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    for (iPoint = 0; iPoint < nSurf_Poin_Par; iPoint++) {
      if (idRecv[ii] == globalP[iPoint]) {
        idRecv[ii] = renumbP[iPoint];
      }
    }
  }
  
  /*--- Now simply reverse the last communication to give the renumbered IDs
   back to the owner of the outlier points. Note everything is flipped. ---*/
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/
  
  send_req = new SU2_MPI::Request[nRecvs];
  recv_req = new SU2_MPI::Request[nSends];
  
  /*--- Launch the non-blocking sends of the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking recv's for the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
#endif
  
  /*--- Copy my own rank's data into the recv buffer directly. ---*/
  
  mm = nElem_Send[rank];
  ll = nElem_Recv[rank];
  kk = nElem_Recv[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) idSend[mm] = idRecv[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  number = nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*--- Add the renumbering for the outliers to the map from before carrying
   the global -> renumber transformation. Note that by construction, 
   nElem_Send[ii] == outliers.size(). We also add in the 1 for viz. here. ---*/
  
  for (int ii = 0; ii < nElem_Send[size]; ii++) {
    Global2Renumber[outliers[ii]] = idSend[ii] + 1;
    Renumber2Global[idSend[ii] + 1] = outliers[ii];
  }
  
  /*--- We can now overwrite the local connectivity for our surface elems
   using our completed map with the new global renumbering. Whew!! Note
   the -1 when accessing the conn from the map. ---*/
  
  for (iElem = 0; iElem < nParallel_Line; iElem++) {
    iNode = (int)iElem*N_POINTS_LINE;
    Conn_BoundLine_Par[iNode+0] = (int)Global2Renumber[Conn_BoundLine_Par[iNode+0]-1];
    Conn_BoundLine_Par[iNode+1] = (int)Global2Renumber[Conn_BoundLine_Par[iNode+1]-1];
  }
  
  for (iElem = 0; iElem < nParallel_BoundTria; iElem++) {
    iNode = (int)iElem*N_POINTS_TRIANGLE;
    Conn_BoundTria_Par[iNode+0] = (int)Global2Renumber[Conn_BoundTria_Par[iNode+0]-1];
    Conn_BoundTria_Par[iNode+1] = (int)Global2Renumber[Conn_BoundTria_Par[iNode+1]-1];
    Conn_BoundTria_Par[iNode+2] = (int)Global2Renumber[Conn_BoundTria_Par[iNode+2]-1];
  }
  
  for (iElem = 0; iElem < nParallel_BoundQuad; iElem++) {
    iNode = (int)iElem*N_POINTS_QUADRILATERAL;
    Conn_BoundQuad_Par[iNode+0] = (int)Global2Renumber[Conn_BoundQuad_Par[iNode+0]-1];
    Conn_BoundQuad_Par[iNode+1] = (int)Global2Renumber[Conn_BoundQuad_Par[iNode+1]-1];
    Conn_BoundQuad_Par[iNode+2] = (int)Global2Renumber[Conn_BoundQuad_Par[iNode+2]-1];
    Conn_BoundQuad_Par[iNode+3] = (int)Global2Renumber[Conn_BoundQuad_Par[iNode+3]-1];
  }
  
  /*--- Free temporary memory ---*/
  
  delete [] idIndex;
  delete [] surfPoint;
  delete [] globalP;
  delete [] renumbP;
  
  delete [] idSend;
  delete [] idRecv;
  delete [] globalSend;
  delete [] globalRecv;
  delete [] renumbSend;
  delete [] renumbRecv;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;
  delete [] Local_Halo;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear_Elems;
  delete [] nPoint_Linear_Nodes;
  delete [] nPoint_Send;
  delete [] nPoint_Recv;
  
}

void COutput::WriteRestart_Parallel_ASCII(CConfig *config, CGeometry *geometry) {
  
  /*--- Local variables ---*/
  
  unsigned short iVar;
  unsigned long iPoint;

  ofstream restart_file;
  string filename;
  
  int iProcessor;
  
  filename = config->GetFilename(RestartFilename, ".dat");
  
  /*--- Only the master node writes the header. ---*/
  
  if (rank == MASTER_NODE) {
    restart_file.open(filename.c_str(), ios::out);
    restart_file.precision(15);
    restart_file << "\"PointID\"";
    for (iVar = 0; iVar < Variable_Names.size()-1; iVar++)
      restart_file << "\t\"" << Variable_Names[iVar] << "\"";
    restart_file << "\t\"" << Variable_Names[Variable_Names.size()-1] << "\"" << endl;
    restart_file.close();
  }
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- All processors open the file. ---*/
  
  restart_file.open(filename.c_str(), ios::out | ios::app);
  restart_file.precision(15);
  
  /*--- Write the restart file in parallel, processor by processor. ---*/
  
  unsigned long myPoint = 0, offset = 0, Global_Index;
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
        
        /*--- Global Index of the current point. (note outer loop over procs) ---*/
        
        Global_Index = iPoint + offset;
        
        /*--- Only write original domain points, i.e., exclude any periodic
         or halo nodes, even if they are output in the viz. files. ---*/
        
        if (Global_Index < geometry->GetGlobal_nPointDomain()) {
          
          /*--- Write global index. (note outer loop over procs) ---*/
          
          restart_file << Global_Index << "\t";
          myPoint++;
          
          /*--- Loop over the variables and write the values to file ---*/
          
          for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
            restart_file << scientific << Parallel_Data[iVar][iPoint] << "\t";
          }
          restart_file << "\n";
        }
      }
    }
    /*--- Flush the file and wait for all processors to arrive. ---*/
    restart_file.flush();
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&myPoint, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
    
  }

  /*--- Write the metadata (master rank alone) ----*/

//  if (rank == MASTER_NODE) {
//    if (dual_time)
//      restart_file <<"EXT_ITER= " << config->GetExtIter() + 1 << endl;
//    else
//      restart_file <<"EXT_ITER= " << config->GetExtIter() + config->GetExtIter_OffSet() + 1 << endl;
//    restart_file <<"AOA= " << config->GetAoA() - config->GetAoA_Offset() << endl;
//    restart_file <<"SIDESLIP_ANGLE= " << config->GetAoS() - config->GetAoS_Offset() << endl;
//    restart_file <<"INITIAL_BCTHRUST= " << config->GetInitial_BCThrust() << endl;
//    restart_file <<"DCD_DCL_VALUE= " << config->GetdCD_dCL() << endl;
//    restart_file <<"DCMX_DCL_VALUE= " << config->GetdCMx_dCL() << endl;
//    restart_file <<"DCMY_DCL_VALUE= " << config->GetdCMy_dCL() << endl;
//    restart_file <<"DCMZ_DCL_VALUE= " << config->GetdCMz_dCL() << endl;

//    if (( config->GetKind_Solver() == DISC_ADJ_EULER ||
//          config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES ||
//          config->GetKind_Solver() == DISC_ADJ_RANS ) && adjoint) {
//      restart_file << "SENS_AOA=" << solver[ADJFLOW_SOL]->GetTotal_Sens_AoA() * PI_NUMBER / 180.0 << endl;
//    }
//  }

  /*--- All processors close the file. ---*/

  restart_file.close();
  
}

void COutput::WriteRestart_Parallel_Binary(CConfig *config, CGeometry *geometry) {

  /*--- Local variables ---*/

  unsigned short iVar;
  unsigned long iPoint;
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool wrt_perf  = config->GetWrt_Performance();
  ofstream restart_file;
  string filename;
  char str_buf[CGNS_STRING_SIZE], fname[100];
  su2double file_size = 0.0, StartTime, StopTime, UsedTime, Bandwidth;
  
  filename = config->GetFilename(RestartFilename, ".dat");

  strcpy(fname, filename.c_str());

  /*--- Prepare the first ints containing the counts. The first is a
   magic number that we can use to check for binary files (it is the hex
   representation for "SU2"). The second two values are number of variables
   and number of points (DoFs). The last two values are for metadata: 
   one int for ExtIter and 8 su2doubles. ---*/

  int var_buf_size = 5;
  int var_buf[5] = {535532, GlobalField_Counter, (int)nGlobalPoint_Sort, 1, 8};

  /*--- Prepare the 1D data buffer on this rank. ---*/

  passivedouble *buf = new passivedouble[nParallel_Poin*GlobalField_Counter];

  /*--- For now, create a temp 1D buffer to load up the data for writing.
   This will be replaced with a derived data type most likely. ---*/

  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++)
    for (iVar = 0; iVar < GlobalField_Counter; iVar++)
      buf[iPoint*GlobalField_Counter+iVar] = SU2_TYPE::GetValue(Parallel_Data[iVar][iPoint]);

  /*--- Prepare metadata. ---*/

  int Restart_ExtIter;
  if (dual_time)
    Restart_ExtIter= (int)config->GetExtIter() + 1;
  else
    Restart_ExtIter = (int)config->GetExtIter() + (int)config->GetExtIter_OffSet() + 1;

  passivedouble Restart_Metadata[8] = {
    SU2_TYPE::GetValue(config->GetAoA() - config->GetAoA_Offset()),
    SU2_TYPE::GetValue(config->GetAoS() - config->GetAoS_Offset()),
    SU2_TYPE::GetValue(config->GetInitial_BCThrust()),
    SU2_TYPE::GetValue(config->GetdCD_dCL()),
    SU2_TYPE::GetValue(config->GetdCMx_dCL()),
    SU2_TYPE::GetValue(config->GetdCMy_dCL()),
    SU2_TYPE::GetValue(config->GetdCMz_dCL()),
    0.0
  };

//  if (( config->GetKind_Solver() == DISC_ADJ_EULER ||
//        config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES ||
//        config->GetKind_Solver() == DISC_ADJ_RANS ) && adjoint) {
//    Restart_Metadata[4] = SU2_TYPE::GetValue(solver[ADJFLOW_SOL]->GetTotal_Sens_AoA() * PI_NUMBER / 180.0);
//  }

  /*--- Set a timer for the binary file writing. ---*/
  
#ifndef HAVE_MPI
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StartTime = MPI_Wtime();
#endif
  
#ifndef HAVE_MPI

  FILE* fhw;
  fhw = fopen(fname, "wb");

  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
  }

  /*--- First, write the number of variables and points. ---*/

  fwrite(var_buf, var_buf_size, sizeof(int), fhw);
  file_size += (su2double)var_buf_size*sizeof(int);
  
  /*--- Write the variable names to the file. Note that we are adopting a
   fixed length of 33 for the string length to match with CGNS. This is 
   needed for when we read the strings later. ---*/

  for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
    strncpy(str_buf, Variable_Names[iVar].c_str(), CGNS_STRING_SIZE);
    fwrite(str_buf, CGNS_STRING_SIZE, sizeof(char), fhw);
    file_size += (su2double)CGNS_STRING_SIZE*sizeof(char);
  }

  /*--- Call to write the entire restart file data in binary in one shot. ---*/

  fwrite(buf, GlobalField_Counter*nParallel_Poin, sizeof(passivedouble), fhw);
  file_size += (su2double)GlobalField_Counter*nParallel_Poin*sizeof(passivedouble);

  /*--- Write the external iteration. ---*/

  fwrite(&Restart_ExtIter, 1, sizeof(int), fhw);
  file_size += (su2double)sizeof(int);

  /*--- Write the metadata. ---*/

  fwrite(Restart_Metadata, 8, sizeof(passivedouble), fhw);
  file_size += (su2double)8*sizeof(passivedouble);

  /*--- Close the file. ---*/

  fclose(fhw);

#else

  /*--- Parallel binary output using MPI I/O. ---*/

  MPI_File fhw;
  SU2_MPI::Status status;
  MPI_Datatype etype, filetype;
  MPI_Offset disp;
  int ierr;

  /*--- We're writing only su2doubles in the data portion of the file. ---*/

  etype = MPI_DOUBLE;

  /*--- Define a derived datatype for this ranks contiguous chunk of data
   that will be placed in the restart (1D array size = num points * num vars). ---*/

  MPI_Type_contiguous(GlobalField_Counter*nParallel_Poin, MPI_DOUBLE, &filetype);
  MPI_Type_commit(&filetype);

  /*--- All ranks open the file using MPI. Here, we try to open the file with
   exclusive so that an error is generated if the file exists. We always want
   to write a fresh restart file, so we delete any existing files and create
   a new one. ---*/

  ierr = MPI_File_open(MPI_COMM_WORLD, fname,
                       MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                       MPI_INFO_NULL, &fhw);
  if (ierr != MPI_SUCCESS)  {
    MPI_File_close(&fhw);
    if (rank == 0)
      MPI_File_delete(fname, MPI_INFO_NULL);
    ierr = MPI_File_open(MPI_COMM_WORLD, fname,
                         MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                         MPI_INFO_NULL, &fhw);
  }

  /*--- Error check opening the file. ---*/

  if (ierr) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
  }

  /*--- First, write the number of variables and points (i.e., cols and rows),
   which we will need in order to read the file later. Also, write the 
   variable string names here. Only the master rank writes the header. ---*/

  if (rank == MASTER_NODE) {
    MPI_File_write(fhw, var_buf, var_buf_size, MPI_INT, MPI_STATUS_IGNORE);
    file_size += (su2double)var_buf_size*sizeof(int);

    /*--- Write the variable names to the file. Note that we are adopting a
     fixed length of 33 for the string length to match with CGNS. This is
     needed for when we read the strings later. ---*/

    for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
      disp = var_buf_size*sizeof(int) + iVar*CGNS_STRING_SIZE*sizeof(char);
      strcpy(str_buf, Variable_Names[iVar].c_str());
      MPI_File_write_at(fhw, disp, str_buf, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);
      file_size += (su2double)CGNS_STRING_SIZE*sizeof(char);
    }
  }

  /*--- Compute the offset for this rank's linear partition of the data in bytes.
   After the calculations above, we have the partition sizes store in nPoint_Linear
   in cumulative storage format. ---*/

  disp = (var_buf_size*sizeof(int) + GlobalField_Counter*CGNS_STRING_SIZE*sizeof(char) +
          GlobalField_Counter*nPoint_Cum[rank]*sizeof(passivedouble));

  /*--- Set the view for the MPI file write, i.e., describe the location in
   the file that this rank "sees" for writing its piece of the restart file. ---*/

  MPI_File_set_view(fhw, disp, etype, filetype, (char*)"native", MPI_INFO_NULL);

  /*--- Collective call for all ranks to write to their view simultaneously. ---*/

  MPI_File_write_all(fhw, buf, GlobalField_Counter*nParallel_Poin, MPI_DOUBLE, &status);
  file_size += (su2double)GlobalField_Counter*nParallel_Poin*sizeof(passivedouble);

  /*--- Free the derived datatype. ---*/

  MPI_Type_free(&filetype);

  /*--- Reset the file view before writing the metadata. ---*/

  MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE, (char*)"native", MPI_INFO_NULL);

  /*--- Finally, the master rank writes the metadata. ---*/

  if (rank == MASTER_NODE) {

    /*--- External iteration. ---*/

    disp = (var_buf_size*sizeof(int) + GlobalField_Counter*CGNS_STRING_SIZE*sizeof(char) +
            GlobalField_Counter*nGlobalPoint_Sort*sizeof(passivedouble));
    MPI_File_write_at(fhw, disp, &Restart_ExtIter, 1, MPI_INT, MPI_STATUS_IGNORE);
    file_size += (su2double)sizeof(int);

    /*--- Additional doubles for AoA, AoS, etc. ---*/

    disp = (var_buf_size*sizeof(int) + GlobalField_Counter*CGNS_STRING_SIZE*sizeof(char) +
            GlobalField_Counter*nGlobalPoint_Sort*sizeof(passivedouble) + 1*sizeof(int));
    MPI_File_write_at(fhw, disp, Restart_Metadata, 8, MPI_DOUBLE, MPI_STATUS_IGNORE);
    file_size += (su2double)8*sizeof(passivedouble);

  }

  /*--- All ranks close the file after writing. ---*/

  MPI_File_close(&fhw);

#endif

  /*--- Compute and store the write time. ---*/
  
#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  UsedTime = StopTime-StartTime;
  
  /*--- Communicate the total file size for the restart ---*/
  
#ifdef HAVE_MPI
  su2double my_file_size = file_size;
  SU2_MPI::Allreduce(&my_file_size, &file_size, 1,
                     MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Compute and store the bandwidth ---*/
  
  Bandwidth = file_size/(1.0e6)/UsedTime;
  config->SetRestart_Bandwidth_Agg(config->GetRestart_Bandwidth_Agg()+Bandwidth);
  
  if ((rank == MASTER_NODE) && (wrt_perf)) {
    cout << "Wrote " << file_size/1.0e6 << " MB to disk in ";
    cout << UsedTime << " s. (" << Bandwidth << " MB/s)." << endl;
  }
  
  /*--- Free temporary data buffer for writing the binary file. ---*/

  delete [] buf;

}

void COutput::WriteCSV_Slice(CConfig *config, CGeometry *geometry,
                             CSolver *FlowSolver, unsigned long iExtIter,
                             unsigned short val_iZone, unsigned short val_direction) {

  /*--- This routine is for exporting slices of field data for 2D cartesian
   grids. It assumes that the grid points lie on lines of constant x-
   or y-coordinates. It is a simple way to export slices or profiles on
   these meshes for use in verification and validation work. It will
   eventually be replaced by general routines for probing/slicing.  ---*/

  int DIRECTION = (int)val_direction;

  su2double coordMin, coordMax;
  coordMin = config->GetStations_Bounds(0);
  coordMax = config->GetStations_Bounds(1);

  unsigned short iVar;
  unsigned long iPoint, iVertex, Global_Index;
  char cstr[200];

  int rank = MASTER_NODE, iProcessor, nProcessor = SINGLE_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

  bool isPeriodic;

  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

      /*--- Checking for less than or equal to the rank, because there may
       be some periodic halo nodes that send info to the same rank. ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                      (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        if (isPeriodic) Local_Halo[iPoint] = false;
      }
    }
  }

  /*--- Sum total number of nodes that belong to the domain ---*/

  unsigned long nTotalPoint;
  unsigned long nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif

  unsigned long *npoint_procs  = new unsigned long[nProcessor];
  unsigned long *nPoint_Linear = new unsigned long[nProcessor+1];

  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < nProcessor; ii++) {
    npoint_procs[ii] = nTotalPoint/nProcessor;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }

  /*--- Get the number of remainder points after the even division. ---*/

  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }

  /*--- Store the point offsets for each rank. ---*/

  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < nProcessor; ii++) {
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[nProcessor] = nTotalPoint;

  unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = NULL;
  unsigned long nLocalVertex_Surface = 0;
  unsigned long MaxLocalVertex_Surface = 0;

  /*--- Find the max number of vertices we will send from among all
   partitions and set up buffers. The master node will handle the
   writing of the CSV file after gathering all of the data. ---*/

  nLocalVertex_Surface = 0;
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {

    /*--- Global Index of the current point. (note outer loop over procs) ---*/

    Global_Index = iPoint + nPoint_Linear[rank];

    /*--- Only write original domain points, i.e., exclude any periodic
     or halo nodes, even if they are output in the viz. files. ---*/

    if (Global_Index < geometry->GetGlobal_nPointDomain()) {
      if ((Parallel_Data[DIRECTION][iPoint] > coordMin) &&
          (Parallel_Data[DIRECTION][iPoint] < coordMax)) {
        nLocalVertex_Surface++;
      }
    }
  }

  /*--- Communicate the number of local vertices on each partition
   to the master node ---*/

  Buffer_Send_nVertex[0] = nLocalVertex_Surface;
  if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long [nProcessor];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalVertex_Surface = nLocalVertex_Surface;
  Buffer_Recv_nVertex[0] = Buffer_Send_nVertex[0];
#endif

  /*--- Send and Recv buffers ---*/

  su2double *Buffer_Send_Data = new su2double [MaxLocalVertex_Surface*GlobalField_Counter];
  su2double *Buffer_Recv_Data = NULL;

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long [MaxLocalVertex_Surface];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Prepare the receive buffers on the master node only. ---*/

  if (rank == MASTER_NODE) {
    Buffer_Recv_Data = new su2double [nProcessor*MaxLocalVertex_Surface*GlobalField_Counter];
    Buffer_Recv_GlobalIndex  = new unsigned long [nProcessor*MaxLocalVertex_Surface];
  }

  /*--- Loop over all vertices in this partition and load the
   data of the specified type into the buffer to be sent to
   the master node. ---*/

  nLocalVertex_Surface = 0;
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {

    /*--- Global Index of the current point. (note outer loop over procs) ---*/

    Global_Index = iPoint + nPoint_Linear[rank];

    /*--- Only write original domain points, i.e., exclude any periodic
     or halo nodes, even if they are output in the viz. files. ---*/

    if (Global_Index < geometry->GetGlobal_nPointDomain()) {
      if ((Parallel_Data[DIRECTION][iPoint] > coordMin) &&
          (Parallel_Data[DIRECTION][iPoint] < coordMax)) {
        Buffer_Send_GlobalIndex[nLocalVertex_Surface] = Global_Index;
        for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
          Buffer_Send_Data[nLocalVertex_Surface*GlobalField_Counter+iVar] = Parallel_Data[iVar][iPoint];
        }
        nLocalVertex_Surface++;
      }
    }
  }

  /*--- Send the information to the master node ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Data, MaxLocalVertex_Surface*GlobalField_Counter, MPI_DOUBLE, Buffer_Recv_Data, MaxLocalVertex_Surface*GlobalField_Counter, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (iVertex = 0; iVertex < Buffer_Recv_nVertex[0]; iVertex++) {
    Buffer_Recv_GlobalIndex[iVertex] = Buffer_Send_GlobalIndex[iVertex];
    for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
      Buffer_Recv_Data[iVertex*GlobalField_Counter+iVar] = Buffer_Send_Data[iVertex*GlobalField_Counter+iVar];
    }
  }
#endif

  /*--- The master node unpacks the data and writes the surface CSV file ---*/

  if (rank == MASTER_NODE) {

    /*--- Write file name with extension if unsteady ---*/
    char buffer[50];
    string filename = "slice";
    if (DIRECTION == 0) {
      SPRINTF (buffer, "_vert.csv");
    } else if (DIRECTION == 1) {
      SPRINTF (buffer, "_hori.csv");
    }
    ofstream SurfFlow_file;

    /*--- Write file name with extension if unsteady ---*/
    strcpy (cstr, filename.c_str());
    strcat (cstr, buffer);
    SurfFlow_file.precision(15);
    SurfFlow_file.open(cstr, ios::out);

    /*--- Global index is first, then the rest of the data. We already have the names. ---*/

    SurfFlow_file << "\"Global_Index\",";
    for (iVar = 0; iVar < Variable_Names.size()-1; iVar++) {
      SurfFlow_file << "\"" << Variable_Names[iVar] << "\",";
    }
    SurfFlow_file << "\"" << Variable_Names[Variable_Names.size()-1] << "\"" << endl;

    /*--- Loop through all of the collected data and write each node's values ---*/

    unsigned long Total_Index;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {

        /*--- Current index position and global index ---*/

        Total_Index  = iProcessor*MaxLocalVertex_Surface+iVertex;
        Global_Index = Buffer_Recv_GlobalIndex[Total_Index];

        /*--- Write the the data ---*/

        SurfFlow_file << scientific << Global_Index;
        Total_Index  = iProcessor*MaxLocalVertex_Surface*GlobalField_Counter+iVertex*GlobalField_Counter;
        for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
          SurfFlow_file << scientific << ", " << Buffer_Recv_Data[Total_Index+iVar];
        }
        SurfFlow_file << endl;

      }
    }

    /*--- Close the CSV file ---*/

    SurfFlow_file.close();

    /*--- Release the recv buffers on the master node ---*/

    delete [] Buffer_Recv_Data;
    delete [] Buffer_Recv_GlobalIndex;
    delete [] Buffer_Recv_nVertex;

  }

  /*--- Release the memory for the remaining buffers and exit ---*/

  delete [] Buffer_Send_Data;
  delete [] Buffer_Send_GlobalIndex;
  
}

void COutput::DeallocateConnectivity_Parallel(bool surf_sol) {
  
  /*--- Deallocate memory for connectivity data on each processor. ---*/
  
  if (surf_sol) {
    if (nParallel_Line > 0      && Conn_BoundLine_Par      != NULL)
      delete [] Conn_BoundLine_Par;
    if (nParallel_BoundTria > 0 && Conn_BoundTria_Par != NULL)
      delete [] Conn_BoundTria_Par;
    if (nParallel_BoundQuad > 0 && Conn_BoundQuad_Par != NULL)
      delete [] Conn_BoundQuad_Par;
  }
  else {
    if (nParallel_Tria > 0 && Conn_Tria_Par != NULL) delete [] Conn_Tria_Par;
    if (nParallel_Quad > 0 && Conn_Quad_Par != NULL) delete [] Conn_Quad_Par;
    if (nParallel_Tetr > 0 && Conn_Tetr_Par != NULL) delete [] Conn_Tetr_Par;
    if (nParallel_Hexa > 0 && Conn_Hexa_Par != NULL) delete [] Conn_Hexa_Par;
    if (nParallel_Pris > 0 && Conn_Pris_Par != NULL) delete [] Conn_Pris_Par;
    if (nParallel_Pyra > 0 && Conn_Pyra_Par != NULL) delete [] Conn_Pyra_Par;
  }
  
}

void COutput::DeallocateData_Parallel() {
  
  /*--- Deallocate memory for solution data ---*/
  
  for (unsigned short iVar = 0; iVar < GlobalField_Counter; iVar++) {
    if (Parallel_Data[iVar] != NULL) delete [] Parallel_Data[iVar];
  }
  if (Parallel_Data != NULL) delete [] Parallel_Data;

}

void COutput::DeallocateSurfaceData_Parallel() {
  
  if (Parallel_Surf_Data != NULL) {
    
    Global2Renumber.clear();
    Renumber2Global.clear();
    /*--- Deallocate memory for surface solution data ---*/
    
    for (unsigned short iVar = 0; iVar < GlobalField_Counter; iVar++) {
      if (Parallel_Surf_Data[iVar] != NULL) delete [] Parallel_Surf_Data[iVar];
    }
    delete [] Parallel_Surf_Data;
  }
  
}

void COutput::MergeInletCoordinates(CConfig *config, CGeometry *geometry) {

  /*--- Local variables needed on all processors ---*/

  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, jPoint, kPoint;

  int iProcessor, nProcessor = size;

  unsigned long iVertex, iMarker;
  unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;

  unsigned long index, iChar;

  char str_buf[MAX_STRING_SIZE];
  vector<string> Marker_Tags;

  unsigned long *nRowCum_Counter = NULL;

  if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];

  /*--- Search all boundaries on the present rank to count the number
   of nodes found on inlet markers. ---*/

  nLocalPoint = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Only communicate owned nodes to avoid duplicates. ---*/

        if (geometry->node[iPoint]->GetDomain())
          nLocalPoint++;

      }
    }
  }
  Buffer_Send_nPoin[0] = nLocalPoint;

  /*--- Communicate the total number of nodes on this domain. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(&Buffer_Send_nPoin, 1, MPI_UNSIGNED_LONG,
                  Buffer_Recv_nPoin, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
#else
  Buffer_Recv_nPoin[0] = Buffer_Send_nPoin[0];
  MaxLocalPoint        = nLocalPoint;
#endif

  /*--- Send and Recv buffers. ---*/

  su2double *Buffer_Send_X = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_X = NULL;

  su2double *Buffer_Send_Y = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Y = NULL;

  su2double *Buffer_Send_Z = NULL, *Buffer_Recv_Z = NULL;
  if (nDim == 3) Buffer_Send_Z = new su2double[MaxLocalPoint];

  char *Buffer_Send_Str = new char[MaxLocalPoint*MAX_STRING_SIZE];
  char *Buffer_Recv_Str = NULL;

  /*--- Prepare the receive buffers in the master node only. ---*/

  if (rank == MASTER_NODE) {

    Buffer_Recv_X = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Y = new su2double[nProcessor*MaxLocalPoint];
    if (nDim == 3) Buffer_Recv_Z = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Str = new char[nProcessor*MaxLocalPoint*MAX_STRING_SIZE];

    /*--- Sum total number of nodes to be written and allocate arrays ---*/

    unsigned long nGlobal_InletPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_InletPoint += Buffer_Recv_nPoin[iProcessor];
    }
    InletCoords = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      InletCoords[iDim] = new su2double[nGlobal_InletPoint];
    }
  }

  /*--- Main communication routine. Loop over each coordinate and perform
   the MPI comm. Temporary 1-D buffers are used to send the coordinates at
   all nodes on each partition to the master node. These are then unpacked
   by the master and sorted by marker tag in one large n-dim. array. ---*/

  /*--- Loop over this partition to collect the coords of the local points. ---*/

  su2double *Coords_Local; jPoint = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Only communicate owned nodes to avoid duplicates. ---*/

        if (geometry->node[iPoint]->GetDomain()) {

          /*--- Retrieve local coordinates at this node. ---*/

          Coords_Local = geometry->node[iPoint]->GetCoord();

          /*--- Load local coords into the temporary send buffer. ---*/

          Buffer_Send_X[jPoint] = Coords_Local[0];
          Buffer_Send_Y[jPoint] = Coords_Local[1];
          if (nDim == 3) Buffer_Send_Z[jPoint] = Coords_Local[2];

          /*--- If US system, the output should be in inches ---*/

          if (config->GetSystemMeasurements() == US) {
            Buffer_Send_X[jPoint] *= 12.0;
            Buffer_Send_Y[jPoint] *= 12.0;
            if (nDim == 3) Buffer_Send_Z[jPoint] *= 12.0;
          }

          /*--- Store the marker tag for this particular node. ---*/

          SPRINTF(&Buffer_Send_Str[jPoint*MAX_STRING_SIZE], "%s",
                  config->GetMarker_All_TagBound(iMarker).c_str());

          /*--- Increment jPoint as the counter. We need this because iPoint
           may include halo nodes that we skip over during this loop. ---*/
          
          jPoint++;
          
        }
      }
    }
  }

  /*--- Gather the coordinate data on the master node using MPI. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_X, MaxLocalPoint, MPI_DOUBLE, Buffer_Recv_X, MaxLocalPoint, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Y, MaxLocalPoint, MPI_DOUBLE, Buffer_Recv_Y, MaxLocalPoint, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) {
    SU2_MPI::Gather(Buffer_Send_Z, MaxLocalPoint, MPI_DOUBLE, Buffer_Recv_Z, MaxLocalPoint, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  }
  SU2_MPI::Gather(Buffer_Send_Str, MaxLocalPoint*MAX_STRING_SIZE, MPI_CHAR, Buffer_Recv_Str, MaxLocalPoint*MAX_STRING_SIZE, MPI_CHAR, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < MaxLocalPoint; iPoint++) {
    Buffer_Recv_X[iPoint] = Buffer_Send_X[iPoint];
    Buffer_Recv_Y[iPoint] = Buffer_Send_Y[iPoint];
    if (nDim == 3) Buffer_Recv_Z[iPoint] = Buffer_Send_Z[iPoint];
    index = iPoint*MAX_STRING_SIZE;
    for (iChar = 0; iChar < MAX_STRING_SIZE; iChar++) {
      Buffer_Recv_Str[index + iChar] = Buffer_Send_Str[index + iChar];
    }
  }
#endif

  /*--- The master node unpacks and sorts this variable by marker tag. ---*/

  if (rank == MASTER_NODE) {

    Marker_Tags_InletFile.clear();

    /*--- First, parse the marker tags to count how many total inlet markers
     we have now on the master. ---*/

    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        index = (iProcessor*MaxLocalPoint + iPoint)*MAX_STRING_SIZE;
        for (iChar = 0; iChar < MAX_STRING_SIZE; iChar++) {
          str_buf[iChar] = Buffer_Recv_Str[index + iChar];
        }
        Marker_Tags.push_back(str_buf);
        Marker_Tags_InletFile.push_back(str_buf);
      }
    }

    /*--- Sort and remove the duplicate inlet marker strings. ---*/

    sort(Marker_Tags_InletFile.begin(), Marker_Tags_InletFile.end());
    Marker_Tags_InletFile.erase(unique(Marker_Tags_InletFile.begin(),
                                       Marker_Tags_InletFile.end()),
                                Marker_Tags_InletFile.end());

    /*--- Store the unique number of markers for writing later. ---*/

    nMarker_InletFile = Marker_Tags_InletFile.size();

    /*--- Count the number of rows (nodes) per marker. ---*/

    nRow_InletFile = new unsigned long[nMarker_InletFile];
    for (iMarker = 0; iMarker < nMarker_InletFile; iMarker++) {
      nRow_InletFile[iMarker] = 0;
    }

    /*--- Now count the number of points per marker. ---*/

    jPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        for (iMarker = 0; iMarker < nMarker_InletFile; iMarker++) {
          if (Marker_Tags_InletFile[iMarker] == Marker_Tags[jPoint]) {
            nRow_InletFile[iMarker]++;
          }
        }
        jPoint++;
      }
    }

    /*--- Now put the number of points per marker into cumulative storage.
     We will also create an extra counter to make sorting easier. ---*/

    nRowCum_InletFile = new unsigned long[nMarker_InletFile+1];
    nRowCum_Counter   = new unsigned long[nMarker_InletFile+1];

    nRowCum_InletFile[0] = 0; nRowCum_Counter[0] = 0;
    for (iMarker = 0; iMarker < nMarker_InletFile; iMarker++) {
      nRowCum_InletFile[iMarker+1] = nRowCum_InletFile[iMarker] + nRow_InletFile[iMarker];
      nRowCum_Counter[iMarker+1]   = nRowCum_Counter[iMarker]   + nRow_InletFile[iMarker];
    }

    /*--- Load up the coordinates, sorted into chunks per marker. ---*/

    jPoint = 0; kPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        for (iMarker = 0; iMarker < nMarker_InletFile; iMarker++) {
          if (Marker_Tags_InletFile[iMarker] == Marker_Tags[kPoint]) {

            /*--- Find our current index for this marker and store coords. ---*/

            index = nRowCum_Counter[iMarker];
            InletCoords[0][index] = Buffer_Recv_X[jPoint];
            InletCoords[1][index] = Buffer_Recv_Y[jPoint];
            if (nDim == 3) InletCoords[2][index] = Buffer_Recv_Z[jPoint];

            /*--- Increment the counter for this marker. ---*/

            nRowCum_Counter[iMarker]++;

          }
        }

        /*--- Increment point counter for marker tags and data. ---*/

        kPoint++;
        jPoint++;
        
      }

      /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

      jPoint = (iProcessor+1)*MaxLocalPoint;

    }
  }

  /*--- Immediately release the temporary data buffers. ---*/

  delete [] Buffer_Send_X;
  delete [] Buffer_Send_Y;
  if (Buffer_Send_Z != NULL) delete [] Buffer_Send_Z;
  delete [] Buffer_Send_Str;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_X;
    delete [] Buffer_Recv_Y;
    if (Buffer_Recv_Z != NULL)  delete [] Buffer_Recv_Z;
    delete [] Buffer_Recv_nPoin;
    delete [] Buffer_Recv_Str;
    delete [] nRowCum_Counter;
  }

}

void COutput::Write_InletFile_Flow(CConfig *config, CGeometry *geometry, CSolver **solver) {

  unsigned short iMarker, iDim, iVar;
  unsigned long iPoint;
  su2double turb_val[2] = {0.0,0.0};

  const unsigned short nDim = geometry->GetnDim();

  bool turbulent = (config->GetKind_Solver() == RANS ||
                    config->GetKind_Solver() == ADJ_RANS ||
                    config->GetKind_Solver() == DISC_ADJ_RANS);

  unsigned short nVar_Turb = 0;
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
        nVar_Turb = 1;
        turb_val[0] = solver[TURB_SOL]->GetNuTilde_Inf();
        break;
      case SST:
        nVar_Turb = 2;
        turb_val[0] = solver[TURB_SOL]->GetTke_Inf();
        turb_val[1] = solver[TURB_SOL]->GetOmega_Inf();
        break;
      default:
        SU2_MPI::Error("Specified turbulence model unavailable or none selected", CURRENT_FUNCTION);
        break;
    }

  /*--- Count the number of columns that we have for this flow case.
   Here, we have nDim entries for node coordinates, 2 entries for the total
   conditions or mass flow, another nDim for the direction vector, and
   finally entries for the number of turbulence variables. ---*/

  unsigned short nCol_InletFile = nDim + 2 + nDim + nVar_Turb;

  /*--- Write the inlet profile file. Note that we have already merged
   all of the information for the markers and coordinates previously
   in the MergeInletCoordinates() routine. ---*/

  ofstream node_file("inlet_example.dat");

  node_file << "NMARK= " << nMarker_InletFile << endl;

  for (iMarker = 0; iMarker < nMarker_InletFile; iMarker++) {

    /*--- Access the default data for this marker. ---*/

    string Marker_Tag   = Marker_Tags_InletFile[iMarker];
    su2double p_total   = config->GetInlet_Ptotal(Marker_Tag);
    su2double t_total   = config->GetInlet_Ttotal(Marker_Tag);
    su2double* flow_dir = config->GetInlet_FlowDir(Marker_Tag);

    /*--- Header information for this marker. ---*/

    node_file << "MARKER_TAG= " << Marker_Tag              << endl;
    node_file << "NROW="        << nRow_InletFile[iMarker] << endl;
    node_file << "NCOL="        << nCol_InletFile          << endl;

    node_file << setprecision(15);
    node_file << std::scientific;

    /*--- Loop over the data structure and write the coords and vars to file. ---*/

    for (iPoint = nRowCum_InletFile[iMarker]; iPoint < nRowCum_InletFile[iMarker+1]; iPoint++) {

      for (iDim = 0; iDim < nDim; iDim++) {
        node_file << InletCoords[iDim][iPoint] << "\t";
      }
      node_file << t_total << "\t" << p_total;
      for (iDim = 0; iDim < nDim; iDim++) {
        node_file << "\t" << flow_dir[iDim];
      }
      for (iVar = 0; iVar < nVar_Turb; iVar++) {
        node_file << "\t" << turb_val[iVar];
      }
      node_file << endl;
    }

  }
  node_file.close();

  /*--- Print a message to inform the user about the template file. ---*/
  
  stringstream err;
  err << endl;
  err << "  Created a template inlet profile file with node coordinates" << endl;
  err << "  and solver variables at `inlet_example.dat`." << endl;
  err << "  You can use this file as a guide for making your own inlet" << endl;
  err << "  specification." << endl << endl;
  SU2_MPI::Error(err.str(), CURRENT_FUNCTION);

}

void COutput::DeallocateInletCoordinates(CConfig *config, CGeometry *geometry) {

  unsigned short iDim, nDim = geometry->GetnDim();

  /*--- The master node alone owns all data found in this routine. ---*/

  if (rank == MASTER_NODE) {

    /*--- Deallocate memory for inlet coordinate data ---*/

    if (nRow_InletFile != NULL)    delete [] nRow_InletFile;
    if (nRowCum_InletFile != NULL) delete [] nRowCum_InletFile;

    Marker_Tags_InletFile.clear();

    for (iDim = 0; iDim < nDim; iDim++) {
      if (InletCoords[iDim] != NULL) delete [] InletCoords[iDim];
    }
    if (InletCoords != NULL) delete [] InletCoords;
  }

}




void COutput::PrepareOffsets(CConfig *config, CGeometry *geometry) {

  unsigned long iPoint;


  /*--- Reset point sorting counters ---*/

  nGlobalPoint_Sort = 0;
  nLocalPoint_Sort  = 0;

  /*--- Prepare the offsets for the FV solver ---*/

  if (!fem_output) {

    /*--- Search all send/recv boundaries on this partition for any periodic
     nodes that were part of the original domain. We want to recover these
     for visualization purposes. ---*/

    unsigned long iVertex;
    bool isPeriodic;

    Local_Halo_Sort = new int[geometry->GetnPoint()];
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      Local_Halo_Sort[iPoint] = !geometry->node[iPoint]->GetDomain();

    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo_Sort[iPoint] = false;
        }
      }
    }

    /*--- Sum total number of nodes that belong to the domain ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo_Sort[iPoint] == false)
        nLocalPoint_Sort++;

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalPoint_Sort, &nGlobalPoint_Sort, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
    nGlobalPoint_Sort = nLocalPoint_Sort;
#endif

    /*--- Set a final variable for the number of points in the restart
     file. We do not write the periodic points for the FV solver, even if
     they show up in the viz. files. ---*/

    nPoint_Restart = geometry->GetGlobal_nPointDomain();

  } else {

    /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
     geometrical information for the FEM DG solver. ---*/

    CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

    unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
    CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

    /*--- Update the solution by looping over the owned volume elements. ---*/

    for(unsigned long l=0; l<nVolElemOwned; ++l) {

      /* Count up the number of local points we have for allocating storage. */

      for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
        nLocalPoint_Sort++;
      }
    }

    Local_Halo_Sort = new int[nLocalPoint_Sort];
    for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++)
      Local_Halo_Sort[iPoint] = false;

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalPoint_Sort, &nGlobalPoint_Sort, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
    nGlobalPoint_Sort = nLocalPoint_Sort;
#endif

    /*--- Set a final variable for the number of points in the restart
     file. We do not write the periodic points for the FV solver, even if
     they show up in the viz. files. ---*/

    nPoint_Restart = nGlobalPoint_Sort;

  }

  /*--- Now that we know the actual number of points we need to output,
   compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/

  beg_node = new unsigned long[size];
  end_node = new unsigned long[size];

  nPoint_Lin = new unsigned long[size];
  nPoint_Cum = new unsigned long[size+1];

  unsigned long total_points = 0;
  for (int ii = 0; ii < size; ii++) {
    nPoint_Lin[ii] = nGlobalPoint_Sort/size;
    total_points  += nPoint_Lin[ii];
  }

  /*--- Get the number of remainder points after the even division. ---*/

  unsigned long remainder = nGlobalPoint_Sort - total_points;
  for (unsigned long ii = 0; ii < remainder; ii++) {
    nPoint_Lin[ii]++;
  }

  /*--- Store the local number of nodes on each proc in the linear
   partitioning, the beginning/end index, and the linear partitioning
   within an array in cumulative storage format. ---*/

  beg_node[0] = 0;
  end_node[0] = beg_node[0] + nPoint_Lin[0];
  nPoint_Cum[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    beg_node[ii]   = end_node[ii-1];
    end_node[ii]   = beg_node[ii] + nPoint_Lin[ii];
    nPoint_Cum[ii] = nPoint_Cum[ii-1] + nPoint_Lin[ii-1];
  }
  nPoint_Cum[size] = nGlobalPoint_Sort;
  
}

void COutput::SortConnectivity_FEM(CConfig *config, CGeometry *geometry, bool surf) {


  /*--- Sort connectivity for each type of element (excluding halos). Note
   In these routines, we sort the connectivity into a linear partitioning
   across all processors based on the global index of the grid nodes. ---*/

  /*--- Sort volumetric grid connectivity. ---*/

  if (!surf) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
      cout <<"Sorting volumetric grid connectivity." << endl;

    SortVolumetricConnectivity_FEM(config, geometry, TRIANGLE     );
    SortVolumetricConnectivity_FEM(config, geometry, QUADRILATERAL);
    SortVolumetricConnectivity_FEM(config, geometry, TETRAHEDRON  );
    SortVolumetricConnectivity_FEM(config, geometry, HEXAHEDRON   );
    SortVolumetricConnectivity_FEM(config, geometry, PRISM        );
    SortVolumetricConnectivity_FEM(config, geometry, PYRAMID      );
    
    /*--- Reduce the total number of cells we will be writing in the output files. ---*/
  
    unsigned long nTotal_Elem = nParallel_Tria + nParallel_Quad + nParallel_Tetr + nParallel_Hexa + nParallel_Pris + nParallel_Pyra;
  #ifndef HAVE_MPI
    nGlobal_Elem_Par = nTotal_Elem;
  #else
    SU2_MPI::Allreduce(&nTotal_Elem, &nGlobal_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  #endif
  }

  /*--- Sort surface grid connectivity. ---*/

  else {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
      cout <<"Sorting surface grid connectivity." << endl;

    SortSurfaceConnectivity_FEM(config, geometry, LINE         );
    SortSurfaceConnectivity_FEM(config, geometry, TRIANGLE     );
    SortSurfaceConnectivity_FEM(config, geometry, QUADRILATERAL);
    
    /*--- Reduce the total number of cells we will be writing in the output files. ---*/
  
    unsigned long nTotal_Surf_Elem = nParallel_Line + nParallel_BoundTria + nParallel_BoundQuad;
  #ifndef HAVE_MPI
    nSurf_Elem_Par   = nTotal_Surf_Elem;
  #else
    SU2_MPI::Allreduce(&nTotal_Surf_Elem, &nSurf_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  #endif
  }
}

void COutput::SortVolumetricConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

}

void COutput::SortSurfaceConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

}

void COutput::SortOutputData_FEM(CConfig *config, CGeometry *geometry) {

  unsigned long iProcessor;
  unsigned long iPoint, Global_Index;

  /* For convenience, set the total number of variables stored at each DOF. */

  int VARS_PER_POINT = GlobalField_Counter;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/

  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

  /*--- Create the map from the global DOF ID to the local index. ---*/

  //map<unsigned long, unsigned long> mapLocal2Global;
  vector<unsigned long> globalID;

  /*--- Update the solution by looping over the owned volume elements. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Loop over the DOFs for this element and store the solution. */
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      globalID.push_back(globalIndex);
    }
  }

  /*--- We start with the grid nodes distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many nodes we must send to each other rank in order to
   have all nodes sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first ~ nGlobalPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/

  int *nPoint_Send = new int[size+1]; nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1]; nPoint_Recv[0] = 0;
  int *nPoint_Flag = new int[size];

  for (int ii=0; ii < size; ii++) {
    nPoint_Send[ii] = 0;
    nPoint_Recv[ii] = 0;
    nPoint_Flag[ii]= -1;
  }
  nPoint_Send[size] = 0; nPoint_Recv[size] = 0;

  for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++ ) {

    /*--- Get the global index of the current point. ---*/

    Global_Index = globalID[iPoint];

    /*--- Search for the processor that owns this point ---*/

    iProcessor = Global_Index/nPoint_Lin[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Cum[iProcessor])
      while(Global_Index >= nPoint_Cum[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Cum[iProcessor])   iProcessor--;

    /*--- If we have not visted this node yet, increment our
     number of elements that must be sent to a particular proc. ---*/

    if (nPoint_Flag[iProcessor] != (int)iPoint) {
      nPoint_Flag[iProcessor] = (int)iPoint;
      nPoint_Send[iProcessor+1]++;
    }

  }

  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nPoint_Recv[1] = nPoint_Send[1];
#endif

  /*--- Prepare to send coordinates. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nPoint_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nPoint_Recv[ii+1] > 0)) nRecvs++;

    nPoint_Send[ii+1] += nPoint_Send[ii];
    nPoint_Recv[ii+1] += nPoint_Recv[ii];
  }

  /*--- Allocate memory to hold the connectivity that we are sending. ---*/

  su2double *connSend = NULL;
  connSend = new su2double[VARS_PER_POINT*nPoint_Send[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Send[size]; ii++)
    connSend[ii] = 0;

  /*--- Allocate arrays for sending the global ID. ---*/

  unsigned long *idSend = new unsigned long[nPoint_Send[size]];
  for (int ii = 0; ii < nPoint_Send[size]; ii++)
    idSend[ii] = 0;

  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = VARS_PER_POINT*nPoint_Send[ii];

  unsigned long *idIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) idIndex[ii] = nPoint_Send[ii];

  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/

  for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++) {

    /*--- Get the index of the current point. ---*/

    Global_Index = globalID[iPoint];

    /*--- Search for the processor that owns this point. ---*/

    iProcessor = Global_Index/nPoint_Lin[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Cum[iProcessor])
      while(Global_Index >= nPoint_Cum[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Cum[iProcessor])   iProcessor--;

    /*--- Load data into the buffer for sending. ---*/

    if (nPoint_Flag[iProcessor] != (int)iPoint) {

      nPoint_Flag[iProcessor] = (int)iPoint;
      unsigned long nn = index[iProcessor];

      /*--- Load the data values. ---*/

      for (unsigned short kk = 0; kk < VARS_PER_POINT; kk++) {
        connSend[nn] = Local_Data[iPoint][kk]; nn++;
      }

      /*--- Load the global ID (minus offset) for sorting the
       points once they all reach the correct processor. ---*/

      nn = idIndex[iProcessor];
      idSend[nn] = Global_Index - beg_node[iProcessor];

      /*--- Increment the index by the message length ---*/

      index[iProcessor]  += VARS_PER_POINT;
      idIndex[iProcessor]++;

    }

  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;
  delete [] idIndex;

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  su2double *connRecv = NULL;
  connRecv = new su2double[VARS_PER_POINT*nPoint_Recv[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Recv[size]; ii++)
    connRecv[ii] = 0;

  unsigned long *idRecv = new unsigned long[nPoint_Recv[size]];
  for (int ii = 0; ii < nPoint_Recv[size]; ii++)
    idRecv[ii] = 0;

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the global IDs. ---*/

  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];

  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = VARS_PER_POINT*nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = VARS_PER_POINT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_DOUBLE, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = VARS_PER_POINT*nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = VARS_PER_POINT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_DOUBLE, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Repeat the process to communicate the global IDs. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the global IDs. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  int mm = VARS_PER_POINT*nPoint_Recv[rank];
  int ll = VARS_PER_POINT*nPoint_Send[rank];
  int kk = VARS_PER_POINT*nPoint_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];

  mm = nPoint_Recv[rank];
  ll = nPoint_Send[rank];
  kk = nPoint_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. First, allocate the
   appropriate amount of memory for this section. ---*/

  Parallel_Data = new su2double*[VARS_PER_POINT];
  for (int jj = 0; jj < VARS_PER_POINT; jj++) {
    Parallel_Data[jj] = new su2double[nPoint_Recv[size]];
    for (int ii = 0; ii < nPoint_Recv[size]; ii++) {
      Parallel_Data[jj][idRecv[ii]] = connRecv[ii*VARS_PER_POINT+jj];
    }
  }

  /*--- Store the total number of local points my rank has for
   the current section after completing the communications. ---*/

  nParallel_Poin = nPoint_Recv[size];

  /*--- Reduce the total number of points we will write in the output files. ---*/

#ifndef HAVE_MPI
  nGlobal_Poin_Par = nParallel_Poin;
#else
  SU2_MPI::Allreduce(&nParallel_Poin, &nGlobal_Poin_Par, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Free temporary memory from communications ---*/

  delete [] connSend;
  delete [] connRecv;
  delete [] idSend;
  delete [] idRecv;
  delete [] nPoint_Recv;
  delete [] nPoint_Send;
  delete [] nPoint_Flag;
  
}

void COutput::SortOutputData_Surface_FEM(CConfig *config, CGeometry *geometry) {

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
  
  cout << "History filename: " << char_histfile << endl;
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
        cout << "  Info: Ignoring the following history output fields/groups:" << endl;
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
  
  
  /*--- First, prepare the offsets needed throughout below. ---*/    
  
  PrepareOffsets(config, geometry);    
  
  /*--- Now that we know the number of fields, create the local data array to temporarily store the volume output 
   * before writing it to file ---*/
   
  Local_Data.resize(nLocalPoint_Sort, std::vector<su2double>(GlobalField_Counter, 0.0));
  
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
      
      if (!Local_Halo_Sort[iPoint] || Wrt_Halo) {
        
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
}

void COutput::Postprocess_HistoryData(CConfig *config){
   
  map<string, su2double> Average;
  map<string, int> Count;
  
  for (unsigned short iField = 0; iField < HistoryOutput_List.size(); iField++){
    HistoryOutputField &currentField = HistoryOutput_Map[HistoryOutput_List[iField]];
    if (currentField.FieldType == TYPE_RESIDUAL){
      if ( SetInit_Residuals(config) ) {
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
      AddHistoryOutput("REL_" + HistoryOutput_List[iField], "rel" + currentField.FieldName, currentField.ScreenFormat, "REL_" + currentField.OutputGroup);
      Average[currentField.OutputGroup] = true;
    }
    if (currentField.FieldType == TYPE_COEFFICIENT){
      AddHistoryOutput("TAVG_"   + HistoryOutput_List[iField], "tavg["  + currentField.FieldName + "]", currentField.ScreenFormat, "TAVG_"   + currentField.OutputGroup);
      AddHistoryOutput("D_"      + HistoryOutput_List[iField], "d["     + currentField.FieldName + "]", currentField.ScreenFormat, "D_"      + currentField.OutputGroup);  
      AddHistoryOutput("D_TAVG_" + HistoryOutput_List[iField], "dtavg[" + currentField.FieldName + "]", currentField.ScreenFormat, "D_TAVG_" + currentField.OutputGroup);  
    }
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

