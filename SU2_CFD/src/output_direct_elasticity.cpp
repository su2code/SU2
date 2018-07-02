/*!
 * \file output_direct_elasticity.cpp
 * \brief Main subroutines for FEA output
 * \author R. Sanchez
 * \version 6.0.1 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include "../include/output_structure.hpp"

CFEAOutput::CFEAOutput(CConfig *config, CGeometry *geometry, unsigned short val_iZone) : COutput(config) {

  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  
  /*--- Initialize number of variables ---*/
  if (linear_analysis) nVar_FEM = nDim;
  if (nonlinear_analysis) nVar_FEM = 3;
  
  nDim = geometry->GetnDim();

}

CFEAOutput::~CFEAOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();

  }

}

void CFEAOutput::LoadOutput_Data(CGeometry ****geometry,
                                     CSolver *****solver_container,
                                     CConfig **config,
                                     CIntegration ****integration,
                                     bool DualTime_Iteration,
                                     su2double timeused,
                                     unsigned short val_iZone,
                                     unsigned short val_iInst) {

  bool fem = ((config[val_iZone]->GetKind_Solver() == FEM_ELASTICITY) ||          // FEM structural solver.
              (config[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM));
  bool linear_analysis = (config[val_iZone]->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config[val_iZone]->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  bool fsi = (config[val_iZone]->GetFSI_Simulation());          // FEM structural solver.
  bool discadj_fem = (config[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM);

  unsigned short iVar;
  unsigned short nDim = geometry[val_iZone][INST_0][MESH_0]->GetnDim();

  SetOutputFieldValue("INT_ITER", config[val_iZone]->GetIntIter());
  SetOutputFieldValue("EXT_ITER", config[val_iZone]->GetExtIter());
  
  SetOutputFieldValue("PHYS_TIME", timeused);
  
  /*--- Residuals: ---*/
  /*--- Linear analysis: RMS of the displacements in the nDim coordinates ---*/
  /*--- Nonlinear analysis: UTOL, RTOL and DTOL (defined in the Postprocessing function) ---*/

  
  if (linear_analysis){
    SetOutputFieldValue("UTOL", log10(solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_RMS(0)));
    SetOutputFieldValue("RTOL", log10(solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_RMS(1)));
    if (nVar_FEM == 3){
      SetOutputFieldValue("ETOL", log10(solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_RMS(2)));    
    }
    SetOutputFieldValue("DISP_X", log10(solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_RMS(0)));
    SetOutputFieldValue("DISP_Y", log10(solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_RMS(1)));
    if (nVar_FEM == 3){
      SetOutputFieldValue("DISP_Z", log10(solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_RMS(2)));    
    }
  } else if (nonlinear_analysis){
    SetOutputFieldValue("UTOL", log10(solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_FEM(0)));
    SetOutputFieldValue("RTOL", log10(solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_FEM(1)));
    if (nVar_FEM == 3){
      SetOutputFieldValue("ETOL", log10(solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_FEM(2)));    
    }
    SetOutputFieldValue("DISP_X", log10(solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_FEM(0)));
    SetOutputFieldValue("DISP_Y", log10(solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_FEM(1)));
    if (nVar_FEM == 3){
      SetOutputFieldValue("DISP_Z", log10(solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetRes_FEM(2)));    
    }
  }
  
  SetOutputFieldValue("VMS", solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetTotal_CFEA());
  SetOutputFieldValue("LOAD_INCREMENT", solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetLoad_Increment());
  SetOutputFieldValue("LOAD_RAMP", solver_container[val_iZone][INST_0][MESH_0][FEA_SOL]->GetForceCoeff());
  
}

void CFEAOutput::SetOutputFields(CConfig *config){
  
  // Iteration numbers
  AddOutputField("INT_ITER",   "Int_Iter",  FORMAT_INTEGER, "INT_ITER");
  AddOutputField("EXT_ITER",   "Ext_Iter",  FORMAT_INTEGER, "EXT_ITER");
  
  // Misc.
  AddOutputField("PHYS_TIME",   "Time(min)", FORMAT_SCIENTIFIC, "PHYS_TIME");
  AddOutputField("LINSOL_ITER", "Linear_Solver_Iterations", FORMAT_INTEGER, "LINSOL_ITER");
  
  // Residuals
  AddOutputField("UTOL",   "Res_FEM[0]", FORMAT_FIXED,  "RESIDUALS");
  AddOutputField("RTOL",   "Res_FEM[1]", FORMAT_FIXED,  "RESIDUALS");
  AddOutputField("ETOL",   "Res_FEM[2]", FORMAT_FIXED,  "RESIDUALS");
  AddOutputField("DISP_X", "Res_FEM[0]", FORMAT_FIXED,  "RESIDUALS");
  AddOutputField("DISP_Y", "Res_FEM[1]", FORMAT_FIXED,  "RESIDUALS");
  AddOutputField("DISP_Z", "Res_FEM[2]", FORMAT_FIXED,  "RESIDUALS");
  
  
  AddOutputField("VMS",            "VonMises_Stress", FORMAT_FIXED, "VMS");
  AddOutputField("LOAD_INCREMENT", "Load_Increment",  FORMAT_FIXED, "LOAD_INCREMENT");
  AddOutputField("LOAD_RAMP",      "Load_Ramp",       FORMAT_FIXED, "LOAD_RAMP");
  
}

inline bool CFEAOutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) { return true;}

inline bool CFEAOutput::WriteScreen_Header(CConfig *config) {  
  
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.

  bool write_header;
  if (nonlinear_analysis) write_header = (config->GetIntIter() == 0);
  else write_header = (((config->GetExtIter() % (config->GetWrt_Con_Freq()*40)) == 0));

  return write_header;
  
}
inline bool CFEAOutput::WriteScreen_Output(CConfig *config, bool write_dualtime) {
  return true;
}


