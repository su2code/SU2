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

void CFEAOutput::LoadHistoryData(CGeometry ****geometry,
                                     CSolver *****solver_container,
                                     CConfig **config,
                                     CIntegration ****integration,
                                     bool DualTime_Iteration,
                                     su2double timeused,
                                     unsigned short val_iZone,
                                     unsigned short val_iInst) {

  CSolver* fea_solver = solver_container[val_iZone][val_iInst][MESH_0][FEA_SOL];
  
  
  bool linear_analysis = (config[val_iZone]->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config[val_iZone]->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.

  SetHistoryOutputValue("INT_ITER", config[val_iZone]->GetIntIter());
  SetHistoryOutputValue("EXT_ITER", config[val_iZone]->GetExtIter());
  
  SetHistoryOutputValue("PHYS_TIME", timeused);
  
  /*--- Residuals: ---*/
  /*--- Linear analysis: RMS of the displacements in the nDim coordinates ---*/
  /*--- Nonlinear analysis: UTOL, RTOL and DTOL (defined in the Postprocessing function) ---*/

  
  if (linear_analysis){
    SetHistoryOutputValue("UTOL", log10(fea_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("RTOL", log10(fea_solver->GetRes_RMS(1)));
    if (nVar_FEM == 3){
      SetHistoryOutputValue("ETOL", log10(fea_solver->GetRes_RMS(2)));    
    }
    SetHistoryOutputValue("DISP_X", log10(fea_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("DISP_Y", log10(fea_solver->GetRes_RMS(1)));
    if (nVar_FEM == 3){
      SetHistoryOutputValue("DISP_Z", log10(fea_solver->GetRes_RMS(2)));    
    }
  } else if (nonlinear_analysis){
    SetHistoryOutputValue("UTOL", log10(fea_solver->GetRes_FEM(0)));
    SetHistoryOutputValue("RTOL", log10(fea_solver->GetRes_FEM(1)));
    if (nVar_FEM == 3){
      SetHistoryOutputValue("ETOL", log10(fea_solver->GetRes_FEM(2)));    
    }
    SetHistoryOutputValue("DISP_X", log10(fea_solver->GetRes_FEM(0)));
    SetHistoryOutputValue("DISP_Y", log10(fea_solver->GetRes_FEM(1)));
    if (nVar_FEM == 3){
      SetHistoryOutputValue("DISP_Z", log10(fea_solver->GetRes_FEM(2)));    
    }
  }
  
  SetHistoryOutputValue("VMS", fea_solver->GetTotal_CFEA());
  SetHistoryOutputValue("LOAD_INCREMENT", fea_solver->GetLoad_Increment());
  SetHistoryOutputValue("LOAD_RAMP", fea_solver->GetForceCoeff());
  
}

void CFEAOutput::SetHistoryOutputFields(CConfig *config){
  
  // Iteration numbers
  AddHistoryOutput("INT_ITER",   "Int_Iter",  FORMAT_INTEGER, "INT_ITER");
  AddHistoryOutput("EXT_ITER",   "Ext_Iter",  FORMAT_INTEGER, "EXT_ITER");
  
  // Misc.
  AddHistoryOutput("PHYS_TIME",   "Time(min)", FORMAT_SCIENTIFIC, "PHYS_TIME");
  AddHistoryOutput("LINSOL_ITER", "Linear_Solver_Iterations", FORMAT_INTEGER, "LINSOL_ITER");
  
  // Residuals
  AddHistoryOutput("UTOL",   "Res[U]", FORMAT_FIXED,  "RESIDUALS");
  AddHistoryOutput("RTOL",   "Res[R]", FORMAT_FIXED,  "RESIDUALS");
  AddHistoryOutput("ETOL",   "Res[E]", FORMAT_FIXED,  "RESIDUALS");
  AddHistoryOutput("DISP_X", "Res[DispX]", FORMAT_FIXED,  "RESIDUALS");
  AddHistoryOutput("DISP_Y", "Res[DispY]", FORMAT_FIXED,  "RESIDUALS");
  AddHistoryOutput("DISP_Z", "Res[DispZ]", FORMAT_FIXED,  "RESIDUALS");
  
  
  AddHistoryOutput("VMS",            "VonMises_Stress", FORMAT_FIXED, "VMS");
  AddHistoryOutput("LOAD_INCREMENT", "Load_Increment",  FORMAT_FIXED, "LOAD_INCREMENT");
  AddHistoryOutput("LOAD_RAMP",      "Load_Ramp",       FORMAT_FIXED, "LOAD_RAMP");
  
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


