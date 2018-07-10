/*!
 * \file output_adjoint_mean.cpp
 * \brief Main subroutines for elasticity discrete adjoint output
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

CDiscAdjFEAOutput::CDiscAdjFEAOutput(CConfig *config, CGeometry *geometry, unsigned short val_iZone) : COutput(config) {
 
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  
  /*--- Initialize number of variables ---*/
  if (linear_analysis) nVar_FEM = nDim;
  if (nonlinear_analysis) nVar_FEM = 3;
  
  nDim = geometry->GetnDim();

}

CDiscAdjFEAOutput::~CDiscAdjFEAOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }


}

inline bool CDiscAdjFEAOutput::WriteHistoryFile_Output(CConfig *config, bool write_dualtime) { return true; }

inline bool CDiscAdjFEAOutput::WriteScreen_Header(CConfig *config) { return true; }

inline bool CDiscAdjFEAOutput::WriteScreen_Output(CConfig *config, bool write_dualtime) { return true; }

void CDiscAdjFEAOutput::SetOutputFields(CConfig *config){
  
  // Iteration numbers
  AddOutputField("INT_ITER",   "Int_Iter",  FORMAT_INTEGER, "INT_ITER");
  AddOutputField("EXT_ITER",   "Ext_Iter",  FORMAT_INTEGER, "EXT_ITER");
  
  // Residuals
  AddOutputField("ADJOINT_DISP_X", "Res[Ux_adj]", FORMAT_FIXED,   "RESIDUALS");
  AddOutputField("ADJOINT_DISP_Y", "Res[Uy_adj]", FORMAT_FIXED,   "RESIDUALS");
  AddOutputField("ADJOINT_DISP_Z", "Res[Uz_adj]", FORMAT_FIXED,   "RESIDUALS");
  
  //Sensitivities
  AddOutputField("SENS_E", "Sens[E]",  FORMAT_FIXED, "SENSITIVITY");
  AddOutputField("SENS_NU","Sens[Nu]", FORMAT_FIXED, "SENSITIVITY");

  
}

inline void CDiscAdjFEAOutput::LoadHistoryData(CGeometry ****geometry, CSolver *****solver_container, CConfig **config,
      CIntegration ****integration, bool DualTime, su2double timeused, unsigned short val_iZone, unsigned short val_iInst) {
  
  SetOutputFieldValue("INT_ITER", config[val_iZone]->GetIntIter());
  SetOutputFieldValue("EXT_ITER", config[val_iZone]->GetExtIter());
  
  SetOutputFieldValue("PHYS_TIME", timeused);
  
  SetOutputFieldValue("ADJOINT_DISP_X", log10(solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetRes_RMS(0)));
  SetOutputFieldValue("ADJOINT_DISP_Y", log10(solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetRes_RMS(1)));
  if (nVar_FEM == 3){
    SetOutputFieldValue("ADJOINT_DISP_Z", log10(solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetRes_RMS(2)));    
  }
  su2double Total_SensE = 0.0; su2double Total_SensNu = 0.0;  
  if (config[val_iZone]->GetnElasticityMod() == 1){
    Total_SensE = solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_E(0);
    Total_SensNu = solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_Nu(0);
  }
  else{
    for (unsigned short iVar = 0; iVar < config[val_iZone]->GetnElasticityMod(); iVar++){
        Total_SensE += solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_E(0)
            *solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_E(0);
        Total_SensNu += solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_Nu(0)
            *solver_container[val_iZone][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_Nu(0);
    }
    Total_SensE = sqrt(Total_SensE);
    Total_SensNu = sqrt(Total_SensNu);

  }
  SetOutputFieldValue("SENS_E", Total_SensE);
  SetOutputFieldValue("SENS_NU", Total_SensNu);
  
}

