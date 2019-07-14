/*!
 * \file output_adj_elasticity.cpp
 * \brief Main subroutines for elasticity discrete adjoint output
 * \author R. Sanchez
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

#include "../../include/output/CAdjElasticityOutput.hpp"

#include "../../../Common/include/geometry_structure.hpp"
#include "../../include/solver_structure.hpp"

CAdjElasticityOutput::CAdjElasticityOutput(CConfig *config, CGeometry *geometry, unsigned short val_iZone) : COutput(config) {
 
  bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  
  /*--- Initialize number of variables ---*/
  if (linear_analysis) nVar_FEM = nDim;
  if (nonlinear_analysis) nVar_FEM = 3;
  
  nDim = geometry->GetnDim();

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    RequestedHistoryFields.push_back("ITER");
    RequestedHistoryFields.push_back("RESIDUALS");
    RequestedHistoryFields.push_back("SENSITIVITY");
    nRequestedHistoryFields = RequestedHistoryFields.size();
  }

  if (nRequestedScreenFields == 0){
    if (multizone) RequestedScreenFields.push_back("OUTER_ITER");
    RequestedScreenFields.push_back("INNER_ITER");
    RequestedScreenFields.push_back("ADJOINT_DISP_X");
    RequestedScreenFields.push_back("ADJOINT_DISP_Y");
    RequestedScreenFields.push_back("SENS_E");
    RequestedScreenFields.push_back("SENS_NU");
    nRequestedScreenFields = RequestedScreenFields.size();
  }

  if (nRequestedVolumeFields == 0){
    RequestedVolumeFields.push_back("COORDINATES");
    RequestedVolumeFields.push_back("SOLUTION");
    nRequestedVolumeFields = RequestedVolumeFields.size();
  }

  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Adj. Comp. Fluid)";
  MultiZoneHeaderString = ss.str();

  /*--- Set the volume filename --- */
  
  VolumeFilename = config->GetAdj_FileName();
  
  /*--- Set the surface filename --- */
  
  SurfaceFilename = config->GetSurfAdjCoeff_FileName();
  
  /*--- Set the restart filename --- */
  
  RestartFilename = config->GetRestart_FileName();

  /*--- Set the default convergence field --- */

  if (Conv_Field.size() == 0 ) Conv_Field = "ADJOINT_DISP_X";
  
}

CAdjElasticityOutput::~CAdjElasticityOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }


}

void CAdjElasticityOutput::SetHistoryOutputFields(CConfig *config){
  
  // Residuals
  AddHistoryOutput("ADJOINT_DISP_X", "Res[Ux_adj]", FORMAT_FIXED,   "RESIDUALS");
  AddHistoryOutput("ADJOINT_DISP_Y", "Res[Uy_adj]", FORMAT_FIXED,   "RESIDUALS");
  AddHistoryOutput("ADJOINT_DISP_Z", "Res[Uz_adj]", FORMAT_FIXED,   "RESIDUALS");
  
  //Sensitivities
  AddHistoryOutput("SENS_E", "Sens[E]",  FORMAT_FIXED, "SENSITIVITY");
  AddHistoryOutput("SENS_NU","Sens[Nu]", FORMAT_FIXED, "SENSITIVITY");

  
}

inline void CAdjElasticityOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {
  
  SetHistoryOutputValue("ADJOINT_DISP_X", log10(solver[ADJFEA_SOL]->GetRes_RMS(0)));
  SetHistoryOutputValue("ADJOINT_DISP_Y", log10(solver[ADJFEA_SOL]->GetRes_RMS(1)));
  if (nVar_FEM == 3){
    SetHistoryOutputValue("ADJOINT_DISP_Z", log10(solver[ADJFEA_SOL]->GetRes_RMS(2)));    
  }
  su2double Total_SensE = 0.0; su2double Total_SensNu = 0.0;  
  if (config->GetnElasticityMod() == 1){
    Total_SensE = solver[ADJFEA_SOL]->GetGlobal_Sens_E(0);
    Total_SensNu = solver[ADJFEA_SOL]->GetGlobal_Sens_Nu(0);
  }
  else{
    // TODO: Update this and change tests
    for (unsigned short iVar = 0; iVar < config->GetnElasticityMod(); iVar++){
      Total_SensE += solver[ADJFEA_SOL]->GetGlobal_Sens_E(0)
                    *solver[ADJFEA_SOL]->GetGlobal_Sens_E(0);
      Total_SensNu += solver[ADJFEA_SOL]->GetGlobal_Sens_Nu(0)
                     *solver[ADJFEA_SOL]->GetGlobal_Sens_Nu(0);
    }
  Total_SensE = sqrt(Total_SensE);
  Total_SensNu = sqrt(Total_SensNu);
  }
  SetHistoryOutputValue("SENS_E", Total_SensE);
  SetHistoryOutputValue("SENS_NU", Total_SensNu);
  
}

void CAdjElasticityOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_Struc = solver[FEA_SOL]->node[iPoint];
  CPoint*    Node_Geo  = geometry->node[iPoint];

  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));

  SetVolumeOutputValue("SENSITIVITY-X", iPoint, Node_Struc->GetSolution(0));
  SetVolumeOutputValue("SENSITIVITY-Y", iPoint, Node_Struc->GetSolution(1));
  if (nDim == 3) SetVolumeOutputValue("SENSITIVITY-Z", iPoint, Node_Struc->GetSolution(2));

}

void CAdjElasticityOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES");

  AddVolumeOutput("SENSITIVITY-X",    "Sensitivity_x", "SOLUTION");
  AddVolumeOutput("SENSITIVITY-Y",    "Sensitivity_y", "SOLUTION");
  if (nDim == 3) AddVolumeOutput("SENSITIVITY-Z", "Sensitivity_z", "SOLUTION");

}
