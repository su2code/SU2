/*!
 * \file output_adj_heat.cpp
 * \brief Main subroutines for flow discrete adjoint output
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

#include "../../include/output/CAdjHeatOutput.hpp"

#include "../../../Common/include/geometry_structure.hpp"
#include "../../include/solver_structure.hpp"

CAdjHeatOutput::CAdjHeatOutput(CConfig *config, CGeometry *geometry, unsigned short val_iZone) : COutput(config) {
  
  nDim = geometry->GetnDim();

  /*--- Set the default history fields if nothing is set in the config file ---*/
  
  if (nRequestedHistoryFields == 0){
    RequestedHistoryFields.push_back("ITER");
    RequestedHistoryFields.push_back("RMS_RES");
    RequestedHistoryFields.push_back("SENSITIVITY");
    nRequestedHistoryFields = RequestedHistoryFields.size();
  }
  
  if (nRequestedScreenFields == 0){
    if (multizone) RequestedScreenFields.push_back("OUTER_ITER");
    RequestedScreenFields.push_back("INNER_ITER");    
    RequestedScreenFields.push_back("RMS_ADJ_TEMPERATURE");
    RequestedScreenFields.push_back("SENS_GEO");
    nRequestedScreenFields = RequestedScreenFields.size();
  }
  
  if (nRequestedVolumeFields == 0){
    RequestedVolumeFields.push_back("COORDINATES");
    RequestedVolumeFields.push_back("SOLUTION");
    RequestedVolumeFields.push_back("SENSITIVITY");
    nRequestedVolumeFields = RequestedVolumeFields.size();
  }
  
  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Adj. Heat)";
  MultiZoneHeaderString = ss.str();
  
  /*--- Set the volume filename --- */
  
  VolumeFilename = config->GetAdj_FileName();
  
  /*--- Set the surface filename --- */
  
  SurfaceFilename = config->GetSurfAdjCoeff_FileName();
  
  /*--- Set the restart filename --- */
  
  RestartFilename = config->GetRestart_AdjFileName();
  
  /*--- Add the obj. function extension --- */
  
  RestartFilename = config->GetObjFunc_Extension(RestartFilename);

  /*--- Set the default convergence field --- */

  if (Conv_Field.size() == 0 ) Conv_Field = "RMS_ADJ_TEMPERATURE";

}

CAdjHeatOutput::~CAdjHeatOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }

}

void CAdjHeatOutput::SetHistoryOutputFields(CConfig *config){
  
  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the conservative variables. 
  /// DESCRIPTION: Root-mean square residual of the adjoint Pressure.
  AddHistoryOutput("RMS_ADJ_TEMPERATURE",    "rms[A_T]",  FORMAT_FIXED, "RMS_RES", TYPE_RESIDUAL); 
  /// END_GROUP
  
  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the conservative variables. 
  /// DESCRIPTION: Maximum residual of the adjoint Pressure.
  AddHistoryOutput("MAX_ADJ_TEMPERATURE",    "max[A_T]",  FORMAT_FIXED, "MAX_RES", TYPE_RESIDUAL);

  
  /// BEGIN_GROUP: SENSITIVITY, DESCRIPTION: Sensitivities of different geometrical or boundary values.   
  /// DESCRIPTION: Sum of the geometrical sensitivities on all markers set in MARKER_MONITORING.
  AddHistoryOutput("SENS_GEO",   "Sens_Geo",   FORMAT_SCIENTIFIC, "SENSITIVITY", TYPE_COEFFICIENT); 
  /// END_GROUP
  
}

void CAdjHeatOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) { 
  
  CSolver* adjheat_solver = solver[ADJHEAT_SOL];

  SetHistoryOutputValue("RMS_ADJ_TEMPERATURE", log10(adjheat_solver->GetRes_RMS(0)));
 
  SetHistoryOutputValue("MAX_ADJ_TEMPERATURE", log10(adjheat_solver->GetRes_Max(0)));
 
  if (multizone)
    SetHistoryOutputValue("MAX_ADJ_TEMPERATURE", log10(adjheat_solver->GetRes_BGS(0)));
    
  SetHistoryOutputValue("SENS_GEO", adjheat_solver->GetTotal_Sens_Geo());
 

}

void CAdjHeatOutput::SetVolumeOutputFields(CConfig *config){
  
  /// BEGIN_GROUP: COORDINATES, DESCRIPTION: Coordinates of the mesh nodes.
  /// DESCRIPTION: x coordinates of the mesh nodes.
  AddVolumeOutput("COORD-X", "x", "COORDINATES"); 
  /// DESCRIPTION: y coordinates of the mesh nodes.
  AddVolumeOutput("COORD-Y", "y", "COORDINATES");
  if (nDim == 3)
    /// DESCRIPTION: z coordinates of the mesh nodes.
    AddVolumeOutput("COORD-Z", "z", "COORDINATES");
  /// END_GROUP
  
  /// BEGIN_GROUP: CONSERVATIVE, DESCRIPTION: The conservative variables of the adjoint solver.
  /// DESCRIPTION: Adjoint Pressure.
  AddVolumeOutput("ADJ_TEMPERATURE",    "Adjoint_Temperature",    "SOLUTION");
  /// END_GROUP
  
  
  /// BEGIN_GROUP: RESIDUAL, DESCRIPTION: Residuals of the conservative variables. 
  /// DESCRIPTION: Residual of the adjoint Pressure.
  AddVolumeOutput("RES_ADJ_TEMPERATURE",    "Residual_Adjoint_Temperature",    "RESIDUAL");
  /// END_GROUP
  
  /// BEGIN_GROUP: SENSITIVITY, DESCRIPTION: Geometrical sensitivities of the current objective function.
  /// DESCRIPTION: Sensitivity x-component.
  AddVolumeOutput("SENSITIVITY_X", "Sensitivity_x", "SENSITIVITY"); 
  /// DESCRIPTION: Sensitivity y-component.
  AddVolumeOutput("SENSITIVITY_Y", "Sensitivity_y", "SENSITIVITY");
  if (nDim == 3)
    /// DESCRIPTION: Sensitivity z-component.
    AddVolumeOutput("SENSITIVITY_Z", "Sensitivity_z", "SENSITIVITY");   
  /// DESCRIPTION: Sensitivity in normal direction.
  AddVolumeOutput("SENSITIVITY", "Surface_Sensitivity", "SENSITIVITY"); 
  /// END_GROUP
 
}

void CAdjHeatOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){
  
  CVariable* Node_AdjHeat = solver[ADJHEAT_SOL]->node[iPoint]; 
  CPoint*    Node_Geo     = geometry->node[iPoint];
  

  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));  
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));
  
  SetVolumeOutputValue("ADJ_TEMPERATURE",    iPoint, Node_AdjHeat->GetSolution(0));
  
  // Residuals
  SetVolumeOutputValue("RES_ADJ_TEMPERATURE", iPoint, Node_AdjHeat->GetSolution(0) - Node_AdjHeat->GetSolution_Old(0));
  
  SetVolumeOutputValue("SENSITIVITY_X", iPoint, Node_AdjHeat->GetSensitivity(0));
  SetVolumeOutputValue("SENSITIVITY_Y", iPoint, Node_AdjHeat->GetSensitivity(1));
  if (nDim == 3)
    SetVolumeOutputValue("SENSITIVITY_Z", iPoint, Node_AdjHeat->GetSensitivity(2));
  
}

void CAdjHeatOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){
  
  SetVolumeOutputValue("SENSITIVITY", iPoint, solver[ADJHEAT_SOL]->GetCSensitivity(iMarker, iVertex));
  
}


