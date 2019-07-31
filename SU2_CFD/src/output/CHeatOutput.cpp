/*!
 * \file output_heat.cpp
 * \brief Main subroutines for the heat solver output
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

#include "../../include/output/CHeatOutput.hpp"
#include "../../../Common/include/geometry_structure.hpp"
#include "../../include/solver_structure.hpp"

CHeatOutput::CHeatOutput(CConfig *config, unsigned short nDim) : COutput(config, nDim) {

  multizone = config->GetMultizone_Problem();

  /*--- Set the default history fields if nothing is set in the config file ---*/
  
  if (nRequestedHistoryFields == 0){
    RequestedHistoryFields.push_back("ITER");
    RequestedHistoryFields.push_back("RMS_RES");
    nRequestedHistoryFields = RequestedHistoryFields.size();
  }
  if (nRequestedScreenFields == 0){
    RequestedScreenFields.push_back("OUTER_ITER");    
    RequestedScreenFields.push_back("INNER_ITER");
    RequestedScreenFields.push_back("RMS_TEMPERATURE");
    nRequestedScreenFields = RequestedScreenFields.size();
  }
  if (nRequestedVolumeFields == 0){
    RequestedVolumeFields.push_back("COORDINATES");
    RequestedVolumeFields.push_back("SOLUTION");
    nRequestedVolumeFields = RequestedVolumeFields.size();
  }
  
  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Solid Heat)";
  MultiZoneHeaderString = ss.str();
  
  /*--- Set the volume filename --- */
  
  VolumeFilename = config->GetVolume_FileName();
  
  /*--- Set the surface filename --- */
  
  SurfaceFilename = config->GetSurfCoeff_FileName();
  
  /*--- Set the restart filename --- */
  
  RestartFilename = config->GetRestart_FileName();

  /*--- Set the default convergence field --- */

  if (Conv_Field.size() == 0 ) Conv_Field = "RMS_TEMPERATURE";


}

CHeatOutput::~CHeatOutput(void) {

  if (rank == MASTER_NODE){
    HistFile.close();
  }

}

void CHeatOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {
  
  CSolver* heat_solver = solver[HEAT_SOL];  

  SetHistoryOutputValue("HEATFLUX",     heat_solver->GetTotal_HeatFlux());
  SetHistoryOutputValue("HEATFLUX_MAX", heat_solver->GetTotal_MaxHeatFlux());
  SetHistoryOutputValue("AVG_TEMPERATURE",  heat_solver->GetTotal_AvgTemperature());
  SetHistoryOutputValue("RMS_TEMPERATURE", log10(heat_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("MAX_TEMPERATURE", log10(heat_solver->GetRes_Max(0)));
  if (multizone)
    SetHistoryOutputValue("BGS_TEMPERATURE", log10(heat_solver->GetRes_BGS(0)));
  
  SetHistoryOutputValue("LINSOL_ITER", heat_solver->GetIterLinSolver());
  
}
  

void CHeatOutput::SetHistoryOutputFields(CConfig *config){
  
  AddHistoryOutput("LINSOL_ITER", "Linear_Solver_Iterations", FORMAT_INTEGER, "LINSOL_ITER", "Linear solver iterations");
  
  AddHistoryOutput("RMS_TEMPERATURE", "rms[T]", FORMAT_FIXED, "RMS_RES", "Root mean square residual of the temperature", TYPE_RESIDUAL);
  AddHistoryOutput("MAX_TEMPERATURE", "max[T]", FORMAT_FIXED, "MAX_RES", "Maximum residual of the temperature", TYPE_RESIDUAL);
  AddHistoryOutput("BGS_TEMPERATURE", "bgs[T]", FORMAT_FIXED, "BGS_RES", "Block-Gauss seidel residual of the temperature", TYPE_RESIDUAL);
  
  AddHistoryOutput("HEATFLUX", "HF",      FORMAT_SCIENTIFIC, "HEAT", "Total heatflux on all surfaces defined in MARKER_MONITORING", TYPE_COEFFICIENT);
  AddHistoryOutput("HEATFLUX_MAX", "MaxHF",    FORMAT_SCIENTIFIC, "HEAT", "Total maximal heatflux on all surfaces defined in MARKER_MONITORING", TYPE_COEFFICIENT);
  AddHistoryOutput("AVG_TEMPERATURE", "AvgTemp", FORMAT_SCIENTIFIC, "HEAT", "Total average temperature on all surfaces defined in MARKER_MONITORING", TYPE_COEFFICIENT);
  
}


void CHeatOutput::SetVolumeOutputFields(CConfig *config){
  
  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES","z-component of the coordinate vector");
  
  // SOLUTION
  AddVolumeOutput("TEMPERATURE", "Temperature", "SOLUTION", "Temperature");

  // Residuals  
  AddVolumeOutput("RES_TEMPERATURE", "Residual_Temperature", "RESIDUAL", "RMS residual of the temperature");
  
}


void CHeatOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){
  
  CVariable* Node_Heat = solver[HEAT_SOL]->node[iPoint]; 
  
  CPoint*    Node_Geo  = geometry->node[iPoint];
  
  // Grid coordinates
  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));  
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));
 
  // SOLUTION
  SetVolumeOutputValue("TEMPERATURE", iPoint, Node_Heat->GetSolution(0));
  
  // Residuals    
  SetVolumeOutputValue("RES_TEMPERATURE", iPoint, solver[HEAT_SOL]->LinSysRes.GetBlock(iPoint, 0));
  
}

