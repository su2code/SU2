/*!
 * \file output_elasticity.cpp
 * \brief Main subroutines for FEA output
 * \author R. Sanchez
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/output/CElasticityOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solver_structure.hpp"

CElasticityOutput::CElasticityOutput(CConfig *config, unsigned short nDim) : COutput(config, nDim, false) {

  linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
  nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
  dynamic = (config->GetTime_Domain());  // Dynamic analysis.

  /*--- Initialize number of variables ---*/
  if (linear_analysis) nVar_FEM = nDim;
  if (nonlinear_analysis) nVar_FEM = 3;

  /*--- Default fields for screen output ---*/
  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.emplace_back("ITER");
    requestedHistoryFields.emplace_back("RMS_RES");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }

  /*--- Default fields for screen output ---*/
  if (nRequestedScreenFields == 0){
    if (dynamic) requestedScreenFields.emplace_back("TIME_ITER");
    if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    if(linear_analysis){
      requestedScreenFields.emplace_back("RMS_DISP_X");
      requestedScreenFields.emplace_back("RMS_DISP_Y");
      requestedScreenFields.emplace_back("RMS_DISP_Z");
    }
    if(nonlinear_analysis){
      requestedScreenFields.emplace_back("RMS_UTOL");
      requestedScreenFields.emplace_back("RMS_RTOL");
      requestedScreenFields.emplace_back("RMS_ETOL");
    }
    requestedScreenFields.emplace_back("VMS");
    nRequestedScreenFields = requestedScreenFields.size();
  }

  /*--- Default fields for volume output ---*/
  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    requestedVolumeFields.emplace_back("STRESS");
    nRequestedVolumeFields = requestedVolumeFields.size();
  }

  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Structure)";
  multiZoneHeaderString = ss.str();

  /*--- Set the volume filename --- */

  volumeFilename = config->GetVolume_FileName();

  /*--- Set the surface filename --- */

  surfaceFilename = config->GetSurfCoeff_FileName();

  /*--- Set the restart filename --- */

  restartFilename = config->GetRestart_FileName();

  /*--- Set the default convergence field --- */

  if (convFields.empty() ) convFields.emplace_back("RMS_DISP_X");

}

CElasticityOutput::~CElasticityOutput(void) {}

void CElasticityOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver)  {

  CSolver* fea_solver = solver[FEA_SOL];

  /*--- Residuals: ---*/
  /*--- Linear analysis: RMS of the displacements in the nDim coordinates ---*/
  /*--- Nonlinear analysis: UTOL, RTOL and DTOL (defined in the Postprocessing function) ---*/


  if (linear_analysis){
    SetHistoryOutputValue("RMS_DISP_X", log10(fea_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("RMS_DISP_Y", log10(fea_solver->GetRes_RMS(1)));
    if (nDim == 3){
      SetHistoryOutputValue("RMS_DISP_Z", log10(fea_solver->GetRes_RMS(2)));
    }
  } else if (nonlinear_analysis){
    SetHistoryOutputValue("RMS_UTOL", log10(fea_solver->LinSysSol.norm()));
    SetHistoryOutputValue("RMS_RTOL", log10(fea_solver->LinSysRes.norm()));
    SetHistoryOutputValue("RMS_ETOL", log10(fea_solver->LinSysSol.dot(fea_solver->LinSysRes)));

  }

  if (multiZone){
    SetHistoryOutputValue("BGS_DISP_X", log10(fea_solver->GetRes_BGS(0)));
    SetHistoryOutputValue("BGS_DISP_Y", log10(fea_solver->GetRes_BGS(1)));
    if (nDim == 3) SetHistoryOutputValue("BGS_DISP_Z", log10(fea_solver->GetRes_BGS(2)));
  }

  SetHistoryOutputValue("VMS", fea_solver->GetTotal_CFEA());
  SetHistoryOutputValue("LOAD_INCREMENT", fea_solver->GetLoad_Increment());
  SetHistoryOutputValue("LOAD_RAMP", fea_solver->GetForceCoeff());

  SetHistoryOutputValue("LINSOL_ITER", fea_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(fea_solver->GetResLinSolver()));
  
} 

void CElasticityOutput::SetHistoryOutputFields(CConfig *config){

  AddHistoryOutput("LINSOL_ITER", "LinSolIter", ScreenOutputFormat::INTEGER, "LINSOL",  "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL",  "Residual of the linear solver.");

  // Residuals

  AddHistoryOutput("RMS_UTOL",   "rms[U]", ScreenOutputFormat::FIXED,  "RMS_RES", "", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("RMS_RTOL",   "rms[R]", ScreenOutputFormat::FIXED,  "RMS_RES", "", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("RMS_ETOL",   "rms[E]", ScreenOutputFormat::FIXED,  "RMS_RES", "", HistoryFieldType::RESIDUAL);

  AddHistoryOutput("RMS_DISP_X", "rms[DispX]", ScreenOutputFormat::FIXED,  "RMS_RES", "", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("RMS_DISP_Y", "rms[DispY]", ScreenOutputFormat::FIXED,  "RMS_RES", "", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("RMS_DISP_Z", "rms[DispZ]", ScreenOutputFormat::FIXED,  "RMS_RES", "", HistoryFieldType::RESIDUAL);

  AddHistoryOutput("BGS_DISP_X", "bgs[DispX]", ScreenOutputFormat::FIXED,  "BGS_RES", "", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("BGS_DISP_Y", "bgs[DispY]", ScreenOutputFormat::FIXED,  "BGS_RES", "", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("BGS_DISP_Z", "bgs[DispZ]", ScreenOutputFormat::FIXED,  "BGS_RES", "", HistoryFieldType::RESIDUAL);

  AddHistoryOutput("VMS",            "VonMises", ScreenOutputFormat::SCIENTIFIC, "", "VMS");
  AddHistoryOutput("LOAD_INCREMENT", "Load_Increment",  ScreenOutputFormat::FIXED, "", "LOAD_INCREMENT");
  AddHistoryOutput("LOAD_RAMP",      "Load_Ramp",       ScreenOutputFormat::FIXED, "", "LOAD_RAMP");

}

void CElasticityOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_Struc = solver[FEA_SOL]->GetNodes();
  CPoint*    Node_Geo  = geometry->node[iPoint];

  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));

  SetVolumeOutputValue("DISPLACEMENT-X", iPoint, Node_Struc->GetSolution(iPoint, 0));
  SetVolumeOutputValue("DISPLACEMENT-Y", iPoint, Node_Struc->GetSolution(iPoint, 1));
  if (nDim == 3) SetVolumeOutputValue("DISPLACEMENT-Z", iPoint, Node_Struc->GetSolution(iPoint, 2));

  if(dynamic){
    SetVolumeOutputValue("VELOCITY-X", iPoint, Node_Struc->GetSolution_Vel(iPoint, 0));
    SetVolumeOutputValue("VELOCITY-Y", iPoint, Node_Struc->GetSolution_Vel(iPoint, 1));
    if (nDim == 3) SetVolumeOutputValue("VELOCITY-Z", iPoint, Node_Struc->GetSolution_Vel(iPoint, 2));

    SetVolumeOutputValue("ACCELERATION-X", iPoint, Node_Struc->GetSolution_Accel(iPoint, 0));
    SetVolumeOutputValue("ACCELERATION-Y", iPoint, Node_Struc->GetSolution_Accel(iPoint, 1));
    if (nDim == 3) SetVolumeOutputValue("ACCELERATION-Z", iPoint, Node_Struc->GetSolution_Accel(iPoint, 2));
  }

  SetVolumeOutputValue("STRESS-XX", iPoint, Node_Struc->GetStress_FEM(iPoint)[0]);
  SetVolumeOutputValue("STRESS-YY", iPoint, Node_Struc->GetStress_FEM(iPoint)[1]);
  SetVolumeOutputValue("STRESS-XY", iPoint, Node_Struc->GetStress_FEM(iPoint)[2]);
  if (nDim == 3){
    SetVolumeOutputValue("STRESS-ZZ", iPoint, Node_Struc->GetStress_FEM(iPoint)[3]);
    SetVolumeOutputValue("STRESS-XZ", iPoint, Node_Struc->GetStress_FEM(iPoint)[4]);
    SetVolumeOutputValue("STRESS-YZ", iPoint, Node_Struc->GetStress_FEM(iPoint)[5]);
  }
  SetVolumeOutputValue("VON_MISES_STRESS", iPoint, Node_Struc->GetVonMises_Stress(iPoint));

}

void CElasticityOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");

  AddVolumeOutput("DISPLACEMENT-X",    "Displacement_x", "SOLUTION", "x-component of the displacement vector");
  AddVolumeOutput("DISPLACEMENT-Y",    "Displacement_y", "SOLUTION", "y-component of the displacement vector");
  if (nDim == 3) AddVolumeOutput("DISPLACEMENT-Z", "Displacement_z", "SOLUTION", "z-component of the displacement vector");

  if(dynamic){
    AddVolumeOutput("VELOCITY-X",    "Velocity_x", "VELOCITY", "x-component of the velocity vector");
    AddVolumeOutput("VELOCITY-Y",    "Velocity_y", "VELOCITY", "y-component of the velocity vector");
    if (nDim == 3) AddVolumeOutput("VELOCITY-Z", "Velocity_z", "VELOCITY", "z-component of the velocity vector");

    AddVolumeOutput("ACCELERATION-X",    "Acceleration_x", "ACCELERATION", "x-component of the acceleration vector");
    AddVolumeOutput("ACCELERATION-Y",    "Acceleration_y", "ACCELERATION", "y-component of the acceleration vector");
    if (nDim == 3) AddVolumeOutput("ACCELERATION-Z", "Acceleration_z", "ACCELERATION", "z-component of the acceleration vector");
  }

  AddVolumeOutput("STRESS-XX",    "Sxx", "STRESS", "x-component of the normal stress vector");
  AddVolumeOutput("STRESS-YY",    "Syy", "STRESS", "y-component of the normal stress vector");
  AddVolumeOutput("STRESS-XY",    "Sxy", "STRESS", "xy shear stress component");

  if (nDim == 3) {
    AddVolumeOutput("STRESS-ZZ",    "Szz", "STRESS", "z-component of the normal stress vector");
    AddVolumeOutput("STRESS-XZ",    "Sxz", "STRESS", "xz shear stress component");
    AddVolumeOutput("STRESS-YZ",    "Syz", "STRESS", "yz shear stress component");
  }

  AddVolumeOutput("VON_MISES_STRESS", "Von_Mises_Stress", "STRESS", "von-Mises stress");

}
bool CElasticityOutput::SetInit_Residuals(CConfig *config){

  return (config->GetTime_Domain() == NO && (curInnerIter  == 0));

}


