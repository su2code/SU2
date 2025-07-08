/*!
 * \file CElasticityOutput.cpp
 * \brief Main subroutines for FEA output
 * \author R. Sanchez
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../include/output/CHeatOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CElasticityOutput::CElasticityOutput(CConfig *config, unsigned short nDim) : COutput(config, nDim, false) {

  linear_analysis = config->GetGeometricConditions() == STRUCT_DEFORMATION::SMALL;
  nonlinear_analysis = config->GetGeometricConditions() == STRUCT_DEFORMATION::LARGE;
  coupled_heat = config->GetWeakly_Coupled_Heat();
  dynamic = config->GetTime_Domain();

  /*--- Initialize number of variables ---*/
  if (linear_analysis) nVar_FEM = nDim;
  if (nonlinear_analysis) nVar_FEM = 3;

  /*--- Default fields for screen output ---*/
  if (nRequestedHistoryFields == 0){
    RequestCommonHistory(dynamic);
    nRequestedHistoryFields = requestedHistoryFields.size();
  }

  /*--- Default fields for screen output ---*/
  if (nRequestedScreenFields == 0) {
    if (dynamic) requestedScreenFields.emplace_back("TIME_ITER");
    if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    if (linear_analysis) {
      requestedScreenFields.emplace_back("RMS_DISP_X");
      requestedScreenFields.emplace_back("RMS_DISP_Y");
      requestedScreenFields.emplace_back("RMS_DISP_Z");
    }
    if (nonlinear_analysis) {
      requestedScreenFields.emplace_back("RMS_UTOL");
      requestedScreenFields.emplace_back("RMS_RTOL");
      requestedScreenFields.emplace_back("RMS_ETOL");
    }
    if (coupled_heat) requestedScreenFields.emplace_back("RMS_TEMPERATURE");
    requestedScreenFields.emplace_back("VMS");
    nRequestedScreenFields = requestedScreenFields.size();
  }

  /*--- Default fields for volume output ---*/
  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    requestedVolumeFields.emplace_back("STRESS");
    if (config->GetTopology_Optimization()) requestedVolumeFields.emplace_back("TOPOLOGY");
    if (coupled_heat) requestedVolumeFields.emplace_back("PRIMITIVE");
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

  if (convFields.empty()) {
    if (linear_analysis) convFields.emplace_back("RMS_DISP_X");
    if (nonlinear_analysis) convFields.emplace_back("RMS_UTOL");
  }
}

void CElasticityOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver)  {

  CSolver* fea_solver = solver[FEA_SOL];
  CSolver* heat_solver = solver[HEAT_SOL];

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
    SetHistoryOutputValue("RMS_UTOL", log10(fea_solver->GetRes_FEM(0)));
    SetHistoryOutputValue("RMS_RTOL", log10(fea_solver->GetRes_FEM(1)));
    SetHistoryOutputValue("RMS_ETOL", log10(fea_solver->GetRes_FEM(2)));
  }

  if (multiZone){
    SetHistoryOutputValue("BGS_DISP_X", log10(fea_solver->GetRes_BGS(0)));
    SetHistoryOutputValue("BGS_DISP_Y", log10(fea_solver->GetRes_BGS(1)));
    if (nDim == 3) SetHistoryOutputValue("BGS_DISP_Z", log10(fea_solver->GetRes_BGS(2)));
  }

  SetHistoryOutputValue("VMS", fea_solver->GetTotal_CFEA());
  SetHistoryOutputValue("LOAD_INCREMENT", fea_solver->GetLoad_Increment()*100);
  SetHistoryOutputValue("LOAD_RAMP", fea_solver->GetForceCoeff());

  SetHistoryOutputValue("LINSOL_ITER", fea_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(fea_solver->GetResLinSolver()));

  SetHistoryOutputValue("REFERENCE_NODE", fea_solver->GetTotal_OFRefNode());
  SetHistoryOutputValue("TOPOL_COMPLIANCE", fea_solver->GetTotal_OFCompliance());
  SetHistoryOutputValue("STRESS_PENALTY", fea_solver->GetTotal_OFStressPenalty());
  if (config->GetRefGeom()) {
    SetHistoryOutputValue("REFERENCE_GEOMETRY", fea_solver->GetTotal_OFRefGeom());
  }
  if (config->GetTopology_Optimization()) {
    SetHistoryOutputValue("VOLUME_FRACTION", fea_solver->GetTotal_OFVolFrac());
    SetHistoryOutputValue("TOPOL_DISCRETENESS", fea_solver->GetTotal_OFDiscreteness());
  }

  /*--- Add heat solver data if available. ---*/
  if (coupled_heat) {
    CHeatOutput::LoadHistoryDataImpl(config, geometry, solver, this);
    SetHistoryOutputValue("LINSOL_ITER_HEAT", heat_solver->GetIterLinSolver());
    SetHistoryOutputValue("LINSOL_RESIDUAL_HEAT", log10(heat_solver->GetResLinSolver()));
  }

  ComputeSimpleCustomOutputs(config);

  /*--- Keep this as last, since it uses the history values that were set. ---*/
  SetCustomAndComboObjectives(FEA_SOL, config, solver);

}

void CElasticityOutput::SetHistoryOutputFields(CConfig *config) {

  AddHistoryOutput("LINSOL_ITER", "LinSolIter", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");

  if (nonlinear_analysis) {
    AddHistoryOutput("RMS_UTOL", "rms[U]", ScreenOutputFormat::FIXED, "RMS_RES", "Norm of displacement increment", HistoryFieldType::RESIDUAL);
    AddHistoryOutput("RMS_RTOL", "rms[R]", ScreenOutputFormat::FIXED, "RMS_RES", "Norm of residual", HistoryFieldType::RESIDUAL);
    AddHistoryOutput("RMS_ETOL", "rms[E]", ScreenOutputFormat::FIXED, "RMS_RES", "Norm of energy/work increment", HistoryFieldType::RESIDUAL);
  } else if (linear_analysis) {
    AddHistoryOutput("RMS_DISP_X", "rms[DispX]", ScreenOutputFormat::FIXED, "RMS_RES", "Residual of X displacement", HistoryFieldType::RESIDUAL);
    AddHistoryOutput("RMS_DISP_Y", "rms[DispY]", ScreenOutputFormat::FIXED, "RMS_RES", "Residual of Y displacement", HistoryFieldType::RESIDUAL);
    AddHistoryOutput("RMS_DISP_Z", "rms[DispZ]", ScreenOutputFormat::FIXED, "RMS_RES", "Residual of Z displacement", HistoryFieldType::RESIDUAL);
  }
  if (multiZone) {
    AddHistoryOutput("BGS_DISP_X", "bgs[DispX]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of X displacement", HistoryFieldType::RESIDUAL);
    AddHistoryOutput("BGS_DISP_Y", "bgs[DispY]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of Y displacement", HistoryFieldType::RESIDUAL);
    AddHistoryOutput("BGS_DISP_Z", "bgs[DispZ]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of Z displacement", HistoryFieldType::RESIDUAL);
  }
  AddHistoryOutput("VMS", "VonMises", ScreenOutputFormat::SCIENTIFIC, "Maximum Von-Misses stress", "VMS");
  AddHistoryOutput("LOAD_RAMP", "Load_Ramp", ScreenOutputFormat::FIXED, "Fraction of total load (ramped)", "LOAD_RAMP");
  AddHistoryOutput("LOAD_INCREMENT", "Load[%]", ScreenOutputFormat::PERCENT, "Percent of total load (incremental)", "LOAD_INCREMENT");

  AddHistoryOutput("REFERENCE_NODE", "RefNode", ScreenOutputFormat::SCIENTIFIC, "STRUCT_COEFF", "Distance to reference node", HistoryFieldType::COEFFICIENT);
  AddHistoryOutput("TOPOL_COMPLIANCE", "TopComp", ScreenOutputFormat::SCIENTIFIC, "STRUCT_COEFF", "Structural compliance", HistoryFieldType::COEFFICIENT);
  AddHistoryOutput("STRESS_PENALTY", "StressPen", ScreenOutputFormat::SCIENTIFIC, "STRUCT_COEFF", "Aggregate stress penalty", HistoryFieldType::COEFFICIENT);
  if (config->GetRefGeom()) {
    AddHistoryOutput("REFERENCE_GEOMETRY", "RefGeom", ScreenOutputFormat::SCIENTIFIC, "STRUCT_COEFF", "L2 norm of difference wrt reference geometry", HistoryFieldType::COEFFICIENT);
  }
  if (config->GetTopology_Optimization()) {
    AddHistoryOutput("VOLUME_FRACTION", "VolFrac", ScreenOutputFormat::SCIENTIFIC, "STRUCT_COEFF", "Fraction of solid material", HistoryFieldType::COEFFICIENT);
    AddHistoryOutput("TOPOL_DISCRETENESS", "TopDisc", ScreenOutputFormat::SCIENTIFIC, "STRUCT_COEFF", "Discreteness of the material distribution", HistoryFieldType::COEFFICIENT);
  }
  AddHistoryOutput("COMBO", "ComboObj", ScreenOutputFormat::SCIENTIFIC, "COMBO", "Combined obj. function value.", HistoryFieldType::COEFFICIENT);

  if (coupled_heat) {
    CHeatOutput::SetHistoryOutputFieldsImpl(config, this);
    AddHistoryOutput("LINSOL_ITER_HEAT", "LinSolIterHeat", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
    AddHistoryOutput("LINSOL_RESIDUAL_HEAT", "LinSolResHeat", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");
  }
}

void CElasticityOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_Struc = solver[FEA_SOL]->GetNodes();
  CPoint*    Node_Geo  = geometry->nodes;

  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(iPoint, 0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(iPoint, 2));

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
  if (coupled_heat) {
    CVariable* Node_Heat = solver[HEAT_SOL]->GetNodes();
    SetVolumeOutputValue("TEMPERATURE", iPoint, Node_Heat->GetSolution(iPoint, 0));
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

  if (config->GetTopology_Optimization()) {
    SetVolumeOutputValue("TOPOL_DENSITY", iPoint, Node_Struc->GetAuxVar(iPoint));
  }

  CSolver* heat_solver = solver[HEAT_SOL];
  if (heat_solver) {
    const auto Node_Heat = heat_solver->GetNodes();
    SetVolumeOutputValue("TEMPERATURE", iPoint, Node_Heat->GetSolution(iPoint, 0));
    SetVolumeOutputValue("RES_TEMPERATURE", iPoint, heat_solver->LinSysRes(iPoint, 0));
  }

}

void CElasticityOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3) AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");

  AddVolumeOutput("DISPLACEMENT-X",    "Displacement_x", "SOLUTION", "x-component of the displacement vector");
  AddVolumeOutput("DISPLACEMENT-Y",    "Displacement_y", "SOLUTION", "y-component of the displacement vector");
  if (nDim == 3) AddVolumeOutput("DISPLACEMENT-Z", "Displacement_z", "SOLUTION", "z-component of the displacement vector");

  if (dynamic) {
    AddVolumeOutput("VELOCITY-X",    "Velocity_x", "SOLUTION", "x-component of the velocity vector");
    AddVolumeOutput("VELOCITY-Y",    "Velocity_y", "SOLUTION", "y-component of the velocity vector");
    if (nDim == 3) AddVolumeOutput("VELOCITY-Z", "Velocity_z", "SOLUTION", "z-component of the velocity vector");

    AddVolumeOutput("ACCELERATION-X",    "Acceleration_x", "SOLUTION", "x-component of the acceleration vector");
    AddVolumeOutput("ACCELERATION-Y",    "Acceleration_y", "SOLUTION", "y-component of the acceleration vector");
    if (nDim == 3) AddVolumeOutput("ACCELERATION-Z", "Acceleration_z", "SOLUTION", "z-component of the acceleration vector");
  }

  if (coupled_heat) {
    AddVolumeOutput("TEMPERATURE", "Temperature", "SOLUTION", "Temperature");
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

  if (config->GetTopology_Optimization()) {
    AddVolumeOutput("TOPOL_DENSITY", "Topology_Density", "TOPOLOGY", "filtered topology density");
  }

  if (coupled_heat) {
    AddVolumeOutput("HEAT_FLUX", "Heat_Flux", "PRIMITIVE", "Heatflux");
    AddVolumeOutput("RES_TEMPERATURE", "Residual_Temperature", "RESIDUAL", "Residual of the temperature");
  }

}

bool CElasticityOutput::SetInitResiduals(const CConfig *config){

  return (config->GetTime_Domain() == NO && (curInnerIter  == 0));

}

void CElasticityOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint,
                                        unsigned short iMarker, unsigned long iVertex) {
  if (!coupled_heat || !config->GetViscous_Wall(iMarker)) return;

  /* Heat flux value at each surface grid node. */
  SetVolumeOutputValue("HEAT_FLUX", iPoint, solver[HEAT_SOL]->GetHeatFlux(iMarker, iVertex));

}
