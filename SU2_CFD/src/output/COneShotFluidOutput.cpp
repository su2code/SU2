/*!
 * \file output_one_shot.cpp
 * \brief Main subroutines for compressible flow output
 * \author B. Munguía
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

#include "../../include/output/COneShotFluidOutput.hpp"

#include "../../../Common/include/geometry_structure.hpp"
#include "../../include/solver_structure.hpp"

COneShotFluidOutput::COneShotFluidOutput(CConfig *config, unsigned short nDim) : CFlowOutput(config, nDim, false) {
  
  turb_model = config->GetKind_Turb_Model();
  
  gridMovement = config->GetGrid_Movement(); 

  /*--- Set the default history fields if nothing is set in the config file ---*/
  
  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.emplace_back("ITER");
    requestedHistoryFields.emplace_back("SEARCH_ITER");
    requestedScreenFields.emplace_back("ONE_SHOT_RHO");
    requestedScreenFields.emplace_back("ONE_SHOT_THETA");
    requestedHistoryFields.emplace_back("GRAD_NORM");
    requestedHistoryFields.emplace_back("RMS_RES");
    requestedHistoryFields.emplace_back("SENSITIVITY");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }
  if (nRequestedScreenFields == 0){
    if (config->GetTime_Domain()) requestedScreenFields.emplace_back("TIME_ITER");
    if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("SEARCH_ITER");
    requestedScreenFields.emplace_back("ONE_SHOT_RHO");
    requestedScreenFields.emplace_back("ONE_SHOT_THETA");
    requestedHistoryFields.emplace_back("GRAD_NORM");
    requestedScreenFields.emplace_back("RMS_DENSITY");
    requestedScreenFields.emplace_back("RMS_MOMENTUM-X");
    requestedScreenFields.emplace_back("RMS_ADJ_DENSITY");
    requestedScreenFields.emplace_back("RMS_ADJ_MOMENTUM-X");
    requestedScreenFields.emplace_back("SENS_GEO");
    requestedScreenFields.emplace_back("SENS_AOA");
    nRequestedScreenFields = requestedScreenFields.size();
  }
  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    requestedVolumeFields.emplace_back("PRIMITIVE");
    requestedVolumeFields.emplace_back("SENSITIVITY");
    nRequestedVolumeFields = requestedVolumeFields.size();
  }
  
  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (One-Shot Comp. Fluid)";
  multiZoneHeaderString = ss.str();
  
  /*--- Set the volume filename --- */
  
  volumeFilename = config->GetVolume_FileName();
  
  /*--- Set the surface filename --- */
  
  surfaceFilename = config->GetSurfCoeff_FileName();
  
  /*--- Set the restart filename --- */
  
  restartFilename = config->GetRestart_FileName();


  /*--- Set the default convergence field --- */

  if (convFields.empty() ) convFields.emplace_back("RMS_DENSITY");
  
}

COneShotFluidOutput::~COneShotFluidOutput(void) {}



void COneShotFluidOutput::SetHistoryOutputFields(CConfig *config){

  /// BEGIN_GROUP: ONE_SHOT, DESCRIPTION: One-shot specific outputs.
  /// DESCRIPTION: Number of Armijo search iterations.
  AddHistoryOutput("SEARCH_ITER",    "Search_Iter", ScreenOutputFormat::INTEGER, "ONE_SHOT", "Number of Armijo search iterations."); 
  /// DESCRIPTION: Number of DVs in the active set.
  AddHistoryOutput("ACTIVE_DV",      "Active_DV",   ScreenOutputFormat::INTEGER, "ONE_SHOT", "Number of DVs in the active set."); 
  /// DESCRIPTION: Primal contraction ratio.
  AddHistoryOutput("ONE_SHOT_RHO",   "Rho_G_u",     ScreenOutputFormat::SCIENTIFIC, "ONE_SHOT", "Primal contraction ratio."); 
  /// DESCRIPTION: Adjoint contraction ratio.
  AddHistoryOutput("ONE_SHOT_THETA", "Theta_N_uu",  ScreenOutputFormat::SCIENTIFIC, "ONE_SHOT", "Adjoint contraction ratio."); 
  /// DESCRIPTION: Norm of the shifted Lagrangian gradient.
  AddHistoryOutput("GRAD_NORM",      "||N_x||",     ScreenOutputFormat::SCIENTIFIC, "ONE_SHOT", "Norm of the shifted Lagrangian gradient."); 
  for(unsigned short iConstr = 0; iConstr < config->GetnConstr(); iConstr++) {
    std::stringstream ss;
    ss << "LAMBDA_" << iConstr;
    std::string LambdaString = ss.str();
    ss.str(std::string());
    ss.clear();
    ss << "LAMBDA[" << iConstr << "]";
    std::string LambdaHeader = ss.str();
    /// DESCRIPTION: Lagrange multiplier.
    AddHistoryOutput(LambdaString, LambdaHeader, ScreenOutputFormat::SCIENTIFIC, "ONE_SHOT", "Lagrange multiplier value.");
  }

  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the SOLUTION variables. 
  /// DESCRIPTION: Root-mean square residual of the density.
  AddHistoryOutput("RMS_DENSITY",    "rms[Rho]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum x-component.
  AddHistoryOutput("RMS_MOMENTUM-X", "rms[RhoU]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum y-component.
  AddHistoryOutput("RMS_MOMENTUM-Y", "rms[RhoV]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the momentum y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum z-component.
  if (nDim == 3) AddHistoryOutput("RMS_MOMENTUM-Z", "rms[RhoW]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the momentum z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the energy.
  AddHistoryOutput("RMS_ENERGY",     "rms[RhoE]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the energy.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint density.
  AddHistoryOutput("RMS_ADJ_DENSITY",    "rms[A_Rho]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint density.", HistoryFieldType::RESIDUAL); 
  /// DESCRIPTION: Root-mean square residual of the adjoint momentum x-component.
  AddHistoryOutput("RMS_ADJ_MOMENTUM-X", "rms[A_RhoU]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint momentum y-component.
  AddHistoryOutput("RMS_ADJ_MOMENTUM-Y", "rms[A_RhoV]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint momentum y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint momentum z-component.
  AddHistoryOutput("RMS_ADJ_MOMENTUM-Z", "rms[A_RhoW]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint momentum z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the adjoint energy.
  AddHistoryOutput("RMS_ADJ_ENERGY",     "rms[A_E]",    ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint energy.", HistoryFieldType::RESIDUAL); 
  
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Root-mean square residual of nu tilde (SA model).  
    AddHistoryOutput("RMS_NU_TILDE", "rms[nu]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of nu tilde (SA model).", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: Root-mean square residual of the adjoint nu tilde.
    AddHistoryOutput("RMS_ADJ_NU_TILDE", "rms[A_nu]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint nu tilde.", HistoryFieldType::RESIDUAL);
    break;  
  case SST: case SST_SUST:
    /// DESCRIPTION: Root-mean square residual of kinetic energy (SST model).    
    AddHistoryOutput("RMS_TKE", "rms[k]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of kinetic energy (SST model).", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: Root-mean square residual of the dissipation (SST model).    
    AddHistoryOutput("RMS_DISSIPATION", "rms[w]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of dissipation (SST model).", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: Root-mean square residual of the adjoint kinetic energy.
    AddHistoryOutput("RMS_ADJ_TKE", "rms[A_k]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint kinetic energy.", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: Root-mean square residual of the adjoint dissipation.
    AddHistoryOutput("RMS_ADJ_DISSIPATION",    "rms[A_w]", ScreenOutputFormat::FIXED, "RMS_RES", " Root-mean square residual of the adjoint dissipation.", HistoryFieldType::RESIDUAL);   
    break;
  default: break;
  }
  /// END_GROUP
   
  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the SOLUTION variables. 
  /// DESCRIPTION: Maximum residual of the density.
  AddHistoryOutput("MAX_DENSITY",    "max[Rho]",  ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum square residual of the density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum x-component. 
  AddHistoryOutput("MAX_MOMENTUM-X", "max[RhoU]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum square residual of the momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum y-component. 
  AddHistoryOutput("MAX_MOMENTUM-Y", "max[RhoV]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum square residual of the momentum y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum z-component. 
  if (nDim == 3) AddHistoryOutput("MAX_MOMENTUM-Z", "max[RhoW]", ScreenOutputFormat::FIXED,"MAX_RES", "Maximum residual of the z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the energy.  
  AddHistoryOutput("MAX_ENERGY",     "max[RhoE]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum residual of the energy.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the adjoint density.
  AddHistoryOutput("MAX_ADJ_DENSITY",    "max[A_Rho]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the adjoint momentum x-component
  AddHistoryOutput("MAX_ADJ_MOMENTUM-X", "max[A_RhoU]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint momentum x-component", HistoryFieldType::RESIDUAL); 
  /// DESCRIPTION: Maximum residual of the adjoint momentum y-component
  AddHistoryOutput("MAX_ADJ_MOMENTUM-Y", "max[A_RhoV]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint momentum y-component", HistoryFieldType::RESIDUAL); 
  /// DESCRIPTION: Maximum residual of the adjoint momentum z-component
  AddHistoryOutput("MAX_ADJ_MOMENTUM-Z", "max[A_RhoW]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint momentum z-component", HistoryFieldType::RESIDUAL); 
  /// DESCRIPTION: Maximum residual of the adjoint energy.
  AddHistoryOutput("MAX_ADJ_ENERGY",     "max[A_E]",    ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint energy.", HistoryFieldType::RESIDUAL);
  
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Maximum residual of nu tilde (SA model).
    AddHistoryOutput("MAX_NU_TILDE",       "max[nu]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of nu tilde (SA model).", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: Maximum residual of the adjoint nu tilde.
    AddHistoryOutput("MAX_ADJ_NU_TILDE", "max[A_nu]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint nu tilde.", HistoryFieldType::RESIDUAL);
    break;  
  case SST: case SST_SUST:
    /// DESCRIPTION: Maximum residual of kinetic energy (SST model). 
    AddHistoryOutput("MAX_TKE", "max[k]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of kinetic energy (SST model).", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: Maximum residual of the dissipation (SST model).   
    AddHistoryOutput("MAX_DISSIPATION",    "max[w]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of dissipation (SST model).", HistoryFieldType::RESIDUAL); 
    /// DESCRIPTION: Maximum residual of the adjoint kinetic energy.
    AddHistoryOutput("MAX_ADJ_TKE", "max[A_k]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint kinetic energy.", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: Maximum residual of the adjoint dissipation.
    AddHistoryOutput("MAX_ADJ_DISSIPATION",    "max[A_w]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint dissipation.", HistoryFieldType::RESIDUAL);
    break;
  default: break;
  }
  /// END_GROUP

  /// DESCRIPTION: Linear solver iterations   
  AddHistoryOutput("LINSOL_ITER", "Linear_Solver_Iterations", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");
 
  /// BEGIN_GROUP: ENGINE_OUTPUT, DESCRIPTION: Engine output
  /// DESCRIPTION: Aero CD drag
  AddHistoryOutput("AEROCDRAG",                  "AeroCDrag",                  ScreenOutputFormat::SCIENTIFIC, "ENGINE_OUTPUT", "Aero CD drag", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Solid CD drag  
  AddHistoryOutput("SOLIDCDRAG",                 "SolidCDrag",                 ScreenOutputFormat::SCIENTIFIC, "ENGINE_OUTPUT", "Solid CD drag ", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Radial distortion 
  AddHistoryOutput("RADIAL_DISTORTION",          "Radial_Distortion",          ScreenOutputFormat::SCIENTIFIC, "ENGINE_OUTPUT", "Radial distortion ", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Circumferential distortion
  AddHistoryOutput("CIRCUMFERENTIAL_DISTORTION", "Circumferential_Distortion", ScreenOutputFormat::SCIENTIFIC, "ENGINE_OUTPUT", "Circumferential distortion", HistoryFieldType::COEFFICIENT);
  /// END_GROUP
  
  /// BEGIN_GROUP: ROTATING_FRAME, DESCRIPTION: Coefficients related to a rotating frame of reference.
  /// DESCRIPTION: Merit  
  AddHistoryOutput("MERIT", "CMerit", ScreenOutputFormat::SCIENTIFIC, "ROTATING_FRAME", "Merit", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: CT 
  AddHistoryOutput("CT",    "CT",     ScreenOutputFormat::SCIENTIFIC, "ROTATING_FRAME", "CT", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: CQ  
  AddHistoryOutput("CQ",    "CQ",     ScreenOutputFormat::SCIENTIFIC, "ROTATING_FRAME", "CQ", HistoryFieldType::COEFFICIENT);
  /// END_GROUP
  
  /// BEGIN_GROUP: EQUIVALENT_AREA, DESCRIPTION: Equivalent area.  
  /// DESCRIPTION: Equivalent area    
  AddHistoryOutput("EQUIV_AREA",   "CEquiv_Area",  ScreenOutputFormat::SCIENTIFIC, "EQUIVALENT_AREA", "Equivalent area", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Nearfield obj. function      
  AddHistoryOutput("NEARFIELD_OF", "CNearFieldOF", ScreenOutputFormat::SCIENTIFIC, "EQUIVALENT_AREA", "Nearfield obj. function ", HistoryFieldType::COEFFICIENT);
  /// END_GROUP

  ///   /// BEGIN_GROUP: HEAT_COEFF, DESCRIPTION: Heat coefficients on all surfaces set with MARKER_MONITORING.
  /// DESCRIPTION: Total heatflux
  AddHistoryOutput("HEATFLUX", "HF",      ScreenOutputFormat::SCIENTIFIC, "HEAT", "Total heatflux on all surfaces set with MARKER_MONITORING.", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Maximal heatflux  
  AddHistoryOutput("HEATFLUX_MAX", "maxHF",    ScreenOutputFormat::SCIENTIFIC, "HEAT", "Total maximum heatflux on all surfaces set with MARKER_MONITORING.", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Temperature
  AddHistoryOutput("TEMPERATURE", "Temp", ScreenOutputFormat::SCIENTIFIC, "HEAT",  "Total avg. temperature on all surfaces set with MARKER_MONITORING.", HistoryFieldType::COEFFICIENT);
  /// END_GROUP
  
  AddHistoryOutput("CFL_NUMBER", "CFL number", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current value of the CFL number");
  
  /*--- Add analyze surface history fields --- */
  
  AddAnalyzeSurfaceOutput(config);
  
  /*--- Add aerodynamic coefficients fields --- */
  
  AddAerodynamicCoefficients(config); 
  
  /*--- Add Cp diff fields ---*/
  
  Add_CpInverseDesignOutput(config);
  
  /*--- Add combo obj value --- */
  
  AddHistoryOutput("COMBO", "ComboObj", ScreenOutputFormat::SCIENTIFIC, "COMBO", "Combined obj. function value.", HistoryFieldType::COEFFICIENT);

  /// BEGIN_GROUP: SENSITIVITY, DESCRIPTION: Sensitivities of different geometrical or boundary values.   
  /// DESCRIPTION: Sum of the geometrical sensitivities on all markers set in MARKER_MONITORING.
  AddHistoryOutput("SENS_GEO",   "Sens_Geo",   ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sum of the geometrical sensitivities on all markers set in MARKER_MONITORING.", HistoryFieldType::COEFFICIENT); 
  /// DESCRIPTION: Sensitivity of the objective function with respect to the angle of attack (only for compressible solver).
  AddHistoryOutput("SENS_AOA",   "Sens_AoA",   ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the angle of attack (only for compressible solver).", HistoryFieldType::COEFFICIENT); 
  /// DESCRIPTION: Sensitivity of the objective function with respect to the Mach number (only of compressible solver).
  AddHistoryOutput("SENS_MACH",  "Sens_Mach",  ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the Mach number (only of compressible solver).", HistoryFieldType::COEFFICIENT); 
  /// DESCRIPTION: Sensitivity of the objective function with respect to the far-field pressure.
  AddHistoryOutput("SENS_PRESS", "Sens_Press", ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the far-field pressure.", HistoryFieldType::COEFFICIENT); 
  /// DESCRIPTION: Sensitivity of the objective function with respect to the far-field temperature.
  AddHistoryOutput("SENS_TEMP",  "Sens_Temp",  ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Sensitivity of the objective function with respect to the far-field temperature.", HistoryFieldType::COEFFICIENT); 
  /// END_GROUP
}

void COneShotFluidOutput::SetVolumeOutputFields(CConfig *config){
  
  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");

  // Solution variables
  AddVolumeOutput("DENSITY",    "Density",    "SOLUTION", "Density");
  AddVolumeOutput("MOMENTUM-X", "Momentum_x", "SOLUTION", "x-component of the momentum vector");
  AddVolumeOutput("MOMENTUM-Y", "Momentum_y", "SOLUTION", "y-component of the momentum vector");
  if (nDim == 3)
    AddVolumeOutput("MOMENTUM-Z", "Momentum_z", "SOLUTION", "z-component of the momentum vector");
  AddVolumeOutput("ENERGY",     "Energy",     "SOLUTION", "Energy");  
  /// DESCRIPTION: Adjoint density.
  AddVolumeOutput("ADJ_DENSITY",    "Adjoint_Density",    "SOLUTION", "Adjoint density"); 
  /// DESCRIPTION: Adjoint momentum x-component.
  AddVolumeOutput("ADJ_MOMENTUM-X", "Adjoint_Momentum_x", "SOLUTION", "x-component of the adjoint momentum vector"); 
  /// DESCRIPTION: Adjoint momentum y-component.
  AddVolumeOutput("ADJ_MOMENTUM-Y", "Adjoint_Momentum_y", "SOLUTION", "y-component of the adjoint momentum vector"); 
  if (nDim == 3)
    /// DESCRIPTION: Adjoint momentum z-component.
    AddVolumeOutput("ADJ_MOMENTUM-Z", "Adjoint_Momentum_z", "SOLUTION", "z-component of the adjoint momentum vector"); 
  /// DESCRIPTION: Adjoint energy.
  AddVolumeOutput("ADJ_ENERGY", "Adjoint_Energy", "SOLUTION", "Adjoint energy");  
  
  // Turbulent Residuals
  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    AddVolumeOutput("TKE", "Turb_Kin_Energy", "SOLUTION", "Turbulent kinetic energy");
    AddVolumeOutput("DISSIPATION", "Omega", "SOLUTION", "Rate of dissipation");
    /// DESCRIPTION: Adjoint kinetic energy.
    AddVolumeOutput("ADJ_TKE", "Adjoint_TKE", "SOLUTION", "Adjoint kinetic energy");
    /// DESCRIPTION: Adjoint dissipation.
    AddVolumeOutput("ADJ_DISSIPATION", "Adjoint_Omega", "SOLUTION", "Adjoint rate of dissipation");
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    AddVolumeOutput("NU_TILDE", "Nu_Tilde", "SOLUTION", "Spalart-Allmaras variable");
    /// DESCRIPTION: Adjoint nu tilde.
    AddVolumeOutput("ADJ_NU_TILDE", "Adjoint_Nu_Tilde", "SOLUTION", "Adjoint Spalart-Allmaras variable"); 
    break;
  case NONE:
    break;
  }
  
  // Grid velocity
  if (config->GetGrid_Movement()){
    AddVolumeOutput("GRID_VELOCITY-X", "Grid_Velocity_x", "GRID_VELOCITY", "x-component of the grid velocity vector");
    AddVolumeOutput("GRID_VELOCITY-Y", "Grid_Velocity_y", "GRID_VELOCITY", "y-component of the grid velocity vector");
    if (nDim == 3 ) 
      AddVolumeOutput("GRID_VELOCITY-Z", "Grid_Velocity_z", "GRID_VELOCITY", "z-component of the grid velocity vector");
  }
  
  // Primitive variables
  AddVolumeOutput("PRESSURE",    "Pressure",                "PRIMITIVE", "Pressure");
  AddVolumeOutput("TEMPERATURE", "Temperature",             "PRIMITIVE", "Temperature");
  AddVolumeOutput("MACH",        "Mach",                    "PRIMITIVE", "Mach number");
  AddVolumeOutput("PRESSURE_COEFF", "Pressure_Coefficient", "PRIMITIVE", "Pressure coefficient");
  
  if (config->GetKind_Solver() == RANS || config->GetKind_Solver() == NAVIER_STOKES){
    AddVolumeOutput("LAMINAR_VISCOSITY", "Laminar_Viscosity", "PRIMITIVE", "Laminar viscosity");
    
    AddVolumeOutput("SKIN_FRICTION-X", "Skin_Friction_Coefficient_x", "PRIMITIVE", "x-component of the skin friction vector");
    AddVolumeOutput("SKIN_FRICTION-Y", "Skin_Friction_Coefficient_y", "PRIMITIVE", "y-component of the skin friction vector");
    if (nDim == 3)
      AddVolumeOutput("SKIN_FRICTION-Z", "Skin_Friction_Coefficient_z", "PRIMITIVE", "z-component of the skin friction vector");
    
    AddVolumeOutput("HEAT_FLUX", "Heat_Flux", "PRIMITIVE", "Heat-flux");
    AddVolumeOutput("Y_PLUS", "Y_Plus", "PRIMITIVE", "Non-dim. wall distance (Y-Plus)");
    
  }
  
  if (config->GetKind_Solver() == RANS) {
    AddVolumeOutput("EDDY_VISCOSITY", "Eddy_Viscosity", "PRIMITIVE", "Turbulent eddy viscosity");
  }
  
  if (config->GetKind_Trans_Model() == BC){
    AddVolumeOutput("INTERMITTENCY", "gamma_BC", "INTERMITTENCY", "Intermittency");
  }

  //Residuals
  AddVolumeOutput("RES_DENSITY", "Residual_Pressure", "RESIDUAL", "Residual of the density");
  AddVolumeOutput("RES_MOMENTUM-X", "Residual_Momentum_x", "RESIDUAL", "Residual of the x-momentum component");
  AddVolumeOutput("RES_MOMENTUM-Y", "Residual_Momentum_y", "RESIDUAL", "Residual of the y-momentum component");
  if (nDim == 3)
    AddVolumeOutput("RES_MOMENTUM-Z", "Residual_Momentum_z", "RESIDUAL", "Residual of the z-momentum component");
  AddVolumeOutput("RES_ENERGY", "Residual_Energy", "RESIDUAL", "Residual of the energy");
  /// DESCRIPTION: Residual of the adjoint density.
  AddVolumeOutput("RES_ADJ_DENSITY",    "Residual_Adjoint_Density",    "RESIDUAL", "Residual of the adjoint density");  
  /// DESCRIPTION: Residual of the adjoint momentum x-component.
  AddVolumeOutput("RES_ADJ_MOMENTUM-X", "Residual_Adjoint_Momentum_x", "RESIDUAL", "Residual of the adjoint x-momentum");  
  /// DESCRIPTION: Residual of the adjoint momentum y-component.
  AddVolumeOutput("RES_ADJ_MOMENTUM-Y", "Residual_Adjoint_Momentum_y", "RESIDUAL", "Residual of the adjoint y-momentum");  
  if (nDim == 3)
    /// DESCRIPTION: Residual of the adjoint momentum z-component.
    AddVolumeOutput("RES_ADJ_MOMENTUM-Z", "Residual_Adjoint_Momentum_z", "RESIDUAL", "Residual of the adjoint z-momentum"); 
  /// DESCRIPTION: Residual of the adjoint energy. 
  AddVolumeOutput("RES_ADJ_ENERGY", "Residual_Adjoint_Energy", "RESIDUAL", "Residual of the adjoint energy");  
  
  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    AddVolumeOutput("RES_TKE", "Residual_TKE", "RESIDUAL", "Residual of turbulent kinetic energy");
    AddVolumeOutput("RES_DISSIPATION", "Residual_Omega", "RESIDUAL", "Residual of the rate of dissipation");
    /// DESCRIPTION: Residual of the adjoint kinetic energy. 
    AddVolumeOutput("RES_ADJ_TKE", "Residual_Adjoint_TKE", "RESIDUAL", "Residual of the turb. kinetic energy");
    /// DESCRIPTION: Residual of the adjoint dissipation.
    AddVolumeOutput("RES_ADJ_DISSIPATION", "Residual_Adjoint_Omega", "RESIDUAL", "Residual of the rate of dissipation");
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    AddVolumeOutput("RES_NU_TILDE", "Residual_Nu_Tilde", "RESIDUAL", "Residual of the Spalart-Allmaras variable");
    /// DESCRIPTION: Residual of the nu tilde. 
    AddVolumeOutput("RES_ADJ_NU_TILDE", "Residual_Adjoint_Nu_Tilde", "RESIDUAL", "Residual of the Spalart-Allmaras variable"); 
    break;
  case NONE:
    break;
  }
  
  // Limiter values
  AddVolumeOutput("LIMITER_DENSITY", "Limiter_Density", "LIMITER", "Limiter value of the density");
  AddVolumeOutput("LIMITER_MOMENTUM-X", "Limiter_Momentum_x", "LIMITER", "Limiter value of the x-momentum");
  AddVolumeOutput("LIMITER_MOMENTUM-Y", "Limiter_Momentum_y", "LIMITER", "Limiter value of the y-momentum");
  if (nDim == 3)
    AddVolumeOutput("LIMITER_MOMENTUM-Z", "Limiter_Momentum_z", "LIMITER", "Limiter value of the z-momentum");
  AddVolumeOutput("LIMITER_ENERGY", "Limiter_Energy", "LIMITER", "Limiter value of the energy");
  
  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    AddVolumeOutput("LIMITER_TKE", "Limiter_TKE", "LIMITER", "Limiter value of turb. kinetic energy");
    AddVolumeOutput("LIMITER_DISSIPATION", "Limiter_Omega", "LIMITER", "Limiter value of dissipation rate");
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    AddVolumeOutput("LIMITER_NU_TILDE", "Limiter_Nu_Tilde", "LIMITER", "Limiter value of the Spalart-Allmaras variable");
    break;
  case NONE:
    break;
  }

  
  // Hybrid RANS-LES
  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
    AddVolumeOutput("DES_LENGTHSCALE", "DES_LengthScale", "DDES", "DES length scale value");
    AddVolumeOutput("WALL_DISTANCE", "Wall_Distance", "DDES", "Wall distance value");
  }
  
  // Roe Low Dissipation
  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    AddVolumeOutput("ROE_DISSIPATION", "Roe_Dissipation", "ROE_DISSIPATION", "Value of the Roe dissipation");
  }
  
  if(config->GetKind_Solver() == RANS || config->GetKind_Solver() == NAVIER_STOKES){
    if (nDim == 3){
      AddVolumeOutput("VORTICITY_X", "Vorticity_x", "VORTEX_IDENTIFICATION", "x-component of the vorticity vector");
      AddVolumeOutput("VORTICITY_Y", "Vorticity_y", "VORTEX_IDENTIFICATION", "y-component of the vorticity vector");
      AddVolumeOutput("Q_CRITERION", "Q_Criterion", "VORTEX_IDENTIFICATION", "Value of the Q-Criterion");      
    }
    AddVolumeOutput("VORTICITY_Z", "Vorticity_z", "VORTEX_IDENTIFICATION", "z-component of the vorticity vector");
  }

  // Anisotropic metric
  if(config->GetError_Estimate() || config->GetKind_SU2() == SU2_MET) {
    if (nDim == 2){
      AddVolumeOutput("ANISO_METR0", "Aniso_Metric[0]", "ANISOTROPIC_METRIC", "x-x-component of the metric");
      AddVolumeOutput("ANISO_METR1", "Aniso_Metric[1]", "ANISOTROPIC_METRIC", "x-y-component of the metric");
      AddVolumeOutput("ANISO_METR2", "Aniso_Metric[2]", "ANISOTROPIC_METRIC", "y-y-component of the metric");
    }
    else{
      AddVolumeOutput("ANISO_METR0", "Aniso_Metric[0]", "ANISOTROPIC_METRIC", "x-x-component of the metric");
      AddVolumeOutput("ANISO_METR1", "Aniso_Metric[1]", "ANISOTROPIC_METRIC", "x-y-component of the metric");
      AddVolumeOutput("ANISO_METR2", "Aniso_Metric[2]", "ANISOTROPIC_METRIC", "x-z-component of the metric");
      AddVolumeOutput("ANISO_METR3", "Aniso_Metric[3]", "ANISOTROPIC_METRIC", "y-y-component of the metric");
      AddVolumeOutput("ANISO_METR4", "Aniso_Metric[4]", "ANISOTROPIC_METRIC", "y-z-component of the metric");
      AddVolumeOutput("ANISO_METR5", "Aniso_Metric[5]", "ANISOTROPIC_METRIC", "z-z-component of the metric");
    }
  }

  /// BEGIN_GROUP: SENSITIVITY, DESCRIPTION: Geometrical sensitivities of the current objective function.
  /// DESCRIPTION: Sensitivity x-component.
  AddVolumeOutput("SENSITIVITY-X", "Sensitivity_x", "SENSITIVITY", "x-component of the sensitivity vector"); 
  /// DESCRIPTION: Sensitivity y-component.
  AddVolumeOutput("SENSITIVITY-Y", "Sensitivity_y", "SENSITIVITY", "y-component of the sensitivity vector");
  if (nDim == 3)
    /// DESCRIPTION: Sensitivity z-component.
    AddVolumeOutput("SENSITIVITY-Z", "Sensitivity_z", "SENSITIVITY", "z-component of the sensitivity vector");   
  /// DESCRIPTION: Sensitivity in normal direction.
  AddVolumeOutput("SENSITIVITY", "Surface_Sensitivity", "SENSITIVITY", "sensitivity in normal direction"); 
  /// END_GROUP
  
  if (config->GetTime_Domain()){
    SetTimeAveragedFields();
  }
}

void COneShotFluidOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){
  
  CVariable* Node_Flow = solver[FLOW_SOL]->GetNodes(); 
  CVariable* Node_Turb = NULL;
  CVariable* Node_AdjFlow = solver[ADJFLOW_SOL]->GetNodes(); 
  CVariable* Node_AdjTurb = NULL;
  CPoint*    Node_Geo     = geometry->node[iPoint];
  
  if (config->GetKind_Turb_Model() != NONE){
    Node_Turb = solver[TURB_SOL]->GetNodes(); 
    Node_AdjTurb = solver[ADJTURB_SOL]->GetNodes(); 
  }
            
  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));  
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));
  
  SetVolumeOutputValue("DENSITY",    iPoint, Node_Flow->GetSolution(iPoint, 0));
  SetVolumeOutputValue("MOMENTUM-X", iPoint, Node_Flow->GetSolution(iPoint, 1));
  SetVolumeOutputValue("MOMENTUM-Y", iPoint, Node_Flow->GetSolution(iPoint, 2));
  if (nDim == 3){
    SetVolumeOutputValue("MOMENTUM-Z", iPoint, Node_Flow->GetSolution(iPoint, 3));
    SetVolumeOutputValue("ENERGY",     iPoint, Node_Flow->GetSolution(iPoint, 4));
  } else {
    SetVolumeOutputValue("ENERGY",     iPoint, Node_Flow->GetSolution(iPoint, 3));    
  }

  SetVolumeOutputValue("ADJ_DENSITY",    iPoint, Node_AdjFlow->GetSolution(iPoint, 0));
  SetVolumeOutputValue("ADJ_MOMENTUM-X", iPoint, Node_AdjFlow->GetSolution(iPoint, 1));
  SetVolumeOutputValue("ADJ_MOMENTUM-Y", iPoint, Node_AdjFlow->GetSolution(iPoint, 2));
  if (nDim == 3){
    SetVolumeOutputValue("ADJ_MOMENTUM-Z", iPoint, Node_AdjFlow->GetSolution(iPoint, 3));
    SetVolumeOutputValue("ADJ_ENERGY",     iPoint, Node_AdjFlow->GetSolution(iPoint, 4));
  } else {
    SetVolumeOutputValue("ADJ_ENERGY",     iPoint, Node_AdjFlow->GetSolution(iPoint, 3));    
  }
  
  // Turbulent Residuals
  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    SetVolumeOutputValue("TKE",         iPoint, Node_Turb->GetSolution(iPoint, 0));
    SetVolumeOutputValue("DISSIPATION", iPoint, Node_Turb->GetSolution(iPoint, 1));
    SetVolumeOutputValue("ADJ_TKE",         iPoint, Node_AdjTurb->GetSolution(iPoint, 0));
    SetVolumeOutputValue("ADJ_DISSIPATION", iPoint, Node_AdjTurb->GetSolution(iPoint, 1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("NU_TILDE", iPoint, Node_Turb->GetSolution(iPoint, 0));
    SetVolumeOutputValue("ADJ_NU_TILDE", iPoint, Node_AdjTurb->GetSolution(iPoint, 0));
    break;
  case NONE:
    break;
  }
  
  if (config->GetGrid_Movement()){
    SetVolumeOutputValue("GRID_VELOCITY-X", iPoint, Node_Geo->GetGridVel()[0]);
    SetVolumeOutputValue("GRID_VELOCITY-Y", iPoint, Node_Geo->GetGridVel()[1]);
    if (nDim == 3)
      SetVolumeOutputValue("GRID_VELOCITY-Z", iPoint, Node_Geo->GetGridVel()[2]);
  }
  
  SetVolumeOutputValue("PRESSURE", iPoint, Node_Flow->GetPressure(iPoint));
  SetVolumeOutputValue("TEMPERATURE", iPoint, Node_Flow->GetTemperature(iPoint));
  SetVolumeOutputValue("MACH", iPoint, sqrt(Node_Flow->GetVelocity2(iPoint))/Node_Flow->GetSoundSpeed(iPoint));

  su2double VelMag = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++){
    VelMag += pow(solver[FLOW_SOL]->GetVelocity_Inf(iDim),2.0); 
  }
  su2double factor = 1.0/(0.5*solver[FLOW_SOL]->GetDensity_Inf()*VelMag); 
  SetVolumeOutputValue("PRESSURE_COEFF", iPoint, (Node_Flow->GetPressure(iPoint) - solver[FLOW_SOL]->GetPressure_Inf())*factor);
  
  if (config->GetKind_Solver() == RANS || config->GetKind_Solver() == NAVIER_STOKES){
    SetVolumeOutputValue("LAMINAR_VISCOSITY", iPoint, Node_Flow->GetLaminarViscosity(iPoint));
  }
  
  if (config->GetKind_Solver() == RANS) {
    SetVolumeOutputValue("EDDY_VISCOSITY", iPoint, Node_Flow->GetEddyViscosity(iPoint));
  }
  
  if (config->GetKind_Trans_Model() == BC){
    SetVolumeOutputValue("INTERMITTENCY", iPoint, Node_Turb->GetGammaBC(iPoint));
  }
  
  SetVolumeOutputValue("RES_DENSITY", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 0));
  SetVolumeOutputValue("RES_MOMENTUM-X", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 1));
  SetVolumeOutputValue("RES_MOMENTUM-Y", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 2));
  if (nDim == 3){
    SetVolumeOutputValue("RES_MOMENTUM-Z", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 3));
    SetVolumeOutputValue("RES_ENERGY", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 4));
  } else {
    SetVolumeOutputValue("RES_ENERGY", iPoint, solver[FLOW_SOL]->LinSysRes.GetBlock(iPoint, 3));
  }

  SetVolumeOutputValue("RES_ADJ_DENSITY",    iPoint, Node_AdjFlow->GetSolution(iPoint, 0) - Node_AdjFlow->GetSolution_Old(iPoint, 0));
  SetVolumeOutputValue("RES_ADJ_MOMENTUM-X", iPoint, Node_AdjFlow->GetSolution(iPoint, 1) - Node_AdjFlow->GetSolution_Old(iPoint, 1));
  SetVolumeOutputValue("RES_ADJ_MOMENTUM-Y", iPoint, Node_AdjFlow->GetSolution(iPoint, 2) - Node_AdjFlow->GetSolution_Old(iPoint, 2));
  if (nDim == 3){
    SetVolumeOutputValue("RES_ADJ_MOMENTUM-Z", iPoint, Node_AdjFlow->GetSolution(iPoint, 3) - Node_AdjFlow->GetSolution_Old(iPoint, 3));
    SetVolumeOutputValue("RES_ADJ_ENERGY",     iPoint, Node_AdjFlow->GetSolution(iPoint, 4) - Node_AdjFlow->GetSolution_Old(iPoint, 4));
  } else {
    SetVolumeOutputValue("RES_ADJ_ENERGY", iPoint, Node_AdjFlow->GetSolution(iPoint, 3) - Node_AdjFlow->GetSolution_Old(iPoint, 3));   
  }
  
  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    SetVolumeOutputValue("RES_TKE", iPoint, solver[TURB_SOL]->LinSysRes.GetBlock(iPoint, 0));
    SetVolumeOutputValue("RES_DISSIPATION", iPoint, solver[TURB_SOL]->LinSysRes.GetBlock(iPoint, 1));
    SetVolumeOutputValue("RES_ADJ_TKE",         iPoint, Node_AdjTurb->GetSolution(iPoint, 0) - Node_AdjTurb->GetSolution_Old(iPoint, 0));
    SetVolumeOutputValue("RES_ADJ_DISSIPATION", iPoint, Node_AdjTurb->GetSolution(iPoint, 1) - Node_AdjTurb->GetSolution_Old(iPoint, 1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("RES_NU_TILDE", iPoint, solver[TURB_SOL]->LinSysRes.GetBlock(iPoint, 0));
    SetVolumeOutputValue("RES_ADJ_NU_TILDE", iPoint, Node_AdjTurb->GetSolution(iPoint, 0) - Node_AdjTurb->GetSolution_Old(iPoint, 0));
    break;
  case NONE:
    break;
  }
  
  SetVolumeOutputValue("LIMITER_DENSITY",    iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 0));
  SetVolumeOutputValue("LIMITER_MOMENTUM-X", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 1));
  SetVolumeOutputValue("LIMITER_MOMENTUM-Y", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 2));
  if (nDim == 3){
    SetVolumeOutputValue("LIMITER_MOMENTUM-Z", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 3));
    SetVolumeOutputValue("LIMITER_ENERGY",     iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 4));
  } else {
    SetVolumeOutputValue("LIMITER_ENERGY", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 3));   
  }
  
  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    SetVolumeOutputValue("LIMITER_TKE",         iPoint, Node_Turb->GetLimiter_Primitive(iPoint, 0));
    SetVolumeOutputValue("LIMITER_DISSIPATION", iPoint, Node_Turb->GetLimiter_Primitive(iPoint, 1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("LIMITER_NU_TILDE", iPoint, Node_Turb->GetLimiter_Primitive(iPoint, 0));
    break;
  case NONE:
    break;
  }
  
  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
    SetVolumeOutputValue("DES_LENGTHSCALE", iPoint, Node_Flow->GetDES_LengthScale(iPoint));
    SetVolumeOutputValue("WALL_DISTANCE", iPoint, Node_Geo->GetWall_Distance());
  }
  
  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    SetVolumeOutputValue("ROE_DISSIPATION", iPoint, Node_Flow->GetRoe_Dissipation(iPoint));
  }
  
  if(config->GetKind_Solver() == RANS || config->GetKind_Solver() == NAVIER_STOKES){
    if (nDim == 3){
      SetVolumeOutputValue("VORTICITY_X", iPoint, Node_Flow->GetVorticity(iPoint)[0]);
      SetVolumeOutputValue("VORTICITY_Y", iPoint, Node_Flow->GetVorticity(iPoint)[1]);
      SetVolumeOutputValue("Q_CRITERION", iPoint, GetQ_Criterion(&(Node_Flow->GetGradient_Primitive(iPoint)[1])));
    } 
    SetVolumeOutputValue("VORTICITY_Z", iPoint, Node_Flow->GetVorticity(iPoint)[2]);      
  }

  if(config->GetError_Estimate() || config->GetKind_SU2() == SU2_MET) {
    if (nDim == 2){
      SetVolumeOutputValue("ANISO_METR0", iPoint, Node_Flow->GetAnisoMetr(iPoint, 0));
      SetVolumeOutputValue("ANISO_METR1", iPoint, Node_Flow->GetAnisoMetr(iPoint, 1));
      SetVolumeOutputValue("ANISO_METR2", iPoint, Node_Flow->GetAnisoMetr(iPoint, 2));
    }
    else{
      SetVolumeOutputValue("ANISO_METR0", iPoint, Node_Flow->GetAnisoMetr(iPoint, 0));
      SetVolumeOutputValue("ANISO_METR1", iPoint, Node_Flow->GetAnisoMetr(iPoint, 1));
      SetVolumeOutputValue("ANISO_METR2", iPoint, Node_Flow->GetAnisoMetr(iPoint, 2));
      SetVolumeOutputValue("ANISO_METR3", iPoint, Node_Flow->GetAnisoMetr(iPoint, 3));
      SetVolumeOutputValue("ANISO_METR4", iPoint, Node_Flow->GetAnisoMetr(iPoint, 4));
      SetVolumeOutputValue("ANISO_METR5", iPoint, Node_Flow->GetAnisoMetr(iPoint, 5));
    }
  }

  SetVolumeOutputValue("SENSITIVITY-X", iPoint, Node_AdjFlow->GetSensitivity(iPoint, 0));
  SetVolumeOutputValue("SENSITIVITY-Y", iPoint, Node_AdjFlow->GetSensitivity(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("SENSITIVITY-Z", iPoint, Node_AdjFlow->GetSensitivity(iPoint, 2));
  
  if (config->GetTime_Domain()){
    LoadTimeAveragedData(iPoint, Node_Flow);
  }
}

void COneShotFluidOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){
  
  if ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver()  == RANS)) {  
    SetVolumeOutputValue("SKIN_FRICTION-X", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 0));
    SetVolumeOutputValue("SKIN_FRICTION-Y", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 1));
    if (nDim == 3)
      SetVolumeOutputValue("SKIN_FRICTION-Z", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 2));
  
    SetVolumeOutputValue("HEAT_FLUX", iPoint, solver[FLOW_SOL]->GetHeatFlux(iMarker, iVertex));
    SetVolumeOutputValue("Y_PLUS", iPoint, solver[FLOW_SOL]->GetYPlus(iMarker, iVertex));
  }

  SetVolumeOutputValue("SENSITIVITY", iPoint, solver[ADJFLOW_SOL]->GetCSensitivity(iMarker, iVertex));

}

void COneShotFluidOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver)  {
  
  CSolver* flow_solver = solver[FLOW_SOL];
  CSolver* turb_solver = solver[TURB_SOL];
  CSolver* adjflow_solver = solver[ADJFLOW_SOL];
  CSolver* adjturb_solver = solver[ADJTURB_SOL];  

  SetHistoryOutputValue("SEARCH_ITER", adjflow_solver->GetArmijoIter());
  SetHistoryOutputValue("ACTIVE_DV", adjflow_solver->GetnActiveDV());
  SetHistoryOutputValue("ONE_SHOT_RHO", adjflow_solver->GetOneShotRho());
  SetHistoryOutputValue("ONE_SHOT_THETA", adjflow_solver->GetOneShotTheta());
  SetHistoryOutputValue("GRAD_NORM", adjflow_solver->GetShiftedLagGradNorm());
  for(unsigned short iConstr = 0; iConstr < config->GetnConstr(); iConstr++) {
    std::stringstream ss;
    ss << "LAMBDA_" << iConstr;
    std::string LambdaString = ss.str();
    SetHistoryOutputValue(LambdaString, adjflow_solver>GetLambdaValue(iConstr));
  }
  
  SetHistoryOutputValue("RMS_DENSITY", log10(flow_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_MOMENTUM-X", log10(flow_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_MOMENTUM-Y", log10(flow_solver->GetRes_RMS(2)));
  if (nDim == 2)
    SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(3)));
  else {
    SetHistoryOutputValue("RMS_MOMENTUM-Z", log10(flow_solver->GetRes_RMS(3)));
    SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(4)));
  }
  SetHistoryOutputValue("RMS_ADJ_DENSITY", log10(adjflow_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_ADJ_MOMENTUM-X", log10(adjflow_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_ADJ_MOMENTUM-Y", log10(adjflow_solver->GetRes_RMS(2)));
  if (geometry->GetnDim() == 3) {
    SetHistoryOutputValue("RMS_ADJ_MOMENTUM-Z", log10(adjflow_solver->GetRes_RMS(3)));
    SetHistoryOutputValue("RMS_ADJ_ENERGY", log10(adjflow_solver->GetRes_RMS(4)));
  } else {
    SetHistoryOutputValue("RMS_ADJ_ENERGY", log10(adjflow_solver->GetRes_RMS(3)));    
  }
  
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetHistoryOutputValue("RMS_NU_TILDE", log10(turb_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("RMS_ADJ_NU_TILDE", log10(adjturb_solver->GetRes_RMS(0)));
    break;  
  case SST: case SST_SUST:
    SetHistoryOutputValue("RMS_TKE", log10(turb_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("RMS_DISSIPATION",    log10(turb_solver->GetRes_RMS(1)));
    SetHistoryOutputValue("RMS_ADJ_TKE", log10(adjturb_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("RMS_ADJ_DISSIPATION",    log10(adjturb_solver->GetRes_RMS(1)));
    break;
  default: break;
  }
  
  SetHistoryOutputValue("MAX_DENSITY", log10(flow_solver->GetRes_Max(0)));
  SetHistoryOutputValue("MAX_MOMENTUM-X", log10(flow_solver->GetRes_Max(1)));
  SetHistoryOutputValue("MAX_MOMENTUM-Y", log10(flow_solver->GetRes_Max(2)));
  if (nDim == 2)
    SetHistoryOutputValue("MAX_ENERGY", log10(flow_solver->GetRes_Max(3)));
  else {
    SetHistoryOutputValue("MAX_MOMENTUM-Z", log10(flow_solver->GetRes_Max(3)));
    SetHistoryOutputValue("MAX_ENERGY", log10(flow_solver->GetRes_Max(4)));
  }
  SetHistoryOutputValue("MAX_ADJ_DENSITY", log10(adjflow_solver->GetRes_Max(0)));
  SetHistoryOutputValue("MAX_ADJ_MOMENTUM-X", log10(adjflow_solver->GetRes_Max(1)));
  SetHistoryOutputValue("MAX_ADJ_MOMENTUM-Y", log10(adjflow_solver->GetRes_Max(2)));
  if (geometry->GetnDim() == 3) {
    SetHistoryOutputValue("MAX_ADJ_MOMENTUM-Z", log10(adjflow_solver->GetRes_Max(3)));
    SetHistoryOutputValue("MAX_ADJ_ENERGY", log10(adjflow_solver->GetRes_Max(4)));
  } else {
    SetHistoryOutputValue("MAX_ADJ_ENERGY", log10(adjflow_solver->GetRes_Max(3)));    
  }
  
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetHistoryOutputValue("MAX_NU_TILDE", log10(turb_solver->GetRes_Max(0)));
    SetHistoryOutputValue("MAX_ADJ_NU_TILDE", log10(adjturb_solver->GetRes_Max(0)));
    break;  
  case SST: case SST_SUST:
    SetHistoryOutputValue("MAX_TKE", log10(turb_solver->GetRes_Max(0)));
    SetHistoryOutputValue("MAX_DISSIPATION",    log10(turb_solver->GetRes_Max(1)));
    SetHistoryOutputValue("MAX_ADJ_TKE", log10(adjturb_solver->GetRes_Max(0)));
    SetHistoryOutputValue("MAX_ADJ_DISSIPATION",    log10(adjturb_solver->GetRes_Max(1)));
    break;
  default: break;
  }
  
  SetHistoryOutputValue("HEATFLUX",     flow_solver->GetTotal_HeatFlux());
  SetHistoryOutputValue("HEATFLUX_MAX", flow_solver->GetTotal_MaxHeatFlux());
  SetHistoryOutputValue("TEMPERATURE",  flow_solver->GetTotal_AvgTemperature());
  
  SetHistoryOutputValue("CFL_NUMBER", config->GetCFL(MESH_0));
  
  SetHistoryOutputValue("LINSOL_ITER", flow_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(flow_solver->GetLinSol_Residual()));
  
  /*--- Set the analyse surface history values --- */
  
  SetAnalyzeSurface(flow_solver, geometry, config, false);
  
  /*--- Set aeroydnamic coefficients --- */
  
  SetAerodynamicCoefficients(config, flow_solver);

  /*--- Set Cp diff fields ---*/
  
  Set_CpInverseDesign(flow_solver, geometry, config);
  
  /*--- Set combo obj value --- */
  
  SetHistoryOutputValue("COMBO", flow_solver->GetTotal_ComboObj());

  SetHistoryOutputValue("SENS_GEO", adjflow_solver->GetTotal_Sens_Geo());
  SetHistoryOutputValue("SENS_AOA", adjflow_solver->GetTotal_Sens_AoA() * PI_NUMBER / 180.0);
  SetHistoryOutputValue("SENS_MACH", adjflow_solver->GetTotal_Sens_Mach());
  SetHistoryOutputValue("SENS_PRESS", adjflow_solver->GetTotal_Sens_Press());
  SetHistoryOutputValue("SENS_TEMP", adjflow_solver->GetTotal_Sens_Temp());
  
}

bool COneShotFluidOutput::SetInit_Residuals(CConfig *config){
  
  return (config->GetTime_Marching() != STEADY && (curInnerIter == 0))||
        (config->GetTime_Marching() == STEADY && (curInnerIter < 2));
  
}

bool COneShotFluidOutput::SetUpdate_Averages(CConfig *config){
  
  return (config->GetTime_Marching() != STEADY && (curInnerIter == config->GetnInner_Iter() - 1 || convergence));
      
}


