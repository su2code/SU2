/*!
 * \file output_flow_comp.cpp
 * \brief Main subroutines for compressible flow output
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

#include "../../include/output/CFlowCompOutput.hpp"

#include "../../../Common/include/geometry_structure.hpp"
#include "../../include/solver_structure.hpp"

CFlowCompOutput::CFlowCompOutput(CConfig *config, unsigned short nDim) : CFlowOutput(config, nDim, false) {
  
  turb_model = config->GetKind_Turb_Model();
  
  gridMovement = config->GetGrid_Movement(); 

  /*--- Set the default history fields if nothing is set in the config file ---*/
  
  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.push_back("ITER");
    requestedHistoryFields.push_back("RMS_RES");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }
  if (nRequestedScreenFields == 0){
    if (config->GetTime_Domain()) requestedScreenFields.push_back("TIME_ITER");
    if (multiZone) requestedScreenFields.push_back("OUTER_ITER");
    requestedScreenFields.push_back("INNER_ITER");
    requestedScreenFields.push_back("RMS_DENSITY");
    requestedScreenFields.push_back("RMS_MOMENTUM-X");
    requestedScreenFields.push_back("RMS_MOMENTUM-Y");
    requestedScreenFields.push_back("RMS_ENERGY");
    nRequestedScreenFields = requestedScreenFields.size();
  }
  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.push_back("COORDINATES");
    requestedVolumeFields.push_back("SOLUTION");
    requestedVolumeFields.push_back("PRIMITIVE");
    nRequestedVolumeFields = requestedVolumeFields.size();
  }
  
  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Comp. Fluid)";
  multiZoneHeaderString = ss.str();
  
  /*--- Set the volume filename --- */
  
  volumeFilename = config->GetVolume_FileName();
  
  /*--- Set the surface filename --- */
  
  surfaceFilename = config->GetSurfCoeff_FileName();
  
  /*--- Set the restart filename --- */
  
  restartFilename = config->GetRestart_FileName();


  /*--- Set the default convergence field --- */

  if (convField.size() == 0 ) convField = "RMS_DENSITY";

  
}

CFlowCompOutput::~CFlowCompOutput(void) {

  if (rank == MASTER_NODE){
    histFile.close();

  }


}



void CFlowCompOutput::SetHistoryOutputFields(CConfig *config){
  

  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the SOLUTION variables. 
  /// DESCRIPTION: Root-mean square residual of the density.
  AddHistoryOutput("RMS_DENSITY",    "rms[Rho]",  FORMAT_FIXED, "RMS_RES", "Root-mean square residual of the density.", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum x-component.
  AddHistoryOutput("RMS_MOMENTUM-X", "rms[RhoU]", FORMAT_FIXED, "RMS_RES", "Root-mean square residual of the momentum x-component.", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum y-component.
  AddHistoryOutput("RMS_MOMENTUM-Y", "rms[RhoV]", FORMAT_FIXED, "RMS_RES", "Root-mean square residual of the momentum y-component.", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum z-component.
  if (nDim == 3) AddHistoryOutput("RMS_MOMENTUM-Z", "rms[RhoW]", FORMAT_FIXED, "RMS_RES", "Root-mean square residual of the momentum z-component.", TYPE_RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the energy.
  AddHistoryOutput("RMS_ENERGY",     "rms[RhoE]", FORMAT_FIXED, "RMS_RES", "Root-mean square residual of the energy.", TYPE_RESIDUAL);
  
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Root-mean square residual of nu tilde (SA model).  
    AddHistoryOutput("RMS_NU_TILDE", "rms[nu]", FORMAT_FIXED, "RMS_RES", "Root-mean square residual of nu tilde (SA model).", TYPE_RESIDUAL);
    break;  
  case SST: case SST_SUST:
    /// DESCRIPTION: Root-mean square residual of kinetic energy (SST model).    
    AddHistoryOutput("RMS_TKE", "rms[k]",  FORMAT_FIXED, "RMS_RES", "Root-mean square residual of kinetic energy (SST model).", TYPE_RESIDUAL);
    /// DESCRIPTION: Root-mean square residual of the dissipation (SST model).    
    AddHistoryOutput("RMS_DISSIPATION", "rms[w]",  FORMAT_FIXED, "RMS_RES", "Root-mean square residual of dissipation (SST model).", TYPE_RESIDUAL);
    break;
  default: break;
  }
  /// END_GROUP
   
  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the SOLUTION variables. 
  /// DESCRIPTION: Maximum residual of the density.
  AddHistoryOutput("MAX_DENSITY",    "max[Rho]",  FORMAT_FIXED,   "MAX_RES", "Maximum square residual of the density.", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum x-component. 
  AddHistoryOutput("MAX_MOMENTUM-X", "max[RhoU]", FORMAT_FIXED,   "MAX_RES", "Maximum square residual of the momentum x-component.", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum y-component. 
  AddHistoryOutput("MAX_MOMENTUM-Y", "max[RhoV]", FORMAT_FIXED,   "MAX_RES", "Maximum square residual of the momentum y-component.", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum z-component. 
  if (nDim == 3) AddHistoryOutput("MAX_MOMENTUM-Z", "max[RhoW]", FORMAT_FIXED,"MAX_RES", "Maximum residual of the z-component.", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the energy.  
  AddHistoryOutput("MAX_ENERGY",     "max[RhoE]", FORMAT_FIXED,   "MAX_RES", "Maximum residual of the energy.", TYPE_RESIDUAL);
  
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Maximum residual of nu tilde (SA model).
    AddHistoryOutput("MAX_NU_TILDE",       "max[nu]", FORMAT_FIXED, "MAX_RES", "Maximum residual of nu tilde (SA model).", TYPE_RESIDUAL);
    break;  
  case SST: case SST_SUST:
    /// DESCRIPTION: Maximum residual of kinetic energy (SST model). 
    AddHistoryOutput("MAX_TKE", "max[k]",  FORMAT_FIXED, "MAX_RES", "Maximum residual of kinetic energy (SST model).", TYPE_RESIDUAL);
    /// DESCRIPTION: Maximum residual of the dissipation (SST model).   
    AddHistoryOutput("MAX_DISSIPATION",    "max[w]",  FORMAT_FIXED, "MAX_RES", "Maximum residual of dissipation (SST model).", TYPE_RESIDUAL); 
    break;
  default: break;
  }
  /// END_GROUP
  
  /// BEGIN_GROUP: BGS_RES, DESCRIPTION: The block Gauss Seidel residuals of the SOLUTION variables. 
  /// DESCRIPTION: Maximum residual of the density.
  AddHistoryOutput("BGS_DENSITY",    "bgs[Rho]",  FORMAT_FIXED,   "BGS_RES", "BGS residual of the density.", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum x-component. 
  AddHistoryOutput("BGS_MOMENTUM-X", "bgs[RhoU]", FORMAT_FIXED,   "BGS_RES", "BGS residual of the momentum x-component.", TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum y-component. 
  AddHistoryOutput("BGS_MOMENTUM-Y", "bgs[RhoV]", FORMAT_FIXED,   "BGS_RES", "BGS residual of the momentum y-component.",  TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum z-component. 
  if (nDim == 3) AddHistoryOutput("BGS_MOMENTUM-Z", "bgs[RhoW]", FORMAT_FIXED, "BGS_RES", "BGS residual of the z-component.",  TYPE_RESIDUAL);
  /// DESCRIPTION: Maximum residual of the energy.  
  AddHistoryOutput("BGS_ENERGY",     "bgs[RhoE]", FORMAT_FIXED,   "BGS_RES", "BGS residual of the energy.",  TYPE_RESIDUAL);
  
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    /// DESCRIPTION: Maximum residual of nu tilde (SA model).
    AddHistoryOutput("BGS_NU_TILDE",       "bgs[nu]", FORMAT_FIXED, "BGS_RES", "BGS residual of nu tilde (SA model).",  TYPE_RESIDUAL);
    break;  
  case SST: case SST_SUST:
    /// DESCRIPTION: Maximum residual of kinetic energy (SST model). 
    AddHistoryOutput("BGS_TKE", "bgs[k]",  FORMAT_FIXED, "BGS_RES", "BGS residual of kinetic energy (SST model).",  TYPE_RESIDUAL);
    /// DESCRIPTION: Maximum residual of the dissipation (SST model).   
    AddHistoryOutput("BGS_DISSIPATION",    "bgs[w]",  FORMAT_FIXED, "BGS_RES", "BGS residual of dissipation (SST model).", TYPE_RESIDUAL); 
    break;
  default: break;
  }
  /// END_GROUP

  vector<string> Marker_Monitoring;
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++){
    Marker_Monitoring.push_back(config->GetMarker_Monitoring_TagBound(iMarker_Monitoring));
  }  
  /// BEGIN_GROUP: AEROELASTIC, DESCRIPTION: Aeroelastic plunge, pitch
  /// DESCRIPTION: Aeroelastic plunge
  AddHistoryOutputPerSurface("PLUNGE", "plunge", FORMAT_FIXED, "AEROELASTIC", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Aeroelastic pitch
  AddHistoryOutputPerSurface("PITCH",  "pitch",  FORMAT_FIXED, "AEROELASTIC", Marker_Monitoring, TYPE_COEFFICIENT);
  /// END_GROUP
   

  /// DESCRIPTION: Linear solver iterations   
  AddHistoryOutput("LINSOL_ITER", "Linear_Solver_Iterations", FORMAT_INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", FORMAT_FIXED, "LINSOL", "Residual of the linear solver.");
 
  /// BEGIN_GROUP: ENGINE_OUTPUT, DESCRIPTION: Engine output
  /// DESCRIPTION: Aero CD drag
  AddHistoryOutput("AEROCDRAG",                  "AeroCDrag",                  FORMAT_SCIENTIFIC, "ENGINE_OUTPUT", "Aero CD drag", TYPE_COEFFICIENT);
  /// DESCRIPTION: Solid CD drag  
  AddHistoryOutput("SOLIDCDRAG",                 "SolidCDrag",                 FORMAT_SCIENTIFIC, "ENGINE_OUTPUT", "Solid CD drag ", TYPE_COEFFICIENT);
  /// DESCRIPTION: Radial distortion 
  AddHistoryOutput("RADIAL_DISTORTION",          "Radial_Distortion",          FORMAT_SCIENTIFIC, "ENGINE_OUTPUT", "Radial distortion ", TYPE_COEFFICIENT);
  /// DESCRIPTION: Circumferential distortion
  AddHistoryOutput("CIRCUMFERENTIAL_DISTORTION", "Circumferential_Distortion", FORMAT_SCIENTIFIC, "ENGINE_OUTPUT", "Circumferential distortion", TYPE_COEFFICIENT);
  /// END_GROUP
  
  /// BEGIN_GROUP: ROTATING_FRAME, DESCRIPTION: Coefficients related to a rotating frame of reference.
  /// DESCRIPTION: Merit  
  AddHistoryOutput("MERIT", "CMerit", FORMAT_SCIENTIFIC, "ROTATING_FRAME", "Merit", TYPE_COEFFICIENT);
  /// DESCRIPTION: CT 
  AddHistoryOutput("CT",    "CT",     FORMAT_SCIENTIFIC, "ROTATING_FRAME", "CT", TYPE_COEFFICIENT);
  /// DESCRIPTION: CQ  
  AddHistoryOutput("CQ",    "CQ",     FORMAT_SCIENTIFIC, "ROTATING_FRAME", "CQ", TYPE_COEFFICIENT);
  /// END_GROUP
  
  /// BEGIN_GROUP: EQUIVALENT_AREA, DESCRIPTION: Equivalent area.  
  /// DESCRIPTION: Equivalent area    
  AddHistoryOutput("EQUIV_AREA",   "CEquiv_Area",  FORMAT_SCIENTIFIC, "EQUIVALENT_AREA", "Equivalent area", TYPE_COEFFICIENT);
  /// DESCRIPTION: Nearfield obj. function      
  AddHistoryOutput("NEARFIELD_OF", "CNearFieldOF", FORMAT_SCIENTIFIC, "EQUIVALENT_AREA", "Nearfield obj. function ", TYPE_COEFFICIENT);
  /// END_GROUP

  ///   /// BEGIN_GROUP: HEAT_COEFF, DESCRIPTION: Heat coefficients on all surfaces set with MARKER_MONITORING.
  /// DESCRIPTION: Total heatflux
  AddHistoryOutput("HEATFLUX", "HF",      FORMAT_SCIENTIFIC, "HEAT", "Total heatflux on all surfaces set with MARKER_MONITORING.", TYPE_COEFFICIENT);
  /// DESCRIPTION: Maximal heatflux  
  AddHistoryOutput("HEATFLUX_MAX", "maxHF",    FORMAT_SCIENTIFIC, "HEAT", "Total maximum heatflux on all surfaces set with MARKER_MONITORING.", TYPE_COEFFICIENT);
  /// DESCRIPTION: Temperature
  AddHistoryOutput("TEMPERATURE", "Temp", FORMAT_SCIENTIFIC, "HEAT",  "Total avg. temperature on all surfaces set with MARKER_MONITORING.", TYPE_COEFFICIENT);
  /// END_GROUP
  
  AddHistoryOutput("CFL_NUMBER", "CFL number", FORMAT_SCIENTIFIC, "CFL_NUMBER", "Current value of the CFL number");
  
  if (config->GetDeform_Mesh()){
    AddHistoryOutput("DEFORM_MIN_VOLUME", "MinVolume", FORMAT_SCIENTIFIC, "DEFORM", "Minimum volume in the mesh");
    AddHistoryOutput("DEFORM_MAX_VOLUME", "MaxVolume", FORMAT_SCIENTIFIC, "DEFORM", "Maximum volume in the mesh");
    AddHistoryOutput("DEFORM_ITER", "DeformIter", FORMAT_INTEGER, "DEFORM", "Linear solver iterations for the mesh deformation");
    AddHistoryOutput("DEFORM_RESIDUAL", "DeformRes", FORMAT_FIXED, "DEFORM", "Residual of the linear solver for the mesh deformation");    
  }
  
  /*--- Add analyze surface history fields --- */
  
  AddAnalyzeSurfaceOutput(config);
  
  /*--- Add aerodynamic coefficients fields --- */
  
  AddAerodynamicCoefficients(config); 
  
  /*--- Add Cp diff fields ---*/
  
  Add_CpInverseDesignOutput(config);
  
  /*--- Add combo obj value --- */
  
  AddHistoryOutput("COMBO", "ComboObj", FORMAT_SCIENTIFIC, "COMBO", "Combined obj. function value.", TYPE_COEFFICIENT);
}

void CFlowCompOutput::SetVolumeOutputFields(CConfig *config){
  
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
  
  // Turbulent Residuals
  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    AddVolumeOutput("TKE", "Turb_Kin_Energy", "SOLUTION", "Turbulent kinetic energy");
    AddVolumeOutput("DISSIPATION", "Omega", "SOLUTION", "Rate of dissipation");
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    AddVolumeOutput("NU_TILDE", "Nu_Tilde", "SOLUTION", "Spalart-Allmaras variable");
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
  
  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    AddVolumeOutput("RES_TKE", "Residual_TKE", "RESIDUAL", "Residual of turbulent kinetic energy");
    AddVolumeOutput("RES_DISSIPATION", "Residual_Omega", "RESIDUAL", "Residual of the rate of dissipation");
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    AddVolumeOutput("RES_NU_TILDE", "Residual_Nu_Tilde", "RESIDUAL", "Residual of the Spalart-Allmaras variable");
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
  
  if (config->GetTime_Domain()){
    AddVolumeOutput("MEAN_DENSITY", "MeanDensity", "TAVG_SOLUTION", "Time-averaged value of the density");
    AddVolumeOutput("MEAN_VELOCITY-X", "MeanVelocity_x", "TAVG_SOLUTION", "Time-averaged value of the density");
    AddVolumeOutput("MEAN_VELOCITY-Y", "MeanVelocity_y", "TAVG_SOLUTION", "Time-averaged value of the density");
    if (nDim == 3)
      AddVolumeOutput("MEAN_VELOCITY-Z", "MeanVelocity_z", "TAVG_SOLUTION", "Time-averaged value of the density");
    if (nDim == 3){
      AddVolumeOutput("MEAN_ENERGY", "MeanEnergy", "TAVG_SOLUTION", "Time-averaged value of the density");
    }
    else {
      AddVolumeOutput("MEAN_ENERGY", "MeanEnergy", "TAVG_SOLUTION", "Time-averaged value of the density");    
    }
    AddVolumeOutput("MEAN_PRESSURE", "MeanPressure", "TAVG_SOLUTION", "Time-averaged value of the density");
    if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
      AddVolumeOutput("MEAN_ROE_DISSIPATION", "MeanRoe_Dissipation", "TAVG_SOLUTION", "Time-averaged value of the density");
    }  
    
    AddVolumeOutput("RMS_U",   "RMSVelocity_x", "TAVG_SOLUTION", "Time-averaged value of the density");
    AddVolumeOutput("RMS_V",   "RMSVelocity_y", "TAVG_SOLUTION", "Time-averaged value of the density");
    AddVolumeOutput("RMS_UV",  "RMSUV",         "TAVG_SOLUTION", "Time-averaged value of the density");    
    AddVolumeOutput("RMS_P",   "RMSPressure",   "TAVG_SOLUTION", "Time-average value pf pressure");
    AddVolumeOutput("UUPRIME", "UUPrimeMean", "TAVG_SOLUTION", "Time-averaged value of the density");
    AddVolumeOutput("VVPRIME", "VVPrimeMean", "TAVG_SOLUTION", "Time-averaged value of the density");
    AddVolumeOutput("UVPRIME", "UVPrimeMean", "TAVG_SOLUTION", "Time-averaged value of the density");
    AddVolumeOutput("PPRIME",  "PPrimeMean",   "TAVG_SOLUTION", "Time-averaged value of the density");
    if (nDim == 3){
      AddVolumeOutput("RMS_W",   "RMSVelocity_z", "TAVG_SOLUTION", "Time-averaged value of the density");
      AddVolumeOutput("RMS_UW", "RMSUW", "TAVG_SOLUTION", "Time-averaged value of the density");
      AddVolumeOutput("RMS_VW", "RMSVW", "TAVG_SOLUTION", "Time-averaged value of the density");
      AddVolumeOutput("WWPRIME", "WWPrimeMean", "TAVG_SOLUTION", "Time-averaged value of the density");
      AddVolumeOutput("UWPRIME", "UWPrimeMean", "TAVG_SOLUTION", "Time-averaged value of the density");
      AddVolumeOutput("VWPRIME", "VWPrimeMean", "TAVG_SOLUTION", "Time-averaged value of the density");
    }
  }
}

void CFlowCompOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){
  
  CVariable* Node_Flow = solver[FLOW_SOL]->node[iPoint]; 
  CVariable* Node_Turb = NULL;
  
  if (config->GetKind_Turb_Model() != NONE){
    Node_Turb = solver[TURB_SOL]->node[iPoint]; 
  }
  
  CPoint*    Node_Geo  = geometry->node[iPoint];
          
  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));  
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));
  
  SetVolumeOutputValue("DENSITY",    iPoint, Node_Flow->GetSolution(0));
  SetVolumeOutputValue("MOMENTUM-X", iPoint, Node_Flow->GetSolution(1));
  SetVolumeOutputValue("MOMENTUM-Y", iPoint, Node_Flow->GetSolution(2));
  if (nDim == 3){
    SetVolumeOutputValue("MOMENTUM-Z", iPoint, Node_Flow->GetSolution(3));
    SetVolumeOutputValue("ENERGY",     iPoint, Node_Flow->GetSolution(4));
  } else {
    SetVolumeOutputValue("ENERGY",     iPoint, Node_Flow->GetSolution(3));    
  }
  
  // Turbulent Residuals
  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    SetVolumeOutputValue("TKE", iPoint, Node_Turb->GetSolution(0));
    SetVolumeOutputValue("DISSIPATION", iPoint, Node_Turb->GetSolution(1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("NU_TILDE", iPoint, Node_Turb->GetSolution(0));
    break;
  case NONE:
    break;
  }
  
  SetVolumeOutputValue("PRESSURE", iPoint, Node_Flow->GetPressure());
  SetVolumeOutputValue("TEMPERATURE", iPoint, Node_Flow->GetTemperature());
  SetVolumeOutputValue("MACH", iPoint, sqrt(Node_Flow->GetVelocity2())/Node_Flow->GetSoundSpeed());
  su2double VelMag = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++){
    VelMag += solver[FLOW_SOL]->GetVelocity_Inf(iDim)*solver[FLOW_SOL]->GetVelocity_Inf(iDim);    
  }
  su2double factor = 1.0/(0.5*solver[FLOW_SOL]->GetDensity_Inf()*VelMag); 
  SetVolumeOutputValue("PRESSURE_COEFF", iPoint, (Node_Flow->GetPressure() - solver[FLOW_SOL]->GetPressure_Inf())*factor);
  
  if (config->GetKind_Solver() == RANS || config->GetKind_Solver() == NAVIER_STOKES){
    SetVolumeOutputValue("LAMINAR_VISCOSITY", iPoint, Node_Flow->GetLaminarViscosity());
  }
  
  if (config->GetKind_Solver() == RANS) {
    SetVolumeOutputValue("EDDY_VISCOSITY", iPoint, Node_Flow->GetEddyViscosity());
  }
  
  if (config->GetKind_Trans_Model() == BC){
    SetVolumeOutputValue("INTERMITTENCY", iPoint, Node_Turb->GetGammaBC());
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
  
  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    SetVolumeOutputValue("RES_TKE", iPoint, solver[TURB_SOL]->LinSysRes.GetBlock(iPoint, 0));
    SetVolumeOutputValue("RES_DISSIPATION", iPoint, solver[TURB_SOL]->LinSysRes.GetBlock(iPoint, 1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("RES_NU_TILDE", iPoint, solver[TURB_SOL]->LinSysRes.GetBlock(iPoint, 0));
    break;
  case NONE:
    break;
  }
  
  SetVolumeOutputValue("LIMITER_DENSITY", iPoint, Node_Flow->GetLimiter_Primitive(0));
  SetVolumeOutputValue("LIMITER_MOMENTUM-X", iPoint, Node_Flow->GetLimiter_Primitive(1));
  SetVolumeOutputValue("LIMITER_MOMENTUM-Y", iPoint, Node_Flow->GetLimiter_Primitive(2));
  if (nDim == 3){
    SetVolumeOutputValue("LIMITER_MOMENTUM-Z", iPoint, Node_Flow->GetLimiter_Primitive(3));
    SetVolumeOutputValue("LIMITER_ENERGY", iPoint, Node_Flow->GetLimiter_Primitive(4));
  } else {
    SetVolumeOutputValue("LIMITER_ENERGY", iPoint, Node_Flow->GetLimiter_Primitive(3));   
  }
  
  switch(config->GetKind_Turb_Model()){
  case SST: case SST_SUST:
    SetVolumeOutputValue("LIMITER_TKE", iPoint, Node_Turb->GetLimiter_Primitive(0));
    SetVolumeOutputValue("LIMITER_DISSIPATION", iPoint, Node_Turb->GetLimiter_Primitive(1));
    break;
  case SA: case SA_COMP: case SA_E: 
  case SA_E_COMP: case SA_NEG: 
    SetVolumeOutputValue("LIMITER_NU_TILDE", iPoint, Node_Turb->GetLimiter_Primitive(0));
    break;
  case NONE:
    break;
  }
  
  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
    SetVolumeOutputValue("DES_LENGTHSCALE", iPoint, Node_Flow->GetDES_LengthScale());
    SetVolumeOutputValue("WALL_DISTANCE", iPoint, Node_Geo->GetWall_Distance());
  }
  
  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    SetVolumeOutputValue("ROE_DISSIPATION", iPoint, Node_Flow->GetRoe_Dissipation());
  }
  
  if(config->GetKind_Solver() == RANS || config->GetKind_Solver() == NAVIER_STOKES){
    if (nDim == 3){
      SetVolumeOutputValue("VORTICITY_X", iPoint, Node_Flow->GetVorticity()[0]);
      SetVolumeOutputValue("VORTICITY_Y", iPoint, Node_Flow->GetVorticity()[1]);
      SetVolumeOutputValue("Q_CRITERION", iPoint, GetQ_Criterion(config, geometry, Node_Flow));            
    } 
    SetVolumeOutputValue("VORTICITY_Z", iPoint, Node_Flow->GetVorticity()[2]);      
  }
  
  if (config->GetTime_Domain()){
    SetAvgVolumeOutputValue("MEAN_DENSITY", iPoint, Node_Flow->GetSolution(0));
    SetAvgVolumeOutputValue("MEAN_VELOCITY-X", iPoint, Node_Flow->GetSolution(1)/Node_Flow->GetSolution(0));
    SetAvgVolumeOutputValue("MEAN_VELOCITY-Y", iPoint, Node_Flow->GetSolution(2)/Node_Flow->GetSolution(0));
    if (nDim == 3)
      SetAvgVolumeOutputValue("MEAN_VELOCITY-Z", iPoint, Node_Flow->GetSolution(3)/Node_Flow->GetSolution(0));
    if (nDim == 3){
      SetAvgVolumeOutputValue("MEAN_ENERGY", iPoint, Node_Flow->GetSolution(4));
    } else {
      SetAvgVolumeOutputValue("MEAN_ENERGY", iPoint, Node_Flow->GetSolution(3));    
    }
    SetAvgVolumeOutputValue("MEAN_PRESSURE", iPoint, Node_Flow->GetPressure());
    if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
      SetAvgVolumeOutputValue("MEAN_ROE_DISSIPATION", iPoint, Node_Flow->GetRoe_Dissipation());
    }  
    
    SetAvgVolumeOutputValue("RMS_U", iPoint, pow(Node_Flow->GetSolution(1)/Node_Flow->GetSolution(0),2));
    SetAvgVolumeOutputValue("RMS_V", iPoint, pow(Node_Flow->GetSolution(2)/Node_Flow->GetSolution(0),2));
    SetAvgVolumeOutputValue("RMS_UV", iPoint, (Node_Flow->GetSolution(1)/Node_Flow->GetSolution(0)) * (Node_Flow->GetSolution(2)/Node_Flow->GetSolution(0)));
    SetAvgVolumeOutputValue("RMS_P", iPoint, pow(Node_Flow->GetPressure(),2));
    if (nDim == 3){
      SetAvgVolumeOutputValue("RMS_W", iPoint, pow(Node_Flow->GetSolution(3)/Node_Flow->GetSolution(0),2));
      SetAvgVolumeOutputValue("RMS_VW", iPoint, (Node_Flow->GetSolution(2)/Node_Flow->GetSolution(0)) * (Node_Flow->GetSolution(3)/Node_Flow->GetSolution(0)));
      SetAvgVolumeOutputValue("RMS_UW", iPoint, (Node_Flow->GetSolution(1)/Node_Flow->GetSolution(0)) * (Node_Flow->GetSolution(3)/Node_Flow->GetSolution(0)));
    }
    
    const su2double umean  = GetVolumeOutputValue("MEAN_VELOCITY-X", iPoint);
    const su2double uumean = GetVolumeOutputValue("RMS_U", iPoint);    
    const su2double vmean  = GetVolumeOutputValue("MEAN_VELOCITY-Y", iPoint);
    const su2double vvmean = GetVolumeOutputValue("RMS_V", iPoint);    
    const su2double uvmean = GetVolumeOutputValue("RMS_UV", iPoint);
    const su2double pmean  = GetVolumeOutputValue("MEAN_PRESSURE", iPoint);
    const su2double ppmean = GetVolumeOutputValue("RMS_P", iPoint);
    
    SetVolumeOutputValue("UUPRIME", iPoint, -(umean*umean - uumean));
    SetVolumeOutputValue("VVPRIME", iPoint, -(vmean*vmean - vvmean));
    SetVolumeOutputValue("UVPRIME", iPoint, -(umean*vmean - uvmean));
    SetVolumeOutputValue("PPRIME",  iPoint, -(pmean*pmean - ppmean));
    if (nDim == 3){
      const su2double wmean  = GetVolumeOutputValue("MEAN_VELOCITY-Z", iPoint);
      const su2double wwmean = GetVolumeOutputValue("RMS_W", iPoint);
      const su2double uwmean = GetVolumeOutputValue("RMS_UW", iPoint);
      const su2double vwmean = GetVolumeOutputValue("RMS_VW", iPoint);
      SetVolumeOutputValue("WWPRIME", iPoint, -(wmean*wmean - wwmean));
      SetVolumeOutputValue("UWPRIME", iPoint, -(umean*wmean - uwmean));
      SetVolumeOutputValue("VWPRIME",  iPoint, -(vmean*wmean - vwmean));
    } 
  }
}

void CFlowCompOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){
  
  if ((config->GetKind_Solver() == NAVIER_STOKES) || (config->GetKind_Solver()  == RANS)) {  
    SetVolumeOutputValue("SKIN_FRICTION-X", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 0));
    SetVolumeOutputValue("SKIN_FRICTION-Y", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 1));
    if (nDim == 3)
      SetVolumeOutputValue("SKIN_FRICTION-Z", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 2));
  
    SetVolumeOutputValue("HEAT_FLUX", iPoint, solver[FLOW_SOL]->GetHeatFlux(iMarker, iVertex));
    SetVolumeOutputValue("Y_PLUS", iPoint, solver[FLOW_SOL]->GetYPlus(iMarker, iVertex));
  }
}

void CFlowCompOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver)  {
  
  CSolver* flow_solver = solver[FLOW_SOL];
  CSolver* turb_solver = solver[TURB_SOL];
  CSolver* mesh_solver = solver[MESH_SOL];
  
  SetHistoryOutputValue("RMS_DENSITY", log10(flow_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_MOMENTUM-X", log10(flow_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_MOMENTUM-Y", log10(flow_solver->GetRes_RMS(2)));
  if (nDim == 2)
    SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(3)));
  else {
    SetHistoryOutputValue("RMS_MOMENTUM-Z", log10(flow_solver->GetRes_RMS(3)));
    SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(4)));
  }
  
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetHistoryOutputValue("RMS_NU_TILDE", log10(turb_solver->GetRes_RMS(0)));
    break;  
  case SST: case SST_SUST:
    SetHistoryOutputValue("RMS_TKE", log10(turb_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("RMS_DISSIPATION",    log10(turb_solver->GetRes_RMS(1)));
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
  
  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    SetHistoryOutputValue("MAX_NU_TILDE", log10(turb_solver->GetRes_Max(0)));
    break;  
  case SST: case SST_SUST:
    SetHistoryOutputValue("MAX_TKE", log10(turb_solver->GetRes_Max(0)));
    SetHistoryOutputValue("MAX_DISSIPATION",    log10(turb_solver->GetRes_Max(1)));
    break;
  default: break;
  }
  
  if (multiZone){
    SetHistoryOutputValue("BGS_DENSITY", log10(flow_solver->GetRes_BGS(0)));
    SetHistoryOutputValue("BGS_MOMENTUM-X", log10(flow_solver->GetRes_BGS(1)));
    SetHistoryOutputValue("BGS_MOMENTUM-Y", log10(flow_solver->GetRes_BGS(2)));
    if (nDim == 2)
      SetHistoryOutputValue("BGS_ENERGY", log10(flow_solver->GetRes_BGS(3)));
    else {
      SetHistoryOutputValue("BGS_MOMENTUM-Z", log10(flow_solver->GetRes_BGS(3)));
      SetHistoryOutputValue("BGS_ENERGY", log10(flow_solver->GetRes_BGS(4)));
    }
    
    
    switch(turb_model){
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      SetHistoryOutputValue("BGS_NU_TILDE", log10(turb_solver->GetRes_BGS(0)));
      break;  
    case SST:
      SetHistoryOutputValue("BGS_TKE", log10(turb_solver->GetRes_BGS(0)));
      SetHistoryOutputValue("BGS_DISSIPATION",    log10(turb_solver->GetRes_BGS(1)));
      break;
    default: break;
    }
  }
  
  SetHistoryOutputValue("HEATFLUX",     flow_solver->GetTotal_HeatFlux());
  SetHistoryOutputValue("HEATFLUX_MAX", flow_solver->GetTotal_MaxHeatFlux());
  SetHistoryOutputValue("TEMPERATURE",  flow_solver->GetTotal_AvgTemperature());
  
  SetHistoryOutputValue("CFL_NUMBER", config->GetCFL(MESH_0));
  
  SetHistoryOutputValue("LINSOL_ITER", flow_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(flow_solver->GetLinSol_Residual()));
  
  if (config->GetDeform_Mesh()){
    SetHistoryOutputValue("DEFORM_MIN_VOLUME", mesh_solver->GetMinimum_Volume());
    SetHistoryOutputValue("DEFORM_MAX_VOLUME", mesh_solver->GetMaximum_Volume());
    SetHistoryOutputValue("DEFORM_ITER", mesh_solver->GetIterLinSolver());
    SetHistoryOutputValue("DEFORM_RESIDUAL", log10(mesh_solver->GetLinSol_Residual()));    
  }
  
  /*--- Set the analyse surface history values --- */
  
  SetAnalyzeSurface(flow_solver, geometry, config, false);
  
  /*--- Set aeroydnamic coefficients --- */
  
  SetAerodynamicCoefficients(config, flow_solver);

  /*--- Set Cp diff fields ---*/
  
  Set_CpInverseDesign(flow_solver, geometry, config);
  
  /*--- Set combo obj value --- */
  
  SetHistoryOutputValue("COMBO", flow_solver->GetTotal_ComboObj());
  
}

su2double CFlowCompOutput::GetQ_Criterion(CConfig *config, CGeometry *geometry, CVariable* node_flow){
  
  unsigned short iDim;
  su2double Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  
  for (iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0 ; jDim < nDim; jDim++) {
      Grad_Vel[iDim][jDim] = node_flow->GetGradient_Primitive(iDim+1, jDim);
    }
  }
  
  su2double s11 = Grad_Vel[0][0];
  su2double s12 = 0.5 * (Grad_Vel[0][1] + Grad_Vel[1][0]);
  su2double s13 = 0.5 * (Grad_Vel[0][2] + Grad_Vel[2][0]);
  su2double s22 = Grad_Vel[1][1];
  su2double s23 = 0.5 * (Grad_Vel[1][2] + Grad_Vel[2][1]);
  su2double s33 = Grad_Vel[2][2];
  su2double omega12 = 0.5 * (Grad_Vel[0][1] - Grad_Vel[1][0]);
  su2double omega13 = 0.5 * (Grad_Vel[0][2] - Grad_Vel[2][0]);
  su2double omega23 = 0.5 * (Grad_Vel[1][2] - Grad_Vel[2][1]);
  
  su2double Q = 2. * pow( omega12, 2.) + 2. * pow( omega13, 2.) + 2. * pow( omega23, 2.) - \
      pow( s11, 2.) - pow( s22, 2.) - pow( s33, 2.0) - 2. * pow( s12, 2.) - 2. * pow( s13, 2.) - 2. * pow( s23, 2.0);
  
  return Q;
}


bool CFlowCompOutput::SetInit_Residuals(CConfig *config){
  
  return (config->GetTime_Marching() != STEADY && (curInnerIter == 0))||
        (config->GetTime_Marching() == STEADY && (curInnerIter < 2));
  
}

bool CFlowCompOutput::SetUpdate_Averages(CConfig *config){
  
  return (config->GetTime_Marching() != STEADY && (curInnerIter == config->GetnInner_Iter() - 1 || convergence));
      
}


