/*!
 * \file CFlowCompFEMOutput.cpp
 * \brief Main subroutines for compressible flow output
 * \author R. Sanchez
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/output/CFlowCompFEMOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CFlowCompFEMOutput::CFlowCompFEMOutput(CConfig *config, unsigned short nDim) : CFlowOutput(config, nDim, true) {

  turb_model = config->GetKind_Turb_Model();

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    requestedHistoryFields.emplace_back("ITER");
    requestedHistoryFields.emplace_back("RMS_RES");
    nRequestedHistoryFields = requestedHistoryFields.size();
  }
  if (nRequestedScreenFields == 0){
    if (config->GetTime_Domain()) requestedScreenFields.emplace_back("TIME_ITER");
    if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("RMS_DENSITY");
    requestedScreenFields.emplace_back("RMS_MOMENTUM-X");
    requestedScreenFields.emplace_back("RMS_MOMENTUM-Y");
    requestedScreenFields.emplace_back("RMS_ENERGY");
    nRequestedScreenFields = requestedScreenFields.size();
  }
  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    requestedVolumeFields.emplace_back("PRIMITIVE");
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

  if (convFields.empty() ) convFields.emplace_back("RMS_DENSITY");

}

CFlowCompFEMOutput::~CFlowCompFEMOutput() = default;



void CFlowCompFEMOutput::SetHistoryOutputFields(CConfig *config){

  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the SOLUTION variables.
  /// DESCRIPTION: Root-mean square residual of the density.
  AddHistoryOutput("RMS_DENSITY",    "rms[Rho]",  ScreenOutputFormat::FIXED,   "RMS_RES", "Root-mean square residual of the density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum x-component.
  AddHistoryOutput("RMS_MOMENTUM-X", "rms[RhoU]", ScreenOutputFormat::FIXED,   "RMS_RES", "Root-mean square residual of the momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum y-component.
  AddHistoryOutput("RMS_MOMENTUM-Y", "rms[RhoV]", ScreenOutputFormat::FIXED,   "RMS_RES", "Root-mean square residual of the momentum y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum z-component.
  if (nDim == 3) AddHistoryOutput("RMS_MOMENTUM-Z", "rms[RhoW]", ScreenOutputFormat::FIXED,   "RMS_RES",  "Root-mean square residual of the momentum z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the energy.
  AddHistoryOutput("RMS_ENERGY",     "rms[RhoE]", ScreenOutputFormat::FIXED,   "RMS_RES", "Root-mean square residual of the energy.", HistoryFieldType::RESIDUAL);
  /// END_GROUP

  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the SOLUTION variables.
  /// DESCRIPTION: Maximum residual of the density.
  AddHistoryOutput("MAX_DENSITY",    "max[Rho]",  ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum square residual of the density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum x-component.
  AddHistoryOutput("MAX_MOMENTUM-X", "max[RhoU]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum square residual of the momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum y-component.
  AddHistoryOutput("MAX_MOMENTUM-Y", "max[RhoV]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum square residual of the momentum y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum z-component.
  if (nDim == 3) AddHistoryOutput("MAX_MOMENTUM-Z", "max[RhoW]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the energy.
  AddHistoryOutput("MAX_ENERGY",     "max[RhoE]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum residual of the energy.", HistoryFieldType::RESIDUAL);
  /// END_GROUP

  AddHistoryOutput("CFL_NUMBER", "CFL number", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current value of the CFL number");

  /*--- Add analyze surface history fields --- */

  AddAnalyzeSurfaceOutput(config);

  /*--- Add aerodynamic coefficients fields --- */

  AddAerodynamicCoefficients(config);

}

void CFlowCompFEMOutput::SetVolumeOutputFields(CConfig *config){

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

  // Primitive variables
  AddVolumeOutput("PRESSURE",    "Pressure",                "PRIMITIVE", "Pressure");
  AddVolumeOutput("TEMPERATURE", "Temperature",             "PRIMITIVE", "Temperature");
  AddVolumeOutput("MACH",        "Mach",                    "PRIMITIVE", "Mach number");
  AddVolumeOutput("PRESSURE_COEFF", "Pressure_Coefficient", "PRIMITIVE", "Pressure coefficient");

  if (config->GetKind_Solver() == MAIN_SOLVER::FEM_NAVIER_STOKES){
    AddVolumeOutput("LAMINAR_VISCOSITY", "Laminar_Viscosity", "PRIMITIVE", "Laminar viscosity");
  }

  if (config->GetKind_Solver() == MAIN_SOLVER::FEM_LES && (config->GetKind_SGS_Model() != TURB_SGS_MODEL::IMPLICIT_LES)) {
    AddVolumeOutput("EDDY_VISCOSITY", "Eddy_Viscosity", "PRIMITIVE", "Turbulent eddy viscosity");
  }
}

void CFlowCompFEMOutput::LoadVolumeDataFEM(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iElem, unsigned long index, unsigned short dof){

  unsigned short iDim;

  unsigned short nVar = solver[FLOW_SOL]->GetnVar();

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/

  auto *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

  /*--- Get a pointer to the fluid model class from the DG-FEM solver
   so that we can access the states below. ---*/

  CFluidModel *DGFluidModel = solver[FLOW_SOL]->GetFluidModel();

  /* Set the pointers for the solution for this element. */

  const unsigned long offset = nVar*volElem[iElem].offsetDOFsSolLocal;
  su2double *solDOFs         = solver[FLOW_SOL]->GetVecSolDOFs() + offset;

  /*--- Get the conservative variables for this particular DOF. ---*/

  const su2double *U = solDOFs+dof*nVar;

  /*--- Load the coordinate values of the solution DOFs. ---*/

  const su2double *coor = volElem[iElem].coorSolDOFs.data() + dof*nDim;

  /*--- Prepare the primitive states. ---*/

  const su2double DensityInv = 1.0/U[0];
  su2double vel[3], Velocity2 = 0.0;
  for(iDim=0; iDim<nDim; ++iDim) {
    vel[iDim] = U[iDim+1]*DensityInv;
    Velocity2 += vel[iDim]*vel[iDim];
  }
  su2double StaticEnergy = U[nDim+1]*DensityInv - 0.5*Velocity2;
  DGFluidModel->SetTDState_rhoe(U[0], StaticEnergy);


  SetVolumeOutputValue("COORD-X",        index, coor[0]);
  SetVolumeOutputValue("COORD-Y",        index, coor[1]);
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z",      index, coor[2]);
  SetVolumeOutputValue("DENSITY",        index, U[0]);
  SetVolumeOutputValue("MOMENTUM-X",     index, U[1]);
  SetVolumeOutputValue("MOMENTUM-Y",     index, U[2]);
  if (nDim == 3){
    SetVolumeOutputValue("MOMENTUM-Z",   index,  U[3]);
    SetVolumeOutputValue("ENERGY",       index,  U[4]);
  } else {
    SetVolumeOutputValue("ENERGY",       index,  U[3]);
  }

  SetVolumeOutputValue("PRESSURE",       index, DGFluidModel->GetPressure());
  SetVolumeOutputValue("TEMPERATURE",    index, DGFluidModel->GetTemperature());
  SetVolumeOutputValue("MACH",           index, sqrt(Velocity2)/DGFluidModel->GetSoundSpeed());
  SetVolumeOutputValue("PRESSURE_COEFF", index, DGFluidModel->GetCp());

  if (config->GetKind_Solver() == MAIN_SOLVER::FEM_NAVIER_STOKES){
    SetVolumeOutputValue("LAMINAR_VISCOSITY", index, DGFluidModel->GetLaminarViscosity());
  }
  if ((config->GetKind_Solver()  == MAIN_SOLVER::FEM_LES) && (config->GetKind_SGS_Model() != TURB_SGS_MODEL::IMPLICIT_LES)){
    // todo: Export Eddy instead of Laminar viscosity
    SetVolumeOutputValue("EDDY_VISCOSITY", index, DGFluidModel->GetLaminarViscosity());
  }
}

void CFlowCompFEMOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {

  CSolver* flow_solver = solver[FLOW_SOL];

  SetHistoryOutputValue("RMS_DENSITY", log10(flow_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("RMS_MOMENTUM-X", log10(flow_solver->GetRes_RMS(1)));
  SetHistoryOutputValue("RMS_MOMENTUM-Y", log10(flow_solver->GetRes_RMS(2)));
  if (nDim == 2)
    SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(3)));
  else {
    SetHistoryOutputValue("RMS_MOMENTUM-Z", log10(flow_solver->GetRes_RMS(3)));
    SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(4)));
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

  SetHistoryOutputValue("AOA", config->GetAoA());
  SetHistoryOutputValue("CFL_NUMBER", config->GetCFL(MESH_0));

  if (config->GetnMarker_Analyze() > 0) {
    SU2_MPI::Error("SetAnalyzeSurface is not implemented for FEM-DG solver.", CURRENT_FUNCTION);
  }

  /*--- Set aeroydnamic coefficients --- */

  SetAerodynamicCoefficients(config, flow_solver);

  ComputeSimpleCustomOutputs(config);
}

bool CFlowCompFEMOutput::SetInitResiduals(const CConfig *config){

  return (config->GetTime_Marching() != TIME_MARCHING::STEADY && (curInnerIter == 0))||
         (config->GetTime_Marching() == TIME_MARCHING::STEADY && (curTimeIter < 2));

}
