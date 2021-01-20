/*!
 * \file CIncEulerSolver.cpp
 * \brief Main subroutines for solving incompressible flow (Euler, Navier-Stokes, etc.).
 * \author F. Palacios, T. Economon
 * \version 7.1.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/solvers/CIncEulerSolver.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/fluid/CConstantDensity.hpp"
#include "../../include/fluid/CIncIdealGas.hpp"
#include "../../include/fluid/CIncIdealGasPolynomial.hpp"
#include "../../include/variables/CIncNSVariable.hpp"


CIncEulerSolver::CIncEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh,
                                 const bool navier_stokes) :
  CFVMFlowSolverBase<CIncEulerVariable, INCOMPRESSIBLE>() {

  /*--- Based on the navier_stokes boolean, determine if this constructor is
   *    being called by itself, or by its derived class CIncNSSolver. ---*/
  const string description = navier_stokes? "Navier-Stokes" : "Euler";

  unsigned short iVar, iMarker, nLineLets;
  ifstream restart_file;
  unsigned short nZone = geometry->GetnZone();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  int Unst_RestartIter;
  unsigned short iZone = config->GetiZone();
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool time_stepping = config->GetTime_Marching() == TIME_STEPPING;
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (!(!restart || (iMesh != MESH_0) || nZone > 1)) {

    /*--- Multizone problems require the number of the zone to be appended. ---*/

    auto filename_ = config->GetSolution_FileName();

    if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone, ".dat");

    /*--- Modify file name for a dual-time unsteady restart ---*/

    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetTime_Marching() == DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter, ".dat");
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/

    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter, ".dat");
    }

    /*--- Read and store the restart metadata. ---*/

//    Read_SU2_Restart_Metadata(geometry, config, false, filename_);

  }

  /*--- Set the gamma value ---*/

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometry constants in the solver structure.
   * Incompressible flow, primitive variables (P, vx, vy, vz, T, rho, beta, lamMu, EddyMu, Kt_eff, Cp, Cv) ---*/

  nDim = geometry->GetnDim();

  nVar = nDim+2; nPrimVar = nDim+9; nPrimVarGrad = nDim+4;

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nPrimVarGrad;

  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Store the number of vertices on each marker for deallocation later ---*/

  nVertex = new unsigned long[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];

  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/

  SetNondimensionalization(config, iMesh);

  /*--- Check if we are executing a verification case. If so, the
   VerificationSolution object will be instantiated for a particular
   option from the available library of verification solutions. Note
   that this is done after SetNondim(), as problem-specific initial
   parameters are needed by the solution constructors. ---*/

  SetVerificationSolution(nDim, nVar, config);

  /// TODO: This type of variables will be replaced.

  AllocateTerribleLegacyTemporaryVariables();

  /*--- Define some auxiliary vectors related to the primitive solution ---*/

  Primitive   = new su2double[nPrimVar] ();
  Primitive_i = new su2double[nPrimVar] ();
  Primitive_j = new su2double[nPrimVar] ();

  /*--- Allocate preconditioning matrix. ---*/

  Preconditioner = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar ++)
    Preconditioner[iVar] = new su2double[nVar];

  /*--- Allocate base class members. ---*/

  Allocate(*config);

  /*--- Jacobians and vector structures for implicit computations ---*/

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {

    if (rank == MASTER_NODE)
      cout << "Initialize Jacobian structure (" << description << "). MG level: " << iMesh <<"." << endl;

    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
  }
  else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (" << description << "). MG level: " << iMesh <<"." << endl;
  }

  /*--- Read farfield conditions ---*/

  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();

  /*--- Initialize the secondary values for direct derivative approxiations ---*/

  switch (config->GetDirectDiff()) {
    case NO_DERIVATIVE:
      /*--- Default ---*/
      break;
    case D_DENSITY:
      SU2_TYPE::SetDerivative(Density_Inf, 1.0);
      break;
    case D_PRESSURE:
      SU2_TYPE::SetDerivative(Pressure_Inf, 1.0);
      break;
    case D_TEMPERATURE:
      SU2_TYPE::SetDerivative(Temperature_Inf, 1.0);
      break;
    case D_MACH: case D_AOA:
    case D_SIDESLIP: case D_REYNOLDS:
    case D_TURB2LAM: case D_DESIGN:
      /*--- Already done in postprocessing of config ---*/
      break;
    default:
      break;
  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  if (navier_stokes) {
    nodes = new CIncNSVariable(Pressure_Inf, Velocity_Inf, Temperature_Inf, nPoint, nDim, nVar, config);
  } else {
    nodes = new CIncEulerVariable(Pressure_Inf, Velocity_Inf, Temperature_Inf, nPoint, nDim, nVar, config);
  }
  SetBaseClassPointerToNodes();

  /*--- Initial comms. ---*/

  CommunicateInitialState(geometry, config);

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "INC.FLOW";

  /*--- Finally, check that the static arrays will be large enough (keep this
   *    check at the bottom to make sure we consider the "final" values). ---*/
  if((nDim > MAXNDIM) || (nPrimVar > MAXNVAR))
    SU2_MPI::Error("Oops! The CIncEulerSolver static array sizes are not large enough.", CURRENT_FUNCTION);
}

CIncEulerSolver::~CIncEulerSolver(void) {

  unsigned short iVar;

  delete [] Primitive;
  delete [] Primitive_i;
  delete [] Primitive_j;

  if (Preconditioner != nullptr) {
    for (iVar = 0; iVar < nVar; iVar ++)
      delete [] Preconditioner[iVar];
    delete [] Preconditioner;
  }

  delete FluidModel;
}

void CIncEulerSolver::SetNondimensionalization(CConfig *config, unsigned short iMesh) {

  su2double Temperature_FreeStream = 0.0,  ModVel_FreeStream = 0.0,Energy_FreeStream = 0.0,
  ModVel_FreeStreamND = 0.0, Omega_FreeStream = 0.0, Omega_FreeStreamND = 0.0, Viscosity_FreeStream = 0.0,
  Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, Pressure_Thermodynamic = 0.0, Tke_FreeStream = 0.0,
  Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Temperature_Ref = 0.0, Velocity_Ref = 0.0, Time_Ref = 0.0,
  Gas_Constant_Ref = 0.0, Omega_Ref = 0.0, Force_Ref = 0.0, Viscosity_Ref = 0.0, Conductivity_Ref = 0.0, Heat_Flux_Ref = 0.0, Energy_Ref= 0.0, Pressure_FreeStreamND = 0.0, Pressure_ThermodynamicND = 0.0, Density_FreeStreamND = 0.0,
  Temperature_FreeStreamND = 0.0, Gas_ConstantND = 0.0, Specific_Heat_CpND = 0.0, Specific_Heat_CvND = 0.0, Thermal_Expansion_CoeffND = 0.0,
  Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStreamND = 0.0,
  Tke_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0,
  Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0;

  unsigned short iDim, iVar;

  /*--- Local variables ---*/

  su2double Mach     = config->GetMach();
  su2double Reynolds = config->GetReynolds();

  bool unsteady      = (config->GetTime_Marching() != NO);
  bool viscous       = config->GetViscous();
  bool turbulent     = ((config->GetKind_Solver() == INC_RANS) ||
                        (config->GetKind_Solver() == DISC_ADJ_INC_RANS));
  bool tkeNeeded     = ((turbulent) && ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)));
  bool energy        = config->GetEnergy_Equation();
  bool boussinesq    = (config->GetKind_DensityModel() == BOUSSINESQ);

  /*--- Compute dimensional free-stream values. ---*/

  Density_FreeStream     = config->GetInc_Density_Init();     config->SetDensity_FreeStream(Density_FreeStream);
  Temperature_FreeStream = config->GetInc_Temperature_Init(); config->SetTemperature_FreeStream(Temperature_FreeStream);
  Pressure_FreeStream    = 0.0; config->SetPressure_FreeStream(Pressure_FreeStream);

  ModVel_FreeStream   = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ModVel_FreeStream += config->GetInc_Velocity_Init()[iDim]*config->GetInc_Velocity_Init()[iDim];
    config->SetVelocity_FreeStream(config->GetInc_Velocity_Init()[iDim],iDim);
  }
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);

  /*--- Depending on the density model chosen, select a fluid model. ---*/

  switch (config->GetKind_FluidModel()) {

    case CONSTANT_DENSITY:

      FluidModel = new CConstantDensity(Density_FreeStream, config->GetSpecific_Heat_Cp());
      FluidModel->SetTDState_T(Temperature_FreeStream);
      break;

    case INC_IDEAL_GAS:

      config->SetGas_Constant(UNIVERSAL_GAS_CONSTANT/(config->GetMolecular_Weight()/1000.0));
      Pressure_Thermodynamic = Density_FreeStream*Temperature_FreeStream*config->GetGas_Constant();
      FluidModel = new CIncIdealGas(config->GetSpecific_Heat_Cp(), config->GetGas_Constant(), Pressure_Thermodynamic);
      FluidModel->SetTDState_T(Temperature_FreeStream);
      Pressure_Thermodynamic = FluidModel->GetPressure();
      config->SetPressure_Thermodynamic(Pressure_Thermodynamic);
      break;

    case INC_IDEAL_GAS_POLY:

      config->SetGas_Constant(UNIVERSAL_GAS_CONSTANT/(config->GetMolecular_Weight()/1000.0));
      Pressure_Thermodynamic = Density_FreeStream*Temperature_FreeStream*config->GetGas_Constant();
      FluidModel = new CIncIdealGasPolynomial<N_POLY_COEFFS>(config->GetGas_Constant(), Pressure_Thermodynamic);
      if (viscous) {
        /*--- Variable Cp model via polynomial. ---*/
        for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++)
          config->SetCp_PolyCoeffND(config->GetCp_PolyCoeff(iVar), iVar);
        FluidModel->SetCpModel(config);
      }
      FluidModel->SetTDState_T(Temperature_FreeStream);
      Pressure_Thermodynamic = FluidModel->GetPressure();
      config->SetPressure_Thermodynamic(Pressure_Thermodynamic);
      break;

    default:

      SU2_MPI::Error("Fluid model not implemented for incompressible solver.", CURRENT_FUNCTION);
      break;
  }

  if (viscous) {

    /*--- The dimensional viscosity is needed to determine the free-stream conditions.
      To accomplish this, simply set the non-dimensional coefficients to the
      dimensional ones. This will be overruled later.---*/

    config->SetMu_RefND(config->GetMu_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref());
    config->SetMu_SND(config->GetMu_S());
    config->SetMu_ConstantND(config->GetMu_Constant());

    for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++)
      config->SetMu_PolyCoeffND(config->GetMu_PolyCoeff(iVar), iVar);

    /*--- Use the fluid model to compute the dimensional viscosity/conductivity. ---*/

    FluidModel->SetLaminarViscosityModel(config);
    Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
    config->SetViscosity_FreeStream(Viscosity_FreeStream);

    Reynolds = Density_FreeStream*ModVel_FreeStream/Viscosity_FreeStream; config->SetReynolds(Reynolds);

    /*--- Turbulence kinetic energy ---*/

    Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());

  }

  /*--- The non-dim. scheme for incompressible flows uses the following ref. values:
     Reference length      = 1 m (fixed by default, grid in meters)
     Reference density     = liquid density or freestream (input)
     Reference velocity    = liquid velocity or freestream (input)
     Reference temperature = liquid temperature or freestream (input)
     Reference pressure    = Reference density * Reference velocity * Reference velocity
     Reference viscosity   = Reference Density * Reference velocity * Reference length
     This is the same non-dim. scheme as in the compressible solver.
     Note that the Re and Re Length are not used as part of initialization. ---*/

  if (config->GetRef_Inc_NonDim() == DIMENSIONAL) {
    Density_Ref     = 1.0;
    Velocity_Ref    = 1.0;
    Temperature_Ref = 1.0;
    Pressure_Ref    = 1.0;
  }
  else if (config->GetRef_Inc_NonDim() == INITIAL_VALUES) {
    Density_Ref     = Density_FreeStream;
    Velocity_Ref    = ModVel_FreeStream;
    Temperature_Ref = Temperature_FreeStream;
    Pressure_Ref    = Density_Ref*Velocity_Ref*Velocity_Ref;
  }
  else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
    Density_Ref     = config->GetInc_Density_Ref();
    Velocity_Ref    = config->GetInc_Velocity_Ref();
    Temperature_Ref = config->GetInc_Temperature_Ref();
    Pressure_Ref    = Density_Ref*Velocity_Ref*Velocity_Ref;
  }
  config->SetDensity_Ref(Density_Ref);
  config->SetVelocity_Ref(Velocity_Ref);
  config->SetTemperature_Ref(Temperature_Ref);
  config->SetPressure_Ref(Pressure_Ref);

  /*--- More derived reference values ---*/

  Length_Ref       = 1.0;                                                config->SetLength_Ref(Length_Ref);
  Time_Ref         = Length_Ref/Velocity_Ref;                            config->SetTime_Ref(Time_Ref);
  Omega_Ref        = Velocity_Ref/Length_Ref;                            config->SetOmega_Ref(Omega_Ref);
  Force_Ref        = Velocity_Ref*Velocity_Ref/Length_Ref;               config->SetForce_Ref(Force_Ref);
  Heat_Flux_Ref    = Density_Ref*Velocity_Ref*Velocity_Ref*Velocity_Ref; config->SetHeat_Flux_Ref(Heat_Flux_Ref);
  Gas_Constant_Ref = Velocity_Ref*Velocity_Ref/Temperature_Ref;          config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref    = Density_Ref*Velocity_Ref*Length_Ref;                config->SetViscosity_Ref(Viscosity_Ref);
  Conductivity_Ref = Viscosity_Ref*Gas_Constant_Ref;                     config->SetConductivity_Ref(Conductivity_Ref);

  /*--- Get the freestream energy. Only useful if energy equation is active. ---*/

  Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
  config->SetEnergy_FreeStream(Energy_FreeStream);
  if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; }; config->SetEnergy_FreeStream(Energy_FreeStream);

  /*--- Compute Mach number ---*/

  if (config->GetKind_FluidModel() == CONSTANT_DENSITY) {
    Mach = ModVel_FreeStream / sqrt(config->GetBulk_Modulus()/Density_FreeStream);
  } else {
    Mach = 0.0;
  }
  config->SetMach(Mach);

  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/

  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref(); config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Pressure_ThermodynamicND = Pressure_Thermodynamic/config->GetPressure_Ref(); config->SetPressure_ThermodynamicND(Pressure_ThermodynamicND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();   config->SetDensity_FreeStreamND(Density_FreeStreamND);

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }

  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);
  Gas_ConstantND      = config->GetGas_Constant()/Gas_Constant_Ref;    config->SetGas_ConstantND(Gas_ConstantND);
  Specific_Heat_CpND  = config->GetSpecific_Heat_Cp()/Gas_Constant_Ref; config->SetSpecific_Heat_CpND(Specific_Heat_CpND);

  /*--- We assume that Cp = Cv for our incompressible fluids. ---*/
  Specific_Heat_CvND  = config->GetSpecific_Heat_Cp()/Gas_Constant_Ref; config->SetSpecific_Heat_CvND(Specific_Heat_CvND);

  Thermal_Expansion_CoeffND = config->GetThermal_Expansion_Coeff()*config->GetTemperature_Ref(); config->SetThermal_Expansion_CoeffND(Thermal_Expansion_CoeffND);

  ModVel_FreeStreamND = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
  ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND); config->SetModVel_FreeStreamND(ModVel_FreeStreamND);

  Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;   config->SetViscosity_FreeStreamND(Viscosity_FreeStreamND);

  Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStream(Tke_FreeStream);

  Tke_FreeStreamND  = 3.0/2.0*(ModVel_FreeStreamND*ModVel_FreeStreamND*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStreamND(Tke_FreeStreamND);

  Omega_FreeStream = Density_FreeStream*Tke_FreeStream/(Viscosity_FreeStream*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStream(Omega_FreeStream);

  Omega_FreeStreamND = Density_FreeStreamND*Tke_FreeStreamND/(Viscosity_FreeStreamND*config->GetTurb2LamViscRatio_FreeStream());
  config->SetOmega_FreeStreamND(Omega_FreeStreamND);

  /*--- Delete the original (dimensional) FluidModel object. No fluid is used for inscompressible cases. ---*/

  delete FluidModel;

  switch (config->GetKind_FluidModel()) {

    case CONSTANT_DENSITY:
      FluidModel = new CConstantDensity(Density_FreeStreamND, Specific_Heat_CpND);
      break;

    case INC_IDEAL_GAS:
      FluidModel = new CIncIdealGas(Specific_Heat_CpND, Gas_ConstantND, Pressure_ThermodynamicND);
      break;

    case INC_IDEAL_GAS_POLY:
      FluidModel = new CIncIdealGasPolynomial<N_POLY_COEFFS>(Gas_ConstantND, Pressure_ThermodynamicND);
      if (viscous) {
        /*--- Variable Cp model via polynomial. ---*/
        config->SetCp_PolyCoeffND(config->GetCp_PolyCoeff(0)/Gas_Constant_Ref, 0);
        for (iVar = 1; iVar < config->GetnPolyCoeffs(); iVar++)
          config->SetCp_PolyCoeffND(config->GetCp_PolyCoeff(iVar)*pow(Temperature_Ref,iVar)/Gas_Constant_Ref, iVar);
        FluidModel->SetCpModel(config);
      }
      break;
      FluidModel->SetTDState_T(Temperature_FreeStreamND);
  }

  Energy_FreeStreamND = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;

  if (viscous) {

    /*--- Constant viscosity model ---*/

    config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);

    /*--- Sutherland's model ---*/

    config->SetMu_RefND(config->GetMu_Ref()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_S()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref()/config->GetTemperature_Ref());

    /*--- Viscosity model via polynomial. ---*/

    config->SetMu_PolyCoeffND(config->GetMu_PolyCoeff(0)/Viscosity_Ref, 0);
    for (iVar = 1; iVar < config->GetnPolyCoeffs(); iVar++)
      config->SetMu_PolyCoeffND(config->GetMu_PolyCoeff(iVar)*pow(Temperature_Ref,iVar)/Viscosity_Ref, iVar);

    /*--- Constant thermal conductivity model ---*/

    config->SetKt_ConstantND(config->GetKt_Constant()/Conductivity_Ref);

    /*--- Conductivity model via polynomial. ---*/

    config->SetKt_PolyCoeffND(config->GetKt_PolyCoeff(0)/Conductivity_Ref, 0);
    for (iVar = 1; iVar < config->GetnPolyCoeffs(); iVar++)
      config->SetKt_PolyCoeffND(config->GetKt_PolyCoeff(iVar)*pow(Temperature_Ref,iVar)/Conductivity_Ref, iVar);

    /*--- Set up the transport property models. ---*/

    FluidModel->SetLaminarViscosityModel(config);
    FluidModel->SetThermalConductivityModel(config);

  }

  if (tkeNeeded) { Energy_FreeStreamND += Tke_FreeStreamND; };  config->SetEnergy_FreeStreamND(Energy_FreeStreamND);

  Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);

  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);

  /*--- Write output to the console if this is the master node and first domain ---*/

  if ((rank == MASTER_NODE) && (iMesh == MESH_0)) {

    cout.precision(6);

    if (config->GetRef_Inc_NonDim() == DIMENSIONAL) {
      cout << "Incompressible flow: rho_ref, vel_ref, temp_ref, p_ref" << endl;
      cout << "are set to 1.0 in order to perform a dimensional calculation." << endl;
      if (dynamic_grid) cout << "Force coefficients computed using MACH_MOTION." << endl;
      else cout << "Force coefficients computed using initial values." << endl;
    }
    else if (config->GetRef_Inc_NonDim() == INITIAL_VALUES) {
      cout << "Incompressible flow: rho_ref, vel_ref, and temp_ref" << endl;
      cout << "are based on the initial values, p_ref = rho_ref*vel_ref^2." << endl;
      if (dynamic_grid) cout << "Force coefficients computed using MACH_MOTION." << endl;
      else cout << "Force coefficients computed using initial values." << endl;
    }
    else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
      cout << "Incompressible flow: rho_ref, vel_ref, and temp_ref" << endl;
      cout << "are user-provided reference values, p_ref = rho_ref*vel_ref^2." << endl;
      if (dynamic_grid) cout << "Force coefficients computed using MACH_MOTION." << endl;
      else cout << "Force coefficients computed using reference values." << endl;
    }
    cout << "The reference area for force coeffs. is " << config->GetRefArea() << " m^2." << endl;
    cout << "The reference length for force coeffs. is " << config->GetRefLength() << " m." << endl;

    cout << "The pressure is decomposed into thermodynamic and dynamic components." << endl;
    cout << "The initial value of the dynamic pressure is 0." << endl;

    cout << "Mach number: "<< config->GetMach();
    if (config->GetKind_FluidModel() == CONSTANT_DENSITY) {
      cout << ", computed using the Bulk modulus." << endl;
    } else {
      cout << ", computed using fluid speed of sound." << endl;
    }

    cout << "For external flows, the initial state is imposed at the far-field." << endl;
    cout << "Angle of attack (deg): "<< config->GetAoA() << ", computed using the initial velocity." << endl;
    cout << "Side slip angle (deg): "<< config->GetAoS() << ", computed using the initial velocity." << endl;

    if (viscous) {
      cout << "Reynolds number per meter: " << config->GetReynolds() << ", computed using initial values."<< endl;
      cout << "Reynolds number is a byproduct of inputs only (not used internally)." << endl;
    }
    cout << "SI units only. The grid should be dimensional (meters)." << endl;

    switch (config->GetKind_DensityModel()) {

      case CONSTANT:
        if (energy) cout << "Energy equation is active and decoupled." << endl;
        else cout << "No energy equation." << endl;
        break;

      case BOUSSINESQ:
        if (energy) cout << "Energy equation is active and coupled through Boussinesq approx." << endl;
        break;

      case VARIABLE:
        if (energy) cout << "Energy equation is active and coupled for variable density." << endl;
        break;

    }

    stringstream NonDimTableOut, ModelTableOut;
    stringstream Unit;

    cout << endl;
    PrintingToolbox::CTablePrinter ModelTable(&ModelTableOut);
    ModelTableOut <<"-- Models:"<< endl;

    ModelTable.AddColumn("Viscosity Model", 25);
    ModelTable.AddColumn("Conductivity Model", 26);
    ModelTable.AddColumn("Fluid Model", 25);
    ModelTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);
    ModelTable.PrintHeader();

    PrintingToolbox::CTablePrinter NonDimTable(&NonDimTableOut);
    NonDimTable.AddColumn("Name", 22);
    NonDimTable.AddColumn("Dim. value", 14);
    NonDimTable.AddColumn("Ref. value", 14);
    NonDimTable.AddColumn("Unit", 10);
    NonDimTable.AddColumn("Non-dim. value", 14);
    NonDimTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);

    NonDimTableOut <<"-- Fluid properties:"<< endl;

    NonDimTable.PrintHeader();

    if (viscous){

      switch(config->GetKind_ViscosityModel()){
      case CONSTANT_VISCOSITY:
        ModelTable << "CONSTANT_VISCOSITY";
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Viscosity" << config->GetMu_Constant() << config->GetMu_Constant()/config->GetMu_ConstantND() << Unit.str() << config->GetMu_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case SUTHERLAND:
        ModelTable << "SUTHERLAND";
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Ref. Viscosity" <<  config->GetMu_Ref() <<  config->GetViscosity_Ref() << Unit.str() << config->GetMu_RefND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "K";
        else if (config->GetSystemMeasurements() == US) Unit << "R";
        NonDimTable << "Sutherland Temp." << config->GetMu_Temperature_Ref() <<  config->GetTemperature_Ref() << Unit.str() << config->GetMu_Temperature_RefND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "K";
        else if (config->GetSystemMeasurements() == US) Unit << "R";
        NonDimTable << "Sutherland Const." << config->GetMu_S() << config->GetTemperature_Ref() << Unit.str() << config->GetMu_SND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case POLYNOMIAL_VISCOSITY:
        ModelTable << "POLYNOMIAL_VISCOSITY";
        for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
          stringstream ss;
          ss << iVar;
          if (config->GetMu_PolyCoeff(iVar) != 0.0)
            NonDimTable << "Mu(T) Poly. Coeff. " + ss.str()  << config->GetMu_PolyCoeff(iVar) << config->GetMu_PolyCoeff(iVar)/config->GetMu_PolyCoeffND(iVar) << "-" << config->GetMu_PolyCoeffND(iVar);
        }
        Unit.str("");
        NonDimTable.PrintFooter();
        break;
      }

      switch(config->GetKind_ConductivityModel()){
      case CONSTANT_PRANDTL:
        ModelTable << "CONSTANT_PRANDTL";
        NonDimTable << "Prandtl (Lam.)"  << "-" << "-" << "-" << config->GetPrandtl_Lam();
        Unit.str("");
        NonDimTable << "Prandtl (Turb.)" << "-" << "-" << "-" << config->GetPrandtl_Turb();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case CONSTANT_CONDUCTIVITY:
        ModelTable << "CONSTANT_CONDUCTIVITY";
        Unit << "W/m^2.K";
        NonDimTable << "Molecular Cond." << config->GetKt_Constant() << config->GetKt_Constant()/config->GetKt_ConstantND() << Unit.str() << config->GetKt_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case POLYNOMIAL_CONDUCTIVITY:
        ModelTable << "POLYNOMIAL_CONDUCTIVITY";
        for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
          stringstream ss;
          ss << iVar;
          if (config->GetKt_PolyCoeff(iVar) != 0.0)
            NonDimTable << "Kt(T) Poly. Coeff. " + ss.str()  << config->GetKt_PolyCoeff(iVar) << config->GetKt_PolyCoeff(iVar)/config->GetKt_PolyCoeffND(iVar) << "-" << config->GetKt_PolyCoeffND(iVar);
        }
        Unit.str("");
        NonDimTable.PrintFooter();
        break;
      }
    } else {
      ModelTable << "-" << "-";
    }

    switch (config->GetKind_FluidModel()){
    case CONSTANT_DENSITY:
      ModelTable << "CONSTANT_DENSITY";
      if (energy){
        Unit << "N.m/kg.K";
        NonDimTable << "Spec. Heat (Cp)" << config->GetSpecific_Heat_Cp() << config->GetSpecific_Heat_Cp()/config->GetSpecific_Heat_CpND() << Unit.str() << config->GetSpecific_Heat_CpND();
        Unit.str("");
      }
      if (boussinesq){
        Unit << "K^-1";
        NonDimTable << "Thermal Exp." << config->GetThermal_Expansion_Coeff() << config->GetThermal_Expansion_Coeff()/config->GetThermal_Expansion_CoeffND() << Unit.str() <<  config->GetThermal_Expansion_CoeffND();
        Unit.str("");
      }
      Unit << "Pa";
      NonDimTable << "Bulk Modulus" << config->GetBulk_Modulus() << 1.0 << Unit.str() <<  config->GetBulk_Modulus();
      Unit.str("");
      NonDimTable.PrintFooter();
      break;

    case INC_IDEAL_GAS:
      ModelTable << "INC_IDEAL_GAS";
      Unit << "N.m/kg.K";
      NonDimTable << "Spec. Heat (Cp)" << config->GetSpecific_Heat_Cp() << config->GetSpecific_Heat_Cp()/config->GetSpecific_Heat_CpND() << Unit.str() << config->GetSpecific_Heat_CpND();
      Unit.str("");
      Unit << "g/mol";
      NonDimTable << "Molecular weight" << config->GetMolecular_Weight()<< 1.0 << Unit.str() << config->GetMolecular_Weight();
      Unit.str("");
      Unit << "N.m/kg.K";
      NonDimTable << "Gas Constant" << config->GetGas_Constant() << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();
      Unit.str("");
      Unit << "Pa";
      NonDimTable << "Therm. Pressure" << config->GetPressure_Thermodynamic() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_ThermodynamicND();
      Unit.str("");
      NonDimTable.PrintFooter();
      break;

    case INC_IDEAL_GAS_POLY:
      ModelTable << "INC_IDEAL_GAS_POLY";
      Unit.str("");
      Unit << "g/mol";
      NonDimTable << "Molecular weight" << config->GetMolecular_Weight()<< 1.0 << Unit.str() << config->GetMolecular_Weight();
      Unit.str("");
      Unit << "N.m/kg.K";
      NonDimTable << "Gas Constant" << config->GetGas_Constant() << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();
      Unit.str("");
      Unit << "Pa";
      NonDimTable << "Therm. Pressure" << config->GetPressure_Thermodynamic() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_ThermodynamicND();
      Unit.str("");
      for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
        stringstream ss;
        ss << iVar;
        if (config->GetCp_PolyCoeff(iVar) != 0.0)
          NonDimTable << "Cp(T) Poly. Coeff. " + ss.str()  << config->GetCp_PolyCoeff(iVar) << config->GetCp_PolyCoeff(iVar)/config->GetCp_PolyCoeffND(iVar) << "-" << config->GetCp_PolyCoeffND(iVar);
      }
      Unit.str("");
      NonDimTable.PrintFooter();
      break;

    }


    NonDimTableOut <<"-- Initial and free-stream conditions:"<< endl;
    NonDimTable.PrintHeader();

    if      (config->GetSystemMeasurements() == SI) Unit << "Pa";
    else if (config->GetSystemMeasurements() == US) Unit << "psf";
    NonDimTable << "Dynamic Pressure" << config->GetPressure_FreeStream() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "Pa";
    else if (config->GetSystemMeasurements() == US) Unit << "psf";
    NonDimTable << "Total Pressure" << config->GetPressure_FreeStream() + 0.5*config->GetDensity_FreeStream()*config->GetModVel_FreeStream()*config->GetModVel_FreeStream()
                << config->GetPressure_Ref() << Unit.str() << config->GetPressure_FreeStreamND() + 0.5*config->GetDensity_FreeStreamND()*config->GetModVel_FreeStreamND()*config->GetModVel_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "kg/m^3";
    else if (config->GetSystemMeasurements() == US) Unit << "slug/ft^3";
    NonDimTable << "Density" << config->GetDensity_FreeStream() << config->GetDensity_Ref() << Unit.str() << config->GetDensity_FreeStreamND();
    Unit.str("");
    if (energy){
      if      (config->GetSystemMeasurements() == SI) Unit << "K";
      else if (config->GetSystemMeasurements() == US) Unit << "R";
      NonDimTable << "Temperature" << config->GetTemperature_FreeStream() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_FreeStreamND();
      Unit.str("");
    }
    if      (config->GetSystemMeasurements() == SI) Unit << "m/s";
    else if (config->GetSystemMeasurements() == US) Unit << "ft/s";
    NonDimTable << "Velocity-X" << config->GetVelocity_FreeStream()[0] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[0];
    NonDimTable << "Velocity-Y" << config->GetVelocity_FreeStream()[1] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[1];
    if (nDim == 3){
      NonDimTable << "Velocity-Z" << config->GetVelocity_FreeStream()[2] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[2];
    }
    NonDimTable << "Velocity Magnitude" << config->GetModVel_FreeStream() << config->GetVelocity_Ref() << Unit.str() << config->GetModVel_FreeStreamND();
    Unit.str("");

    if (viscous){
      NonDimTable.PrintFooter();
      if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
      else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
      NonDimTable << "Viscosity" << config->GetViscosity_FreeStream() << config->GetViscosity_Ref() << Unit.str() << config->GetViscosity_FreeStreamND();
      Unit.str("");
      if      (config->GetSystemMeasurements() == SI) Unit << "W/m^2.K";
      else if (config->GetSystemMeasurements() == US) Unit << "lbf/ft.s.R";
      NonDimTable << "Conductivity" << "-" << config->GetConductivity_Ref() << Unit.str() << "-";
      Unit.str("");
      if (turbulent){
        if      (config->GetSystemMeasurements() == SI) Unit << "m^2/s^2";
        else if (config->GetSystemMeasurements() == US) Unit << "ft^2/s^2";
        NonDimTable << "Turb. Kin. Energy" << config->GetTke_FreeStream() << config->GetTke_FreeStream()/config->GetTke_FreeStreamND() << Unit.str() << config->GetTke_FreeStreamND();
        Unit.str("");
        if      (config->GetSystemMeasurements() == SI) Unit << "1/s";
        else if (config->GetSystemMeasurements() == US) Unit << "1/s";
        NonDimTable << "Spec. Dissipation" << config->GetOmega_FreeStream() << config->GetOmega_FreeStream()/config->GetOmega_FreeStreamND() << Unit.str() << config->GetOmega_FreeStreamND();
        Unit.str("");
      }
    }

    NonDimTable.PrintFooter();
    NonDimTable << "Mach Number" << "-" << "-" << "-" << config->GetMach();
    if (viscous){
      NonDimTable << "Reynolds Number" << "-" << "-" << "-" << config->GetReynolds();
    }

    NonDimTable.PrintFooter();
    ModelTable.PrintFooter();

    if (unsteady){
      NonDimTableOut << "-- Unsteady conditions" << endl;
      NonDimTable.PrintHeader();
      NonDimTable << "Total Time" << config->GetMax_Time() << config->GetTime_Ref() << "s" << config->GetMax_Time()/config->GetTime_Ref();
      Unit.str("");
      NonDimTable << "Time Step" << config->GetTime_Step() << config->GetTime_Ref() << "s" << config->GetDelta_UnstTimeND();
      Unit.str("");
      NonDimTable.PrintFooter();
    }


    cout << ModelTableOut.str();
    cout << NonDimTableOut.str();
  }



}

void CIncEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  unsigned long iPoint, Point_Fine;
  unsigned short iMesh, iChildren, iVar;
  su2double Area_Children, Area_Parent, *Solution_Fine, *Solution;

  const bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  const bool rans = (config->GetKind_Turb_Model() != NONE);
  const bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                          (config->GetTime_Marching() == DT_STEPPING_2ND));

  /*--- Check if a verification solution is to be computed. ---*/
  if ((VerificationSolution) && (TimeIter == 0) && !restart) {

    /*--- Loop over the multigrid levels. ---*/
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {

      /*--- Loop over all grid points. ---*/
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {

        /* Set the pointers to the coordinates and solution of this DOF. */
        const su2double *coor = geometry[iMesh]->nodes->GetCoord(iPoint);
        su2double *solDOF     = solver_container[iMesh][FLOW_SOL]->GetNodes()->GetSolution(iPoint);

        /* Set the solution in this DOF to the initial condition provided by
           the verification solution class. This can be the exact solution,
           but this is not necessary. */
        VerificationSolution->GetInitialCondition(coor, solDOF);
      }
    }
  }

  /*--- If restart solution, then interpolate the flow solution to
   all the multigrid levels, this is important with the dual time strategy ---*/

  if (restart && (TimeIter == 0)) {

    Solution = new su2double[nVar];
    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        Area_Parent = geometry[iMesh]->nodes->GetVolume(iPoint);
        for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
        for (iChildren = 0; iChildren < geometry[iMesh]->nodes->GetnChildren_CV(iPoint); iChildren++) {
          Point_Fine = geometry[iMesh]->nodes->GetChildren_CV(iPoint, iChildren);
          Area_Children = geometry[iMesh-1]->nodes->GetVolume(Point_Fine);
          Solution_Fine = solver_container[iMesh-1][FLOW_SOL]->GetNodes()->GetSolution(Point_Fine);
          for (iVar = 0; iVar < nVar; iVar++) {
            Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
          }
        }
        solver_container[iMesh][FLOW_SOL]->GetNodes()->SetSolution(iPoint,Solution);
      }
      solver_container[iMesh][FLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
      solver_container[iMesh][FLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    }
    delete [] Solution;

    /*--- Interpolate the turblence variable also, if needed ---*/

    if (rans) {

      unsigned short nVar_Turb = solver_container[MESH_0][TURB_SOL]->GetnVar();
      Solution = new su2double[nVar_Turb];
      for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          Area_Parent = geometry[iMesh]->nodes->GetVolume(iPoint);
          for (iVar = 0; iVar < nVar_Turb; iVar++) Solution[iVar] = 0.0;
          for (iChildren = 0; iChildren < geometry[iMesh]->nodes->GetnChildren_CV(iPoint); iChildren++) {
            Point_Fine = geometry[iMesh]->nodes->GetChildren_CV(iPoint, iChildren);
            Area_Children = geometry[iMesh-1]->nodes->GetVolume(Point_Fine);
            Solution_Fine = solver_container[iMesh-1][TURB_SOL]->GetNodes()->GetSolution(Point_Fine);
            for (iVar = 0; iVar < nVar_Turb; iVar++) {
              Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
            }
          }
          solver_container[iMesh][TURB_SOL]->GetNodes()->SetSolution(iPoint,Solution);
        }
        solver_container[iMesh][TURB_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION_EDDY);
        solver_container[iMesh][TURB_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION_EDDY);
        solver_container[iMesh][TURB_SOL]->Postprocessing(geometry[iMesh], solver_container[iMesh], config, iMesh);
      }
      delete [] Solution;
    }

  }

  /*--- The value of the solution for the first iteration of the dual time ---*/

  if (dual_time && (TimeIter == 0 || (restart && TimeIter == config->GetRestart_Iter()))) {
    PushSolutionBackInTime(TimeIter, restart, rans, solver_container, geometry, config);
  }
}

void CIncEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long ErrorCounter = 0;

  unsigned long InnerIter = config->GetInnerIter();
  bool cont_adjoint     = config->GetContinuous_Adjoint();
  bool implicit         = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool muscl            = (config->GetMUSCL_Flow() || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == ROE));
  bool limiter          = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  bool center           = ((config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (cont_adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED));
  bool center_jst       = center && (config->GetKind_Centered_Flow() == JST);
  bool van_albada       = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;
  bool outlet           = ((config->GetnMarker_Outlet() != 0));

  /*--- Set the primitive variables ---*/

  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

  /*--- Upwind second order reconstruction ---*/

  if ((muscl && !center) && (iMesh == MESH_0) && !Output) {

    /*--- Gradient computation for MUSCL reconstruction. ---*/

    if (config->GetKind_Gradient_Method_Recon() == GREEN_GAUSS)
      SetPrimitive_Gradient_GG(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == LEAST_SQUARES)
      SetPrimitive_Gradient_LS(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES)
      SetPrimitive_Gradient_LS(geometry, config, true);

    /*--- Limiter computation ---*/

    if ((limiter) && (iMesh == MESH_0) && !Output && !van_albada) {
      SetPrimitive_Limiter(geometry, config);
    }

  }

  /*--- Artificial dissipation ---*/

  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }

  /*--- Update the beta value based on the maximum velocity. ---*/

  SetBeta_Parameter(geometry, solver_container, config, iMesh);

  /*--- Compute properties needed for mass flow BCs. ---*/

  if (outlet) GetOutlet_Properties(geometry, config, iMesh, Output);

  /*--- Initialize the Jacobian matrices ---*/

  if (implicit && !Output) Jacobian.SetValZero();

  /*--- Error message ---*/

  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
  }

}

void CIncEulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                  unsigned short iMesh) { }

unsigned long CIncEulerSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {

  unsigned long iPoint, nonPhysicalPoints = 0;
  bool physical = true;

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Incompressible flow, primitive variables ---*/

    physical = nodes->SetPrimVar(iPoint,FluidModel);

    /* Check for non-realizable states for reporting. */

    if (!physical) nonPhysicalPoints++;

    /*--- Initialize the convective, source and viscous residual vector ---*/

    if (!Output) LinSysRes.SetBlock_Zero(iPoint);

  }

  return nonPhysicalPoints;
}

void CIncEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh, unsigned long Iteration) {

  su2double Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0,
  Mean_BetaInc2, Lambda, Local_Delta_Time,
  Global_Delta_Time = 1E6, Global_Delta_UnstTimeND, ProjVel, ProjVel_i, ProjVel_j;
  const su2double* Normal;

  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;

  bool implicit      = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool time_stepping = config->GetTime_Marching() == TIME_STEPPING;
  bool dual_time     = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));

  Min_Delta_Time = 1.E30; Max_Delta_Time = 0.0;

  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    nodes->SetMax_Lambda_Inv(iPoint,0.0);

  /*--- Loop interior edges ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Point identification, Normal vector and area ---*/

    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    Normal = geometry->edges->GetNormal(iEdge);

    Area = GeometryToolbox::Norm(nDim, Normal);

    /*--- Mean Values ---*/

    Mean_ProjVel    = 0.5 * (nodes->GetProjVel(iPoint,Normal) + nodes->GetProjVel(jPoint,Normal));
    Mean_BetaInc2   = 0.5 * (nodes->GetBetaInc2(iPoint)      + nodes->GetBetaInc2(jPoint));
    Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);

    /*--- Adjustment for grid movement ---*/

    if (dynamic_grid) {
      su2double *GridVel_i = geometry->nodes->GetGridVel(iPoint);
      su2double *GridVel_j = geometry->nodes->GetGridVel(jPoint);
      ProjVel_i = 0.0; ProjVel_j = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }

    /*--- Inviscid contribution ---*/

    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->nodes->GetDomain(iPoint)) nodes->AddMax_Lambda_Inv(iPoint,Lambda);
    if (geometry->nodes->GetDomain(jPoint)) nodes->AddMax_Lambda_Inv(jPoint,Lambda);

  }

  /*--- Loop boundary edges ---*/

  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

      /*--- Point identification, Normal vector and area ---*/

      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

      Area = GeometryToolbox::Norm(nDim, Normal);

      /*--- Mean Values ---*/

      Mean_ProjVel    = nodes->GetProjVel(iPoint,Normal);
      Mean_BetaInc2   = nodes->GetBetaInc2(iPoint);
      Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);

      /*--- Adjustment for grid movement ---*/

      if (dynamic_grid) {
        su2double *GridVel = geometry->nodes->GetGridVel(iPoint);
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }

      /*--- Inviscid contribution ---*/

      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->nodes->GetDomain(iPoint)) {
        nodes->AddMax_Lambda_Inv(iPoint,Lambda);
      }

    }
    }
  }

  /*--- Local time-stepping: each element uses their own speed for steady state
   simulations or for pseudo time steps in a dual time simulation. ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    Vol = geometry->nodes->GetVolume(iPoint);

    if (Vol != 0.0) {
      Local_Delta_Time  = nodes->GetLocalCFL(iPoint)*Vol / nodes->GetMax_Lambda_Inv(iPoint);
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time    = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time    = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      nodes->SetDelta_Time(iPoint,Local_Delta_Time);
    }
    else {
      nodes->SetDelta_Time(iPoint,0.0);
    }

  }

  /*--- Compute the max and the min dt (in parallel) ---*/

  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Min_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Min_Delta_Time = rbuf_time;

    sbuf_time = Max_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Max_Delta_Time = rbuf_time;
#endif
  }

  /*--- For time-accurate simulations use the minimum delta time of the whole mesh (global) ---*/

  if (time_stepping) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    /*--- If the unsteady CFL is set to zero, it uses the defined
     unsteady time step, otherwise it computes the time step based
     on the unsteady CFL ---*/

    if (config->GetUnst_CFL() == 0.0) {
      Global_Delta_Time = config->GetDelta_UnstTime();
    }
    config->SetDelta_UnstTimeND(Global_Delta_Time);
    for (iPoint = 0; iPoint < nPointDomain; iPoint++){

      /*--- Sets the regular CFL equal to the unsteady CFL ---*/

      nodes->SetLocalCFL(iPoint, config->GetUnst_CFL());
      nodes->SetDelta_Time(iPoint, Global_Delta_Time);
      Min_Delta_Time = Global_Delta_Time;
      Max_Delta_Time = Global_Delta_Time;

    }
  }

  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/

  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {

    Global_Delta_UnstTimeND = 1e30;
    for (iPoint = 0; iPoint < nPointDomain; iPoint++){
      Global_Delta_UnstTimeND = min(Global_Delta_UnstTimeND,config->GetUnst_CFL()*Global_Delta_Time/nodes->GetLocalCFL(iPoint));
    }

#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_UnstTimeND;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_UnstTimeND = rbuf_time;
#endif
    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }

  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/

  if (dual_time)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!implicit) {
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), nodes->GetDelta_Time(iPoint));
        nodes->SetDelta_Time(iPoint,Local_Delta_Time);
      }
    }

}

void CIncEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  CNumerics* numerics = numerics_container[CONV_TERM];

  unsigned long iEdge, iPoint, jPoint;

  bool implicit    = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool jst_scheme  = ((config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0));

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge, set normal vectors, and number of neighbors ---*/

    iPoint = geometry->edges->GetNode(iEdge,0); jPoint = geometry->edges->GetNode(iEdge,1);
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));
    numerics->SetNeighbor(geometry->nodes->GetnNeighbor(iPoint), geometry->nodes->GetnNeighbor(jPoint));

    /*--- Set primitive variables w/o reconstruction ---*/

    numerics->SetPrimitive(nodes->GetPrimitive(iPoint), nodes->GetPrimitive(jPoint));

    /*--- Set the largest convective eigenvalue ---*/

    numerics->SetLambda(nodes->GetLambda(iPoint), nodes->GetLambda(jPoint));

    /*--- Set undivided laplacian and pressure-based sensor ---*/

    if (jst_scheme) {
      numerics->SetUndivided_Laplacian(nodes->GetUndivided_Laplacian(iPoint), nodes->GetUndivided_Laplacian(jPoint));
      numerics->SetSensor(nodes->GetSensor(iPoint), nodes->GetSensor(jPoint));
    }

    /*--- Grid movement ---*/

    if (dynamic_grid) {
      numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(jPoint));
    }

    /*--- Compute residuals, and Jacobians ---*/

    auto residual = numerics->ComputeResidual(config);

    /*--- Update convective and artificial dissipation residuals ---*/

    LinSysRes.AddBlock(iPoint, residual);
    LinSysRes.SubtractBlock(jPoint, residual);

    /*--- Store implicit contributions from the residual calculation. ---*/

    if (implicit) {
      Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
    }
  }

}

void CIncEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                                      CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  CNumerics* numerics = numerics_container[CONV_TERM];

  su2double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j,
  *V_i, *V_j, *S_i, *S_j, *Limiter_i = nullptr, *Limiter_j = nullptr;

  unsigned long iEdge, iPoint, jPoint, counter_local = 0, counter_global = 0;
  unsigned short iDim, iVar;

  unsigned long InnerIter = config->GetInnerIter();
  bool implicit         = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool muscl            = (config->GetMUSCL_Flow() && (iMesh == MESH_0));
  bool limiter          = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  bool van_albada       = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;

  /*--- Loop over all the edges ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge and normal vectors ---*/

    iPoint = geometry->edges->GetNode(iEdge,0); jPoint = geometry->edges->GetNode(iEdge,1);
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    /*--- Grid movement ---*/

    if (dynamic_grid)
      numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(jPoint));

    /*--- Get primitive variables ---*/

    V_i = nodes->GetPrimitive(iPoint); V_j = nodes->GetPrimitive(jPoint);
    S_i = nodes->GetSecondary(iPoint); S_j = nodes->GetSecondary(jPoint);

    /*--- High order reconstruction using MUSCL strategy ---*/

    if (muscl) {

      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->nodes->GetCoord(jPoint, iDim) - geometry->nodes->GetCoord(iPoint, iDim));
        Vector_j[iDim] = 0.5*(geometry->nodes->GetCoord(iPoint, iDim) - geometry->nodes->GetCoord(jPoint, iDim));
      }

      Gradient_i = nodes->GetGradient_Reconstruction(iPoint);
      Gradient_j = nodes->GetGradient_Reconstruction(jPoint);

      if (limiter) {
        Limiter_i = nodes->GetLimiter_Primitive(iPoint);
        Limiter_j = nodes->GetLimiter_Primitive(jPoint);
      }

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        if (limiter) {
          if (van_albada){
            Limiter_i[iVar] = (V_j[iVar]-V_i[iVar])*(2.0*Project_Grad_i + V_j[iVar]-V_i[iVar])/(4*Project_Grad_i*Project_Grad_i+(V_j[iVar]-V_i[iVar])*(V_j[iVar]-V_i[iVar])+EPS);
            Limiter_j[iVar] = (V_j[iVar]-V_i[iVar])*(-2.0*Project_Grad_j + V_j[iVar]-V_i[iVar])/(4*Project_Grad_j*Project_Grad_j+(V_j[iVar]-V_i[iVar])*(V_j[iVar]-V_i[iVar])+EPS);
          }
          Primitive_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          Primitive_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          Primitive_i[iVar] = V_i[iVar] + Project_Grad_i;
          Primitive_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }

      for (iVar = nPrimVarGrad; iVar < nPrimVar; iVar++) {
        Primitive_i[iVar] = V_i[iVar];
        Primitive_j[iVar] = V_j[iVar];
      }

      /*--- Check for non-physical solutions after reconstruction. If found,
       use the cell-average value of the solution. This results in a locally
       first-order approximation, but this is typically only active
       during the start-up of a calculation or difficult transients. For
       incompressible flow, only the temperature and density need to be
       checked. Pressure is the dynamic pressure (can be negative). ---*/

      if (config->GetEnergy_Equation()) {
        bool neg_temperature_i = (Primitive_i[nDim+1] < 0.0);
        bool neg_temperature_j = (Primitive_j[nDim+1] < 0.0);

        bool neg_density_i  = (Primitive_i[nDim+2] < 0.0);
        bool neg_density_j  = (Primitive_j[nDim+2] < 0.0);

        if (neg_density_i || neg_temperature_i) {
          nodes->SetNon_Physical(iPoint, true);
        } else {
          nodes->SetNon_Physical(iPoint, false);
        }

        if (neg_density_j || neg_temperature_j) {
          nodes->SetNon_Physical(jPoint, true);
        } else {
          nodes->SetNon_Physical(jPoint, false);
        }

        /* Lastly, check for existing first-order points still active
         from previous iterations. */

        if (nodes->GetNon_Physical(iPoint)) {
          counter_local++;
          for (iVar = 0; iVar < nPrimVar; iVar++)
            Primitive_i[iVar] = V_i[iVar];
        }
        if (nodes->GetNon_Physical(jPoint)) {
          counter_local++;
          for (iVar = 0; iVar < nPrimVar; iVar++)
            Primitive_j[iVar] = V_j[iVar];
        }
      }

      numerics->SetPrimitive(Primitive_i, Primitive_j);

    } else {

      /*--- Set conservative variables without reconstruction ---*/

      numerics->SetPrimitive(V_i, V_j);
      numerics->SetSecondary(S_i, S_j);

    }

    /*--- Compute the residual ---*/

    auto residual = numerics->ComputeResidual(config);

    /*--- Update residual value ---*/

    LinSysRes.AddBlock(iPoint, residual);
    LinSysRes.SubtractBlock(jPoint, residual);

    /*--- Set implicit Jacobians ---*/

    if (implicit) {
      Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
    }
  }

  /*--- Warning message about non-physical reconstructions. ---*/

  if (config->GetComm_Level() == COMM_FULL) {
    if (iMesh == MESH_0) {
      SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
      config->SetNonphysical_Reconstr(counter_global);
    }
  }

}

void CIncEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                      CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM];

  unsigned short iVar;
  unsigned long iPoint;

  const bool implicit       = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool rotating_frame = config->GetRotating_Frame();
  const bool axisymmetric   = config->GetAxisymmetric();
  const bool body_force     = config->GetBody_Force();
  const bool boussinesq     = (config->GetKind_DensityModel() == BOUSSINESQ);
  const bool viscous        = config->GetViscous();
  const bool radiation      = config->AddRadiation();
  const bool vol_heat       = config->GetHeatSource();

  if (body_force) {

    /*--- Loop over all points ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the conservative variables ---*/

      numerics->SetConservative(nodes->GetSolution(iPoint),
                                nodes->GetSolution(iPoint));

      /*--- Set incompressible density  ---*/

      numerics->SetDensity(nodes->GetDensity(iPoint),
                           nodes->GetDensity(iPoint));

      /*--- Load the volume of the dual mesh cell ---*/

      numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Compute the body force source residual ---*/

      auto residual = numerics->ComputeResidual(config);

      /*--- Add the source residual to the total ---*/

      LinSysRes.AddBlock(iPoint, residual);

    }
  }

  if (boussinesq) {

    /*--- Loop over all points ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the conservative variables ---*/

      numerics->SetConservative(nodes->GetSolution(iPoint),
                                nodes->GetSolution(iPoint));

      /*--- Set incompressible density  ---*/

      numerics->SetDensity(nodes->GetDensity(iPoint),
                           nodes->GetDensity(iPoint));

      /*--- Load the volume of the dual mesh cell ---*/

      numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Compute the boussinesq source residual ---*/

      auto residual = numerics->ComputeResidual(config);

      /*--- Add the source residual to the total ---*/

      LinSysRes.AddBlock(iPoint, residual);

    }
  }

  if (rotating_frame) {

    /*--- Loop over all points ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the primitive variables ---*/

      numerics->SetPrimitive(nodes->GetPrimitive(iPoint), nullptr);

      /*--- Set incompressible density ---*/

      numerics->SetDensity(nodes->GetDensity(iPoint), 0.0);

      /*--- Load the volume of the dual mesh cell ---*/

      numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Compute the rotating frame source residual ---*/

      auto residual = numerics->ComputeResidual(config);

      /*--- Add the source residual to the total ---*/

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Add the implicit Jacobian contribution ---*/

      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }

  if (axisymmetric) {

    /*--- For viscous problems, we need an additional gradient. ---*/

    if (viscous) {

      for (iPoint = 0; iPoint < nPoint; iPoint++) {

        su2double yCoord          = geometry->nodes->GetCoord(iPoint, 1);
        su2double yVelocity       = nodes->GetVelocity(iPoint,1);
        su2double Total_Viscosity = (nodes->GetLaminarViscosity(iPoint) +
                                     nodes->GetEddyViscosity(iPoint));
        su2double AuxVar = 0.0;
        if (yCoord > EPS)
          AuxVar = Total_Viscosity*yVelocity/yCoord;

        /*--- Set the auxilairy variable for this node. ---*/

        nodes->SetAuxVar(iPoint, 0, AuxVar);

      }

      /*--- Compute the auxiliary variable gradient with GG or WLS. ---*/

      if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
        SetAuxVar_Gradient_GG(geometry, config);
      }
      if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
        SetAuxVar_Gradient_LS(geometry, config);
      }

    }

    /*--- loop over points ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Conservative variables w/o reconstruction ---*/

      numerics->SetPrimitive(nodes->GetPrimitive(iPoint), nullptr);

      /*--- Set incompressible density  ---*/

      numerics->SetDensity(nodes->GetDensity(iPoint),
                           nodes->GetDensity(iPoint));

      /*--- Set control volume ---*/

      numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Set y coordinate ---*/

      numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                         geometry->nodes->GetCoord(iPoint));

      /*--- If viscous, we need gradients for extra terms. ---*/

      if (viscous) {

        /*--- Gradient of the primitive variables ---*/

        numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint), nullptr);

        /*--- Load the aux variable gradient that we already computed. ---*/

        numerics->SetAuxVarGrad(nodes->GetAuxVarGradient(iPoint), nullptr);

      }

      /*--- Compute Source term Residual ---*/

      auto residual = numerics->ComputeResidual(config);

      /*--- Add Residual ---*/

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Implicit part ---*/

      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }

  if (radiation) {

    CNumerics* second_numerics = numerics_container[SOURCE_SECOND_TERM];

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Store the radiation source term ---*/

      second_numerics->SetRadVarSource(solver_container[RAD_SOL]->GetNodes()->GetRadiative_SourceTerm(iPoint));

      /*--- Set control volume ---*/

      second_numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Compute the residual ---*/

      auto residual = second_numerics->ComputeResidual(config);

      /*--- Add Residual ---*/

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Implicit part ---*/

      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      if (vol_heat) {

        if(solver_container[RAD_SOL]->GetNodes()->GetVol_HeatSource(iPoint)) {

          auto Volume = geometry->nodes->GetVolume(iPoint);

          /*--- Subtract integrated source from the residual. ---*/
          LinSysRes(iPoint, nDim+1) -= config->GetHeatSource_Val()*Volume;
        }

      }

    }

  }

  /*--- Check if a verification solution is to be computed. ---*/

  if (VerificationSolution) {
    if ( VerificationSolution->IsManufacturedSolution() ) {

      /*--- Get the physical time. ---*/
      su2double time = 0.0;
      if (config->GetTime_Marching()) time = config->GetPhysicalTime();

      /*--- Loop over points ---*/
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

        /*--- Get control volume size. ---*/
        su2double Volume = geometry->nodes->GetVolume(iPoint);

        /*--- Get the current point coordinates. ---*/
        const su2double *coor = geometry->nodes->GetCoord(iPoint);

        /*--- Get the MMS source term. ---*/
        vector<su2double> sourceMan(nVar,0.0);
        VerificationSolution->GetMMSSourceTerm(coor, time, sourceMan.data());

        /*--- Compute the residual for this control volume and subtract. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          LinSysRes[iPoint*nVar+iVar] -= sourceMan[iVar]*Volume;
        }

      }
    }
  }

}

void CIncEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {

  /* This method should be used to call any new source terms for a particular problem*/
  /* This method calls the new child class in CNumerics, where the new source term should be implemented.  */

  /* Next we describe how to get access to some important quanties for this method */
  /* Access to all points in the current geometric mesh by saying: nPointDomain */
  /* Get the vector of conservative variables at some point iPoint = nodes->GetSolution(iPoint) */
  /* Get the volume (or area in 2D) associated with iPoint = nodes->GetVolume(iPoint) */
  /* Get the vector of geometric coordinates of point iPoint = nodes->GetCoord(iPoint) */

}

void CIncEulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {

  su2double Area, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0,
  Mean_BetaInc2, Lambda, ProjVel, ProjVel_i, ProjVel_j, *GridVel, *GridVel_i, *GridVel_j;
  const su2double* Normal;

  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;

  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    nodes->SetLambda(iPoint,0.0);
  }

  /*--- Loop interior edges ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Point identification, Normal vector and area ---*/

    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    Normal = geometry->edges->GetNormal(iEdge);
    Area = GeometryToolbox::Norm(nDim, Normal);

    /*--- Mean Values ---*/

    Mean_ProjVel    = 0.5 * (nodes->GetProjVel(iPoint,Normal) + nodes->GetProjVel(jPoint,Normal));
    Mean_BetaInc2   = 0.5 * (nodes->GetBetaInc2(iPoint)      + nodes->GetBetaInc2(jPoint));
    Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);

    /*--- Adjustment for grid movement ---*/

    if (dynamic_grid) {
      GridVel_i = geometry->nodes->GetGridVel(iPoint);
      GridVel_j = geometry->nodes->GetGridVel(jPoint);
      ProjVel_i = 0.0; ProjVel_j =0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }

    /*--- Inviscid contribution ---*/

    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->nodes->GetDomain(iPoint)) nodes->AddLambda(iPoint,Lambda);
    if (geometry->nodes->GetDomain(jPoint)) nodes->AddLambda(jPoint,Lambda);

  }

  /*--- Loop boundary edges ---*/

  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

      /*--- Point identification, Normal vector and area ---*/

      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = GeometryToolbox::Norm(nDim, Normal);

      /*--- Mean Values ---*/

      Mean_ProjVel    = nodes->GetProjVel(iPoint,Normal);
      Mean_BetaInc2   = nodes->GetBetaInc2(iPoint);
      Mean_SoundSpeed = sqrt(Mean_BetaInc2*Area*Area);

      /*--- Adjustment for grid movement ---*/

      if (dynamic_grid) {
        GridVel = geometry->nodes->GetGridVel(iPoint);
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }

      /*--- Inviscid contribution ---*/

      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->nodes->GetDomain(iPoint)) {
        nodes->AddLambda(iPoint,Lambda);
      }

    }
    }
  }

  /*--- Correct the eigenvalue values across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_MAX_EIG);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_MAX_EIG);
  }

  /*--- MPI parallelization ---*/

  InitiateComms(geometry, config, MAX_EIGENVALUE);
  CompleteComms(geometry, config, MAX_EIGENVALUE);

}

void CIncEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint, jPoint, iEdge;
  su2double *Diff;
  unsigned short iVar;
  bool boundary_i, boundary_j;

  Diff = new su2double[nVar];

  nodes->SetUnd_LaplZero();

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    /*--- Solution differences ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      Diff[iVar] = nodes->GetSolution(iPoint,iVar) - nodes->GetSolution(jPoint,iVar);

    boundary_i = geometry->nodes->GetPhysicalBoundary(iPoint);
    boundary_j = geometry->nodes->GetPhysicalBoundary(jPoint);

    /*--- Both points inside the domain, or both in the boundary ---*/

    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->nodes->GetDomain(iPoint)) nodes->SubtractUnd_Lapl(iPoint,Diff);
      if (geometry->nodes->GetDomain(jPoint)) nodes->AddUnd_Lapl(jPoint,Diff);
    }

    /*--- iPoint inside the domain, jPoint on the boundary ---*/

    if (!boundary_i && boundary_j)
      if (geometry->nodes->GetDomain(iPoint)) nodes->SubtractUnd_Lapl(iPoint,Diff);

    /*--- jPoint inside the domain, iPoint on the boundary ---*/

    if (boundary_i && !boundary_j)
      if (geometry->nodes->GetDomain(jPoint)) nodes->AddUnd_Lapl(jPoint,Diff);

  }

  /*--- Correct the Laplacian values across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_LAPLACIAN);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_LAPLACIAN);
  }

  /*--- MPI parallelization ---*/

  InitiateComms(geometry, config, UNDIVIDED_LAPLACIAN);
  CompleteComms(geometry, config, UNDIVIDED_LAPLACIAN);

  delete [] Diff;

}

void CIncEulerSolver::SetCentered_Dissipation_Sensor(CGeometry *geometry, CConfig *config) {

  unsigned long iEdge, iPoint, jPoint;
  su2double Pressure_i = 0.0, Pressure_j = 0.0;
  bool boundary_i, boundary_j;

  /*--- Reset variables to store the undivided pressure ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    iPoint_UndLapl[iPoint] = 0.0;
    jPoint_UndLapl[iPoint] = 0.0;
  }

  /*--- Evaluate the pressure sensor ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    /*--- Get the pressure, or density for incompressible solvers ---*/

    Pressure_i = nodes->GetDensity(iPoint);
    Pressure_j = nodes->GetDensity(jPoint);

    boundary_i = geometry->nodes->GetPhysicalBoundary(iPoint);
    boundary_j = geometry->nodes->GetPhysicalBoundary(jPoint);

    /*--- Both points inside the domain, or both on the boundary ---*/

    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {

      if (geometry->nodes->GetDomain(iPoint)) {
        iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i);
        jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j);
      }

      if (geometry->nodes->GetDomain(jPoint)) {
        iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j);
        jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j);
      }

    }

    /*--- iPoint inside the domain, jPoint on the boundary ---*/

    if (!boundary_i && boundary_j)
      if (geometry->nodes->GetDomain(iPoint)) {
        iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i);
        jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j);
      }

    /*--- jPoint inside the domain, iPoint on the boundary ---*/

    if (boundary_i && !boundary_j)
      if (geometry->nodes->GetDomain(jPoint)) {
        iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j);
        jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j);
      }

  }

  /*--- Correct the sensor values across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_SENSOR);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_SENSOR);
  }

  /*--- Set pressure switch for each point ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    nodes->SetSensor(iPoint,fabs(iPoint_UndLapl[iPoint]) / jPoint_UndLapl[iPoint]);

  /*--- MPI parallelization ---*/

  InitiateComms(geometry, config, SENSOR);
  CompleteComms(geometry, config, SENSOR);

}

void CIncEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {

  su2double *Residual, *Res_TruncError, Vol, Delta, Res;
  unsigned short iVar, jVar;
  unsigned long iPoint;

  su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  bool adjoint = config->GetContinuous_Adjoint();

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = (geometry->nodes->GetVolume(iPoint) +
           geometry->nodes->GetPeriodicVolume(iPoint));
    Delta = nodes->GetDelta_Time(iPoint) / Vol;

    Res_TruncError = nodes->GetResTruncError(iPoint);
    Residual = LinSysRes.GetBlock(iPoint);

    if (!adjoint) {
      SetPreconditioner(config, iPoint);
      for (iVar = 0; iVar < nVar; iVar ++ ) {
        Res = 0.0;
        for (jVar = 0; jVar < nVar; jVar ++ )
          Res += Preconditioner[iVar][jVar]*(Residual[jVar] + Res_TruncError[jVar]);
        nodes->AddSolution(iPoint,iVar, -Res*Delta*RK_AlphaCoeff);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint));
      }
    }
  }

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

  /*--- For verification cases, compute the global error metrics. ---*/

  ComputeVerificationError(geometry, config);

}

void CIncEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  unsigned short iVar, jVar;
  unsigned long iPoint;

  bool adjoint = config->GetContinuous_Adjoint();

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = (geometry->nodes->GetVolume(iPoint) +
           geometry->nodes->GetPeriodicVolume(iPoint));
    Delta = nodes->GetDelta_Time(iPoint) / Vol;

    local_Res_TruncError = nodes->GetResTruncError(iPoint);
    local_Residual = LinSysRes.GetBlock(iPoint);


    if (!adjoint) {
      SetPreconditioner(config, iPoint);
      for (iVar = 0; iVar < nVar; iVar ++ ) {
        Res = 0.0;
        for (jVar = 0; jVar < nVar; jVar ++ )
          Res += Preconditioner[iVar][jVar]*(local_Residual[jVar] + local_Res_TruncError[jVar]);
        nodes->AddSolution(iPoint,iVar, -Res*Delta);
        AddRes_RMS(iVar, Res*Res);
        AddRes_Max(iVar, fabs(Res), geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint));
      }
    }
  }

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

  /*--- For verification cases, compute the global error metrics. ---*/

  ComputeVerificationError(geometry, config);

}

void CIncEulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iVar, jVar;
  unsigned long iPoint, total_index, IterLinSol = 0;
  su2double Delta, *local_Res_TruncError, Vol;

  bool adjoint = config->GetContinuous_Adjoint();

  /*--- Set maximum residual to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Build implicit system ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Read the residual ---*/

    local_Res_TruncError = nodes->GetResTruncError(iPoint);

    /*--- Read the volume ---*/

    Vol = (geometry->nodes->GetVolume(iPoint) +
           geometry->nodes->GetPeriodicVolume(iPoint));

    /*--- Apply the preconditioner and add to the diagonal. ---*/

    if (nodes->GetDelta_Time(iPoint) != 0.0) {
      Delta = Vol / nodes->GetDelta_Time(iPoint);
      SetPreconditioner(config, iPoint);
      for (iVar = 0; iVar < nVar; iVar ++ ) {
        for (jVar = 0; jVar < nVar; jVar ++ ) {
          Preconditioner[iVar][jVar] = Delta*Preconditioner[iVar][jVar];
        }
      }
      Jacobian.AddBlock2Diag(iPoint, Preconditioner);
    } else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/

    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint));
    }

  }

  /*--- Initialize residual and solution at the ghost points ---*/

  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }

  /*--- Solve or smooth the linear system ---*/

  IterLinSol = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  /*--- Store the value of the residual. ---*/

  SetResLinSolver(System.GetResidual());

  /*--- The the number of iterations of the linear solver ---*/

  SetIterLinSolver(IterLinSol);

  /*--- Update solution (system written in terms of increments) ---*/

  if (!adjoint) {
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        nodes->AddSolution(iPoint, iVar, nodes->GetUnderRelaxation(iPoint)*LinSysSol[iPoint*nVar+iVar]);
      }
    }
  }

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

  /*--- For verification cases, compute the global error metrics. ---*/

  ComputeVerificationError(geometry, config);

}

void CIncEulerSolver::SetBeta_Parameter(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh) {

  su2double epsilon2  = config->GetBeta_Factor();
  su2double epsilon2_default = 4.1;
  su2double maxVel2 = 0.0;
  su2double Beta = 1.0;

  unsigned long iPoint;

  /*--- For now, only the finest mesh level stores the Beta for all levels. ---*/

  if (iMesh == MESH_0) {

    for (iPoint = 0; iPoint < nPoint; iPoint++) {

      /*--- Store the local maximum of the squared velocity in the field. ---*/

      if (nodes->GetVelocity2(iPoint) > maxVel2)
        maxVel2 = nodes->GetVelocity2(iPoint);

    }

    /*--- Communicate the max globally to give a conservative estimate. ---*/

#ifdef HAVE_MPI
    su2double myMaxVel2 = maxVel2; maxVel2 = 0.0;
    SU2_MPI::Allreduce(&myMaxVel2, &maxVel2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    Beta = max(1e-10,maxVel2);
    config->SetMax_Vel2(Beta);

  }

  /*--- Allow an override if user supplies a large epsilon^2. ---*/

  epsilon2 = max(epsilon2_default,epsilon2);

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    nodes->SetBetaInc2(iPoint,epsilon2*config->GetMax_Vel2());

}

void CIncEulerSolver::SetPreconditioner(CConfig *config, unsigned long iPoint) {

  unsigned short iDim, jDim;

  su2double  BetaInc2, Density, dRhodT, Temperature, oneOverCp, Cp;
  su2double  Velocity[3] = {0.0,0.0,0.0};

  bool variable_density = (config->GetKind_DensityModel() == VARIABLE);
  bool implicit         = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool energy           = config->GetEnergy_Equation();

  /*--- Access the primitive variables at this node. ---*/

  Density     = nodes->GetDensity(iPoint);
  BetaInc2    = nodes->GetBetaInc2(iPoint);
  Cp          = nodes->GetSpecificHeatCp(iPoint);
  oneOverCp   = 1.0/Cp;
  Temperature = nodes->GetTemperature(iPoint);

  for (iDim = 0; iDim < nDim; iDim++)
    Velocity[iDim] = nodes->GetVelocity(iPoint,iDim);

  /*--- We need the derivative of the equation of state to build the
   preconditioning matrix. For now, the only option is the ideal gas
   law, but in the future, dRhodT should be in the fluid model. ---*/

  if (variable_density) {
    dRhodT = -Density/Temperature;
  } else {
    dRhodT = 0.0;
  }

  /*--- Calculating the inverse of the preconditioning matrix
   that multiplies the time derivative during time integration. ---*/

  if (implicit) {

    /*--- For implicit calculations, we multiply the preconditioner
     by the cell volume over the time step and add to the Jac diagonal. ---*/

    Preconditioner[0][0] = 1.0/BetaInc2;
    for (iDim = 0; iDim < nDim; iDim++)
      Preconditioner[iDim+1][0] = Velocity[iDim]/BetaInc2;

    if (energy) Preconditioner[nDim+1][0] = Cp*Temperature/BetaInc2;
    else        Preconditioner[nDim+1][0] = 0.0;

    for (jDim = 0; jDim < nDim; jDim++) {
      Preconditioner[0][jDim+1] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        if (iDim == jDim) Preconditioner[iDim+1][jDim+1] = Density;
        else Preconditioner[iDim+1][jDim+1] = 0.0;
      }
      Preconditioner[nDim+1][jDim+1] = 0.0;
    }

    Preconditioner[0][nDim+1] = dRhodT;
    for (iDim = 0; iDim < nDim; iDim++)
      Preconditioner[iDim+1][nDim+1] = Velocity[iDim]*dRhodT;

    if (energy) Preconditioner[nDim+1][nDim+1] = Cp*(dRhodT*Temperature + Density);
    else        Preconditioner[nDim+1][nDim+1] = 1.0;

  } else {

    /*--- For explicit calculations, we move the residual to the
     right-hand side and pre-multiply by the preconditioner inverse.
     Therefore, we build inv(Precon) here and multiply by the residual
     later in the R-K and Euler Explicit time integration schemes. ---*/

    Preconditioner[0][0] = Temperature*BetaInc2*dRhodT/Density + BetaInc2;
    for (iDim = 0; iDim < nDim; iDim ++)
      Preconditioner[iDim+1][0] = -1.0*Velocity[iDim]/Density;

    if (energy) Preconditioner[nDim+1][0] = -1.0*Temperature/Density;
    else        Preconditioner[nDim+1][0] = 0.0;


    for (jDim = 0; jDim < nDim; jDim++) {
      Preconditioner[0][jDim+1] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        if (iDim == jDim) Preconditioner[iDim+1][jDim+1] = 1.0/Density;
        else Preconditioner[iDim+1][jDim+1] = 0.0;
      }
      Preconditioner[nDim+1][jDim+1] = 0.0;
    }

    Preconditioner[0][nDim+1] = -1.0*BetaInc2*dRhodT*oneOverCp/Density;
    for (iDim = 0; iDim < nDim; iDim ++)
      Preconditioner[iDim+1][nDim+1] = 0.0;

    if (energy) Preconditioner[nDim+1][nDim+1] = oneOverCp/Density;
    else        Preconditioner[nDim+1][nDim+1] = 0.0;

  }

}

void CIncEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;

  su2double *V_infty, *V_domain;

  bool implicit      = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;
  bool viscous       = config->GetViscous();

  su2double *Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Allocate the value at the infinity ---*/

    V_infty = GetCharacPrimVar(val_marker, iVertex);

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Index of the closest interior node ---*/

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Recompute and store the velocity in the primitive variable vector. ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        V_infty[iDim+1] = GetVelocity_Inf(iDim);

      /*--- Far-field pressure set to static pressure (0.0). ---*/

      V_infty[0] = GetPressure_Inf();

      /*--- Dirichlet condition for temperature at far-field (if energy is active). ---*/

      V_infty[nDim+1] = GetTemperature_Inf();

      /*--- Store the density.  ---*/

      V_infty[nDim+2] = GetDensity_Inf();

      /*--- Beta coefficient stored at the node ---*/

      V_infty[nDim+3] = nodes->GetBetaInc2(iPoint);

      /*--- Cp is needed for Temperature equation. ---*/

      V_infty[nDim+7] = nodes->GetSpecificHeatCp(iPoint);

      /*--- Set various quantities in the numerics class ---*/

      conv_numerics->SetPrimitive(V_domain, V_infty);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the convective residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Update residual value ---*/

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Convective Jacobian contribution for implicit integration ---*/

      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      /*--- Viscous residual contribution ---*/

      if (viscous) {

        /*--- Set transport properties at infinity. ---*/

        V_infty[nDim+4] = nodes->GetLaminarViscosity(iPoint);
        V_infty[nDim+5] = nodes->GetEddyViscosity(iPoint);
        V_infty[nDim+6] = nodes->GetThermalConductivity(iPoint);

        /*--- Set the normal vector and the coordinates ---*/

        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                                geometry->nodes->GetCoord(Point_Normal));

        /*--- Primitive variables, and gradient ---*/

        visc_numerics->SetPrimitive(V_domain, V_infty);
        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                          nodes->GetGradient_Primitive(iPoint));

        /*--- Turbulent kinetic energy ---*/

        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));

        /*--- Compute and update viscous residual ---*/

        auto residual = visc_numerics->ComputeResidual(config);
        LinSysRes.SubtractBlock(iPoint, residual);

        /*--- Viscous Jacobian contribution for implicit integration ---*/

        if (implicit)
          Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

      }

    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] Normal;

}

void CIncEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  unsigned long Point_Normal;
  su2double *Flow_Dir, Flow_Dir_Mag, Vel_Mag, Area, P_total, P_domain, Vn;
  su2double *V_inlet, *V_domain;
  su2double UnitFlowDir[3] = {0.0,0.0,0.0};
  su2double dV[3] = {0.0,0.0,0.0};
  su2double Damping = config->GetInc_Inlet_Damping();

  bool implicit      = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool viscous       = config->GetViscous();

  string Marker_Tag  = config->GetMarker_All_TagBound(val_marker);

  unsigned short Kind_Inlet = config->GetKind_Inc_Inlet(Marker_Tag);

  su2double *Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Allocate the value at the inlet ---*/

    V_inlet = GetCharacPrimVar(val_marker, iVertex);

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Index of the closest interior node ---*/

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = GeometryToolbox::Norm(nDim, Normal);

      /*--- Both types of inlets may use the prescribed flow direction.
       Ensure that the flow direction is a unit vector. ---*/

      Flow_Dir = Inlet_FlowDir[val_marker][iVertex];
      Flow_Dir_Mag = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Flow_Dir_Mag += Flow_Dir[iDim]*Flow_Dir[iDim];
      Flow_Dir_Mag = sqrt(Flow_Dir_Mag);

      /*--- Store the unit flow direction vector. ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        UnitFlowDir[iDim] = Flow_Dir[iDim]/Flow_Dir_Mag;

      /*--- Retrieve solution at this boundary node. ---*/

      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Neumann condition for dynamic pressure ---*/

      V_inlet[0] = nodes->GetPressure(iPoint);

      /*--- The velocity is either prescribed or computed from total pressure. ---*/

      switch (Kind_Inlet) {

          /*--- Velocity and temperature (if required) been specified at the inlet. ---*/

        case VELOCITY_INLET:

          /*--- Retrieve the specified velocity and temperature for the inlet. ---*/

          Vel_Mag  = Inlet_Ptotal[val_marker][iVertex]/config->GetVelocity_Ref();

          /*--- Store the velocity in the primitive variable vector. ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Vel_Mag*UnitFlowDir[iDim];

          /*--- Dirichlet condition for temperature (if energy is active) ---*/

          V_inlet[nDim+1] = Inlet_Ttotal[val_marker][iVertex]/config->GetTemperature_Ref();

          break;

          /*--- Stagnation pressure has been specified at the inlet. ---*/

        case PRESSURE_INLET:

          /*--- Retrieve the specified total pressure for the inlet. ---*/

          P_total = Inlet_Ptotal[val_marker][iVertex]/config->GetPressure_Ref();

          /*--- Store the current static pressure for clarity. ---*/

          P_domain = nodes->GetPressure(iPoint);

          /*--- Check for back flow through the inlet. ---*/

          Vn = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Vn += V_domain[iDim+1]*(-1.0*Normal[iDim]/Area);
          }

          /*--- If the local static pressure is larger than the specified
           total pressure or the velocity is directed upstream, we have a
           back flow situation. The specified total pressure should be used
           as a static pressure condition and the velocity from the domain
           is used for the BC. ---*/

          if ((P_domain > P_total) || (Vn < 0.0)) {

            /*--- Back flow: use the prescribed P_total as static pressure. ---*/

            V_inlet[0] = Inlet_Ptotal[val_marker][iVertex]/config->GetPressure_Ref();

            /*--- Neumann condition for velocity. ---*/

            for (iDim = 0; iDim < nDim; iDim++)
              V_inlet[iDim+1] = V_domain[iDim+1];

            /*--- Neumann condition for the temperature. ---*/

            V_inlet[nDim+1] = nodes->GetTemperature(iPoint);

          } else {

            /*--- Update the velocity magnitude using the total pressure. ---*/

            Vel_Mag = sqrt((P_total - P_domain)/(0.5*nodes->GetDensity(iPoint)));

            /*--- If requested, use the local boundary normal (negative),
             instead of the prescribed flow direction in the config. ---*/

            if (config->GetInc_Inlet_UseNormal()) {
              for (iDim = 0; iDim < nDim; iDim++)
                UnitFlowDir[iDim] = -Normal[iDim]/Area;
            }

            /*--- Compute the delta change in velocity in each direction. ---*/

            for (iDim = 0; iDim < nDim; iDim++)
              dV[iDim] = Vel_Mag*UnitFlowDir[iDim] - V_domain[iDim+1];

            /*--- Update the velocity in the primitive variable vector.
             Note we use damping here to improve stability/convergence. ---*/

            for (iDim = 0; iDim < nDim; iDim++)
              V_inlet[iDim+1] = V_domain[iDim+1] + Damping*dV[iDim];

            /*--- Dirichlet condition for temperature (if energy is active) ---*/

            V_inlet[nDim+1] = Inlet_Ttotal[val_marker][iVertex]/config->GetTemperature_Ref();

          }

          break;

      }

      /*--- Access density at the node. This is either constant by
        construction, or will be set fixed implicitly by the temperature
        and equation of state. ---*/

      V_inlet[nDim+2] = nodes->GetDensity(iPoint);

      /*--- Beta coefficient from the config file ---*/

      V_inlet[nDim+3] = nodes->GetBetaInc2(iPoint);

      /*--- Cp is needed for Temperature equation. ---*/

      V_inlet[nDim+7] = nodes->GetSpecificHeatCp(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Update residual value ---*/

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

      /*--- Viscous contribution, commented out because serious convergence problems ---*/

      if (viscous) {

        /*--- Set transport properties at the inlet ---*/

        V_inlet[nDim+4] = nodes->GetLaminarViscosity(iPoint);
        V_inlet[nDim+5] = nodes->GetEddyViscosity(iPoint);
        V_inlet[nDim+6] = nodes->GetThermalConductivity(iPoint);

        /*--- Set the normal vector and the coordinates ---*/

        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                                geometry->nodes->GetCoord(Point_Normal));

        /*--- Primitive variables, and gradient ---*/

        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                          nodes->GetGradient_Primitive(iPoint));

        /*--- Turbulent kinetic energy ---*/

        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));

        /*--- Compute and update residual ---*/

        auto residual = visc_numerics->ComputeResidual(config);

        LinSysRes.SubtractBlock(iPoint, residual);

        /*--- Jacobian contribution for implicit integration ---*/

        if (implicit)
          Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

      }

    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] Normal;

}

void CIncEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_outlet, *V_domain, P_Outlet = 0.0, P_domain;
  su2double mDot_Target, mDot_Old, dP, Density_Avg, Area_Outlet;
  su2double Damping = config->GetInc_Outlet_Damping();

  bool implicit      = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool viscous       = config->GetViscous();
  string Marker_Tag  = config->GetMarker_All_TagBound(val_marker);

  su2double *Normal = new su2double[nDim];

  unsigned short Kind_Outlet = config->GetKind_Inc_Outlet(Marker_Tag);

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Allocate the value at the outlet ---*/

    V_outlet = GetCharacPrimVar(val_marker, iVertex);

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Index of the closest interior node ---*/

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Current solution at this boundary node ---*/

      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Store the current static pressure for clarity. ---*/

      P_domain = nodes->GetPressure(iPoint);

      /*--- Compute a boundary value for the pressure depending on whether
       we are prescribing a back pressure or a mass flow target. ---*/

      switch (Kind_Outlet) {

          /*--- Velocity and temperature (if required) been specified at the inlet. ---*/

        case PRESSURE_OUTLET:

          /*--- Retrieve the specified back pressure for this outlet. ---*/

          P_Outlet = config->GetOutlet_Pressure(Marker_Tag)/config->GetPressure_Ref();

          /*--- The pressure is prescribed at the outlet. ---*/

          V_outlet[0] = P_Outlet;

          /*--- Neumann condition for the velocity. ---*/

          for (iDim = 0; iDim < nDim; iDim++) {
            V_outlet[iDim+1] = nodes->GetVelocity(iPoint,iDim);
          }

          break;

          /*--- A mass flow target has been specified for the outlet. ---*/

        case MASS_FLOW_OUTLET:

          /*--- Retrieve the specified target mass flow at the outlet. ---*/

          mDot_Target = config->GetOutlet_Pressure(Marker_Tag)/(config->GetDensity_Ref() * config->GetVelocity_Ref());

          /*--- Retrieve the old mass flow, density, and area of the outlet,
           which has been computed in a preprocessing step. These values
           were stored in non-dim. form in the config container. ---*/

          mDot_Old    = config->GetOutlet_MassFlow(Marker_Tag);
          Density_Avg = config->GetOutlet_Density(Marker_Tag);
          Area_Outlet = config->GetOutlet_Area(Marker_Tag);

          /*--- Compute the pressure increment based on the difference
           between the current and target mass flow. Note that increasing
           pressure decreases flow speed. ---*/

          dP = 0.5*Density_Avg*(mDot_Old*mDot_Old - mDot_Target*mDot_Target)/((Density_Avg*Area_Outlet)*(Density_Avg*Area_Outlet));

          /*--- Update the new outlet pressure. Note that we use damping
           here to improve stability/convergence. ---*/

          P_Outlet = P_domain + Damping*dP;

          /*--- The pressure is prescribed at the outlet. ---*/

          V_outlet[0] = P_Outlet;

          /*--- Neumann condition for the velocity ---*/

          for (iDim = 0; iDim < nDim; iDim++) {
            V_outlet[iDim+1] = nodes->GetVelocity(iPoint,iDim);
          }

          break;

      }

      /*--- Neumann condition for the temperature. ---*/

      V_outlet[nDim+1] = nodes->GetTemperature(iPoint);

      /*--- Access density at the interior node. This is either constant by
        construction, or will be set fixed implicitly by the temperature
        and equation of state. ---*/

      V_outlet[nDim+2] = nodes->GetDensity(iPoint);

      /*--- Beta coefficient from the config file ---*/

      V_outlet[nDim+3] = nodes->GetBetaInc2(iPoint);

      /*--- Cp is needed for Temperature equation. ---*/

      V_outlet[nDim+7] = nodes->GetSpecificHeatCp(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_outlet);

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Update residual value ---*/

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit) {
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
      }

      /*--- Viscous contribution, commented out because serious convergence problems ---*/

      if (viscous) {

        /*--- Set transport properties at the outlet. ---*/

        V_outlet[nDim+4] = nodes->GetLaminarViscosity(iPoint);
        V_outlet[nDim+5] = nodes->GetEddyViscosity(iPoint);
        V_outlet[nDim+6] = nodes->GetThermalConductivity(iPoint);

        /*--- Set the normal vector and the coordinates ---*/

        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                                geometry->nodes->GetCoord(Point_Normal));

        /*--- Primitive variables, and gradient ---*/

        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                          nodes->GetGradient_Primitive(iPoint));

        /*--- Turbulent kinetic energy ---*/

        if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST))
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                              solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));

        /*--- Compute and update residual ---*/

        auto residual = visc_numerics->ComputeResidual(config);

        LinSysRes.SubtractBlock(iPoint, residual);

        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

      }

    }
  }

  /*--- Free locally allocated memory ---*/
  delete [] Normal;

}

void CIncEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                        unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {

  /*--- Local variables ---*/

  unsigned short iVar, jVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;

  su2double Density, Cp;
  su2double *V_time_nM1, *V_time_n, *V_time_nP1;
  su2double U_time_nM1[5], U_time_n[5], U_time_nP1[5];
  su2double Volume_nM1, Volume_nP1, TimeStep;
  su2double *GridVel_i = nullptr, *GridVel_j = nullptr, Residual_GCL;
  const su2double* Normal;

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool energy   = config->GetEnergy_Equation();

  /*--- Store the physical time step ---*/

  TimeStep = config->GetDelta_UnstTimeND();

  /*--- Compute the dual time-stepping source term for static meshes ---*/

  if (!dynamic_grid) {

    /*--- Loop over all nodes (excluding halos) ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Initialize the Residual / Jacobian container to zero. ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
        if (implicit) {
        for (jVar = 0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;
        }
      }

      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. These are actually
       the primitive values, but we will convert to conservatives. ---*/

      V_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      V_time_n   = nodes->GetSolution_time_n(iPoint);
      V_time_nP1 = nodes->GetSolution(iPoint);

      /*--- Access the density and Cp at this node (constant for now). ---*/

      Density     = nodes->GetDensity(iPoint);
      Cp          = nodes->GetSpecificHeatCp(iPoint);

      /*--- Compute the conservative variable vector for all time levels. ---*/

      U_time_nM1[0] = Density;
      U_time_n[0]   = Density;
      U_time_nP1[0] = Density;

      for (iDim = 0; iDim < nDim; iDim++) {
        U_time_nM1[iDim+1] = Density*V_time_nM1[iDim+1];
        U_time_n[iDim+1]   = Density*V_time_n[iDim+1];
        U_time_nP1[iDim+1] = Density*V_time_nP1[iDim+1];
      }

      U_time_nM1[nDim+1] = Density*Cp*V_time_nM1[nDim+1];
      U_time_n[nDim+1]   = Density*Cp*V_time_n[nDim+1];
      U_time_nP1[nDim+1] = Density*Cp*V_time_nP1[nDim+1];

      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/

      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order). Note that for an
       incompressible problem, the pressure equation does not have a
       contribution, as the time derivative should always be zero. ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetTime_Marching() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (config->GetTime_Marching() == DT_STEPPING_2ND)
          Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                            +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }

      if (!energy) Residual[nDim+1] = 0.0;

      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/

      LinSysRes.AddBlock(iPoint, Residual);

      if (implicit) {
        for (iVar = 1; iVar < nVar; iVar++) {
          if (config->GetTime_Marching() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetTime_Marching() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_i[iDim+1][iDim+1] = Density*Jacobian_i[iDim+1][iDim+1];
        if (energy) Jacobian_i[nDim+1][nDim+1] = Density*Cp*Jacobian_i[nDim+1][nDim+1];

        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
      }

    }

  }

  else {

    /*--- For unsteady flows on dynamic meshes (rigidly transforming or
     dynamically deforming), the Geometric Conservation Law (GCL) should be
     satisfied in conjunction with the ALE formulation of the governing
     equations. The GCL prevents accuracy issues caused by grid motion, i.e.
     a uniform free-stream should be preserved through a moving grid. First,
     we will loop over the edges and boundaries to compute the GCL component
     of the dual time source term that depends on grid velocities. ---*/

    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

      /*--- Initialize the Residual / Jacobian container to zero. ---*/

      for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;

      /*--- Get indices for nodes i & j plus the face normal ---*/

      iPoint = geometry->edges->GetNode(iEdge,0);
      jPoint = geometry->edges->GetNode(iEdge,1);
      Normal = geometry->edges->GetNormal(iEdge);

      /*--- Grid velocities stored at nodes i & j ---*/

      GridVel_i = geometry->nodes->GetGridVel(iPoint);
      GridVel_j = geometry->nodes->GetGridVel(jPoint);

      /*--- Compute the GCL term by averaging the grid velocities at the
       edge mid-point and dotting with the face normal. ---*/

      Residual_GCL = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual_GCL += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

      /*--- Compute the GCL component of the source term for node i ---*/

      V_time_n = nodes->GetSolution_time_n(iPoint);

      /*--- Access the density and Cp at this node (constant for now). ---*/

      Density     = nodes->GetDensity(iPoint);
      Cp          = nodes->GetSpecificHeatCp(iPoint);

      /*--- Compute the conservative variable vector for all time levels. ---*/

      U_time_n[0] = Density;
      for (iDim = 0; iDim < nDim; iDim++) {
        U_time_n[iDim+1] = Density*V_time_n[iDim+1];
      }
      U_time_n[nDim+1] = Density*Cp*V_time_n[nDim+1];

      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;

      if (!energy) Residual[nDim+1] = 0.0;
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Compute the GCL component of the source term for node j ---*/

      V_time_n = nodes->GetSolution_time_n(jPoint);

      U_time_n[0] = Density;
      for (iDim = 0; iDim < nDim; iDim++) {
        U_time_n[iDim+1] = Density*V_time_n[iDim+1];
      }
      U_time_n[nDim+1] = Density*Cp*V_time_n[nDim+1];

      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;

      if (!energy) Residual[nDim+1] = 0.0;
      LinSysRes.SubtractBlock(jPoint, Residual);

    }

    /*---  Loop over the boundary edges ---*/

    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
          (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        /*--- Initialize the Residual / Jacobian container to zero. ---*/

        for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;

        /*--- Get the index for node i plus the boundary face normal ---*/

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

        /*--- Grid velocities stored at boundary node i ---*/

        GridVel_i = geometry->nodes->GetGridVel(iPoint);

        /*--- Compute the GCL term by dotting the grid velocity with the face
         normal. The normal is negated to match the boundary convention. ---*/

        Residual_GCL = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Residual_GCL -= 0.5*(GridVel_i[iDim]+GridVel_i[iDim])*Normal[iDim];

        /*--- Compute the GCL component of the source term for node i ---*/

        V_time_n = nodes->GetSolution_time_n(iPoint);

        /*--- Access the density and Cp at this node (constant for now). ---*/

        Density     = nodes->GetDensity(iPoint);
        Cp          = nodes->GetSpecificHeatCp(iPoint);

        U_time_n[0] = Density;
        for (iDim = 0; iDim < nDim; iDim++) {
          U_time_n[iDim+1] = Density*V_time_n[iDim+1];
        }
        U_time_n[nDim+1] = Density*Cp*V_time_n[nDim+1];

        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = U_time_n[iVar]*Residual_GCL;

        if (!energy) Residual[nDim+1] = 0.0;
        LinSysRes.AddBlock(iPoint, Residual);

      }
      }
    }

    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Initialize the Residual / Jacobian container to zero. ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
      }

      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/

      V_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      V_time_n   = nodes->GetSolution_time_n(iPoint);
      V_time_nP1 = nodes->GetSolution(iPoint);

      /*--- Access the density and Cp at this node (constant for now). ---*/

      Density     = nodes->GetDensity(iPoint);
      Cp          = nodes->GetSpecificHeatCp(iPoint);

      /*--- Compute the conservative variable vector for all time levels. ---*/

      U_time_nM1[0] = Density;
      U_time_n[0]   = Density;
      U_time_nP1[0] = Density;

      for (iDim = 0; iDim < nDim; iDim++) {
        U_time_nM1[iDim+1] = Density*V_time_nM1[iDim+1];
        U_time_n[iDim+1]   = Density*V_time_n[iDim+1];
        U_time_nP1[iDim+1] = Density*V_time_nP1[iDim+1];
      }

      U_time_nM1[nDim+1] = Density*Cp*V_time_nM1[nDim+1];
      U_time_n[nDim+1]   = Density*Cp*V_time_n[nDim+1];
      U_time_nP1[nDim+1] = Density*Cp*V_time_nP1[nDim+1];

      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/

      Volume_nM1 = geometry->nodes->GetVolume_nM1(iPoint);
      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetTime_Marching() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
        if (config->GetTime_Marching() == DT_STEPPING_2ND)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
          + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      }

      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      if (!energy) Residual[nDim+1] = 0.0;
      LinSysRes.AddBlock(iPoint, Residual);

      if (implicit) {
        for (iVar = 1; iVar < nVar; iVar++) {
          if (config->GetTime_Marching() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetTime_Marching() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_i[iDim+1][iDim+1] = Density*Jacobian_i[iDim+1][iDim+1];
        if (energy) Jacobian_i[nDim+1][nDim+1] = Density*Cp*Jacobian_i[nDim+1][nDim+1];

        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
      }

    }
  }

}

void CIncEulerSolver::GetOutlet_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) {

  unsigned short iDim, iMarker;
  unsigned long iVertex, iPoint;
  su2double *V_outlet = nullptr, Velocity[3], MassFlow,
  Velocity2, Density, Area, AxiFactor;
  unsigned short iMarker_Outlet, nMarker_Outlet;
  string Inlet_TagBound, Outlet_TagBound;

  bool axisymmetric = config->GetAxisymmetric();

  bool write_heads = ((((config->GetInnerIter() % (config->GetScreen_Wrt_Freq(2)*40)) == 0)
                       && (config->GetInnerIter()!= 0))
                      || (config->GetInnerIter() == 1));

  /*--- Get the number of outlet markers and check for any mass flow BCs. ---*/

  nMarker_Outlet = config->GetnMarker_Outlet();
  bool Evaluate_BC = false;
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
    Outlet_TagBound = config->GetMarker_Outlet_TagBound(iMarker_Outlet);
    if (config->GetKind_Inc_Outlet(Outlet_TagBound) == MASS_FLOW_OUTLET)
      Evaluate_BC = true;
  }

  /*--- If we have a massflow outlet BC, then we need to compute and
   communicate the total massflow, density, and area through each outlet
   boundary, so that it can be used in the iterative procedure to update
   the back pressure until we converge to the desired mass flow. This
   routine is called only once per iteration as a preprocessing and the
   values for all outlets are stored and retrieved later in the BC_Outlet
   routines. ---*/

  if (Evaluate_BC) {

    su2double *Outlet_MassFlow = new su2double[config->GetnMarker_All()];
    su2double *Outlet_Density  = new su2double[config->GetnMarker_All()];
    su2double *Outlet_Area     = new su2double[config->GetnMarker_All()];

    /*--- Comute MassFlow, average temp, press, etc. ---*/

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

      Outlet_MassFlow[iMarker] = 0.0;
      Outlet_Density[iMarker]  = 0.0;
      Outlet_Area[iMarker]     = 0.0;

      if ((config->GetMarker_All_KindBC(iMarker) == OUTLET_FLOW) ) {

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->nodes->GetDomain(iPoint)) {

            V_outlet = nodes->GetPrimitive(iPoint);

            geometry->vertex[iMarker][iVertex]->GetNormal(Vector);

            if (axisymmetric) {
              if (geometry->nodes->GetCoord(iPoint, 1) != 0.0)
                AxiFactor = 2.0*PI_NUMBER*geometry->nodes->GetCoord(iPoint, 1);
              else
                AxiFactor = 1.0;
            } else {
              AxiFactor = 1.0;
            }

            Density      = V_outlet[nDim+2];

            Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0;

            for (iDim = 0; iDim < nDim; iDim++) {
              Area += (Vector[iDim] * AxiFactor) * (Vector[iDim] * AxiFactor);
              Velocity[iDim] = V_outlet[iDim+1];
              Velocity2 += Velocity[iDim] * Velocity[iDim];
              MassFlow += Vector[iDim] * AxiFactor * Density * Velocity[iDim];
            }
            Area = sqrt (Area);

            Outlet_MassFlow[iMarker] += MassFlow;
            Outlet_Density[iMarker]  += Density*Area;
            Outlet_Area[iMarker]     += Area;

          }
        }
      }
    }

    /*--- Copy to the appropriate structure ---*/

    su2double *Outlet_MassFlow_Local = new su2double[nMarker_Outlet];
    su2double *Outlet_Density_Local  = new su2double[nMarker_Outlet];
    su2double *Outlet_Area_Local     = new su2double[nMarker_Outlet];

    su2double *Outlet_MassFlow_Total = new su2double[nMarker_Outlet];
    su2double *Outlet_Density_Total  = new su2double[nMarker_Outlet];
    su2double *Outlet_Area_Total     = new su2double[nMarker_Outlet];

    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      Outlet_MassFlow_Local[iMarker_Outlet] = 0.0;
      Outlet_Density_Local[iMarker_Outlet]  = 0.0;
      Outlet_Area_Local[iMarker_Outlet]     = 0.0;

      Outlet_MassFlow_Total[iMarker_Outlet] = 0.0;
      Outlet_Density_Total[iMarker_Outlet]  = 0.0;
      Outlet_Area_Total[iMarker_Outlet]     = 0.0;
    }

    /*--- Copy the values to the local array for MPI ---*/

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == OUTLET_FLOW)) {
        for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
          Outlet_TagBound = config->GetMarker_Outlet_TagBound(iMarker_Outlet);
          if (config->GetMarker_All_TagBound(iMarker) == Outlet_TagBound) {
            Outlet_MassFlow_Local[iMarker_Outlet] += Outlet_MassFlow[iMarker];
            Outlet_Density_Local[iMarker_Outlet]  += Outlet_Density[iMarker];
            Outlet_Area_Local[iMarker_Outlet]     += Outlet_Area[iMarker];
          }
        }
      }
    }

    /*--- All the ranks to compute the total value ---*/

#ifdef HAVE_MPI

    SU2_MPI::Allreduce(Outlet_MassFlow_Local, Outlet_MassFlow_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Density_Local, Outlet_Density_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Outlet_Area_Local, Outlet_Area_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else

    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      Outlet_MassFlow_Total[iMarker_Outlet] = Outlet_MassFlow_Local[iMarker_Outlet];
      Outlet_Density_Total[iMarker_Outlet]  = Outlet_Density_Local[iMarker_Outlet];
      Outlet_Area_Total[iMarker_Outlet]     = Outlet_Area_Local[iMarker_Outlet];
    }

#endif

    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      if (Outlet_Area_Total[iMarker_Outlet] != 0.0) {
        Outlet_Density_Total[iMarker_Outlet] /= Outlet_Area_Total[iMarker_Outlet];
      }
      else {
        Outlet_Density_Total[iMarker_Outlet] = 0.0;
      }

      if (iMesh == MESH_0) {
        config->SetOutlet_MassFlow(iMarker_Outlet, Outlet_MassFlow_Total[iMarker_Outlet]);
        config->SetOutlet_Density(iMarker_Outlet, Outlet_Density_Total[iMarker_Outlet]);
        config->SetOutlet_Area(iMarker_Outlet, Outlet_Area_Total[iMarker_Outlet]);
      }
    }

    /*--- Screen output using the values already stored in the config container ---*/

    if ((rank == MASTER_NODE) && (iMesh == MESH_0) ) {

      cout.precision(5);
      cout.setf(ios::fixed, ios::floatfield);

      if (write_heads && Output && !config->GetDiscrete_Adjoint()) {
        cout << endl   << "---------------------------- Outlet properties --------------------------" << endl;
      }

      for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
        Outlet_TagBound = config->GetMarker_Outlet_TagBound(iMarker_Outlet);
        if (write_heads && Output && !config->GetDiscrete_Adjoint()) {

          /*--- Geometry defintion ---*/

          cout <<"Outlet surface: " << Outlet_TagBound << "." << endl;

          if ((nDim ==3) || axisymmetric) {
            cout <<"Area (m^2): " << config->GetOutlet_Area(Outlet_TagBound) << endl;
          }
          if (nDim == 2) {
            cout <<"Length (m): " << config->GetOutlet_Area(Outlet_TagBound) << "." << endl;
          }

          cout << setprecision(5) << "Outlet Avg. Density (kg/m^3): " <<  config->GetOutlet_Density(Outlet_TagBound) * config->GetDensity_Ref() << endl;
          su2double Outlet_mDot = fabs(config->GetOutlet_MassFlow(Outlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
          cout << "Outlet mass flow (kg/s): "; cout << setprecision(5) << Outlet_mDot;

        }
      }

      if (write_heads && Output && !config->GetDiscrete_Adjoint()) {cout << endl;
        cout << "-------------------------------------------------------------------------" << endl << endl;
      }

      cout.unsetf(ios_base::floatfield);

    }

    delete [] Outlet_MassFlow_Local;
    delete [] Outlet_Density_Local;
    delete [] Outlet_Area_Local;

    delete [] Outlet_MassFlow_Total;
    delete [] Outlet_Density_Total;
    delete [] Outlet_Area_Total;

    delete [] Outlet_MassFlow;
    delete [] Outlet_Density;
    delete [] Outlet_Area;

  }

}

void CIncEulerSolver::PrintVerificationError(const CConfig *config) const {

  if ((rank != MASTER_NODE) || (MGLevel != MESH_0)) return;

  if (config && !config->GetDiscrete_Adjoint()) {

    cout.precision(6);
    cout.setf(ios::scientific, ios::floatfield);

    cout << endl   << "------------------------ Global Error Analysis --------------------------" << endl;

    cout << setw(20) << "RMS Error [P]: " << setw(12) << VerificationSolution->GetError_RMS(0) << "     | ";
    cout << setw(20) << "Max Error [P]: " << setw(12) << VerificationSolution->GetError_Max(0);
    cout << endl;

    cout << setw(20) << "RMS Error [U]: " << setw(12) << VerificationSolution->GetError_RMS(1) << "     | ";
    cout << setw(20) << "Max Error [U]: " << setw(12) << VerificationSolution->GetError_Max(1);
    cout << endl;

    cout << setw(20) << "RMS Error [V]: " << setw(12) << VerificationSolution->GetError_RMS(2) << "     | ";
    cout << setw(20) << "Max Error [V]: " << setw(12) << VerificationSolution->GetError_Max(2);
    cout << endl;

    if (nDim == 3) {
      cout << setw(20) << "RMS Error [W]: " << setw(12) << VerificationSolution->GetError_RMS(3) << "     | ";
      cout << setw(20) << "Max Error [W]: " << setw(12) << VerificationSolution->GetError_Max(3);
      cout << endl;
    }

    if (config->GetEnergy_Equation()) {
      cout << setw(20) << "RMS Error [T]: " << setw(12) << VerificationSolution->GetError_RMS(nDim+1) << "     | ";
      cout << setw(20) << "Max Error [T]: " << setw(12) << VerificationSolution->GetError_Max(nDim+1);
      cout << endl;
    }

    cout << "-------------------------------------------------------------------------" << endl << endl;
    cout.unsetf(ios_base::floatfield);
  }
}

void CIncEulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh, iMeshFine;
  unsigned long iPoint, index, iChildren, Point_Fine;
  unsigned short turb_model = config->GetKind_Turb_Model();
  su2double Area_Children, Area_Parent, Coord[3] = {0.0}, *Solution_Fine;
  bool static_fsi = ((config->GetTime_Marching() == STEADY) && config->GetFSI_Simulation());
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool steady_restart = config->GetSteadyRestart();
  bool turbulent = (config->GetKind_Solver() == INC_RANS) || (config->GetKind_Solver() == DISC_ADJ_INC_RANS);

  string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Store the number of variables for the turbulence model
   (that could appear in the restart file before the grid velocities). ---*/
  unsigned short turbVars = 0;
  if (turbulent){
    if ((turb_model == SST) || (turb_model == SST_SUST)) turbVars = 2;
    else turbVars = 1;
  }

  /*--- Adjust the number of solution variables in the restart. We always
   carry a space in nVar for the energy equation in the solver, but we only
   write it to the restart if it is active. Therefore, we must reduce nVar
   here if energy is inactive so that the restart is read correctly. ---*/

  bool energy               = config->GetEnergy_Equation();
  bool weakly_coupled_heat  = config->GetWeakly_Coupled_Heat();

  unsigned short nVar_Restart = nVar;
  if ((!energy) && (!weakly_coupled_heat)) nVar_Restart--;
  Solution[nVar-1] = GetTemperature_Inf();

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar_Restart; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      nodes->SetSolution(iPoint_Local,Solution);
      iPoint_Global_Local++;

      /*--- For dynamic meshes, read in and store the
       grid coordinates and grid velocities for each node. ---*/

      if (dynamic_grid && val_update_geo) {

        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        /*--- If we are restarting the solution from a previously computed static calculation (no grid movement) ---*/
        /*--- the grid velocities are set to 0. This is useful for FSI computations ---*/

        /*--- Rewind the index to retrieve the Coords. ---*/
        index = counter*Restart_Vars[1];
        for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim]; }

        su2double GridVel[3] = {0.0,0.0,0.0};
        if (!steady_restart) {
          /*--- Move the index forward to get the grid velocities. ---*/
          index = counter*Restart_Vars[1] + skipVars + nVar_Restart + turbVars;
          for (iDim = 0; iDim < nDim; iDim++) { GridVel[iDim] = Restart_Data[index+iDim]; }
        }

        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->nodes->SetCoord(iPoint_Local, iDim, Coord[iDim]);
          geometry[MESH_0]->nodes->SetGridVel(iPoint_Local, iDim, GridVel[iDim]);
        }
      }

      /*--- For static FSI problems, grid_movement is 0 but we need to read in and store the
       grid coordinates for each node (but not the grid velocities, as there are none). ---*/

      if (static_fsi && val_update_geo) {
       /*--- Rewind the index to retrieve the Coords. ---*/
        index = counter*Restart_Vars[1];
        for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim];}

        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->nodes->SetCoord(iPoint_Local, iDim, Coord[iDim]);
        }
      }

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;

    }
  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < nPointDomain) {
    SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Update the geometry for flows on deforming meshes ---*/

  if ((dynamic_grid || static_fsi) && val_update_geo) {

    /*--- Communicate the new coordinates and grid velocities at the halos ---*/

    geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, COORDINATES);
    geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, COORDINATES);

    if (dynamic_grid) {
      geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, GRID_VELOCITY);
      geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, GRID_VELOCITY);
    }

    /*--- Recompute the edges and  dual mesh control volumes in the
     domain and on the boundaries. ---*/

    geometry[MESH_0]->SetControlVolume(config, UPDATE);
    geometry[MESH_0]->SetBoundControlVolume(config, UPDATE);
    geometry[MESH_0]->SetMaxLength(config);

    /*--- Update the multigrid structure after setting up the finest grid,
     including computing the grid velocities on the coarser levels. ---*/

    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      iMeshFine = iMesh-1;
      geometry[iMesh]->SetControlVolume(config, geometry[iMeshFine], UPDATE);
      geometry[iMesh]->SetBoundControlVolume(config, geometry[iMeshFine],UPDATE);
      geometry[iMesh]->SetCoord(geometry[iMeshFine]);
      if (dynamic_grid) {
        geometry[iMesh]->SetRestricted_GridVelocity(geometry[iMeshFine], config);
      }
      geometry[iMesh]->SetMaxLength(config);
    }
  }

  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. We alo call the preprocessing routine
   on the fine level in order to have all necessary quantities updated,
   especially if this is a turbulent simulation (eddy viscosity). ---*/

  solver[MESH_0][FLOW_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][FLOW_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);

  /*--- For turbulent simulations the flow preprocessing is done by the turbulence solver
   *    after it loads its variables (they are needed to compute flow primitives). ---*/
  if (!turbulent) {
    solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  }

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->nodes->GetVolume(iPoint);
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->nodes->GetnChildren_CV(iPoint); iChildren++) {
        Point_Fine = geometry[iMesh]->nodes->GetChildren_CV(iPoint, iChildren);
        Area_Children = geometry[iMesh-1]->nodes->GetVolume(Point_Fine);
        Solution_Fine = solver[iMesh-1][FLOW_SOL]->GetNodes()->GetSolution(Point_Fine);
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][FLOW_SOL]->GetNodes()->SetSolution(iPoint,Solution);
    }
    solver[iMesh][FLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][FLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);

    if (!turbulent) {
      solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    }
  }

  /*--- Update the old geometry (coordinates n and n-1) in dual time-stepping strategy ---*/
  if (dual_time && config->GetGrid_Movement() && !config->GetDeform_Mesh() &&
      (config->GetKind_GridMovement() != RIGID_MOTION)) {
    Restart_OldGeometry(geometry[MESH_0], config);
  }

  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars; Restart_Vars = nullptr;
  delete [] Restart_Data; Restart_Data = nullptr;

}

void CIncEulerSolver::SetFreeStream_Solution(CConfig *config){

  unsigned long iPoint;
  unsigned short iDim;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    nodes->SetSolution(iPoint,0, Pressure_Inf);
    for (iDim = 0; iDim < nDim; iDim++){
      nodes->SetSolution(iPoint,iDim+1, Velocity_Inf[iDim]);
    }
    nodes->SetSolution(iPoint,nDim+1, Temperature_Inf);
  }
}
