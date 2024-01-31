/*!
 * \file CIncEulerSolver.cpp
 * \brief Main subroutines for solving incompressible flow (Euler, Navier-Stokes, etc.).
 * \author F. Palacios, T. Economon
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

#include "../../include/solvers/CIncEulerSolver.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/fluid/CConstantDensity.hpp"
#include "../../include/fluid/CIncIdealGas.hpp"
#include "../../include/fluid/CIncIdealGasPolynomial.hpp"
#include "../../include/variables/CIncNSVariable.hpp"
#include "../../include/limiters/CLimiterDetails.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/fluid/CFluidScalar.hpp"
#include "../../include/fluid/CFluidFlamelet.hpp"
#include "../../include/fluid/CFluidModel.hpp"


CIncEulerSolver::CIncEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh,
                                 const bool navier_stokes) :
  CFVMFlowSolverBase<CIncEulerVariable, ENUM_REGIME::INCOMPRESSIBLE>(*geometry, *config) {

  /*--- Based on the navier_stokes boolean, determine if this constructor is
   *    being called by itself, or by its derived class CIncNSSolver. ---*/
  const string description = navier_stokes? "Navier-Stokes" : "Euler";

  unsigned short iMarker;
  ifstream restart_file;
  unsigned short nZone = geometry->GetnZone();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  int Unst_RestartIter = 0;
  unsigned short iZone = config->GetiZone();
  bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));
  bool time_stepping = config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING;
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (restart && (iMesh == MESH_0) && nZone <= 1) {

    /*--- Multizone problems require the number of the zone to be appended. ---*/

    auto filename_ = config->GetSolution_FileName();

    if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone, ".dat");

    /*--- Modify file name for a dual-time unsteady restart ---*/

    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST)
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

    filename_ = "flow";
    filename_ = config->GetFilename(filename_, ".meta", Unst_RestartIter);
    Read_SU2_Restart_Metadata(geometry, config, adjoint, filename_);

  }
  if (restart && (config->GetKind_Streamwise_Periodic() == ENUM_STREAMWISE_PERIODIC::MASSFLOW)) {
    string filename_ = "flow";
    filename_ = config->GetFilename(filename_, ".meta", Unst_RestartIter);
    Read_SU2_Restart_Metadata(geometry, config, adjoint, filename_);
    if (rank==MASTER_NODE) cout << "Setting streamwise periodic pressure drop from restart metadata file." << endl;
  }

  /*--- Set the gamma value ---*/

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometry constants in the solver structure.
   * Incompressible flow, primitive variables (P, vx, vy, vz, T, rho, beta, lamMu, EddyMu, Kt_eff, Cp, Cv) ---*/

  nDim = geometry->GetnDim();

  /*--- Make sure to align the sizes with the constructor of CIncEulerVariable. ---*/
  nVar = nDim+2; nPrimVar = nDim+9; nPrimVarGrad = nDim+4;

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nPrimVarGrad;

  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Store the number of vertices on each marker for deallocation later ---*/

  nVertex.resize(nMarker);
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

  /*--- Allocate base class members. ---*/

  Allocate(*config);

  /*--- MPI + OpenMP initialization. ---*/

  HybridParallelInitialization(*config, *geometry);

  /*--- Jacobians and vector structures for implicit computations ---*/

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {

    if (rank == MASTER_NODE)
      cout << "Initialize Jacobian structure (" << description << "). MG level: " << iMesh <<"." << endl;

    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
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

  SetReferenceValues(*config);

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  if (navier_stokes) {
    nodes = new CIncNSVariable(Pressure_Inf, Velocity_Inf, Temperature_Inf, nPoint, nDim, nVar, config);
  } else {
    nodes = new CIncEulerVariable(Pressure_Inf, Velocity_Inf, Temperature_Inf, nPoint, nDim, nVar, config);
  }
  SetBaseClassPointerToNodes();

  if (iMesh == MESH_0) {
    nodes->NonPhysicalEdgeCounter.resize(geometry->GetnEdge()) = 0;
  }

  /*--- Initial comms. ---*/

  CommunicateInitialState(geometry, config);

  /*--- Sizing edge mass flux array ---*/
  if (config->GetBounded_Scalar())
    EdgeMassFluxes.resize(geometry->GetnEdge()) = su2double(0.0);

  /*--- Add the solver name. ---*/
  SolverName = "INC.FLOW";

  /*--- Finally, check that the static arrays will be large enough (keep this
   *    check at the bottom to make sure we consider the "final" values). ---*/
  if((nDim > MAXNDIM) || (nPrimVar > MAXNVAR))
    SU2_MPI::Error("Oops! The CIncEulerSolver static array sizes are not large enough.", CURRENT_FUNCTION);
}

CIncEulerSolver::~CIncEulerSolver() {

  for(auto& model : FluidModel) delete model;
}

void CIncEulerSolver::SetNondimensionalization(CConfig *config, unsigned short iMesh) {

  su2double Temperature_FreeStream = 0.0,  ModVel_FreeStream = 0.0,Energy_FreeStream = 0.0,
  ModVel_FreeStreamND = 0.0, Omega_FreeStream = 0.0, Omega_FreeStreamND = 0.0, Viscosity_FreeStream = 0.0,
  Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, Pressure_Thermodynamic = 0.0, Tke_FreeStream = 0.0,
  Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Temperature_Ref = 0.0, Velocity_Ref = 0.0, Time_Ref = 0.0,
  Gas_Constant_Ref = 0.0, Omega_Ref = 0.0, Force_Ref = 0.0, Viscosity_Ref = 0.0, Conductivity_Ref = 0.0, Heat_Flux_Ref = 0.0, Energy_Ref= 0.0, Pressure_FreeStreamND = 0.0, Pressure_ThermodynamicND = 0.0, Density_FreeStreamND = 0.0,
  Temperature_FreeStreamND = 0.0, Gas_ConstantND = 0.0, Specific_Heat_CpND = 0.0, Thermal_Expansion_CoeffND = 0.0,
  Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStreamND = 0.0,
  Tke_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0,
  Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0;

  unsigned short iDim, iVar;

  /*--- Local variables ---*/

  su2double Mach     = config->GetMach();
  su2double Reynolds = config->GetReynolds();

  bool unsteady      = (config->GetTime_Marching() != TIME_MARCHING::STEADY);
  bool viscous       = config->GetViscous();
  bool turbulent     = ((config->GetKind_Solver() == MAIN_SOLVER::INC_RANS) ||
                        (config->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_INC_RANS));
  bool tkeNeeded     = ((turbulent) && ((config->GetKind_Turb_Model() == TURB_MODEL::SST)));
  bool energy        = config->GetEnergy_Equation();
  bool boussinesq    = (config->GetKind_DensityModel() == INC_DENSITYMODEL::BOUSSINESQ);

  /*--- Compute dimensional free-stream values. ---*/

  Density_FreeStream     = config->GetInc_Density_Init();     config->SetDensity_FreeStream(Density_FreeStream);
  Temperature_FreeStream = config->GetInc_Temperature_Init(); config->SetTemperature_FreeStream(Temperature_FreeStream);
  Pressure_FreeStream    = 0.0; config->SetPressure_FreeStream(Pressure_FreeStream);

  /*--- The dimensional viscosity is needed to determine the free-stream conditions.
    To accomplish this, simply set the non-dimensional coefficients to the
    dimensional ones. This will be overruled later.---*/

  config->SetTemperature_Ref(1.0);
  config->SetViscosity_Ref(1.0);
  config->SetConductivity_Ref(1.0);
  config->SetGas_Constant_Ref(1.0);

  ModVel_FreeStream   = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ModVel_FreeStream += config->GetInc_Velocity_Init()[iDim]*config->GetInc_Velocity_Init()[iDim];
    config->SetVelocity_FreeStream(config->GetInc_Velocity_Init()[iDim],iDim);
  }
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);

  CFluidModel* auxFluidModel = nullptr;

  switch (config->GetKind_FluidModel()) {

    case CONSTANT_DENSITY:

      auxFluidModel = new CConstantDensity(Density_FreeStream, config->GetSpecific_Heat_Cp());
      auxFluidModel->SetTDState_T(Temperature_FreeStream);
      break;

    case INC_IDEAL_GAS:

      config->SetGas_Constant(UNIVERSAL_GAS_CONSTANT/(config->GetMolecular_Weight()/1000.0));
      Pressure_Thermodynamic = Density_FreeStream*Temperature_FreeStream*config->GetGas_Constant();
      auxFluidModel = new CIncIdealGas(config->GetSpecific_Heat_Cp(), config->GetGas_Constant(), Pressure_Thermodynamic);
      auxFluidModel->SetTDState_T(Temperature_FreeStream);
      Pressure_Thermodynamic = auxFluidModel->GetPressure();
      config->SetPressure_Thermodynamic(Pressure_Thermodynamic);
      break;

    case INC_IDEAL_GAS_POLY:

      config->SetGas_Constant(UNIVERSAL_GAS_CONSTANT/(config->GetMolecular_Weight()/1000.0));
      Pressure_Thermodynamic = Density_FreeStream*Temperature_FreeStream*config->GetGas_Constant();
      auxFluidModel = new CIncIdealGasPolynomial<N_POLY_COEFFS>(config->GetGas_Constant(), Pressure_Thermodynamic);
      if (viscous) {
        /*--- Variable Cp model via polynomial. ---*/
        for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++)
          config->SetCp_PolyCoeffND(config->GetCp_PolyCoeff(iVar), iVar);
        auxFluidModel->SetCpModel(config);
      }
      auxFluidModel->SetTDState_T(Temperature_FreeStream);
      Pressure_Thermodynamic = auxFluidModel->GetPressure();
      config->SetPressure_Thermodynamic(Pressure_Thermodynamic);
      break;

    case FLUID_MIXTURE:

      config->SetGas_Constant(UNIVERSAL_GAS_CONSTANT / (config->GetMolecular_Weight() / 1000.0));
      Pressure_Thermodynamic = config->GetPressure_Thermodynamic();
      auxFluidModel = new CFluidScalar(config->GetSpecific_Heat_Cp(), config->GetGas_Constant(), Pressure_Thermodynamic, config);
      auxFluidModel->SetTDState_T(Temperature_FreeStream, config->GetSpecies_Init());
      break;

    case FLUID_FLAMELET:

      config->SetGas_Constant(UNIVERSAL_GAS_CONSTANT / (config->GetMolecular_Weight() / 1000.0));
      Pressure_Thermodynamic = config->GetPressure_Thermodynamic();
      auxFluidModel = new CFluidFlamelet(config, Pressure_Thermodynamic);
      config->SetPressure_Thermodynamic(Pressure_Thermodynamic);
      auxFluidModel->SetTDState_T(Temperature_FreeStream, config->GetSpecies_Init());
      break;

    default:

      SU2_MPI::Error("Fluid model not implemented for incompressible solver.", CURRENT_FUNCTION);
      break;
  }

  if (viscous) {

    for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++)
      config->SetMu_PolyCoeffND(config->GetMu_PolyCoeff(iVar), iVar);

    /*--- Use the fluid model to compute the dimensional viscosity/conductivity. ---*/

    auxFluidModel->SetLaminarViscosityModel(config);
    Viscosity_FreeStream = auxFluidModel->GetLaminarViscosity();
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

  Energy_FreeStream = auxFluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;
  if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; };
  config->SetEnergy_FreeStream(Energy_FreeStream);

  /*--- Auxilary (dimensional) FluidModel no longer needed. ---*/
  delete auxFluidModel;

  /*--- Compute Mach number ---*/

  if (config->GetKind_FluidModel() == CONSTANT_DENSITY) {
    Mach = ModVel_FreeStream / sqrt(config->GetBulk_Modulus()/Density_FreeStream);
  } else {
    Mach = 0.0;
  }
  config->SetMach(Mach);

  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/

  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();       config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Pressure_ThermodynamicND = Pressure_Thermodynamic/config->GetPressure_Ref(); config->SetPressure_ThermodynamicND(Pressure_ThermodynamicND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();         config->SetDensity_FreeStreamND(Density_FreeStreamND);

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }

  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);
  Gas_ConstantND      = config->GetGas_Constant()/Gas_Constant_Ref;               config->SetGas_ConstantND(Gas_ConstantND);
  Specific_Heat_CpND  = config->GetSpecific_Heat_CpND();

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

  const su2double MassDiffusivityND = config->GetDiffusivity_Constant() / (Velocity_Ref * Length_Ref);
  config->SetDiffusivity_ConstantND(MassDiffusivityND);

  /*--- Create one final fluid model object per OpenMP thread to be able to use them in parallel.
   *    GetFluidModel() should be used to automatically access the "right" object of each thread. ---*/

  assert(FluidModel.empty() && "Potential memory leak!");
  FluidModel.resize(omp_get_max_threads());

  for (auto& fluidModel : FluidModel) {

    switch (config->GetKind_FluidModel()) {

      case CONSTANT_DENSITY:
        fluidModel = new CConstantDensity(Density_FreeStreamND, Specific_Heat_CpND);
        break;

      case INC_IDEAL_GAS:
        fluidModel = new CIncIdealGas(Specific_Heat_CpND, Gas_ConstantND, Pressure_ThermodynamicND);
        fluidModel->SetTDState_T(Temperature_FreeStreamND);
        break;

      case FLUID_MIXTURE:
        fluidModel = new CFluidScalar(Specific_Heat_CpND, Gas_ConstantND, Pressure_ThermodynamicND, config);
        fluidModel->SetTDState_T(Temperature_FreeStreamND, config->GetSpecies_Init());
        break;

      case FLUID_FLAMELET:
        fluidModel = new CFluidFlamelet(config, Pressure_Thermodynamic);
        fluidModel->SetTDState_T(Temperature_FreeStreamND, config->GetSpecies_Init());
        break;

      case INC_IDEAL_GAS_POLY:
        fluidModel = new CIncIdealGasPolynomial<N_POLY_COEFFS>(Gas_ConstantND, Pressure_ThermodynamicND);
        if (viscous) {
          /*--- Variable Cp model via polynomial. ---*/
          config->SetCp_PolyCoeffND(config->GetCp_PolyCoeff(0)/Gas_Constant_Ref, 0);
          for (iVar = 1; iVar < config->GetnPolyCoeffs(); iVar++)
            config->SetCp_PolyCoeffND(config->GetCp_PolyCoeff(iVar)*pow(Temperature_Ref,iVar)/Gas_Constant_Ref, iVar);
          fluidModel->SetCpModel(config);
        }
        fluidModel->SetTDState_T(Temperature_FreeStreamND);
        break;

    }

    if (viscous) {

      /*--- Viscosity model via polynomial. ---*/

      config->SetMu_PolyCoeffND(config->GetMu_PolyCoeff(0)/Viscosity_Ref, 0);
      for (iVar = 1; iVar < config->GetnPolyCoeffs(); iVar++)
        config->SetMu_PolyCoeffND(config->GetMu_PolyCoeff(iVar)*pow(Temperature_Ref,iVar)/Viscosity_Ref, iVar);

      /*--- Conductivity model via polynomial. ---*/

      config->SetKt_PolyCoeffND(config->GetKt_PolyCoeff(0)/Conductivity_Ref, 0);
      for (iVar = 1; iVar < config->GetnPolyCoeffs(); iVar++)
        config->SetKt_PolyCoeffND(config->GetKt_PolyCoeff(iVar)*pow(Temperature_Ref,iVar)/Conductivity_Ref, iVar);

      /*--- Set up the transport property models. ---*/

      fluidModel->SetLaminarViscosityModel(config);
      fluidModel->SetThermalConductivityModel(config);
      fluidModel->SetMassDiffusivityModel(config);
    }

  }

  Energy_FreeStreamND = GetFluidModel()->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;

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

      case INC_DENSITYMODEL::CONSTANT:
        if (energy) cout << "Energy equation is active and decoupled." << endl;
        else cout << "No energy equation." << endl;
        break;

      case INC_DENSITYMODEL::BOUSSINESQ:
        if (energy) cout << "Energy equation is active and coupled through Boussinesq approx." << endl;
        break;

      case INC_DENSITYMODEL::VARIABLE:
        if (energy) cout << "Energy equation is active and coupled for variable density." << endl;
        break;

      case INC_DENSITYMODEL::FLAMELET:
        cout << "Energy equation is disabled and density is obtained through flamelet manifold." << endl;
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
      case VISCOSITYMODEL::CONSTANT:
        ModelTable << "CONSTANT_VISCOSITY";
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Viscosity" << config->GetMu_Constant() << config->GetMu_Constant()/config->GetMu_ConstantND() << Unit.str() << config->GetMu_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case VISCOSITYMODEL::FLAMELET:
        ModelTable << "FLAMELET";
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Viscosity" << "--" << "--" << Unit.str() << config->GetMu_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case VISCOSITYMODEL::COOLPROP:
        ModelTable << "COOLPROP_VISCOSITY";
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Viscosity" << "--" << "--" << Unit.str() << config->GetMu_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case VISCOSITYMODEL::SUTHERLAND:
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

      case VISCOSITYMODEL::POLYNOMIAL:
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
      case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
        ModelTable << "CONSTANT_PRANDTL";
        NonDimTable << "Prandtl (Lam.)"  << "-" << "-" << "-" << config->GetPrandtl_Lam();
        Unit.str("");
        NonDimTable << "Prandtl (Turb.)" << "-" << "-" << "-" << config->GetPrandtl_Turb();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case CONDUCTIVITYMODEL::CONSTANT:
        ModelTable << "CONSTANT";
        Unit << "W/m^2.K";
        NonDimTable << "Molecular Cond." << config->GetThermal_Conductivity_Constant() << config->GetThermal_Conductivity_Constant()/config->GetThermal_Conductivity_ConstantND() << Unit.str() << config->GetThermal_Conductivity_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case CONDUCTIVITYMODEL::FLAMELET:
        ModelTable << "FLAMELET";
        Unit << "W/m^2.K";
        NonDimTable << "Molecular Cond." << "--" << "--" << Unit.str() << config->GetThermal_Conductivity_ConstantND();
        Unit.str("");
        break;

      case CONDUCTIVITYMODEL::COOLPROP:
        ModelTable << "COOLPROP";
        Unit << "W/m^2.K";
        NonDimTable << "Molecular Cond." << "--" << "--" << Unit.str() << config->GetThermal_Conductivity_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case CONDUCTIVITYMODEL::POLYNOMIAL:
        ModelTable << "POLYNOMIAL";
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

    case FLUID_MIXTURE:
      ModelTable << "FLUID_MIXTURE";
      Unit << "N.m/kg.K";
      NonDimTable << "Spec. Heat (Cp)" << config->GetSpecific_Heat_Cp() << config->GetSpecific_Heat_Cp() / config->GetSpecific_Heat_CpND() << Unit.str() << config->GetSpecific_Heat_CpND();
      Unit.str("");
      Unit << "g/mol";
      NonDimTable << "Molecular weight" << config->GetMolecular_Weight() << 1.0 << Unit.str() << config->GetMolecular_Weight();
      Unit.str("");
      Unit << "N.m/kg.K";
      NonDimTable << "Gas Constant" << config->GetGas_Constant() << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();
      Unit.str("");
      Unit << "Pa";
      NonDimTable << "Therm. Pressure" << config->GetPressure_Thermodynamic() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_ThermodynamicND();
      Unit.str("");
      NonDimTable.PrintFooter();
      break;

    case FLUID_FLAMELET:
      ModelTable << "FLAMELET";
      Unit << "N.m/kg.K";
      NonDimTable << "Spec. Heat (Cp)" << "--" << "--" << Unit.str() << config->GetSpecific_Heat_CpND();
      Unit.str("");
      Unit << "g/mol";
      NonDimTable << "Molecular weight" << "--" << "--" << Unit.str() << config->GetMolecular_Weight();
      Unit.str("");
      Unit << "N.m/kg.K";
      NonDimTable << "Gas Constant" << "--" << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();
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
      NonDimTable << "Conductivity" << "-" << config->GetThermal_Conductivity_Ref() << Unit.str() << "-";
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
      if (config->GetKind_Species_Model() != SPECIES_MODEL::NONE) {
        if      (config->GetSystemMeasurements() == SI) Unit << "m^2/s";
        else if (config->GetSystemMeasurements() == US) Unit << "ft^2/s";
        NonDimTable << "Mass Diffusivity" << config->GetDiffusivity_Constant() << config->GetDiffusivity_Constant()/config->GetDiffusivity_ConstantND() << Unit.str() << config->GetDiffusivity_ConstantND();
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

void CIncEulerSolver::SetReferenceValues(const CConfig& config) {

  /*--- Evaluate reference values for non-dimensionalization. For dimensional or non-dim
   based on initial values, use the far-field state (inf). For a custom non-dim based
   on user-provided reference values, use the ref values to compute the forces. ---*/

  su2double RefDensity, RefVel2;

  if ((config.GetRef_Inc_NonDim() == DIMENSIONAL) ||
      (config.GetRef_Inc_NonDim() == INITIAL_VALUES)) {
    RefDensity = Density_Inf;
    RefVel2 = GeometryToolbox::SquaredNorm(nDim, Velocity_Inf);
  }
  else {
    RefDensity = config.GetInc_Density_Ref();
    RefVel2 = pow(config.GetInc_Velocity_Ref(), 2);
  }

  DynamicPressureRef = 0.5 * RefDensity * RefVel2;
  AeroCoeffForceRef =  DynamicPressureRef * config.GetRefArea();

}

void CIncEulerSolver::CommonPreprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                          unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  const bool implicit   = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool center     = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  const bool center_jst = (config->GetKind_Centered_Flow() == CENTERED::JST) && (iMesh == MESH_0);
  const bool outlet     = (config->GetnMarker_Outlet() != 0);

  /*--- Set the primitive variables ---*/

  ompMasterAssignBarrier(ErrorCounter, 0);

  SU2_OMP_ATOMIC
  ErrorCounter += SetPrimitive_Variables(solver_container, config);

  if ((iMesh == MESH_0) && (config->GetComm_Level() == COMM_FULL)) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    {
      unsigned long tmp = ErrorCounter;
      SU2_MPI::Allreduce(&tmp, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      config->SetNonphysical_Points(ErrorCounter);
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Artificial dissipation ---*/

  if (center && !Output) {
    SetMax_Eigenvalue(geometry, config);
    if (center_jst) {
      SetCentered_Dissipation_Sensor(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }

  /*--- Update the beta value based on the maximum velocity. ---*/

  SetBeta_Parameter(geometry, solver_container, config, iMesh);

  /*--- Compute properties needed for mass flow BCs. ---*/

  if (outlet) {
    SU2_OMP_SAFE_GLOBAL_ACCESS(GetOutlet_Properties(geometry, config, iMesh, Output);)
  }

  /*--- Initialize the Jacobian matrix and residual, not needed for the reducer strategy
   *    as we set blocks (including diagonal ones) and completely overwrite. ---*/

  if(!ReducerStrategy && !Output) {
    LinSysRes.SetValZero();
    if (implicit) Jacobian.SetValZero();
    else {SU2_OMP_BARRIER} // because of "nowait" in LinSysRes
  }
}

void CIncEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                    unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  const auto InnerIter = config->GetInnerIter();
  const bool muscl = config->GetMUSCL_Flow() && (iMesh == MESH_0);
  const bool center = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  const bool limiter = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE) && (InnerIter <= config->GetLimiterIter());
  const bool van_albada = (config->GetKind_SlopeLimit_Flow() == LIMITER::VAN_ALBADA_EDGE);

  /*--- Common preprocessing steps. ---*/

  CommonPreprocessing(geometry, solver_container, config, iMesh, iRKStep, RunTime_EqSystem, Output);

  /*--- Upwind second order reconstruction ---*/

  if (!Output && muscl && !center) {

    /*--- Gradient computation for MUSCL reconstruction. ---*/

    switch (config->GetKind_Gradient_Method_Recon()) {
      case GREEN_GAUSS:
        SetPrimitive_Gradient_GG(geometry, config, true); break;
      case LEAST_SQUARES:
      case WEIGHTED_LEAST_SQUARES:
        SetPrimitive_Gradient_LS(geometry, config, true); break;
      default: break;
    }

    /*--- Limiter computation ---*/

    if (limiter && !van_albada) SetPrimitive_Limiter(geometry, config);
  }
}

unsigned long CIncEulerSolver::SetPrimitive_Variables(CSolver **solver_container, const CConfig *config) {

  unsigned long iPoint, nonPhysicalPoints = 0;

  AD::StartNoSharedReading();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Incompressible flow, primitive variables ---*/

    auto physical = nodes->SetPrimVar(iPoint,GetFluidModel());

    /* Check for non-realizable states for reporting. */

    if (!physical) nonPhysicalPoints++;
  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

  return nonPhysicalPoints;
}

void CIncEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                   unsigned short iMesh, unsigned long Iteration) {

  /*--- Define an object to compute the speed of sound. ---*/
  struct SoundSpeed {
    FORCEINLINE su2double operator() (const CIncEulerVariable& nodes, unsigned long iPoint, unsigned long jPoint) const {
      return sqrt(0.5 * (nodes.GetBetaInc2(iPoint) + nodes.GetBetaInc2(jPoint)));
    }

    FORCEINLINE su2double operator() (const CIncEulerVariable& nodes, unsigned long iPoint) const {
      return sqrt(nodes.GetBetaInc2(iPoint));
    }

  } soundSpeed;

  /*--- Define an object to compute the viscous eigenvalue. ---*/
  struct LambdaVisc {
    const bool energy;

    LambdaVisc(bool e) : energy(e) {}

    FORCEINLINE su2double lambda(su2double lamVisc, su2double eddyVisc, su2double rho, su2double k, su2double cv) const {
      su2double Lambda_1 = (4.0/3.0)*(lamVisc + eddyVisc);
      su2double Lambda_2 = 0.0;
      if (energy) Lambda_2 = k / cv;
      return (Lambda_1 + Lambda_2) / rho;
    }

    FORCEINLINE su2double operator() (const CIncEulerVariable& nodes, unsigned long iPoint, unsigned long jPoint) const {
      su2double thermalCond = 0.5*(nodes.GetThermalConductivity(iPoint) + nodes.GetThermalConductivity(jPoint));
      su2double laminarVisc = 0.5*(nodes.GetLaminarViscosity(iPoint) + nodes.GetLaminarViscosity(jPoint));
      su2double eddyVisc = 0.5*(nodes.GetEddyViscosity(iPoint) + nodes.GetEddyViscosity(jPoint));
      su2double density = 0.5*(nodes.GetDensity(iPoint) + nodes.GetDensity(jPoint));
      su2double cv = 0.5*(nodes.GetSpecificHeatCv(iPoint) + nodes.GetSpecificHeatCv(jPoint));
      return lambda(laminarVisc, eddyVisc, density, thermalCond, cv);
    }

    FORCEINLINE su2double operator() (const CIncEulerVariable& nodes, unsigned long iPoint) const {
      su2double thermalCond = nodes.GetThermalConductivity(iPoint);
      su2double laminarVisc = nodes.GetLaminarViscosity(iPoint);
      su2double eddyVisc = nodes.GetEddyViscosity(iPoint);
      su2double density = nodes.GetDensity(iPoint);
      su2double cv = nodes.GetSpecificHeatCv(iPoint);
      return lambda(laminarVisc, eddyVisc, density, thermalCond, cv);
    }

  } lambdaVisc(config->GetEnergy_Equation());

  /*--- Now instantiate the generic implementation with the two functors above. ---*/

  SetTime_Step_impl(soundSpeed, lambdaVisc, geometry, solver_container, config, iMesh, Iteration);

}

void CIncEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  CNumerics* numerics = numerics_container[CONV_TERM + omp_get_thread_num()*MAX_TERMS];

  unsigned long iPoint, jPoint;

  const bool implicit    = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool jst_scheme  = ((config->GetKind_Centered_Flow() == CENTERED::JST) && (iMesh == MESH_0));
  const bool bounded_scalar = config->GetBounded_Scalar();

  /*--- For hybrid parallel AD, pause preaccumulation if there is shared reading of
  * variables, otherwise switch to the faster adjoint evaluation mode. ---*/
  bool pausePreacc = false;
  if (ReducerStrategy) pausePreacc = AD::PausePreaccumulation();
  else AD::StartNoSharedReading();

  /*--- Loop over edge colors. ---*/
  for (auto color : EdgeColoring)
  {
  /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
  SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
  for(auto k = 0ul; k < color.size; ++k) {

    auto iEdge = color.indices[k];

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

    if (bounded_scalar) EdgeMassFluxes[iEdge] = residual[0];

    /*--- Update residual value ---*/

    if (ReducerStrategy) {
      EdgeFluxes.SetBlock(iEdge, residual);
      if (implicit)
        Jacobian.SetBlocks(iEdge, residual.jacobian_i, residual.jacobian_j);
    }
    else {
      LinSysRes.AddBlock(iPoint, residual);
      LinSysRes.SubtractBlock(jPoint, residual);

      /*--- Set implicit computation ---*/
      if (implicit)
        Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
    }

    /*--- Viscous contribution. ---*/

    Viscous_Residual(iEdge, geometry, solver_container,
                     numerics_container[VISC_TERM + omp_get_thread_num()*MAX_TERMS], config);
  }
  END_SU2_OMP_FOR
  } // end color loop

  /*--- Restore preaccumulation and adjoint evaluation state. ---*/
  AD::ResumePreaccumulation(pausePreacc);
  if (!ReducerStrategy) AD::EndNoSharedReading();

  if (ReducerStrategy) {
    SumEdgeFluxes(geometry);
    if (implicit)
      Jacobian.SetDiagonalAsColumnSum();
  }

}

void CIncEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                                      CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  CNumerics* numerics = numerics_container[CONV_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Static arrays of MUSCL-reconstructed primitives and secondaries (thread safety). ---*/
  su2double Primitive_i[MAXNVAR] = {0.0}, Primitive_j[MAXNVAR] = {0.0};

  unsigned long iPoint, jPoint, counter_local = 0;
  unsigned short iDim, iVar;

  SU2_OMP_MASTER
  ErrorCounter = 0;
  END_SU2_OMP_MASTER

  const bool implicit   = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool muscl      = (config->GetMUSCL_Flow() && (iMesh == MESH_0));
  const bool limiter    = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE);
  const bool van_albada = (config->GetKind_SlopeLimit_Flow() == LIMITER::VAN_ALBADA_EDGE);
  const bool bounded_scalar = config->GetBounded_Scalar();

  /*--- For hybrid parallel AD, pause preaccumulation if there is shared reading of
  * variables, otherwise switch to the faster adjoint evaluation mode. ---*/
  bool pausePreacc = false;
  if (ReducerStrategy) pausePreacc = AD::PausePreaccumulation();
  else AD::StartNoSharedReading();

  /*--- Loop over edge colors. ---*/
  for (auto color : EdgeColoring)
  {
  /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
  SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
  for(auto k = 0ul; k < color.size; ++k) {

    auto iEdge = color.indices[k];

    /*--- Points in edge and normal vectors ---*/

    iPoint = geometry->edges->GetNode(iEdge,0); jPoint = geometry->edges->GetNode(iEdge,1);
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    /*--- Grid movement ---*/

    if (dynamic_grid)
      numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(jPoint));

    /*--- Get primitive variables ---*/

    auto V_i = nodes->GetPrimitive(iPoint);
    auto V_j = nodes->GetPrimitive(jPoint);

    /*--- High order reconstruction using MUSCL strategy ---*/

    if (muscl) {

      auto Coord_i = geometry->nodes->GetCoord(iPoint);
      auto Coord_j = geometry->nodes->GetCoord(jPoint);

      su2double Vector_ij[MAXNDIM] = {0.0};
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_ij[iDim] = 0.5*(Coord_j[iDim] - Coord_i[iDim]);
      }

      auto Gradient_i = nodes->GetGradient_Reconstruction(iPoint);
      auto Gradient_j = nodes->GetGradient_Reconstruction(jPoint);

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {

        su2double Project_Grad_i = 0.0;
        su2double Project_Grad_j = 0.0;

        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_ij[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j -= Vector_ij[iDim]*Gradient_j[iVar][iDim];
        }

        su2double lim_i = 1.0;
        su2double lim_j = 1.0;

        if (van_albada) {
          su2double V_ij = V_j[iVar] - V_i[iVar];
          lim_i = LimiterHelpers<>::vanAlbadaFunction(Project_Grad_i, V_ij, EPS);
          lim_j = LimiterHelpers<>::vanAlbadaFunction(-Project_Grad_j, V_ij, EPS);
        }
        else if (limiter) {
          lim_i = nodes->GetLimiter_Primitive(iPoint, iVar);
          lim_j = nodes->GetLimiter_Primitive(jPoint, iVar);
        }

        Primitive_i[iVar] = V_i[iVar] + lim_i * Project_Grad_i;
        Primitive_j[iVar] = V_j[iVar] + lim_j * Project_Grad_j;
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
        const bool neg_temperature_i = (Primitive_i[prim_idx.Temperature()] < 0.0);
        const bool neg_temperature_j = (Primitive_j[prim_idx.Temperature()] < 0.0);

        const bool neg_density_i  = (Primitive_i[prim_idx.Density()] < 0.0);
        const bool neg_density_j  = (Primitive_j[prim_idx.Density()] < 0.0);

        bool bad_recon = neg_temperature_i || neg_temperature_j || neg_density_i || neg_density_j;
        bad_recon = nodes->UpdateNonPhysicalEdgeCounter(iEdge, bad_recon);
        counter_local += bad_recon;

        if (bad_recon) {
          for (iVar = 0; iVar < nPrimVar; iVar++) {
            Primitive_i[iVar] = V_i[iVar];
            Primitive_j[iVar] = V_j[iVar];
          }
        }
      }

      numerics->SetPrimitive(Primitive_i, Primitive_j);

    } else {

      /*--- Set primitive variables without reconstruction ---*/

      numerics->SetPrimitive(V_i, V_j);

    }

    /*--- Compute the residual ---*/

    auto residual = numerics->ComputeResidual(config);

    if (bounded_scalar) EdgeMassFluxes[iEdge] = residual[0];

    /*--- Update residual value ---*/

    if (ReducerStrategy) {
      EdgeFluxes.SetBlock(iEdge, residual);
      if (implicit)
        Jacobian.SetBlocks(iEdge, residual.jacobian_i, residual.jacobian_j);
    }
    else {
      LinSysRes.AddBlock(iPoint, residual);
      LinSysRes.SubtractBlock(jPoint, residual);

      /*--- Set implicit computation ---*/
      if (implicit)
        Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
    }

    /*--- Viscous contribution. ---*/

    Viscous_Residual(iEdge, geometry, solver_container,
                     numerics_container[VISC_TERM + omp_get_thread_num()*MAX_TERMS], config);

  }
  END_SU2_OMP_FOR
  } // end color loop

  FinalizeResidualComputation(geometry, pausePreacc, counter_local, config);
}

void CIncEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                      CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  unsigned short iVar;
  unsigned long iPoint;

  const bool implicit       = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool rotating_frame = config->GetRotating_Frame();
  const bool axisymmetric   = config->GetAxisymmetric();
  const bool body_force     = config->GetBody_Force();
  const bool boussinesq     = (config->GetKind_DensityModel() == INC_DENSITYMODEL::BOUSSINESQ);
  const bool viscous        = config->GetViscous();
  const bool radiation      = config->AddRadiation();
  const bool vol_heat       = config->GetHeatSource();
  const bool turbulent      = (config->GetKind_Turb_Model() != TURB_MODEL::NONE);
  const bool energy         = config->GetEnergy_Equation();
  const bool streamwise_periodic             = (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE);
  const bool streamwise_periodic_temperature = config->GetStreamwise_Periodic_Temperature();

  AD::StartNoSharedReading();

  if (body_force) {

    /*--- Loop over all points ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
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
    END_SU2_OMP_FOR
  }

  if (boussinesq) {

    /*--- Loop over all points ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
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
    END_SU2_OMP_FOR
  }

  if (rotating_frame) {

    /*--- Loop over all points ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
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
    END_SU2_OMP_FOR
  }

  AD::EndNoSharedReading();

  if (axisymmetric) {

    /*--- For viscous problems, we need an additional gradient. ---*/

    if (viscous) {

      AD::StartNoSharedReading();

      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < nPoint; iPoint++) {

        su2double yCoord          = geometry->nodes->GetCoord(iPoint, 1);
        su2double yVelocity       = nodes->GetVelocity(iPoint,1);
        su2double Total_Viscosity = (nodes->GetLaminarViscosity(iPoint) +
                                     nodes->GetEddyViscosity(iPoint));
        su2double AuxVar = 0.0;
        if (yCoord > EPS)
          AuxVar = Total_Viscosity*yVelocity/yCoord;

        /*--- Set the auxiliary variable for this node. ---*/

        nodes->SetAuxVar(iPoint, 0, AuxVar);

      }
      END_SU2_OMP_FOR

      AD::EndNoSharedReading();

      /*--- Compute the auxiliary variable gradient with GG or WLS. ---*/

      if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
        SetAuxVar_Gradient_GG(geometry, config);
      }
      if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
        SetAuxVar_Gradient_LS(geometry, config);
      }

    }

    /*--- loop over points ---*/

    AD::StartNoSharedReading();

    SU2_OMP_FOR_STAT(omp_chunk_size)
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
    END_SU2_OMP_FOR

    AD::EndNoSharedReading();
  }

  if (radiation) {

    AD::StartNoSharedReading();

    CNumerics* second_numerics = numerics_container[SOURCE_SECOND_TERM + omp_get_thread_num()*MAX_TERMS];

    SU2_OMP_FOR_STAT(omp_chunk_size)
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
    END_SU2_OMP_FOR

    AD::EndNoSharedReading();
  }

  if (streamwise_periodic) {

    /*--- For turbulent streamwise periodic problems w/ energy eq, we need an additional gradient of Eddy viscosity. ---*/
    if (streamwise_periodic_temperature && turbulent) {

      AD::StartNoSharedReading();

      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        /*--- Set the auxiliary variable, Eddy viscosity mu_t, for this node. ---*/
        nodes->SetAuxVar(iPoint, 0, nodes->GetEddyViscosity(iPoint));
      }
      END_SU2_OMP_FOR

      AD::EndNoSharedReading();

      /*--- Compute the auxiliary variable gradient with GG or WLS. ---*/
      if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
        SetAuxVar_Gradient_GG(geometry, config);
      }
      if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
        SetAuxVar_Gradient_LS(geometry, config);
      }

    } // if turbulent

    if (config->GetKind_Streamwise_Periodic() == ENUM_STREAMWISE_PERIODIC::MASSFLOW) {
      /*---------------------------------------------------------------------------------------------*/
      /*--- Update the Pressure Drop [Pa] for the Momentum source term if Massflow is prescribed. ---*/
      /*--- The Pressure drop is iteratively adapted to result in the prescribed Target-Massflow. ---*/
      /*---------------------------------------------------------------------------------------------*/

      /*--- Compute update to Delta p based on massflow-difference ---*/
      const su2double Average_Density_Global = SPvals.Streamwise_Periodic_AvgDensity;
      const su2double Area_Global = SPvals.Streamwise_Periodic_BoundaryArea;
      const su2double TargetMassFlow = config->GetStreamwise_Periodic_TargetMassFlow() / (config->GetDensity_Ref() * config->GetVelocity_Ref());
      const su2double MassFlow_Global = SPvals.Streamwise_Periodic_MassFlow;
      const su2double ddP = 0.5 / ( Average_Density_Global * pow(Area_Global, 2)) * (pow(TargetMassFlow, 2) - pow(MassFlow_Global, 2));

      /*--- Store updated pressure difference ---*/
      const su2double damping_factor = config->GetInc_Outlet_Damping();
      SPvalsUpdated = SPvals;
      SPvalsUpdated.Streamwise_Periodic_PressureDrop += damping_factor*ddP;
      if (!config->GetDiscrete_Adjoint())
        SPvals = SPvalsUpdated;

      /*--- Set delta_p, m_dot, inlet_T, integrated_heat ---*/
      numerics->SetStreamwisePeriodicValues(SPvalsUpdated);
    }
    else {
      numerics->SetStreamwisePeriodicValues(SPvals);
    }

    AD::StartNoSharedReading();

    /*--- Loop over all points ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the primitive variables ---*/
      numerics->SetPrimitive(nodes->GetPrimitive(iPoint), nullptr);

      /*--- Set incompressible density ---*/
      numerics->SetDensity(nodes->GetDensity(iPoint), 0.0);

      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Load the aux variable gradient that we already computed. ---*/
      if(streamwise_periodic_temperature && turbulent)
        numerics->SetAuxVarGrad(nodes->GetAuxVarGradient(iPoint), nullptr);

      /*--- Compute the streamwise periodic source residual and add to the total ---*/
      auto residual = numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    } // for iPoint
    END_SU2_OMP_FOR

    AD::EndNoSharedReading();

    if(!streamwise_periodic_temperature && energy) {

      CNumerics* second_numerics = numerics_container[SOURCE_SECOND_TERM + omp_get_thread_num()*MAX_TERMS];

      /*--- Set delta_p, m_dot, inlet_T, integrated_heat ---*/
      second_numerics->SetStreamwisePeriodicValues(SPvals);

      /*--- This bit acts as a boundary condition rather than a source term. But logically it fits better here. ---*/
      for (auto iMarker = 0ul; iMarker < config->GetnMarker_All(); iMarker++) {

        /*--- Only "inlet"/donor periodic marker ---*/
        if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY &&
            config->GetMarker_All_PerBound(iMarker) == 1) {

          SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
          for (auto iVertex = 0ul; iVertex < nVertex[iMarker]; iVertex++) {

            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

            if (!geometry->nodes->GetDomain(iPoint)) continue;

            /*--- Load the primitive variables ---*/
            second_numerics->SetPrimitive(nodes->GetPrimitive(iPoint), nullptr);

            /*--- Set incompressible density ---*/
            second_numerics->SetDensity(nodes->GetDensity(iPoint), 0.0);

            /*--- Set the specific heat ---*/
            second_numerics->SetSpecificHeat(nodes->GetSpecificHeatCp(iPoint), 0.0);

            /*--- Set the area normal ---*/
            second_numerics->SetNormal(geometry->vertex[iMarker][iVertex]->GetNormal());

            /*--- Compute the streamwise periodic source residual and add to the total ---*/
            auto residual = second_numerics->ComputeResidual(config);
            LinSysRes.AddBlock(iPoint, residual);

          }// for iVertex
          END_SU2_OMP_FOR
        }// if periodic inlet boundary
      }// for iMarker
    }// if !streamwise_periodic_temperature
  }// if streamwise_periodic

  /*--- Check if a verification solution is to be computed. ---*/

  if (VerificationSolution) {
    if ( VerificationSolution->IsManufacturedSolution() ) {

      /*--- Get the physical time. ---*/
      su2double time = 0.0;
      if (config->GetTime_Marching() != TIME_MARCHING::STEADY) time = config->GetPhysicalTime();

      AD::StartNoSharedReading();

      /*--- Loop over points ---*/
      SU2_OMP_FOR_STAT(omp_chunk_size)
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
      END_SU2_OMP_FOR

      AD::EndNoSharedReading();
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

void CIncEulerSolver::SetMax_Eigenvalue(CGeometry *geometry, const CConfig *config) {

  /*--- Define an object to compute the speed of sound. ---*/
  struct SoundSpeed {
    FORCEINLINE su2double operator() (const CIncEulerVariable& nodes, unsigned long iPoint, unsigned long jPoint) const {
      return sqrt(0.5 * (nodes.GetBetaInc2(iPoint) + nodes.GetBetaInc2(jPoint)));
    }

    FORCEINLINE su2double operator() (const CIncEulerVariable& nodes, unsigned long iPoint) const {
      return sqrt(nodes.GetBetaInc2(iPoint));
    }

  } soundSpeed;

  /*--- Instantiate generic implementation. ---*/

  SetMax_Eigenvalue_impl(soundSpeed, geometry, config);

}

void CIncEulerSolver::SetCentered_Dissipation_Sensor(CGeometry *geometry, const CConfig *config) {

  /*--- Define an object for the sensor variable, density. ---*/
  struct SensVar {
    FORCEINLINE su2double operator() (const CIncEulerVariable& nodes, unsigned long iPoint) const {
      return nodes.GetDensity(iPoint);
    }
  } sensVar;

  /*--- Instantiate generic implementation. ---*/
  SetCentered_Dissipation_Sensor_impl(sensVar, geometry, config);

}

template<ENUM_TIME_INT IntegrationType>
FORCEINLINE void CIncEulerSolver::Explicit_Iteration(CGeometry *geometry, CSolver **solver_container,
                                                     CConfig *config, unsigned short iRKStep) {
  struct Precond {
    const CIncEulerSolver* solver;
    su2activematrix matrix;
    unsigned short nVar;

    Precond(const CIncEulerSolver* s, unsigned short n) : solver(s), nVar(n) {
      matrix.resize(nVar,nVar);
    }

    FORCEINLINE void compute(const CConfig* config, unsigned long iPoint) {
      solver->SetPreconditioner(config, iPoint, 1.0, matrix);
    }

    FORCEINLINE su2double apply(unsigned short iVar, const su2double* res, const su2double* resTrunc) const {
      su2double resPrec = 0.0;
      for (unsigned short jVar = 0; jVar < nVar; ++jVar)
        resPrec += matrix(iVar,jVar) * (res[jVar] + resTrunc[jVar]);
      return resPrec;
    }
  } precond(this, nVar);

  Explicit_Iteration_impl<IntegrationType>(precond, geometry, solver_container, config, iRKStep);
}

void CIncEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                           CConfig *config, unsigned short iRKStep) {

  Explicit_Iteration<RUNGE_KUTTA_EXPLICIT>(geometry, solver_container, config, iRKStep);
}

void CIncEulerSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                             CConfig *config, unsigned short iRKStep) {

  Explicit_Iteration<CLASSICAL_RK4_EXPLICIT>(geometry, solver_container, config, iRKStep);
}

void CIncEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  Explicit_Iteration<EULER_EXPLICIT>(geometry, solver_container, config, 0);
}

void CIncEulerSolver::PrepareImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) {

  struct IncPrec {
    const CIncEulerSolver* solver;
    const bool active = true;
    su2activematrix matrix;

    IncPrec(const CIncEulerSolver* s, unsigned short nVar) : solver(s) { matrix.resize(nVar,nVar); }

    FORCEINLINE const su2activematrix& operator() (const CConfig* config, unsigned long iPoint, su2double delta) {
      solver->SetPreconditioner(config, iPoint, delta, matrix);
      return matrix;
    }

  } precond(this, nVar);

  PrepareImplicitIteration_impl(precond, geometry, config);
}

void CIncEulerSolver::CompleteImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) {

  CompleteImplicitIteration_impl<false>(geometry, config);
}

void CIncEulerSolver::SetBeta_Parameter(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iMesh) {
  static su2double MaxVel2;
  const su2double epsilon2_default = 4.1;

  /*--- For now, only the finest mesh level stores the Beta for all levels. ---*/

  if (iMesh == MESH_0) {
    SU2_OMP_MASTER
    MaxVel2 = 0.0;
    END_SU2_OMP_MASTER
    su2double maxVel2 = 0.0;

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++)
      maxVel2 = max(maxVel2, nodes->GetVelocity2(iPoint));
    END_SU2_OMP_FOR

    SU2_OMP_CRITICAL
    MaxVel2 = max(MaxVel2, maxVel2);
    END_SU2_OMP_CRITICAL

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
    {
      maxVel2 = MaxVel2;
      SU2_MPI::Allreduce(&maxVel2, &MaxVel2, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());

      config->SetMax_Vel2(max(1e-10, MaxVel2));
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Allow an override if user supplies a large epsilon^2. ---*/

  su2double BetaInc2 = max(epsilon2_default, config->GetBeta_Factor()) * config->GetMax_Vel2();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++)
    nodes->SetBetaInc2(iPoint, BetaInc2);
  END_SU2_OMP_FOR

}

void CIncEulerSolver::SetPreconditioner(const CConfig *config, unsigned long iPoint,
                                        su2double delta, su2activematrix& Preconditioner) const {

  unsigned short iDim, jDim, iVar, jVar;

  su2double  BetaInc2, Density, dRhodT, Temperature, oneOverCp, Cp;
  su2double  Velocity[MAXNDIM] = {0.0};

  bool variable_density = (config->GetVariable_Density_Model());
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

    for (iVar = 0; iVar < nVar; iVar ++ )
      for (jVar = 0; jVar < nVar; jVar ++ )
        Preconditioner[iVar][jVar] = delta*Preconditioner[iVar][jVar];

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

  const bool implicit = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;
  const bool viscous = config->GetViscous();

  su2double Normal[MAXNDIM] = {0.0};

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- Allocate the value at the infinity ---*/

    auto V_infty = GetCharacPrimVar(val_marker, iVertex);

    /*--- Index of the closest interior node ---*/

    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

    /*--- Normal vector for this vertex (negate for outward convention) ---*/

    geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
    for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
    conv_numerics->SetNormal(Normal);

    /*--- Retrieve solution at the farfield boundary node ---*/

    auto V_domain = nodes->GetPrimitive(iPoint);

    /*--- Recompute and store the velocity in the primitive variable vector. ---*/

    for (iDim = 0; iDim < nDim; iDim++)
      V_infty[iDim+prim_idx.Velocity()] = GetVelocity_Inf(iDim);

    /*--- Far-field pressure set to static pressure (0.0). ---*/

    V_infty[prim_idx.Pressure()] = GetPressure_Inf();

    /*--- Dirichlet condition for temperature at far-field (if energy is active). ---*/

    V_infty[prim_idx.Temperature()] = GetTemperature_Inf();

    /*--- Store the density.  ---*/

    V_infty[prim_idx.Density()] = GetDensity_Inf();

    /*--- Beta coefficient stored at the node ---*/

    V_infty[prim_idx.Beta()] = nodes->GetBetaInc2(iPoint);

    /*--- Cp is needed for Temperature equation. ---*/

    V_infty[prim_idx.CpTotal()] = nodes->GetSpecificHeatCp(iPoint);

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

    if (!viscous) continue;

    /*--- Set transport properties at infinity. ---*/

    V_infty[prim_idx.LaminarViscosity()] = nodes->GetLaminarViscosity(iPoint);
    V_infty[prim_idx.EddyViscosity()] = nodes->GetEddyViscosity(iPoint);
    V_infty[prim_idx.ThermalConductivity()] = nodes->GetThermalConductivity(iPoint);

    /*--- Set the normal vector and the coordinates ---*/

    visc_numerics->SetNormal(Normal);
    su2double Coord_Reflected[MAXNDIM];
    GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                             geometry->nodes->GetCoord(iPoint), Coord_Reflected);
    visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);

    /*--- Primitive variables, and gradient ---*/

    visc_numerics->SetPrimitive(V_domain, V_infty);
    visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                      nodes->GetGradient_Primitive(iPoint));

    /*--- Turbulent kinetic energy ---*/

    if (config->GetKind_Turb_Model() == TURB_MODEL::SST)
      visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                          solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));

    /*--- Compute and update viscous residual ---*/

    auto residual_v = visc_numerics->ComputeResidual(config);
    LinSysRes.SubtractBlock(iPoint, residual_v);

    /*--- Viscous Jacobian contribution for implicit integration ---*/

    if (implicit)
      Jacobian.SubtractBlock2Diag(iPoint, residual_v.jacobian_i);

  }
  END_SU2_OMP_FOR

}

void CIncEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  unsigned long Point_Normal;
  su2double *Flow_Dir, Flow_Dir_Mag, Vel_Mag, Area, P_total, P_domain, Vn;
  su2double *V_inlet, *V_domain;
  su2double UnitFlowDir[MAXNDIM] = {0.0}, dV[MAXNDIM] = {0.0};
  su2double Damping = config->GetInc_Inlet_Damping();

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool viscous = config->GetViscous();

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  INLET_TYPE Kind_Inlet = config->GetKind_Inc_Inlet(Marker_Tag);

  su2double Normal[MAXNDIM] = {0.0};

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- Allocate the value at the inlet ---*/

    V_inlet = GetCharacPrimVar(val_marker, iVertex);

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
    Flow_Dir_Mag = GeometryToolbox::Norm(nDim, Flow_Dir);

    /*--- Store the unit flow direction vector.
     If requested, use the local boundary normal (negative),
     instead of the prescribed flow direction in the config. ---*/

    if (config->GetInc_Inlet_UseNormal()) {
      for (iDim = 0; iDim < nDim; iDim++)
        UnitFlowDir[iDim] = -Normal[iDim]/Area;
    } else {
      for (iDim = 0; iDim < nDim; iDim++)
        UnitFlowDir[iDim] = Flow_Dir[iDim]/Flow_Dir_Mag;
    }

    /*--- Retrieve solution at this boundary node. ---*/

    V_domain = nodes->GetPrimitive(iPoint);

    /*--- Neumann condition for dynamic pressure ---*/

    V_inlet[prim_idx.Pressure()] = nodes->GetPressure(iPoint);

    /*--- The velocity is either prescribed or computed from total pressure. ---*/

    switch (Kind_Inlet) {

        /*--- Velocity and temperature (if required) been specified at the inlet. ---*/

      case INLET_TYPE::VELOCITY_INLET:

        /*--- Retrieve the specified velocity and temperature for the inlet. ---*/

        Vel_Mag  = Inlet_Ptotal[val_marker][iVertex]/config->GetVelocity_Ref();

        /*--- Store the velocity in the primitive variable vector. ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          V_inlet[iDim+prim_idx.Velocity()] = Vel_Mag*UnitFlowDir[iDim];

        /*--- Dirichlet condition for temperature (if energy is active) ---*/

        V_inlet[prim_idx.Temperature()] = Inlet_Ttotal[val_marker][iVertex]/config->GetTemperature_Ref();

        break;

        /*--- Stagnation pressure has been specified at the inlet. ---*/

      case INLET_TYPE::PRESSURE_INLET:

        /*--- Retrieve the specified total pressure for the inlet. ---*/

        P_total = Inlet_Ptotal[val_marker][iVertex]/config->GetPressure_Ref();

        /*--- Store the current static pressure for clarity. ---*/

        P_domain = nodes->GetPressure(iPoint);

        /*--- Check for back flow through the inlet. ---*/

        Vn = -GeometryToolbox::DotProduct(nDim, &V_domain[prim_idx.Velocity()], Normal) / Area;

        /*--- If the local static pressure is larger than the specified
         total pressure or the velocity is directed upstream, we have a
         back flow situation. The specified total pressure should be used
         as a static pressure condition and the velocity from the domain
         is used for the BC. ---*/

        if ((P_domain > P_total) || (Vn < 0.0)) {

          /*--- Back flow: use the prescribed P_total as static pressure. ---*/

          V_inlet[prim_idx.Pressure()] = Inlet_Ptotal[val_marker][iVertex]/config->GetPressure_Ref();

          /*--- Neumann condition for velocity. ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+prim_idx.Velocity()] = V_domain[iDim+prim_idx.Velocity()];

          /*--- Neumann condition for the temperature. ---*/

          V_inlet[prim_idx.Temperature()] = nodes->GetTemperature(iPoint);

        } else {

          /*--- Update the velocity magnitude using the total pressure. ---*/

          Vel_Mag = sqrt((P_total - P_domain)/(0.5*nodes->GetDensity(iPoint)));

          /*--- Compute the delta change in velocity in each direction. ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            dV[iDim] = Vel_Mag*UnitFlowDir[iDim] - V_domain[iDim+prim_idx.Velocity()];

          /*--- Update the velocity in the primitive variable vector.
           Note we use damping here to improve stability/convergence. ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+prim_idx.Velocity()] = V_domain[iDim+prim_idx.Velocity()] + Damping*dV[iDim];

          /*--- Dirichlet condition for temperature (if energy is active) ---*/

          V_inlet[prim_idx.Temperature()] = Inlet_Ttotal[val_marker][iVertex]/config->GetTemperature_Ref();

        }

        break;

      default:
        SU2_MPI::Error("Unsupported INC_INLET_TYPE.", CURRENT_FUNCTION);
        break;
    }

    /*--- check if the inlet node is shared with a viscous wall ---*/

    if (geometry->nodes->GetViscousBoundary(iPoint)) {

      /*--- match the velocity and pressure for the viscous wall---*/

      for (iDim = 0; iDim < nDim; iDim++)
        V_inlet[iDim+prim_idx.Velocity()] = nodes->GetVelocity(iPoint,iDim);

      /*--- pressure obtained from interior ---*/

      V_inlet[prim_idx.Pressure()] = nodes->GetPressure(iPoint);
    }

    /*--- Access density at the node. This is either constant by
      construction, or will be set fixed implicitly by the temperature
      and equation of state. ---*/

    V_inlet[prim_idx.Density()] = nodes->GetDensity(iPoint);

    /*--- Beta coefficient from the config file ---*/

    V_inlet[prim_idx.Beta()] = nodes->GetBetaInc2(iPoint);

    /*--- Cp is needed for Temperature equation. ---*/

    V_inlet[prim_idx.CpTotal()] = nodes->GetSpecificHeatCp(iPoint);

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

    if (!viscous) continue;

    /*--- Set transport properties at the inlet ---*/

    V_inlet[prim_idx.LaminarViscosity()] = nodes->GetLaminarViscosity(iPoint);
    V_inlet[prim_idx.EddyViscosity()] = nodes->GetEddyViscosity(iPoint);
    V_inlet[prim_idx.ThermalConductivity()] = nodes->GetThermalConductivity(iPoint);

    /*--- Set the normal vector and the coordinates ---*/

    visc_numerics->SetNormal(Normal);
    su2double Coord_Reflected[MAXNDIM];
    GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                             geometry->nodes->GetCoord(iPoint), Coord_Reflected);
    visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);

    /*--- Primitive variables, and gradient ---*/

    visc_numerics->SetPrimitive(V_domain, V_inlet);
    visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                      nodes->GetGradient_Primitive(iPoint));

    /*--- Turbulent kinetic energy ---*/

    if (config->GetKind_Turb_Model() == TURB_MODEL::SST)
      visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                          solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));

    /*--- Compute and update residual ---*/

    auto residual_v = visc_numerics->ComputeResidual(config);

    LinSysRes.SubtractBlock(iPoint, residual_v);

    /*--- Jacobian contribution for implicit integration ---*/

    if (implicit)
      Jacobian.SubtractBlock2Diag(iPoint, residual_v.jacobian_i);
  }
  END_SU2_OMP_FOR
}

void CIncEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_outlet, *V_domain, P_Outlet = 0.0, P_domain;
  su2double mDot_Target, mDot_Old, dP, Density_Avg, Area_Outlet;
  su2double Damping = config->GetInc_Outlet_Damping();

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool viscous = config->GetViscous();
  string Marker_Tag  = config->GetMarker_All_TagBound(val_marker);

  su2double Normal[MAXNDIM] = {0.0};

  INC_OUTLET_TYPE Kind_Outlet = config->GetKind_Inc_Outlet(Marker_Tag);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- Allocate the value at the outlet ---*/

    V_outlet = GetCharacPrimVar(val_marker, iVertex);

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

      case INC_OUTLET_TYPE::PRESSURE_OUTLET:

        /*--- Retrieve the specified back pressure for this outlet. ---*/

        P_Outlet = config->GetOutlet_Pressure(Marker_Tag)/config->GetPressure_Ref();

        /*--- The pressure is prescribed at the outlet. ---*/

        V_outlet[prim_idx.Pressure()] = P_Outlet;

        /*--- Neumann condition for the velocity. ---*/

        for (iDim = 0; iDim < nDim; iDim++) {
          V_outlet[iDim+prim_idx.Velocity()] = nodes->GetVelocity(iPoint,iDim);
        }

        break;

        /*--- A mass flow target has been specified for the outlet. ---*/

      case INC_OUTLET_TYPE::MASS_FLOW_OUTLET:

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

        V_outlet[prim_idx.Pressure()] = P_Outlet;

        /*--- Neumann condition for the velocity ---*/

        for (iDim = 0; iDim < nDim; iDim++) {
          V_outlet[iDim+prim_idx.Velocity()] = nodes->GetVelocity(iPoint,iDim);
        }

        break;

    }

    /*--- Neumann condition for the temperature. ---*/

    V_outlet[prim_idx.Temperature()] = nodes->GetTemperature(iPoint);

    /*--- Access density at the interior node. This is either constant by
      construction, or will be set fixed implicitly by the temperature
      and equation of state. ---*/

    V_outlet[prim_idx.Density()] = nodes->GetDensity(iPoint);

    /*--- Beta coefficient from the config file ---*/

    V_outlet[prim_idx.Beta()] = nodes->GetBetaInc2(iPoint);

    /*--- Cp is needed for Temperature equation. ---*/

    V_outlet[prim_idx.CpTotal()] = nodes->GetSpecificHeatCp(iPoint);

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

    if (!viscous) continue;

    /*--- Set transport properties at the outlet. ---*/

    V_outlet[prim_idx.LaminarViscosity()] = nodes->GetLaminarViscosity(iPoint);
    V_outlet[prim_idx.EddyViscosity()] = nodes->GetEddyViscosity(iPoint);
    V_outlet[prim_idx.ThermalConductivity()] = nodes->GetThermalConductivity(iPoint);

    /*--- Set the normal vector and the coordinates ---*/

    visc_numerics->SetNormal(Normal);
    su2double Coord_Reflected[MAXNDIM];
    GeometryToolbox::PointPointReflect(nDim, geometry->nodes->GetCoord(Point_Normal),
                                             geometry->nodes->GetCoord(iPoint), Coord_Reflected);
    visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint), Coord_Reflected);

    /*--- Primitive variables, and gradient ---*/

    visc_numerics->SetPrimitive(V_domain, V_outlet);
    visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                      nodes->GetGradient_Primitive(iPoint));

    /*--- Turbulent kinetic energy ---*/

    if (config->GetKind_Turb_Model() == TURB_MODEL::SST)
      visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0),
                                          solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0));

    /*--- Compute and update residual ---*/

    auto residual_v = visc_numerics->ComputeResidual(config);

    LinSysRes.SubtractBlock(iPoint, residual_v);

    /*--- Jacobian contribution for implicit integration ---*/
    if (implicit)
      Jacobian.SubtractBlock2Diag(iPoint, residual_v.jacobian_i);

  }
  END_SU2_OMP_FOR

}

void CIncEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                        unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {

  /*--- Local variables ---*/

  unsigned short iVar, iMarker, iDim, iNeigh;
  unsigned long iPoint, jPoint, iEdge, iVertex;

  const su2double *V_time_nM1 = nullptr, *V_time_n = nullptr, *V_time_nP1 = nullptr;
  su2double U_time_nM1[MAXNVAR], U_time_n[MAXNVAR], U_time_nP1[MAXNVAR];
  su2double Volume_nM1, Volume_nP1, TimeStep;
  const su2double *Normal = nullptr, *GridVel_i = nullptr, *GridVel_j = nullptr;
  su2double Density, Cp;

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool first_order = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST);
  const bool second_order = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool energy = config->GetEnergy_Equation();

  const int ndim = nDim;
  auto V2U = [ndim](su2double Density, su2double Cp, const su2double* V, su2double* U) {
    U[0] = Density;
    for (int iDim = 0; iDim < ndim; iDim++) U[iDim+1] = Density*V[iDim+1];
    U[ndim+1] = Density*Cp*V[ndim+1];
  };

  /*--- Store the physical time step ---*/

  TimeStep = config->GetDelta_UnstTimeND();

  /*--- Compute the dual time-stepping source term for static meshes ---*/

  if (!dynamic_grid) {

    /*--- Loop over all nodes (excluding halos) ---*/

    AD::StartNoSharedReading();

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. These are actually
       the primitive values, but we will convert to conservatives. ---*/

      V_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      V_time_n   = nodes->GetSolution_time_n(iPoint);
      V_time_nP1 = nodes->GetSolution(iPoint);

      /*--- Access the density and Cp at this node (constant for now). ---*/

      Density = nodes->GetDensity(iPoint);
      Cp = nodes->GetSpecificHeatCp(iPoint);

      /*--- Compute the conservative variable vector for all time levels. ---*/

      V2U(Density, Cp, V_time_nM1, U_time_nM1);
      V2U(Density, Cp, V_time_n, U_time_n);
      V2U(Density, Cp, V_time_nP1, U_time_nP1);

      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/

      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/

      for (iVar = 0; iVar < nVar-!energy; iVar++) {
        if (first_order)
          LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (second_order)
          LinSysRes(iPoint,iVar) += ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                                     +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }

      /*--- Compute the Jacobian contribution due to the dual time source term. ---*/

      if (implicit) {
        su2double delta = (second_order? 1.5 : 1.0) * Volume_nP1 * Density / TimeStep;

        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian.AddVal2Diag(iPoint, iDim+1, delta);

        if (energy) delta *= Cp;
        Jacobian.AddVal2Diag(iPoint, nDim+1, delta);
      }
    }
    END_SU2_OMP_FOR

    AD::EndNoSharedReading();
  }

  else {

    /*--- For unsteady flows on dynamic meshes (rigidly transforming or
     dynamically deforming), the Geometric Conservation Law (GCL) should be
     satisfied in conjunction with the ALE formulation of the governing
     equations. The GCL prevents accuracy issues caused by grid motion, i.e.
     a uniform free-stream should be preserved through a moving grid. First,
     we will loop over the edges and boundaries to compute the GCL component
     of the dual time source term that depends on grid velocities. ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

      /*--- Compute the conservative variables. ---*/

      V_time_n = nodes->GetSolution_time_n(iPoint);
      Density = nodes->GetDensity(iPoint);
      Cp = nodes->GetSpecificHeatCp(iPoint);
      V2U(Density, Cp, V_time_n, U_time_n);

      GridVel_i = geometry->nodes->GetGridVel(iPoint);

      for (iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); iNeigh++) {

        iEdge = geometry->nodes->GetEdge(iPoint, iNeigh);
        Normal = geometry->edges->GetNormal(iEdge);

        jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
        GridVel_j = geometry->nodes->GetGridVel(jPoint);

        /*--- Determine whether to consider the normal outward or inward. ---*/
        su2double dir = (iPoint < jPoint)? 0.5 : -0.5;

        su2double Residual_GCL = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Residual_GCL += dir*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

        for (iVar = 0; iVar < nVar-!energy; iVar++)
          LinSysRes(iPoint,iVar) += U_time_n[iVar]*Residual_GCL;
      }
    }
    END_SU2_OMP_FOR

    /*--- Loop over the boundary edges ---*/

    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
          (config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY) &&
          (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {

        SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

          /*--- Get the index for node i plus the boundary face normal ---*/

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          /*--- Grid velocities stored at boundary node i ---*/

          GridVel_i = geometry->nodes->GetGridVel(iPoint);

          /*--- Compute the GCL term by dotting the grid velocity with the face
           normal. The normal is negated to match the boundary convention. ---*/

          su2double Residual_GCL = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Residual_GCL -= 0.5*(GridVel_i[iDim]+GridVel_i[iDim])*Normal[iDim];

          /*--- Compute the GCL component of the source term for node i ---*/

          V_time_n = nodes->GetSolution_time_n(iPoint);
          Density = nodes->GetDensity(iPoint);
          Cp = nodes->GetSpecificHeatCp(iPoint);
          V2U(Density, Cp, V_time_n, U_time_n);

          for (iVar = 0; iVar < nVar-!energy; iVar++)
            LinSysRes(iPoint,iVar) += U_time_n[iVar]*Residual_GCL;
        }
        END_SU2_OMP_FOR
      }
    }

    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/

    AD::StartNoSharedReading();

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. These are actually
       the primitive values, but we will convert to conservatives. ---*/

      V_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      V_time_n   = nodes->GetSolution_time_n(iPoint);
      V_time_nP1 = nodes->GetSolution(iPoint);

      /*--- Access the density and Cp at this node (constant for now). ---*/

      Density = nodes->GetDensity(iPoint);
      Cp = nodes->GetSpecificHeatCp(iPoint);

      /*--- Compute the conservative variable vector for all time levels. ---*/

      V2U(Density, Cp, V_time_nM1, U_time_nM1);
      V2U(Density, Cp, V_time_n, U_time_n);
      V2U(Density, Cp, V_time_nP1, U_time_nP1);

      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/

      Volume_nM1 = geometry->nodes->GetVolume_nM1(iPoint);
      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/

      for (iVar = 0; iVar < nVar-!energy; iVar++) {
        if (first_order)
          LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
        if (second_order)
          LinSysRes(iPoint,iVar) += (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
                                     + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      }

      /*--- Compute the Jacobian contribution due to the dual time source term. ---*/

      if (implicit) {
        su2double delta = (second_order? 1.5 : 1.0) * Volume_nP1 * Density / TimeStep;

        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian.AddVal2Diag(iPoint, iDim+1, delta);

        if (energy) delta *= Cp;
        Jacobian.AddVal2Diag(iPoint, nDim+1, delta);
      }
    }
    END_SU2_OMP_FOR

    AD::EndNoSharedReading();
  }

}

void CIncEulerSolver::GetOutlet_Properties(CGeometry *geometry, CConfig *config, unsigned short iMesh, bool Output) {

  unsigned short iDim, iMarker;
  unsigned long iVertex, iPoint;
  su2double *V_outlet = nullptr, Velocity[3], MassFlow,
  Velocity2, Density, Area;
  unsigned short iMarker_Outlet, nMarker_Outlet;
  string Inlet_TagBound, Outlet_TagBound;
  su2double Vector[MAXNDIM] = {0.0};

  bool axisymmetric = config->GetAxisymmetric();

  bool write_heads = ((((config->GetInnerIter() % (config->GetScreen_Wrt_Freq(2)*40)) == 0)
                       && (config->GetInnerIter()!= 0))
                      || (config->GetInnerIter() == 1));

  /*--- Get the number of outlet markers and check for any mass flow BCs. ---*/

  nMarker_Outlet = config->GetnMarker_Outlet();
  bool Evaluate_BC = false;
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
    Outlet_TagBound = config->GetMarker_Outlet_TagBound(iMarker_Outlet);
    if (config->GetKind_Inc_Outlet(Outlet_TagBound) == INC_OUTLET_TYPE::MASS_FLOW_OUTLET)
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

    auto *Outlet_MassFlow = new su2double[config->GetnMarker_All()];
    auto *Outlet_Density  = new su2double[config->GetnMarker_All()];
    auto *Outlet_Area     = new su2double[config->GetnMarker_All()];

    /*--- Comute MassFlow, average temp, press, etc. ---*/

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

      Outlet_MassFlow[iMarker] = 0.0;
      Outlet_Density[iMarker]  = 0.0;
      Outlet_Area[iMarker]     = 0.0;

      if ((config->GetMarker_All_KindBC(iMarker) == OUTLET_FLOW) ) {

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->nodes->GetDomain(iPoint)) {

            geometry->vertex[iMarker][iVertex]->GetNormal(Vector);

            su2double AxiFactor = 1.0;
            if (axisymmetric) {
              if (geometry->nodes->GetCoord(iPoint, 1) > EPS)
                AxiFactor = 2.0*PI_NUMBER*geometry->nodes->GetCoord(iPoint, 1);
              else {
                for (const auto jPoint : geometry->nodes->GetPoints(iPoint)) {
                  if (geometry->nodes->GetVertex(jPoint,iMarker) >= 0) {
                    AxiFactor = PI_NUMBER * geometry->nodes->GetCoord(jPoint, 1);
                  }
                }
              }
            }

            V_outlet = nodes->GetPrimitive(iPoint);

            Density = V_outlet[prim_idx.Density()];

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

    auto *Outlet_MassFlow_Local = new su2double[nMarker_Outlet];
    auto *Outlet_Density_Local  = new su2double[nMarker_Outlet];
    auto *Outlet_Area_Local     = new su2double[nMarker_Outlet];

    auto *Outlet_MassFlow_Total = new su2double[nMarker_Outlet];
    auto *Outlet_Density_Total  = new su2double[nMarker_Outlet];
    auto *Outlet_Area_Total     = new su2double[nMarker_Outlet];

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

    SU2_MPI::Allreduce(Outlet_MassFlow_Local, Outlet_MassFlow_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Outlet_Density_Local, Outlet_Density_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Outlet_Area_Local, Outlet_Area_Total, nMarker_Outlet, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

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

          cout << setprecision(5) << "Outlet Avg. Density (kg/m^3): " << config->GetOutlet_Density(Outlet_TagBound) * config->GetDensity_Ref() << endl;
          su2double Outlet_mDot = fabs(config->GetOutlet_MassFlow(Outlet_TagBound)) * config->GetDensity_Ref() * config->GetVelocity_Ref();
          su2double Outlet_mDot_Target = fabs(config->GetOutlet_Pressure(Outlet_TagBound)) / (config->GetDensity_Ref() * config->GetVelocity_Ref());
          cout << "Outlet mass flow (kg/s): " << setprecision(5) << Outlet_mDot << endl;
          cout << "target mass flow (kg/s): " << setprecision(5) << Outlet_mDot_Target << endl;
          su2double goal = 100.0*Outlet_mDot/Outlet_mDot_Target;
          cout << "Target achieved:" << setprecision(5) << goal << " % "<< endl;
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

  /*--- Adjust the number of solution variables in the restart. We always
   carry a space in nVar for the energy equation in the solver, but we only
   write it to the restart if it is active. Therefore, we must reduce nVar
   here if energy is inactive so that the restart is read correctly. ---*/

  su2double Solution[MAXNVAR] = {0.0};

  auto nVar_Restart = nVar;
  if (!(config->GetEnergy_Equation() || config->GetWeakly_Coupled_Heat())) {
    nVar_Restart--;
    Solution[nVar-1] = GetTemperature_Inf();
  }

  LoadRestart_impl(geometry, solver, config, val_iter, val_update_geo, Solution, nVar_Restart);

}

void CIncEulerSolver::SetFreeStream_Solution(const CConfig *config){
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
    nodes->SetSolution(iPoint,0, Pressure_Inf);
    for (unsigned short iDim = 0; iDim < nDim; iDim++){
      nodes->SetSolution(iPoint,iDim+1, Velocity_Inf[iDim]);
    }
    nodes->SetSolution(iPoint,nDim+1, Temperature_Inf);
  }
  END_SU2_OMP_FOR

}

unsigned long CIncEulerSolver::RegisterSolutionExtra(bool input, const CConfig* config) {
  if (config->GetKind_Streamwise_Periodic() == ENUM_STREAMWISE_PERIODIC::MASSFLOW) {
    if (input) AD::RegisterInput(SPvals.Streamwise_Periodic_PressureDrop);
    else AD::RegisterOutput(SPvalsUpdated.Streamwise_Periodic_PressureDrop);
    return 1;
  }
  return 0;
}

void CIncEulerSolver::SetAdjoint_SolutionExtra(const su2activevector& adj_sol, const CConfig* config) {
  if (config->GetKind_Streamwise_Periodic() == ENUM_STREAMWISE_PERIODIC::MASSFLOW) {
    SU2_TYPE::SetDerivative(SPvalsUpdated.Streamwise_Periodic_PressureDrop, SU2_TYPE::GetValue(adj_sol[0]));
  }
}

void CIncEulerSolver::ExtractAdjoint_SolutionExtra(su2activevector& adj_sol, const CConfig* config) {
  if (config->GetKind_Streamwise_Periodic() == ENUM_STREAMWISE_PERIODIC::MASSFLOW) {
    adj_sol[0] = SU2_TYPE::GetDerivative(SPvals.Streamwise_Periodic_PressureDrop);
  }
}
