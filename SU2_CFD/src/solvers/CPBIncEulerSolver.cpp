/*!
 * \file CPBIncEulerSolver.cpp
 * \brief Main subroutines for solving pressure based incompressible flow (Euler, Navier-Stokes, etc.).
 * \author F. Palacios, T. Economon
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/solvers/CPBIncEulerSolver.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/fluid/CConstantDensity.hpp"
#include "../../include/fluid/CIncIdealGas.hpp"
#include "../../include/fluid/CIncIdealGasPolynomial.hpp"
#include "../../include/variables/CPBIncNSVariable.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/fluid/CFluidScalar.hpp"
#include "../../include/fluid/CFluidFlamelet.hpp"
#include "../../include/fluid/CFluidModel.hpp"

// TODO: Temporarily directly copied from CIncEulerSolver's constructor to adhere to necessary
// initializations, only real difference is the definition of the nVar,nprimVar etc
CPBIncEulerSolver::CPBIncEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh,
                                 const bool navier_stokes) :
  CFVMFlowSolverBase<CPBIncEulerVariable, ENUM_REGIME::INCOMPRESSIBLE>(*geometry, *config) {
  SU2_ZONE_SCOPED


  /*--- Based on the navier_stokes boolean, determine if this constructor is
   *    being called by itself, or by its derived class CIncNSSolver. ---*/
  const string description = navier_stokes? "Navier-Stokes" : "Euler";

  unsigned long iMarker;
  ifstream restart_file;
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  int Unst_RestartIter = 0;
  bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));
  bool time_stepping = config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING;
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  const bool centered = config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED;
  const su2double* scalar_init = nullptr;

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (restart) {

    if (iMesh == MESH_0) {

      /*--- Modify file name for a dual-time unsteady restart ---*/

      if (dual_time) {
        if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
        else if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST)
          Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
        else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
      }

      /*--- Modify file name for a time stepping unsteady restart ---*/

      if (time_stepping) {
        if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
        else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      }

    }

    if (config->GetKind_Streamwise_Periodic() == ENUM_STREAMWISE_PERIODIC::MASSFLOW) {
      if (rank==MASTER_NODE) cout << "Setting streamwise periodic pressure drop from restart metadata file." << endl;
    }

    auto filename_ = config->GetFilename("flow", ".meta", Unst_RestartIter);
    Read_SU2_Restart_Metadata(geometry, config, adjoint, filename_);
  }


  /*--- Set the gamma value ---*/

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometry constants in the solver structure.
   * Incompressible flow, primitive variables (P, vx, vy, vz, T, rho, beta, lamMu, EddyMu, Kt_eff, Cp, Cv) ---*/

  nDim = geometry->GetnDim();

  /*--- Make sure to align the sizes with the constructor of CIncEulerVariable. ---*/
  nVar = nDim;
  nPrimVar = nDim + 10;
  /*--- Centered schemes only need gradients for viscous fluxes (T and v, but we need also to include P). ---*/
  nPrimVarGrad = nDim + (centered ? 2 : 4);

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

    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy, false, true);
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
  if (config->GetKind_Species_Model() != SPECIES_MODEL::NONE) scalar_init = config->GetSpecies_Init();
  GetFluidModel()->SetTDState_T(Temperature_Inf, scalar_init);
  Enthalpy_Inf = GetFluidModel()->GetEnthalpy();

  /*--- Initialize the secondary values for direct derivative approximations ---*/

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
    nodes = new CPBIncNSVariable(Density_Inf, Pressure_Inf, Velocity_Inf, Temperature_Inf, nPoint, nDim, nVar, config);
  } else {
    nodes = new CPBIncEulerVariable(Density_Inf, Pressure_Inf, Velocity_Inf, Temperature_Inf, nPoint, nDim, nVar, config);
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

  /*--- Sizing edge velocity arrays ---*/
  EdgeVelocity.resize(geometry->GetnEdge(), nDim) = su2double(0.0);
  for (auto color : EdgeColoring) {
    /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
    SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
    for(auto k = 0ul; k < color.size; ++k) {
      auto iEdge = color.indices[k];
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) 
        EdgeVelocity[iEdge][iDim] = Velocity_Inf[iDim];
    }
    END_SU2_OMP_FOR
  } // end color loop

  /*--- Sizing of correction arrays ---*/
  // MomentumCorrection.resize(nPoint, nDim) = su2double(0.0);
  // MomentumEdgeCorrection.resize(geometry->GetnEdge(), nDim) = su2double(0.0);
  // PressureCorrection.resize(nPoint) = su2double(0.0);

  /*--- Add the solver name. ---*/
  SolverName = "INC.FLOW";

  /*--- Finally, check that the static arrays will be large enough (keep this
   *    check at the bottom to make sure we consider the "final" values). ---*/
  if((nDim > MAXNDIM) || (nPrimVar > MAXNVAR))
    SU2_MPI::Error("Oops! The CPBIncEulerSolver static array sizes are not large enough.", CURRENT_FUNCTION);

}

CPBIncEulerSolver::~CPBIncEulerSolver() {
  SU2_ZONE_SCOPED

  for (auto& model : FluidModel) delete model;

}

//TODO: this funciton is largely if not entirely the same as the one from cinceulersolver
void CPBIncEulerSolver::SetNondimensionalization(CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED
  
  su2double Temperature_FreeStream = 0.0,  ModVel_FreeStream = 0.0,Energy_FreeStream = 0.0,
  ModVel_FreeStreamND = 0.0, Omega_FreeStream = 0.0, Omega_FreeStreamND = 0.0, Viscosity_FreeStream = 0.0,
  Thermal_Conductivity_FreeStream = 0.0, SpecificHeat_Cp_FreeStream = 0.0,
  Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, Pressure_Thermodynamic = 0.0, Tke_FreeStream = 0.0,
  Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Temperature_Ref = 0.0, Velocity_Ref = 0.0, Time_Ref = 0.0,
  Gas_Constant_Ref = 0.0, Omega_Ref = 0.0, Force_Ref = 0.0, Viscosity_Ref = 0.0, Conductivity_Ref = 0.0, Heat_Flux_Ref = 0.0, Energy_Ref= 0.0, Pressure_FreeStreamND = 0.0, Pressure_ThermodynamicND = 0.0, Density_FreeStreamND = 0.0,
  Temperature_FreeStreamND = 0.0, Gas_ConstantND = 0.0, Specific_Heat_CpND = 0.0, Specific_Heat_CvND = 0.0, Thermal_Expansion_CoeffND = 0.0,
  Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0}, Viscosity_FreeStreamND = 0.0,
  Thermal_Conductivity_FreeStreamND = 0.0, SpecificHeat_Cp_FreeStreamND = 0.0,
  Tke_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0,
  Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0;

  unsigned short iDim, iVar;

  /*--- Local variables ---*/

  su2double Mach     = config->GetMach();
  su2double Reynolds = config->GetReynolds();

  bool unsteady      = (config->GetTime_Marching() != TIME_MARCHING::STEADY);
  bool viscous       = config->GetViscous();
  bool turbulent     = (config->GetKind_Solver() == MAIN_SOLVER::RANS);

  bool tkeNeeded     = ((turbulent) && (config->GetKind_Turb_Model() == TURB_MODEL::SST));
  bool energy        = config->GetEnergy_Equation();
  bool boussinesq    = (config->GetKind_DensityModel() == INC_DENSITYMODEL::BOUSSINESQ);


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

  /*--- Only constant density model is currently available. ---*/

  CFluidModel* auxFluidModel = nullptr;

  switch (config->GetKind_FluidModel()) {

    case CONSTANT_DENSITY:
      auxFluidModel = new CConstantDensity(Density_FreeStream, config->GetSpecific_Heat_Cp(), Temperature_FreeStream);
      auxFluidModel->SetTDState_T(Temperature_FreeStream);
      break;
    default:
      SU2_MPI::Error("Currently, only a constant density fluid model is implemented for the pressure based incompressible solver.", CURRENT_FUNCTION);
      break;
  }

  if (viscous) {

    for (iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++)
      config->SetMu_PolyCoeffND(config->GetMu_PolyCoeff(iVar), iVar);

    /*--- Use the fluid model to compute the dimensional viscosity/conductivity. ---*/

    auxFluidModel->SetLaminarViscosityModel(config);
    auxFluidModel->SetThermalConductivityModel(config);
    Viscosity_FreeStream = auxFluidModel->GetLaminarViscosity();
    config->SetViscosity_FreeStream(Viscosity_FreeStream);
    Thermal_Conductivity_FreeStream = auxFluidModel->GetThermalConductivity();
    config->SetThermalConductivity_FreeStream(Thermal_Conductivity_FreeStream);
    SpecificHeat_Cp_FreeStream = auxFluidModel->GetCp();
    config->SetSpecificHeatCp_FreeStream(SpecificHeat_Cp_FreeStream);

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
  Gas_Constant_Ref = Velocity_Ref*Velocity_Ref/Temperature_Ref;          config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref    = Density_Ref*Velocity_Ref*Length_Ref;                config->SetViscosity_Ref(Viscosity_Ref);
  Conductivity_Ref = Viscosity_Ref*Gas_Constant_Ref;                     config->SetConductivity_Ref(Conductivity_Ref);
  //Heat_Flux_Ref    = Density_Ref*Velocity_Ref*Velocity_Ref*Velocity_Ref; config->SetHeat_Flux_Ref(Heat_Flux_Ref); //To avoid error in Friction Force routine only

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

  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref(); config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Pressure_ThermodynamicND = Pressure_Thermodynamic/config->GetPressure_Ref(); config->SetPressure_ThermodynamicND(Pressure_ThermodynamicND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();   config->SetDensity_FreeStreamND(Density_FreeStreamND);

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
  Thermal_Conductivity_FreeStreamND = Thermal_Conductivity_FreeStream / Conductivity_Ref;
  config->SetThermalConductivity_FreeStreamND(Thermal_Conductivity_FreeStreamND);
  SpecificHeat_Cp_FreeStreamND = SpecificHeat_Cp_FreeStream / Gas_Constant_Ref;
  config->SetSpecificHeatCp_FreeStreamND(SpecificHeat_Cp_FreeStreamND);

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
        fluidModel = new CConstantDensity(Density_FreeStreamND, Specific_Heat_CpND, STD_REF_TEMP / config->GetTemperature_Ref());
        break;
      
      default:
        SU2_MPI::Error("The requested fluid model is not (yet) available for the pressure based incompressible solver.", CURRENT_FUNCTION);
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

  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);

  /*--- Write output to the console if this is the master node and first domain ---*/

  if ((rank == MASTER_NODE) && (iMesh == MESH_0)) {

    cout.precision(6);

    if (config->GetRef_Inc_NonDim() == DIMENSIONAL) {
      cout << "Pressure based incompressible flow: rho_ref, vel_ref, temp_ref, p_ref" << endl;
      cout << "are set to 1.0 in order to perform a dimensional calculation." << endl;
      if (dynamic_grid) cout << "Force coefficients computed using MACH_MOTION." << endl;
      else cout << "Force coefficients computed using initial values." << endl;
    }
    else if (config->GetRef_Inc_NonDim() == INITIAL_VALUES) {
      cout << "Pressure based incompressible flow: rho_ref, vel_ref, and temp_ref" << endl;
      cout << "are based on the initial values, p_ref = rho_ref*vel_ref^2." << endl;
      if (dynamic_grid) cout << "Force coefficients computed using MACH_MOTION." << endl;
      else cout << "Force coefficients computed using initial values." << endl;
    }
    else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
      cout << "Pressure based incompressible flow: rho_ref, vel_ref, and temp_ref" << endl;
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
        cout << "No energy equation." << endl;
        break;

      default:
        SU2_MPI::Error("No other density models are availabel for the pressure based solver as of yet.", CURRENT_FUNCTION);
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
      default:
        SU2_MPI::Error("No other viscosity models are availabel for the pressure based solver as of yet.", CURRENT_FUNCTION);
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

    switch (config->GetKind_FluidModel()) {
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
      NonDimTable.PrintHeader();
      NonDimTableOut << "-- Unsteady conditions" << endl;
      NonDimTable << "Total Time" << config->GetTotal_UnstTime() << config->GetTime_Ref() << "s" << config->GetTotal_UnstTimeND();
      Unit.str("");
      NonDimTable << "Time Step" << config->GetDelta_UnstTime() << config->GetTime_Ref() << "s" << config->GetDelta_UnstTimeND();
      Unit.str("");
      NonDimTable.PrintFooter();
    }


    cout << ModelTableOut.str();
    cout << NonDimTableOut.str();
  }

}


void CPBIncEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_ZONE_SCOPED

  unsigned long ErrorCounter = 0;

  unsigned long InnerIter = config->GetInnerIter();
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool muscl            = config->GetMUSCL_Flow();
  bool limiter          = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE) && (InnerIter <= config->GetLimiterIter());
  unsigned long         iVertex;
  unsigned short        iVar, iMarker;

  /*--- Set the primitive variables ---*/

  ErrorCounter = SetPrimitive_Variables(solver_container, config);

  /*--- Compute Primitive gradients to be used in Mass flux computation and pressure residual ---*/

  if ((iMesh == MESH_0) && !Output) {

    /*--- Gradient computation for MUSCL reconstruction. ---*/

    if (config->GetKind_Gradient_Method_Recon() == GREEN_GAUSS)
      SetPrimitive_Gradient_GG(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == LEAST_SQUARES)
      SetPrimitive_Gradient_LS(geometry, config, true);
    if (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES)
      SetPrimitive_Gradient_LS(geometry, config, true);

    /*--- Limiter computation ---*/

    if ((limiter) && (iMesh == MESH_0) && !Output) {
      SetPrimitive_Limiter(geometry, config);
    }

    /*--- Compute gradient of the primitive variables for pressure source term ---*/

    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
      SetPrimitive_Gradient_GG(geometry, config, false);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      SetPrimitive_Gradient_LS(geometry, config, false);
    }
  }

  /*--- Reset flag for strong BCs. ---*/
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
    nodes->ResetStrongBC(iPoint);

  if(!ReducerStrategy && !Output) {
    LinSysRes.SetValZero();
    if (implicit) Jacobian.SetValDiagonalZero();
    else {SU2_OMP_BARRIER} // because of "nowait" in LinSysRes
  }

  /*--- Error message ---*/

  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
  }

}

void CPBIncEulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                  unsigned short iMesh) {
  SU2_ZONE_SCOPED
  
  unsigned long iPoint, ErrorCounter = 0;

  /*--- Set the current estimate of velocity as a primitive variable, needed for momentum interpolation.
   * -- This does not change the pressure, it remains the same as the old value i.e. previous (pseudo) time step. ---*/
  if (iMesh == MESH_0) ErrorCounter = SetPrimitive_Variables(solver_container, config);

  /*--- Compute gradients to be used in Rhie Chow interpolation by poisson solver ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
  }
}

//TODO: exact copy of CIncEulerSolver
unsigned long CPBIncEulerSolver::SetPrimitive_Variables(CSolver **solver_container, const CConfig *config) {

  SU2_ZONE_SCOPED

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


// TODO: exact copy of CIncEulerSolver
void CPBIncEulerSolver::SetReferenceValues(const CConfig &config) {
  SU2_ZONE_SCOPED
  
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

  if (DynamicPressureRef < EPS) {
    DynamicPressureRef = 1.0;
  }

  AeroCoeffForceRef =  DynamicPressureRef * config.GetRefArea();
  
}

void CPBIncEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                           CConfig *config, unsigned short iRKStep) {
  SU2_ZONE_SCOPED

  Explicit_Iteration<RUNGE_KUTTA_EXPLICIT>(geometry, solver_container, config, iRKStep);
}

void CPBIncEulerSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                             CConfig *config, unsigned short iRKStep) {
  SU2_ZONE_SCOPED

  Explicit_Iteration<CLASSICAL_RK4_EXPLICIT>(geometry, solver_container, config, iRKStep);
}

void CPBIncEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  SU2_ZONE_SCOPED

  Explicit_Iteration<EULER_EXPLICIT>(geometry, solver_container, config, 0);
}

void CPBIncEulerSolver::PrepareImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) {
  SU2_ZONE_SCOPED

  /*--- There is no preconditioning required for the pressure-based momentum solve, thus active = false. */
  struct IncPrec {
    const CPBIncEulerSolver* solver;
    const bool active = false;
    su2activematrix matrix;

    IncPrec(const CPBIncEulerSolver* s, unsigned short nVar) : solver(s) { matrix.resize(nVar,nVar); }

    FORCEINLINE const su2activematrix& operator() (const CConfig* config, unsigned long iPoint, su2double delta) {
      return matrix;
    }

  } precond(this, nVar);

  PrepareImplicitIteration_impl(precond, geometry, config);
}

void CPBIncEulerSolver::Centered_Residual(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics **numerics_container,
                        CConfig *config,
                        unsigned short iMesh,
                        unsigned short iRKStep) {
  SU2_ZONE_SCOPED

  CNumerics* numerics = numerics_container[CONV_TERM + omp_get_thread_num()*MAX_TERMS];

  unsigned long iPoint, jPoint;

  const bool bounded_scalar = config->GetBounded_Scalar();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

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

    /*--- Grid movement ---*/
    if (dynamic_grid)
      numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(jPoint));

    /*--- Set the edge velocity (computed / corrected by Poisson solver!) ---*/

    numerics->SetEdgeVelocity(EdgeVelocity[iEdge]);

    /*--- Compute residuals, and Jacobians ---*/

    auto conv_residual = numerics->ComputeResidual(config);

    if (bounded_scalar) {
      su2double MeanDensity = 0.5*(nodes->GetPrimitive(iPoint)[nDim+2] + nodes->GetPrimitive(jPoint)[nDim+2]);
      auto Normal = geometry->edges->GetNormal(iEdge);

      /*--- Find mass flux (note that edgevelocity itself already included grid movement) ---*/

      EdgeMassFluxes[iEdge] = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        EdgeMassFluxes[iEdge] += MeanDensity * EdgeVelocity[iEdge][iDim] * Normal[iDim];
    }

    /*--- Update residual value ---*/

    if (ReducerStrategy) {
      EdgeFluxes.SetBlock(iEdge, conv_residual);
    } else {
      LinSysRes.AddBlock(iPoint, conv_residual);
      LinSysRes.SubtractBlock(jPoint, conv_residual);
    }

    /*--- Viscous contribution, returns its Jacobians so that the matrix is updated once. ---*/

    const auto visc_residual = Viscous_Residual(
        iEdge, geometry, solver_container, numerics_container[VISC_TERM + omp_get_thread_num()*MAX_TERMS], config);

    if (implicit) UpdateJacobian(iEdge, iPoint, jPoint, conv_residual, visc_residual);
  }
  END_SU2_OMP_FOR
  } // end color loop
}

void CPBIncEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                                      CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED
  
  CNumerics* numerics = numerics_container[CONV_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Static arrays of MUSCL-reconstructed primitives and secondaries (thread safety). ---*/
  su2double Primitive_i[MAXNVAR] = {0.0}, Primitive_j[MAXNVAR] = {0.0};

  unsigned long iPoint, jPoint, counter_local = 0;

  SU2_OMP_MASTER
  ErrorCounter = 0;
  END_SU2_OMP_MASTER

  const bool implicit   = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool muscl      = (config->GetMUSCL_Flow() && (iMesh == MESH_0));
  const bool limiter    = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE);
  const bool van_albada = (config->GetKind_SlopeLimit_Flow() == LIMITER::VAN_ALBADA_EDGE);
  const bool bounded_scalar = config->GetBounded_Scalar();
  const bool multicomponent = (config->GetKind_FluidModel() == FLUID_MIXTURE);

  const su2double kappa = config->GetMUSCL_Kappa_Flow();
  const su2double musclRamp = config->GetMUSCLRampValue() * config->GetNewtonKrylovRelaxation();

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
      GeometryToolbox::Distance(nDim, Coord_j, Coord_i, Vector_ij);

      auto Gradient_i = nodes->GetGradient_Reconstruction(iPoint);
      auto Gradient_j = nodes->GetGradient_Reconstruction(jPoint);

      for (auto iVar = 0u; iVar < nPrimVarGrad; iVar++) {
        const su2double V_ij = V_j[iVar] - V_i[iVar];

        const su2double Project_Grad_i = MUSCL_Reconstruction(Gradient_i[iVar], Vector_ij, V_ij, kappa, musclRamp);
        const su2double Project_Grad_j = MUSCL_Reconstruction(Gradient_j[iVar], Vector_ij, V_ij, kappa, musclRamp);

        su2double lim_i = 1.0;
        su2double lim_j = 1.0;
        if (van_albada) {
          lim_i = LimiterHelpers<>::vanAlbadaFunction(Project_Grad_i, V_ij, EPS);
          lim_j = LimiterHelpers<>::vanAlbadaFunction(Project_Grad_j, V_ij, EPS);
        }
        else if (limiter) {
          lim_i = nodes->GetLimiter_Primitive(iPoint, iVar);
          lim_j = nodes->GetLimiter_Primitive(jPoint, iVar);
        }

        Primitive_i[iVar] = V_i[iVar] + 0.5 * lim_i * Project_Grad_i;
        Primitive_j[iVar] = V_j[iVar] - 0.5 * lim_j * Project_Grad_j;
      }

      for (auto iVar = nPrimVarGrad; iVar < nPrimVar; iVar++) {
        Primitive_i[iVar] = V_i[iVar];
        Primitive_j[iVar] = V_j[iVar];
      }

      numerics->SetPrimitive(Primitive_i, Primitive_j);

    } else {

      /*--- Set primitive variables without reconstruction ---*/

      numerics->SetPrimitive(V_i, V_j);

    }

    /*--- Set the edge velocity (computed / corrected by Poisson solver!) ---*/

    numerics->SetEdgeVelocity(EdgeVelocity[iEdge]);
    
    /*--- Compute residuals, and Jacobians ---*/

    auto conv_residual = numerics->ComputeResidual(config);

    if (bounded_scalar) {
      su2double MeanDensity = 0.5*(V_i[nDim+2] + V_j[nDim+2]);
      auto Normal = geometry->edges->GetNormal(iEdge);

      /*--- Find mass flux (note that edgevelocity itself already included grid movement) ---*/

      EdgeMassFluxes[iEdge] = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        EdgeMassFluxes[iEdge] += MeanDensity * EdgeVelocity[iEdge][iDim] * Normal[iDim];
    }

    /*--- Update residual value ---*/

    if (ReducerStrategy) {
      EdgeFluxes.SetBlock(iEdge, conv_residual);
    } else {
      LinSysRes.AddBlock(iPoint, conv_residual);
      LinSysRes.SubtractBlock(jPoint, conv_residual);
    }

    /*--- Viscous contribution, returns its Jacobians so that the matrix is updated once. ---*/

    const auto visc_residual = Viscous_Residual(
        iEdge, geometry, solver_container, numerics_container[VISC_TERM + omp_get_thread_num()*MAX_TERMS], config);

    if (implicit) UpdateJacobian(iEdge, iPoint, jPoint, conv_residual, visc_residual);

  }
  END_SU2_OMP_FOR
  } // end color loop

  FinalizeResidualComputation(geometry, pausePreacc, counter_local, config);
}

void CPBIncEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                      CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

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
  const bool multicomponent = (config->GetKind_FluidModel() == FLUID_MIXTURE);

  /*--- Add pressure source term ---*/
  unsigned long iPoint;
  unsigned short iVar;
  
  // TODO: ideas on how this could be implemented nicer would be appreciated.
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    su2double Vdpdx[MAXNDIM];

    /*--- Compute the residual based on the pressure gradient. ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Vdpdx[iVar] = geometry->nodes->GetVolume(iPoint)*nodes->GetGradient_Primitive(iPoint,0,iVar);

    auto residual = CNumerics::ResidualType<>(Vdpdx, nullptr, nullptr);

    /*--- Add Residual to the total ---*/
    LinSysRes.AddBlock(iPoint, residual);
  }
  END_SU2_OMP_FOR

  // TODO: Write all other possible source terms e.g. body force, axisymmetry etc. etc.

}

void CPBIncEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh, unsigned long Iteration) {
  SU2_ZONE_SCOPED

  /* Define an object to compute the speed of sound, as the speed of sound is theoretically infinite, 
  this makes no sense. However to be able to reuse the time step routine we artificially define the speed of sound
  to be zero such that a regular advective time step is computed */
  struct SoundSpeed {
    FORCEINLINE su2double operator() (const CPBIncEulerVariable& nodes, unsigned long iPoint, unsigned long jPoint) const {
      return 0.0;
    }

    FORCEINLINE su2double operator() (const CPBIncEulerVariable& nodes, unsigned long iPoint) const {
      return 0.0;
    }

  } soundSpeed;

  /*--- Define an object to compute the viscous eigenvalue. ---*/
  struct LambdaVisc {

    FORCEINLINE su2double lambda(su2double lamVisc, su2double eddyVisc, su2double rho) const {
      return (4.0/3.0) * (lamVisc + eddyVisc) / rho;
    }

    FORCEINLINE su2double operator() (const CPBIncEulerVariable& nodes, unsigned long iPoint, unsigned long jPoint) const {
      su2double laminarVisc = 0.5*(nodes.GetLaminarViscosity(iPoint) + nodes.GetLaminarViscosity(jPoint));
      su2double eddyVisc = 0.5*(nodes.GetEddyViscosity(iPoint) + nodes.GetEddyViscosity(jPoint));
      su2double density = 0.5*(nodes.GetDensity(iPoint) + nodes.GetDensity(jPoint));
      return lambda(laminarVisc, eddyVisc, density);
    }

    FORCEINLINE su2double operator() (const CPBIncEulerVariable& nodes, unsigned long iPoint) const {
      su2double thermalCond = nodes.GetThermalConductivity(iPoint);
      su2double laminarVisc = nodes.GetLaminarViscosity(iPoint);
      su2double eddyVisc = nodes.GetEddyViscosity(iPoint);
      su2double density = nodes.GetDensity(iPoint);
      return lambda(laminarVisc, eddyVisc, density);
    }

  } lambdaVisc;

  /*--- Now instantiate the generic implementation with the two functors above. ---*/

  SetTime_Step_impl(soundSpeed, lambdaVisc, geometry, solver_container, config, iMesh, Iteration);
}


void CPBIncEulerSolver::ComputeRhieChowVelocities(CGeometry *geometry, CSolver **solver_container) {
  SU2_ZONE_SCOPED

  unsigned short iDim;
  unsigned long iPoint, jPoint;
  su2double Edge_Vector[MAXNDIM], dist_ij_2, Grad_Avg;
  su2double *Coord_i, *Coord_j, *GridVel_i,*GridVel_j;
  su2double GradP_f[MAXNDIM], GradP_in[MAXNDIM], GradP_proj, Coeff_Mom;

  CSolver* poisson_solver = solver_container[POISSON_SOL];
  CVariable* poisson_nodes = poisson_solver->GetNodes();

  /*--- Mass flux is computed over all edges ---*/
  for (auto color : EdgeColoring) {
    SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
    for (auto k = 0ul; k < color.size; ++k) {
      auto iEdge = color.indices[k];

      iPoint = geometry->edges->GetNode(iEdge,0); jPoint = geometry->edges->GetNode(iEdge,1);
      
      if (dynamic_grid) {
        GridVel_i = geometry->nodes->GetGridVel(iPoint);
        GridVel_j = geometry->nodes->GetGridVel(jPoint);
      }

      /*--- Rhie Chow interpolation ---*/
      Coord_i = geometry->nodes->GetCoord(iPoint);
      Coord_j = geometry->nodes->GetCoord(jPoint);
      dist_ij_2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
        dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
      }
      /*--- 1. Interpolate the pressure gradient based on node values ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        Grad_Avg = 0.5*(nodes->GetGradient_Primitive(iPoint,0,iDim) + nodes->GetGradient_Primitive(jPoint,0,iDim));
        GradP_in[iDim] = Grad_Avg;
      }

      /*--- 2. Compute pressure gradient at the face ---*/
      /*--- Eq 15.62 F Moukalled, L Mangani M. Darwish OpenFOAM and uFVM book. ---*/
      GradP_proj = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        GradP_proj += GradP_in[iDim]*Edge_Vector[iDim];
      }
      if (dist_ij_2 != 0.0) {
        for (iDim = 0; iDim < nDim; iDim++) {
          GradP_f[iDim] = GradP_in[iDim] - (GradP_proj - (nodes->GetPressure(jPoint) - nodes->GetPressure(iPoint)))*Edge_Vector[iDim]/ dist_ij_2;
        }
      }

      /*--- Linearly interpolated coefficient. ---*/

      Coeff_Mom = 0.5*(poisson_nodes->GetMomCoeff(iPoint) + poisson_nodes->GetMomCoeff(jPoint));

      for (iDim = 0; iDim < nDim; iDim++) {

        /*--- Face average velocity. ---*/

        su2double meanVelocity = 0.5*(nodes->GetVelocity(iPoint,iDim) + nodes->GetVelocity(jPoint,iDim));
        
        if (dynamic_grid) {
          meanVelocity -= 0.5*(GridVel_i[iDim]+GridVel_j[iDim]);
        }

        /*--- Correction based on Rhie-Chow. ---*/

        su2double RhieChowCorrection = Coeff_Mom*(GradP_f[iDim] - GradP_in[iDim]);

        su2double CorrectedEdgeVelocity = meanVelocity - RhieChowCorrection;

        /*--- Update edge velocity. ---*/

        SetEdgeVelocity(iEdge, iDim, CorrectedEdgeVelocity);
      }

    }
    END_SU2_OMP_FOR
  }
}


void CPBIncEulerSolver::ApplyPressureVelocityCorrection(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  SU2_ZONE_SCOPED

  /*--- Start of computing the corrections ---*/
  unsigned long iEdge, iPoint, jPoint, iMarker, iVertex;
  unsigned short iDim, iVar, KindBC;
  su2double Vel, Current_Pressure, factor, PCorr_Ref, Vol, delT, Density;
  string Marker_Tag;
  su2activevector pressureCorrection, alpha_p;
  su2activematrix velocityCorrection;
  su2activematrix velocityEdgeCorrection;
  long Pref_local;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool piso = (config->GetPISO_corrections() > 1);

  // CSolver* flow_solution = solver_container[FLOW_SOL];
  // CVariable* flow_nodes = flow_solution->GetNodes();

  CSolver* poisson_solver = solver_container[POISSON_SOL];
  CVariable* poisson_nodes = poisson_solver->GetNodes();

  /*--- Allocate corrections and relaxation ---*/
  // TODO: can be moved into constructor for performance.
  pressureCorrection.resize(nPointDomain);
  velocityCorrection.resize(nPointDomain,nDim);
  velocityEdgeCorrection.resize(geometry->GetnEdge(),nDim);
  alpha_p.resize(nPointDomain);

  /*--- Combine all pressure corrections into a vector for easy access ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    pressureCorrection[iPoint] = poisson_nodes->GetSolution(iPoint,0);

  /*--- Define a reference pressure ---*/
  // TODO: look at this, currently copied (but working?) logic from old solver.
  unsigned long PRef_Point = 1;
  Pref_local = geometry->GetGlobal_to_Local_Point(PRef_Point);
  PCorr_Ref = 0.0;
  if (Pref_local >= 0)
    if(geometry->nodes->GetDomain(Pref_local))
      PCorr_Ref = 0.0;//Pressure_Correc[Pref_local];

  /*--- Compute Velocity Corrections and under relaxation factor for the pressure. ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    factor = 0.0;
    Vol = geometry->nodes->GetVolume(iPoint);
    delT = nodes->GetDelta_Time(iPoint);
    const auto view = Jacobian.GetBlockView(iPoint, iPoint);
    for (iDim = 0; iDim < nDim; iDim++) {
      velocityCorrection[iPoint][iDim] = -poisson_nodes->GetMomCoeff(iPoint)*(poisson_nodes->GetGradient(iPoint,0,iDim));
      if (implicit && !piso) factor += view(iDim, iDim);
    }

    /*--- The PISO algorithm should not underrelax pressure, for explicit iterations we simply do not have the necessary jacobian information to compute this ---*/
    if (piso) alpha_p[iPoint] = 1.0;
    else {
      alpha_p[iPoint] = config->GetRelaxation_Factor_PBFlow();
      // if (implicit) alpha_p[iPoint] *= (Vol/delT) / (factor+(Vol/delT));     
    } 
  }

  // TODO: optional addition of HbyA 
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      velocityCorrection[iPoint][iDim] += poisson_nodes->GetHbyACorrection(iPoint, iDim);
    }
  }

  /*--- Compute the edge corrections based on the average of the momentum coefficients and the average of the p' gradient. ---*/
  su2double* Coord_i,* Coord_j;
  su2double GradP_f[MAXNDIM], GradP_in[MAXNDIM];
  for (auto color : EdgeColoring) {
    SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
    for (auto k = 0ul; k < color.size; ++k) {
      auto iEdge = color.indices[k];

      iPoint = geometry->edges->GetNode(iEdge,0); jPoint = geometry->edges->GetNode(iEdge,1);

      Coord_i = geometry->nodes->GetCoord(iPoint);
      Coord_j = geometry->nodes->GetCoord(jPoint);
      su2double dist_ij_2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        su2double dist = Coord_j[iDim]-Coord_i[iDim];
        dist_ij_2 += dist*dist;
      }
      /*--- 1. Interpolate the pressure gradient based on node values ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        su2double Grad_Avg = 0.5*(poisson_nodes->GetGradient(iPoint,0,iDim) + poisson_nodes->GetGradient(jPoint,0,iDim));
        GradP_in[iDim] = Grad_Avg;
      }

      /*--- 2. Compute pressure gradient at the face ---*/
      /*--- Eq 15.62 F Moukalled, L Mangani M. Darwish OpenFOAM and uFVM book. ---*/
      su2double GradP_proj = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        GradP_proj += GradP_in[iDim]*(Coord_j[iDim]-Coord_i[iDim]);
      }
      if (dist_ij_2 != 0.0) {
        for (iDim = 0; iDim < nDim; iDim++) {
          GradP_f[iDim] = GradP_in[iDim] - (GradP_proj - (pressureCorrection[jPoint] - pressureCorrection[iPoint]))*(Coord_j[iDim]-Coord_i[iDim])/ dist_ij_2;
        }
      }

      // pressure-correction face gradient  == the Poisson stencil
      // su2double edgeVec[MAXNDIM], gradAvg[MAXNDIM];
      // GeometryToolbox::Distance(nDim, Coord_j, Coord_i, edgeVec);
      // const su2double dist2 = GeometryToolbox::SquaredNorm(nDim, edgeVec);

      // su2double proj = 0.0;
      // for (iDim = 0; iDim < nDim; iDim++) {
      //   gradAvg[iDim] = 0.5*(poisson_nodes->GetGradient(iPoint,0,iDim) +
      //                       poisson_nodes->GetGradient(jPoint,0,iDim));
      //   proj += gradAvg[iDim]*edgeVec[iDim];
      // }

      // const su2double dp   = pressureCorrection[jPoint] - pressureCorrection[iPoint];
      // const su2double d_f  = 0.5*(poisson_nodes->GetMomCoeff(iPoint) +
      //                             poisson_nodes->GetMomCoeff(jPoint));

      // for (iDim = 0; iDim < nDim; iDim++) {
      //   const su2double gradFace = gradAvg[iDim] + (dp - proj)*edgeVec[iDim]/dist2;
      //   velocityEdgeCorrection[iEdge][iDim] = -d_f*gradFace;
      //                                      // + 0.5*(HbyA_i[iDim] + HbyA_j[iDim]);
      // }
      
      for (iDim = 0; iDim < nDim; iDim++) {

        velocityEdgeCorrection[iEdge][iDim] = -0.5 * (poisson_nodes->GetMomCoeff(iPoint)+poisson_nodes->GetMomCoeff(jPoint))
                                                   * GradP_f[iDim];
        // velocityEdgeCorrection[iEdge][iDim] = -0.25 * (poisson_nodes->GetMomCoeff(iPoint)+poisson_nodes->GetMomCoeff(jPoint))
        //                                            * (poisson_nodes->GetGradient(iPoint,0,iDim)+poisson_nodes->GetGradient(jPoint,0,iDim));


        // // 2nd piso correction term.
        // velocityEdgeCorrection[iEdge][iDim] += 0.5*(poisson_nodes->GetHbyACorrection(iPoint, iDim)+poisson_nodes->GetHbyACorrection(jPoint, iDim));
      }
    }
    END_SU2_OMP_FOR
  }

  /*--- Reassign strong boundary conditions ---*/
  /*--- For now I only have velocity inlet and fully developed outlet. Will need to add other types of inlet/outlet conditions
   *  where different treatment of pressure might be needed. Symmetry and Euler wall are weak BCs. ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    KindBC = config->GetMarker_All_KindBC(iMarker);
    Marker_Tag  = config->GetMarker_All_TagBound(iMarker);
    switch (KindBC) {
      case EULER_WALL: case SYMMETRY_PLANE:
        break;

      /*--- Nothing at MPI boundaries ---*/
      case SEND_RECEIVE:
        break;

      /*--- Only a fully developed outlet is implemented. For pressure, a dirichlet
            BC has to be applied and no correction is necessary. Velocity has a neumann BC. ---*/
      case OUTLET_FLOW:{
        auto Kind_Outlet = config->GetKind_Inc_Outlet(Marker_Tag);
        switch (Kind_Outlet) {
          case INC_OUTLET_TYPE::PRESSURE_OUTLET:{
            for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              if (geometry->nodes->GetDomain(iPoint))
                pressureCorrection[iPoint] = PCorr_Ref;
            }
            break;
          }
          //TODO: other outlet types
          default: 
            SU2_MPI::Error("The requested outflow boundary condition has not yet been implemented for the pressure based poisson solver", CURRENT_FUNCTION);
            break;
        }
        break;
      }

      /*--- Only a fixed velocity inlet is implemented now. Along with the wall boundaries,
        * the velocity is known and thus no correction is necessary.---*/
      case ISOTHERMAL: case HEAT_FLUX: case INLET_FLOW: {
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (geometry->nodes->GetDomain(iPoint)) {
            for (iDim = 0; iDim < nDim; iDim++)
              velocityCorrection[iPoint][iDim] = 0.0;
            alpha_p[iPoint] = 1.0;
            }
        }
        break;
      }

      /*--- Farfield is treated as a fully developed flow for pressure and a fixed pressure is
      * used, thus no correction is necessary. The treatment for velocity depends on whether the
      * flow is into the domain or out. If flow is in, a dirichlet bc is applied and no correction
      * is made, otherwise a Neumann BC is used and velocity is adjusted. ---*/

      case FAR_FIELD:
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (geometry->nodes->GetDomain(iPoint)) {
            // Check if the boundary condition is an inlet or not
            if (nodes->GetStrongBC(iPoint)) {
              for (iDim = 0; iDim < nDim; iDim++)
                velocityCorrection[iPoint][iDim] = 0.0;
            }
            pressureCorrection[iPoint] = PCorr_Ref;
          }
        }
        
        break;

      default: 
        SU2_MPI::Error("The requested boundary condition has not yet been implemented for the pressure based poisson solver", CURRENT_FUNCTION);
        break;
    }
  }

  /*--- Apply corrections to the nodal solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Velocity corrections ---*/

    for (iVar = 0; iVar < nDim; iVar++) {
      Vel = nodes->GetVelocity(iPoint,iVar);
      Vel += velocityCorrection[iPoint][iVar];
      Density = nodes->GetDensity(iPoint);
      nodes->SetSolution(iPoint,iVar,Density*Vel);
      poisson_nodes->SetMomCorrection(iPoint,iVar,velocityCorrection[iPoint][iVar]); // TODO: TEMPORARY JUST A TEST
    }
    nodes->SetVelocity(iPoint);

    /*--- Pressure corrections ---*/

    Current_Pressure = nodes->GetPressure(iPoint);
    Current_Pressure += alpha_p[iPoint] * (pressureCorrection[iPoint] - PCorr_Ref);
    nodes->SetPrimitive(iPoint,0,Current_Pressure);
  }

  /*--- Add corrections to the edge velocities ---*/

  for (auto color : EdgeColoring) {
    SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
    for (auto k = 0ul; k < color.size; ++k) {
      auto iEdge = color.indices[k];

      for (iDim = 0; iDim < nDim; iDim++) 
        AddEdgeVelocity(iEdge, iDim, velocityEdgeCorrection[iEdge][iDim]);
    }
    END_SU2_OMP_FOR
  }

  /*--- Reset HbyA for next iteration ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      poisson_nodes->SetHbyACorrection(iPoint, iDim, 0.0);
  }


  /*--- periodic communication for both the momentum and the poisson equations as both are now updated ---*/
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
   InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
   CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);

   InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_PRESSURE);
   CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_PRESSURE);
  }

  /*--- Communicate updated velocities and pressure ---*/
  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

  InitiateComms(geometry, config, MPI_QUANTITIES::PRESSURE_VAR);
  CompleteComms(geometry, config, MPI_QUANTITIES::PRESSURE_VAR);
}

void CPBIncEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;

  const bool implicit = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;
  const bool viscous = config->GetViscous();
  bool inflow = false;

  su2double Normal[MAXNDIM] = {0.0};
  su2double Face_Flux;
  su2double *GridVel_i;

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

    /*--- Store the density.  ---*/

    V_infty[prim_idx.Density()] = GetDensity_Inf();

    /*--- Set various quantities in the numerics class ---*/

    conv_numerics->SetPrimitive(V_domain, V_infty);

    if (dynamic_grid) {
      GridVel_i = geometry->nodes->GetGridVel(iPoint);
      conv_numerics->SetGridVel(GridVel_i, GridVel_i);
    }

    /*--- Decide if the boundary should be an inlet or an outlet ---*/

    Face_Flux = 0.0;
    if (dynamic_grid)
      for (iDim = 0; iDim < nDim; iDim++) 
        Face_Flux += nodes->GetDensity(iPoint)*(V_domain[iDim+1]-GridVel_i[iDim])*Normal[iDim];
    else
      for (iDim = 0; iDim < nDim; iDim++)
        Face_Flux += nodes->GetDensity(iPoint)*V_domain[iDim+1]*Normal[iDim];

    inflow = false;
    if ((Face_Flux < 0.0) && (fabs(Face_Flux) > EPS)) inflow = true;

    if (inflow) {
      /*--- Set this face as an inlet. ---*/
      LinSysRes.SetBlock_Zero(iPoint);

      /*--- Mark as a strong bc ---*/
      nodes->SetStrongBC(iPoint);

      if (implicit)
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian.DeleteValsRowi(iPoint, iDim);
    } else {
      /*--- Set this face as an outlet. ---*/
      
      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->SetPrimitive(V_domain, V_domain);
      
      /*--- Set the edge velocity (computed / corrected by Poisson solver!) ---*/

      conv_numerics->SetEdgeVelocity(&V_domain[1]);

      auto residual = conv_numerics->ComputeResidual(config);

      LinSysRes.AddBlock(iPoint, residual);
      nodes->SetPressure(iPoint,GetPressure_Inf());

      if (implicit) 
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
    }

  }
  END_SU2_OMP_FOR

}

void CPBIncEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { 
  SU2_ZONE_SCOPED
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  unsigned long Point_Normal;
  su2double *Flow_Dir, Flow_Dir_Mag, Vel_Mag, Area, P_total, P_domain, Vn;
  su2double *V_inlet, *V_domain;
  su2double UnitFlowDir[MAXNDIM] = {0.0}, dV[MAXNDIM] = {0.0};
  su2double Damping = config->GetInc_Inlet_Damping();

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool viscous = config->GetViscous();
  const bool energy_multicomponent = config->GetKind_FluidModel() == FLUID_MIXTURE && config->GetEnergy_Equation();
  const bool species_model = config->GetKind_Species_Model() != SPECIES_MODEL::NONE;

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

    if (config->GetInletUseNormal()) {
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


      default:
        SU2_MPI::Error("Unsupported INC_INLET_TYPE for pressure based solver.", CURRENT_FUNCTION);
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

    /*--- Enforce a strong boundary condition by directly applying the velocity ---*/

    nodes->SetVelocity_Old(iPoint,V_inlet+1);

    LinSysRes.SetBlock_Zero(iPoint);

    /*--- Jacobian contribution for implicit integration ---*/

    if (implicit)
      for (iDim = 0; iDim < nDim; iDim++) 
        Jacobian.DeleteValsRowi(iPoint, iDim);

  }
  END_SU2_OMP_FOR
}

void CPBIncEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { 
  SU2_ZONE_SCOPED
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_outlet, *V_domain, P_Outlet = 0.0, P_domain;
  su2double mDot_Target, mDot_Old, dP, Density_Avg, Area_Outlet;
  su2double Damping = config->GetInc_Outlet_Damping();

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool viscous = config->GetViscous();
  const bool energy_multicomponent = config->GetKind_FluidModel() == FLUID_MIXTURE && config->GetEnergy_Equation();
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

      default:
        SU2_MPI::Error("Unsupported INC_OUTLET_TYPE for pressure based solver.", CURRENT_FUNCTION);
        break;  

    }

    /*--- Neumann condition for the temperature. ---*/

    V_outlet[prim_idx.Temperature()] = nodes->GetTemperature(iPoint);

    /*--- Access density at the interior node. This is either constant by
      construction, or will be set fixed implicitly by the temperature
      and equation of state. ---*/

    V_outlet[prim_idx.Density()] = nodes->GetDensity(iPoint);

    /*--- Cp is needed for Temperature equation. ---*/

    V_outlet[prim_idx.CpTotal()] = nodes->GetSpecificHeatCp(iPoint);

    /*-- Neumann condition for Enthalpy in energy equation. ---*/
    V_outlet[prim_idx.Enthalpy()] = nodes->GetEnthalpy(iPoint);

    /*--- Set various quantities in the solver class ---*/

    conv_numerics->SetPrimitive(V_domain, V_outlet);

    if (dynamic_grid)
      conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                geometry->nodes->GetGridVel(iPoint));

    /*--- Set the edge velocity (computed / corrected by Poisson solver!) ---*/

    conv_numerics->SetEdgeVelocity(&V_domain[1]);

    /*--- Compute the residual using an upwind scheme ---*/

    auto residual = conv_numerics->ComputeResidual(config);

    /*--- Update residual value ---*/

    LinSysRes.AddBlock(iPoint, residual);

    /*--- Jacobian contribution for implicit integration ---*/

    if (implicit) {
      Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
    }

    /*--- Viscous contribution, commented out because serious convergence problems ---*/

    if (!viscous || energy_multicomponent) continue;

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