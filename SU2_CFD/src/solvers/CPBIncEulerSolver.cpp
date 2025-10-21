/*!
 * \file CPBIncEulerSolver.cpp
 * \brief Main subroutines for solving  pressure based incompressible flow (Euler, Navier-Stokes, etc.).
 * \author F. Palacios, T. Economon, A. Koodly
 * \version 8.0.0 "Harrier"
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

#include "../../include/solvers/CPBIncEulerSolver.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/fluid/CConstantDensity.hpp"
#include "../../include/fluid/CIncIdealGas.hpp"
#include "../../include/fluid/CIncIdealGasPolynomial.hpp"
#include "../../include/variables/CPBIncNSVariable.hpp"
#include "../../include/solvers/CFVMFlowSolverBase.inl"

template class CFVMFlowSolverBase<CPBIncEulerVariable, ENUM_REGIME::INCOMPRESSIBLE>;

CPBIncEulerSolver::CPBIncEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh,
                                 const bool navier_stokes) :
  CFVMFlowSolverBase<CPBIncEulerVariable, ENUM_REGIME::INCOMPRESSIBLE>(*geometry, *config) {

  /*--- Based on the navier_stokes boolean, determine if this constructor is
   *    being called by itself, or by the Navier Stokes solver. ---*/
  const string description = navier_stokes? "Navier-Stokes" : "Euler";

  unsigned long iPoint, iVertex, iEdge, jPoint;
  unsigned short iVar, iDim, iMarker;
  ifstream restart_file;
  unsigned short nZone = geometry->GetnZone();
  bool restart   = (config->GetRestart() || config->GetRestart_Flow());
  int Unst_RestartIter;
  unsigned short iZone = config->GetiZone();
  bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));
  bool time_stepping = config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING;
  string filename_ = config->GetSolution_FileName();

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/
  if (!(!restart || (iMesh != MESH_0) || nZone > 1)) {

    /*--- Multizone problems require the number of the zone to be appended. ---*/

    if (nZone > 1) filename_ = config->GetMultizone_FileName(filename_, iZone, ".dat");

    /*--- Modify file name for a dual-time unsteady restart ---*/

    if (dual_time) {
      if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter, ".dat");
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/

    if (time_stepping) {
      Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      filename_ = config->GetUnsteady_FileName(filename_, Unst_RestartIter, ".dat");
    }
  }

  /*--- Define geometry constants in the solver structure.
   * Incompressible flow, primitive variables (P, vx, vy, vz, rho, mu_lam, mu_t) ---*/

  nDim = geometry->GetnDim();

  nVar = nDim; nPrimVar = nDim+4; nPrimVarGrad = nDim+2;

  Residual_RMS.resize(nVar,0.0); // TODO: Added PBFlow
  Residual_Max.resize(nVar,0.0);
  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nPrimVarGrad;

  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nEdge        = geometry->GetnEdge();

  /*--- Store the number of vertices on each marker for deallocation later ---*/

  nVertex.resize(nMarker);
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];

  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/

  SetNondimensionalization(config, iMesh);

  AllocateTerribleLegacyTemporaryVariables();

  /*--- Define some auxiliary vectors related to the primitive solution ---*/

  Primitive   = new su2double[nPrimVar] ();
  Primitive_i = new su2double[nPrimVar] ();
  Primitive_j = new su2double[nPrimVar] ();

  /*--- Allocate base class members. ---*/

  Allocate(*config);

  /*--- Jacobians and vector structures for implicit computations ---*/

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure ("<<description<<"). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

  }

  else {
    if (rank == MASTER_NODE) cout << "Explicit scheme. No Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
  }

  /*--- Read farfield conditions ---*/

  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();


  /*--- Initialize the solution to the far-field state everywhere. ---*/
  if (navier_stokes) {
    nodes = new CPBIncNSVariable(Density_Inf, Pressure_Inf, Velocity_Inf, nPoint, nDim, nVar, config);
  } else {
    nodes = new CPBIncEulerVariable(Density_Inf, Pressure_Inf, Velocity_Inf, nPoint, nDim, nVar, config);
  }
  SetBaseClassPointerToNodes();

  /*--- Allocate velocities at every face. This is for momentum interpolation. ---*/

  FaceVelocity.resize(nEdge,nDim);
  FaceVelocityCorrec.resize(nEdge,nDim);
  su2double Normal[MAXNDIM];
  for (iEdge = 0; iEdge < nEdge; iEdge++) {
    geometry->edges->GetNormal(iEdge, Normal);
    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);
    for (iDim = 0; iDim < nDim; iDim++) {
      FaceVelocity[iEdge][iDim]= 0.5*(nodes->GetSolution(iPoint,iDim)+nodes->GetSolution(jPoint,iDim));
      FaceVelocityCorrec[iEdge][iDim] = 0.0;
     }
   }
   PseudoTimeCorr.resize(nEdge,nDim) = su2double(0.0);
   TimeMarchingCorr_n.resize(nEdge,nDim) = su2double(0.0);
   TimeMarchingCorr_n1.resize(nEdge,nDim) = su2double(0.0);

  /*--- Initial comms. ---*/

  CommunicateInitialState(geometry, config);
  InitiateComms(geometry, config, PRESSURE_VAR);
  CompleteComms(geometry, config, PRESSURE_VAR);

  /*--- Add the solver name (max 8 characters) ---*/
  SolverName = "INC.FLOW";

}

CPBIncEulerSolver::~CPBIncEulerSolver() {

  unsigned short iMarker;
  unsigned long iVertex;
  unsigned long iEdge;

  delete [] Primitive;
  delete [] Primitive_i;
  delete [] Primitive_j;

  delete FluidModel;

}

void CPBIncEulerSolver::SetNondimensionalization(CConfig *config, unsigned short iMesh) {

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
  auto turb_model = config->GetKind_Turb_Model();

  bool unsteady      = (config->GetTime_Marching() != TIME_MARCHING::STEADY);
  bool viscous       = config->GetViscous();
  bool turbulent     = (config->GetKind_Solver() == MAIN_SOLVER::RANS);
  bool tkeNeeded     = ((turbulent) && (turb_model == TURB_MODEL::SST));

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

  /*--- Only constant density model is available. ---*/

  switch (config->GetKind_FluidModel()) {

    case CONSTANT_DENSITY:

      FluidModel = new CConstantDensity(Density_FreeStream, config->GetSpecific_Heat_Cp());
      FluidModel->SetTDState_T(Temperature_FreeStream);
      break;
    default:
      SU2_MPI::Error("Fluid model not implemented for pressure based incompressible solver.", CURRENT_FUNCTION);
      break;
  }

  if (viscous) {

    /*--- The dimensional viscosity is needed to determine the free-stream conditions.
      To accomplish this, simply set the non-dimensional coefficients to the
      dimensional ones. This will be overruled later.---*/

    // config->SetMu_RefND(config->GetMu_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref()); 
    config->SetMu_SND(config->GetMu_S());
    config->SetMu_ConstantND(config->GetMu_Constant());

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
  Gas_Constant_Ref = Velocity_Ref*Velocity_Ref/Temperature_Ref;          config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref    = Density_Ref*Velocity_Ref*Length_Ref;                config->SetViscosity_Ref(Viscosity_Ref);
  Heat_Flux_Ref    = Density_Ref*Velocity_Ref*Velocity_Ref*Velocity_Ref; config->SetHeat_Flux_Ref(Heat_Flux_Ref); //To avoid error in Friction Force routine only

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
      FluidModel->SetTDState_T(Temperature_FreeStreamND);
      break;

  }

  if (viscous) {

    /*--- Constant viscosity model ---*/

    config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);

    /*--- Set up the transport property models. ---*/

    FluidModel->SetLaminarViscosityModel(config);

  }

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
        cout << "No energy equation." << endl;
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
      }
    } else {
      ModelTable << "-" << "-";
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

void CPBIncEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

}

void CPBIncEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long ErrorCounter = 0;

  unsigned long InnerIter = config->GetInnerIter();
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool muscl            = config->GetMUSCL_Flow();
  bool limiter          = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE) && (InnerIter <= config->GetLimiterIter());
  unsigned long         iVertex;
  unsigned short        iVar, iMarker;

  /*--- Set the primitive variables ---*/

  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

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

  /*--- Reset mass flux residual. ---*/
  SetResMassFluxZero();

  /*--- Reset flag for strong BCs. ---*/
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
      nodes->ResetStrongBC(iPoint);


  /*--- Initialize the Jacobian matrices ---*/

  if (implicit) Jacobian.SetValZero();

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
bool RightSol = true;
unsigned long iPoint, ErrorCounter = 0;

  /*--- Set the current estimate of velocity as a primitive variable, needed for momentum interpolation.
   * -- This does not change the pressure, it remains the same as the old value .i.e. previous (pseudo)time step. ---*/
  if (iMesh == MESH_0) ErrorCounter = SetPrimitive_Variables(solver_container, config, true);

  /*--- Compute gradients to be used in Rhie Chow interpolation ---*/

    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
      SetPrimitive_Gradient_GG(geometry, config);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      SetPrimitive_Gradient_LS(geometry, config);
    }
}

unsigned long CPBIncEulerSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {

  unsigned long iPoint, nonPhysicalPoints = 0;
  unsigned short iVar;
  su2double pressure_val;

  bool physical = true;

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- PB Incompressible flow, primitive variables nDim+4, (P, vx, vy, vz, rho, lam_mu, eddy_visc) ---*/

    physical = nodes->SetPrimVar(iPoint, Density_Inf, config);

    /* Check for non-realizable states for reporting. */

    if (!physical) nonPhysicalPoints++;

    /*--- Initialize the convective, source and viscous residual vector ---*/

    if (!Output) LinSysRes.SetBlock_Zero(iPoint);

  }

  return nonPhysicalPoints;
}


void CPBIncEulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh, iMeshFine;
  unsigned long iPoint, index, iChildren, Point_Fine;
  auto turb_model = config->GetKind_Turb_Model();
  su2double Area_Children, Area_Parent, Coord[3] = {0.0}, *Solution_Fine;

  bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));
  bool steady_restart = config->GetSteadyRestart();
  bool time_stepping = config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING;
  bool turbulent = (config->GetKind_Solver() == MAIN_SOLVER::INC_RANS) || (config->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_INC_RANS);

  string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- This variable is only used to skip upto grid velocities.
        Pressure + ndim velocities need to be read from the restart file.
        Mass flux is also stored in restart file, need to skip it when
        reading grid velocities.---*/
  unsigned short nVar_Restart = nVar+2;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Store the number of variables for the turbulence model
   (that could appear in the restart file before the grid velocities). ---*/
  unsigned short turbVars = 0;
  if (turbulent){
    if (turb_model == TURB_MODEL::SST) turbVars = 2;
    else turbVars = 1;
  }

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
       offset in the buffer of data from the restart file and load it.
       Pressure is not part of the solution vector, so store it in the
       primitive variable vector---*/

      index = counter*Restart_Vars[1] + skipVars;
      nodes->SetPressure(iPoint_Local, Restart_Data[index]);
      for (iVar = 1; iVar <= nVar; iVar++) Solution[iVar-1] = Restart_Data[index+iVar];
      nodes->SetSolution(iPoint_Local, Solution);
      iPoint_Global_Local++;

      /*--- Remove mass flux. ---*/
      index++;

      /*--- For dynamic meshes, read in and store the
       grid coordinates and grid velocities for each node. ---*/

      if (dynamic_grid) {

        /*--- First, remove any variables for the turbulence model that
         appear in the restart file before the grid velocities. ---*/

        if (turb_model == TURB_MODEL::SA) {
          index++;
        } else if (turb_model == TURB_MODEL::SST) {
          index+=2;
        }

        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        /*--- If we are restarting the solution from a previously computed static calculation (no grid movement) ---*/
        /*--- the grid velocities are set to 0. This is useful for FSI computations ---*/

        su2double GridVel[3] = {0.0,0.0,0.0};
        if (!steady_restart) {

          /*--- Rewind the index to retrieve the Coords. ---*/
          index = counter*Restart_Vars[1];
          for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim]; }

          /*--- Move the index forward to get the grid velocities. ---*/
          index = counter*Restart_Vars[1] + skipVars + nVar_Restart + turbVars;
          for (iDim = 0; iDim < nDim; iDim++) { GridVel[iDim] = Restart_Data[index+iDim]; }
        }

        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->nodes->SetCoord(iPoint_Local, iDim, Coord[iDim]);
          geometry[MESH_0]->nodes->SetGridVel(iPoint_Local, iDim, GridVel[iDim]);
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

  if (dynamic_grid) {

    /*--- Communicate the new coordinates and grid velocities at the halos ---*/

    geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, COORDINATES);
    geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, COORDINATES);

    if (dynamic_grid) {
      geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, GRID_VELOCITY);
      geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, GRID_VELOCITY);
    }

    /*--- Recompute the edges and  dual mesh control volumes in the
     domain and on the boundaries. ---*/

    geometry[MESH_0]->SetCoord_CG();
    geometry[MESH_0]->SetControlVolume(config, UPDATE);
    geometry[MESH_0]->SetBoundControlVolume(config, UPDATE);
    geometry[MESH_0]->SetMaxLength(config);

    /*--- Update the multigrid structure after setting up the finest grid,
     including computing the grid velocities on the coarser levels. ---*/

    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      iMeshFine = iMesh-1;
      geometry[iMesh]->SetControlVolume(config, UPDATE);
      geometry[iMesh]->SetBoundControlVolume(config, UPDATE);
      geometry[iMesh]->SetCoord(geometry[iMeshFine]);
      if (dynamic_grid) {
        geometry[iMesh]->SetRestricted_GridVelocity(geometry[iMeshFine]);
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
  solver[MESH_0][FLOW_SOL]->InitiateComms(geometry[MESH_0], config, PRESSURE_VAR);
  solver[MESH_0][FLOW_SOL]->CompleteComms(geometry[MESH_0], config, PRESSURE_VAR);

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
      solver[iMesh][FLOW_SOL]->GetNodes()->SetSolution(iPoint, Solution);
    }

    solver[iMesh][FLOW_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][FLOW_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    solver[MESH_0][FLOW_SOL]->InitiateComms(geometry[iMesh], config, PRESSURE_VAR);
    solver[MESH_0][FLOW_SOL]->CompleteComms(geometry[iMesh], config, PRESSURE_VAR);
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




void CPBIncEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  CNumerics* numerics = numerics_container[CONV_TERM];
  su2double *V_i, *V_j;
  unsigned long iEdge, iPoint, jPoint, counter_local = 0, counter_global = 0;
  unsigned short iDim, iVar;
  unsigned long InnerIter = config->GetInnerIter();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Loop over all the edges ---*/

  for (iEdge = 0; iEdge < nEdge; iEdge++) {

    /*--- Points in edge and normal vectors ---*/
    iPoint = geometry->edges->GetNode(iEdge,0); jPoint = geometry->edges->GetNode(iEdge,1);
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));
    numerics->SetNeighbor(geometry->nodes->GetnNeighbor(iPoint), geometry->nodes->GetnNeighbor(jPoint));

    /*--- Grid movement ---*/
    if (dynamic_grid)
      numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(jPoint));

    /*--- Get primitive variables ---*/
    V_i = nodes->GetPrimitive(iPoint); V_j = nodes->GetPrimitive(jPoint);
    numerics->SetPrimitive(V_i, V_j);

    /*--- Compute residuals, and Jacobians ---*/
    auto residual = numerics->ComputeResidual(config);

    /*--- Update residual value ---*/
    LinSysRes.AddBlock(iPoint, residual);
    LinSysRes.SubtractBlock(jPoint, residual);


    /*--- Set implicit Jacobians ---*/
    if (implicit)
      Jacobian.UpdateBlocks(iEdge,iPoint,jPoint,residual.jacobian_i, residual.jacobian_j);
  }
}

void CPBIncEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                   CConfig *config, unsigned short iMesh) {

  CNumerics* numerics = numerics_container[CONV_TERM];
  CMatrixView<su2double>  Gradient_i, Gradient_j;
  su2double  Project_Grad_i, Project_Grad_j, Normal[3], *GV,
  *V_i, *V_j, *S_i, *S_j, *Limiter_i = NULL, *Limiter_j = NULL, Non_Physical = 1.0;
  unsigned long iEdge, iPoint, jPoint, counter_local = 0, counter_global = 0;
  unsigned short iDim, iVar;
  unsigned long InnerIter = config->GetInnerIter();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool muscl    = (config->GetMUSCL_Flow() && (iMesh == MESH_0));
  bool limiter  = ((config->GetKind_SlopeLimit_Flow() != LIMITER::NONE) && (InnerIter <= config->GetLimiterIter()));
  bool van_albada = config->GetKind_SlopeLimit_Flow() == LIMITER::VAN_ALBADA_EDGE;

  /*--- Loop over all the edges ---*/

  for (iEdge = 0; iEdge < nEdge; iEdge++) {

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
      numerics->SetPrimitive(Primitive_i, Primitive_j);
    }  else {
      /*--- Set conservative variables without reconstruction ---*/
      numerics->SetPrimitive(V_i, V_j);
      numerics->SetSecondary(S_i, S_j);

    }
    //numerics->SetFaceVel(FaceVelocity[iEdge]);

    auto residual = numerics->ComputeResidual(config);

    /*--- Update residual value ---*/
    LinSysRes.AddBlock(iPoint, residual);
    LinSysRes.SubtractBlock(jPoint, residual);


    /*--- Set implicit Jacobians ---*/
    if (implicit)
      Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);

  }

  /*--- Warning message about non-physical reconstructions. ---*/

  if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
#else
    counter_global = counter_local;
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Reconstr(counter_global);
  }

}



void CPBIncEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                   CConfig *config, unsigned short iMesh) {

  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM];

  unsigned short iVar,iDim,jVar,jDim;
  unsigned long iEdge, iPoint, jPoint;

  const bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  const bool rotating_frame = config->GetRotating_Frame();
  const bool axisymmetric   = config->GetAxisymmetric();
  const bool gravity        = (config->GetGravityForce() == YES);
  const bool body_force     = config->GetBody_Force();

  su2double **Jacobian_Temp, *Residual_Temp;

  /*--- Add pressure contribution. ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Initialize residual to zero. ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Residual[iVar] = 0.0;

    /*--- Assign the pressure gradient to the residual. ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Residual[iVar] = geometry->nodes->GetVolume(iPoint)*nodes->GetGradient_Primitive(iPoint,0,iVar);

    /*--- Add Residual ---*/
    LinSysRes.AddBlock(iPoint, Residual);
  }

  /*--- Other source terms for body force, rotation etc ---*/
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
      /*--- Compute the rotating frame source residual ---*/
      auto residual = numerics->ComputeResidual(config);
      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, residual);
    }
  }

  if (rotating_frame) {

    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the primitive variables ---*/
      numerics->SetPrimitive(nodes->GetPrimitive(iPoint), NULL);

      /*--- Set incompressible density ---*/
      numerics->SetDensity(nodes->GetDensity(iPoint), 0.0);

      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Compute the rotating frame source residual (CSourceIncRotatingFrame_Flow class) ---*/
      auto residual = numerics->ComputeResidual(config);

      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }

  if (axisymmetric) {

  }
}

void CPBIncEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res,alfa,Mom_Coeff[3];
  unsigned short iVar, jVar;
  unsigned long iPoint;

  for (iVar = 0; iVar < nVar; iVar++) { //TODO: TO be fixed for PBFlow was replaced in preprocess
    // SetResidual_RMS(iVar, 0.0); // TODO: PBFlow
    Residual_RMS[iVar] = 0.0; // TODO: PBFlow: Updated with current std in the code.
    Residual_Max[iVar] = 0.0;
  }
  alfa = 1.0;
  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = (geometry->nodes->GetVolume(iPoint) +
           geometry->nodes->GetPeriodicVolume(iPoint));
    Delta = nodes->GetDelta_Time(iPoint) / Vol;

    local_Res_TruncError = nodes->GetResTruncError(iPoint);
    local_Residual = LinSysRes.GetBlock(iPoint);

    /*--- Workaround to deal with nodes that are part of multiple boundaries and where
     *    one face might be strong BC and another weak BC (mostly for farfield boundary
     *    where the boundary face is strong or weak depending on local flux. ---*/
    if (nodes->GetStrongBC(iPoint)) {
       LinSysRes.SetBlock_Zero(iPoint);
    }

    for (iVar = 0; iVar < nVar; iVar ++ ) { // TODO: PBFlow
      Res = local_Residual[iVar] + local_Res_TruncError[iVar];
      nodes->AddSolution(iPoint, iVar, -alfa*Res*Delta);
      // AddRes_RMS(iVar, Res*Res);  // Old COde
      Residual_RMS[iVar] += Res*Res; // TODO: PBFlow: Updated with current std in the code.
      AddRes_Max(iVar, fabs(Res), geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint));
      
        //   AddRes_Max(iVar, fabs(Res), geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint));
    }
  }
  SetIterLinSolver(1);

  /*-- Note here that there is an assumption that solution[0] is pressure/density and velocities start from 1 ---*/
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}

void CPBIncEulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iVar, iDim;
  unsigned long iPoint, total_index, IterLinSol = 0;
  su2double Delta, *local_Res_TruncError, Vol,Mom_Coeff[3];

  /*--- Set maximum residual to zero ---*/

//   for (iVar = 0; iVar < nVar; iVar++) { // TODO: PBFlow
    // SetRes_RMS(iVar, 0.0);
    SetResToZero(); // TODO: PBFlow: Updated with current std in the code.
//     Residual_Max[iVar] = 0.0;
//   }

//   config->SetLinear_Solver_Iter(config->GetPoisson_Linear_Solver_Iter()); // TODO:PBFlow

  /*--- Build implicit system ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Workaround to deal with nodes that are part of multiple boundaries and where
     *    one face might be strong BC and another weak BC (mostly for farfield boundary
     *    where the boundary face is strong or weak depending on local flux. ---*/
    if (nodes->GetStrongBC(iPoint)) {
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
      LinSysRes.SetBlock_Zero(iPoint);
    }

    /*--- Read the volume ---*/
    Vol = (geometry->nodes->GetVolume(iPoint) +
           geometry->nodes->GetPeriodicVolume(iPoint));

    if (nodes->GetDelta_Time(iPoint) != 0.0) {
      Delta = Vol / nodes->GetDelta_Time(iPoint);
      Jacobian.AddVal2Diag(iPoint, Delta);
    } else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }

    /*--- Read the residual ---*/

    local_Res_TruncError = nodes->GetResTruncError(iPoint);
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      Residual_RMS[iVar] += LinSysRes[total_index]*LinSysRes[total_index]; // TODO:PBFlow
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
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      nodes->AddSolution(iPoint, iVar, config->GetRelaxation_Factor_PBFlow()*LinSysSol[iPoint*nVar+iVar]);
    }
  }

  /*-- Note here that there is an assumption that solution[0] is pressure/density and velocities start from 1 ---*/
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  /*--- MPI solution ---*/
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);

}

void CPBIncEulerSolver::SetMomCoeff(CGeometry *geometry, CSolver **solver_container, CConfig *config, bool periodic, unsigned short iMesh) {

  unsigned short iVar, jVar, iDim, jDim;
  unsigned long iPoint, jPoint, iNeigh;
  su2double Mom_Coeff[MAXNDIM], Mom_Coeff_nb[MAXNDIM], Vol, delT;
  bool simplec = (config->GetKind_PBIter() == ENUM_PBITER::SIMPLEC);
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  if (!periodic)  {

    if (implicit) {
      /* First sum up the momentum coefficient using the jacobian from given point and it's neighbors. ---*/
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

        /*--- Self contribution. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          Mom_Coeff[iVar] = Jacobian.GetBlock(iPoint,iPoint,iVar,iVar);
        }

        /*--- Contribution from neighbors. ---*/
        nodes->Set_Mom_Coeff_nbZero(iPoint);

        if (simplec) {
          for (iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); iNeigh++) {
            jPoint = geometry->nodes->GetPoint(iPoint,iNeigh);
            for (iVar = 0; iVar < nVar; iVar++) {
              nodes->Add_Mom_Coeff_nb(iPoint, Jacobian.GetBlock(iPoint,jPoint,iVar,iVar),iVar);
            }
          }
        }

        Vol = geometry->nodes->GetVolume(iPoint); delT = nodes->GetDelta_Time(iPoint);

        for (iVar = 0; iVar < nVar; iVar++) {
          Mom_Coeff[iVar] = Mom_Coeff[iVar] - nodes->Get_Mom_Coeff_nb(iPoint, iVar) - config->GetRCFactor()*(Vol/delT);
          //Mom_Coeff[iVar] = nodes->GetDensity(iPoint)*Vol/Mom_Coeff[iVar];
          Mom_Coeff[iVar] = Vol/Mom_Coeff[iVar];
        }

        nodes->Set_Mom_Coeff(iPoint, Mom_Coeff);
      }
    }
    else {
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

        Vol = geometry->nodes->GetVolume(iPoint); delT = nodes->GetDelta_Time(iPoint);

        for(iVar = 0; iVar < nVar; iVar++)
          Mom_Coeff[iVar] = delT/Vol;

        nodes->Set_Mom_Coeff(iPoint, Mom_Coeff);
      }
    }
  /*--- Insert MPI call here. ---*/
  InitiateComms(geometry, config, MOM_COEFF);
  CompleteComms(geometry, config, MOM_COEFF);
  }

}


void CPBIncEulerSolver::SetMomCoeffPer(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iVar, jVar, iDim, jDim;
  unsigned long iPoint, jPoint, iNeigh;
  su2double Mom_Coeff[MAXNDIM], Mom_Coeff_nb[MAXNDIM], Vol, delT;
  bool simplec = (config->GetKind_PBIter() == ENUM_PBITER::SIMPLEC);

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Self contribution. ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      Mom_Coeff[iVar] = Jacobian.GetBlock(iPoint,iPoint,iVar,iVar);
    }

    /*--- Contribution from neighbors. ---*/
    nodes->Set_Mom_Coeff_nbZero(iPoint);
    if (simplec) {
      for (iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); iNeigh++) {
        jPoint = geometry->nodes->GetPoint(iPoint,iNeigh);
        for (iVar = 0; iVar < nVar; iVar++) {
          nodes->Add_Mom_Coeff_nb(iPoint, Jacobian.GetBlock(iPoint,jPoint,iVar,iVar),iVar);
        }
      }
    }

    Vol = geometry->nodes->GetVolume(iPoint);   delT = nodes->GetDelta_Time(iPoint);

    for (iVar = 0; iVar < nVar; iVar++) {
      Mom_Coeff[iVar] = Mom_Coeff[iVar] + (1.0 - config->GetRCFactor())*(Vol/delT) - nodes->Get_Mom_Coeff_nb(iPoint, iVar);
      //Mom_Coeff[iVar] = nodes->GetDensity(iPoint)*Vol/Mom_Coeff[iVar];
      Mom_Coeff[iVar] = Vol/Mom_Coeff[iVar];
    }
    nodes->Set_Mom_Coeff(iPoint, Mom_Coeff);
  }

  /*--- Insert MPI call here. ---*/
  InitiateComms(geometry, config, MOM_COEFF); // TODO: PBFlow: FIXED
  CompleteComms(geometry, config, MOM_COEFF); // TODO: PBFlow: FIXED

}

void CPBIncEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh, unsigned long Iteration) {

  su2double *Normal, Area, Vol, Length, Mean_ProjVel = 0.0, Mean_Vel, Mean_Density,
  Mean_BetaInc2 = 4.1, Lambda, Local_Delta_Time, Mean_DensityInc, GradVel, Lambda1,
  Global_Delta_Time = 1E6, Global_Delta_UnstTimeND, ProjVel, RefProjFlux, MinRefProjFlux, ProjVel_i, ProjVel_j;

  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker, iVar;

  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool time_steping  = config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING;
  bool dual_time     = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));

  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0; MinRefProjFlux = 0.0;
  Normal = new su2double[nDim];

  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    nodes->SetMax_Lambda_Inv(iPoint, 0.0);

  /*--- Loop interior edges ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Point identification, Normal vector and area ---*/

    iPoint = geometry->edges->GetNode(iEdge,0);
    jPoint = geometry->edges->GetNode(iEdge,1);

    geometry->edges->GetNormal(iEdge, Normal);


    Area = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
    Area = sqrt(Area);

    /*--- Compute flux across the cell ---*/

    Mean_Density = 0.5*(nodes->GetDensity(iPoint) + nodes->GetDensity(jPoint));
    Mean_ProjVel = 0.0;
    for (iVar = 0; iVar < nVar; iVar++) {
        Mean_Vel = 0.5*(nodes->GetVelocity(iPoint, iVar) + nodes->GetVelocity(jPoint, iVar));
        Mean_ProjVel += (Mean_Vel*Normal[iVar]);
    }

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

    RefProjFlux = fabs(config->GetInc_Velocity_Ref()*Area);
    MinRefProjFlux = max(RefProjFlux, MinRefProjFlux);

    Lambda = fabs(Mean_ProjVel) + RefProjFlux;

    /*--- Inviscid contribution ---*/

    if (geometry->nodes->GetDomain(iPoint)) nodes->AddMax_Lambda_Inv(iPoint, Lambda+EPS);
    if (geometry->nodes->GetDomain(jPoint)) nodes->AddMax_Lambda_Inv(jPoint, Lambda+EPS);

  }

  /*--- Loop boundary edges ---*/

  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

      /*--- Point identification, Normal vector and area ---*/

      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      geometry->vertex[iMarker][iVertex]->GetNormal(Normal);

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt(Area);

      /*--- Compute flux across the cell ---*/
      Mean_Density = nodes->GetDensity(iPoint);
      Mean_ProjVel = 0.0;
      for (iVar = 0; iVar < nVar; iVar++) {
          Mean_Vel = nodes->GetVelocity(iPoint, iVar);
          Mean_ProjVel += Mean_Vel*Normal[iVar];
      }

      /*--- Adjustment for grid movement ---*/

      if (dynamic_grid) {
        su2double *GridVel = geometry->nodes->GetGridVel(iPoint);
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }

      RefProjFlux = fabs(config->GetInc_Velocity_Ref()*Area);
      MinRefProjFlux = max(RefProjFlux, MinRefProjFlux);

      Lambda = fabs(Mean_ProjVel) + RefProjFlux;

      if (geometry->nodes->GetDomain(iPoint)) {
        nodes->AddMax_Lambda_Inv(iPoint, Lambda+EPS);
      }
     }
    }
  }

  /*--- Local time-stepping: each element uses their own speed for steady state
   simulations or for pseudo time steps in a dual time simulation. ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    Vol = geometry->nodes->GetVolume(iPoint);

    if (Vol != 0.0) {
      Local_Delta_Time  = config->GetCFL(iMesh)*Vol/nodes->GetMax_Lambda_Inv(iPoint);
      if (!implicit) Local_Delta_Time  = 0.5*config->GetCFL(iMesh)*Vol/nodes->GetMax_Lambda_Inv(iPoint);
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time    = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time    = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      nodes->SetDelta_Time(iPoint, Local_Delta_Time);
     }
     else {
       nodes->SetDelta_Time(iPoint, 0.0);
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

  if (time_steping) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++){

      /*--- Sets the regular CFL equal to the unsteady CFL ---*/

      config->SetCFL(iMesh,config->GetUnst_CFL());

      /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step, otherwise
       it computes the time step based on the unsteady CFL ---*/

      if (config->GetCFL(iMesh) == 0.0){
        nodes->SetDelta_Time(iPoint, config->GetDelta_UnstTime());
      } else {
        nodes->SetDelta_Time(iPoint, Global_Delta_Time);
      }
    }
  }

  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/

  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);

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
        nodes->SetDelta_Time(iPoint, Local_Delta_Time);
      }
    }

}

void CPBIncEulerSolver::SetPoissonSourceTerm(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {


  unsigned short iVar, jVar, iDim, jDim, KindBC;
  unsigned long iPoint, jPoint, iEdge, iMarker, iVertex, iNeigh, n_inlet;
  su2double Edge_Vector[MAXNDIM], dist_ij_2;
  su2double *Coord_i, *Coord_j;
  su2double MassFlux_Part, MassFlux_Avg, Mom_Coeff[MAXNDIM], *Normal,Vel_Avg, Grad_Avg;
  su2double Area, MeanDensity, Vol , TimeStep;
  su2double GradP_f[MAXNDIM], GradP_in[MAXNDIM], GradP_proj, RhieChowInterp, Coeff_Mom, PsCorr[MAXNDIM], PsCorrFace;
  su2double *Flow_Dir, Flow_Dir_Mag, Vel_Mag, Adj_Mass,*GridVel_i,*GridVel_j;
  su2double Net_Mass, alfa, Mass_In, Mass_Out, Mass_Free_In, Mass_Free_Out, Mass_Corr, Area_out;
  string Marker_Tag;
  su2double ProjGridVelFlux, *MeshVel_i, *MeshVel_j;
  Normal = new su2double [nDim];
  bool unsteady = (config->GetTime_Marching() != TIME_MARCHING::STEADY);

  /*--- Initialize mass flux to zero ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    nodes->SetMassFluxZero(iPoint);

  Net_Mass = 0.0;

  for (iEdge = 0; iEdge < nEdge; iEdge++) {

    iPoint = geometry->edges->GetNode(iEdge,0); jPoint = geometry->edges->GetNode(iEdge,1);
    geometry->edges->GetNormal(iEdge, Normal);

    Area = GeometryToolbox::Norm(nDim, Normal);

    MeanDensity = 0.5*(nodes->GetDensity(iPoint) + nodes->GetDensity(jPoint));
    
    if (dynamic_grid) {
      GridVel_i = geometry->nodes->GetGridVel(iPoint);
      GridVel_j = geometry->nodes->GetGridVel(jPoint);
    }

    /*--- Face average mass flux. ---*/
    MassFlux_Avg = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Vel_Avg = 0.5*(nodes->GetVelocity(iPoint,iDim) + nodes->GetVelocity(jPoint,iDim));
      MassFlux_Avg += MeanDensity*Vel_Avg*Normal[iDim];
    }
    
    if (dynamic_grid)
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Avg = 0.5*(GridVel_i[iDim]+GridVel_j[iDim]);
        MassFlux_Avg -= MeanDensity*Vel_Avg*Normal[iDim];
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

    /*--- 2. Compute pressure gradient at the face ---*
         Eq 15.62 F Moukalled, L Mangani M. Darwish OpenFOAM and uFVM book. ---*/
    GradP_proj = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      GradP_proj += GradP_in[iDim]*Edge_Vector[iDim];
    }
    if (dist_ij_2 != 0.0) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradP_f[iDim] = GradP_in[iDim] - (GradP_proj - (nodes->GetPressure(jPoint) - nodes->GetPressure(iPoint)))*Edge_Vector[iDim]/ dist_ij_2;
      }
    }

    /*--- Correct the massflux by adding the pressure terms.
     * --- GradP_f is the gradient computed directly at the face and GradP_in is the
     * --- gradient linearly interpolated based on node values. This effectively adds a third
     * --- order derivative of pressure to remove odd-even decoupling of pressure and velocities.
     * --- GradP_f = (p_F^n - p_P^n)/ds , GradP_in = 0.5*(GradP_P^n + GradP_F^n)---*/
    RhieChowInterp = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      /*--- Linearly interpolated coefficient. ---*/
      Coeff_Mom = 0.5*(nodes->Get_Mom_Coeff(iPoint,iDim) + nodes->Get_Mom_Coeff(jPoint,iDim));
      /*--- Difference of pressure gradients. ---*/
      RhieChowInterp += Coeff_Mom*(GradP_f[iDim] - GradP_in[iDim])*Normal[iDim]*MeanDensity;
      /*--- Save the pressure gradient contribution for the correction term used in the next iteration. ---*/
      PsCorr[iDim] = -Coeff_Mom*(GradP_f[iDim] - GradP_in[iDim]);
    }

    /*--- Rhie Chow correction for time step must go here ---*/
    su2double beta = 0.0, beta_n = 0.0, beta_n1 = 0.0;
    su2double den_i,den_j,num_i,num_n_i,num_n1_i,num_j,num_n1_j,num_n_j;
    if (unsteady) {
      TimeStep = config->GetDelta_UnstTimeND();
      
      Vol = 0.5*(geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetVolume(jPoint));
      
      for (iDim = 0; iDim < nDim; iDim++) {
          den_i = geometry->nodes->GetVolume(iPoint)/nodes->Get_Mom_Coeff(iPoint,iDim);
          den_j = geometry->nodes->GetVolume(jPoint)/nodes->Get_Mom_Coeff(jPoint,iDim);
          Coeff_Mom = 0.5*(nodes->Get_Mom_Coeff(iPoint,iDim) + nodes->Get_Mom_Coeff(jPoint,iDim));
      }
      PsCorrFace = 0.0;
      if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) {
        for (iDim = 0; iDim < nDim; iDim++) {
          // Coefficient of correction for pseudo time iteration
          num_i = den_i - geometry->nodes->GetVolume(iPoint)/TimeStep;
          num_j = den_j - geometry->nodes->GetVolume(jPoint)/TimeStep;
          beta = 0.5*(num_i/den_i + num_j/den_j);
          // Coefficient of correction for unsteady time level n
          num_n_i = -geometry->nodes->GetVolume(iPoint)/TimeStep;
          num_n_j = -geometry->nodes->GetVolume(jPoint)/TimeStep;
          beta_n = -0.5*(num_n_i/den_i + num_n_j/den_j);

          PsCorrFace += (beta*PseudoTimeCorr[iEdge][iDim] + beta_n*TimeMarchingCorr_n[iEdge][iDim])*Normal[iDim]*MeanDensity;
        }
      }
      if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND) {
        for (iDim = 0; iDim < nDim; iDim++) {
          // Coefficient of correction for pseudo time iteration
          num_i = den_i - 3.0*geometry->nodes->GetVolume(iPoint)/(2.0*TimeStep);
          num_j = den_j - 3.0*geometry->nodes->GetVolume(jPoint)/(2.0*TimeStep);
          beta = 0.5*(num_i/den_i + num_j/den_j);
          // Coefficient of correction for unsteady time level n
          num_n_i = -4.0*geometry->nodes->GetVolume(iPoint)/(2.0*TimeStep);
          num_n_j = -4.0*geometry->nodes->GetVolume(jPoint)/(2.0*TimeStep);
          beta_n = -0.5*(num_n_i/den_i + num_n_j/den_j);
          // Coefficient of correction for unsteady time level n-1
          num_n1_i = geometry->nodes->GetVolume(iPoint)/(2.0*TimeStep);
          num_n1_j = geometry->nodes->GetVolume(jPoint)/(2.0*TimeStep);
          beta_n1 = 0.5*(num_n1_i/den_i + num_n1_j/den_j);

          PsCorrFace += (beta*PseudoTimeCorr[iEdge][iDim] + beta_n*TimeMarchingCorr_n[iEdge][iDim] + beta_n1*TimeMarchingCorr_n1[iEdge][iDim])*Normal[iDim]*MeanDensity;
        }
      }
    }
    else {
      beta = 1.0; PsCorrFace = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        PsCorrFace += PseudoTimeCorr[iEdge][iDim]*Normal[iDim]*MeanDensity;
      PsCorrFace = beta*PsCorrFace;
    }

    /*--- Calculate the mass flux at the face including the linearly interpolated velocities, pressure 
     *    gradient difference contribution and correction for time stepping (both pseudo and dual time). ---*/ 
    MassFlux_Part = MassFlux_Avg - RhieChowInterp + PsCorrFace;


    /*--- Update correction of face velocity for the next iteration. ---*/
    for (iDim = 0; iDim < nDim; iDim++) 
      PseudoTimeCorr[iEdge][iDim] = PsCorr[iDim] ;//+ PseudoTimeCorr[iEdge][iDim];

    if (geometry->nodes->GetDomain(iPoint)) nodes->AddMassFlux(iPoint,MassFlux_Part);
    if (geometry->nodes->GetDomain(jPoint)) nodes->SubtractMassFlux(jPoint,MassFlux_Part);
  }

  /*--- Mass flux correction for outflow ---*/
  Mass_In = 0.0; Mass_Out = 0.0; Mass_Free_In = 0.0; Mass_Free_Out = 0.0;
  Area_out = 0.0; Adj_Mass = 0.0; n_inlet = 0;
  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    KindBC = config->GetMarker_All_KindBC(iMarker);
    Marker_Tag  = config->GetMarker_All_TagBound(iMarker);

    switch (KindBC) {
    /*--- Wall boundaries have zero mass flux (irrespective of grid movement) ---*/
      case EULER_WALL: case ISOTHERMAL: case HEAT_FLUX: case SYMMETRY_PLANE:
      break;

      case INLET_FLOW:
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->nodes->GetDomain(iPoint)) {

            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);

            if (dynamic_grid)
              GridVel_i = geometry->nodes->GetGridVel(iPoint);
                
            MassFlux_Part = 0.0;
            if (dynamic_grid)
              for (iDim = 0; iDim < nDim; iDim++)
                MassFlux_Part -= nodes->GetDensity(iPoint)*(nodes->GetVelocity(iPoint, iDim)-GridVel_i[iDim])*Normal[iDim];
            else
             for (iDim = 0; iDim < nDim; iDim++)
              MassFlux_Part -= nodes->GetDensity(iPoint)*(nodes->GetVelocity(iPoint, iDim))*Normal[iDim];

            if (geometry->nodes->GetDomain(iPoint))
              nodes->AddMassFlux(iPoint, MassFlux_Part);

             /*--- Sum up the mass flux entering to be used for mass flow correction at outflow ---*/
             Mass_In += fabs(MassFlux_Part);
          }
        }
        break;

      case FAR_FIELD:
      /*--- Treat the farfield as a fully developed outlet for pressure. I still compute the mass fluxes
       * to use when dealing with truncated open boundaries (not implemented yet). ---*/
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->nodes->GetDomain(iPoint)) {
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);

            if (dynamic_grid)
              GridVel_i = geometry->nodes->GetGridVel(iPoint);
                
            MassFlux_Part = 0.0;
            if (dynamic_grid)
              for (iDim = 0; iDim < nDim; iDim++)
                MassFlux_Part -= nodes->GetDensity(iPoint)*(nodes->GetVelocity(iPoint, iDim)-GridVel_i[iDim])*Normal[iDim];
            else
             for (iDim = 0; iDim < nDim; iDim++)
              MassFlux_Part -= nodes->GetDensity(iPoint)*(nodes->GetVelocity(iPoint, iDim))*Normal[iDim];

            if ((MassFlux_Part < 0.0) && (fabs(MassFlux_Part) > EPS)) {
              Mass_Free_In += fabs(MassFlux_Part);
              nodes->AddMassFlux(iPoint, MassFlux_Part);
            }
            else {
              Mass_Free_Out += fabs(MassFlux_Part);
              nodes->AddMassFlux(iPoint, MassFlux_Part);
            }

            nodes->SetMassFluxZero(iPoint);
          }
        }
        break;

      case OUTLET_FLOW:{
      /*--- Note I am assuming a fully developed outlet, thus the pressure value is prescribed
       * -- and a dirichlet bc has to be applied along outlet faces. The Massflux, which forms the RHS
       * -- of the equation, is set to zero to enforce the dirichlet bc. ---*/

        auto Kind_Outlet = config->GetKind_Inc_Outlet(Marker_Tag);

        switch (Kind_Outlet) {
          case INC_OUTLET_TYPE::PRESSURE_OUTLET:
            for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

              if (geometry->nodes->GetDomain(iPoint)) {
                geometry->vertex[iMarker][iVertex]->GetNormal(Normal);

                if (dynamic_grid)
                  GridVel_i = geometry->nodes->GetGridVel(iPoint);
                
                MassFlux_Part = 0.0;
                if (dynamic_grid)
                  for (iDim = 0; iDim < nDim; iDim++)
                    MassFlux_Part -= nodes->GetDensity(iPoint)*(nodes->GetVelocity(iPoint, iDim)-GridVel_i[iDim])*Normal[iDim];
                else
                  for (iDim = 0; iDim < nDim; iDim++)
                    MassFlux_Part -= nodes->GetDensity(iPoint)*(nodes->GetVelocity(iPoint, iDim))*Normal[iDim];

                /*--- Sum up the mass flux leaving to be used for mass flow correction at outflow ---*/
                Mass_Out += fabs(MassFlux_Part);

                nodes->AddMassFlux(iPoint, MassFlux_Part);
              }
            }
            break;
          // Not working properly - Also removed the mass flux correction outside the loop which was summing up mass in and out
          case INC_OUTLET_TYPE::OPEN: {
            for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

              if (geometry->nodes->GetDomain(iPoint)) {
                geometry->vertex[iMarker][iVertex]->GetNormal(Normal);

                MassFlux_Part = 0.0;
                for (iDim = 0; iDim < nDim; iDim++)
                  MassFlux_Part -= nodes->GetDensity(iPoint)*nodes->GetVelocity(iPoint,iDim)*Normal[iDim];

                /*--- Sum up the mass flux leaving to be used for mass flow correction at outflow ---*/
                if (MassFlux_Part > 0.0)
                  Mass_Free_Out += fabs(MassFlux_Part);
                else
                  Mass_Free_In += fabs(MassFlux_Part);
                }
              }
            }
            break;
          }
        }
        break;

        default:{
            cout<<"Invalid option!"<<endl;
            break;
        }
    }
  }

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    AddResMassFlux(nodes->GetMassFlux(iPoint)*nodes->GetMassFlux(iPoint));
    //cout<<iPoint<<"\t"<<nodes->GetMassFlux(iPoint)<<endl;
  }

  InitiateComms(geometry, config, MASS_FLUX);
  CompleteComms(geometry, config, MASS_FLUX);

  SetResMassFluxRMS(geometry, config);
  delete [] Normal;
}

void CPBIncEulerSolver::SetMomentumCorrection_DualTime() {

unsigned long  iEdge;
unsigned short iDim;

for (iEdge = 0; iEdge < nEdge; iEdge++) {
  for (iDim = 0; iDim < nDim; iDim++) {
    TimeMarchingCorr_n1[iEdge][iDim] = TimeMarchingCorr_n[iEdge][iDim];
    TimeMarchingCorr_n[iEdge][iDim] = PseudoTimeCorr[iEdge][iDim];
  }
}

}


void CPBIncEulerSolver::SetResMassFluxRMS(CGeometry *geometry, CConfig *config) {
  unsigned short iVar;

#ifndef HAVE_MPI

  if (GetResMassFlux() != GetResMassFlux()) {
      SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);
  }

  ResMassFlux = sqrt(ResMassFlux/geometry->GetnPoint());

#else

  int nProcessor = size, iProcessor;

  su2double sbuf_residual, rbuf_residual;
  unsigned long  Global_nPointDomain;
  unsigned short iDim;

  /*--- Set the L2 Norm residual in all the processors ---*/

  sbuf_residual  = 0.0;
  rbuf_residual  = 0.0;

  sbuf_residual = GetResMassFlux();

  if (config->GetComm_Level() == COMM_FULL) {

    unsigned long Local_nPointDomain = geometry->GetnPointDomain();
    SU2_MPI::Allreduce(&sbuf_residual, &rbuf_residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  } else {

    /*--- Reduced MPI comms have been requested. Use a local residual only. ---*/

    rbuf_residual = sbuf_residual;
    Global_nPointDomain = geometry->GetnPointDomain();

  }


  if (rbuf_residual != rbuf_residual) {
    SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);
  }

  SetResMassFlux(max(EPS*EPS, sqrt(rbuf_residual/Global_nPointDomain)));
#endif

}


void CPBIncEulerSolver:: Flow_Correction(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned long iEdge, iPoint, jPoint, iMarker, iVertex;
  unsigned short iDim, iVar, KindBC;
  su2double Vel, Current_Pressure, factor, PCorr_Ref, Vol, delT, Density;
  string Marker_Tag;
  su2activevector Pressure_Correc, alpha_p;
  su2activematrix vel_corr;
  long Pref_local;
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Allocate corrections and relaxation ---*/
  Pressure_Correc.resize(nPointDomain);
  vel_corr.resize(nPointDomain,nVar);
  alpha_p.resize(nPointDomain);

  /*--- Pressure Corrections ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    Pressure_Correc[iPoint] = solver_container[POISSON_SOL]->GetNodes()->GetSolution(iPoint,0);

  Pref_local = geometry->GetGlobal_to_Local_Point(PRef_Point);
  PCorr_Ref = 0.0;
  if (Pref_local >= 0)
    if(geometry->nodes->GetDomain(Pref_local))
    PCorr_Ref = 0.0;//Pressure_Correc[Pref_local];

  /*--- Velocity Corrections. ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    factor = 0.0;
    Vol = geometry->nodes->GetVolume(iPoint);
    delT = nodes->GetDelta_Time(iPoint);
    for (iVar = 0; iVar < nVar; iVar++) {
      vel_corr[iPoint][iVar] = nodes->Get_Mom_Coeff(iPoint,iVar)*(solver_container[POISSON_SOL]->GetNodes()->GetGradient(iPoint,0,iVar));
      if (implicit) factor += Jacobian.GetBlock(iPoint, iPoint, iVar, iVar);
    }
    if (implicit)
      alpha_p[iPoint] = config->GetRelaxation_Factor_PBFlow()*(Vol/delT) / (factor+(Vol/delT));
    else
      alpha_p[iPoint] = config->GetRelaxation_Factor_PBFlow();
  }

  /*--- Reassign strong boundary conditions ---*/
  /*--- For now I only have velocity inlet and fully developed outlet. Will need to add other types of inlet/outlet conditions
   *  where different treatment of pressure might be needed. Symmetry and Euler wall are weak BCs. ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    KindBC = config->GetMarker_All_KindBC(iMarker);
    Marker_Tag  = config->GetMarker_All_TagBound(iMarker);
    switch (KindBC) {
      /*--- Only a fully developed outlet is implemented. For pressure, a dirichlet
            BC has to be applied and no correction is necessary. Velocity has a neumann BC. ---*/
      case OUTLET_FLOW:{
        auto Kind_Outlet = config->GetKind_Inc_Outlet(Marker_Tag);
        switch (Kind_Outlet) {
          case INC_OUTLET_TYPE::PRESSURE_OUTLET:{
            for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              if (geometry->nodes->GetDomain(iPoint))
                Pressure_Correc[iPoint] = PCorr_Ref;
            }
          break;
          }
          // Not working yet
          case INC_OUTLET_TYPE::OPEN:{
            break;
          }
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
              vel_corr[iPoint][iDim] = 0.0;
            alpha_p[iPoint] = 1.0;
           }
        }
        break;
        }

        /*--- Farfield is treated as a fully developed flow for pressure and a fixed pressure is
         * used, thus no correction is necessary. The treatment for velocity depends on whether the
         * flow is into the domain or out. If flow is in, a dirichlet bc is applied and no correction
         * is made, otherwise a Neumann BC is used and velocity is adjusted. The fixed value of velocity
         * is the one from the previous iteration. ---*/

         case FAR_FIELD:{
           for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
             iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
             if (geometry->nodes->GetDomain(iPoint)) {
               // Check if inlet or not
               if (nodes->GetStrongBC(iPoint)) {
                 for (iDim = 0; iDim < nDim; iDim++)
                   vel_corr[iPoint][iDim] = 0.0;
                }
               Pressure_Correc[iPoint] = PCorr_Ref;
             }
            }
        break;
        }
        default:{
        break;
        }
    }
  }

  /*--- Velocity corrections ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      Vel = nodes->GetVelocity(iPoint,iVar);
      Vel = Vel - vel_corr[iPoint][iVar];
      Density = nodes->GetDensity(iPoint);
      nodes->SetSolution(iPoint,iVar,Density*Vel);
    }
    nodes->SetVelocity(iPoint);
    /*--- Pressure corrections ---*/
    Current_Pressure = nodes->GetPressure(iPoint);
    Current_Pressure += alpha_p[iPoint]*(Pressure_Correc[iPoint] - PCorr_Ref);
    nodes->SetPressure(iPoint,Current_Pressure);
  }
  
  /*--- Correct face velocity. ---*/
  /*su2double Area, Vel_Mag,rho, GradP_in[MAXNDIM],GradP_f[MAXNDIM], GradP_proj;
  su2double *Coord_i, *Coord_j, dist_ij, delP, Pressure_j, Pressure_i;
  su2double Edge_Vector[MAXNDIM], dist_ij_2, proj_vector_ij;
  for (iEdge = 0; iEdge < nEdge; iEdge++) {

    iPoint = geometry->edges->GetNode(iEdge,0); jPoint = geometry->edges->GetNode(iEdge,1);
    geometry->edges->GetNormal(iEdge, Normal);
    dist_ij_2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Edge_Vector[iDim] = (geometry->nodes->GetCoord(jPoint, iDim) -
                           geometry->nodes->GetCoord(iPoint, iDim));
      dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    }

    for (iDim = 0; iDim < nDim; iDim++)
      GradP_in[iDim] = ( solver_container[POISSON_SOL]->GetNodes()->GetGradient(iPoint,0,iDim) +
                         solver_container[POISSON_SOL]->GetNodes()->GetGradient(jPoint,0,iDim) );

    GradP_proj = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      GradP_proj += GradP_in[iDim]*Edge_Vector[iDim];
    }
    if (dist_ij_2 != 0.0) {
      for (iDim = 0; iDim < nDim; iDim++) {
        GradP_f[iDim] = GradP_in[iDim] - (GradP_proj - (Pressure_Correc[jPoint] - Pressure_Correc[iPoint]))*Edge_Vector[iDim]/ dist_ij_2;
      }
    }

    for (iDim =0; iDim < nDim; iDim++) {
      Coeff_Mom = 0.5*(nodes->Get_Mom_Coeff(iPoint,iDim) + nodes->Get_Mom_Coeff(jPoint,iDim));
      FaceVelocity[iEdge][iDim] = FaceVelocity[iEdge][iDim] - Coeff_Mom*(GradP_f[iDim]);
    }
  }*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);

    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_PRESSURE);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_PRESSURE);
  }

  /*--- Communicate updated velocities and pressure ---*/
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  InitiateComms(geometry, config, PRESSURE_VAR);
  CompleteComms(geometry, config, PRESSURE_VAR);

  /*--- Reset pressure corrections to zero for next iteration. ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    solver_container[POISSON_SOL]->GetNodes()->SetSolution(iPoint,0,0.0);
  }

  /*--- Communicate updated Poisson solution ---*/
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    solver_container[POISSON_SOL]->InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    solver_container[POISSON_SOL]->CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  solver_container[POISSON_SOL]->InitiateComms(geometry, config, SOLUTION);
  solver_container[POISSON_SOL]->CompleteComms(geometry, config, SOLUTION);

}


void CPBIncEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                        unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {

  /*--- Local variables ---*/

  unsigned short iVar, jVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;

  su2double *U_time_nM1, *U_time_n, *U_time_nP1;
  su2double Volume_nM1, Volume_nP1, TimeStep;
  su2double *GridVel_i = nullptr, *GridVel_j = nullptr, Residual_GCL;
  const su2double* Normal;

  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Store the physical time step ---*/

  TimeStep = config->GetDelta_UnstTimeND();

  /*--- Compute the dual time-stepping source term for static meshes ---*/

  if (!dynamic_grid) {

    /*--- Loop over all nodes (excluding halos) ---*/

    for (iPoint = 0; iPoint < nPoint; iPoint++) {

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

      U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      U_time_n   = nodes->GetSolution_time_n(iPoint);
      U_time_nP1 = nodes->GetSolution(iPoint);

      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/

      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order). ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)
          Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                            +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }

      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/

      LinSysRes.AddBlock(iPoint, Residual);

      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }
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

      U_time_n = nodes->GetSolution_time_n(iPoint);

      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;

      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Compute the GCL component of the source term for node j ---*/

      U_time_n = nodes->GetSolution_time_n(jPoint);

      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;

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

        U_time_n = nodes->GetSolution_time_n(iPoint);

        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = U_time_n[iVar]*Residual_GCL;

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

      U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      U_time_n   = nodes->GetSolution_time_n(iPoint);
      U_time_nP1 = nodes->GetSolution(iPoint);

      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/

      Volume_nM1 = geometry->nodes->GetVolume_nM1(iPoint);
      Volume_nP1 = geometry->nodes->GetVolume(iPoint);

      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
        if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
          + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      }

      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }
        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
      }
    }
  }
}

void CPBIncEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim, jDim, iVar;
  unsigned long iVertex, iPoint, Point_Normal, total_index;

  su2double *V_infty, *V_domain;
  su2double Face_Flux, Flux0, Flux1, MeanDensity, proj_vel;
  su2double *Coord_i, *Coord_j, dist_ij, delP, Pressure_j, Pressure_i;
  bool implicit       = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  bool viscous        = config->GetViscous();
  bool inflow         = false;

  su2double *Normal = new su2double[nDim];
  auto turb_model = config->GetKind_Turb_Model();

  string Marker_Tag  = config->GetMarker_All_TagBound(val_marker);
  su2double *GridVel_i;

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Set farfield soultion. ---*/
      V_infty = GetCharacPrimVar(val_marker, iVertex);
      V_infty[0] = GetPressure_Inf();
      for (iDim = 0; iDim < nDim; iDim++)
        V_infty[iDim+1] = GetVelocity_Inf(iDim);
      V_infty[nDim+1] = nodes->GetDensity(iPoint);

      if (dynamic_grid) {
        GridVel_i = geometry->nodes->GetGridVel(iPoint);   
        conv_numerics->SetGridVel(GridVel_i,GridVel_i);
      }

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

        nodes->SetStrongBC(iPoint);

        if (implicit) {
          for (iDim = 0; iDim < nDim; iDim++) {
            total_index = iPoint*nVar+iDim;
            Jacobian.DeleteValsRowi(total_index);
          }
        }
      }
      else {

        if (dynamic_grid)
          conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));
        /*--- Compute the residual using an upwind scheme ---*/

        conv_numerics->SetPrimitive(V_domain, V_domain);
        auto residual = conv_numerics->ComputeResidual(config);

        LinSysRes.AddBlock(iPoint, residual);
        nodes->SetPressure(iPoint,GetPressure_Inf());

        if (implicit) 
         Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
      }

      /*--- Set transport properties at the outlet (used by turb solver also). ---*/
      if (viscous) {
        V_domain[nDim+2] = nodes->GetLaminarViscosity(iPoint);
        V_domain[nDim+3] = nodes->GetEddyViscosity(iPoint);
        V_infty[nDim+2] = nodes->GetLaminarViscosity(iPoint);
        V_infty[nDim+3] = nodes->GetEddyViscosity(iPoint);
      }
      if (viscous && !inflow) {

        /*--- Set the normal vector and the coordinates ---*/
        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                                geometry->nodes->GetCoord(Point_Normal));

        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_domain);
        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                          nodes->GetGradient_Primitive(iPoint));

        /*--- Turbulent kinetic energy ---*/
        if (turb_model == TURB_MODEL::SST)
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


void CPBIncEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim, iVar, jDim;
  unsigned long iVertex, iPoint,total_index, Point_Normal;
  su2double *Flow_Dir, Vel_Mag, Area, Flow_Dir_Mag;
  su2double UnitFlowDir[3] = {0.0,0.0,0.0};
  su2double Face_Flux, proj_vel;
  su2double *V_inlet = new su2double[nDim];
  su2double *V_Charac;
  su2double *V_domain;
  su2double *Coord_i, *Coord_j, dist_ij, delP, Pressure_j, Pressure_i;
  string Marker_Tag  = config->GetMarker_All_TagBound(val_marker);
  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool viscous        = config->GetViscous();

  su2double *Normal = new su2double[nDim];
  su2double *val_normal = new su2double[nDim];
  su2double *GridVel_i;

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

      /*--- Retrieve solution at the farfield boundary node ---*/

      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Retrieve the specified velocity for the inlet. ---*/

      Vel_Mag  = Inlet_Ptotal[val_marker][iVertex]/config->GetVelocity_Ref();

      Flow_Dir = Inlet_FlowDir[val_marker][iVertex];
      Flow_Dir_Mag = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Flow_Dir_Mag += Flow_Dir[iDim]*Flow_Dir[iDim];
      Flow_Dir_Mag = sqrt(Flow_Dir_Mag);

      /*--- Store the unit flow direction vector. ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        UnitFlowDir[iDim] = Flow_Dir[iDim]/Flow_Dir_Mag;

      /*--- Store the velocity in the primitive variable vector. ---*/

      for (iDim = 0; iDim < nDim; iDim++)
          V_inlet[iDim] = Vel_Mag*UnitFlowDir[iDim];

      /*--- Adding the grid velocity because SetVelocityOld sets the solution which is U = U_rel + U_grid. ---*/
      /*if (dynamic_grid) {
        GridVel_i = geometry->nodes->GetGridVel(iPoint); 
         for (iDim = 0; iDim < nDim; iDim++)
          V_inlet[iDim] += GridVel_i[iDim];
      }*/

      /*--- Update the CharacPrimVar for this vertex on inlet marker. ---*/
      /*--- This is necessary for the turbulent solver. ---*/
      V_Charac = GetCharacPrimVar(val_marker, iVertex);

      V_Charac[0] = nodes->GetPressure(iPoint);
      for (iDim = 0; iDim < nDim; iDim++)
          V_Charac[iDim+1] = V_domain[iDim];
      V_Charac[nDim+1] = nodes->GetDensity(iPoint);

      if (viscous) {
          V_Charac[nDim+2] = nodes->GetLaminarViscosity(iPoint);
          V_Charac[nDim+3] = nodes->GetEddyViscosity(iPoint);
      }

      /*--- Impose the value of the velocity as a strong boundary condition (Dirichlet).
       * Fix the velocity and remove any contribution to the residual at this node. 
       * Note that SetVelocityOld sets the solution which needs to be multiplied by rho. ---*/

      nodes->SetVelocity_Old(iPoint,V_inlet);

      nodes->SetStrongBC(iPoint);

      LinSysRes.SetBlock_Zero(iPoint);
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit) {
        for (iDim = 0; iDim < nDim; iDim++) {
          total_index = iPoint*nVar+iDim;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] Normal;
  delete [] V_inlet;
  delete [] val_normal;

}

void CPBIncEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim,iVar,jDim;
  unsigned long iVertex, iPoint, Point_Normal, total_index;
  su2double Area, yCoordRef, yCoord,proj_vel;
  su2double *V_outlet, *V_domain, P_Outlet, Face_Flux, Flux0;
  su2double *Coord_i, *Coord_j, dist_ij, delP, Pressure_j, Pressure_i;

  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool viscous       = config->GetViscous();
  string Marker_Tag  = config->GetMarker_All_TagBound(val_marker);

  su2double *Normal = new su2double[nDim];
  auto Kind_Outlet = config->GetKind_Inc_Outlet(Marker_Tag);
  auto turb_model = config->GetKind_Turb_Model();
  su2double *GridVel_j = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

      Area = GeometryToolbox::Norm(nDim, Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Current solution at this boundary node ---*/

      V_domain = nodes->GetPrimitive(iPoint);

      V_outlet = GetCharacPrimVar(val_marker, iVertex);

      for (iDim = 0; iDim < nDim; iDim++) GridVel_j[iDim] = 0.0;
      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), GridVel_j);

      switch (Kind_Outlet) {
        case INC_OUTLET_TYPE::PRESSURE_OUTLET: {

          /*--- Retrieve the specified back pressure for this outlet. ---*/
          P_Outlet = config->GetOutlet_Pressure(Marker_Tag)/config->GetPressure_Ref();

          nodes->SetPressure(iPoint, P_Outlet);
          V_outlet[0] = P_Outlet;
          for (iDim = 0; iDim < nDim; iDim++)
            V_outlet[iDim+1] = 0.0;
          V_outlet[nDim+1] = nodes->GetDensity(iPoint);

          conv_numerics->SetPrimitive(V_domain, V_outlet);

          /*--- Compute the residual using an upwind scheme ---*/

          auto residual = conv_numerics->ComputeResidual(config);

          for (iVar = 0; iVar < nVar; iVar++) {
            Residual[iVar] = 2.0*residual.residual[iVar];
          }

          for (iDim = 0; iDim < nDim; iDim++)
            for (jDim = 0; jDim < nDim; jDim++)
              Jacobian_i[iDim][jDim] = 2.0*residual.jacobian_i[iDim][jDim];

          /*--- Update residual value ---*/
          LinSysRes.AddBlock(iPoint, Residual);

          if (implicit)
            Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

          for (iDim = 0; iDim < nDim; iDim++)
            V_outlet[iDim+1] = V_domain[iDim+1];

          break;
          }

          case INC_OUTLET_TYPE::OPEN:
                   /* Not working yet */

          break;
      }

      if (viscous) {
        /*--- Set transport properties at the outlet. ---*/
        V_domain[nDim+2] = nodes->GetLaminarViscosity(iPoint);
        V_domain[nDim+3] = nodes->GetEddyViscosity(iPoint);
        V_outlet[nDim+2] = nodes->GetLaminarViscosity(iPoint);
        V_outlet[nDim+3] = nodes->GetEddyViscosity(iPoint);

        /*--- Set the normal vector and the coordinates ---*/
        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                                geometry->nodes->GetCoord(Point_Normal));

        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                          nodes->GetGradient_Primitive(iPoint));

        /*--- Turbulent kinetic energy ---*/
        if (turb_model == TURB_MODEL::SST)
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
  delete [] GridVel_j;

}

/*--- Note that velocity indices in residual are hard coded in solver_structure. Need to be careful. ---*/
// TODO: PBFlow: Already initialized in FVMFlowSolverBase
void CPBIncEulerSolver::BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                               CNumerics *numerics, CConfig *config) {

  /*--- Complete residuals for periodic boundary conditions. We loop over
   the periodic BCs in matching pairs so that, in the event that there are
   adjacent periodic markers, the repeated points will have their residuals
   accumulated correctly during the communications. For implicit calculations,
   the Jacobians and linear system are also correctly adjusted here. ---*/

  SetMomCoeffPer(geometry, solver_container, config);

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
  }
}

void CPBIncEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) { }


void CPBIncEulerSolver::PrintVerificationError(const CConfig *config) const { }
