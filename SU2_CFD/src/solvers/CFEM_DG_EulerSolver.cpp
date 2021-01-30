/*!
 * \file CFEM_DG_EulerSolver.cpp
 * \brief Main subroutines for solving finite element Euler flow problems
 * \author J. Alonso, E. van der Weide, T. Economon
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


#include "../../include/solvers/CFEM_DG_EulerSolver.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/fluid/CIdealGas.hpp"
#include "../../include/fluid/CVanDerWaalsGas.hpp"
#include "../../include/fluid/CPengRobinson.hpp"

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(void)
  : CFEM_DG_SolverBase() {}

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(CConfig        *config,
                                         unsigned short val_nDim,
                                         unsigned short iMesh)
  : CFEM_DG_SolverBase() {

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Dummy solver constructor that calls the SetNondim. routine in
        order to load the flow non-dim. information into the config class.
        This is needed to complete a partitioning for time-accurate local
        time stepping that depends on the flow state. ---*/
  nDim = val_nDim;
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  SetNondimensionalization(config, iMesh, false);
}

CFEM_DG_EulerSolver::CFEM_DG_EulerSolver(CGeometry      *geometry,
                                         CConfig        *config,
                                         unsigned short iMesh)
  : CFEM_DG_SolverBase(geometry, config, iMesh) {

  /*--- Set the gamma value ---*/
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

CFEM_DG_EulerSolver::~CFEM_DG_EulerSolver(void) {
}

void CFEM_DG_EulerSolver::SetNondimensionalization(CConfig        *config,
                                                   unsigned short iMesh,
                                                   const bool     writeOutput) {

  su2double Temperature_FreeStream   = 0.0, Mach2Vel_FreeStream   = 0.0, ModVel_FreeStream      = 0.0,
            Energy_FreeStream        = 0.0, ModVel_FreeStreamND   = 0.0, Velocity_Reynolds      = 0.0,
            Omega_FreeStream         = 0.0, Omega_FreeStreamND    = 0.0, Viscosity_FreeStream   = 0.0,
            Density_FreeStream       = 0.0, Pressure_FreeStream   = 0.0, Tke_FreeStream         = 0.0,
            Length_Ref               = 0.0, Density_Ref           = 0.0, Pressure_Ref           = 0.0,
            Velocity_Ref             = 0.0, Temperature_Ref       = 0.0, Time_Ref               = 0.0,
            Omega_Ref                = 0.0, Force_Ref             = 0.0, Gas_Constant_Ref       = 0.0,
            Viscosity_Ref            = 0.0, Conductivity_Ref      = 0.0, Energy_Ref             = 0.0,
            Froude                   = 0.0, Pressure_FreeStreamND = 0.0, Density_FreeStreamND   = 0.0,
            Temperature_FreeStreamND = 0.0, Gas_ConstantND        = 0.0, Viscosity_FreeStreamND = 0.0,
            Tke_FreeStreamND         = 0.0, Energy_FreeStreamND   = 0.0, Total_UnstTimeND       = 0.0,
            Delta_UnstTimeND         = 0.0;

  su2double Velocity_FreeStreamND[MAXNDIM] = {0.0};

  unsigned short iDim;

  /*--- Local variables ---*/
  su2double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  su2double Mach             = config->GetMach();
  su2double Reynolds         = config->GetReynolds();
  bool unsteady           = (config->GetTime_Marching() != NO);
  bool viscous            = config->GetViscous();
  bool grid_movement      = config->GetGrid_Movement();
  bool turbulent          = (config->GetKind_Solver() == FEM_RANS) || (config->GetKind_Solver() == FEM_LES);
  bool tkeNeeded          = ((turbulent) && ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)));
  bool free_stream_temp   = (config->GetKind_FreeStreamOption() == TEMPERATURE_FS);
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);

  /*--- Compute the Free Stream velocity, using the Mach number ---*/
  Pressure_FreeStream = config->GetPressure_FreeStream();
  Density_FreeStream  = config->GetDensity_FreeStream();
  Temperature_FreeStream  = config->GetTemperature_FreeStream();

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:

      if (config->GetSystemMeasurements() == SI) config->SetGas_Constant(287.058);
      else if (config->GetSystemMeasurements() == US) config->SetGas_Constant(1716.49);

      FluidModel = new CIdealGas(1.4, config->GetGas_Constant(), config->GetCompute_Entropy());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case IDEAL_GAS:

      FluidModel = new CIdealGas(Gamma, config->GetGas_Constant(), config->GetCompute_Entropy());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case VW_GAS:

      FluidModel = new CVanDerWaalsGas(Gamma, config->GetGas_Constant(),
                                       config->GetPressure_Critical(), config->GetTemperature_Critical());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case PR_GAS:

      FluidModel = new CPengRobinson(Gamma, config->GetGas_Constant(), config->GetPressure_Critical(),
                                     config->GetTemperature_Critical(), config->GetAcentric_Factor());
      if (free_stream_temp) {
        FluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = FluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        FluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = FluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

  }

  Mach2Vel_FreeStream = FluidModel->GetSoundSpeed();

  /*--- Compute the Free Stream velocity, using the Mach number ---*/
  if (nDim == 2) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
  }
  if (nDim == 3) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
  }

  /*--- Compute the modulus of the free stream velocity ---*/
  ModVel_FreeStream = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);

  /*--- Viscous initialization ---*/
  if (viscous) {

    /*--- The dimensional viscosity is needed to determine the free-stream conditions.
          To accomplish this, simply set the non-dimensional coefficients to the
          dimensional ones. This will be overruled later.---*/
    config->SetMu_RefND(config->GetMu_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref());
    config->SetMu_SND(config->GetMu_S());

    config->SetMu_ConstantND(config->GetMu_Constant());

    /*--- Reynolds based initialization ---*/
    if (reynolds_init) {

      /*--- First, check if there is mesh motion. If yes, use the Mach
       number relative to the body to initialize the flow. ---*/
      if (grid_movement) Velocity_Reynolds = config->GetMach_Motion()*Mach2Vel_FreeStream;
      else Velocity_Reynolds = ModVel_FreeStream;

      /*--- For viscous flows, pressure will be computed from a density
            that is found from the Reynolds number. The viscosity is computed
            from the dimensional version of Sutherland's law or the constant
            viscosity, depending on the input option.---*/
      FluidModel->SetLaminarViscosityModel(config);

      Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);

      Density_FreeStream = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds*config->GetLength_Reynolds());
      config->SetDensity_FreeStream(Density_FreeStream);
      FluidModel->SetTDState_rhoT(Density_FreeStream, Temperature_FreeStream);
      Pressure_FreeStream = FluidModel->GetPressure();
      config->SetPressure_FreeStream(Pressure_FreeStream);
      Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*--- Thermodynamics quantities based initialization ---*/
    else {

      FluidModel->SetLaminarViscosityModel(config);
      Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);
      Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*--- Turbulence kinetic energy ---*/
    Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());

  }
  else {

    /*--- For inviscid flow, energy is calculated from the specified
     FreeStream quantities using the proper gas law. ---*/
    Energy_FreeStream = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

  }

  /*-- Compute the freestream energy. ---*/
  if (tkeNeeded) { Energy_FreeStream += Tke_FreeStream; }; config->SetEnergy_FreeStream(Energy_FreeStream);

  /*--- Compute non dimensional quantities. By definition,
   Lref is one because we have converted the grid to meters. ---*/
  if (config->GetRef_NonDim() == DIMENSIONAL) {
    Pressure_Ref      = 1.0;
    Density_Ref       = 1.0;
    Temperature_Ref   = 1.0;
  }
  else if (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
    Pressure_Ref      = Pressure_FreeStream;     // Pressure_FreeStream = 1.0
    Density_Ref       = Density_FreeStream;      // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;  // Temperature_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
    Pressure_Ref      = Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
    Pressure_Ref      = Mach*Mach*Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/(Gamma*(M_inf)^2)
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  config->SetPressure_Ref(Pressure_Ref);
  config->SetDensity_Ref(Density_Ref);
  config->SetTemperature_Ref(Temperature_Ref);

  Length_Ref        = 1.0;                                                         config->SetLength_Ref(Length_Ref);
  Velocity_Ref      = sqrt(config->GetPressure_Ref()/config->GetDensity_Ref());    config->SetVelocity_Ref(Velocity_Ref);
  Time_Ref          = Length_Ref/Velocity_Ref;                                     config->SetTime_Ref(Time_Ref);
  Omega_Ref         = Velocity_Ref/Length_Ref;                                     config->SetOmega_Ref(Omega_Ref);
  Force_Ref         = Velocity_Ref*Velocity_Ref/Length_Ref;                        config->SetForce_Ref(Force_Ref);
  Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/config->GetTemperature_Ref();      config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref     = config->GetDensity_Ref()*Velocity_Ref*Length_Ref;            config->SetViscosity_Ref(Viscosity_Ref);
  Conductivity_Ref  = Viscosity_Ref*Gas_Constant_Ref;                              config->SetConductivity_Ref(Conductivity_Ref);
  Froude            = ModVel_FreeStream/sqrt(STANDARD_GRAVITY*Length_Ref);         config->SetFroude(Froude);

  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/
  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();  config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();    config->SetDensity_FreeStreamND(Density_FreeStreamND);

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }

  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);

  Gas_ConstantND = config->GetGas_Constant()/Gas_Constant_Ref;    config->SetGas_ConstantND(Gas_ConstantND);


  ModVel_FreeStreamND = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
  ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND); config->SetModVel_FreeStreamND(ModVel_FreeStreamND);

  Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;   config->SetViscosity_FreeStreamND(Viscosity_FreeStreamND);

  Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStream(Tke_FreeStream);

  Tke_FreeStreamND  = 3.0/2.0*(ModVel_FreeStreamND*ModVel_FreeStreamND*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());
  config->SetTke_FreeStreamND(Tke_FreeStreamND);

  Omega_FreeStream = Density_FreeStream*Tke_FreeStream/max((Viscosity_FreeStream*config->GetTurb2LamViscRatio_FreeStream()), 1.e-25);
  config->SetOmega_FreeStream(Omega_FreeStream);

  Omega_FreeStreamND = Density_FreeStreamND*Tke_FreeStreamND/max((Viscosity_FreeStreamND*config->GetTurb2LamViscRatio_FreeStream()), 1.e-25);
  config->SetOmega_FreeStreamND(Omega_FreeStreamND);

  /*--- Initialize the dimensionless Fluid Model that will be used to solve the dimensionless problem ---*/

  /*--- Delete the original (dimensional) FluidModel object before replacing. ---*/
  delete FluidModel;

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:
      FluidModel = new CIdealGas(1.4, Gas_ConstantND, config->GetCompute_Entropy());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

    case IDEAL_GAS:
      FluidModel = new CIdealGas(Gamma, Gas_ConstantND, config->GetCompute_Entropy());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

    case VW_GAS:
      FluidModel = new CVanDerWaalsGas(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                       config->GetTemperature_Critical()/config->GetTemperature_Ref());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

    case PR_GAS:
      FluidModel = new CPengRobinson(Gamma, Gas_ConstantND, config->GetPressure_Critical() /config->GetPressure_Ref(),
                                     config->GetTemperature_Critical()/config->GetTemperature_Ref(), config->GetAcentric_Factor());
      FluidModel->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
      break;

  }

  Energy_FreeStreamND = FluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;

  if (viscous) {

    /*--- Constant viscosity model ---*/
    config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);

    /*--- Sutherland's model ---*/

    config->SetMu_RefND(config->GetMu_Ref()/Viscosity_Ref);
    config->SetMu_SND(config->GetMu_S()/config->GetTemperature_Ref());
    config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref()/config->GetTemperature_Ref());

    /* constant thermal conductivity model */
    config->SetKt_ConstantND(config->GetKt_Constant()/Conductivity_Ref);

    FluidModel->SetLaminarViscosityModel(config);
    FluidModel->SetThermalConductivityModel(config);

  }

  if (tkeNeeded) { Energy_FreeStreamND += Tke_FreeStreamND; };  config->SetEnergy_FreeStreamND(Energy_FreeStreamND);

  Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);

  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);

  /*--- Write output to the console if this is required and if this is the master node and first domain ---*/
  if ((rank == MASTER_NODE) && (iMesh == MESH_0) && writeOutput) {

    cout.precision(6);

    if (viscous) {
      cout << "Viscous flow: Computing pressure using the ideal gas law" << endl;
      cout << "based on the free-stream temperature and a density computed" << endl;
      cout << "from the Reynolds number." << endl;
    } else {
      cout << "Inviscid flow: Computing density based on free-stream" << endl;
      cout << "temperature and pressure using the ideal gas law." << endl;
    }

    if (grid_movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
    else cout << "Force coefficients computed using free-stream values." << endl;

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

    if (viscous) {

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

      }
    } else {
      ModelTable << "-" << "-";
    }

    if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
    else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
    NonDimTable << "Gas Constant" << config->GetGas_Constant() << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
    else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
    NonDimTable << "Spec. Heat Ratio" << "-" << "-" << "-" << Gamma;
    Unit.str("");

    switch(config->GetKind_FluidModel()){
    case STANDARD_AIR:
      ModelTable << "STANDARD_AIR";
      break;
    case IDEAL_GAS:
      ModelTable << "IDEAL_GAS";
      break;
    case VW_GAS:
      ModelTable << "VW_GAS";
      break;
    case PR_GAS:
      ModelTable << "PR_GAS";
      break;
    }

    if (config->GetKind_FluidModel() == VW_GAS || config->GetKind_FluidModel() == PR_GAS){
        NonDimTable << "Critical Pressure" << config->GetPressure_Critical() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_Critical() /config->GetPressure_Ref();
        Unit.str("");
        Unit << "K";
        NonDimTable << "Critical Temperature" << config->GetTemperature_Critical() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_Critical() /config->GetTemperature_Ref();
        Unit.str("");
    }
    NonDimTable.PrintFooter();

    NonDimTableOut <<"-- Initial and free-stream conditions:"<< endl;

    NonDimTable.PrintHeader();

    if      (config->GetSystemMeasurements() == SI) Unit << "Pa";
    else if (config->GetSystemMeasurements() == US) Unit << "psf";
    NonDimTable << "Static Pressure" << config->GetPressure_FreeStream() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "kg/m^3";
    else if (config->GetSystemMeasurements() == US) Unit << "slug/ft^3";
    NonDimTable << "Density" << config->GetDensity_FreeStream() << config->GetDensity_Ref() << Unit.str() << config->GetDensity_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "K";
    else if (config->GetSystemMeasurements() == US) Unit << "R";
    NonDimTable << "Temperature" << config->GetTemperature_FreeStream() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "m^2/s^2";
    else if (config->GetSystemMeasurements() == US) Unit << "ft^2/s^2";
    NonDimTable << "Total Energy" << config->GetEnergy_FreeStream() << config->GetEnergy_Ref() << Unit.str() << config->GetEnergy_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "m/s";
    else if (config->GetSystemMeasurements() == US) Unit << "ft/s";
    NonDimTable << "Velocity-X" << config->GetVelocity_FreeStream()[0] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[0];
    NonDimTable << "Velocity-Y" << config->GetVelocity_FreeStream()[1] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[1];
    if (nDim == 3) {
      NonDimTable << "Velocity-Z" << config->GetVelocity_FreeStream()[2] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[2];
    }
    NonDimTable << "Velocity Magnitude" << config->GetModVel_FreeStream() << config->GetVelocity_Ref() << Unit.str() << config->GetModVel_FreeStreamND();
    Unit.str("");

    if (viscous) {
      NonDimTable.PrintFooter();
      if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
      else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
      NonDimTable << "Viscosity" << config->GetViscosity_FreeStream() << config->GetViscosity_Ref() << Unit.str() << config->GetViscosity_FreeStreamND();
      Unit.str("");
      if      (config->GetSystemMeasurements() == SI) Unit << "W/m^2.K";
      else if (config->GetSystemMeasurements() == US) Unit << "lbf/ft.s.R";
      NonDimTable << "Conductivity" << "-" << config->GetConductivity_Ref() << Unit.str() << "-";
      Unit.str("");
      if (turbulent) {
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
    if (viscous) {
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

void CFEM_DG_EulerSolver::DetermineGraphDOFs(const CMeshFEM_DG *FEMGeometry,
                                             CConfig           *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::MetaDataJacobianComputation(const CMeshFEM_DG *FEMGeometry,
                                                      const vector<int> &colorLocalDOFs) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::SetUpTaskList(CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Prepare_MPI_Communication(const CMeshFEM_DG *FEMGeometry,
                                                    CConfig           *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Initiate_MPI_Communication(CConfig *config,
                                                     const unsigned short timeLevel) {
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

bool CFEM_DG_EulerSolver::Complete_MPI_Communication(CConfig *config,
                                                     const unsigned short timeLevel,
                                                     const bool commMustBeCompleted) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  return false;
}

void CFEM_DG_EulerSolver::Initiate_MPI_ReverseCommunication(CConfig *config,
                                                            const unsigned short timeLevel) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

bool CFEM_DG_EulerSolver::Complete_MPI_ReverseCommunication(CConfig *config,
                                                            const unsigned short timeLevel,
                                                            const bool commMustBeCompleted) {
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  return false;
}

void CFEM_DG_EulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iStep, unsigned short RunTime_EqSystem, bool Output) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                      unsigned short iMesh) { }

void CFEM_DG_EulerSolver::ComputeSpatialJacobian(CGeometry *geometry,  CSolver **solver_container,
                                                 CNumerics **numerics, CConfig *config,
                                                 unsigned short iMesh, unsigned short RunTime_EqSystem) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Set_OldSolution() {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Set_NewSolution() {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::CheckTimeSynchronization(CConfig         *config,
                                                   const su2double TimeSync,
                                                   su2double       &timeEvolved,
                                                   bool            &syncTimeReached) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ProcessTaskList_DG(CGeometry *geometry,  CSolver **solver_container,
                                             CNumerics **numerics, CConfig *config,
                                             unsigned short iMesh) {
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_SpaceTimeIntegration(CGeometry *geometry,  CSolver **solver_container,
                                                    CNumerics **numerics, CConfig *config,
                                                    unsigned short iMesh, unsigned short RunTime_EqSystem) {
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::TolerancesADERPredictorStep(void) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_PredictorStep(CConfig             *config,
                                                const unsigned long elemBeg,
                                                const unsigned long elemEnd,
                                                su2double           *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_AliasedPredictorResidual_2D(CConfig              *config,
                                                              CVolumeElementFEM_DG *elem,
                                                              const su2double      *sol,
                                                              const unsigned short nSimul,
                                                              const unsigned short NPad,
                                                              su2double            *res,
                                                              su2double            *work) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_AliasedPredictorResidual_3D(CConfig              *config,
                                                              CVolumeElementFEM_DG *elem,
                                                              const su2double      *sol,
                                                              const unsigned short nSimul,
                                                              const unsigned short NPad,
                                                              su2double            *res,
                                                              su2double            *work) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_NonAliasedPredictorResidual_2D(CConfig              *config,
                                                                 CVolumeElementFEM_DG *elem,
                                                                 const su2double      *sol,
                                                                 const unsigned short nSimul,
                                                                 const unsigned short NPad,
                                                                 su2double            *res,
                                                                 su2double            *work) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_NonAliasedPredictorResidual_3D(CConfig              *config,
                                                                 CVolumeElementFEM_DG *elem,
                                                                 const su2double      *sol,
                                                                 const unsigned short nSimul,
                                                                 const unsigned short NPad,
                                                                 su2double            *res,
                                                                 su2double            *work) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_TimeInterpolatePredictorSol(CConfig             *config,
                                                              const unsigned short iTime,
                                                              const unsigned long  elemBeg,
                                                              const unsigned long  elemEnd,
                                                              const unsigned long  nAdjElem,
                                                              const unsigned long  *adjElem,
                                                              const bool           secondPartTimeInt,
                                                              su2double            *solTimeLevel) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Shock_Capturing_DG(CConfig             *config,
                                             const unsigned long elemBeg,
                                             const unsigned long elemEnd,
                                             su2double           *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Volume_Residual(CConfig             *config,
                                          const unsigned long elemBeg,
                                          const unsigned long elemEnd,
                                          su2double           *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Boundary_Conditions(const unsigned short timeLevel,
                                              CConfig              *config,
                                              CNumerics            **numerics,
                                              const bool           haloInfoNeededForBC,
                                              su2double            *workArray){

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ResidualFaces(CConfig             *config,
                                        const unsigned long indFaceBeg,
                                        const unsigned long indFaceEnd,
                                        unsigned long       &indResFaces,
                                        CNumerics           *numerics,
                                        su2double           *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::InviscidFluxesInternalMatchingFace(
                                              CConfig              *config,
                                              const unsigned long  lBeg,
                                              const unsigned long  lEnd,
                                              const unsigned short NPad,
                                              su2double            *solIntL,
                                              su2double            *solIntR,
                                              su2double            *fluxes,
                                              CNumerics            *numerics) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::AccumulateSpaceTimeResidualADEROwnedElem(
                                                     CConfig             *config,
                                                     const unsigned short timeLevel,
                                                     const unsigned short intPoint) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::AccumulateSpaceTimeResidualADERHaloElem(
                                                     CConfig             *config,
                                                     const unsigned short timeLevel,
                                                     const unsigned short intPoint) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::CreateFinalResidual(const unsigned short timeLevel,
                                              const bool ownedElements) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::MultiplyResidualByInverseMassMatrix(
                                              CConfig            *config,
                                              const bool          useADER,
                                              const unsigned long elemBeg,
                                              const unsigned long elemEnd,
                                              su2double           *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::Pressure_Forces(const CGeometry *geometry, const CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                               CConfig *config, unsigned short iRKStep) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                                 CConfig *config, unsigned short iRKStep) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::SetResidual_RMS_FEM(CGeometry *geometry,
                                              CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ComputeVerificationError(CGeometry *geometry,
                                                   CConfig   *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ADER_DG_Iteration(const unsigned long elemBeg,
                                            const unsigned long elemEnd) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BoundaryStates_Euler_Wall(CConfig                  *config,
                                                    const unsigned short     nFaceSimul,
                                                    const unsigned short     NPad,
                                                    const CSurfaceElementFEM *surfElem,
                                                    const su2double          *solIntL,
                                                    su2double                *solIntR) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BoundaryStates_Inlet(CConfig                  *config,
                                               const unsigned short     nFaceSimul,
                                               const unsigned short     NPad,
                                               const CSurfaceElementFEM *surfElem,
                                               unsigned short           val_marker,
                                               const su2double          *solIntL,
                                               su2double                *solIntR) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BoundaryStates_Outlet(CConfig                  *config,
                                                const unsigned short     nFaceSimul,
                                                const unsigned short     NPad,
                                                const CSurfaceElementFEM *surfElem,
                                                unsigned short           val_marker,
                                                const su2double          *solIntL,
                                                su2double                *solIntR) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BoundaryStates_Riemann(CConfig                  *config,
                                                 const unsigned short     nFaceSimul,
                                                 const unsigned short     NPad,
                                                 const CSurfaceElementFEM *surfElem,
                                                 unsigned short           val_marker,
                                                 const su2double          *solIntL,
                                                 su2double                *solIntR) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Euler_Wall(CConfig                  *config,
                                        const unsigned long      surfElemBeg,
                                        const unsigned long      surfElemEnd,
                                        const CSurfaceElementFEM *surfElem,
                                        su2double                *resFaces,
                                        CNumerics                *conv_numerics,
                                        su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Far_Field(CConfig                  *config,
                                       const unsigned long      surfElemBeg,
                                       const unsigned long      surfElemEnd,
                                       const CSurfaceElementFEM *surfElem,
                                       su2double                *resFaces,
                                       CNumerics                *conv_numerics,
                                       su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Sym_Plane(CConfig                  *config,
                                       const unsigned long      surfElemBeg,
                                       const unsigned long      surfElemEnd,
                                       const CSurfaceElementFEM *surfElem,
                                       su2double                *resFaces,
                                       CNumerics                *conv_numerics,
                                       su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Supersonic_Outlet(CConfig                  *config,
                                               const unsigned long      surfElemBeg,
                                               const unsigned long      surfElemEnd,
                                               const CSurfaceElementFEM *surfElem,
                                               su2double                *resFaces,
                                               CNumerics                *conv_numerics,
                                               su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Inlet(CConfig                  *config,
                                   const unsigned long      surfElemBeg,
                                   const unsigned long      surfElemEnd,
                                   const CSurfaceElementFEM *surfElem,
                                   su2double                *resFaces,
                                   CNumerics                *conv_numerics,
                                   unsigned short           val_marker,
                                   su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Outlet(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    unsigned short           val_marker,
                                    su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Riemann(CConfig                  *config,
                                     const unsigned long      surfElemBeg,
                                     const unsigned long      surfElemEnd,
                                     const CSurfaceElementFEM *surfElem,
                                     su2double                *resFaces,
                                     CNumerics                *conv_numerics,
                                     unsigned short           val_marker,
                                     su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::BC_Custom(CConfig                  *config,
                                    const unsigned long      surfElemBeg,
                                    const unsigned long      surfElemEnd,
                                    const CSurfaceElementFEM *surfElem,
                                    su2double                *resFaces,
                                    CNumerics                *conv_numerics,
                                    su2double                *workArray) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ResidualInviscidBoundaryFace(
                                      CConfig                  *config,
                                      const unsigned short     nFaceSimul,
                                      const unsigned short     NPad,
                                      CNumerics                *conv_numerics,
                                      const CSurfaceElementFEM *surfElem,
                                      su2double                *solInt0,
                                      su2double                *solInt1,
                                      su2double                *fluxes,
                                      su2double                *resFaces,
                                      unsigned long            &indResFaces) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::LeftStatesIntegrationPointsBoundaryFace(
                                             CConfig                  *config,
                                             const unsigned short     nFaceSimul,
                                             const unsigned short     NPad,
                                             const CSurfaceElementFEM *surfElem,
                                             su2double                *solFace,
                                             su2double                *solIntL) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::ComputeInviscidFluxesFace(CConfig              *config,
                                                    const unsigned short nFaceSimul,
                                                    const unsigned short NPad,
                                                    const unsigned long  nPoints,
                                                    const su2double      *normalsFace[],
                                                    const su2double      *gridVelsFace[],
                                                    const su2double      *solL,
                                                    const su2double      *solR,
                                                    su2double            *fluxes,
                                                    CNumerics            *numerics) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CFEM_DG_EulerSolver::LoadRestart(CGeometry **geometry,
                                      CSolver   ***solver,
                                      CConfig   *config,
                                      int       val_iter,
                                      bool      val_update_geo) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}
