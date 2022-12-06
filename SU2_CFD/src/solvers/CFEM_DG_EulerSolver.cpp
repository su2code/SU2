/*!
 * \file CFEM_DG_EulerSolver.cpp
 * \brief Main subroutines for solving finite element Euler flow problems
 * \author J. Alonso, E. van der Weide, T. Economon
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../include/fluid/CIdealGas.hpp"
#include "../../include/fluid/CVanDerWaalsGas.hpp"
#include "../../include/fluid/CPengRobinson.hpp"
#include "../../include/fluid/CCoolProp.hpp"

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

  /*--- Set the gamma value and the number of variables. ---*/
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  nVar = nDim + 2;

  /*--- Define some auxiliary vectors related to the residual ---*/
  Residual_RMS.resize(nVar,1.e-35);
  Residual_Max.resize(nVar,1.e-35);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

  /*--- Perform the non-dimensionalization for the flow equations using the
        specified reference values. ---*/
  SetNondimensionalization(config, iMesh, true);

  /*--- Check if we are executing a verification case. If so, the
        VerificationSolution object will be instantiated for a particular
        option from the available library of verification solutions. Note
        that this is done after SetNondim(), as problem-specific initial
        parameters are needed by the solution constructors. ---*/
  SetVerificationSolution(nDim, nVar, config);

  /*--- Read farfield conditions ---*/
  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Energy_Inf      = config->GetEnergy_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Mach_Inf        = config->GetMach();

  /*--- Set the conservative variables of the free-stream. ---*/
  ConsVarFreeStream.resize(nVar);

  ConsVarFreeStream[0]      = Density_Inf;
  ConsVarFreeStream[nVar-1] = Density_Inf*Energy_Inf;

  for(unsigned short iDim=0; iDim<nDim; ++iDim)
    ConsVarFreeStream[iDim+1] = Density_Inf*Velocity_Inf[iDim];

  /*--- Set the entropy variables of the free-stream. ---*/
  EntropyVarFreeStream.resize(nVar);

  const su2double ovgm1 = 1.0/Gamma_Minus_One;
  const su2double s     = log(Pressure_Inf/pow(Density_Inf,Gamma));
  const su2double pInv  = 1.0/Pressure_Inf;

  su2double V2_Inf = 0.0;
  for(unsigned short iDim=0; iDim<nDim; ++iDim) {
    EntropyVarFreeStream[iDim+1] = Density_Inf*Velocity_Inf[iDim]*pInv;
    V2_Inf                      += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }

  EntropyVarFreeStream[0]      =  (Gamma-s)*ovgm1 - 0.5*pInv*Density_Inf*V2_Inf;
  EntropyVarFreeStream[nVar-1] = -pInv*Density_Inf;

  /*--- A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain ---*/
  dynamic_grid = config->GetDynamic_Grid();

  /*--- Start of the parallel region. ---*/
  SU2_OMP_PARALLEL
  {
    /*--- Parallel loop over the volume elements. ---*/
#ifdef HAVE_OMP
    const size_t omp_chunk_size_vol = computeStaticChunkSize(nVolElemTot, omp_get_num_threads(), 64);
#endif
    SU2_OMP_FOR_STAT(omp_chunk_size_vol)
    for(unsigned long i=0; i<nVolElemTot; ++i) {

      /*--- Allocate the memory to store the volume solution(s)
            and the residuals. ---*/
      volElem[i].AllocateCompressibleFlowVar(config, nVar);
      volElem[i].AllocateResiduals(config, nVar);

      /*--- Initialize the solution to a uniform flow. This is overruled
            when a restart takes place. ---*/
      volElem[i].SetConstantSolution(EntropyVarFreeStream.data(), nVar, 0);
    }
    END_SU2_OMP_FOR

    /*--- Parallel loop over the interior faces. ---*/
    const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
    const unsigned long  nFaces      = nMatchingInternalFacesWithHaloElem[nTimeLevels];
#ifdef HAVE_OMP
    const size_t omp_chunk_size_face = computeStaticChunkSize(nFaces, omp_get_num_threads(), 64);
#endif
    SU2_OMP_FOR_STAT(omp_chunk_size_face)
    for(unsigned long i=0; i<nFaces; ++i) {

      /*--- Allocate the memory to store the residuals. ---*/
      matchingInternalFaces[i].AllocateResiduals(config, nVar);
    }
    END_SU2_OMP_FOR

    /*--- Loop over the physical boundaries. ---*/
    for(unsigned short iMarker=0; iMarker<config->GetnMarker_All(); iMarker++) {
      if(config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY) {
        SU2_OMP_FOR_STAT(omp_chunk_size_face)
        for(unsigned long i=0; i<boundaries[iMarker].surfElem.size(); ++i) {
          boundaries[iMarker].surfElem[i].AllocateResiduals(config, nVar);
        }
        END_SU2_OMP_FOR
      }
    }
  }
  END_SU2_OMP_PARALLEL

  /*--- Determine the communication pattern. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);
  if( !DGGeometry) SU2_MPI::Error(string("Dynamic cast failed"), CURRENT_FUNCTION);
  Prepare_MPI_Communication(DGGeometry, config);

  /*--- Set up the list of tasks to be carried out in the computational expensive
        part of the code. For the Runge-Kutta schemes this is typically the
        computation of the spatial residual, while for ADER this list contains
        the tasks to be done for one space time step. ---*/
  SetUpTaskList(config);
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
  bool unsteady           = (config->GetTime_Marching() != TIME_MARCHING::STEADY);
  bool viscous            = config->GetViscous();
  bool grid_movement      = config->GetGrid_Movement();
  bool turbulent          = (config->GetKind_Solver() == MAIN_SOLVER::FEM_RANS) || (config->GetKind_Solver() == MAIN_SOLVER::FEM_LES);
  bool tkeNeeded          = ((turbulent) && ((config->GetKind_Turb_Model() == TURB_MODEL::SST) ));
  bool free_stream_temp   = (config->GetKind_FreeStreamOption() == FREESTREAM_OPTION::TEMPERATURE_FS);
  bool reynolds_init      = (config->GetKind_InitOption() == REYNOLDS);

  /*--- Compute the Free Stream velocity, using the Mach number ---*/
  Pressure_FreeStream = config->GetPressure_FreeStream();
  Density_FreeStream  = config->GetDensity_FreeStream();
  Temperature_FreeStream  = config->GetTemperature_FreeStream();

  CFluidModel* auxFluidModel = nullptr;

  /*--- The dimensional viscosity is needed to determine the free-stream conditions.
        To accomplish this, simply set the non-dimensional coefficients to the
        dimensional ones. This will be overruled later.---*/

  config->SetTemperature_Ref(1.0);
  config->SetViscosity_Ref(1.0);
  config->SetConductivity_Ref(1.0);

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:

      if (config->GetSystemMeasurements() == SI) config->SetGas_Constant(287.058);
      else if (config->GetSystemMeasurements() == US) config->SetGas_Constant(1716.49);

      auxFluidModel = new CIdealGas(1.4, config->GetGas_Constant(), config->GetCompute_Entropy());
      if (free_stream_temp) {
        auxFluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = auxFluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        auxFluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = auxFluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case IDEAL_GAS:

      auxFluidModel = new CIdealGas(Gamma, config->GetGas_Constant(), config->GetCompute_Entropy());
      if (free_stream_temp) {
        auxFluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = auxFluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        auxFluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = auxFluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case VW_GAS:

      auxFluidModel = new CVanDerWaalsGas(Gamma, config->GetGas_Constant(),
                                          config->GetPressure_Critical(), config->GetTemperature_Critical());
      if (free_stream_temp) {
        auxFluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = auxFluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        auxFluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = auxFluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case PR_GAS:

      auxFluidModel = new CPengRobinson(Gamma, config->GetGas_Constant(), config->GetPressure_Critical(),
                                        config->GetTemperature_Critical(), config->GetAcentric_Factor());
      if (free_stream_temp) {
        auxFluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = auxFluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        auxFluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = auxFluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

    case COOLPROP:

      auxFluidModel = new CCoolProp(config->GetFluid_Name());
      if (free_stream_temp) {
        auxFluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
        Density_FreeStream = auxFluidModel->GetDensity();
        config->SetDensity_FreeStream(Density_FreeStream);
      }
      else {
        auxFluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
        Temperature_FreeStream = auxFluidModel->GetTemperature();
        config->SetTemperature_FreeStream(Temperature_FreeStream);
      }
      break;

  }

  Mach2Vel_FreeStream = auxFluidModel->GetSoundSpeed();

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
      auxFluidModel->SetLaminarViscosityModel(config);

      Viscosity_FreeStream = auxFluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);

      Density_FreeStream = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds*config->GetLength_Reynolds());
      config->SetDensity_FreeStream(Density_FreeStream);
      auxFluidModel->SetTDState_rhoT(Density_FreeStream, Temperature_FreeStream);
      Pressure_FreeStream = auxFluidModel->GetPressure();
      config->SetPressure_FreeStream(Pressure_FreeStream);
      Energy_FreeStream = auxFluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*--- Thermodynamics quantities based initialization ---*/
    else {

      auxFluidModel->SetLaminarViscosityModel(config);
      Viscosity_FreeStream = auxFluidModel->GetLaminarViscosity();
      config->SetViscosity_FreeStream(Viscosity_FreeStream);
      Energy_FreeStream = auxFluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

    }

    /*--- Turbulence kinetic energy ---*/
    Tke_FreeStream  = 3.0/2.0*(ModVel_FreeStream*ModVel_FreeStream*config->GetTurbulenceIntensity_FreeStream()*config->GetTurbulenceIntensity_FreeStream());

  }
  else {

    /*--- For inviscid flow, energy is calculated from the specified
     FreeStream quantities using the proper gas law. ---*/
    Energy_FreeStream = auxFluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

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

  /*--- Auxilary (dimensional) FluidModel no longer needed. ---*/
  delete auxFluidModel;

  /*--- Create one final fluid model object per OpenMP thread to be able to use them in parallel.
   *    GetFluidModel() should be used to automatically access the "right" object of each thread. ---*/

  assert(FluidModel.empty() && "Potential memory leak!");
  FluidModel.resize(omp_get_max_threads());

  SU2_OMP_PARALLEL
  {
    const int thread = omp_get_thread_num();

    switch (config->GetKind_FluidModel()) {

      case STANDARD_AIR:
        FluidModel[thread] = new CIdealGas(1.4, Gas_ConstantND, config->GetCompute_Entropy());
        break;

      case IDEAL_GAS:
        FluidModel[thread] = new CIdealGas(Gamma, Gas_ConstantND, config->GetCompute_Entropy());
        break;

      case VW_GAS:
        FluidModel[thread] = new CVanDerWaalsGas(Gamma, Gas_ConstantND,
                                                 config->GetPressure_Critical() /config->GetPressure_Ref(),
                                                 config->GetTemperature_Critical()/config->GetTemperature_Ref());
        break;

      case PR_GAS:
        FluidModel[thread] = new CPengRobinson(Gamma, Gas_ConstantND,
                                               config->GetPressure_Critical() /config->GetPressure_Ref(),
                                               config->GetTemperature_Critical()/config->GetTemperature_Ref(),
                                               config->GetAcentric_Factor());
        break;

      case COOLPROP:
        FluidModel[thread] = new CCoolProp(config->GetFluid_Name());
        break;
    }

    GetFluidModel()->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);
    if (viscous) {
      GetFluidModel()->SetLaminarViscosityModel(config);
      GetFluidModel()->SetThermalConductivityModel(config);
    }

  }
  END_SU2_OMP_PARALLEL

  Energy_FreeStreamND = GetFluidModel()->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;

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
      case VISCOSITYMODEL::CONSTANT:
        ModelTable << "CONSTANT_VISCOSITY";
        if      (config->GetSystemMeasurements() == SI) Unit << "N.s/m^2";
        else if (config->GetSystemMeasurements() == US) Unit << "lbf.s/ft^2";
        NonDimTable << "Viscosity" << config->GetMu_Constant() << config->GetMu_Constant()/config->GetMu_ConstantND() << Unit.str() << config->GetMu_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      case VISCOSITYMODEL::COOLPROP:
        ModelTable << "COOLPROP";
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

      default:
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

      case CONDUCTIVITYMODEL::COOLPROP:
        ModelTable << "COOLPROP";
        Unit << "W/m^2.K";
        NonDimTable << "Molecular Cond." << "--" << "--" << Unit.str() << config->GetThermal_Conductivity_ConstantND();
        Unit.str("");
        NonDimTable.PrintFooter();
        break;

      default:
        break;

      }
    } else {
      ModelTable << "-" << "-";
    }

    if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
    else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
    if (config->GetKind_FluidModel() == COOLPROP) {
      CCoolProp auxFluidModel(config->GetFluid_Name());
      NonDimTable << "Gas Constant" << auxFluidModel.GetGas_Constant() << config->GetGas_Constant_Ref()
                  << Unit.str() << auxFluidModel.GetGas_Constant()/config->GetGas_Constant_Ref();
    }
    else {
        NonDimTable << "Gas Constant" << config->GetGas_Constant() << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();
    }
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
    else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
    if (config->GetKind_FluidModel() == COOLPROP) {
      NonDimTable << "Spec. Heat Ratio" << "-" << "-" << "-" << "-";
    }
    else {
      NonDimTable << "Spec. Heat Ratio" << "-" << "-" << "-" << Gamma;
    }
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
    case COOLPROP:
      ModelTable << "CoolProp library";
      break;
    }

    if (config->GetKind_FluidModel() == VW_GAS || config->GetKind_FluidModel() == PR_GAS){
        NonDimTable << "Critical Pressure" << config->GetPressure_Critical() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_Critical() /config->GetPressure_Ref();
        Unit.str("");
        Unit << "K";
        NonDimTable << "Critical Temperature" << config->GetTemperature_Critical() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_Critical() /config->GetTemperature_Ref();
        Unit.str("");
    }
    if (config->GetKind_FluidModel() == COOLPROP) {
        CCoolProp auxFluidModel(config->GetFluid_Name());
        NonDimTable << "Critical Pressure" << auxFluidModel.GetPressure_Critical() << config->GetPressure_Ref() << Unit.str() << auxFluidModel.GetPressure_Critical() /config->GetPressure_Ref();
        Unit.str("");
        Unit << "K";
        NonDimTable << "Critical Temperature" << auxFluidModel.GetTemperature_Critical() << config->GetTemperature_Ref() << Unit.str() << auxFluidModel.GetTemperature_Critical() /config->GetTemperature_Ref();
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
      NonDimTable << "Conductivity" << "-" << config->GetThermal_Conductivity_Ref() << Unit.str() << "-";
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

void CFEM_DG_EulerSolver::SetUpTaskList(CConfig *config) {

  /* Check whether an ADER space-time step must be carried out.
     When only a spatial Jacobian is computed this is false per definition.  */
  if( (config->GetKind_TimeIntScheme_Flow() == ADER_DG) &&
     !(config->GetJacobian_Spatial_Discretization_Only()) ) {

    /*------------------------------------------------------------------------*/
    /* ADER time integration with local time stepping. There are 2^(M-1)      */
    /* subtime steps, where M is the number of time levels. For each of       */
    /* the subtime steps a number of tasks must be carried out, hence the     */
    /* total list can become rather lengthy.                                  */
    /*------------------------------------------------------------------------*/

    /*--- Determine whether or not the predictor solution must be interpolated
          for the time levels on this rank. Make a distinction between owned
          and halo elements. Also determine whether or not boundary conditions
          are present and whether or not these depend on halo data. ---*/
    const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();
    vector<bool> interpolOwnedElem(nTimeLevels, false);
    vector<bool> interpolHaloElem(nTimeLevels, false);
    vector<bool> BCPresent(nTimeLevels, false);
    vector<bool> BCDependOnHalos(nTimeLevels, false);

    for(unsigned short level=0; level<nTimeLevels; ++level) {

      /* First check the boundary conditions. */
      for(unsigned short iMarker=0; iMarker<config->GetnMarker_All(); iMarker++) {

        const unsigned long surfElemBeg = boundaries[iMarker].nSurfElem[level];
        const unsigned long surfElemEnd = boundaries[iMarker].nSurfElem[level+1];
        if(surfElemEnd > surfElemBeg) {
          BCPresent[level] = true;
          if( boundaries[iMarker].haloInfoNeededForBC ) BCDependOnHalos[level] = true;
        }
      }

      /* Determine whether or not the owned elements must be interpolated. */
      if(nVolElemOwnedPerTimeLevel[level+1] > nVolElemOwnedPerTimeLevel[level])
        interpolOwnedElem[level] = true;
      if(nMatchingInternalFacesLocalElem[level+1] > nMatchingInternalFacesLocalElem[level])
       interpolOwnedElem[level] = true;
      if(nMatchingInternalFacesWithHaloElem[level+1] > nMatchingInternalFacesWithHaloElem[level])
        interpolOwnedElem[level] = true;
      if( BCPresent[level] )
        interpolOwnedElem[level] = true;

      /* Determine whether or not the halo elements must be interpolated. */
      if(nMatchingInternalFacesWithHaloElem[level+1] > nMatchingInternalFacesWithHaloElem[level])
        interpolHaloElem[level] = true;
      if( BCDependOnHalos[level] )
        interpolHaloElem[level] = true;
    }

    /* Determine the number of subtime steps and abbreviate the number of
       time integration points a bit easier.  */
    const unsigned int   nSubTimeSteps          = pow(2,nTimeLevels-1);
    const unsigned short nTimeIntegrationPoints = config->GetnTimeIntegrationADER_DG();

    /* Define the two dimensional vector to store the latest index for a
       certain task for every time level. */
    vector<vector<int> > indexInList(CTaskDefinition::ADER_UPDATE_SOLUTION+1,
                                     vector<int>(nTimeLevels, -1));

    /* Loop over the number of subtime steps in the algorithm. */
    for(unsigned int step=0; step<nSubTimeSteps; ++step) {

      /* Determine the time level for which an update must be carried out
         for this step. */
      unsigned int ii = step+1;
      unsigned short timeLevel = 0;
      while( !(ii%2) ) {ii/=2; ++timeLevel;}

      /* Definition of the variable to store previous indices of
         tasks that must have been completed. */
      int prevInd[5];

      /* Carry out the predictor step of the communication elements of level 0
         if these elements are present on this rank. */
      unsigned long elemBeg = nVolElemOwnedPerTimeLevel[0]
                            + nVolElemInternalPerTimeLevel[0];
      unsigned long elemEnd = nVolElemOwnedPerTimeLevel[1];
      if(elemEnd > elemBeg) {
        prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][0];
        indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][0] = (int)tasksList.size();
        tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS, 0, prevInd[0]));
      }

      /* Initiate the communication of elements of level 0, if there is
         something to be communicated. */
#ifdef HAVE_MPI
      if( commRequests[0].size() ) {
        if(elemEnd > elemBeg)   // Data needs to be computed before sending.
          prevInd[0] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][0];
        else                    // Data only needs to be received.
          prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][0];

        /* Make sure that any previous communication has been completed. */
        prevInd[1] = indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][0];

        /* Create the task. */
        indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][0] = tasksList.size();
        tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_MPI_COMMUNICATION, 0,
                                            prevInd[0], prevInd[1]));
      }
#endif

      /* Check if the time level is not nTimeLevels-1. In that case carry out
         the predictor step of the communication elements for the next time
         level and initiate their communication, if appropriate on this rank. */
      if(timeLevel < (nTimeLevels-1)) {
        const unsigned short nL = timeLevel+1;

        /* Carry out the predictor step of the communication elements of level nL
           if these elements are present on this rank. */
        elemBeg = nVolElemOwnedPerTimeLevel[nL] + nVolElemInternalPerTimeLevel[nL];
        elemEnd = nVolElemOwnedPerTimeLevel[nL+1];

        if(elemEnd > elemBeg) {
          prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][nL];
          indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][nL] = (int)tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS, nL, prevInd[0]));
        }

        /* Initiate the communication of elements of level nL, if there is
           something to be communicated. */
#ifdef HAVE_MPI
        if( commRequests[nL].size() ) {
          if(elemEnd > elemBeg)   // Data needs to be computed before sending.
            prevInd[0] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][nL];
          else                    // Data only needs to be received.
            prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][nL];

          /* Make sure that any previous communication has been completed. */
          prevInd[1] = indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][nL];

          /* Create the actual task. */
          indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][nL] = tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_MPI_COMMUNICATION, nL,
                                              prevInd[0], prevInd[1]));
        }
#endif
      }

      /* Carry out the predictor step of the internal elements of time level 0,
         if these elements are present on this rank. */
      elemBeg = nVolElemOwnedPerTimeLevel[0];
      elemEnd = nVolElemOwnedPerTimeLevel[0] + nVolElemInternalPerTimeLevel[0];

      if(elemEnd > elemBeg) {
        prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][0];
        indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][0] = (int)tasksList.size();
        tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS, 0, prevInd[0]));
      }

      /* Determine the tasks to be completed before the communication of time
         level 0 can be completed. */
      prevInd[0] = prevInd[1] = prevInd[2] = -1;
#ifdef HAVE_MPI
      if( commRequests[0].size() )
        prevInd[0] = indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][0];
#endif

      if( elementsSendSelfComm[0].size() ) {
        prevInd[1] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][0];
        prevInd[2] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][0];
      }

      /* Make sure that the -1 in prevInd, if any, are numbered last. */
      sort(prevInd, prevInd+3, greater<int>());

      /* Complete the communication of time level 0, if there is something
         to be completed. */
      if(prevInd[0] > -1) {
        indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][0] = (int)tasksList.size();
        tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_MPI_COMMUNICATION, 0,
                                            prevInd[0], prevInd[1], prevInd[2]));
      }

      /* Check if the time level is not nTimeLevels-1. In that case carry out
         the predictor step of the internal elements for the next time
         level and complete the communication for this time level,
         if appropriate on this rank. */
      if(timeLevel < (nTimeLevels-1)) {
        const unsigned short nL = timeLevel+1;

        /* Carry out the predictor step of the internal elements of level nL
           if these elements are present on this rank. */
        elemBeg = nVolElemOwnedPerTimeLevel[nL];
        elemEnd = nVolElemOwnedPerTimeLevel[nL] + nVolElemInternalPerTimeLevel[nL];

        if(elemEnd > elemBeg) {
          prevInd[0] = indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][nL];
          indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][nL] = (int)tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS, nL, prevInd[0]));
        }

        /* Determine the tasks to be completed before the communication of time
           level nL can be completed. */
        prevInd[0] = prevInd[1] = prevInd[2] = -1;
#ifdef HAVE_MPI
        if( commRequests[nL].size() )
          prevInd[0] = indexInList[CTaskDefinition::INITIATE_MPI_COMMUNICATION][nL];
#endif

        if( elementsSendSelfComm[nL].size() ) {
          prevInd[1] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][nL];
          prevInd[2] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][nL];
        }

        /* Make sure that the -1 in prevInd, if any, are numbered last. */
        sort(prevInd, prevInd+3, greater<int>());

        /* Complete the communication of time level nL, if there is something
           to be completed. */
        if(prevInd[0] > -1) {
          indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][nL] = (int)tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_MPI_COMMUNICATION, nL,
                                              prevInd[0], prevInd[1], prevInd[2]));
        }
      }

      /* Loop over the time integration points to compute the space time
         integral for the corrector step. */
      for(unsigned short intPoint=0; intPoint<nTimeIntegrationPoints; ++intPoint) {

        /* Loop over the time levels to be treated. */
        for(unsigned short level=0; level<=timeLevel; ++level) {

          /* Easier storage of the number of owned elements for this level
             and the number of owned elements of the adjacent time level
             that are needed for the computation of the surface integral. */
          unsigned long nOwnedElem = nVolElemOwnedPerTimeLevel[level+1]
                                   - nVolElemOwnedPerTimeLevel[level];

          unsigned long nAdjOwnedElem = 0;
          if(level < (nTimeLevels-1))
            nAdjOwnedElem = ownedElemAdjLowTimeLevel[level+1].size();

          /* Check if the solution of the owned elements of the current time
             level and the adjacent owned elements of the next time level must
             be interpolated to the integration point intPoint. */
          if( interpolOwnedElem[level] ) {

            /* Determine the tasks that should have been completed before this
               task can be carried out. This depends on the time integration
               point. */
            if(intPoint == 0) {

              /* First time integration point. The predictor steps of the
                 owned elements must have been completed before this task
                 can be carried out. */
              if( nOwnedElem ) {
                prevInd[0] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][level];
                prevInd[1] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][level];
              }
              else {
                prevInd[0] = prevInd[1] = -1;
              }

              if( nAdjOwnedElem ) {
                prevInd[2] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS][level+1];
                prevInd[3] = indexInList[CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS][level+1];
              }
              else {
                prevInd[2] = prevInd[3] = -1;
              }

              prevInd[4] = -1;
            }
            else {

              /* Not the first integration point. The residual computations of
                 the previous integration point must have been completed, because
                 the solution in the work vectors is overwritten. */
              prevInd[0] = indexInList[CTaskDefinition::VOLUME_RESIDUAL][level];
              prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED][level];
              prevInd[2] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO][level];
              prevInd[3] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];
              prevInd[4] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS][level];
            }

            /* Create the task. */
            indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS,
                                                level, prevInd[0], prevInd[1], prevInd[2], prevInd[3], prevInd[4]));

            /* The info on the integration point and whether or not this
               time integration corresponds to the second part for the
               adjacent elements must be added to the task just created. */
            tasksList.back().intPointADER          = intPoint;
            tasksList.back().secondPartTimeIntADER = level < timeLevel;

            /* If artificial viscosity is used for the shock capturing
               terms, these terms must be computed for the owned elements,
               including the adjacent ones. */
            prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS][level];
            indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS,
                                                level, prevInd[0]));
          }

          /* The solution of the halo elements of the current time level and
             the adjacent halo elements of the next time level must be
             interpolated to the integration point intPoint. Check if these
             elements are present. */
          if( interpolHaloElem[level] ) {

            unsigned long nAdjHaloElem = 0;
            if(level < (nTimeLevels-1))
              nAdjHaloElem = haloElemAdjLowTimeLevel[level+1].size();

            unsigned long nHaloElem = nVolElemHaloPerTimeLevel[level+1]
                                    - nVolElemHaloPerTimeLevel[level];

            /* Determine the tasks that should have been completed before this
               task can be carried out. This depends on the time integration
               point. */
            if(intPoint == 0) {

              /* First time integration point. The communication of the
                 halo elements must have been completed before this task
                 can be carried out. */
             if( nHaloElem )
               prevInd[0] = indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][level];
             else
               prevInd[0] = -1;

             if( nAdjHaloElem )
               prevInd[1] = indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][level+1];
             else
               prevInd[1] = -1;
            }
            else {

              /* Not the first integration point. The residual computation of
                 the previous integration point must have been completed, because
                 the solution in the work vectors is overwritten. */
              prevInd[0] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];
              prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO][level];
            }

            /* Create the task. */
            indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS,
                                                level, prevInd[0], prevInd[1]));

            /* The info on the integration point and whether or not this
               time integration corresponds to the second part for the
               adjacent elements must be added to the task just created. */
            tasksList.back().intPointADER          = intPoint;
            tasksList.back().secondPartTimeIntADER = level < timeLevel;

            /* If artificial viscosity is used for the shock capturing
               terms, these terms must be computed for the halo elements,
               including the adjacent ones. */
            prevInd[0] = indexInList[CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS][level];
            indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS,
                                                level, prevInd[0]));
          }

          /* Check whether there are boundary conditions that involve halo elements. */
          if( BCDependOnHalos[level] ) {

            /* Create the dependency list for this task. */
            prevInd[0] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level];
            prevInd[1] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS][level];

            /* Create the task for the boundary conditions that involve halo elements. */
            indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO,
                                                level, prevInd[0], prevInd[1]));
          }

          /* Compute the surface residuals for this time level that involve
             halo elements, if present. Afterwards, accumulate the space time
             residuals for the halo elements. */
          if(nMatchingInternalFacesWithHaloElem[level+1] >
             nMatchingInternalFacesWithHaloElem[level]) {

            /* Create the dependencies for the surface residual part that involve
               halo elements. For all but the first integration point, make sure
               that the previous residual is already accumulated, because this
               task will overwrite that residual. */
            prevInd[0] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level];
            prevInd[1] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS][level];

            if(intPoint == 0)
              prevInd[2] = -1;
            else
              prevInd[2] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level];

            /* Create the task for the surface residual. */
            indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS,
                                                level, prevInd[0], prevInd[1], prevInd[2]));

            /* Create the task to accumulate the surface residuals of the halo
               elements. Make sure to set the integration point for this task. */
            prevInd[0] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];
            indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS,
                                                level, prevInd[0]));
            tasksList.back().intPointADER = intPoint;
          }

          /* If this is the last integration point, initiate the reverse
             communication for the residuals of this time level, if there
             is something to communicate. */
          if(intPoint == (nTimeIntegrationPoints-1)) {

            /* Check if there is something to communicate. Despite the fact
               that self communication takes place when the reverse communication
               is completed, a check for self communication is still needed,
               because of the periodic transformations. */
            bool commData = false;

#ifdef HAVE_MPI
            if( commRequests[level].size() ) commData = true;
#endif
            if(commData || elementsSendSelfComm[level].size()) {

              /* Determine the dependencies before the reverse communication
                 can start. */
              prevInd[0] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level];
              if( haloElemAdjLowTimeLevel[level].size() )
                prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS][level-1];
              else
                prevInd[1] = -1;

              prevInd[2] = indexInList[CTaskDefinition::COMPLETE_MPI_COMMUNICATION][level];

              /* Create the task. */
              indexInList[CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION][level] = (int)tasksList.size();
              tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION,
                                                  level, prevInd[0], prevInd[1], prevInd[2]));
            }
          }

          /* Compute the contribution of the volume integral, boundary conditions that only
             involve owned elements and surface integral between owned elements for this
             time level, if present. */
          if( nOwnedElem ) {

            /* Create the dependencies for this task. The shock capturing viscosity of
               the owned elements must be completed and for all integration points but
               the first, the computed residuals must have been accumulated in the space
               time residual, because the tasks below will overwrite the residuals. */
            prevInd[0] = indexInList[CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS][level];
            if(intPoint == 0)
              prevInd[1] = -1;
            else
              prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level];

            /* Create the tasks. */
            indexInList[CTaskDefinition::VOLUME_RESIDUAL][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::VOLUME_RESIDUAL, level,
                                                prevInd[0], prevInd[1]));

            indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED, level,
                                                prevInd[0], prevInd[1]));

            if(nMatchingInternalFacesLocalElem[level+1] >
               nMatchingInternalFacesLocalElem[level]) {
              indexInList[CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS][level] = (int)tasksList.size();
              tasksList.push_back(CTaskDefinition(CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS,
                                                  level, prevInd[0], prevInd[1]));
            }
          }

          /* Accumulate the space time residuals for the owned elements of this
             level, if these elements are present. */
          if(nAdjOwnedElem || nOwnedElem) {

            /* Create the dependencies for this task. */
            prevInd[0] = indexInList[CTaskDefinition::VOLUME_RESIDUAL][level];
            prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED][level];
            prevInd[1] = indexInList[CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO][level];
            prevInd[3] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS][level];
            prevInd[4] = indexInList[CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS][level];

            /* Create the task. */
            indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level] = (int)tasksList.size();
            tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS,
                                                level, prevInd[0], prevInd[1], prevInd[2], prevInd[3], prevInd[4]));
            tasksList.back().intPointADER = intPoint;
          }

          /* If this is the last integration point, complete the reverse
             communication for the residuals of this time level, if there
             is something to communicate. */
          if(intPoint == (nTimeIntegrationPoints-1)) {

            bool commData = false;

#ifdef HAVE_MPI
            if( commRequests[level].size() ) commData = true;
#endif
            if(commData || elementsSendSelfComm[level].size()) {

              /* Determine the dependencies before the reverse communication
                 can be completed. */
              prevInd[0] = indexInList[CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION][level];
              prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level];
              if( ownedElemAdjLowTimeLevel[level].size() )
                prevInd[2] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level-1];
              else
                prevInd[2] = -1;

              /* Create the task. */
              indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][level] = (int)tasksList.size();
              tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION,
                                                  level, prevInd[0], prevInd[1], prevInd[2]));
            }
          }

        } /* End loop time level. */

      } /* End loop time integration points. */

      /* Loop again over the number of active time levels for this subtime
         step and compute the update for the DOFs of the owned elements. */
      for(unsigned short level=0; level<=timeLevel; ++level) {
        if(nVolElemOwnedPerTimeLevel[level+1] > nVolElemOwnedPerTimeLevel[level]) {

          /* Updates must be carried out for this time level. First multiply the
             residuals by the inverse of the mass matrix. This task can be
             completed if the communication of this level has been completed.
             However, it is possible that no communication is needed and
             therefore the accumulation of the space time residuals is also
             added to the dependency list. */
          prevInd[0] = indexInList[CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION][level];
          prevInd[1] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level];
          if( ownedElemAdjLowTimeLevel[level].size() )
            prevInd[2] = indexInList[CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS][level-1];
          else
            prevInd[2] = -1;

          indexInList[CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX][level] = (int)tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX,
                                               level, prevInd[0], prevInd[1], prevInd[2]));

          /* Compute the new state vector for this time level. */
          prevInd[0] = indexInList[CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX][level];
          indexInList[CTaskDefinition::ADER_UPDATE_SOLUTION][level] = (int)tasksList.size();
          tasksList.push_back(CTaskDefinition(CTaskDefinition::ADER_UPDATE_SOLUTION,
                                               level, prevInd[0]));
        }
      }

    } /* End loop subtime steps. */


    // EXTRA FOR DEBUGGING
/*  SU2_OMP_SINGLE
    {
      if(rank == MASTER_NODE) {
        cout << endl;
        cout << "Task list for rank " << rank << endl;
        cout << "------------------------------------------------" << endl;
        cout << "Number of tasks: " << tasksList.size() << endl;
        for(unsigned long j=0; j<tasksList.size(); ++j) {

          cout << "Task " << j << ": ";
          switch( tasksList[j].task ) {
            case CTaskDefinition::NO_TASK: cout << "NO_TASK" << endl; break;
            case CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS: cout << "ADER_PREDICTOR_STEP_COMM_ELEMENTS" << endl; break;
            case CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS: cout << "ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS" << endl; break;
            case CTaskDefinition::INITIATE_MPI_COMMUNICATION: cout << "INITIATE_MPI_COMMUNICATION" << endl; break;
            case CTaskDefinition::COMPLETE_MPI_COMMUNICATION: cout << "COMPLETE_MPI_COMMUNICATION" << endl; break;
            case CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION: cout << "INITIATE_REVERSE_MPI_COMMUNICATION" << endl; break;
            case CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION: cout << "COMPLETE_REVERSE_MPI_COMMUNICATION" << endl; break;
            case CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS: cout << "ADER_TIME_INTERPOLATE_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS: cout << "ADER_TIME_INTERPOLATE_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS: cout << "SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS: cout << "SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::VOLUME_RESIDUAL: cout << "VOLUME_RESIDUAL" << endl; break;
            case CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS: cout << "SURFACE_RESIDUAL_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS: cout << "SURFACE_RESIDUAL_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED: cout << "BOUNDARY_CONDITIONS_DEPEND_ON_OWNED" << endl; break;
            case CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO: cout << "BOUNDARY_CONDITIONS_DEPEND_ON_HALO" << endl; break;
            case CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS: cout << "SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS: cout << "SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS: cout << "ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS" << endl; break;
            case CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS: cout << "ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS" << endl; break;
            case CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX: cout << "MULTIPLY_INVERSE_MASS_MATRIX" << endl; break;
            case CTaskDefinition::ADER_UPDATE_SOLUTION: cout << "ADER_UPDATE_SOLUTION" << endl; break;
            default: cout << "This cannot happen" << endl;
          }
          cout << " Time level: " << tasksList[j].timeLevel
               << " Integration point: " << tasksList[j].intPointADER
               << " Second part: " << tasksList[j].secondPartTimeIntADER << endl;
          cout << " Depends on tasks:";
          for(unsigned short k=0; k<tasksList[j].nIndMustBeCompleted; ++k)
            cout << " " << tasksList[j].indMustBeCompleted[k];
          cout << endl << endl;
        }

        cout << "CFEM_DG_EulerSolver::SetUpTaskList: ADER tasklist printed" << endl;
      }
#ifdef HAVE_MPI
      SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif
    }
    END_SU2_OMP_SINGLE

    exit(1); */

    // END EXTRA FOR DEBUGGING.
  }
  else {

    /*------------------------------------------------------------------------*/
    /* Standard time integration scheme for which the spatial residual must   */
    /* be computed for the DOFS of the owned elements. This results in a      */
    /* relatively short tasks list, which can be set easily.                  */
    /*------------------------------------------------------------------------*/

    tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_MPI_COMMUNICATION,                   0));                   // Task  0
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS,     0));                   // Task  1
    tasksList.push_back(CTaskDefinition(CTaskDefinition::VOLUME_RESIDUAL,                              0,  1));               // Task  2
    tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_MPI_COMMUNICATION,                   0,  0));               // Task  3
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS,      0,  3));               // Task  4
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS,               0,  1, 4));            // Task  5
    tasksList.push_back(CTaskDefinition(CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO,           0,  1, 3));            // Task  6
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS,  0,  5));               // Task  7
    tasksList.push_back(CTaskDefinition(CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION,           0,  7));               // Task  8
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS,              0,  1));               // Task  9
    tasksList.push_back(CTaskDefinition(CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED,          0,  1));               // Task 10
    tasksList.push_back(CTaskDefinition(CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION,           0,  2, 8));            // Task 11
    tasksList.push_back(CTaskDefinition(CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS, 0,  2, 6, 9, 10, 11)); // Task 12
    tasksList.push_back(CTaskDefinition(CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX,                 0, 12));               // Task 13
  }
}

void CFEM_DG_EulerSolver::ConservativeToEntropyVariables(ColMajorMatrix<su2double> &sol) {

  /*--- Easier storage of the number of items to be converted and
        the inverse of gamma-1. ---*/
  const unsigned short nItems = sol.rows();
  const su2double ovgm1       = 1.0/Gamma_Minus_One;

  /*--- Make a distinction between two and three space dimensions
        in order to have the most efficient code. ---*/
  switch( nDim ) {

    case 2: {

      /*--- Two dimensional simulation. Loop over the number of entities
            and carry out the conversion. ---*/
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nItems; ++i) {
        const su2double rho    = sol(i,0);
        const su2double rhoInv = 1.0/rho;
        const su2double u      = rhoInv*sol(i,1);
        const su2double v      = rhoInv*sol(i,2);
        const su2double kin    = 0.5*(u*u + v*v);
        const su2double p      = Gamma_Minus_One*(sol(i,3) - rho*kin);
        const su2double pInv   = 1.0/p;
        const su2double s      = log(p*pow(rhoInv,Gamma));
        const su2double abv    = pInv*rho;

        sol(i,0) =  (Gamma-s)*ovgm1 - abv*kin;
        sol(i,1) =  abv*u;
        sol(i,2) =  abv*v;
        sol(i,3) = -abv;
      }

      break;
    }

    case 3: {

      /*--- Three dimensional simulation. Loop over the number of entities
            and carry out the conversion. ---*/
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nItems; ++i) {
        const su2double rho    = sol(i,0);
        const su2double rhoInv = 1.0/rho;
        const su2double u      = rhoInv*sol(i,1);
        const su2double v      = rhoInv*sol(i,2);
        const su2double w      = rhoInv*sol(i,3);
        const su2double kin    = 0.5*(u*u + v*v + w*w);
        const su2double p      = Gamma_Minus_One*(sol(i,4) - rho*kin);
        const su2double pInv   = 1.0/p;
        const su2double s      = log(p*pow(rhoInv,Gamma));
        const su2double abv    = pInv*rho;

        sol(i,0) =  (Gamma-s)*ovgm1 - pInv*rho*kin;
        sol(i,1) =  abv*u;
        sol(i,2) =  abv*v;
        sol(i,3) =  abv*w;
        sol(i,4) = -abv;
      }

      break;
    }
  }
}

void CFEM_DG_EulerSolver::EntropyToPrimitiveVariables(ColMajorMatrix<su2double> &sol) {

  /*--- Easier storage of the number of items to be converted and
        the inverse of 1-gamma. ---*/
  const unsigned short nItems =  sol.rows();
  const su2double gm1         =  Gamma_Minus_One;
  const su2double ov1mg       = -1.0/gm1;

  /*--- Make a distinction between two and three space dimensions
        in order to have the most efficient code. ---*/
  switch( nDim ) {

    case 2: {

      /*--- Two dimensional simulation. Loop over the number of entities
            and carry out the conversion. ---*/
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nItems; ++i) {
        const su2double V3Inv =  1.0/sol(i,3);
        const su2double u     = -V3Inv*sol(i,1);
        const su2double v     = -V3Inv*sol(i,2);
        const su2double eKin  =  0.5*(u*u + v*v);
        const su2double s     =  Gamma - gm1*(sol(i,0) - sol(i,3)*eKin);
        const su2double tmp   = -sol(i,3)*exp(s);
        const su2double rho   =  pow(tmp, ov1mg);
        const su2double p     = -rho*V3Inv;

        sol(i,0) = rho;
        sol(i,1) = u;
        sol(i,2) = v;
        sol(i,3) = p;
      }

      break;
    }

    case 3: {

      /*--- Three dimensional simulation. Loop over the number of entities
            and carry out the conversion. ---*/
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nItems; ++i) {
        const su2double V4Inv =  1.0/sol(i,4);
        const su2double u     = -V4Inv*sol(i,1);
        const su2double v     = -V4Inv*sol(i,2);
        const su2double w     = -V4Inv*sol(i,3);
        const su2double eKin  =  0.5*(u*u + v*v + w*w);
        const su2double s     =  Gamma - gm1*(sol(i,0) - sol(i,4)*eKin);
        const su2double tmp   = -sol(i,4)*exp(s);
        const su2double rho   =  pow(tmp, ov1mg);
        const su2double p     = -rho*V4Inv;

        sol(i,0) = rho;
        sol(i,1) = u;
        sol(i,2) = v;
        sol(i,3) = w;
        sol(i,4) = p;
      }

      break;
    }
  }
}

void CFEM_DG_EulerSolver::Initiate_MPI_Communication(CConfig *config,
                                                     const unsigned short timeLevel) {
#ifdef HAVE_MPI

  /*--- Check if there is anything to communicate. ---*/
  if( commRequests[timeLevel].size() ) {

    /*--- Easier storage whether or not an ADER simulation is performed
          and the number of time DOFs to be communicated. ---*/
    const bool ADER = (config->GetKind_TimeIntScheme_Flow() == ADER_DG);
    const unsigned short nTimeDOFs = config->GetnTimeDOFsADER_DG();

    /*--- OpenMP loop over the number of ranks to which
          data must be sent. ---*/
    SU2_OMP_FOR_(schedule(static,1) SU2_NOWAIT)
    for(unsigned long i=0; i<ranksSendMPI[timeLevel].size(); ++i) {
      unsigned long ii = 0;

      /*--- Loop over the elements to be sent and copy the solution data
            into the send buffer. ---*/
      su2double *sendBuf = commSendBuf[timeLevel][i].data();
      for(unsigned long j=0; j<elementsSendMPIComm[timeLevel][i].size(); ++j) {
        const unsigned long jj = elementsSendMPIComm[timeLevel][i][j];
        const unsigned long nItems = volElem[jj].nDOFsSol;

        /*--- Make a distinction between ADER and a standard time
              integration scheme. ---*/
        if( ADER ) {
          for(unsigned short k=0; k<nTimeDOFs; ++k) {
            for(unsigned short iVar=0; iVar<nVar; ++iVar) {
              SU2_OMP_SIMD
              for(unsigned long mm=0; mm<nItems; ++mm)
                sendBuf[mm+ii] = volElem[jj].solDOFsADERPredictor[k](mm,iVar);
              ii += nItems;
            }
          }
        }
        else {
          for(unsigned short iVar=0; iVar<nVar; ++iVar) {
            SU2_OMP_SIMD
            for(unsigned long mm=0; mm<nItems; ++mm)
              sendBuf[mm+ii] = volElem[jj].solDOFsWork(mm,iVar);
            ii += nItems;
          }
        }
      }

      /*--- Send the data using non-blocking sends. ---*/
      int dest = ranksSendMPI[timeLevel][i];
      int tag  = dest + timeLevel;
      SU2_MPI::Isend(sendBuf, ii, MPI_DOUBLE, dest, tag, SU2_MPI::GetComm(),
                     &commRequests[timeLevel][i]);
    }
    END_SU2_OMP_FOR

    /*--- OpenMP loop over the number of ranks from which data is received. ---*/
    SU2_OMP_FOR_STAT(1)
    for(unsigned long i=0; i<ranksRecvMPI[timeLevel].size(); ++i) {

      /*--- Post the non-blocking receive. ---*/
      unsigned long indComm = i + ranksSendMPI[timeLevel].size();
      int source = ranksRecvMPI[timeLevel][i];
      int tag    = rank + timeLevel;
      SU2_MPI::Irecv(commRecvBuf[timeLevel][i].data(),
                     commRecvBuf[timeLevel][i].size(),
                     MPI_DOUBLE, source, tag, SU2_MPI::GetComm(),
                     &commRequests[timeLevel][indComm]);
    }
    END_SU2_OMP_FOR
  }
#endif
}

bool CFEM_DG_EulerSolver::Complete_MPI_Communication(CConfig *config,
                                                     const unsigned short timeLevel,
                                                     const bool commMustBeCompleted) {

  /*--- Easier storage whether or not an ADER simulation is performed
        and the number of time DOFs to be communicated. ---*/
  const bool ADER = (config->GetKind_TimeIntScheme_Flow() == ADER_DG);
  const unsigned short nTimeDOFs = config->GetnTimeDOFsADER_DG();

  /*-----------------------------------------------------------------------*/
  /*--- Complete the MPI communication, if needed and if possible, and  ---*/
  /*--- copy the data from the receive buffers into the correct         ---*/
  /*--- location in volElem.                                            ---*/
  /*-----------------------------------------------------------------------*/

#ifdef HAVE_MPI
  if( commRequests[timeLevel].size() ) {

    /*--- There are communication requests to be completed. Check if these
          requests must be completed. In that case Waitall is used.
          Otherwise, Testall is used to check if all the requests have
          been completed. If not, false is returned. Use the member variable
          counter as a flag for all threads. ---*/
    SU2_OMP_SINGLE
    {
      counter = 1;
      if( commMustBeCompleted ) {
        SU2_MPI::Waitall(commRequests[timeLevel].size(),
                         commRequests[timeLevel].data(), MPI_STATUSES_IGNORE);
      }
      else {
        int flag;
        SU2_MPI::Testall(commRequests[timeLevel].size(),
                         commRequests[timeLevel].data(), &flag, MPI_STATUSES_IGNORE);
        if( !flag ) counter = 0;
      }
    }
    END_SU2_OMP_SINGLE

    if( !counter ) return false;

    /*--- OpenMP loop over the number of ranks from which
          this rank has received data. ---*/
    SU2_OMP_FOR_(schedule(static,1) SU2_NOWAIT)
    for(unsigned long i=0; i<ranksRecvMPI[timeLevel].size(); ++i) {
      unsigned long ii = 0;

      /*--- Loop over the elements to be received and copy the solution data
            into the appropriate locations in volElem. ---*/
      su2double *recvBuf = commRecvBuf[timeLevel][i].data();
      for(unsigned long j=0; j<elementsRecvMPIComm[timeLevel][i].size(); ++j) {
        const unsigned long jj = elementsRecvMPIComm[timeLevel][i][j];
        const unsigned long nItems = volElem[jj].nDOFsSol;

        /*--- Make a distinction between ADER and a standard time
              integration scheme. ---*/
        if( ADER ) {
          for(unsigned short k=0; k<nTimeDOFs; ++k) {
            for(unsigned short iVar=0; iVar<nVar; ++iVar) {
              SU2_OMP_SIMD
              for(unsigned long mm=0; mm<nItems; ++mm)
                volElem[jj].solDOFsADERPredictor[k](mm,iVar) = recvBuf[mm+ii];
              ii += nItems;
            }
          }
        }
        else {
          for(unsigned short iVar=0; iVar<nVar; ++iVar) {
            SU2_OMP_SIMD
            for(unsigned long mm=0; mm<nItems; ++mm)
              volElem[jj].solDOFsWork(mm,iVar) = recvBuf[mm+ii];
            ii += nItems;
          }
        }
      }
    }
    END_SU2_OMP_FOR
  }

#endif

  /*-----------------------------------------------------------------------*/
  /*---               Carry out the self communication.                 ---*/
  /*-----------------------------------------------------------------------*/

  /*--- Determine the chunk size for the OMP loops, if supported. ---*/
  const unsigned long nSelfCom = elementsSendSelfComm[timeLevel].size();
#ifdef HAVE_OMP
  const size_t omp_chunk_size_elem = computeStaticChunkSize(nSelfCom, omp_get_num_threads(), 64);
#endif

  /*--- Loop over the elements for self-communication. ---*/
  SU2_OMP_FOR_STAT(omp_chunk_size_elem)
  for(unsigned long i=0; i<nSelfCom; ++i) {
    const unsigned long elemS  = elementsSendSelfComm[timeLevel][i];
    const unsigned long elemR  = elementsRecvSelfComm[timeLevel][i];
    const unsigned long nItems = volElem[elemS].nDOFsSol;

    /*--- Make a distinction between ADER and a standard time
          integration scheme. ---*/
    if( ADER ) {
      for(unsigned short k=0; k<nTimeDOFs; ++k) {
        for(unsigned short iVar=0; iVar<nVar; ++iVar) {
          SU2_OMP_SIMD
          for(unsigned long mm=0; mm<nItems; ++mm)
            volElem[elemR].solDOFsADERPredictor[k](mm,iVar) = volElem[elemS].solDOFsADERPredictor[k](mm,iVar);
        }
      }
    }
    else {
      for(unsigned short iVar=0; iVar<nVar; ++iVar) {
        SU2_OMP_SIMD
        for(unsigned long mm=0; mm<nItems; ++mm)
          volElem[elemR].solDOFsWork(mm,iVar) = volElem[elemS].solDOFsWork(mm,iVar);
      }
    }
  }
  END_SU2_OMP_FOR

  /*------------------------------------------------------------------------*/
  /*--- Correct the vector quantities in the rotational periodic halo's. ---*/
  /*------------------------------------------------------------------------*/

  /*--- Loop over the markers for which a rotational periodic
        correction must be applied to the momentum variables. ---*/
  unsigned int ii = 0;
  const unsigned long nRotPerMarkers = halosRotationalPeriodicity[timeLevel].size();
  for(unsigned long l=0; l<nRotPerMarkers; ++l) {

    /*--- Easier storage of the rotational matrix. ---*/
    su2double rotMatrix[3][3];

    rotMatrix[0][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[0][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[0][2] = rotationMatricesPeriodicity[ii++];

    rotMatrix[1][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[1][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[1][2] = rotationMatricesPeriodicity[ii++];

    rotMatrix[2][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[2][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[2][2] = rotationMatricesPeriodicity[ii++];

    /*--- Determine the number of elements for this transformation and
          loop over them. ---*/
    const unsigned long nHaloElem = halosRotationalPeriodicity[timeLevel][l].size();
#ifdef HAVE_OMP
    const size_t omp_chunk_size_halo = computeStaticChunkSize(nHaloElem, omp_get_num_threads(), 64);
#endif
    SU2_OMP_FOR_STAT(omp_chunk_size_halo)
    for(unsigned long j=0; j<nHaloElem; ++j) {
      const unsigned long ind = halosRotationalPeriodicity[timeLevel][l][j];
      const unsigned long nItems = volElem[ind].nDOFsSol;

      /*--- Make a distinction between ADER and a standard time
            integration scheme. ---*/
      if( ADER ) {
        for(unsigned short k=0; k<nTimeDOFs; ++k) {

          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned long mm=0; mm<nItems; ++mm) {
            const su2double u = volElem[ind].solDOFsADERPredictor[k](mm,1),
                            v = volElem[ind].solDOFsADERPredictor[k](mm,2),
                            w = volElem[ind].solDOFsADERPredictor[k](mm,3);

            volElem[ind].solDOFsADERPredictor[k](mm,1) = rotMatrix[0][0]*u + rotMatrix[0][1]*v + rotMatrix[0][2]*w;
            volElem[ind].solDOFsADERPredictor[k](mm,2) = rotMatrix[1][0]*u + rotMatrix[1][1]*v + rotMatrix[1][2]*w;
            volElem[ind].solDOFsADERPredictor[k](mm,3) = rotMatrix[2][0]*u + rotMatrix[2][1]*v + rotMatrix[2][2]*w;
          }
        }
      }
      else {
        SU2_OMP_SIMD_IF_NOT_AD
        for(unsigned long mm=0; mm<nItems; ++mm) {
          const su2double u = volElem[ind].solDOFsWork(mm,1),
                          v = volElem[ind].solDOFsWork(mm,2),
                          w = volElem[ind].solDOFsWork(mm,3);

          volElem[ind].solDOFsWork(mm,1) = rotMatrix[0][0]*u + rotMatrix[0][1]*v + rotMatrix[0][2]*w;
          volElem[ind].solDOFsWork(mm,2) = rotMatrix[1][0]*u + rotMatrix[1][1]*v + rotMatrix[1][2]*w;
          volElem[ind].solDOFsWork(mm,3) = rotMatrix[2][0]*u + rotMatrix[2][1]*v + rotMatrix[2][2]*w;
        }
      }
    }
    END_SU2_OMP_FOR
  }

  /*--- Return true to indicate that the communication has been completed. ---*/
  return true;
}

void CFEM_DG_EulerSolver::Initiate_MPI_ReverseCommunication(CConfig *config,
                                                            const unsigned short timeLevel) {

  /*--- Easier storage whether or not an ADER simulation is performed. ---*/
  const bool ADER = (config->GetKind_TimeIntScheme_Flow() == ADER_DG);

  /*------------------------------------------------------------------------*/
  /*--- Correct the vector residuals in the rotational periodic halo's.  ---*/
  /*------------------------------------------------------------------------*/

  /*--- Loop over the markers for which a rotational periodic
        correction must be applied to the momentum variables. ---*/
  unsigned int ii = 0;
  const unsigned short nRotPerMarkers = halosRotationalPeriodicity[0].size();
  for(unsigned short l=0; l<nRotPerMarkers; ++l) {

    /* Easier storage of the transpose of the rotational matrix. */
    su2double rotMatrix[3][3];

    rotMatrix[0][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[1][0] = rotationMatricesPeriodicity[ii++];
    rotMatrix[2][0] = rotationMatricesPeriodicity[ii++];

    rotMatrix[0][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[1][1] = rotationMatricesPeriodicity[ii++];
    rotMatrix[2][1] = rotationMatricesPeriodicity[ii++];

    rotMatrix[0][2] = rotationMatricesPeriodicity[ii++];
    rotMatrix[1][2] = rotationMatricesPeriodicity[ii++];
    rotMatrix[2][2] = rotationMatricesPeriodicity[ii++];

    /*--- Determine the number of elements for this transformation and
          loop over them. ---*/
    const unsigned long nHaloElem = halosRotationalPeriodicity[timeLevel][l].size();
#ifdef HAVE_OMP
    const size_t omp_chunk_size_halo = computeStaticChunkSize(nHaloElem, omp_get_num_threads(), 64);
#endif
    SU2_OMP_FOR_STAT(omp_chunk_size_halo)
    for(unsigned long j=0; j<nHaloElem; ++j) {
      const unsigned long ind = halosRotationalPeriodicity[timeLevel][l][j];
      const unsigned long nItems = volElem[ind].nDOFsSol;

      /*--- Set the reference for the residual, depending whether
            or not an ADER simulation is carried out. ---*/
      ColMajorMatrix<su2double> &res = ADER ? volElem[ind].resTotDOFsADER
                                            : volElem[ind].resDOFs;

      /*--- Loop over the items and correct the momentum residuals. ---*/
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned long mm=0; mm<nItems; ++mm) {
        const su2double ru = res(mm,1), rv = res(mm,2), rw = res(mm,3);

        res(mm,1) = rotMatrix[0][0]*ru + rotMatrix[0][1]*rv + rotMatrix[0][2]*rw;
        res(mm,2) = rotMatrix[1][0]*ru + rotMatrix[1][1]*rv + rotMatrix[1][2]*rw;
        res(mm,3) = rotMatrix[2][0]*ru + rotMatrix[2][1]*rv + rotMatrix[2][2]*rw;
      }
    }
    END_SU2_OMP_FOR
  }

#ifdef HAVE_MPI

  /*--- Check if there is anything to communicate. ---*/
  if( commRequests[timeLevel].size() ) {

    /*--- Loop over the number of ranks from which this rank receives data in
          the original communication pattern. In the reverse pattern, data
          has to be sent. ---*/
    SU2_OMP_FOR_(schedule(static,1) SU2_NOWAIT)
    for(unsigned long i=0; i<ranksRecvMPI[timeLevel].size(); ++i) {
      unsigned long ii = 0;

      /*--- Loop over the elements copy the residual data into recvBuf. ---*/
      su2double *recvBuf = commRecvBuf[timeLevel][i].data();
      for(unsigned long j=0; j<elementsRecvMPIComm[timeLevel][i].size(); ++j) {
        const unsigned long jj = elementsRecvMPIComm[timeLevel][i][j];
        const unsigned long nItems = volElem[jj].nDOFsSol;

        /*--- Set the reference to the correct residual to be communicated. ---*/
        ColMajorMatrix<su2double> &res = ADER ? volElem[jj].resTotDOFsADER
                                              : volElem[jj].resDOFs;

        /*--- Copy the data into recvBuf. ---*/
        for(unsigned short iVar=0; iVar<nVar; ++iVar) {
          SU2_OMP_SIMD
          for(unsigned long mm=0; mm<nItems; ++mm)
            recvBuf[mm+ii] = res(mm,iVar);
           ii += nItems;
        }
      }

      /*--- Send the data using non-blocking sends. ---*/
      int dest = ranksRecvMPI[timeLevel][i];
      int tag  = dest + timeLevel + 20;
      SU2_MPI::Isend(recvBuf, ii, MPI_DOUBLE, dest, tag, SU2_MPI::GetComm(),
                     &commRequests[timeLevel][i]);
    }
    END_SU2_OMP_FOR

    /*--- Post the non-blocking receives. As this is the reverse communication,
          a loop over the sending ranks must be carried out. ---*/
    SU2_OMP_FOR_STAT(1)
    for(unsigned long i=0; i<ranksSendMPI[timeLevel].size(); ++i) {

      unsigned long indComm = i + ranksRecvMPI[timeLevel].size();
      int source = ranksSendMPI[timeLevel][i];
      int tag    = rank + timeLevel + 20;
      SU2_MPI::Irecv(commSendBuf[timeLevel][i].data(),
                     commSendBuf[timeLevel][i].size(),
                     MPI_DOUBLE, source, tag, SU2_MPI::GetComm(),
                     &commRequests[timeLevel][indComm]);
    }
    END_SU2_OMP_FOR
  }
#endif
}

bool CFEM_DG_EulerSolver::Complete_MPI_ReverseCommunication(CConfig *config,
                                                            const unsigned short timeLevel,
                                                            const bool commMustBeCompleted) {

  /*--- Easier storage whether or not an ADER simulation is performed. ---*/
  const bool ADER = (config->GetKind_TimeIntScheme_Flow() == ADER_DG);

  /*-----------------------------------------------------------------------*/
  /*--- Complete the MPI communication, if needed and if possible, and  ---*/
  /*--- copy the data from the receive buffers into the correct         ---*/
  /*--- location in volElem.                                            ---*/
  /*-----------------------------------------------------------------------*/

#ifdef HAVE_MPI
  if( commRequests[timeLevel].size() ) {

    /*--- There are communication requests to be completed. Check if these
          requests must be completed. In that case Waitall is used.
          Otherwise, Testall is used to check if all the requests have
          been completed. If not, false is returned. Use the member variable
          counter as a flag for all threads. ---*/
    SU2_OMP_SINGLE
    {
      counter = 1;
      if( commMustBeCompleted ) {
        SU2_MPI::Waitall(commRequests[timeLevel].size(),
                         commRequests[timeLevel].data(), MPI_STATUSES_IGNORE);
      }
      else {
        int flag;
        SU2_MPI::Testall(commRequests[timeLevel].size(),
                         commRequests[timeLevel].data(), &flag, MPI_STATUSES_IGNORE);
        if( !flag ) counter = 0;
      }
    }
    END_SU2_OMP_SINGLE

    if( !counter ) return false;

    /*-------------------------------------------------------------------------*/
    /*---    Update the residuals of the owned DOFs with the data received. ---*/
    /*-------------------------------------------------------------------------*/

    /*--- Loop over the received residual data from all ranks and update the
          residual of the DOFs of the corresponding elements. Note that in
          reverse mode the send communication data must be used.
          To avoid a race condition in the update, the loop below must
          be carried out by a single thread. ---*/
    SU2_OMP_SINGLE
    {
      for(unsigned long i=0; i<ranksSendMPI[timeLevel].size(); ++i) {
        unsigned long ii = 0;

        /*--- Loop over the elements that must be updated. ---*/
        su2double *sendBuf = commSendBuf[timeLevel][i].data();
        for(unsigned long j=0; j<elementsSendMPIComm[timeLevel][i].size(); ++j) {
          const unsigned long jj = elementsSendMPIComm[timeLevel][i][j];
          const unsigned long nItems = volElem[jj].nDOFsSol;

          /*--- Set the reference for the residual, depending on the situation. ---*/
          ColMajorMatrix<su2double> &res = ADER ? volElem[jj].resTotDOFsADER
                                                : volElem[jj].resDOFs;

          /*--- Update the residuals of the DOFs. ---*/
          for(unsigned short iVar=0; iVar<nVar; ++iVar) {
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned long mm=0; mm<nItems; ++mm)
              res(mm,iVar) += sendBuf[mm+ii];
            ii += nItems;
          }
        }
      }
    }
    END_SU2_OMP_SINGLE
  }
#endif

  /*-----------------------------------------------------------------------*/
  /*---               Carry out the self communication.                 ---*/
  /*-----------------------------------------------------------------------*/

  /*--- Loop over the elements for self-communication. Also this loop must
        be carried out by a single thread to avoid a race condition. ---*/
  SU2_OMP_SINGLE
  for(unsigned long i=0; i<elementsSendSelfComm[timeLevel].size(); ++i) {
    const unsigned long elemS  = elementsSendSelfComm[timeLevel][i];
    const unsigned long elemR  = elementsRecvSelfComm[timeLevel][i];
    const unsigned long nItems = volElem[elemS].nDOFsSol;

    /*--- Set the references for the residual, depending on the situation. ---*/
    ColMajorMatrix<su2double> &resS = ADER ? volElem[elemS].resTotDOFsADER
                                           : volElem[elemS].resDOFs;
    ColMajorMatrix<su2double> &resR = ADER ? volElem[elemR].resTotDOFsADER
                                           : volElem[elemR].resDOFs;

    /*--- Update the send residual. ---*/
    for(unsigned short iVar=0; iVar<nVar; ++iVar) {
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned long j=0; j<nItems; ++j)
        resS(j,iVar) += resR(j,iVar);
    }
  }
  END_SU2_OMP_SINGLE

  /*-----------------------------------------------------------------------*/
  /*---   Initialize the halo residuals for this time level for ADER.   ---*/
  /*-----------------------------------------------------------------------*/

  if( ADER ) {

#ifdef HAVE_OMP
    const unsigned long nHaloElem = nVolElemHaloPerTimeLevel[timeLevel+1]
                                  - nVolElemHaloPerTimeLevel[timeLevel];
    const size_t omp_chunk_size_halo = computeStaticChunkSize(nHaloElem, omp_get_num_threads(), 64);
#endif
    SU2_OMP_FOR_STAT(omp_chunk_size_halo)
    for(unsigned long i=nVolElemHaloPerTimeLevel[timeLevel];
                      i<nVolElemHaloPerTimeLevel[timeLevel+1]; ++i)
      volElem[i].resTotDOFsADER.setConstant(0.0);
    END_SU2_OMP_FOR
  }

  /*--- Return true to indicate that the communication has been completed. ---*/
  return true;
}

void CFEM_DG_EulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  /*--- Check if a verification solution is to be computed. ---*/
  if ((VerificationSolution)  && (TimeIter == 0)) {

    /*--- Start OpenMP parallel region. ---*/
    SU2_OMP_PARALLEL {

      /*--- Loop over the owned elements. ---*/
#ifdef HAVE_OMP
      const size_t omp_chunk_size_vol = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif
      SU2_OMP_FOR_STAT(omp_chunk_size_vol)
      for(unsigned long i=0; i<nVolElemOwned; ++i) {

        /*--- Define the help variables to store the coordinates
              and the solution in the nodal DOFs. ---*/
        su2double coor[3] = {0.0}, solDOF[5] = {0.0};

        /*--- Loop over the DOFs of this element. ---*/
        const unsigned short nDOFs = volElem[i].standardElemFlow->GetNDOFs();
        for(unsigned short j=0; j<nDOFs; ++j) {

          /*--- Store the coordinates of the DOFs in coor. ---*/
          for(unsigned short iDim=0; iDim<nDim; ++iDim)
            coor[iDim] = volElem[i].coorSolDOFs(j,iDim);

          /*--- Compute the initial condition provided by the verification solution
                class. This can be the exact solution, but this is not necessary.
                Note that the conservative variables are computed. ---*/
          VerificationSolution->GetInitialCondition(coor, solDOF);

          /*--- Store the conservative variables in the DOFs for the element
                for the time being. ---*/
          for(unsigned short iVar=0; iVar<nVar; ++iVar)
            volElem[i].solDOFs(j,iVar) = solDOF[iVar];
        }

        /*--- Convert the conservative variables to entropy variables. ---*/
        ConservativeToEntropyVariables(volElem[i].solDOFs);

        /*--- Convert the nodal solution to the modal solution. ---*/
        volElem[i].NodalToModalFlow();
      }
      END_SU2_OMP_FOR
    }
    END_SU2_OMP_PARALLEL
  }
}

void CFEM_DG_EulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iStep, unsigned short RunTime_EqSystem, bool Output) {

  /*--- Initialize the norms for the residual monotoring to zero. ---*/
  SetResToZero();

/*--- Determine the chunk size for the OMP loops, if supported. ---*/
#ifdef HAVE_OMP
  const size_t omp_chunk_size_elem = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif

  /*--- Initialize the counter for the number of bad elements to zero. ---*/
  SU2_OMP_SINGLE
  counter = 0;
  END_SU2_OMP_SINGLE

  /*--- Define the local counter for the number of bad elements. ---*/
  unsigned long ErrorCounter = 0;

  /*--- Loop over the number of owned elements and check for non-physical
        solutions in the integration points. When found the solution in the
        element is overwritten with a constant free-stream solution. This
        does not work usually. ---*/
  SU2_OMP_FOR_STAT(omp_chunk_size_elem)
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /*--- Determine the solution in the integration points and
          convert it to primitive variables. ---*/
    ColMajorMatrix<su2double> &solInt = volElem[l].ComputeSolIntPoints(volElem[l].solDOFs);
    EntropyToPrimitiveVariables(solInt);

    /*--- Check for negative density and pressure in the integration
          points. If found the element is flag as bad. ---*/
    bool badElement = false;
    const unsigned short nInt = volElem[l].standardElemFlow->GetNIntegration();
    for(unsigned short i=0; i<nInt; ++i)
      if((solInt(i,0) <= 0.0) || (solInt(i,nDim+1) <= 0.0)) badElement = true;

    /*--- If the element is bad, reset the state in the element and
          update the counter ErrorCounter. ---*/
    if( badElement ) {
      volElem[l].SetConstantSolution(EntropyVarFreeStream.data(), nVar, 0);
      ++ErrorCounter;
    }
  }
  END_SU2_OMP_FOR

  /*--- Determine the total number of bad elements in this partition. ---*/
  SU2_OMP_ATOMIC
  counter += ErrorCounter;

  SU2_OMP_BARRIER

  /*--- Collect the number of non-physical points for this iteration. ---*/
  SU2_OMP_SINGLE
  {
    if (config->GetComm_Level() == COMM_FULL) {
#ifdef HAVE_MPI
      ErrorCounter = counter;
      SU2_MPI::Allreduce(&ErrorCounter, &counter, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
#endif
      if(iMesh == MESH_0) config->SetNonphysical_Points(counter);
    }
  }
  END_SU2_OMP_SINGLE

  /*-----------------------------------------------------------------------------*/
  /*                       Check for grid motion.                                */
  /*-----------------------------------------------------------------------------*/

  const bool harmonic_balance = config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE;
  if(config->GetGrid_Movement() && !harmonic_balance) {

    /*--- Determine the type of grid motion. ---*/
    switch( config->GetKind_GridMovement() ) {

      case ROTATING_FRAME:
      case STEADY_TRANSLATION: {

        /*--- Steady motion. This is accounted for in the general preprocessing,
              where grid velocities are computed. ---*/
        break;
      }

      case RIGID_MOTION: {

        /*--- Rigid body motion described. At the moment this is only
              possible for the Classical Runge Kutta scheme. ---*/
        if(config->GetKind_TimeIntScheme() != CLASSICAL_RK4_EXPLICIT) {
          SU2_OMP_SINGLE
          SU2_MPI::Error("Rigid body motion only possible for CLASSICAL_RK4_EXPLICIT.",
                         CURRENT_FUNCTION);
          END_SU2_OMP_SINGLE
        }

        /*--- Determine whether or not it is needed to compute the motion data. ---*/
        const unsigned long TimeIter = config->GetTimeIter();

        bool computeMotion = false, firstTime = false;
        if(TimeIter == 0 && iStep == 0) computeMotion = firstTime = true;
        if(iStep == 1 || iStep == 3)    computeMotion = true;

        if( computeMotion ) {

          /*--- Determine the number of time levels. For the classical Runge-Kutta
                scheme this should be 1. ---*/
          const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();

          /*--- Determine the time for which the motion data must be determined. ---*/
          const su2double deltaT = config->GetDelta_UnstTimeND();

          su2double tNew = TimeIter*deltaT;
          if( iStep ) tNew += 0.25*(iStep+1)*deltaT;

          /*--- Determine the time for which the currently stored position was
                calculated. ---*/
          const su2double tOld = tNew - 0.5*deltaT;

          /*--- Hard code the plunging motion in y-direction. ---*/
          const su2double b2New    = 0.25*tNew*tNew*(3.0-tNew);
          const su2double db2Newdt = 0.75*tNew*(2.0-tNew);

          const su2double b2Old = 0.25*tOld*tOld*(3.0-tOld);
          const su2double dB2   = b2New - b2Old;

          /*--- Compute the chunk sizes for the mesh points and internal faces. ---*/
#ifdef HAVE_OMP
          const size_t omp_chunk_size_poin = computeStaticChunkSize(nMeshPoints, omp_get_num_threads(), 64);
          const size_t omp_chunk_size_face = computeStaticChunkSize(nMatchingInternalFacesWithHaloElem[nTimeLevels],
                                                                    omp_get_num_threads(), 64);
#endif

          /*-------------------------------------------------------------------*/
          /*--- Update the coordinates, if this is not the first time step. ---*/
          /*-------------------------------------------------------------------*/

          if( !firstTime ) {

            /*--- Update the grid points. ---*/
            SU2_OMP_FOR_STAT(omp_chunk_size_poin)
            for(unsigned long l=0; l<nMeshPoints; ++l)
              meshPoints[l].coor[1] += dB2;
	    END_SU2_OMP_FOR

            /*--- Update the coordinates of the volume elements. ---*/
            SU2_OMP_FOR_STAT(omp_chunk_size_elem)
            for(unsigned long l=0; l<nVolElemOwned; ++l) {

              /*--- Adapt the coordinates of the volume integration points. ---*/
              unsigned short nItems = volElem[l].coorIntegrationPoints.rows();
              SU2_OMP_SIMD_IF_NOT_AD
              for(unsigned short i=0; i<nItems; ++i)
                volElem[l].coorIntegrationPoints(i,1) += dB2;

              /*--- Update the coordinates of the grid DOFs. ---*/
              nItems = volElem[l].coorGridDOFs.rows();
              SU2_OMP_SIMD_IF_NOT_AD
              for(unsigned short i=0; i<nItems; ++i)
                volElem[l].coorGridDOFs(i,1) += dB2;

              /*--- Update the coordinates of the solution DOFs. ---*/
              nItems = volElem[l].coorSolDOFs.rows();
              SU2_OMP_SIMD_IF_NOT_AD
              for(unsigned short i=0; i<nItems; ++i)
                volElem[l].coorSolDOFs(i,1) += dB2;
            }
            END_SU2_OMP_FOR

            /*--- Loop over the internal matching faces. ---*/
            SU2_OMP_FOR_STAT(omp_chunk_size_face)
            for(unsigned long l=0; l<nMatchingInternalFacesWithHaloElem[nTimeLevels]; ++l) {

              const unsigned short nItems = matchingInternalFaces[l].coorIntegrationPoints.rows();
              SU2_OMP_SIMD_IF_NOT_AD
              for(unsigned short i=0; i<nItems; ++i)
                matchingInternalFaces[l].coorIntegrationPoints(i,1) += dB2;
            }
            END_SU2_OMP_FOR

            /*--- The physical boundary faces. Exclude the periodic boundaries,
                  because these are not physical boundaries. ---*/
            for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
              if( !(boundaries[iMarker].periodicBoundary) ) {

                /*--- Determine the chunk size for the surface elements. ---*/
                const unsigned long nSurfElem = boundaries[iMarker].surfElem.size();
#ifdef HAVE_OMP
                const size_t omp_chunk_size_surf = computeStaticChunkSize(nSurfElem, omp_get_num_threads(), 64);
#endif
                /*--- Loop over the number of surface elements. ---*/
                SU2_OMP_FOR_STAT(nSurfElem)
                for(unsigned long l=0; l<nSurfElem; ++l) {

                  const unsigned short nItems = boundaries[iMarker].surfElem[l].coorIntegrationPoints.rows();
                  SU2_OMP_SIMD_IF_NOT_AD
                  for(unsigned short i=0; i<nItems; ++i)
                    boundaries[iMarker].surfElem[l].coorIntegrationPoints(i,1) += dB2;
                }
                END_SU2_OMP_FOR
              }
            }
          }

          /*-------------------------------------------------------------------*/
          /*---                Compute all the grid velocities.             ---*/
          /*-------------------------------------------------------------------*/

          /*--- Loop over the volume elements. ---*/
          SU2_OMP_FOR_STAT(omp_chunk_size_elem)
          for(unsigned long l=0; l<nVolElemOwned; ++l) {

            /*--- Grid velocities of the integration points. ---*/
            unsigned short nItems = volElem[l].gridVelocitiesInt.rows();
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nItems; ++i)
              volElem[l].gridVelocitiesInt(i,1) = db2Newdt;

            /*--- Grid velocities of the solution DOFs. ---*/
            nItems = volElem[l].gridVelocitiesSolDOFs.rows();
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nItems; ++i)
              volElem[l].gridVelocitiesSolDOFs(i,1) = db2Newdt;
          }
          END_SU2_OMP_FOR

          /*--- Loop over the internal matching faces. ---*/
          SU2_OMP_FOR_STAT(omp_chunk_size_face)
          for(unsigned long l=0; l<nMatchingInternalFacesWithHaloElem[nTimeLevels]; ++l) {

            /*--- Grid velocities of the integration points. ---*/
            const unsigned short nItems = matchingInternalFaces[l].gridVelocities.rows();
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nItems; ++i)
              matchingInternalFaces[l].gridVelocities(i,1) = db2Newdt;
          }
          END_SU2_OMP_FOR

          /*--- The physical boundary faces. Exclude the periodic boundaries,
                because these are not physical boundaries. ---*/
          for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
            if( !(boundaries[iMarker].periodicBoundary) ) {

              /*--- Determine the chunk size for the surface elements. ---*/
              const unsigned long nSurfElem = boundaries[iMarker].surfElem.size();
#ifdef HAVE_OMP
              const size_t omp_chunk_size_surf = computeStaticChunkSize(nSurfElem, omp_get_num_threads(), 64);
#endif
              /*--- Loop over the number of surface elements. ---*/
              SU2_OMP_FOR_STAT(nSurfElem)
              for(unsigned long l=0; l<nSurfElem; ++l) {

                const unsigned short nItems = boundaries[iMarker].surfElem[l].gridVelocities.rows();
                SU2_OMP_SIMD_IF_NOT_AD
                for(unsigned short i=0; i<nItems; ++i)
                  boundaries[iMarker].surfElem[l].gridVelocities(i,1) = db2Newdt;
              }
              END_SU2_OMP_FOR
            }
          }
        }

        break;
      }

      default: {
        SU2_OMP_SINGLE
        SU2_MPI::Error("Only rigid body motion possible for DG-FEM solver at the moment.",
                       CURRENT_FUNCTION);
        END_SU2_OMP_SINGLE
      }
    }
  }
}

void CFEM_DG_EulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                      unsigned short iMesh) { }

void CFEM_DG_EulerSolver::ComputeSpatialJacobian(CGeometry *geometry,  CSolver **solver_container,
                                                 CNumerics **numerics, CConfig *config,
                                                 unsigned short iMesh, unsigned short RunTime_EqSystem) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::Set_OldSolution() {

  /*--- Loop over owned elements. ---*/
#ifdef HAVE_OMP
  const size_t omp_chunk_size_elem = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif
  SU2_OMP_FOR_STAT(omp_chunk_size_elem)
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /*--- Determine the total number of items to be copied. ---*/
    const unsigned int nItems = volElem[l].solDOFs.rows()*volElem[l].solDOFs.cols();

    /*--- Set the pointers and copy the data. ---*/
    const su2double *solDOFs = volElem[l].solDOFs.data();
    su2double *solDOFsWork   = volElem[l].solDOFsWork.data();

    SU2_OMP_SIMD
    for(unsigned int i=0; i<nItems; ++i)
      solDOFsWork[i] = solDOFs[i];
  }
  END_SU2_OMP_FOR
}

void CFEM_DG_EulerSolver::Set_NewSolution() {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh, unsigned long Iteration) {

  /*--- Initialize the minimum and maximum time step and determine the
        active polynomial degree when grid sequencing is used. ---*/
  SU2_OMP_SINGLE
  {
    DetermineCurrentPInPSequencing(config);
    Min_Delta_Time = 1.e25;
    Max_Delta_Time = 0.0;
  }
  END_SU2_OMP_SINGLE

  /*--- Determine the chunk size for the OMP loops, if supported. ---*/
#ifdef HAVE_OMP
  const size_t omp_chunk_size_elem = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif

  /*--- Check whether or not a time stepping scheme is used and store
        the CFL number a bit easier. Note that if we are using explicit
        time stepping, the regular CFL condition has been overwritten with the
        unsteady CFL condition in the config post-processing (if non-zero). ---*/
  const bool time_stepping = config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING;
  const su2double CFL = config->GetCFL(iMesh);

  /*--- Check for explicit time stepping with imposed time step. If the unsteady
        CFL is set to zero (default), it uses the defined unsteady time step,
        otherwise it computes the time step based on the provided unsteady CFL.
        Note that the regular CFL option in the config is always ignored with
        time stepping. ---*/
  if(time_stepping && (config->GetUnst_CFL() == 0.0)) {

    /*--- Loop over the owned volume elements and set the fixed dt. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size_elem)
    for(unsigned long i=0; i<nVolElemOwned; ++i)
      volElem[i].deltaTime = config->GetDelta_UnstTimeND();
    END_SU2_OMP_FOR

  } else {

    /*--- Define the thread local variables for the minimum and maximum time step. ---*/
    su2double MaxDeltaT = 0.0, MinDeltaT = 1.e25;

    /*--- Loop over the owned volume elements. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size_elem)
    for(unsigned long l=0; l<nVolElemOwned; ++l) {

      /*--- Determine the solution in the integration points. ---*/
      ColMajorMatrix<su2double> &solInt = volElem[l].ComputeSolIntPoints(volElem[l].solDOFs);

      /*--- Abbreviate the number of integration points and its padded version. ---*/
      const unsigned short nInt    = volElem[l].standardElemFlow->GetNIntegration();
      const unsigned short nIntPad = volElem[l].standardElemFlow->GetNIntegrationPad();

      /*--- Make a distinction between two and three space dimensions
            in order to have the most efficient code. ---*/
      switch( nDim ) {

        case 2: {

          /*--- Two dimensional simulation. Loop over the padded integration
                points and compute the inviscid spectral radius squared, which
                is stored in the first entry of solInt. Note that the formulation
                used is a rather conservative estimate. ---*/
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nIntPad; ++i) {

            const su2double pOvRho = -1.0/solInt(i,3);
            const su2double u      =  pOvRho*solInt(i,1) - volElem[l].gridVelocitiesInt(i,0);
            const su2double v      =  pOvRho*solInt(i,2) - volElem[l].gridVelocitiesInt(i,1);
            const su2double a      =  sqrt(fabs(Gamma*pOvRho));

            const su2double radx = fabs(u) + a;
            const su2double rady = fabs(v) + a;

            solInt(i,0) = radx*radx + rady*rady;
          }

          break;
        }

        case 3: {

          /*--- Three dimensional simulation. Loop over the padded integration
                points and compute the inviscid spectral radius squared, which
                is stored in the first entry of solInt. Note that the formulation
                used is a rather conservative estimate. ---*/
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nIntPad; ++i) {

            const su2double pOvRho = -1.0/solInt(i,4);
            const su2double u      =  pOvRho*solInt(i,1) - volElem[l].gridVelocitiesInt(i,0);
            const su2double v      =  pOvRho*solInt(i,2) - volElem[l].gridVelocitiesInt(i,1);
            const su2double w      =  pOvRho*solInt(i,3) - volElem[l].gridVelocitiesInt(i,2);
            const su2double a      =  sqrt(fabs(Gamma*pOvRho));

            const su2double radx = fabs(u) + a;
            const su2double rady = fabs(v) + a;
            const su2double radz = fabs(w) + a;

            solInt(i,0) = radx*radx + rady*rady + radz*radz;
          }

          break;
        }
      }

      /*--- Determine the maximum value of the inviscid spectral radius
            squared for the integration points. ---*/
      su2double charVel2Max = 0.0;
      for(unsigned short i=0; i<nInt; ++i)
        charVel2Max = max(charVel2Max, solInt(i,0));

      /*--- Compute the time step for the element. Note that for the spectral
            radius a correction factor, which is a function of the polynomial degree
            and the element type, must be taken into account. ---*/
      const unsigned short pCurrent = min(currentPInPSequencing, volElem[l].standardElemFlow->GetPolyDegree());
      const passivedouble factInv = volElem[l].standardElemFlow->GetFactorInviscidSpectralRadius(pCurrent);
      volElem[l].deltaTime = CFL*volElem[l].lenScale/(factInv*sqrt(charVel2Max));

      /*--- Update the minimum and maximum value, for which the factor for
            time accurate local time stepping must be taken into account. ---*/
      const su2double dtEff = volElem[l].factTimeLevel*volElem[l].deltaTime;
      MinDeltaT = min(MinDeltaT, dtEff);
      MaxDeltaT = max(MaxDeltaT, dtEff);
    }
    END_SU2_OMP_FOR

    /*--- Update the shared variables Min_Delta_Time and Max_Delta_Time. ---*/
    SU2_OMP_CRITICAL
    {
      Min_Delta_Time = min(Min_Delta_Time, MinDeltaT);
      Max_Delta_Time = max(Max_Delta_Time, MaxDeltaT);
    }
    END_SU2_OMP_CRITICAL

    /*--- Compute the max and the min dt (in parallel). Note that we only
          do this for steady calculations if the high verbosity is set, but we
          always perform the reduction for unsteady calculations where the CFL
          limit is used to set the global time step. ---*/
#ifdef HAVE_MPI
    SU2_OMP_SINGLE
    {
      if ((config->GetComm_Level() == COMM_FULL) || time_stepping) {
        su2double rbuf_time = Min_Delta_Time;
        SU2_MPI::Allreduce(&rbuf_time, &Min_Delta_Time, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());

        rbuf_time = Max_Delta_Time;
        SU2_MPI::Allreduce(&rbuf_time, &Max_Delta_Time, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
      }
    }
    END_SU2_OMP_SINGLE
#endif

    /*--- For explicit time stepping with an unsteady CFL imposed, use the
          minimum delta time of the entire mesh. As Min_Delta_Time is scaled to
          the time step of the largest time level, a correction must be used
          for the time level when time accurate local time stepping is used. ---*/
    if (time_stepping) {
      SU2_OMP_FOR_STAT(omp_chunk_size_elem)
      for(unsigned long l=0; l<nVolElemOwned; ++l)
        volElem[l].deltaTime = Min_Delta_Time/volElem[l].factTimeLevel;
      END_SU2_OMP_FOR

      SU2_OMP_SINGLE
      config->SetDelta_UnstTimeND(Min_Delta_Time);
      END_SU2_OMP_SINGLE
    }
  }
}

void CFEM_DG_EulerSolver::CheckTimeSynchronization(CConfig         *config,
                                                   const su2double TimeSync,
                                                   su2double       &timeEvolved,
                                                   bool            &syncTimeReached) {

  /*--- Only one thread needs to carry out the first part of this function. ---*/
  SU2_OMP_SINGLE
  {
    /*--- Check if this is the first time this check is carried out
          and determine the new time evolved. ---*/
    const bool firstTime = timeEvolved == 0.0;
    timeEvolved         += Min_Delta_Time;

    /*--- Check for a (too) small a value for the synchronization time and
          print a warning if this happens. ---*/
    if(firstTime && timeEvolved >= 1.5*TimeSync) {

      if(rank == MASTER_NODE) {
        cout << endl << "              WARNING" << endl;
        cout << "The specified synchronization time is " << timeEvolved/TimeSync
             << " times smaller than the time step for stability" << endl;
        cout << "This is inefficient!!!!!" << endl << endl;
      }
    }
  }
  END_SU2_OMP_SINGLE

  /*--- If the current value of timeEvolved is larger or equal than the
        synchronization time, syncTimeReached is set to true and a correction
        to the time step is carried out. The factor for the time accurate local
        time stepping must be taken into account. If the synchronization time
        has not been reached yet, syncTimeReached is set to false. ---*/
  if(timeEvolved >= TimeSync) {

    const su2double newDeltaTime = Min_Delta_Time + (TimeSync - timeEvolved);

#ifdef HAVE_OMP
    const size_t omp_chunk_size_vol = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif
    SU2_OMP_FOR_STAT(omp_chunk_size_vol)
    for(unsigned long i=0; i<nVolElemOwned; ++i)
      volElem[i].deltaTime = newDeltaTime/volElem[i].factTimeLevel;
    END_SU2_OMP_FOR

    SU2_OMP_SINGLE
    syncTimeReached = true;
    END_SU2_OMP_SINGLE
  }
  else {
    SU2_OMP_SINGLE
    syncTimeReached = false;
    END_SU2_OMP_SINGLE
  }
}

void CFEM_DG_EulerSolver::ProcessTaskList_DG(CGeometry *geometry,  CSolver **solver_container,
                                             CNumerics **numerics, CConfig *config,
                                             unsigned short iMesh) {

  /*--- Initialize the bool vector, that indicates whether or
        not the tasks from the list have been completed. ---*/
  SU2_OMP_SINGLE
  taskCompleted.assign(tasksList.size(), false);
  END_SU2_OMP_SINGLE

  /*--- Easier storage of the number of time levels. ---*/
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();

  /*--- While loop to carry out all the tasks in tasksList. ---*/
  unsigned long lowestIndexInList = 0;
  while(lowestIndexInList < tasksList.size()) {

    /*--- Find the next task that can be carried out. The outer loop is there
          to make sure that a communication is completed in case there are no
          other tasks ---*/
    for(unsigned short j=0; j<2; ++j) {

      bool taskCarriedOut = false;
      for(unsigned long i=lowestIndexInList; i<tasksList.size(); ++i) {

        /*--- Determine whether or not it can be attempted to carry out
              this task. ---*/
        bool taskCanBeCarriedOut = !taskCompleted[i];
        for(unsigned short ind=0; ind<tasksList[i].nIndMustBeCompleted; ++ind) {
          if( !taskCompleted[tasksList[i].indMustBeCompleted[ind]] )
            taskCanBeCarriedOut = false;
        }
        SU2_OMP_BARRIER

        if( taskCanBeCarriedOut ) {

          /*--- Determine the actual task to be carried out and do so. The
                only tasks that may fail are the completion of the non-blocking
                communication. If that is the case the next task needs to be
                found. ---*/
          switch( tasksList[i].task ) {

            case CTaskDefinition::ADER_PREDICTOR_STEP_COMM_ELEMENTS: {

              /*--- Carry out the ADER predictor step for the elements whose
                    solution must be communicated for this time level. ---*/
              const unsigned short level   = tasksList[i].timeLevel;
              const unsigned long  elemBeg = nVolElemOwnedPerTimeLevel[level]
                                           + nVolElemInternalPerTimeLevel[level];
              const unsigned long  elemEnd = nVolElemOwnedPerTimeLevel[level+1];

              ADER_DG_PredictorStep(config, elemBeg, elemEnd);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::ADER_PREDICTOR_STEP_INTERNAL_ELEMENTS: {

              /*--- Carry out the ADER predictor step for the elements whose
                    solution must not be communicated for this time level. ---*/
              const unsigned short level   = tasksList[i].timeLevel;
              const unsigned long  elemBeg = nVolElemOwnedPerTimeLevel[level];
              const unsigned long  elemEnd = nVolElemOwnedPerTimeLevel[level]
                                           + nVolElemInternalPerTimeLevel[level];
              ADER_DG_PredictorStep(config, elemBeg, elemEnd);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::INITIATE_MPI_COMMUNICATION: {

              /*--- Start the MPI communication of the solution in the halo elements. ---*/
              Initiate_MPI_Communication(config, tasksList[i].timeLevel);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::COMPLETE_MPI_COMMUNICATION: {

              /*--- Attempt to complete the MPI communication of the solution data.
                    For j==0, SU2_MPI::Testall will be used, which returns false if
                    not all requests can be completed. In that case the next task on
                    the list is carried out. If j==1, this means that the next
                    tasks are waiting for this communication to be completed and
                    hence MPI_Waitall is used. ---*/
              if( Complete_MPI_Communication(config, tasksList[i].timeLevel,
                                             j==1) ) {
                taskCarriedOut = true;
                SU2_OMP_SINGLE
                taskCompleted[i] = true;
                END_SU2_OMP_SINGLE
              }
              break;
            }

            case CTaskDefinition::INITIATE_REVERSE_MPI_COMMUNICATION: {

              /*--- Start the communication of the residuals, for which the
                    reverse communication must be used. ---*/
              Initiate_MPI_ReverseCommunication(config, tasksList[i].timeLevel);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::COMPLETE_REVERSE_MPI_COMMUNICATION: {

              /*--- Attempt to complete the MPI communication of the residual data.
                    For j==0, SU2_MPI::Testall will be used, which returns false if
                    not all requests can be completed. In that case the next task on
                    the list is carried out. If j==1, this means that the next
                    tasks are waiting for this communication to be completed and
                    hence MPI_Waitall is used. ---*/
              if( Complete_MPI_ReverseCommunication(config, tasksList[i].timeLevel,
                                                    j==1) ) {
                taskCarriedOut = true;
                SU2_OMP_SINGLE
                taskCompleted[i] = true;
                END_SU2_OMP_SINGLE
              }
              break;
            }

            case CTaskDefinition::ADER_TIME_INTERPOLATE_OWNED_ELEMENTS: {

              /*--- Interpolate the predictor solution of the owned elements
                    in time to the given time integration point for the
                    given time level. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              unsigned long nAdjElem = 0, *adjElem = nullptr;
              if(level < (nTimeLevels-1)) {
                nAdjElem = ownedElemAdjLowTimeLevel[level+1].size();
                adjElem  = ownedElemAdjLowTimeLevel[level+1].data();
              }

              ADER_DG_TimeInterpolatePredictorSol(config, tasksList[i].intPointADER,
                                                  nVolElemOwnedPerTimeLevel[level],
                                                  nVolElemOwnedPerTimeLevel[level+1],
                                                  nAdjElem, adjElem,
                                                  tasksList[i].secondPartTimeIntADER);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::ADER_TIME_INTERPOLATE_HALO_ELEMENTS: {

              /*--- Interpolate the predictor solution of the halo elements
                    in time to the given time integration point for the
                    given time level. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              unsigned long nAdjElem = 0, *adjElem = nullptr;
              if(level < (nTimeLevels-1)) {
                nAdjElem = haloElemAdjLowTimeLevel[level+1].size();
                adjElem  = haloElemAdjLowTimeLevel[level+1].data();
              }

              ADER_DG_TimeInterpolatePredictorSol(config, tasksList[i].intPointADER,
                                                  nVolElemHaloPerTimeLevel[level],
                                                  nVolElemHaloPerTimeLevel[level+1],
                                                  nAdjElem, adjElem,
                                                  tasksList[i].secondPartTimeIntADER);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_OWNED_ELEMENTS: {

              /*--- Compute the artificial viscosity for shock capturing in DG. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              Shock_Capturing_DG(config, nVolElemOwnedPerTimeLevel[level],
                                 nVolElemOwnedPerTimeLevel[level+1]);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::SHOCK_CAPTURING_VISCOSITY_HALO_ELEMENTS: {

              /*--- Compute the artificial viscosity for shock capturing in DG. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              Shock_Capturing_DG(config, nVolElemHaloPerTimeLevel[level],
                                 nVolElemHaloPerTimeLevel[level+1]);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::VOLUME_RESIDUAL: {

              /*--- Compute the volume portion of the residual. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              Volume_Residual(config, nVolElemOwnedPerTimeLevel[level],
                              nVolElemOwnedPerTimeLevel[level+1]);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::SURFACE_RESIDUAL_OWNED_ELEMENTS: {

              /*--- Compute the residual of the faces that only involve owned elements. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              ResidualFaces(config, nMatchingInternalFacesLocalElem[level],
                            nMatchingInternalFacesLocalElem[level+1],
                            numerics[CONV_TERM]);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::SURFACE_RESIDUAL_HALO_ELEMENTS: {

              /*--- Compute the residual of the faces that involve a halo element. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              ResidualFaces(config, nMatchingInternalFacesWithHaloElem[level],
                            nMatchingInternalFacesWithHaloElem[level+1],
                            numerics[CONV_TERM]);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_OWNED: {

              /*--- Apply the boundary conditions that only depend on data
                    of owned elements. ---*/
              Boundary_Conditions(tasksList[i].timeLevel, config, numerics, false);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::BOUNDARY_CONDITIONS_DEPEND_ON_HALO: {

              /*--- Apply the boundary conditions that also depend on data
                    of halo elements. ---*/
              Boundary_Conditions(tasksList[i].timeLevel, config, numerics, true);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_OWNED_ELEMENTS: {

              /*--- Create the final residual by summing up all contributions. ---*/
              CreateFinalResidual(tasksList[i].timeLevel, true);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::SUM_UP_RESIDUAL_CONTRIBUTIONS_HALO_ELEMENTS: {

              /*--- Create the final residual by summing up all contributions. ---*/
              CreateFinalResidual(tasksList[i].timeLevel, false);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_OWNED_ELEMENTS: {

              /*--- Accumulate the space time residuals for the owned elements
                    for ADER-DG. ---*/
              AccumulateSpaceTimeResidualADEROwnedElem(config, tasksList[i].timeLevel,
                                                       tasksList[i].intPointADER);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::ADER_ACCUMULATE_SPACETIME_RESIDUAL_HALO_ELEMENTS: {

              /*--- Accumulate the space time residuals for the halo elements
                    for ADER-DG. ---*/
              AccumulateSpaceTimeResidualADERHaloElem(config, tasksList[i].timeLevel,
                                                      tasksList[i].intPointADER);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::MULTIPLY_INVERSE_MASS_MATRIX: {

              /*--- Multiply the residual by the (lumped) mass matrix, to obtain the final value. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              const bool useADER = config->GetKind_TimeIntScheme() == ADER_DG;
              MultiplyResidualByInverseMassMatrix(config, useADER,
                                                  nVolElemOwnedPerTimeLevel[level],
                                                  nVolElemOwnedPerTimeLevel[level+1]);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            case CTaskDefinition::ADER_UPDATE_SOLUTION: {

              /*--- Perform the update step for ADER-DG. ---*/
              const unsigned short level = tasksList[i].timeLevel;
              ADER_DG_Iteration(nVolElemOwnedPerTimeLevel[level],
                                nVolElemOwnedPerTimeLevel[level+1]);
              taskCarriedOut = true;
              SU2_OMP_SINGLE
              taskCompleted[i] = true;
              END_SU2_OMP_SINGLE
              break;
            }

            default: {

              cout << "Task not defined. This should not happen." << endl;
              exit(1);
            }
          }
        }

        /* Break the inner loop if a task has been carried out. */
        if( taskCarriedOut ) break;
      }

      /* Break the outer loop if a task has been carried out. */
      if( taskCarriedOut ) break;
    }

    /* Update the value of lowestIndexInList. */
    for(; lowestIndexInList < tasksList.size(); ++lowestIndexInList)
      if( !taskCompleted[lowestIndexInList] ) break;
  }
}

void CFEM_DG_EulerSolver::ADER_SpaceTimeIntegration(CGeometry *geometry,  CSolver **solver_container,
                                                    CNumerics **numerics, CConfig *config,
                                                    unsigned short iMesh, unsigned short RunTime_EqSystem) {

  /*--- Preprocessing. ---*/
  Preprocessing(geometry, solver_container, config, iMesh, 0, RunTime_EqSystem, false);
  TolerancesADERPredictorStep();

  /*--- Process the tasks list to carry out one ADER space time integration step. ---*/
  ProcessTaskList_DG(geometry, solver_container, numerics, config, iMesh);

  /*--- Postprocessing. ---*/
  Postprocessing(geometry, solver_container, config, iMesh);
}

void CFEM_DG_EulerSolver::TolerancesADERPredictorStep(void) {

  /*--- Initialize the tolerance vector. ---*/
  SU2_OMP_SINGLE
  TolSolADER.assign(nVar, 0.0);
  END_SU2_OMP_SINGLE

  /*--- Parallel loop over the owned elements to determine
        the local values of the tolerances. ---*/
  su2double WRef[5] = {0.0};
#ifdef HAVE_OMP
  const size_t omp_chunk_size_vol = computeStaticChunkSize(nVolElemTot, omp_get_num_threads(), 64);
#endif
  SU2_OMP_FOR_STAT(omp_chunk_size_vol)
  for(unsigned long i=0; i<nVolElemTot; ++i) {

    /*--- Determine the value of the first basis function
          for this element. ---*/
    const passivedouble basis0 = volElem[i].standardElemFlow->ValBasis0();

    /*--- Loop over the number of variables and update WRef. ---*/
    for(unsigned short iVar=0; iVar<nVar; ++iVar)
      WRef[iVar] = max(WRef[iVar], basis0*volElem[i].solDOFs(0,iVar));
  }
  END_SU2_OMP_FOR

  /*--- Determine the maximum values of all threads. ---*/
  SU2_OMP_CRITICAL
  {
    for(unsigned short iVar=0; iVar<nVar; ++iVar)
      TolSolADER[iVar] = max(TolSolADER[iVar], WRef[iVar]);
  }
  END_SU2_OMP_CRITICAL
  SU2_OMP_BARRIER

  /*--- Determine the global maximum when MPI is used. ---*/
#ifdef HAVE_MPI
  SU2_OMP_SINGLE
  {
    for(unsigned short iVar=0; iVar<nVar; ++iVar)
      WRef[iVar] = TolSolADER[iVar];

    SU2_MPI::Allreduce(WRef, TolSolADER.data(), nVar, MPI_DOUBLE, MPI_MAX,
                       SU2_MPI::GetComm());
  }
  END_SU2_OMP_SINGLE
#endif
}

void CFEM_DG_EulerSolver::ADER_DG_PredictorStep(CConfig             *config,
                                                const unsigned long elemBeg,
                                                const unsigned long elemEnd) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::ADER_DG_AliasedPredictorResidual_2D(CConfig              *config,
                                                              CVolumeElementFEM_DG *elem,
                                                              const su2double      *sol,
                                                              const unsigned short nSimul,
                                                              const unsigned short NPad,
                                                              su2double            *res,
                                                              su2double            *work) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::ADER_DG_AliasedPredictorResidual_3D(CConfig              *config,
                                                              CVolumeElementFEM_DG *elem,
                                                              const su2double      *sol,
                                                              const unsigned short nSimul,
                                                              const unsigned short NPad,
                                                              su2double            *res,
                                                              su2double            *work) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::ADER_DG_NonAliasedPredictorResidual_2D(CConfig              *config,
                                                                 CVolumeElementFEM_DG *elem,
                                                                 const su2double      *sol,
                                                                 const unsigned short nSimul,
                                                                 const unsigned short NPad,
                                                                 su2double            *res,
                                                                 su2double            *work) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::ADER_DG_NonAliasedPredictorResidual_3D(CConfig              *config,
                                                                 CVolumeElementFEM_DG *elem,
                                                                 const su2double      *sol,
                                                                 const unsigned short nSimul,
                                                                 const unsigned short NPad,
                                                                 su2double            *res,
                                                                 su2double            *work) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::ADER_DG_TimeInterpolatePredictorSol(CConfig             *config,
                                                              const unsigned short iTime,
                                                              const unsigned long  elemBeg,
                                                              const unsigned long  elemEnd,
                                                              const unsigned long  nAdjElem,
                                                              const unsigned long  *adjElem,
                                                              const bool           secondPartTimeInt) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::Shock_Capturing_DG(CConfig             *config,
                                             const unsigned long elemBeg,
                                             const unsigned long elemEnd) {

  /*--- Run shock capturing algorithm ---*/
  switch( config->GetKind_FEM_DG_Shock() ) {
    case FEM_SHOCK_CAPTURING_DG::NONE:
      break;
    case FEM_SHOCK_CAPTURING_DG::PERSSON:
      // to be done
      //Shock_Capturing_DG_Persson(elemBeg, elemEnd, workArray);
      break;
  }
}

void CFEM_DG_EulerSolver::Volume_Residual(CConfig             *config,
                                          const unsigned long elemBeg,
                                          const unsigned long elemEnd) {

  /*--- Abbreviation of 1/(Gamma-1). ---*/
  const su2double ovgm1 = 1.0/Gamma_Minus_One;

  /*--- Determine the chunk size for the OpenMP parallelization. ---*/
#ifdef HAVE_OMP
    const unsigned long nElem = elemEnd - elemBeg;
    const size_t omp_chunk_size = computeStaticChunkSize(nElem, omp_get_num_threads(), 64);
#endif

  /*--- Loop over the given element range. ---*/
  SU2_OMP_FOR_DYN(omp_chunk_size)
  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Determine the primitive solution in the integration      ---*/
    /*---         points and store the reference to memory for the fluxes. ---*/
    /*------------------------------------------------------------------------*/

    ColMajorMatrix<su2double>          &solInt = volElem[l].ComputeSolIntPoints(volElem[l].solDOFsWork);
    vector<ColMajorMatrix<su2double> > &fluxes = volElem[l].standardElemFlow->workGradSolInt[omp_get_thread_num()];

    EntropyToPrimitiveVariables(solInt);

    /*--- Compute the transformation matrix between conservative and
          entropy variables in the integration poins. ---*/
    VolumeTransformationMatrix(&volElem[l], solInt);

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the inviscid fluxes, multiplied by minus the     ---*/
    /*---         integration weight, in the integration points.           ---*/
    /*---         If needed, also the source terms are computed.           ---*/
    /*------------------------------------------------------------------------*/

    const unsigned short nIntPad = volElem[l].standardElemFlow->GetNIntegrationPad();
    const passivedouble *weights = volElem[l].standardElemFlow->GetIntegrationWeights();

    /*--- Make a distinction between two and three space dimensions
          in order to have the most efficient code. ---*/
    switch( nDim ) {

      case 2: {

        /*--- Two dimensional simulation. Easier storage of the metric terms. ---*/
        ColMajorMatrix<su2double> &dParDx = volElem[l].metricTermsInt[0];
        ColMajorMatrix<su2double> &dParDy = volElem[l].metricTermsInt[1];
        
        /*--- Loop over the padded number of integration points, such that vectorization
              is more efficient, and compute the fluxes in r- and s-direction. ---*/
        SU2_OMP_SIMD_IF_NOT_AD
        for(unsigned short i=0; i<nIntPad; ++i) {

          /*--- Easier storage of the primitive variables. ---*/
          const su2double rho = solInt(i,0);
          const su2double u   = solInt(i,1);
          const su2double v   = solInt(i,2);
          const su2double p   = solInt(i,3);

          /*--- Compute the total energy per unit volume and the
                relative velocities. ---*/
          const su2double rE   = ovgm1*p + 0.5*rho*(u*u + v*v);
          const su2double uRel = u - volElem[l].gridVelocitiesInt(i,0);
          const su2double vRel = v - volElem[l].gridVelocitiesInt(i,1);

          /*--- Compute the metric terms multiplied by minus the integration weight.
                The minus sign comes from the integration by parts in the weak
                formulation. ---*/
          const su2double drdx = -weights[i]*dParDx(i,0);
          const su2double dsdx = -weights[i]*dParDx(i,1);

          const su2double drdy = -weights[i]*dParDy(i,0);
          const su2double dsdy = -weights[i]*dParDy(i,1);

          /*--- Fluxes in r-direction. */
          const su2double Ur = uRel*drdx + vRel*drdy;

          fluxes[0](i,0) = Ur*rho;
          fluxes[0](i,1) = Ur*rho*u + p*drdx;
          fluxes[0](i,2) = Ur*rho*v + p*drdy;
          fluxes[0](i,3) = Ur*rE + p*(u*drdx + v*drdy);

          /*--- Fluxes in s-direction. */
          const su2double Us = uRel*dsdx + vRel*dsdy;

          fluxes[1](i,0) = Us*rho;
          fluxes[1](i,1) = Us*rho*u + p*dsdx;
          fluxes[1](i,2) = Us*rho*v + p*dsdy;
          fluxes[1](i,3) = Us*rE + p*(u*dsdx + v*dsdy);
        }

        break;
      }

      /*----------------------------------------------------------------------*/

      case 3: {

        /*--- Three dimensional simulation. Easier storage of the metric terms. ---*/
        ColMajorMatrix<su2double> &dParDx = volElem[l].metricTermsInt[0];
        ColMajorMatrix<su2double> &dParDy = volElem[l].metricTermsInt[1];
        ColMajorMatrix<su2double> &dParDz = volElem[l].metricTermsInt[2];

        /*--- Loop over the padded number of integration points, such that vectorization
              is more efficient, and compute the fluxes in r-, s- and t-direction. ---*/
        SU2_OMP_SIMD_IF_NOT_AD
        for(unsigned short i=0; i<nIntPad; ++i) {

          /*--- Easier storage of the primitive variables. ---*/
          const su2double rho = solInt(i,0);
          const su2double u   = solInt(i,1);
          const su2double v   = solInt(i,2);
          const su2double w   = solInt(i,3);
          const su2double p   = solInt(i,4);

          /*--- Compute the total energy per unit volume and the
                relative velocities. ---*/
          const su2double rE   = ovgm1*p + 0.5*rho*(u*u + v*v + w*w);
          const su2double uRel = u - volElem[l].gridVelocitiesInt(i,0);
          const su2double vRel = v - volElem[l].gridVelocitiesInt(i,1);
          const su2double wRel = w - volElem[l].gridVelocitiesInt(i,2);

          /*--- Compute the metric terms multiplied by minus the integration weight.
                The minus sign comes from the integration by parts in the weak
                formulation. ---*/
          const su2double drdx = -weights[i]*dParDx(i,0);
          const su2double dsdx = -weights[i]*dParDx(i,1);
          const su2double dtdx = -weights[i]*dParDx(i,2);

          const su2double drdy = -weights[i]*dParDy(i,0);
          const su2double dsdy = -weights[i]*dParDy(i,1);
          const su2double dtdy = -weights[i]*dParDy(i,2);

          const su2double drdz = -weights[i]*dParDz(i,0);
          const su2double dsdz = -weights[i]*dParDz(i,1);
          const su2double dtdz = -weights[i]*dParDz(i,2);

          /*--- Fluxes in r-direction. */
          const su2double Ur = uRel*drdx + vRel*drdy + wRel*drdz;

          fluxes[0](i,0) = Ur*rho;
          fluxes[0](i,1) = Ur*rho*u + p*drdx;
          fluxes[0](i,2) = Ur*rho*v + p*drdy;
          fluxes[0](i,3) = Ur*rho*w + p*drdz;
          fluxes[0](i,4) = Ur*rE + p*(u*drdx + v*drdy + w*drdz);

          /*--- Fluxes in s-direction. */
          const su2double Us = uRel*dsdx + vRel*dsdy + wRel*dsdz;

          fluxes[1](i,0) = Us*rho;
          fluxes[1](i,1) = Us*rho*u + p*dsdx;
          fluxes[1](i,2) = Us*rho*v + p*dsdy;
          fluxes[1](i,3) = Us*rho*w + p*dsdz;
          fluxes[1](i,4) = Us*rE + p*(u*dsdx + v*dsdy + w*dsdz);

          /*--- Fluxes in t-direction. */
          const su2double Ut = uRel*dtdx + vRel*dtdy + wRel*dtdz;

          fluxes[2](i,0) = Ut*rho;
          fluxes[2](i,1) = Ut*rho*u + p*dtdx;
          fluxes[2](i,2) = Ut*rho*v + p*dtdy;
          fluxes[2](i,3) = Ut*rho*w + p*dtdz;
          fluxes[2](i,4) = Ut*rE + p*(u*dtdx + v*dtdy + w*dtdz);
        }

        break;
      }
    }

    /*--- Compute the volume source terms, if needed. ---*/
    const bool addSourceTerms = VolumeSourceTerms(config, &volElem[l], solInt);

    /*------------------------------------------------------------------------*/
    /*--- Step 3: Compute the contribution to the residuals from the       ---*/
    /*---         integration over the volume element.                     ---*/
    /*------------------------------------------------------------------------*/

    /*--- Initialize the residual to zero. ---*/
    volElem[l].resDOFs.setConstant(0.0);

    /*--- Add the contribution from the fluxes. ---*/
    volElem[l].ResidualGradientBasisFunctions(fluxes);

    /*--- Add the contribution from the source terms, if needed.
          The source terms are stored in solInt. ---*/
    if( addSourceTerms ) volElem[l].ResidualBasisFunctions(solInt);
  }
  END_SU2_OMP_FOR
}

bool CFEM_DG_EulerSolver::VolumeSourceTerms(CConfig                  *config,
                                            CVolumeElementFEM_DG     *elem,
                                            ColMajorMatrix<su2double> &sol) {

  /*--- Determine whether a body force term is present and. ---*/
  const bool body_force = config->GetBody_Force();

  if( body_force ) {
    const su2double *body_force_vector = config->GetBody_Force_Vector();
    const unsigned short nIntPad = elem->standardElemFlow->GetNIntegrationPad();
    const passivedouble *weights = elem->standardElemFlow->GetIntegrationWeights();

    /*--- Loop over the padded number of integration points to initialize
          the source term for the density and energy. ---*/
    const unsigned short indE = nDim+1;
    SU2_OMP_SIMD
    for(unsigned short i=0; i<nIntPad; ++i) {
      sol(i,0)    = 0.0;
      sol(i,indE) = 0.0;
    }

    /*--- Loop over the number of dimensions to set the momentum source terms
          and to update the energy source term. Note that the source terms are
          multiplied with minus the integration weight in order to be consistent
          with the formulation of the residual. Also note that for the energy
          source term the absolute velocity must be taken and not the relative. ---*/
    for(unsigned short iDim=0; iDim<nDim; ++iDim) {
      const unsigned short jDim = iDim+1;
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nIntPad; ++i) {
        const su2double abv = -weights[i]*elem->JacobiansInt(i)*body_force_vector[iDim];
        sol(i,indE) += abv*sol(i,jDim);
        sol(i,jDim)  = abv;
      }
    }
  }

  /*--- Initialize addSourceTerms to body_force. The value of addSourceTerms
        will be set to true when a manufactured solution is computed. ---*/
  bool addSourceTerms = body_force;

  /*--- Check whether or not a manufactured solution is used. ---*/
  if( VerificationSolution ) {
    if( VerificationSolution->IsManufacturedSolution() ) {

      /*--- Get the physical time if necessary. ---*/
      su2double time = 0.0;
      if (config->GetTime_Marching() != TIME_MARCHING::STEADY) time = config->GetPhysicalTime();

      /*--- Store the number of integration points and the integration weights. ---*/
      const unsigned short nIntPad = elem->standardElemFlow->GetNIntegrationPad();
      const unsigned short nInt    = elem->standardElemFlow->GetNIntegration();
      const passivedouble *weights = elem->standardElemFlow->GetIntegrationWeights();

      /*--- For the manufactured solutions a source term must be added. If a
            standard source term has not been specified, initialize the source
            terms, stored in sol, to zero and set addSourceTerms to true. ---*/
      addSourceTerms = true;
      if( !body_force ) {
        for(unsigned short iVar=0; iVar<nVar; ++iVar) {
          SU2_OMP_SIMD
          for(unsigned short i=0; i<nIntPad; ++i)
            sol(i,iVar) = 0.0;
        }
      }

      /*--- Loop over the integration points. ---*/
      for(unsigned short i=0; i<nInt; ++i) {

        /*--- Copy the coordinates of this integration point and determine
              the source terms for the manufactured solution. ---*/
        su2double coor[3] = {0.0};
        for(unsigned iDim=0; iDim<nDim; ++iDim)
          coor[iDim] = elem->coorIntegrationPoints(i,iDim);

        su2double sourceMan[5] = {0.0};
        VerificationSolution->GetMMSSourceTerm(coor, time, sourceMan);

        /*--- Subtract the source term of the manufactured solution, multiplied
              by the appropriate weight, from the possibly earlier computed
              source term. It is subtracted in order to be consistent with
              the definition of the residual used in this code. ---*/
        const su2double weightJac = weights[i]*elem->JacobiansInt(i);
        for(unsigned short iVar=0; iVar<nVar; ++iVar)
          sol(i,iVar) -= weightJac*sourceMan[iVar];
      }
    }
  }

  /*--- Return the value of addSourceTerms. ---*/
  return addSourceTerms;
}

void CFEM_DG_EulerSolver::VolumeTransformationMatrix(CVolumeElementFEM_DG     *elem,
                                                     ColMajorMatrix<su2double> &sol) {

  /*--- Easier storage of the padded number of integration points
        and the inverse of Gamma-1. ---*/
  const unsigned short nIntPad = elem->standardElemFlow->GetNIntegrationPad();
  const su2double      ovgm1   = 1.0/Gamma_Minus_One;

  /*--- Make a distinction between 2D and 3D. ---*/
  switch( nDim ) {
    case 2: {

      /*--- Two dimensional simulation. Loop over the integration points. ---*/
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nIntPad; ++i) {

        /*--- Compute the necessary variables that appear in the
              transformation matrix. ---*/
        const su2double rho = sol(i,0);
        const su2double u   = sol(i,1);
        const su2double v   = sol(i,2);
        const su2double p   = sol(i,3);

        const su2double ru   = rho*u;
        const su2double rv   = rho*v;
        const su2double eKin = 0.5*(u*u + v*v);
        const su2double rE   = p*ovgm1 + rho*eKin;
        const su2double rH   = rE + p;

        /*--- Store the elements of the transformation matrix dUdV in this
              integration point. Note that this matrix is symmetric and hence
              only the upper-diagonal part (or lower diagonal part) is stored. ---*/
        elem->dUdVInt(i,0) = rho;                // dUdV(0,0)
        elem->dUdVInt(i,1) = ru;                 // dUdV(0,1) = dUdV(1,0)
        elem->dUdVInt(i,2) = rv;                 // dUdV(0,2) = dUdV(2,0)
        elem->dUdVInt(i,3) = rE;                 // dUdV(0,3) = dUdV(3,0)
        elem->dUdVInt(i,4) = ru*u + p;           // dUdV(1,1)
        elem->dUdVInt(i,5) = ru*v;               // dUdV(1,2) = dUdV(2,1)
        elem->dUdVInt(i,6) = rH*u;               // dUdV(1,3) = dUdV(3,1)
        elem->dUdVInt(i,7) = rv*v + p;           // dUdV(2,2)
        elem->dUdVInt(i,8) = rH*v;               // dUdV(2,3) = dUdV(3,2)
        elem->dUdVInt(i,9) = rE*rH/rho + p*eKin; // dUdV(3,3)
      }

      break;
    }

    case 3: {

      /*--- Two dimensional simulation. Loop over the integration points. ---*/
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nIntPad; ++i) {

        /*--- Compute the necessary variables that appear in the
              transformation matrix. ---*/
        const su2double rho = sol(i,0);
        const su2double u   = sol(i,1);
        const su2double v   = sol(i,2);
        const su2double w   = sol(i,3);
        const su2double p   = sol(i,4);

        const su2double ru   = rho*u;
        const su2double rv   = rho*v;
        const su2double rw   = rho*w;
        const su2double eKin = 0.5*(u*u + v*v + w*w);
        const su2double rE   = p*ovgm1 + rho*eKin;
        const su2double rH   = rE + p;

        /*--- Store the elements of the transformation matrix dUdV in this
              integration point. Note that this matrix is symmetric and hence
              only the upper-diagonal part (or lower diagonal part) is stored. ---*/
        elem->dUdVInt(i, 0) = rho;                // dUdV(0,0)
        elem->dUdVInt(i, 1) = ru;                 // dUdV(0,1) = dUdV(1,0)
        elem->dUdVInt(i, 2) = rv;                 // dUdV(0,2) = dUdV(2,0)
        elem->dUdVInt(i, 3) = rw;                 // dUdV(0,3) = dUdV(3,0)
        elem->dUdVInt(i, 4) = rE;                 // dUdV(0,4) = dUdV(4,0)
        elem->dUdVInt(i, 5) = ru*u + p;           // dUdV(1,1)
        elem->dUdVInt(i, 6) = ru*v;               // dUdV(1,2) = dUdV(2,1)
        elem->dUdVInt(i, 7) = ru*w;               // dUdV(1,3) = dUdV(3,1)
        elem->dUdVInt(i, 8) = rH*u;               // dUdV(1,4) = dUdV(4,1)
        elem->dUdVInt(i, 9) = rv*v + p;           // dUdV(2,2)
        elem->dUdVInt(i,10) = rv*w;               // dUdV(2,3) = dUdV(3,2)
        elem->dUdVInt(i,11) = rH*v;               // dUdV(2,4) = dUdV(4,2)
        elem->dUdVInt(i,12) = rw*w + p;           // dUdV(3,3)
        elem->dUdVInt(i,13) = rH*w;               // dUdV(3,4) = dUdV(4,3)
        elem->dUdVInt(i,14) = rE*rH/rho + p*eKin; // dUdV(4,4)
      }

      break;
    }
  }
}

void CFEM_DG_EulerSolver::Boundary_Conditions(const unsigned short timeLevel,
                                              CConfig              *config,
                                              CNumerics            **numerics,
                                              const bool           haloInfoNeededForBC){

  /* Loop over all boundaries. */
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    /* Check if this boundary marker must be treated at all. */
    if(boundaries[iMarker].haloInfoNeededForBC == haloInfoNeededForBC) {

      /* Determine the range of faces for this time level and test if any
         surface element for this marker must be treated at all. */
      const unsigned long surfElemBeg = boundaries[iMarker].nSurfElem[timeLevel];
      const unsigned long surfElemEnd = boundaries[iMarker].nSurfElem[timeLevel+1];

      if(surfElemEnd > surfElemBeg) {

        /* Set the pointer for the boundary faces for readability. */
        CSurfaceElementFEM *surfElem = boundaries[iMarker].surfElem.data();

        /* Apply the appropriate boundary condition. */
        switch (config->GetMarker_All_KindBC(iMarker)) {
          case EULER_WALL:
            BC_Euler_Wall(config, surfElemBeg, surfElemEnd, surfElem,
                          numerics[CONV_BOUND_TERM]);
            break;
          case FAR_FIELD:
            BC_Far_Field(config, surfElemBeg, surfElemEnd, surfElem,
                         numerics[CONV_BOUND_TERM]);
            break;
          case SYMMETRY_PLANE:
            BC_Sym_Plane(config, surfElemBeg, surfElemEnd, surfElem,
                         numerics[CONV_BOUND_TERM]);
            break;
          case SUPERSONIC_INLET: /* Use far field for this. When a more detailed state
                                    needs to be specified, use a Riemann boundary. */
            BC_Far_Field(config, surfElemBeg, surfElemEnd, surfElem,
                         numerics[CONV_BOUND_TERM]);
            break;
          case SUPERSONIC_OUTLET:
            BC_Supersonic_Outlet(config, surfElemBeg, surfElemEnd, surfElem,
                                 numerics[CONV_BOUND_TERM]);
            break;
          case INLET_FLOW:
            BC_Inlet(config, surfElemBeg, surfElemEnd, surfElem,
                     numerics[CONV_BOUND_TERM], iMarker);
            break;
          case OUTLET_FLOW:
            BC_Outlet(config, surfElemBeg, surfElemEnd, surfElem,
                      numerics[CONV_BOUND_TERM], iMarker);
            break;
          case ISOTHERMAL:
            BC_Isothermal_Wall(config, surfElemBeg, surfElemEnd, surfElem,
                               numerics[CONV_BOUND_TERM], iMarker);
            break;
          case HEAT_FLUX:
            BC_HeatFlux_Wall(config, surfElemBeg, surfElemEnd, surfElem,
                             numerics[CONV_BOUND_TERM], iMarker);
            break;
          case RIEMANN_BOUNDARY:
            BC_Riemann(config, surfElemBeg, surfElemEnd, surfElem,
                       numerics[CONV_BOUND_TERM], iMarker);
            break;
          case CUSTOM_BOUNDARY:
            BC_Custom(config, surfElemBeg, surfElemEnd, surfElem,
                      numerics[CONV_BOUND_TERM]);
            break;
          case PERIODIC_BOUNDARY:  // Nothing to be done for a periodic boundary.
            break;
          default:
            SU2_MPI::Error("BC not implemented.", CURRENT_FUNCTION);
        }
      }
    }
  }
}

void CFEM_DG_EulerSolver::ResidualFaces(CConfig             *config,
                                        const unsigned long indFaceBeg,
                                        const unsigned long indFaceEnd,
                                        CNumerics           *numerics) {

  /*--- Determine the chunk size for the OpenMP parallelization. ---*/
#ifdef HAVE_OMP
    const unsigned long nFaces = indFaceEnd - indFaceBeg;
    const size_t omp_chunk_size = computeStaticChunkSize(nFaces, omp_get_num_threads(), 64);
#endif

  /*--- Loop over the given face range. ---*/
  SU2_OMP_FOR_DYN(omp_chunk_size)
  for(unsigned long l=indFaceBeg; l<indFaceEnd; ++l) {

    /*--- Abbreviate the number of padded integration points
          and the integration weights. ---*/
    const unsigned short nIntPad = matchingInternalFaces[l].standardElemFlow->GetNIntegrationPad();
    const passivedouble *weights = matchingInternalFaces[l].standardElemFlow->GetIntegrationWeights();

    /*------------------------------------------------------------------------*/
    /*--- Step 1: Compute the inviscid fluxes in the integration points of ---*/
    /*---         this matching face multiplied by the integration weight. ---*/
    /*------------------------------------------------------------------------*/

    /*--- Compute the primitive variables of the left and right state
          in the integration point of the face. ---*/
    ColMajorMatrix<su2double> &solIntLeft  = matchingInternalFaces[l].ComputeSolSide0IntPoints(volElem);
    ColMajorMatrix<su2double> &solIntRight = matchingInternalFaces[l].ComputeSolSide1IntPoints(volElem);

    EntropyToPrimitiveVariables(solIntLeft);
    EntropyToPrimitiveVariables(solIntRight);

    /*--- Compute the invisid fluxes in the integration points of the face. ---*/
    const unsigned int indFlux = omp_get_num_threads() + omp_get_thread_num();
    ColMajorMatrix<su2double> &fluxes = matchingInternalFaces[l].standardElemFlow->elem0->workSolInt[indFlux];
    ComputeInviscidFluxesFace(config, solIntLeft, solIntRight, matchingInternalFaces[l].JacobiansFace,
                              matchingInternalFaces[l].metricNormalsFace,
                              matchingInternalFaces[l].gridVelocities, numerics, fluxes);

    /*--- Multiply the fluxes with the integration weight of the
          corresponding integration point. ---*/
    for(unsigned short j=0; j<nVar; ++j) {
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nIntPad; ++i)
        fluxes(i,j) *= weights[i];
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 2: Compute the contribution to the residuals of the DOFs of ---*/
    /*---         the elements on the left and right side of the face.     ---*/
    /*------------------------------------------------------------------------*/

    /*--- Initialize the residual to zero. ---*/
    matchingInternalFaces[l].resDOFsSide0.setConstant(0.0);
    matchingInternalFaces[l].resDOFsSide1.setConstant(0.0);

    /*--- Add the contribution from the fluxes. Note that for side 1 
          the fluxes must be negated because for side 1 the normal
          is inward pointing. ---*/
    matchingInternalFaces[l].ResidualBasisFunctionsSide0(fluxes);
    
    for(unsigned short j=0; j<nVar; ++j) {
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nIntPad; ++i)
        fluxes(i,j) = -fluxes(i,j);
    }

    matchingInternalFaces[l].ResidualBasisFunctionsSide1(fluxes);
  }
  END_SU2_OMP_FOR
}

void CFEM_DG_EulerSolver::AccumulateSpaceTimeResidualADEROwnedElem(
                                                     CConfig             *config,
                                                     const unsigned short timeLevel,
                                                     const unsigned short intPoint) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::AccumulateSpaceTimeResidualADERHaloElem(
                                                     CConfig             *config,
                                                     const unsigned short timeLevel,
                                                     const unsigned short intPoint) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::CreateFinalResidual(const unsigned short timeLevel,
                                              const bool ownedElements) {


  /* Determine the element range for which the final residual must
     be created. */
  unsigned long elemBeg, elemEnd;
  if( ownedElements ) {
    elemBeg = nVolElemOwnedPerTimeLevel[0];
    elemEnd = nVolElemOwnedPerTimeLevel[timeLevel+1];
  }
  else {
    elemBeg = nVolElemHaloPerTimeLevel[0];
    elemEnd = nVolElemHaloPerTimeLevel[timeLevel+1];
  }

  /*--- Determine the chunk size for the OpenMP parallelization. ---*/
#ifdef HAVE_OMP
    const unsigned long nElem = elemEnd - elemBeg;
    const size_t omp_chunk_size = computeStaticChunkSize(nElem, omp_get_num_threads(), 64);
#endif

  /*--- Loop over the given element range. ---*/
  SU2_OMP_FOR_DYN(omp_chunk_size)
  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

    /* For a halo element the residual is initialized to zero. */
    if( !ownedElements ) volElem[l].resDOFs.setConstant(0.0);

    /* Determine the total number of entries in the residual and set
       the pointer to the actual data. */
    const unsigned int nResTot = volElem[l].resDOFs.cols() * volElem[l].resDOFs.rows();
    su2double *resVol = volElem[l].resDOFs.data();

    /* Loop over the boundary faces of this element and add the residual. */
    for(unsigned long j=0; j<volElem[l].boundaryFaceIDs.size(); ++j) {
      const unsigned long iMarker = volElem[l].boundaryFaceIDs[j].long0;
      const unsigned long jj      = volElem[l].boundaryFaceIDs[j].long1;

      const su2double *resSurf = boundaries[iMarker].surfElem[jj].resDOFsElem.data();

      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned int i=0; i<nResTot; ++i) resVol[i] += resSurf[i];
    }

    /* Loop over the internal faces of this element and add the residual. */
    for(unsigned long j=0; j<volElem[l].internalFaceIDs.size(); ++j) {
      const unsigned long jj = volElem[l].internalFaceIDs[j];
      const su2double *resSurf;
      if(matchingInternalFaces[jj].elemID0 == l) resSurf = matchingInternalFaces[jj].resDOFsSide0.data();
      else                                       resSurf = matchingInternalFaces[jj].resDOFsSide1.data();

      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned int i=0; i<nResTot; ++i) resVol[i] += resSurf[i];
    }
  }
  END_SU2_OMP_FOR
}

void CFEM_DG_EulerSolver::MultiplyResidualByInverseMassMatrix(
                                              CConfig            *config,
                                              const bool          useADER,
                                              const unsigned long elemBeg,
                                              const unsigned long elemEnd) {

  /*--- Determine the chunk size for the OpenMP parallelization. ---*/
#ifdef HAVE_OMP
    const unsigned long nElem = elemEnd - elemBeg;
    const size_t omp_chunk_size = computeStaticChunkSize(nElem, omp_get_num_threads(), 64);
#endif

  /*--- Definition of the object for the blas functionalities. ---*/
  CBlasStructure blas;

  /*--- Define the local variables for the reduction of the residuals. ---*/
  su2double L2Res[MAXNVAR] = {0.0}, LInfRes[MAXNVAR] = {0.0}, coorResMax[MAXNVAR][MAXNDIM] = {0.0};
  unsigned long indResMax[MAXNVAR] = {0};

  /*--- Loop over the given element range. ---*/
  SU2_OMP_FOR_DYN(omp_chunk_size)
  for(unsigned long l=elemBeg; l<elemEnd; ++l) {

    /*--- Easier storage of the number of DOFs and integration points
          and the integration weights. ---*/
    const unsigned short nDOFs    = volElem[l].standardElemFlow->GetNDOFs();
    const unsigned short nDOFsPad = volElem[l].standardElemFlow->GetNDOFsPad();
    const unsigned short nIntPad  = volElem[l].standardElemFlow->GetNIntegrationPad();
    const passivedouble *weights  = volElem[l].standardElemFlow->GetIntegrationWeights();

    /*--- Set the reference to the correct residual. This depends
          whether or not the ADER scheme is used. ---*/
    ColMajorMatrix<su2double> &rVec = useADER ? volElem[l].resTotDOFsADER :
                                                volElem[l].resDOFs;

    /*--- The multiplication by the inverse of mass matrix is carried out by solving
          the linear system of equations M x = res, where M is the generalized mass
          matrix  and res the residual. The generalized mass matrix for the entropy
          variables is a positive definite matrix, but, as it depends on the solution,
          it changes during the simulation. Therefore a matrix free preconditioned
          flexible conjugate gradient method is used to solve this system efficiently.
          Abbreviate the required vectors for the algorithm to improve readability. ---*/
    const int thread = omp_get_thread_num();

    ColMajorMatrix<su2double> &zVec   = volElem[l].standardElemFlow->workDOFs[thread][0];
    ColMajorMatrix<su2double> &pVec   = volElem[l].standardElemFlow->workDOFs[thread][1];
    ColMajorMatrix<su2double> &ApVec  = volElem[l].standardElemFlow->workDOFs[thread][2];
    ColMajorMatrix<su2double> &solVec = volElem[l].standardElemFlow->workDOFs[thread][3];

    ColMajorMatrix<su2double> &pInt = volElem[l].standardElemFlow->workSolInt[thread];

    /*--- Easier storage of the transformation matrix dUdV in the
          integration points of the element. ---*/
    const ColMajorMatrix<su2double> &dUdV = volElem[l].dUdVInt;

    /*--- Determine the inverse of the average Jacobian. ---*/
    const su2double JacInv = 1.0/volElem[l].avgJacobian;

    /*--- Determine the maximum value of rVec. This serves as a scaling value
          to determine the convergence. Update the norms of residuals as well. ---*/
    su2double rVecInitMax = 0.0;
    for(unsigned short j=0; j<nVar; ++j) {
      for(unsigned short i=0; i<nDOFs; ++i) {
        rVecInitMax = max(rVecInitMax, fabs(rVec(i,j)));
        const su2double val = fabs(JacInv*rVec(i,j));
        L2Res[j] += val*val;
        if(val > LInfRes[j]) {
          LInfRes[j] = val;
          indResMax[j] = volElem[l].elemIDGlobal;
          for(unsigned short iDim=0; iDim<nDim; ++iDim)
            coorResMax[j][iDim] = volElem[l].coorSolDOFs(i,iDim);
        }
      }
    }

    /*--- Initialize the solution vector to zero. ---*/
    for(unsigned short j=0; j<nVar; ++j) {
      SU2_OMP_SIMD
      for(unsigned short i=0; i<nDOFsPad; ++i)
        solVec(i,j) = 0.0;
    }

    /*--- If the value of rVecInitMax is very small, it does not make sense to
          solve the linear system, because the solution is zero.
          Test if this is not the case. ---*/
    if(rVecInitMax*JacInv > ((su2double) 1.e-15)) {

      /*--- Abbreviations involving gamma. ---*/
      const su2double gm1    =  Gamma_Minus_One;
      const su2double ov1mg  = -1.0/Gamma_Minus_One;
      const su2double govgm1 =  Gamma/gm1;

      /*--- Determine the value of basis function zero, which can be used
            to compute the averaged state of the element. ---*/
      const passivedouble b0 = volElem[l].standardElemFlow->ValBasis0();

      /*--- Make a distinction between two and three space dimensions
            in order to have the most efficient code. ---*/
      switch( nDim ) {

        case 2: {

          /*------------------------------------------------------------------*/
          /*--- Compute the elements of the transformation matrix dVdU     ---*/
          /*--- based on the average solution in the element.              ---*/
          /*------------------------------------------------------------------*/

          /*--- Compute the entropy variables in the center of the element. ---*/
          const su2double V0 = b0*volElem[l].solDOFs(0,0);
          const su2double V1 = b0*volElem[l].solDOFs(0,1);
          const su2double V2 = b0*volElem[l].solDOFs(0,2);
          const su2double V3 = b0*volElem[l].solDOFs(0,3);

          /*--- Compute the primitive variables from the entropy ones. ---*/
          const su2double V3Inv =  1.0/V3;
          const su2double u     = -V3Inv*V1;
          const su2double v     = -V3Inv*V2;
          const su2double eKin  =  0.5*(u*u + v*v);
          const su2double s     =  Gamma - gm1*(V0 - V3*eKin);
          const su2double tmp   = -V3*exp(s);
          const su2double rho   =  pow(tmp, ov1mg);
          const su2double p     = -rho*V3Inv;
          const su2double pInv  =  1.0/p;

          /*--- Two abbreviations that appear in the elements of dVdU. ---*/
          const su2double abv1 = gm1*rho*pInv*pInv;
          const su2double abv2 = eKin*abv1;

          /*--- Store the elements of the transformation matrix dVdU of the
                averaged state. Note that this matrix is symmetric and hence
                only the upper-diagonal part (or lower diagonal part) is stored.
                Multiply the transformation matrix by the inverse of the Jacobian
                to account for the transformation to the standard element. ---*/
          su2double dVdU[10];

          dVdU[0] =  JacInv*(govgm1/rho + abv2*eKin);       // dVdU(0,0)
          dVdU[1] = -JacInv*u*abv2;                         // dVdU(0,1) = dVdU(1,0)
          dVdU[2] = -JacInv*v*abv2;                         // dVdU(0,2) = dVdU(2,0)
          dVdU[3] =  JacInv*(abv2 - pInv);                  // dVdU(0,3) = dVdU(3,0)
          dVdU[4] =  JacInv*(abv1*u*u + pInv);              // dVdU(1,1)
          dVdU[5] =  JacInv*abv1*u*v;                       // dVdU(1,2) = dVdU(2,1)
          dVdU[6] = -JacInv*abv1*u;                         // dVdU(1,3) = dVdU(3,1)
          dVdU[7] =  JacInv*(abv1*v*v + pInv);              // dVdU(2,2)
          dVdU[8] = -JacInv*abv1*v;                         // dVdU(2,3) = dVdU(3,2)
          dVdU[9] =  JacInv*abv1;                           // dVdU(3,3)

          /*------------------------------------------------------------------*/
          /*--- Initialization phase of the preconditioned conjugate       ---*/
          /*--- gradient method.                                           ---*/
          /*------------------------------------------------------------------*/

          /*--- Preconditioning step. Multiply rVec with the dVdU and store the
                result in zVec. Loop over the padded value of nDOFs for
                performance reasons. ---*/
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nDOFsPad; ++i)
            zVec(i,0) = dVdU[0]*rVec(i,0) + dVdU[1]*rVec(i,1)
                      + dVdU[2]*rVec(i,2) + dVdU[3]*rVec(i,3);

          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nDOFsPad; ++i)
            zVec(i,1) = dVdU[1]*rVec(i,0) + dVdU[4]*rVec(i,1)
                      + dVdU[5]*rVec(i,2) + dVdU[6]*rVec(i,3);

          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nDOFsPad; ++i)
            zVec(i,2) = dVdU[2]*rVec(i,0) + dVdU[5]*rVec(i,1)
                      + dVdU[7]*rVec(i,2) + dVdU[8]*rVec(i,3);

          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nDOFsPad; ++i)
            zVec(i,3) = dVdU[3]*rVec(i,0) + dVdU[6]*rVec(i,1)
                      + dVdU[8]*rVec(i,2) + dVdU[9]*rVec(i,3);

          /*--- Copy zVec into pVec. ---*/
          for(unsigned short j=0; j<nVar; ++j) {
          SU2_OMP_SIMD
            for(unsigned short i=0; i<nDOFsPad; ++i)
              pVec(i,j) = zVec(i,j);
          }

          /*--- Loop over the number of iterations. ---*/
          for(int iter=0;;++iter)
          {
            /*--- Safeguard to avoid an infinite loop. ---*/
            if(iter == 50)
              SU2_MPI::Error(string("Convergence not reached for CG algorithm"),
                             CURRENT_FUNCTION);

            /*--- Compute the dot product of rVec and zVec. ---*/
            su2double dotrz = 0.0;
            for(unsigned short j=0; j<nVar; ++j)
              dotrz += blas.dot(nDOFsPad, &rVec(0,j), &zVec(0,j));

            /*--- The matrix vector product of the modified mass matrix and the
                  vector pVec must be determined. This is done in three steps.
                  Step 1 is to interpolate the data of pVec to the integration
                  points. The result will be stored in pInt. ---*/
            volElem[l].standardElemFlow->SolIntPointsDOFsPadded(pVec, pInt);

            /*--- Step 2 of the matrix vector product. Multiply pInt with dUdV,
                  the integration weight and the Jacobian. The loop is carried
                  out over the padded number of integration points for
                  performance reasons. ---*/
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nIntPad; ++i) {
              const su2double intWeight = weights[i]*volElem[l].JacobiansInt(i);

              const su2double p0 = intWeight*pInt(i,0);
              const su2double p1 = intWeight*pInt(i,1);
              const su2double p2 = intWeight*pInt(i,2);
              const su2double p3 = intWeight*pInt(i,3);

              pInt(i,0) = dUdV(i,0)*p0 + dUdV(i,1)*p1 + dUdV(i,2)*p2 + dUdV(i,3)*p3;
              pInt(i,1) = dUdV(i,1)*p0 + dUdV(i,4)*p1 + dUdV(i,5)*p2 + dUdV(i,6)*p3;
              pInt(i,2) = dUdV(i,2)*p0 + dUdV(i,5)*p1 + dUdV(i,7)*p2 + dUdV(i,8)*p3;
              pInt(i,3) = dUdV(i,3)*p0 + dUdV(i,6)*p1 + dUdV(i,8)*p2 + dUdV(i,9)*p3;
            }

            /*--- Step 3 of the matrix vector product. Scatter the results of pInt
                  back to the DOFs. This is the final result, which is stored in
                  ApVec. This array must be initialized to zero, because the
                  function ResidualBasisFunctions accumulates the data.  ---*/
            ApVec.setConstant(0.0);
            volElem[l].standardElemFlow->ResidualBasisFunctions(pInt, ApVec);

            /*--- Determine the dot product between pVec and ApVec. ---*/
            su2double dotpAp = 0.0;
            for(unsigned short j=0; j<nVar; ++j)
              dotpAp += blas.dot(nDOFsPad, &pVec(0,j), &ApVec(0,j));

            /*--- Determine the coefficient alpha, which is the contribution of
                  pVec to the solution and -ApVec to the right hand side. ---*/
            const su2double alpha = dotrz/dotpAp;

            /*--- Compute the new solution and right hand side. ---*/
            for(unsigned short j=0; j<nVar; ++j) {
              blas.axpy(nDOFsPad,  alpha,  &pVec(0,j), &solVec(0,j));
              blas.axpy(nDOFsPad, -alpha, &ApVec(0,j),   &rVec(0,j));
            }

            /*--- Determine the Linf norm of rVec. Needed to check the convergence
                  of the iterative algorithm. ---*/
            su2double rVecMax = 0.0;
            for(unsigned short j=0; j<nVar; ++j)
              for(unsigned short i=0; i<nDOFs; ++i)
                rVecMax = max(rVecMax, fabs(rVec(i,j)));

            /*--- Convergence criterion. At the moment the convergence criterion is
                  hard coded, but could be replaced by a user defined parameter. ---*/
            if(rVecMax < ((su2double) 1.e-10*rVecInitMax)) break;

            /*--- Preconditioning step. Multiply rVec with the dVdU and store the
                  result in zVec. Loop over the padded value of nDOFs for
                  performance reasons. ---*/
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nDOFsPad; ++i)
              zVec(i,0) = dVdU[0]*rVec(i,0) + dVdU[1]*rVec(i,1)
                        + dVdU[2]*rVec(i,2) + dVdU[3]*rVec(i,3);
  
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nDOFsPad; ++i)
              zVec(i,1) = dVdU[1]*rVec(i,0) + dVdU[4]*rVec(i,1)
                        + dVdU[5]*rVec(i,2) + dVdU[6]*rVec(i,3);

            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nDOFsPad; ++i)
              zVec(i,2) = dVdU[2]*rVec(i,0) + dVdU[5]*rVec(i,1)
                        + dVdU[7]*rVec(i,2) + dVdU[8]*rVec(i,3);

            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nDOFsPad; ++i)
              zVec(i,3) = dVdU[3]*rVec(i,0) + dVdU[6]*rVec(i,1)
                        + dVdU[8]*rVec(i,2) + dVdU[9]*rVec(i,3);

            /*--- Compute the value of the coefficient beta. This is a flexible CG,
                  hence the Polak-Ribiere formula must be used. ---*/
            su2double dotzAp = 0.0;
            for(unsigned short j=0; j<nVar; ++j)
              dotzAp += blas.dot(nDOFsPad, &zVec(0,j), &ApVec(0,j));

            const su2double beta = -alpha*dotzAp/dotrz;

            /*--- Compute the new p-vector. ---*/
            for(unsigned short j=0; j<nVar; ++j) {
              SU2_OMP_SIMD_IF_NOT_AD
              for(unsigned short i=0; i<nDOFsPad; ++i)
                pVec(i,j) = beta*pVec(i,j) + zVec(i,j);
            }
          }

          break;
        }

        /*--------------------------------------------------------------------*/

        case 3: {

          /*------------------------------------------------------------------*/
          /*--- Compute the elements of the transformation matrix dVdU     ---*/
          /*--- based on the average solution in the element.              ---*/
          /*------------------------------------------------------------------*/

          /*--- Compute the entropy variables in the center of the element. ---*/
          const su2double V0 = b0*volElem[l].solDOFs(0,0);
          const su2double V1 = b0*volElem[l].solDOFs(0,1);
          const su2double V2 = b0*volElem[l].solDOFs(0,2);
          const su2double V3 = b0*volElem[l].solDOFs(0,3);
          const su2double V4 = b0*volElem[l].solDOFs(0,4);

          /*--- Compute the primitive variables from the entropy ones. ---*/
          const su2double V4Inv =  1.0/V4;
          const su2double u     = -V4Inv*V1;
          const su2double v     = -V4Inv*V2;
          const su2double w     = -V4Inv*V3;
          const su2double eKin  =  0.5*(u*u + v*v + w*w);
          const su2double s     =  Gamma - gm1*(V0 - V4*eKin);
          const su2double tmp   = -V4*exp(s);
          const su2double rho   =  pow(tmp, ov1mg);
          const su2double p     = -rho*V4Inv;
          const su2double pInv  =  1.0/p;

          /*--- Two abbreviations that appear in the elements of dVdU. ---*/
          const su2double abv1 = gm1*rho*pInv*pInv;
          const su2double abv2 = eKin*abv1;

          /*--- Store the elements of the transformation matrix dVdU of the
                averaged state. Note that this matrix is symmetric and hence
                only the upper-diagonal part (or lower diagonal part) is stored.
                Multiply the transformation matrix by the inverse of the Jacobian
                to account for the transformation to the standard element. ---*/
          su2double dVdU[15];

          dVdU[0]  =  JacInv*(govgm1/rho + abv2*eKin);      // dVdU(0,0)
          dVdU[1]  = -JacInv*u*abv2;                        // dVdU(0,1) = dVdU(1,0)
          dVdU[2]  = -JacInv*v*abv2;                        // dVdU(0,2) = dVdU(2,0)
          dVdU[3]  = -JacInv*w*abv2;                        // dVdU(0,3) = dVdU(3,0)
          dVdU[4]  =  JacInv*(abv2 - pInv);                 // dVdU(0,4) = dVdU(4,0)
          dVdU[5]  =  JacInv*(abv1*u*u + pInv);             // dVdU(1,1)
          dVdU[6]  =  JacInv*abv1*u*v;                      // dVdU(1,2) = dVdU(2,1)
          dVdU[7]  =  JacInv*abv1*u*w;                      // dVdU(1,3) = dVdU(3,1)
          dVdU[8]  = -JacInv*abv1*u;                        // dVdU(1,4) = dVdU(4,1)
          dVdU[9]  =  JacInv*(abv1*v*v + pInv);             // dVdU(2,2)
          dVdU[10] =  JacInv*abv1*v*w;                      // dVdU(2,3) = dVdU(3,2)
          dVdU[11] = -JacInv*abv1*v;                        // dVdU(2,4) = dVdU(4,2)
          dVdU[12] =  JacInv*(abv1*w*w + pInv);             // dVdU(3,3)
          dVdU[13] = -JacInv*abv1*w;                        // dVdU(3,4) = dVdU(4,3)
          dVdU[14] =  JacInv*abv1;                          // dVdU(4,4)

          /*------------------------------------------------------------------*/
          /*--- Initialization phase of the preconditioned conjugate       ---*/
          /*--- gradient method.                                           ---*/
          /*------------------------------------------------------------------*/

          /*--- Preconditioning step. Multiply rVec with the dVdU and store the
                result in zVec. Loop over the padded value of nDOFs for
                performance reasons. ---*/
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nDOFsPad; ++i)
            zVec(i,0) = dVdU[0] *rVec(i,0) + dVdU[1] *rVec(i,1) + dVdU[2] *rVec(i,2)
                      + dVdU[3] *rVec(i,3) + dVdU[4] *rVec(i,4);

          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nDOFsPad; ++i)
            zVec(i,1) = dVdU[1] *rVec(i,0) + dVdU[5] *rVec(i,1) + dVdU[6] *rVec(i,2)
                      + dVdU[7] *rVec(i,3) + dVdU[8] *rVec(i,4);

          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nDOFsPad; ++i)
            zVec(i,2) = dVdU[2] *rVec(i,0) + dVdU[6] *rVec(i,1) + dVdU[9] *rVec(i,2)
                      + dVdU[10]*rVec(i,3) + dVdU[11]*rVec(i,4);

          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nDOFsPad; ++i)
            zVec(i,3) = dVdU[3] *rVec(i,0) + dVdU[7] *rVec(i,1) + dVdU[10]*rVec(i,2)
                      + dVdU[12]*rVec(i,3) + dVdU[13]*rVec(i,4);

          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned short i=0; i<nDOFsPad; ++i)
            zVec(i,4) = dVdU[4] *rVec(i,0) + dVdU[8] *rVec(i,1) + dVdU[11]*rVec(i,2)
                      + dVdU[13]*rVec(i,3) + dVdU[14]*rVec(i,4);

          /*--- Copy zVec into pVec. ---*/
          for(unsigned short j=0; j<nVar; ++j) {
          SU2_OMP_SIMD
            for(unsigned short i=0; i<nDOFsPad; ++i)
              pVec(i,j) = zVec(i,j);
          }

          /*------------------------------------------------------------------*/
          /*--- Iterative phase of the preconditioned conjugate gradient   ---*/
          /*--- method.                                                    ---*/
          /*------------------------------------------------------------------*/

          /*--- Loop over the number of iterations. ---*/
          for(int iter=0;;++iter)
          {
            /*--- Safeguard to avoid an infinite loop. ---*/
            if(iter == 50)
              SU2_MPI::Error(string("Convergence not reached for CG algorithm"),
                             CURRENT_FUNCTION);

            /*--- Compute the dot product of rVec and zVec. ---*/
            su2double dotrz = 0.0;
            for(unsigned short j=0; j<nVar; ++j)
              dotrz += blas.dot(nDOFsPad, &rVec(0,j), &zVec(0,j));

            /*--- The matrix vector product of the modified mass matrix and the
                  vector pVec must be determined. This is done in three steps.
                  Step 1 is to interpolate the data of pVec to the integration
                  points. The result will be stored in pInt. ---*/
            volElem[l].standardElemFlow->SolIntPointsDOFsPadded(pVec, pInt);

            /*--- Step 2 of the matrix vector product. Multiply pInt with dUdV,
                  the integration weight and the Jacobian. The loop is carried
                  out over the padded number of integration points for
                  performance reasons. ---*/
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nIntPad; ++i) {
              const su2double intWeight = weights[i]*volElem[l].JacobiansInt(i);

              const su2double p0 = intWeight*pInt(i,0);
              const su2double p1 = intWeight*pInt(i,1);
              const su2double p2 = intWeight*pInt(i,2);
              const su2double p3 = intWeight*pInt(i,3);
              const su2double p4 = intWeight*pInt(i,4);

              pInt(i,0) = dUdV(i,0)*p0 + dUdV(i,1)*p1 + dUdV(i, 2)*p2 + dUdV(i, 3)*p3 + dUdV(i, 4)*p4;
              pInt(i,1) = dUdV(i,1)*p0 + dUdV(i,5)*p1 + dUdV(i, 6)*p2 + dUdV(i, 7)*p3 + dUdV(i, 8)*p4;
              pInt(i,2) = dUdV(i,2)*p0 + dUdV(i,6)*p1 + dUdV(i, 9)*p2 + dUdV(i,10)*p3 + dUdV(i,11)*p4;
              pInt(i,3) = dUdV(i,3)*p0 + dUdV(i,7)*p1 + dUdV(i,10)*p2 + dUdV(i,12)*p3 + dUdV(i,13)*p4;
              pInt(i,4) = dUdV(i,4)*p0 + dUdV(i,8)*p1 + dUdV(i,11)*p2 + dUdV(i,13)*p3 + dUdV(i,14)*p4;
            }

            /*--- Step 3 of the matrix vector product. Scatter the results of pInt
                  back to the DOFs. This is the final result, which is stored in
                  ApVec. This array must be initialized to zero, because the
                  function ResidualBasisFunctions accumulates the data.  ---*/
            ApVec.setConstant(0.0);
            volElem[l].standardElemFlow->ResidualBasisFunctions(pInt, ApVec);

            /*--- Determine the dot product between pVec and ApVec. ---*/
            su2double dotpAp = 0.0;
            for(unsigned short j=0; j<nVar; ++j)
              dotpAp += blas.dot(nDOFsPad, &pVec(0,j), &ApVec(0,j));

            /*--- Determine the coefficient alpha, which is the contribution of
                  pVec to the solution and -ApVec to the right hand side. ---*/
            const su2double alpha = dotrz/dotpAp;
            /*--- Compute the new solution and right hand side. ---*/
            for(unsigned short j=0; j<nVar; ++j) {
              blas.axpy(nDOFsPad,  alpha,  &pVec(0,j), &solVec(0,j));
              blas.axpy(nDOFsPad, -alpha, &ApVec(0,j),   &rVec(0,j));
            }

            /*--- Determine the Linf norm of rVec. Needed to check the convergence
                  of the iterative algorithm. ---*/
            su2double rVecMax = 0.0;
            for(unsigned short j=0; j<nVar; ++j)
              for(unsigned short i=0; i<nDOFs; ++i)
                rVecMax = max(rVecMax, fabs(rVec(i,j)));

            /*--- Convergence criterion. At the moment the convergence criterion is
                  hard coded, but could be replaced by a user defined parameter. ---*/
            if(rVecMax < ((su2double) 1.e-10*rVecInitMax)) break;

            /*--- Preconditioning step. Multiply rVec with the dVdU and store the
                  result in zVec. Loop over the padded value of nDOFs for
                  performance reasons. ---*/
            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nDOFsPad; ++i)
              zVec(i,0) = dVdU[0] *rVec(i,0) + dVdU[1] *rVec(i,1) + dVdU[2] *rVec(i,2)
                        + dVdU[3] *rVec(i,3) + dVdU[4] *rVec(i,4);

            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nDOFsPad; ++i)
              zVec(i,1) = dVdU[1] *rVec(i,0) + dVdU[5] *rVec(i,1) + dVdU[6] *rVec(i,2)
                        + dVdU[7] *rVec(i,3) + dVdU[8] *rVec(i,4);

            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nDOFsPad; ++i)
              zVec(i,2) = dVdU[2] *rVec(i,0) + dVdU[6] *rVec(i,1) + dVdU[9] *rVec(i,2)
                        + dVdU[10]*rVec(i,3) + dVdU[11]*rVec(i,4);

            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nDOFsPad; ++i)
              zVec(i,3) = dVdU[3] *rVec(i,0) + dVdU[7] *rVec(i,1) + dVdU[10]*rVec(i,2)
                        + dVdU[12]*rVec(i,3) + dVdU[13]*rVec(i,4);

            SU2_OMP_SIMD_IF_NOT_AD
            for(unsigned short i=0; i<nDOFsPad; ++i)
              zVec(i,4) = dVdU[4] *rVec(i,0) + dVdU[8] *rVec(i,1) + dVdU[11]*rVec(i,2)
                        + dVdU[13]*rVec(i,3) + dVdU[14]*rVec(i,4);

            /*--- Compute the value of the coefficient beta. This is a flexible CG,
                  hence the Polak-Ribiere formula must be used. ---*/
            su2double dotzAp = 0.0;
            for(unsigned short j=0; j<nVar; ++j)
              dotzAp += blas.dot(nDOFsPad, &zVec(0,j), &ApVec(0,j));

            const su2double beta = -alpha*dotzAp/dotrz;

            /*--- Compute the new p-vector. ---*/
            for(unsigned short j=0; j<nVar; ++j) {
              SU2_OMP_SIMD_IF_NOT_AD
              for(unsigned short i=0; i<nDOFsPad; ++i)
                pVec(i,j) = beta*pVec(i,j) + zVec(i,j);
            }
          }

          break;
        }
      }
    }

    /*------------------------------------------------------------------------*/
    /*--- Copy the data from solVec back into resVec.                      ---*/
    /*------------------------------------------------------------------------*/

    for(unsigned short j=0; j<nVar; ++j) {
      SU2_OMP_SIMD
      for(unsigned short i=0; i<nDOFsPad; ++i)
        rVec(i,j) = solVec(i,j);
    }

    /*--- If the current polynomial degree in the P sequencing is less than
          the polynomial degree of the element, set the residuals of the higher
          DOFs to zero to improve robustness. ---*/
    if(currentPInPSequencing < volElem[l].standardElemFlow->GetPolyDegree()) {
      const unsigned short VTK_Type = volElem[l].standardElemFlow->GetVTK_Type();
      const unsigned short nDOFsUse = CFEMStandardElementBase::GetNDOFsStatic(VTK_Type, currentPInPSequencing);

      for(unsigned short j=0; j<nVar; ++j) {
        for(unsigned short i=nDOFsUse; i<nDOFsPad; ++i)
          rVec(i,j) = 0.0;
      }
    }
  }
  END_SU2_OMP_FOR

  /*--- Reduce residual information over all threads in this rank. ---*/
  SU2_OMP_CRITICAL
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Residual_RMS[iVar] += L2Res[iVar];
    AddRes_Max(iVar, LInfRes[iVar], indResMax[iVar], coorResMax[iVar]);
  }
  END_SU2_OMP_CRITICAL
  SU2_OMP_BARRIER
}

void CFEM_DG_EulerSolver::SetReferenceValues(const CConfig& config) {

  /*--- Evaluate reference values for non-dimensionalization. For dynamic meshes,
   use the motion Mach number as a reference value for computing the force coefficients.
   Otherwise, use the freestream values, which is the standard convention. ---*/
  su2double RefVel2;

  if (dynamic_grid && !config.GetFSI_Simulation()) {
    su2double Gas_Constant = config.GetGas_ConstantND();
    su2double Mach2Vel = sqrt(Gamma * Gas_Constant * Temperature_Inf);
    su2double Mach_Motion = config.GetMach_Motion();
    RefVel2 = pow(Mach_Motion * Mach2Vel, 2);
  }
  else {
    RefVel2 = GeometryToolbox::SquaredNorm(nDim, Velocity_Inf);
  }

  DynamicPressureRef = 0.5 * Density_Inf * RefVel2;
  AeroCoeffForceRef =  DynamicPressureRef * config.GetRefArea();
}

void CFEM_DG_EulerSolver::Pressure_Forces(const CGeometry *geometry, const CConfig *config) {

  /*--- Get some information from config. ---*/
  su2double Alpha = config->GetAoA() * PI_NUMBER / 180.0;
  su2double Beta = config->GetAoS() * PI_NUMBER / 180.0;
  su2double RefArea = config->GetRefArea();
  su2double RefLength = config->GetRefLength();
  auto Origin = config->GetRefOriginMoment(0);

  /*--- Reference pressure is always the far-field value. ---*/
  const su2double RefPressure = Pressure_Inf;

  /*--- Set the reference values and initialization. ---*/
  SU2_OMP_SINGLE
  {
    SetReferenceValues(*config);
    TotalCoeff.setZero();
    AllBoundInvCoeff.setZero();
    SurfaceInvCoeff.setZero();
    SurfaceCoeff.setZero();
  }
  END_SU2_OMP_SINGLE

  const su2double factor = 1.0 / AeroCoeffForceRef;

  /*--- Loop over the Euler and Navier-Stokes markers ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    /* Check if this boundary must be monitored. */
    const unsigned short Monitoring = config->GetMarker_All_Monitoring(iMarker);
    if(Monitoring == YES) {

      /* Easier storage of the boundary condition. */
      const unsigned short Boundary = config->GetMarker_All_KindBC(iMarker);

      /*--- Obtain the origin for the moment computation for a particular marker ---*/
      for(unsigned short iMarker_Monitoring=0; iMarker_Monitoring<config->GetnMarker_Monitoring();
                       ++iMarker_Monitoring) {
        string Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        string Marker_Tag     = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }

      /* Check for a boundary for which the forces must be computed. */
      if((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) ||
         (Boundary == ISOTHERMAL)) {

        /*--- Forces initialization at each Marker ---*/
        SU2_OMP_SINGLE
        {
          InvCoeff.setZero(iMarker);
        }
        END_SU2_OMP_SINGLE

        /* Variables to store the contributions from the individual threads. */
        su2double ForceInviscid[MAXNDIM] = {0.0}, MomentInviscid[MAXNDIM] = {0.0};
        su2double MomentX_Force[MAXNDIM] = {0.0}, MomentY_Force[MAXNDIM] = {0.0},
                  MomentZ_Force[MAXNDIM] = {0.0};

        /* OpenMP parallel loop over the number of faces of this boundary. */
        const unsigned long nFaces   = boundaries[iMarker].surfElem.size();
        CSurfaceElementFEM *surfElem = boundaries[iMarker].surfElem.data();
#ifdef HAVE_OMP
        const size_t omp_chunk_size = computeStaticChunkSize(nFaces, omp_get_num_threads(), 64);
#endif
        SU2_OMP_FOR_DYN(omp_chunk_size)
        for(unsigned long l=0; l<nFaces; ++l) {

          /*--- Compute the variables of in the integration point of the face
                and convert them to primitive variables. ---*/
          ColMajorMatrix<su2double> &solInt = surfElem[l].ComputeSolSide0IntPoints(volElem);
          EntropyToPrimitiveVariables(solInt);

          /* Easier storage of the number of integration points and the weights,
             unit normals, Jacobians and coordinates of the integration points. */
          const unsigned short nInt = surfElem[l].standardElemFlow->GetNIntegration();
          const passivedouble *weights = surfElem[l].standardElemFlow->GetIntegrationWeights();
          const ColMajorMatrix<su2double> &normals = surfElem[l].metricNormalsFace;
          const su2activevector &Jacobians = surfElem[l].JacobiansFace;
          const ColMajorMatrix<su2double> &Coord = surfElem[l].coorIntegrationPoints;

          /*--- Make a distinction between 2D and 3D for performance reasons. ---*/
          switch( nDim ) {
            case 2: {

              /* Two dimensional simulation. Loop over the integration points. */
              for(unsigned short i=0; i<nInt; ++i) {

                /* Easier storage of the pressure and the distance to the origin
                   for the moment computation. */
                const su2double Pressure = solInt(i,3);
                const su2double DistX = Coord(i,0) - Origin[0];
                const su2double DistY = Coord(i,1) - Origin[1];

                /* Update the pressure contribution to the forces. Note that the
                   normal points into the geometry, hence no minus sign. */
                const su2double ForceMag = (Pressure - Pressure_Inf)*weights[i]
                                         * Jacobians[i]*factor;
                const su2double ForceX   = ForceMag*normals(i,0);
                const su2double ForceY   = ForceMag*normals(i,1);

                ForceInviscid[0] += ForceX;
                ForceInviscid[1] += ForceY;

                /* Update the pressure contribution to the moments. */
                MomentInviscid[2] += (ForceY*DistX - ForceX*DistY)/RefLength;
                MomentZ_Force[0] += (-ForceX * Coord(i,1));
                MomentZ_Force[1] += (ForceY * Coord(i,0));
              }
              break;
            }

            case 3: {

                /* Three dimensional simulation. Loop over the integration points. */
              for(unsigned short i=0; i<nInt; ++i) {

                /* Easier storage of the pressure and the distance to the origin
                   for the moment computation. */
                const su2double Pressure = solInt(i,4);
                const su2double DistX = Coord(i,0) - Origin[0];
                const su2double DistY = Coord(i,1) - Origin[1];
                const su2double DistZ = Coord(i,2) - Origin[2];

                /* Update the pressure contribution to the forces. Note that the
                   normal points into the geometry, hence no minus sign. */
                const su2double ForceMag = (Pressure - Pressure_Inf)*weights[i]
                                         * Jacobians[i]*factor;
                const su2double ForceX   = ForceMag*normals(i,0);
                const su2double ForceY   = ForceMag*normals(i,1);
                const su2double ForceZ   = ForceMag*normals(i,2);

                ForceInviscid[0] += ForceX;
                ForceInviscid[1] += ForceY;
                ForceInviscid[2] += ForceZ;

                /* Update the pressure contribution to the moments. */
                MomentInviscid[0] += (ForceZ*DistY - ForceY*DistZ)/RefLength;
                MomentInviscid[1] += (ForceX*DistZ - ForceZ*DistX)/RefLength;
                MomentInviscid[2] += (ForceY*DistX - ForceX*DistY)/RefLength;

                MomentX_Force[1] += (-ForceY * Coord(i,2));
                MomentX_Force[2] += (ForceZ * Coord(i,1));

                MomentY_Force[2] += (-ForceZ * Coord(i,0));
                MomentY_Force[0] += (ForceX * Coord(i,2));

                MomentZ_Force[0] += (-ForceX * Coord(i,1));
                MomentZ_Force[1] += (ForceY * Coord(i,0));
              }

              break;
            }
          }
        }
        END_SU2_OMP_FOR

        /*--- Reduce the forces and moments over all threads in this rank. ---*/
        if (nDim == 2) {
          SU2_OMP_CRITICAL
          {
            InvCoeff.CD[iMarker] += ForceInviscid[0] * cos(Alpha) + ForceInviscid[1] * sin(Alpha);
            InvCoeff.CL[iMarker] += -ForceInviscid[0] * sin(Alpha) + ForceInviscid[1] * cos(Alpha);
            InvCoeff.CMz[iMarker] += MomentInviscid[2];
            InvCoeff.CoPx[iMarker] += MomentZ_Force[1];
            InvCoeff.CoPy[iMarker] += -MomentZ_Force[0];
            InvCoeff.CFx[iMarker] += ForceInviscid[0];
            InvCoeff.CFy[iMarker] += ForceInviscid[1];
          }
          END_SU2_OMP_CRITICAL
          SU2_OMP_BARRIER

          SU2_OMP_SINGLE
          {
            InvCoeff.CEff[iMarker] = InvCoeff.CL[iMarker] / (InvCoeff.CD[iMarker] + EPS);
            InvCoeff.CT[iMarker] = -InvCoeff.CFx[iMarker];
            InvCoeff.CQ[iMarker] = -InvCoeff.CMz[iMarker];
            InvCoeff.CMerit[iMarker] = InvCoeff.CT[iMarker] / (InvCoeff.CQ[iMarker] + EPS);
          }
          END_SU2_OMP_SINGLE
        }

        if (nDim == 3) {
          SU2_OMP_CRITICAL
          {
            InvCoeff.CD[iMarker] += ForceInviscid[0] * cos(Alpha) * cos(Beta) + ForceInviscid[1] * sin(Beta) +
                                    ForceInviscid[2] * sin(Alpha) * cos(Beta);
            InvCoeff.CL[iMarker] += -ForceInviscid[0] * sin(Alpha) + ForceInviscid[2] * cos(Alpha);
            InvCoeff.CSF[iMarker] += -ForceInviscid[0] * sin(Beta) * cos(Alpha) + ForceInviscid[1] * cos(Beta) -
                                     ForceInviscid[2] * sin(Beta) * sin(Alpha);
            InvCoeff.CMx[iMarker] += MomentInviscid[0];
            InvCoeff.CMy[iMarker] += MomentInviscid[1];
            InvCoeff.CMz[iMarker] += MomentInviscid[2];
            InvCoeff.CoPx[iMarker] += -MomentY_Force[0];
            InvCoeff.CoPz[iMarker] += MomentY_Force[2];
            InvCoeff.CFx[iMarker] += ForceInviscid[0];
            InvCoeff.CFy[iMarker] += ForceInviscid[1];
            InvCoeff.CFz[iMarker] += ForceInviscid[2];
          }
          END_SU2_OMP_CRITICAL
          SU2_OMP_BARRIER

          SU2_OMP_SINGLE
          {
            InvCoeff.CEff[iMarker] = InvCoeff.CL[iMarker] / (InvCoeff.CD[iMarker] + EPS);
            InvCoeff.CT[iMarker] = -InvCoeff.CFz[iMarker];
            InvCoeff.CQ[iMarker] = -InvCoeff.CMz[iMarker];
            InvCoeff.CMerit[iMarker] = InvCoeff.CT[iMarker] / (InvCoeff.CQ[iMarker] + EPS);
          }
          END_SU2_OMP_SINGLE
        }

        /*--- Update the coefficients for all boundaries. ---*/
        SU2_OMP_SINGLE
        {
          AllBoundInvCoeff.CD += InvCoeff.CD[iMarker];
          AllBoundInvCoeff.CL += InvCoeff.CL[iMarker];
          AllBoundInvCoeff.CSF += InvCoeff.CSF[iMarker];
          AllBoundInvCoeff.CEff = AllBoundInvCoeff.CL / (AllBoundInvCoeff.CD + EPS);
          AllBoundInvCoeff.CMx += InvCoeff.CMx[iMarker];
          AllBoundInvCoeff.CMy += InvCoeff.CMy[iMarker];
          AllBoundInvCoeff.CMz += InvCoeff.CMz[iMarker];
          AllBoundInvCoeff.CoPx += InvCoeff.CoPx[iMarker];
          AllBoundInvCoeff.CoPy += InvCoeff.CoPy[iMarker];
          AllBoundInvCoeff.CoPz += InvCoeff.CoPz[iMarker];
          AllBoundInvCoeff.CFx += InvCoeff.CFx[iMarker];
          AllBoundInvCoeff.CFy += InvCoeff.CFy[iMarker];
          AllBoundInvCoeff.CFz += InvCoeff.CFz[iMarker];
          AllBoundInvCoeff.CT += InvCoeff.CT[iMarker];
          AllBoundInvCoeff.CQ += InvCoeff.CQ[iMarker];
          AllBoundInvCoeff.CMerit = AllBoundInvCoeff.CT / (AllBoundInvCoeff.CQ + EPS);

          /*--- Compute the coefficients per surface ---*/
          for(unsigned short iMarker_Monitoring = 0;
               iMarker_Monitoring < config->GetnMarker_Monitoring();
               ++iMarker_Monitoring) {
            string Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
            string Marker_Tag = config->GetMarker_All_TagBound(iMarker);

            if(Marker_Tag == Monitoring_Tag) {
              SurfaceInvCoeff.CL[iMarker_Monitoring] += InvCoeff.CL[iMarker];
              SurfaceInvCoeff.CD[iMarker_Monitoring] += InvCoeff.CD[iMarker];
              SurfaceInvCoeff.CSF[iMarker_Monitoring] += InvCoeff.CSF[iMarker];
              SurfaceInvCoeff.CEff[iMarker_Monitoring] = SurfaceInvCoeff.CL[iMarker_Monitoring] / (SurfaceInvCoeff.CD[iMarker_Monitoring] + EPS);
              SurfaceInvCoeff.CFx[iMarker_Monitoring] += InvCoeff.CFx[iMarker];
              SurfaceInvCoeff.CFy[iMarker_Monitoring] += InvCoeff.CFy[iMarker];
              SurfaceInvCoeff.CFz[iMarker_Monitoring] += InvCoeff.CFz[iMarker];
              SurfaceInvCoeff.CMx[iMarker_Monitoring] += InvCoeff.CMx[iMarker];
              SurfaceInvCoeff.CMy[iMarker_Monitoring] += InvCoeff.CMy[iMarker];
              SurfaceInvCoeff.CMz[iMarker_Monitoring] += InvCoeff.CMz[iMarker];
            }
          }
        }
        END_SU2_OMP_SINGLE
      }
    }
  }

  SU2_OMP_SINGLE
  {
#ifdef HAVE_MPI
    /*--- Accumulate the boundary contribution from all nodes. ---*/
    if (config->GetComm_Level() == COMM_FULL) {

      /* Lambda function to carry out the reduction. */
      auto Allreduce = [](su2double x) {
        su2double tmp = x;
        x = 0.0;
        SU2_MPI::Allreduce(&tmp, &x, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
        return x;
      };

      /* The reduction of all elements of AllBoundInvCoeff. */
      AllBoundInvCoeff.CD = Allreduce(AllBoundInvCoeff.CD);
      AllBoundInvCoeff.CL = Allreduce(AllBoundInvCoeff.CL);
      AllBoundInvCoeff.CSF = Allreduce(AllBoundInvCoeff.CSF);
      AllBoundInvCoeff.CEff = AllBoundInvCoeff.CL / (AllBoundInvCoeff.CD + EPS);

      AllBoundInvCoeff.CMx = Allreduce(AllBoundInvCoeff.CMx);
      AllBoundInvCoeff.CMy = Allreduce(AllBoundInvCoeff.CMy);
      AllBoundInvCoeff.CMz = Allreduce(AllBoundInvCoeff.CMz);

      AllBoundInvCoeff.CoPx = Allreduce(AllBoundInvCoeff.CoPx);
      AllBoundInvCoeff.CoPy = Allreduce(AllBoundInvCoeff.CoPy);
      AllBoundInvCoeff.CoPz = Allreduce(AllBoundInvCoeff.CoPz);

      AllBoundInvCoeff.CFx = Allreduce(AllBoundInvCoeff.CFx);
      AllBoundInvCoeff.CFy = Allreduce(AllBoundInvCoeff.CFy);
      AllBoundInvCoeff.CFz = Allreduce(AllBoundInvCoeff.CFz);

      AllBoundInvCoeff.CT = Allreduce(AllBoundInvCoeff.CT);
      AllBoundInvCoeff.CQ = Allreduce(AllBoundInvCoeff.CQ);
      AllBoundInvCoeff.CMerit = AllBoundInvCoeff.CT / (AllBoundInvCoeff.CQ + EPS);

      /*--- Lambda function for the reduction of multiple entities. ---*/
      int nMarkerMon = config->GetnMarker_Monitoring();
      su2double* buffer = new su2double[nMarkerMon];

      auto Allreduce_inplace = [buffer](int size, su2double* x) {
        SU2_MPI::Allreduce(x, buffer, size, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
        for (int i = 0; i < size; ++i) x[i] = buffer[i];
      };

      /*--- The reduction of the forces and moments of the surfaces. ---*/
      Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CL);
      Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CD);
      Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CSF);

      for (int iMarker_Monitoring = 0; iMarker_Monitoring < nMarkerMon; ++iMarker_Monitoring)
        SurfaceInvCoeff.CEff[iMarker_Monitoring] =
            SurfaceInvCoeff.CL[iMarker_Monitoring] / (SurfaceInvCoeff.CD[iMarker_Monitoring] + EPS);

      Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CFx);
      Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CFy);
      Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CFz);

      Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CMx);
      Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CMy);
      Allreduce_inplace(nMarkerMon, SurfaceInvCoeff.CMz);

      delete[] buffer;
    }
#endif

    /*--- Update the total coefficients (note that all the MPI ranks have the same value) ---*/
    TotalCoeff.CD = AllBoundInvCoeff.CD;
    TotalCoeff.CL = AllBoundInvCoeff.CL;
    TotalCoeff.CSF = AllBoundInvCoeff.CSF;
    TotalCoeff.CEff = TotalCoeff.CL / (TotalCoeff.CD + EPS);
    TotalCoeff.CFx = AllBoundInvCoeff.CFx;
    TotalCoeff.CFy = AllBoundInvCoeff.CFy;
    TotalCoeff.CFz = AllBoundInvCoeff.CFz;
    TotalCoeff.CMx = AllBoundInvCoeff.CMx;
    TotalCoeff.CMy = AllBoundInvCoeff.CMy;
    TotalCoeff.CMz = AllBoundInvCoeff.CMz;
    TotalCoeff.CoPx = AllBoundInvCoeff.CoPx;
    TotalCoeff.CoPy = AllBoundInvCoeff.CoPy;
    TotalCoeff.CoPz = AllBoundInvCoeff.CoPz;
    TotalCoeff.CT = AllBoundInvCoeff.CT;
    TotalCoeff.CQ = AllBoundInvCoeff.CQ;
    TotalCoeff.CMerit = TotalCoeff.CT / (TotalCoeff.CQ + EPS);

    /*--- And the same for the surfaces. ---*/
    for (int iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
      SurfaceCoeff.CL[iMarker_Monitoring] = SurfaceInvCoeff.CL[iMarker_Monitoring];
      SurfaceCoeff.CD[iMarker_Monitoring] = SurfaceInvCoeff.CD[iMarker_Monitoring];
      SurfaceCoeff.CSF[iMarker_Monitoring] = SurfaceInvCoeff.CSF[iMarker_Monitoring];
      SurfaceCoeff.CEff[iMarker_Monitoring] =
          SurfaceCoeff.CL[iMarker_Monitoring] / (SurfaceCoeff.CD[iMarker_Monitoring] + EPS);
      SurfaceCoeff.CFx[iMarker_Monitoring] = SurfaceInvCoeff.CFx[iMarker_Monitoring];
      SurfaceCoeff.CFy[iMarker_Monitoring] = SurfaceInvCoeff.CFy[iMarker_Monitoring];
      SurfaceCoeff.CFz[iMarker_Monitoring] = SurfaceInvCoeff.CFz[iMarker_Monitoring];
      SurfaceCoeff.CMx[iMarker_Monitoring] = SurfaceInvCoeff.CMx[iMarker_Monitoring];
      SurfaceCoeff.CMy[iMarker_Monitoring] = SurfaceInvCoeff.CMy[iMarker_Monitoring];
      SurfaceCoeff.CMz[iMarker_Monitoring] = SurfaceInvCoeff.CMz[iMarker_Monitoring];
    }
  }
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                               CConfig *config, unsigned short iRKStep) {

  /*--- Get the required parameters for the Runge-Kutta time integration. ---*/
  const su2double      RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  const unsigned short nRKStages     = config->GetnRKStep();

  /*--- Loop over owned elements. ---*/
#ifdef HAVE_OMP
  const size_t omp_chunk_size_elem = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif
  SU2_OMP_FOR_STAT(omp_chunk_size_elem)
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /*--- Set the reference for the variable where the new solution must be stored. ---*/
    ColMajorMatrix<su2double> &solNew = (iRKStep == (nRKStages-1)) ? volElem[l].solDOFs :
                                                                     volElem[l].solDOFsWork;

    /*--- Determine the update coefficient for this time step. ---*/
    const su2double tmp = RK_AlphaCoeff*volElem[l].deltaTime;

    /*--- Compute the new state. ---*/
    const unsigned short nDOFs = volElem[l].standardElemFlow->GetNDOFs();

    for(unsigned short j=0; j<nVar; ++j) {
      SU2_OMP_SIMD
      for(unsigned short i=0; i<nDOFs; ++i)
        solNew(i,j) = volElem[l].solDOFs(i,j) - tmp*volElem[l].resDOFs(i,j);
    }
  }
  END_SU2_OMP_FOR

  /*--- Test for the last RK step. ---*/
  if(iRKStep == (nRKStages-1)) {

    /*--- Compute the root mean square residual. Note that the SetResidual_RMS
          function of CSolver cannot be used, because that is for the FV solver. ---*/
    SetResidual_RMS_FEM(geometry, config);

    /*--- For verification cases, compute the global error metrics. ---*/
    ComputeVerificationError(geometry, config);
  }
}

void CFEM_DG_EulerSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                                 CConfig *config, unsigned short iRKStep) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::SetResidual_RMS_FEM(CGeometry *geometry,
                                              CConfig *config) {

  SU2_OMP_SINGLE
  {
#ifdef HAVE_MPI
    /* Parallel mode. Disable the reduce for the residual to avoid overhead if requested. */
    if (config->GetComm_Level() == COMM_FULL) {

      /*--- The local L2 norms must be added to obtain the
            global value. Also check for divergence. ---*/
      vector<su2double> rbufRes(nVar);
      SU2_MPI::Allreduce(Residual_RMS.data(), rbufRes.data(), nVar, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

      for(unsigned short iVar=0; iVar<nVar; ++iVar) {
        if (rbufRes[iVar] != rbufRes[iVar])
          SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);

        Residual_RMS[iVar] = max(EPS*EPS, sqrt(rbufRes[iVar]/nDOFsGlobal));
      }

      /*--- The global maximum norms must be obtained. ---*/
      rbufRes.resize(nVar*size);
      SU2_MPI::Allgather(Residual_Max.data(), nVar, MPI_DOUBLE, rbufRes.data(),
                         nVar, MPI_DOUBLE, SU2_MPI::GetComm());

      vector<unsigned long> rbufPoint(nVar*size);
      SU2_MPI::Allgather(Point_Max.data(), nVar, MPI_UNSIGNED_LONG, rbufPoint.data(),
                         nVar, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

      vector<su2double> sbufCoor(nDim*nVar);
      for(unsigned short iVar=0; iVar<nVar; ++iVar) {
        for(unsigned short iDim=0; iDim<nDim; ++iDim)
          sbufCoor[iVar*nDim+iDim] = Point_Max_Coord[iVar][iDim];
      }

      vector<su2double> rbufCoor(nDim*nVar*size);
      SU2_MPI::Allgather(sbufCoor.data(), nVar*nDim, MPI_DOUBLE, rbufCoor.data(),
                         nVar*nDim, MPI_DOUBLE, SU2_MPI::GetComm());

      for(unsigned short iVar=0; iVar<nVar; ++iVar) {
        for(int proc=0; proc<size; ++proc)
          AddRes_Max(iVar, rbufRes[proc*nVar+iVar], rbufPoint[proc*nVar+iVar],
                     &rbufCoor[proc*nVar*nDim+iVar*nDim]);
      }
    }

#else
    /*--- Sequential mode. Check for a divergence of the solver and compute
          the L2-norm of the residuals. ---*/
    for(unsigned short iVar=0; iVar<nVar; ++iVar) {

      if(GetRes_RMS(iVar) != GetRes_RMS(iVar))
        SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);

      Residual_RMS[iVar] = max(EPS*EPS, sqrt(GetRes_RMS(iVar)/nDOFsGlobal));
    }
#endif
  }
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::ComputeVerificationError(CGeometry *geometry,
                                                   CConfig   *config) {

  /*--- The errors only need to be computed on the finest grid. ---*/
  if(MGLevel != MESH_0) return;

  /*--- If this is a verification case, we can compute the global
   error metrics by using the difference between the local error
   and the known solution at each DOF. This is then collected into
   RMS (L2) and maximum (Linf) global error norms. From these
   global measures, one can compute the order of accuracy. ---*/

  bool write_heads = ((((config->GetTimeIter() % (config->GetScreen_Wrt_Freq(2)*40)) == 0)
                       && (config->GetTimeIter()!= 0))
                      || (config->GetTimeIter() == 1));
  if( !write_heads ) return;

  /*--- Check if there actually is an exact solution for this
        verification case, if computed at all. ---*/
  if (VerificationSolution) {
    if (VerificationSolution->ExactSolutionKnown()) {

      /*--- Get the physical time if necessary. ---*/
      su2double time = 0.0;
      if (config->GetTime_Marching() != TIME_MARCHING::STEADY) time = config->GetPhysicalTime();

      /*--- Reset the global error measures to zero. ---*/
      SU2_OMP_SINGLE
      {
        for (unsigned short iVar = 0; iVar < nVar; iVar++) {
          VerificationSolution->SetError_RMS(iVar, 0.0);
          VerificationSolution->SetError_Max(iVar, 0.0, 0);
        }
      }
      END_SU2_OMP_SINGLE


      for(int i=0; i<size; ++i) {

        if(i == rank) {

          const int thread = omp_get_thread_num();
          for(int j=0; j<omp_get_num_threads(); ++j) {
            if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
            SU2_OMP_BARRIER
          }
        }

        SU2_OMP_SINGLE
        SU2_MPI::Barrier(SU2_MPI::GetComm());
        END_SU2_OMP_SINGLE
      }

      SU2_OMP_SINGLE
      SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
      END_SU2_OMP_SINGLE
    }
  }
}

void CFEM_DG_EulerSolver::ADER_DG_Iteration(const unsigned long elemBeg,
                                            const unsigned long elemEnd) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::BoundaryStates_Euler_Wall(const CSurfaceElementFEM        *surfElem,
                                                    const ColMajorMatrix<su2double> &solIntL,
                                                    ColMajorMatrix<su2double>       &solIntR,
                                                    const su2double                 factNorm) {

  /*--- Easier storage of the number of items for which the BC must be applied.
        Usually the (padded) number of integration points. ---*/
  const unsigned short nItems = solIntR.rows();

  /*--- Make a distinction between 2D and 3D for optimal performance. ---*/
  switch( nDim ) {
    case 2: {
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nItems; ++i) {

        /*--- Compute the normal velocity relative to the prescribed grid velocity. ---*/
        const su2double nx = surfElem->metricNormalsFace(i,0);
        const su2double ny = surfElem->metricNormalsFace(i,1);

        const su2double Vn = (solIntL(i,1) - surfElem->gridVelocities(i,0))*nx
                           + (solIntL(i,2) - surfElem->gridVelocities(i,1))*ny;

        /*--- Compute the boundary state by removing or mirroring the normal component,
              factNorm is 1.0 is removing, factNorm is 2.0 is mirroring. ---*/
        solIntR(i,0) = solIntL(i,0);
        solIntR(i,1) = solIntL(i,1) - factNorm*Vn*nx;
        solIntR(i,2) = solIntL(i,2) - factNorm*Vn*ny;
        solIntR(i,3) = solIntL(i,3);
      }
      break;
    }

    case 3: {
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nItems; ++i) {

        /*--- Compute the normal velocity relative to the prescribed grid velocity. ---*/
        const su2double nx = surfElem->metricNormalsFace(i,0);
        const su2double ny = surfElem->metricNormalsFace(i,1);
        const su2double nz = surfElem->metricNormalsFace(i,2);

        const su2double Vn = (solIntL(i,1) - surfElem->gridVelocities(i,0))*nx
                           + (solIntL(i,2) - surfElem->gridVelocities(i,1))*ny
                           + (solIntL(i,3) - surfElem->gridVelocities(i,2))*nz;

        /*--- Compute the boundary state by removing or mirroring the normal component,
              factNorm is 1.0 is removing, factNorm is 2.0 is mirroring. ---*/
        solIntR(i,0) = solIntL(i,0);
        solIntR(i,1) = solIntL(i,1) - factNorm*Vn*nx;
        solIntR(i,2) = solIntL(i,2) - factNorm*Vn*ny;
        solIntR(i,3) = solIntL(i,3) - factNorm*Vn*nz;
        solIntR(i,4) = solIntL(i,4);
      }
      break;
    }
  }
}

void CFEM_DG_EulerSolver::BoundaryStates_Inlet(CConfig                  *config,
                                               const unsigned short     nFaceSimul,
                                               const unsigned short     NPad,
                                               const CSurfaceElementFEM *surfElem,
                                               unsigned short           val_marker,
                                               const su2double          *solIntL,
                                               su2double                *solIntR) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::BoundaryStates_Outlet(CConfig                  *config,
                                                const unsigned short     nFaceSimul,
                                                const unsigned short     NPad,
                                                const CSurfaceElementFEM *surfElem,
                                                unsigned short           val_marker,
                                                const su2double          *solIntL,
                                                su2double                *solIntR) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::BoundaryStates_Riemann(CConfig                  *config,
                                                 const unsigned short     nFaceSimul,
                                                 const unsigned short     NPad,
                                                 const CSurfaceElementFEM *surfElem,
                                                 unsigned short           val_marker,
                                                 const su2double          *solIntL,
                                                 su2double                *solIntR) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::BC_Euler_Wall(CConfig             *config,
                                        const unsigned long surfElemBeg,
                                        const unsigned long surfElemEnd,
                                        CSurfaceElementFEM  *surfElem,
                                        CNumerics           *conv_numerics) {

  /*--- Determine the chunk size for the OpenMP parallelization. ---*/
#ifdef HAVE_OMP
    const unsigned long nFaces = surfElemEnd - surfElemBeg;
    const size_t omp_chunk_size = computeStaticChunkSize(nFaces, omp_get_num_threads(), 64);
#endif

  /*--- Determine the index in the work arrays where the right solution must be stored. ---*/
  const unsigned int indRight = omp_get_num_threads() + omp_get_thread_num();

  /*--- Loop over the requested range of surface faces. ---*/
  SU2_OMP_FOR_DYN(omp_chunk_size)
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /*--- Compute the variables of the left state in the integration point of the face
          and convert them to primitive variables. ---*/
    ColMajorMatrix<su2double> &solIntLeft = surfElem[l].ComputeSolSide0IntPoints(volElem);
    EntropyToPrimitiveVariables(solIntLeft);

    /*--- Compute the right state by applying the inviscid wall BC's. ---*/
    ColMajorMatrix<su2double> &solIntRight = surfElem[l].standardElemFlow->workSolInt[indRight];
    BoundaryStates_Euler_Wall(&surfElem[l], solIntLeft, solIntRight, 1.0);

    /*--- The remainder of the contribution of this boundary face to the residual
          is the same for all boundary conditions. Hence a generic function can
          be used to carry out this task. ---*/
    ResidualInviscidBoundaryFace(config, conv_numerics, &surfElem[l], solIntLeft, solIntRight);
  }
  END_SU2_OMP_FOR
}

void CFEM_DG_EulerSolver::BC_Far_Field(CConfig             *config,
                                       const unsigned long surfElemBeg,
                                       const unsigned long surfElemEnd,
                                       CSurfaceElementFEM  *surfElem,
                                       CNumerics           *conv_numerics) {

  /*--- Determine the chunk size for the OpenMP parallelization. ---*/
#ifdef HAVE_OMP
    const unsigned long nFaces = surfElemEnd - surfElemBeg;
    const size_t omp_chunk_size = computeStaticChunkSize(nFaces, omp_get_num_threads(), 64);
#endif

  /*--- Determine the index in the work arrays where the left solution must be stored. ---*/
  const unsigned int indRight = omp_get_num_threads() + omp_get_thread_num();

  /*--- Loop over the requested range of surface faces. ---*/
  SU2_OMP_FOR_DYN(omp_chunk_size)
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /*--- Compute the variables of the left state in the integration point of the face. ---*/
    ColMajorMatrix<su2double> &solIntLeft = surfElem[l].ComputeSolSide0IntPoints(volElem);

    /*--- Set the right state in the integration points to the free stream value. ---*/
    ColMajorMatrix<su2double> &solIntRight = surfElem[l].standardElemFlow->workSolInt[indRight];

    const unsigned short nRows = solIntRight.rows();
    for(unsigned short j=0; j<nVar; ++j) {
      SU2_OMP_SIMD
      for(unsigned short i=0; i<nRows; ++i)
        solIntRight(i,j) = EntropyVarFreeStream[j];
    }

    /*--- Convert the entropy variables to the primitive variables. ---*/
    EntropyToPrimitiveVariables(solIntLeft);
    EntropyToPrimitiveVariables(solIntRight);

    /*--- The remainder of the contribution of this boundary face to the residual
          is the same for all boundary conditions. Hence a generic function can
          be used to carry out this task. ---*/
    ResidualInviscidBoundaryFace(config, conv_numerics, &surfElem[l], solIntLeft, solIntRight);
  }
  END_SU2_OMP_FOR
}

void CFEM_DG_EulerSolver::BC_Sym_Plane(CConfig             *config,
                                       const unsigned long surfElemBeg,
                                       const unsigned long surfElemEnd,
                                       CSurfaceElementFEM  *surfElem,
                                       CNumerics           *conv_numerics) {

  /*--- Determine the chunk size for the OpenMP parallelization. ---*/
#ifdef HAVE_OMP
    const unsigned long nFaces = surfElemEnd - surfElemBeg;
    const size_t omp_chunk_size = computeStaticChunkSize(nFaces, omp_get_num_threads(), 64);
#endif

  /*--- Determine the index in the work arrays where the right solution must be stored. ---*/
  const unsigned int indRight = omp_get_num_threads() + omp_get_thread_num();

  /*--- Loop over the requested range of surface faces. ---*/
  SU2_OMP_FOR_DYN(omp_chunk_size)
  for(unsigned long l=surfElemBeg; l<surfElemEnd; ++l) {

    /*--- Compute the variables of the left state in the integration point of the face
          and convert them to primitive variables. ---*/
    ColMajorMatrix<su2double> &solIntLeft = surfElem[l].ComputeSolSide0IntPoints(volElem);
    EntropyToPrimitiveVariables(solIntLeft);

    /*--- Compute the right state by applying the mirror BC's. This is almost the same
          as the Euler wall, except that the normal velocity is mirrored. ---*/
    ColMajorMatrix<su2double> &solIntRight = surfElem[l].standardElemFlow->workSolInt[indRight];
    BoundaryStates_Euler_Wall(&surfElem[l], solIntLeft, solIntRight, 2.0);

    /*--- The remainder of the contribution of this boundary face to the residual
          is the same for all boundary conditions. Hence a generic function can
          be used to carry out this task. ---*/
    ResidualInviscidBoundaryFace(config, conv_numerics, &surfElem[l], solIntLeft, solIntRight);
  }
  END_SU2_OMP_FOR
}

void CFEM_DG_EulerSolver::BC_Supersonic_Outlet(CConfig             *config,
                                               const unsigned long surfElemBeg,
                                               const unsigned long surfElemEnd,
                                               CSurfaceElementFEM  *surfElem,
                                               CNumerics           *conv_numerics) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::BC_Inlet(CConfig             *config,
                                   const unsigned long surfElemBeg,
                                   const unsigned long surfElemEnd,
                                   CSurfaceElementFEM  *surfElem,
                                   CNumerics           *conv_numerics,
                                   unsigned short      val_marker) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::BC_Outlet(CConfig             *config,
                                    const unsigned long surfElemBeg,
                                    const unsigned long surfElemEnd,
                                    CSurfaceElementFEM  *surfElem,
                                    CNumerics           *conv_numerics,
                                    unsigned short      val_marker) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::BC_HeatFlux_Wall(CConfig            *config,
                                           const unsigned long surfElemBeg,
                                           const unsigned long surfElemEnd,
                                           CSurfaceElementFEM *surfElem,
                                           CNumerics           *conv_numerics,
                                           unsigned short      val_marker) {
  /*--- This is a viscous wall boundary condition, which cannot be applied for
        an inviscid analysis. Instead an inviscid wall boundary condition is
        applied. ---*/
  BC_Euler_Wall(config, surfElemBeg, surfElemEnd, surfElem, conv_numerics);
}

void CFEM_DG_EulerSolver::BC_Isothermal_Wall(CConfig            *config,
                                             const unsigned long surfElemBeg,
                                             const unsigned long surfElemEnd,
                                             CSurfaceElementFEM *surfElem,
                                             CNumerics           *conv_numerics,
                                             unsigned short      val_marker) {
  /*--- This is a viscous wall boundary condition, which cannot be applied for
        an inviscid analysis. Instead an inviscid wall boundary condition is
        applied. ---*/
  BC_Euler_Wall(config, surfElemBeg, surfElemEnd, surfElem, conv_numerics);
}

void CFEM_DG_EulerSolver::BC_Riemann(CConfig             *config,
                                     const unsigned long surfElemBeg,
                                     const unsigned long surfElemEnd,
                                     CSurfaceElementFEM  *surfElem,
                                     CNumerics           *conv_numerics,
                                     unsigned short      val_marker) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::BC_Custom(CConfig             *config,
                                    const unsigned long surfElemBeg,
                                    const unsigned long surfElemEnd,
                                    CSurfaceElementFEM  *surfElem,
                                    CNumerics           *conv_numerics) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::ResidualInviscidBoundaryFace(
                                      CConfig                  *config,
                                      CNumerics                 *conv_numerics,
                                      CSurfaceElementFEM        *surfElem,
                                      ColMajorMatrix<su2double> &solInt0,
                                      ColMajorMatrix<su2double> &solInt1) {

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Compute the fluxes in the integration points using the     ---*/
  /*---         approximate Riemann solver.                                ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Abbreviate the number of padded integration points
        and the integration weights. ---*/
  const unsigned short nIntPad = surfElem->standardElemFlow->GetNIntegrationPad();
  const passivedouble *weights = surfElem->standardElemFlow->GetIntegrationWeights();

  /*--- Compute the invisid fluxes in the integration points of the face.
        Note that solInt0 is used to store the fluxes. ---*/
  ComputeInviscidFluxesFace(config, solInt0, solInt1, surfElem->JacobiansFace,
                            surfElem->metricNormalsFace, surfElem->gridVelocities,
                            conv_numerics, solInt0);

  /*--- Multiply the fluxes with the integration weight of the
        corresponding integration point. ---*/
  for(unsigned short j=0; j<nVar; ++j) {
    SU2_OMP_SIMD_IF_NOT_AD
    for(unsigned short i=0; i<nIntPad; ++i)
      solInt0(i,j) *= weights[i];
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Compute the contribution to the residuals of the DOFs of   ---*/
  /*---         the elements on the left and right side of the face.       ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Initialize the residual to zero and add the contribution
        from the fluxes. ---*/
  surfElem->resDOFsElem.setConstant(0.0);
  surfElem->ResidualBasisFunctions(solInt0);
}

void CFEM_DG_EulerSolver::LeftStatesIntegrationPointsBoundaryFace(
                                             CConfig                  *config,
                                             const unsigned short     nFaceSimul,
                                             const unsigned short     NPad,
                                             const CSurfaceElementFEM *surfElem,
                                             su2double                *solFace,
                                             su2double                *solIntL) {

  for(int i=0; i<size; ++i) {

    if(i == rank) {

      const int thread = omp_get_thread_num();
      for(int j=0; j<omp_get_num_threads(); ++j) {
        if(j == thread) cout << "Rank: " << i << ", thread: " << j << endl << flush;
        SU2_OMP_BARRIER
      }
    }

    SU2_OMP_SINGLE
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    END_SU2_OMP_SINGLE
  }

  SU2_OMP_SINGLE
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
  END_SU2_OMP_SINGLE
}

void CFEM_DG_EulerSolver::ComputeInviscidFluxesFace(CConfig                   *config,
                                                    ColMajorMatrix<su2double> &solLeft,
                                                    ColMajorMatrix<su2double> &solRight,
                                                    su2activevector           &JacobiansFace,
                                                    ColMajorMatrix<su2double> &normalsFace,
                                                    ColMajorMatrix<su2double> &gridVelocities,
                                                    CNumerics                 *numerics,
                                                    ColMajorMatrix<su2double> &fluxes) {

  /*--- Some abbreviations for terms involving the specific heat ratio and
        the number of items for which the flux must be computed. ---*/
  const su2double gm1 = Gamma_Minus_One;
  const su2double ovgm1 = 1.0/gm1;
  const unsigned int nItems = gridVelocities.rows();

  /*--- Make a distinction between the several Riemann solvers. ---*/
  switch( config->GetRiemann_Solver_FEM() ) {

    case UPWIND::ISMAIL_ROE: {

      /*--- Entropy satisfying Riemann flux of Ismail and Roe.
            Values to scale the acoustic eigenvalues to obtain an adequate amount
            of dissipation to be entropy satisfying in the Ismail_Roe flux. Note
            that only beta is taken, assuming that the jump in Mach number over the
            interface is less than 0.5. For alphaMax = 0 the EC1 flux is obtained. ---*/
      const su2double beta     = 1.0/6.0;
      //const su2double alphaMax = 2.0;
      const su2double alphaMax = 0.0;

      /*--- Some more constants involving Gamma. ---*/
      const su2double gp1Ovg = (Gamma+1.0)/Gamma;
      const su2double gm1Ovg = gm1/Gamma;

      /*--- Make a distinction between two and three space dimensions
            in order to have the most efficient code. ---*/
      switch( nDim ) {

        case 2: {

          /*--- Loop over the items. ---*/
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned int i=0; i<nItems; ++i) {

            /*--- Easier storage of the components of the unit normal and
                  the Jacobian (which is the area of the face). ---*/
            const su2double nx = normalsFace(i,0);
            const su2double ny = normalsFace(i,1);
            const su2double Area = JacobiansFace(i);

            /*--- Compute the normal grid velocity. ---*/
            const su2double gridVelNorm = gridVelocities(i,0)*nx + gridVelocities(i,1)*ny;

            /*--- Easier storage of the primitive variables of
                  the left and right state. ---*/
            const su2double rhoL = solLeft(i,0), rhoR = solRight(i,0);
            const su2double uL   = solLeft(i,1), uR   = solRight(i,1);
            const su2double vL   = solLeft(i,2), vR   = solRight(i,2);
            const su2double pL   = solLeft(i,3), pR   = solRight(i,3);

            /*--- Compute the total energy of the left and right state. ---*/
            const su2double rEL = ovgm1*pL + 0.5*rhoL*(uL*uL + vL*vL);
            const su2double rER = ovgm1*pR + 0.5*rhoR*(uR*uR + vR*vR);

            /*--- Compute the entropy variables of the left and right state. ---*/
            const su2double sL    = log(pL/pow(rhoL,Gamma));
            const su2double pInvL = 1.0/pL;

            const su2double V3L = -pInvL*rhoL;
            const su2double V2L = -V3L*vL;
            const su2double V1L = -V3L*uL;
            const su2double V0L = (Gamma-sL)*ovgm1 + 0.5*V3L*(uL*uL + vL*vL);

            const su2double sR    = log(pR/pow(rhoR,Gamma));
            const su2double pInvR = 1.0/pR;

            const su2double V3R = -pInvR*rhoR;
            const su2double V2R = -V3R*vR;
            const su2double V1R = -V3R*uR;
            const su2double V0R = (Gamma-sR)*ovgm1 + 0.5*V3R*(uR*uR + vR*vR);

            /*--- Compute the z-variables of the left and right states. ---*/
            const su2double z0L = sqrt(rhoL/pL), z0R = sqrt(rhoR/pR);
            const su2double z1L = uL*z0L,        z1R = uR*z0R;
            const su2double z2L = vL*z0L,        z2R = vR*z0R;
            const su2double z3L = sqrt(rhoL*pL), z3R = sqrt(rhoR*pR);

            /*--- Compute the arithmetic average of the z-variables. ---*/
            const su2double z0Avg = 0.5*(z0L + z0R);
            const su2double z1Avg = 0.5*(z1L + z1R);
            const su2double z2Avg = 0.5*(z2L + z2R);
            const su2double z3Avg = 0.5*(z3L + z3R);

            /*--- Compute the logarithmic mean of z0. ---*/
            su2double zeta, f, u, F;

            zeta = z0L/z0R;
            f    = (zeta-1.0)/(zeta+1.0);
            u    = f*f;
            if(u < (su2double) 0.01)
              F = 1.0 + u/3.0 + u*u/5.0 + u*u*u/7.0 + u*u*u*u/9.0;
            else
              F = log(zeta)/(2.0*f);

            const su2double z0LogAvg = z0Avg/F;

            /*--- Compute the logarithmic mean of z3. ---*/
            zeta = z3L/z3R;
            f    = (zeta-1.0)/(zeta+1.0);
            u    = f*f;
            if(u < (su2double) 0.01)
              F = 1.0 + u/3.0 + u*u/5.0 + u*u*u/7.0 + u*u*u*u/9.0;
            else
              F = log(zeta)/(2.0*f);

            const su2double z3LogAvg = z3Avg/F;

            /*--- Compute the other averaged quantities that are necessary. ---*/
            const su2double oneOvz0Avg = 1.0/z0Avg;
            const su2double rhoAvg = z0Avg*z3LogAvg;
            const su2double p1Avg  = oneOvz0Avg*z3Avg;
            const su2double p2Avg  = 0.5*(gp1Ovg*z3LogAvg/z0LogAvg + gm1Ovg*p1Avg);
            const su2double uAvg   = oneOvz0Avg*z1Avg;
            const su2double vAvg   = oneOvz0Avg*z2Avg;

            const su2double vnAvg  = uAvg*nx + vAvg*ny;
            const su2double unAvg  = vnAvg - gridVelNorm;
            const su2double kinAvg = 0.5*(uAvg*uAvg + vAvg*vAvg);
            const su2double a2Avg  = Gamma*p2Avg/rhoAvg;
            const su2double aAvg   = sqrt(a2Avg);
            const su2double HAvg   = a2Avg*ovgm1 + kinAvg;
            const su2double EAvg   = HAvg - p2Avg/rhoAvg;

            const su2double ovaAvg  = 1.0/aAvg;
            const su2double ova2Avg = 1.0/a2Avg;

            /*--- Compute the difference in entropy variables. ---*/
            const su2double dV0 = V0R - V0L, dV1 = V1R - V1L, dV2 = V2R - V2L,
                            dV3 = V3R - V3L;

            /*--- Define the difference in conservative variables as dU/dV deltaV, where the
                  transformation matrix dU/dV must be evaluated at the averaged state. ---*/
            su2double abv1 = rhoAvg*(uAvg*dV1 + vAvg*dV2);
            su2double abv2 = abv1 + rhoAvg*(dV0 + HAvg*dV3);

            const su2double dr  = abv1 + rhoAvg*(dV0 + EAvg*dV3);
            const su2double dru = uAvg*abv2 + p2Avg*dV1;
            const su2double drv = vAvg*abv2 + p2Avg*dV2;
            const su2double drE = HAvg*abv1 + rhoAvg*EAvg*(dV0 + HAvg*dV3)
                                + p2Avg*kinAvg*dV3;

            /*--- Compute the absolute values of the eigenvalues of the flux Jacobian. ---*/
            su2double lam1 = fabs(unAvg + aAvg);
            su2double lam2 = fabs(unAvg - aAvg);
            su2double lam3 = fabs(unAvg);

            /*--- Scale the acoustic eigenvalue, such that the EC2 (or EC1) flux of Ismail
                  and Roe is obtained. Also multiply all eigenvalues by half to obtain
                  the correct scaling. ---*/
            lam1 *= 0.5*(1.0 + beta + alphaMax);
            lam2 *= 0.5*(1.0 + beta + alphaMax);
            lam3 *= 0.5;

            /*--- Some abbreviations, which occur quite often in the dissipation terms. ---*/
            abv1 = 0.5*(lam1 + lam2);
            abv2 = 0.5*(lam1 - lam2);

            const su2double abv3 = abv1 - lam3;
            const su2double abv4 = gm1*(kinAvg*dr - uAvg*dru - vAvg*drv + drE);
            const su2double abv5 = nx*dru + ny*drv - vnAvg*dr;
            const su2double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
            const su2double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

            /*--- Compute the fluxes. ---*/
            fluxes(i,0) = Area*(rhoAvg*unAvg - lam3*dr -abv6);
            fluxes(i,1) = Area*(rhoAvg*unAvg*uAvg + p1Avg*nx - lam3*dru - uAvg*abv6 - nx*abv7);
            fluxes(i,2) = Area*(rhoAvg*unAvg*vAvg + p1Avg*ny - lam3*drv - vAvg*abv6 - ny*abv7);
            fluxes(i,3) = Area*(rhoAvg*unAvg*HAvg + p2Avg*gridVelNorm  - lam3*drE - HAvg*abv6 - vnAvg*abv7);
          }

          break;
        }

        /*--------------------------------------------------------------------*/

        case 3: {

          /*--- Loop over the items. ---*/
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned int i=0; i<nItems; ++i) {

            /*--- Easier storage of the components of the unit normal and
                  the Jacobian (which is the area of the face). ---*/
            const su2double nx = normalsFace(i,0);
            const su2double ny = normalsFace(i,1);
            const su2double nz = normalsFace(i,2);
            const su2double Area = JacobiansFace(i);

            /*--- Compute the normal grid velocity. ---*/
            const su2double gridVelNorm = gridVelocities(i,0)*nx + gridVelocities(i,1)*ny
                                        + gridVelocities(i,2)*nz;

            /*--- Easier storage of the primitive variables of
                  the left and right state. ---*/
            const su2double rhoL = solLeft(i,0), rhoR = solRight(i,0);
            const su2double uL   = solLeft(i,1), uR   = solRight(i,1);
            const su2double vL   = solLeft(i,2), vR   = solRight(i,2);
            const su2double wL   = solLeft(i,3), wR   = solRight(i,3);
            const su2double pL   = solLeft(i,4), pR   = solRight(i,4);

            /*--- Compute the total energy of the left and right state. ---*/
            const su2double rEL = ovgm1*pL + 0.5*rhoL*(uL*uL + vL*vL + wL*wL);
            const su2double rER = ovgm1*pR + 0.5*rhoR*(uR*uR + vR*vR + wR*wR);

            /*--- Compute the entropy variables of the left and right state. ---*/
            const su2double sL    = log(pL/pow(rhoL,Gamma));
            const su2double pInvL = 1.0/pL;

            const su2double V4L = -pInvL*rhoL;
            const su2double V3L = -V4L*wL;
            const su2double V2L = -V4L*vL;
            const su2double V1L = -V4L*uL;
            const su2double V0L = (Gamma-sL)*ovgm1 + 0.5*V4L*(uL*uL + vL*vL + wL*wL);

            const su2double sR    = log(pR/pow(rhoR,Gamma));
            const su2double pInvR = 1.0/pR;

            const su2double V4R = -pInvR*rhoR;
            const su2double V3R = -V4R*wR;
            const su2double V2R = -V4R*vR;
            const su2double V1R = -V4R*uR;
            const su2double V0R = (Gamma-sR)*ovgm1 + 0.5*V4R*(uR*uR + vR*vR + wR*wR);

            /*--- Compute the z-variables of the left and right states. ---*/
            const su2double z0L = sqrt(rhoL/pL), z0R = sqrt(rhoR/pR);
            const su2double z1L = uL*z0L,        z1R = uR*z0R;
            const su2double z2L = vL*z0L,        z2R = vR*z0R;
            const su2double z3L = wL*z0L,        z3R = wR*z0R;
            const su2double z4L = sqrt(rhoL*pL), z4R = sqrt(rhoR*pR);

            /*--- Compute the arithmetic average of the z-variables. ---*/
            const su2double z0Avg = 0.5*(z0L + z0R);
            const su2double z1Avg = 0.5*(z1L + z1R);
            const su2double z2Avg = 0.5*(z2L + z2R);
            const su2double z3Avg = 0.5*(z3L + z3R);
            const su2double z4Avg = 0.5*(z4L + z4R);

            /*--- Compute the logarithmic mean of z0. ---*/
            su2double zeta, f, u, F;

            zeta = z0L/z0R;
            f    = (zeta-1.0)/(zeta+1.0);
            u    = f*f;
            if(u < (su2double) 0.01)
              F = 1.0 + u/3.0 + u*u/5.0 + u*u*u/7.0 + u*u*u*u/9.0;
            else
              F = log(zeta)/(2.0*f);

            const su2double z0LogAvg = z0Avg/F;

            /*--- Compute the logarithmic mean of z4. ---*/
            zeta = z4L/z4R;
            f    = (zeta-1.0)/(zeta+1.0);
            u    = f*f;
            if(u < (su2double) 0.01)
              F = 1.0 + u/3.0 + u*u/5.0 + u*u*u/7.0 + u*u*u*u/9.0;
            else
              F = log(zeta)/(2.0*f);

            const su2double z4LogAvg = z4Avg/F;

            /*--- Compute the other averaged quantities that are necessary. ---*/
            const su2double oneOvz0Avg = 1.0/z0Avg;
            const su2double rhoAvg = z0Avg*z4LogAvg;
            const su2double p1Avg  = oneOvz0Avg*z4Avg;
            const su2double p2Avg  = 0.5*(gp1Ovg*z4LogAvg/z0LogAvg + gm1Ovg*p1Avg);
            const su2double uAvg   = oneOvz0Avg*z1Avg;
            const su2double vAvg   = oneOvz0Avg*z2Avg;
            const su2double wAvg   = oneOvz0Avg*z3Avg;

            const su2double vnAvg  = uAvg*nx + vAvg*ny + wAvg*nz;
            const su2double unAvg  = vnAvg - gridVelNorm;
            const su2double kinAvg = 0.5*(uAvg*uAvg + vAvg*vAvg + wAvg*wAvg);
            const su2double a2Avg  = Gamma*p2Avg/rhoAvg;
            const su2double aAvg   = sqrt(a2Avg);
            const su2double HAvg   = a2Avg*ovgm1 + kinAvg;
            const su2double EAvg   = HAvg - p2Avg/rhoAvg;

            const su2double ovaAvg  = 1.0/aAvg;
            const su2double ova2Avg = 1.0/a2Avg;

            /*--- Compute the difference in entropy variables. ---*/
            const su2double dV0 = V0R - V0L, dV1 = V1R - V1L, dV2 = V2R - V2L,
                            dV3 = V3R - V3L, dV4 = V4R - V4L;

            /*--- Define the difference in conservative variables as dU/dV deltaV, where the
                  transformation matrix dU/dV must be evaluated at the averaged state. ---*/
            su2double abv1 = rhoAvg*(uAvg*dV1 + vAvg*dV2 + wAvg*dV3);
            su2double abv2 = abv1 + rhoAvg*(dV0 + HAvg*dV4);

            const su2double dr  = abv1 + rhoAvg*(dV0 + EAvg*dV4);
            const su2double dru = uAvg*abv2 + p2Avg*dV1;
            const su2double drv = vAvg*abv2 + p2Avg*dV2;
            const su2double drw = wAvg*abv2 + p2Avg*dV3;
            const su2double drE = HAvg*abv1 + rhoAvg*EAvg*(dV0 + HAvg*dV4)
                                + p2Avg*kinAvg*dV4;

            /*--- Compute the absolute values of the eigenvalues of the flux Jacobian. ---*/
            su2double lam1 = fabs(unAvg + aAvg);
            su2double lam2 = fabs(unAvg - aAvg);
            su2double lam3 = fabs(unAvg);

            /*--- Scale the acoustic eigenvalue, such that the EC2 (or EC1) flux of Ismail
                  and Roe is obtained. Also multiply all eigenvalues by half to obtain
                  the correct scaling. ---*/
            lam1 *= 0.5*(1.0 + beta + alphaMax);
            lam2 *= 0.5*(1.0 + beta + alphaMax);
            lam3 *= 0.5;

            /*--- Some abbreviations, which occur quite often in the dissipation terms. ---*/
            abv1 = 0.5*(lam1 + lam2);
            abv2 = 0.5*(lam1 - lam2);

            const su2double abv3 = abv1 - lam3;
            const su2double abv4 = gm1*(kinAvg*dr - uAvg*dru - vAvg*drv - wAvg*drw + drE);
            const su2double abv5 = nx*dru + ny*drv + nz*drw - vnAvg*dr;
            const su2double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
            const su2double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

            /*--- Compute the fluxes. ---*/
            fluxes(i,0) = Area*(rhoAvg*unAvg - lam3*dr -abv6);
            fluxes(i,1) = Area*(rhoAvg*unAvg*uAvg + p1Avg*nx - lam3*dru - uAvg*abv6 - nx*abv7);
            fluxes(i,2) = Area*(rhoAvg*unAvg*vAvg + p1Avg*ny - lam3*drv - vAvg*abv6 - ny*abv7);
            fluxes(i,3) = Area*(rhoAvg*unAvg*wAvg + p1Avg*nz - lam3*drw - wAvg*abv6 - nz*abv7);
            fluxes(i,4) = Area*(rhoAvg*unAvg*HAvg + p2Avg*gridVelNorm  - lam3*drE - HAvg*abv6 - vnAvg*abv7);
          }

          break;
        }
      }

      break;
    }

    case UPWIND::ROE: {

      /*--- Roe's approximate Riemann solver. Easier storage of the cut off
            value for the entropy correction. ---*/
      const su2double Delta = config->GetEntropyFix_Coeff();

      /*--- Make a distinction between two and three space dimensions
            in order to have the most efficient code. ---*/
      switch( nDim ) {

        case 2: {

          /*--- Loop over the items. ---*/
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned int i=0; i<nItems; ++i) {

            /*--- Easier storage of the components of the unit normal and
                  half the Jacobian (which is half the area of the face). ---*/
            const su2double nx = normalsFace(i,0);
            const su2double ny = normalsFace(i,1);
            const su2double halfArea = 0.5*JacobiansFace(i);

            /*--- Compute the normal grid velocity. ---*/
            const su2double gridVelNorm = gridVelocities(i,0)*nx + gridVelocities(i,1)*ny;

            /*--- Easier storage of the primitive variables of
                  the left and right state. ---*/
            const su2double rhoL = solLeft(i,0), rhoR = solRight(i,0);
            const su2double uL   = solLeft(i,1), uR   = solRight(i,1);
            const su2double vL   = solLeft(i,2), vR   = solRight(i,2);
            const su2double pL   = solLeft(i,3), pR   = solRight(i,3);

            /*--- Compute the total energy of the left and right state. ---*/
            const su2double rEL = ovgm1*pL + 0.5*rhoL*(uL*uL + vL*vL);
            const su2double rER = ovgm1*pR + 0.5*rhoR*(uR*uR + vR*vR);

            /*--- Compute the difference of the conservative mean flow variables. ---*/
            const su2double dr  = rhoR    - rhoL;
            const su2double dru = rhoR*uR - rhoL*uL;
            const su2double drv = rhoR*vR - rhoL*vL;
            const su2double drE = rER     - rEL;

            /*--- Compute the Roe average state. ---*/
            const su2double zL = sqrt(rhoL);
            const su2double zR = sqrt(rhoR);
            su2double tmp      = 1.0/(zL + zR);

            const su2double rHL = rEL + pL;
            const su2double rHR = rER + pR;

            const su2double uAvg = tmp*(zL*uL + zR*uR);
            const su2double vAvg = tmp*(zL*vL + zR*vR);
            const su2double HAvg = tmp*(rHL/zL + rHR/zR);

            /*--- Compute from the Roe average state some variables, which occur
                  quite often in the matrix vector product to be computed. ---*/
            const su2double alphaAvg = 0.5*(uAvg*uAvg + vAvg*vAvg);
            tmp                      = gm1*(HAvg - alphaAvg);
            const su2double a2Avg    = fabs(tmp);
            const su2double aAvg     = sqrt(a2Avg);
            const su2double vnAvg    = uAvg*nx + vAvg*ny;
            const su2double unAvg    = vnAvg - gridVelNorm;
            const su2double ovaAvg   = 1.0/aAvg;
            const su2double ova2Avg  = 1.0/a2Avg;

            /*--- Compute the absolute values of the three eigenvalues and
                  apply the entropy correction. ---*/
            su2double lam1 = fabs(unAvg + aAvg);
            su2double lam2 = fabs(unAvg - aAvg);
            su2double lam3 = fabs(unAvg);

            tmp  = Delta*max(lam1, lam2);
            lam1 = max(lam1, tmp);
            lam2 = max(lam2, tmp);
            lam3 = max(lam3, tmp);

            /*--- Some abbreviations, which occur quite often in the dissipation terms. ---*/
            const su2double abv1 = 0.5*(lam1 + lam2);
            const su2double abv2 = 0.5*(lam1 - lam2);
            const su2double abv3 = abv1 - lam3;

            const su2double abv4 = gm1*(alphaAvg*dr - uAvg*dru - vAvg*drv + drE);
            const su2double abv5 = nx*dru + ny*drv - vnAvg*dr;
            const su2double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
            const su2double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

            /*--- Compute the Roe flux vector, which is 0.5*(FL + FR - |A|(UR-UL)). ---*/
            const su2double vnL = uL*nx + vL*ny;
            const su2double vnR = uR*nx + vR*ny;
            const su2double unL = vnL - gridVelNorm;
            const su2double unR = vnR - gridVelNorm;
            const su2double pa  = pL + pR;

            fluxes(i,0) = halfArea*(rhoL*unL + rhoR*unR - (lam3*dr + abv6));
            fluxes(i,1) = halfArea*(rhoL*uL*unL + rhoR*uR*unR + pa*nx
                        -           (lam3*dru + uAvg*abv6 + nx*abv7));
            fluxes(i,2) = halfArea*(rhoL*vL*unL + rhoR*vR*unR + pa*ny
                        -           (lam3*drv + vAvg*abv6 + ny*abv7));
            fluxes(i,3) = halfArea*(rEL*unL + rER*unR + pL*vnL + pR*vnR
                        -           (lam3*drE + HAvg*abv6 + vnAvg*abv7));
          }

          break;
        }

        /*--------------------------------------------------------------------*/

        case 3: {

          /*--- Loop over the items. ---*/
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned int i=0; i<nItems; ++i) {

            /*--- Easier storage of the components of the unit normal and
                  half the Jacobian (which is half the area of the face). ---*/
            const su2double nx = normalsFace(i,0);
            const su2double ny = normalsFace(i,1);
            const su2double nz = normalsFace(i,2);
            const su2double halfArea = 0.5*JacobiansFace(i);

            /*--- Compute the normal grid velocity. ---*/
            const su2double gridVelNorm = gridVelocities(i,0)*nx + gridVelocities(i,1)*ny
                                        + gridVelocities(i,2)*nz;

            /*--- Easier storage of the primitive variables of
                  the left and right state. ---*/
            const su2double rhoL = solLeft(i,0), rhoR = solRight(i,0);
            const su2double uL   = solLeft(i,1), uR   = solRight(i,1);
            const su2double vL   = solLeft(i,2), vR   = solRight(i,2);
            const su2double wL   = solLeft(i,3), wR   = solRight(i,3);
            const su2double pL   = solLeft(i,4), pR   = solRight(i,4);

            /*--- Compute the total energy of the left and right state. ---*/
            const su2double rEL = ovgm1*pL + 0.5*rhoL*(uL*uL + vL*vL + wL*wL);
            const su2double rER = ovgm1*pR + 0.5*rhoR*(uR*uR + vR*vR + wR*wR);

            /*--- Compute the difference of the conservative mean flow variables. ---*/
            const su2double dr  = rhoR    - rhoL;
            const su2double dru = rhoR*uR - rhoL*uL;
            const su2double drv = rhoR*vR - rhoL*vL;
            const su2double drw = rhoR*wR - rhoL*wL;
            const su2double drE = rER     - rEL;

            /*--- Compute the Roe average state. ---*/
            const su2double zL = sqrt(rhoL);
            const su2double zR = sqrt(rhoR);
            su2double tmp      = 1.0/(zL + zR);

            const su2double rHL = rEL + pL;
            const su2double rHR = rER + pR;

            const su2double uAvg = tmp*(zL*uL + zR*uR);
            const su2double vAvg = tmp*(zL*vL + zR*vR);
            const su2double wAvg = tmp*(zL*wL + zR*wR);
            const su2double HAvg = tmp*(rHL/zL + rHR/zR);

            /*--- Compute from the Roe average state some variables, which occur
                  quite often in the matrix vector product to be computed. ---*/
            const su2double alphaAvg = 0.5*(uAvg*uAvg + vAvg*vAvg + wAvg*wAvg);
            tmp                      = gm1*(HAvg - alphaAvg);
            const su2double a2Avg    = fabs(tmp);
            const su2double aAvg     = sqrt(a2Avg);
            const su2double vnAvg    = uAvg*nx + vAvg*ny + wAvg*nz;
            const su2double unAvg    = vnAvg - gridVelNorm;
            const su2double ovaAvg   = 1.0/aAvg;
            const su2double ova2Avg  = 1.0/a2Avg;

            /*--- Compute the absolute values of the three eigenvalues and
                  apply the entropy correction. ---*/
            su2double lam1 = fabs(unAvg + aAvg);
            su2double lam2 = fabs(unAvg - aAvg);
            su2double lam3 = fabs(unAvg);

            tmp  = Delta*max(lam1, lam2);
            lam1 = max(lam1, tmp);
            lam2 = max(lam2, tmp);
            lam3 = max(lam3, tmp);

            /*--- Some abbreviations, which occur quite often in the dissipation terms. ---*/
            const su2double abv1 = 0.5*(lam1 + lam2);
            const su2double abv2 = 0.5*(lam1 - lam2);
            const su2double abv3 = abv1 - lam3;

            const su2double abv4 = gm1*(alphaAvg*dr - uAvg*dru - vAvg*drv -wAvg*drw + drE);
            const su2double abv5 = nx*dru + ny*drv + nz*drw - vnAvg*dr;
            const su2double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
            const su2double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

            /*--- Compute the Roe flux vector, which is 0.5*(FL + FR - |A|(UR-UL)). ---*/
            const su2double vnL = uL*nx + vL*ny + wL*nz;
            const su2double vnR = uR*nx + vR*ny + wR*nz;
            const su2double unL = vnL - gridVelNorm;
            const su2double unR = vnR - gridVelNorm;
            const su2double pa  = pL + pR;

            fluxes(i,0) = halfArea*(rhoL*unL + rhoR*unR - (lam3*dr + abv6));
            fluxes(i,1) = halfArea*(rhoL*uL*unL + rhoR*uR*unR + pa*nx
                        -           (lam3*dru + uAvg*abv6 + nx*abv7));
            fluxes(i,2) = halfArea*(rhoL*vL*unL + rhoR*vR*unR + pa*ny
                        -           (lam3*drv + vAvg*abv6 + ny*abv7));
            fluxes(i,3) = halfArea*(rhoL*wL*unL + rhoR*wR*unR + pa*nz
                        -           (lam3*drw + wAvg*abv6 + nz*abv7));
            fluxes(i,4) = halfArea*(rEL*unL + rER*unR + pL*vnL + pR*vnR
                        -           (lam3*drE + HAvg*abv6 + vnAvg*abv7));
          }

          break;
        }
      }

      break;
    }

    /*------------------------------------------------------------------------*/

    case UPWIND::LAX_FRIEDRICH: {

      /*--- Local Lax-Friedrich (Rusanov) flux Make a distinction between two and
            three space dimensions in order to have the most efficient code. ---*/
      switch( nDim ) {

        case 2: {

          /*--- Loop over the items. ---*/
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned int i=0; i<nItems; ++i) {

            /*--- Easier storage of the components of the unit normal and
                  half the Jacobian (which is half the area of the face). ---*/
            const su2double nx = normalsFace(i,0);
            const su2double ny = normalsFace(i,1);
            const su2double halfArea = 0.5*JacobiansFace(i);

            /*--- Compute the normal grid velocity. ---*/
            const su2double gridVelNorm = gridVelocities(i,0)*nx + gridVelocities(i,1)*ny;

            /*--- Easier storage of the primitive variables of
                  the left and right state. ---*/
            const su2double rhoL = solLeft(i,0),  rhoR = solRight(i,0);
            const su2double uL   = solLeft(i,1),  uR   = solRight(i,1);
            const su2double vL   = solLeft(i,2),  vR   = solRight(i,2);
            const su2double pL   = solLeft(i,3),  pR   = solRight(i,3);

            /*--- Compute the speed of sound squared and the total energy of
                  the left and right state. ---*/
            const su2double a2L = Gamma*pL/rhoL;
            const su2double a2R = Gamma*pR/rhoR;
            const su2double rEL = ovgm1*pL + 0.5*rhoL*(uL*uL + vL*vL);
            const su2double rER = ovgm1*pR + 0.5*rhoR*(uR*uR + vR*vR);

            /*--- Compute the difference of the conservative mean flow variables. ---*/
            const su2double dr  = rhoR    - rhoL;
            const su2double dru = rhoR*uR - rhoL*uL;
            const su2double drv = rhoR*vR - rhoL*vL;
            const su2double drE = rER     - rEL;

            /*--- Compute the spectral radii of the left and right state
                  and take the maximum for the dissipation terms. ---*/
            const su2double vnL = uL*nx + vL*ny;
            const su2double vnR = uR*nx + vR*ny;
            const su2double unL = vnL - gridVelNorm;
            const su2double unR = vnR - gridVelNorm;

            const su2double radL = fabs(unL) + sqrt(fabs(a2L));
            const su2double radR = fabs(unR) + sqrt(fabs(a2R));
            const su2double rad  = max(radL, radR);

            /*--- Compute the flux vector, which is 0.5*(FL + FR - rad(UR-UL)). ---*/
            const su2double pa = pL + pR;

            fluxes(i,0) = halfArea*(rhoL*unL + rhoR*unR - rad*dr);
            fluxes(i,1) = halfArea*(rhoL*uL*unL + rhoR*uR*unR + pa*nx - rad*dru);
            fluxes(i,2) = halfArea*(rhoL*vL*unL + rhoR*vR*unR + pa*ny - rad*drv);
            fluxes(i,3) = halfArea*(rEL*unL + rER*unR + pL*vnL + pR*vnR - rad*drE);
          }

          break;
        }

        /*--------------------------------------------------------------------*/

        case 3: {

          /*--- Loop over the items. ---*/
          SU2_OMP_SIMD_IF_NOT_AD
          for(unsigned int i=0; i<nItems; ++i) {

            /*--- Easier storage of the components of the unit normal and
                  half the Jacobian (which is half the area of the face). ---*/
            const su2double nx = normalsFace(i,0);
            const su2double ny = normalsFace(i,1);
            const su2double nz = normalsFace(i,2);
            const su2double halfArea = 0.5*JacobiansFace(i);

            /*--- Compute the normal grid velocity. ---*/
            const su2double gridVelNorm = gridVelocities(i,0)*nx + gridVelocities(i,1)*ny
                                        + gridVelocities(i,2)*nz;

            /*--- Easier storage of the primitive variables of
                  the left and right state. ---*/
            const su2double rhoL = solLeft(i,0), rhoR = solRight(i,0);
            const su2double uL   = solLeft(i,1), uR   = solRight(i,1);
            const su2double vL   = solLeft(i,2), vR   = solRight(i,2);
            const su2double wL   = solLeft(i,3), wR   = solRight(i,3);
            const su2double pL   = solLeft(i,4), pR   = solRight(i,4);

            /*--- Compute the speed of sound squared and the total energy of
                  the left and right state. ---*/
            const su2double a2L = Gamma*pL/rhoL;
            const su2double a2R = Gamma*pR/rhoR;
            const su2double rEL = ovgm1*pL + 0.5*rhoL*(uL*uL + vL*vL + wL*wL);
            const su2double rER = ovgm1*pR + 0.5*rhoR*(uR*uR + vR*vR + wR*wR);

            /*--- Compute the difference of the conservative mean flow variables. ---*/
            const su2double dr  = rhoR    - rhoL;
            const su2double dru = rhoR*uR - rhoL*uL;
            const su2double drv = rhoR*vR - rhoL*vL;
            const su2double drw = rhoR*wR - rhoL*wL;
            const su2double drE = rER     - rEL;

            /*--- Compute the spectral radii of the left and right state
                  and take the maximum for the dissipation terms. ---*/
            const su2double vnL = uL*nx + vL*ny + wL*nz;
            const su2double vnR = uR*nx + vR*ny + wR*nz;
            const su2double unL = vnL - gridVelNorm;
            const su2double unR = vnR - gridVelNorm;

            const su2double radL = fabs(unL) + sqrt(fabs(a2L));
            const su2double radR = fabs(unR) + sqrt(fabs(a2R));
            const su2double rad  = max(radL, radR);

            /*--- Compute the flux vector, which is 0.5*(FL + FR - rad(UR-UL)). ---*/
            const su2double pa  = pL + pR;

            fluxes(i,0) = halfArea*(rhoL*unL + rhoR*unR - rad*dr);
            fluxes(i,1) = halfArea*(rhoL*uL*unL + rhoR*uR*unR + pa*nx - rad*dru);
            fluxes(i,2) = halfArea*(rhoL*vL*unL + rhoR*vR*unR + pa*ny - rad*drv);
            fluxes(i,3) = halfArea*(rhoL*wL*unL + rhoR*wR*unR + pa*nz - rad*drw);
            fluxes(i,4) = halfArea*(rEL*unL + rER*unR + pL*vnL + pR*vnR - rad*drE);
          }

          break;
        }
      }

      break;
    }

    /*------------------------------------------------------------------------*/

    default: {

      if(rank == MASTER_NODE) {
        SU2_OMP_SINGLE
        SU2_MPI::Error(string("Riemann solver not implemented yet"), CURRENT_FUNCTION);
        END_SU2_OMP_SINGLE
      }

      SU2_OMP_SINGLE
      SU2_MPI::Barrier(SU2_MPI::GetComm());
      END_SU2_OMP_SINGLE
    }
  }
}

void CFEM_DG_EulerSolver::LoadRestart(CGeometry **geometry,
                                      CSolver   ***solver,
                                      CConfig   *config,
                                      int       val_iter,
                                      bool      val_update_geo) {

  /*--- The actual reading is done by a single thread. ---*/
  SU2_OMP_SINGLE
  {
    /*--- Read the restarRestart_Varst data from either an ASCII or binary SU2 file. ---*/
    string restart_filename = config->GetSolution_FileName();
    restart_filename = config->GetFilename(restart_filename, "", val_iter);

    if (config->GetRead_Binary_Restart()) {
      Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
    } else {
      Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
    }

    /*--- Determine a map from the ID of the global DOF to the local
          element and DOF of the element. ---*/
    map<unsigned long, CUnsignedLong2T>  mapGlobalDOFToLocal;
    for(unsigned long i=0; i<nVolElemOwned; ++i) {
      for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {
        unsigned long GlobalDOF = volElem[i].offsetDOFsSolGlobal+j;
        mapGlobalDOFToLocal[GlobalDOF] = CUnsignedLong2T(i,j);
      }
    }

    /*--- Skip coordinates ---*/
    unsigned short skipVars = nDim;

    /*--- Loop over the global DOFs and determine whether they are
          stored locally. ---*/
    unsigned long counter = 0;
    for (unsigned long iPoint_Global=0; iPoint_Global<geometry[MESH_0]->GetGlobal_nPointDomain(); ++iPoint_Global) {

      /*--- Check if the DOF is stored on this rank. ---*/
      auto it = mapGlobalDOFToLocal.find(iPoint_Global);
      if(it != mapGlobalDOFToLocal.cend()) {

        /*--- Retrieve the element and DOF inside the element. ---*/
        const unsigned long i = it->second.long0;
        const unsigned long j = it->second.long1;

        /*--- Determine the start index in the read buffer and
              update the counter afterwards. ---*/
        const unsigned long index = counter*Restart_Vars[1] + skipVars;
        ++counter;

        /*--- Store the conservative variables in the local data structures. ---*/
        for(unsigned short iVar=0; iVar<nVar; ++iVar)
          volElem[i].solDOFs(j,iVar) = Restart_Data[index+iVar];
      }
    }
  
    /*--- Detect a wrong solution file ---*/
    unsigned short rbuf_NotMatching = 0;
    if(counter < nDOFsLocOwned) rbuf_NotMatching = 1;

#ifdef HAVE_MPI
    unsigned short sbuf_NotMatching = rbuf_NotMatching;
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_MAX, SU2_MPI::GetComm());
#endif

    if (rbuf_NotMatching != 0)
      SU2_MPI::Error(string("The solution file ") + restart_filename.data() +
                     string(" doesn't match with the mesh file!\n") +
                     string("It could be empty lines at the end of the file."),
                     CURRENT_FUNCTION);

    /*--- Delete the class memory that is used to load the restart. ---*/
    delete[] Restart_Vars;
    delete[] Restart_Data;
    Restart_Vars = nullptr; Restart_Data = nullptr;

    /* Initialize the counter for the number of bad DOFs. */
    counter = 0;
  }
  END_SU2_OMP_SINGLE

  /*--- Definition of the number of bad elements for this thread.
        The reduction is handled manually to avoid complications
        with CODIPACK. ---*/
  unsigned long nDOFsBad = 0;

  /*--- Loop over the owned elements. ---*/
#ifdef HAVE_OMP
  const size_t omp_chunk_size_vol = computeStaticChunkSize(nVolElemOwned, omp_get_num_threads(), 64);
#endif
  SU2_OMP_FOR_STAT(omp_chunk_size_vol)
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /*--- Loop over the DOFs of this element, which currently
          stores the conservative variables. ---*/
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j) {

      /*--- Compute the primitive variables. ---*/
      const su2double rho    = volElem[i].solDOFs(j,0);
      const su2double rhoInv = 1.0/rho;

      su2double Velocity2 = 0.0;
      for(unsigned short iDim=1; iDim<=nDim; ++iDim) {
        const su2double vel = volElem[i].solDOFs(j,iDim)*rhoInv;
        Velocity2 += vel*vel;
      }

      const su2double StaticEnergy = volElem[i].solDOFs(j,nDim+1)*rhoInv - 0.5*Velocity2;

      GetFluidModel()->SetTDState_rhoe(rho, StaticEnergy);
      const su2double Pressure = GetFluidModel()->GetPressure();
      const su2double Temperature = GetFluidModel()->GetTemperature();

      /*--- Check for negative pressure, density or temperature. ---*/
      if((Pressure < 0.0) || (rho < 0.0) || (Temperature < 0.0)) {

        /*--- Reset the state to the infinity state and update nDOFsBad. ---*/
        for(unsigned short k=0; k<nVar; ++k)
          volElem[i].solDOFs(j,k) = ConsVarFreeStream[k];
        ++nDOFsBad;
      }
    }

    /*--- Convert the conservative variables to entropy variables. ---*/
    ConservativeToEntropyVariables(volElem[i].solDOFs);

    /*--- Convert the nodal solution to the modal solution. ---*/
    volElem[i].NodalToModalFlow();
  }
  END_SU2_OMP_FOR

  /*--- Carry out the reduction over the threads. ---*/
  if (config->GetComm_Level() == COMM_FULL) {
    SU2_OMP_ATOMIC
    counter += nDOFsBad;

    SU2_OMP_BARRIER
  }

  /*--- Warning message about non-physical points ---*/
  if (config->GetComm_Level() == COMM_FULL) {

    SU2_OMP_SINGLE
    {
#ifdef HAVE_MPI
      unsigned long nBadDOFsLoc = counter;
      SU2_MPI::Reduce(&nBadDOFsLoc, &nDOFsBad, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
#endif

      if((rank == MASTER_NODE) && (nDOFsBad != 0))
        cout << "Warning. The initial solution contains "<< nDOFsBad << " DOFs that are not physical." << endl;
    }
    END_SU2_OMP_SINGLE
  }
}
