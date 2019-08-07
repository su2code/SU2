/*!
 * \file CEulerVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, T. Economon
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

#include "../../include/variables/CEulerVariable.hpp"

CEulerVariable::CEulerVariable(void) : CVariable() {

  /*--- Array initialization ---*/

  HB_Source = NULL;
  Primitive = NULL;
  Secondary = NULL;

  Gradient_Primitive = NULL;
  Gradient_Secondary = NULL;

  Limiter_Primitive = NULL;
  Limiter_Secondary = NULL;

  WindGust    = NULL;
  WindGustDer = NULL;

  nPrimVar     = 0;
  nPrimVarGrad = 0;

  nSecondaryVar     = 0;
  nSecondaryVarGrad = 0;

  Solution_New = NULL;

  Solution_BGS_k = NULL;
}

CEulerVariable::CEulerVariable(su2double val_density, su2double *val_velocity, su2double val_energy, unsigned short val_nDim,
                               unsigned short val_nvar, CConfig *config) : CVariable(val_nDim, val_nvar, config) {
    unsigned short iVar, iDim, iMesh, nMGSmooth = 0;

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool viscous = config->GetViscous();
  bool windgust = config->GetWind_Gust();
  bool classical_rk4 = (config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT);
  bool fsi = config->GetFSI_Simulation();
  bool multizone = config->GetMultizone_Problem();

  /*--- Array initialization ---*/

  HB_Source = NULL;
  Primitive = NULL;
  Secondary = NULL;

  Gradient_Primitive = NULL;
  Gradient_Secondary = NULL;

  Limiter_Primitive = NULL;
  Limiter_Secondary = NULL;

  WindGust    = NULL;
  WindGustDer = NULL;

  nPrimVar     = 0;
  nPrimVarGrad = 0;

  nSecondaryVar     = 0;
  nSecondaryVarGrad = 0;

  Solution_New = NULL;

  /*--- Allocate and initialize the primitive variables and gradients ---*/
  nPrimVar = nDim+9; nPrimVarGrad = nDim+4;
  if (viscous) { nSecondaryVar = 8; nSecondaryVarGrad = 2; }
  else { nSecondaryVar = 2; nSecondaryVarGrad = 2; }


  /*--- Allocate residual structures ---*/

  Res_TruncError = new su2double [nVar];

  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }

  /*--- Only for residual smoothing (multigrid) ---*/

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);

  if (nMGSmooth > 0) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }

  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    Undivided_Laplacian = new su2double [nVar];
  }

  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/

  Limiter_Primitive = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Limiter_Primitive[iVar] = 0.0;

  Limiter_Secondary = new su2double [nSecondaryVarGrad];
  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
    Limiter_Secondary[iVar] = 0.0;

  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;

  Solution_Max = new su2double [nPrimVarGrad];
  Solution_Min = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }

  /*--- Solution and old solution initialization ---*/

  Solution[0] = val_density;
  Solution_Old[0] = val_density;
  for (iDim = 0; iDim < nDim; iDim++) {
    Solution[iDim+1] = val_density*val_velocity[iDim];
    Solution_Old[iDim+1] = val_density*val_velocity[iDim];
  }
  Solution[nVar-1] = val_density*val_energy;
  Solution_Old[nVar-1] = val_density*val_energy;

  /*--- New solution initialization for Classical RK4 ---*/

  if (classical_rk4) {
    Solution_New = new su2double[nVar];
    Solution_New[0] = val_density;
    for (iDim = 0; iDim < nDim; iDim++) {
      Solution_New[iDim+1] = val_density*val_velocity[iDim];
    }
    Solution_New[nVar-1] = val_density*val_energy;
  }

    /*--- Allocate and initialize solution for dual time strategy ---*/

  if (dual_time) {
    Solution_time_n[0] = val_density;
    Solution_time_n1[0] = val_density;
    for (iDim = 0; iDim < nDim; iDim++) {
      Solution_time_n[iDim+1] = val_density*val_velocity[iDim];
      Solution_time_n1[iDim+1] = val_density*val_velocity[iDim];
    }
    Solution_time_n[nVar-1] = val_density*val_energy;
    Solution_time_n1[nVar-1] = val_density*val_energy;
  }


  /*--- Allocate space for the harmonic balance source terms ---*/

  if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    HB_Source = new su2double[nVar];
    for (iVar = 0; iVar < nVar; iVar++) HB_Source[iVar] = 0.0;
  }

  /*--- Allocate vector for wind gust and wind gust derivative field ---*/

  if (windgust) {
    WindGust = new su2double [nDim];
    WindGustDer = new su2double [nDim+1];
  }

  /*--- Incompressible flow, primitive variables nDim+3, (P, vx, vy, vz, rho, beta) ---*/

  Primitive = new su2double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;

  Secondary = new su2double [nSecondaryVar];
  for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar] = 0.0;

  /*--- Compressible flow, gradients primitive variables nDim+4, (T, vx, vy, vz, P, rho, h)
        We need P, and rho for running the adjoint problem ---*/

  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }

  Gradient_Secondary = new su2double* [nSecondaryVarGrad];
  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
    Gradient_Secondary[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Secondary[iVar][iDim] = 0.0;
  }

  Solution_BGS_k = NULL;
  if (fsi || multizone){
      Solution_BGS_k  = new su2double [nVar];
      Solution[0] = val_density;
      for (iDim = 0; iDim < nDim; iDim++) {
        Solution_BGS_k[iDim+1] = val_density*val_velocity[iDim];
      }
      Solution_BGS_k[nVar-1] = val_density*val_energy;
  }

}

CEulerVariable::CEulerVariable(su2double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config) : CVariable(val_nDim, val_nvar, config) {
    unsigned short iVar, iDim, iMesh, nMGSmooth = 0;

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool viscous = config->GetViscous();
  bool windgust = config->GetWind_Gust();
  bool classical_rk4 = (config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT);
  bool fsi = config->GetFSI_Simulation();
  bool multizone = config->GetMultizone_Problem();

  /*--- Array initialization ---*/

  HB_Source = NULL;
  Primitive = NULL;
  Secondary = NULL;

  Gradient_Primitive = NULL;
  Gradient_Secondary = NULL;

  Limiter_Primitive = NULL;
  Limiter_Secondary = NULL;

  WindGust    = NULL;
  WindGustDer = NULL;

  nPrimVar     = 0;
  nPrimVarGrad = 0;

  nSecondaryVar     = 0;
  nSecondaryVarGrad = 0;

  Solution_New = NULL;

  /*--- Allocate and initialize the primitive variables and gradients ---*/

  nPrimVar = nDim+9; nPrimVarGrad = nDim+4;
  if (viscous) { nSecondaryVar = 8; nSecondaryVarGrad = 2; }
  else { nSecondaryVar = 2; nSecondaryVarGrad = 2; }


  /*--- Allocate residual structures ---*/

  Res_TruncError = new su2double [nVar];

  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }

  /*--- Only for residual smoothing (multigrid) ---*/

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);

  if (nMGSmooth > 0) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }

  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
    Undivided_Laplacian = new su2double [nVar];

  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/

  Limiter_Primitive = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Limiter_Primitive[iVar] = 0.0;

  Limiter_Secondary = new su2double [nSecondaryVarGrad];
  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
    Limiter_Secondary[iVar] = 0.0;

  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;

  Solution_Max = new su2double [nPrimVarGrad];
  Solution_Min = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }

  /*--- Solution initialization ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Solution[iVar] = val_solution[iVar];
    Solution_Old[iVar] = val_solution[iVar];
  }

  /*--- New solution initialization for Classical RK4 ---*/

  if (classical_rk4) {
    Solution_New = new su2double[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Solution_New[iVar] = val_solution[iVar];
    }
  }

  /*--- Allocate and initializate solution for dual time strategy ---*/

  if (dual_time) {
    Solution_time_n = new su2double [nVar];
    Solution_time_n1 = new su2double [nVar];

    for (iVar = 0; iVar < nVar; iVar++) {
      Solution_time_n[iVar] = val_solution[iVar];
      Solution_time_n1[iVar] = val_solution[iVar];
    }
  }

  /*--- Allocate space for the harmonic balance source terms ---*/

  if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    HB_Source = new su2double[nVar];
    for (iVar = 0; iVar < nVar; iVar++) HB_Source[iVar] = 0.0;
  }

  /*--- Allocate vector for wind gust and wind gust derivative field ---*/

  if (windgust) {
    WindGust = new su2double [nDim];
    WindGustDer = new su2double [nDim+1];
  }

  /*--- Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c) ---*/

  Primitive = new su2double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;

  Secondary = new su2double [nSecondaryVar];
  for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar] = 0.0;


  /*--- Compressible flow, gradients primitive variables nDim+4, (T, vx, vy, vz, P, rho, h)
        We need P, and rho for running the adjoint problem ---*/

  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }

  Gradient_Secondary = new su2double* [nSecondaryVarGrad];
  for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) {
    Gradient_Secondary[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Secondary[iVar][iDim] = 0.0;
  }

  Solution_BGS_k = NULL;
  if (fsi || multizone){
      Solution_BGS_k  = new su2double [nVar];
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_BGS_k[iVar] = val_solution[iVar];
      }
  }

}

CEulerVariable::~CEulerVariable(void) {
    unsigned short iVar;

  if (HB_Source         != NULL) delete [] HB_Source;
  if (Primitive         != NULL) delete [] Primitive;
  if (Secondary         != NULL) delete [] Secondary;
  if (Limiter_Primitive != NULL) delete [] Limiter_Primitive;
  if (Limiter_Secondary != NULL) delete [] Limiter_Secondary;
  if (WindGust          != NULL) delete [] WindGust;
  if (WindGustDer       != NULL) delete [] WindGustDer;

  if (Gradient_Primitive != NULL) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      if (Gradient_Primitive[iVar] != NULL) delete [] Gradient_Primitive[iVar];
    delete [] Gradient_Primitive;
  }
  if (Gradient_Secondary != NULL) {
    for (iVar = 0; iVar < nSecondaryVarGrad; iVar++)
      if (Gradient_Secondary[iVar] != NULL) delete [] Gradient_Secondary[iVar];
    delete [] Gradient_Secondary;
  }

  if (Solution_New != NULL) delete [] Solution_New;

  if (Solution_BGS_k  != NULL) delete [] Solution_BGS_k;

}

void CEulerVariable::SetGradient_PrimitiveZero(unsigned short val_primvar) {
    unsigned short iVar, iDim;

    for (iVar = 0; iVar < val_primvar; iVar++)
        for (iDim = 0; iDim < nDim; iDim++)
            Gradient_Primitive[iVar][iDim] = 0.0;
}

void CEulerVariable::SetGradient_SecondaryZero(unsigned short val_secondaryvar) {
    unsigned short iVar, iDim;

    for (iVar = 0; iVar < val_secondaryvar; iVar++)
        for (iDim = 0; iDim < nDim; iDim++)
            Gradient_Secondary[iVar][iDim] = 0.0;
}

su2double CEulerVariable::GetProjVel(su2double *val_vector) {
    su2double ProjVel;
    unsigned short iDim;

    ProjVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
        ProjVel += Primitive[iDim+1]*val_vector[iDim];

    return ProjVel;
}

bool CEulerVariable::SetPrimVar(CFluidModel *FluidModel) {
  unsigned short iVar;
  bool check_dens = false, check_press = false, check_sos = false, check_temp = false, RightVol = true;


  SetVelocity();   // Computes velocity and velocity^2
  su2double density = GetDensity();
  su2double staticEnergy = GetEnergy()-0.5*Velocity2;

  /*--- Check will be moved inside fluid model plus error description strings ---*/

  FluidModel->SetTDState_rhoe(density, staticEnergy);

  check_dens = SetDensity();
  check_press = SetPressure(FluidModel->GetPressure());
  check_sos = SetSoundSpeed(FluidModel->GetSoundSpeed2());
  check_temp = SetTemperature(FluidModel->GetTemperature());

  /*--- Check that the solution has a physical meaning ---*/

  if (check_dens || check_press || check_sos || check_temp) {

    /*--- Copy the old solution ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];

    /*--- Recompute the primitive variables ---*/

    SetVelocity();   // Computes velocity and velocity^2
    su2double density = GetDensity();
    su2double staticEnergy = GetEnergy()-0.5*Velocity2;
    /* check will be moved inside fluid model plus error description strings*/
    FluidModel->SetTDState_rhoe(density, staticEnergy);

    SetDensity();
    SetPressure(FluidModel->GetPressure());
    SetSoundSpeed(FluidModel->GetSoundSpeed2());
    SetTemperature(FluidModel->GetTemperature());

    RightVol = false;

  }

  /*--- Set enthalpy ---*/

  SetEnthalpy();                                // Requires pressure computation.

  return RightVol;

}

void CEulerVariable::SetSecondaryVar(CFluidModel *FluidModel) {

   /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/

   SetdPdrho_e(FluidModel->GetdPdrho_e());
   SetdPde_rho(FluidModel->GetdPde_rho());

}

