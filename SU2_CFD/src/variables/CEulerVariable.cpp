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


CEulerVariable::CEulerVariable(su2double density, const su2double *velocity, su2double energy, unsigned long npoint,
                               unsigned long ndim, unsigned long nvar, CConfig *config) : CVariable(npoint, ndim, nvar, config) {

  bool dual_time = (config->GetTime_Marching() == DT_STEPPING_1ST) ||
                   (config->GetTime_Marching() == DT_STEPPING_2ND);
  bool viscous   = config->GetViscous();
  bool windgust  = config->GetWind_Gust();
  bool classical_rk4 = (config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT);

  /*--- Allocate and initialize the primitive variables and gradients ---*/

  nPrimVar          = nDim+9;
  nPrimVarGrad      = nDim+4;
  nSecondaryVar     = viscous? 8 : 2;
  nSecondaryVarGrad = 2;

  /*--- Allocate residual structures ---*/

  Res_TruncError.resize(nPoint,nVar) = su2double(0.0);

  /*--- Only for residual smoothing (multigrid) ---*/

  for (unsigned long iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    if (config->GetMG_CorrecSmooth(iMesh) > 0) {
      Residual_Sum.resize(nPoint,nVar);
      Residual_Old.resize(nPoint,nVar);
      break;
    }
  }

  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
    Undivided_Laplacian.resize(nPoint,nVar);

  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/

  Limiter_Primitive.resize(nPoint,nPrimVarGrad) = su2double(0.0);
  Limiter.resize(nPoint,nVar) = su2double(0.0);

  Solution_Max.resize(nPoint,nPrimVarGrad) = su2double(0.0);
  Solution_Min.resize(nPoint,nPrimVarGrad) = su2double(0.0);

  /*--- Solution initialization ---*/

  su2double val_solution[5] = {su2double(1.0), velocity[0], velocity[1], energy, energy};
  if(nDim==3) val_solution[3] = velocity[2];

  for(unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; ++iVar)
      Solution(iPoint,iVar) = density*val_solution[iVar];

  Solution_Old = Solution;

  /*--- New solution initialization for Classical RK4 ---*/

  if (classical_rk4) Solution_New = Solution;

  /*--- Allocate and initializate solution for dual time strategy ---*/

  if (dual_time) {
    Solution_time_n = Solution;
    Solution_time_n1 = Solution;
  }

  /*--- Allocate space for the harmonic balance source terms ---*/

  if (config->GetTime_Marching() == HARMONIC_BALANCE)
    HB_Source.resize(nPoint,nVar) = su2double(0.0);

  /*--- Allocate vector for wind gust and wind gust derivative field ---*/

  if (windgust) {
    WindGust.resize(nPoint,nDim);
    WindGustDer.resize(nPoint,nDim+1);
  }

  /*--- Compressible flow, primitive variables nDim+5, (T, vx, vy, vz, P, rho, h, c) ---*/

  Primitive.resize(nPoint,nPrimVar) = su2double(0.0);
  Secondary.resize(nPoint,nSecondaryVar) = su2double(0.0);

  /*--- Compressible flow, gradients primitive variables nDim+4, (T, vx, vy, vz, P, rho, h)
        We need P, and rho for running the adjoint problem ---*/

  Gradient_Primitive.resize(nPoint,nPrimVarGrad,nDim,0.0);

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    Rmatrix.resize(nPoint,nDim,nDim,0.0);
  }

  if (config->GetMultizone_Problem())
    Set_BGSSolution_k();

  Velocity2.resize(nPoint) = su2double(0.0);
  Max_Lambda_Inv.resize(nPoint) = su2double(0.0);
  Delta_Time.resize(nPoint) = su2double(0.0);
  Lambda.resize(nPoint) = su2double(0.0);
  Sensor.resize(nPoint) = su2double(0.0);

}

void CEulerVariable::SetGradient_PrimitiveZero() {
  Gradient_Primitive.storage.setConstant(0.0);
}

bool CEulerVariable::SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) {

  bool RightVol = true;

  SetVelocity(iPoint);   // Computes velocity and velocity^2
  su2double density      = GetDensity(iPoint);
  su2double staticEnergy = GetEnergy(iPoint)-0.5*Velocity2(iPoint);

  /*--- Check will be moved inside fluid model plus error description strings ---*/

  FluidModel->SetTDState_rhoe(density, staticEnergy);

  bool check_dens  = SetDensity(iPoint);
  bool check_press = SetPressure(iPoint, FluidModel->GetPressure());
  bool check_sos   = SetSoundSpeed(iPoint, FluidModel->GetSoundSpeed2());
  bool check_temp  = SetTemperature(iPoint, FluidModel->GetTemperature());

  /*--- Check that the solution has a physical meaning ---*/

  if (check_dens || check_press || check_sos || check_temp) {

    /*--- Copy the old solution ---*/

    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint, iVar) = Solution_Old(iPoint, iVar);

    /*--- Recompute the primitive variables ---*/

    SetVelocity(iPoint);   // Computes velocity and velocity^2
    su2double density = GetDensity(iPoint);
    su2double staticEnergy = GetEnergy(iPoint)-0.5*Velocity2(iPoint);
    /* check will be moved inside fluid model plus error description strings*/
    FluidModel->SetTDState_rhoe(density, staticEnergy);

    SetDensity(iPoint);
    SetPressure(iPoint, FluidModel->GetPressure());
    SetSoundSpeed(iPoint, FluidModel->GetSoundSpeed2());
    SetTemperature(iPoint, FluidModel->GetTemperature());

    RightVol = false;

  }

  SetEnthalpy(iPoint); // Requires pressure computation.

  return RightVol;
}

void CEulerVariable::SetSecondaryVar(unsigned long iPoint, CFluidModel *FluidModel) {

   /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/

   SetdPdrho_e(iPoint, FluidModel->GetdPdrho_e());
   SetdPde_rho(iPoint, FluidModel->GetdPde_rho());

}

void CEulerVariable::SetSolution_New() { Solution_New = Solution; }
