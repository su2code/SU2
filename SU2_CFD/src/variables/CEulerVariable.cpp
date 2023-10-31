/*!
 * \file CEulerVariable.cpp
 * \brief Definition of the solution fields.
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

#include "../../include/variables/CEulerVariable.hpp"
#include "../../include/fluid/CFluidModel.hpp"

CEulerVariable::CEulerVariable(su2double density, const su2double *velocity, su2double energy, unsigned long npoint,
                               unsigned long ndim, unsigned long nvar, const CConfig *config)
  : CFlowVariable(npoint, ndim, nvar, ndim + 9, ndim + 4, config),
    indices(ndim, 0) {

  const bool dual_time = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                         (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool classical_rk4 = (config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT);

  nSecondaryVar = config->GetViscous() ? 8 : 2,
  nSecondaryVarGrad = 2;

  /*--- Solution initialization ---*/

  su2double val_solution[5] = {su2double(1.0), velocity[0], velocity[1], energy, energy};
  if(nDim==3) val_solution[3] = velocity[2];

  for(unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; ++iVar)
      Solution(iPoint,iVar) = density*val_solution[iVar];

  Solution_Old = Solution;

  if (classical_rk4) Solution_New = Solution;

  /*--- Allocate and initializate solution for dual time strategy ---*/

  if (dual_time) {
    Solution_time_n = Solution;
    Solution_time_n1 = Solution;
  }

  Secondary.resize(nPoint,nSecondaryVar) = su2double(0.0);

  if (config->GetAxisymmetric()){
    nAuxVar = 3;
    Grad_AuxVar.resize(nPoint,nAuxVar,nDim,0.0);
    AuxVar.resize(nPoint,nAuxVar) = su2double(0.0);
  }

  if (config->GetWind_Gust()) {
    WindGust.resize(nPoint,nDim);
  }

  if (config->GetVorticityConfinement()) {
    nAuxVar = 1;
    Grad_AuxVar.resize(nPoint, nAuxVar, nDim, 0.0);
    AuxVar.resize(nPoint, nAuxVar) = su2double(0.0);
  }
  
  if (config->GetKind_FluidModel() == ENUM_FLUIDMODEL::DATADRIVEN_FLUID){
    DataDrivenFluid = true;
    DatasetExtrapolation.resize(nPoint) = 0;
    NIterNewtonsolver.resize(nPoint) = 0;
    FluidEntropy.resize(nPoint) = su2double(0.0);
  }
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

  /*--- Set look-up variables in case of data-driven fluid model ---*/
  if (DataDrivenFluid) {
    SetDataExtrapolation(iPoint, FluidModel->GetExtrapolation());
    SetEntropy(iPoint, FluidModel->GetEntropy());
  }

  return RightVol;
}

void CEulerVariable::SetSecondaryVar(unsigned long iPoint, CFluidModel *FluidModel) {

   /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/

   SetdPdrho_e(iPoint, FluidModel->GetdPdrho_e());
   SetdPde_rho(iPoint, FluidModel->GetdPde_rho());

}
