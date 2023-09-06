/*!
 * \file CIncEulerVariable.cpp
 * \brief Definition of the variable classes for incompressible flow.
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

#include "../../include/variables/CIncEulerVariable.hpp"
#include "../../include/fluid/CFluidModel.hpp"

CIncEulerVariable::CIncEulerVariable(su2double pressure, const su2double *velocity, su2double temperature,
                                     unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config)
  : CFlowVariable(npoint, ndim, nvar, ndim + 9, ndim + 4, config),
    indices(ndim, 0) {

  const bool dual_time = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                         (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool classical_rk4 = (config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT);

  /*--- Solution initialization ---*/

  su2double val_solution[5] = {pressure, velocity[0], velocity[1], temperature, temperature};
  if(nDim==3) val_solution[3] = velocity[2];

  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = val_solution[iVar];

  Solution_Old = Solution;

  if (classical_rk4) Solution_New = Solution;

  /*--- Allocate and initialize solution for dual time strategy ---*/

  if (dual_time) {
    Solution_time_n = Solution;
    Solution_time_n1 = Solution;
  }

  if (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE) {
    Streamwise_Periodic_RecoveredPressure.resize(nPoint) = su2double(0.0);
    if (config->GetStreamwise_Periodic_Temperature())
      Streamwise_Periodic_RecoveredTemperature.resize(nPoint) = su2double(0.0);
  }
}

bool CIncEulerVariable::SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) {

  bool physical = true;

  /*--- Set the value of the pressure ---*/

  SetPressure(iPoint);

  /*--- Set the value of the temperature directly ---*/

  su2double Temperature = Solution(iPoint, nDim+1);
  const auto check_temp = SetTemperature(iPoint, Temperature);

  /*--- Use the fluid model to compute the new value of density.
  Note that the thermodynamic pressure is constant and decoupled
  from the dynamic pressure being iterated. ---*/

  /*--- Use the fluid model to compute the new value of density. ---*/

  FluidModel->SetTDState_T(Temperature);

  /*--- Set the value of the density ---*/

  const auto check_dens = SetDensity(iPoint, FluidModel->GetDensity());

  /*--- Non-physical solution found. Revert to old values. ---*/

  if (check_dens || check_temp) {

    /*--- Copy the old solution ---*/

    for (auto iVar = 0ul; iVar < nVar; iVar++)
      Solution(iPoint, iVar) = Solution_Old(iPoint, iVar);

    /*--- Recompute the primitive variables ---*/

    Temperature = Solution(iPoint, nDim+1);
    SetTemperature(iPoint, Temperature);
    FluidModel->SetTDState_T(Temperature);
    SetDensity(iPoint, FluidModel->GetDensity());

    /*--- Flag this point as non-physical. ---*/

    physical = false;

  }

  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/

  SetVelocity(iPoint);

  /*--- Set specific heats (only necessary for consistency with preconditioning). ---*/

  SetSpecificHeatCp(iPoint, FluidModel->GetCp());
  SetSpecificHeatCv(iPoint, FluidModel->GetCv());

  return physical;

}
