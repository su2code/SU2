/*!
 * \file CPBIncEulerVariable.cpp
 * \brief Definition of the variable classes for pressure based incompressible flow.
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

#include "../../include/variables/CPBIncEulerVariable.hpp"
#include "../../include/fluid/CFluidModel.hpp"

CPBIncEulerVariable::CPBIncEulerVariable(su2double density, su2double pressure, const su2double *velocity, su2double temperature,
                                     unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config)
  : CFlowVariable(npoint, ndim, nvar, ndim + 10,
                  ndim + (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED ? 2 : 4), config),
    indices(ndim, 0) {

  const bool dual_time = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                         (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool classical_rk4 = (config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT);
  TemperatureLimits[0]= config->GetTemperatureLimits(0);
  TemperatureLimits[1]= config->GetTemperatureLimits(1);
  
  /*--- Solution initialization ---*/
  for(unsigned long iPoint=0; iPoint < nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = density*velocity[iVar];

  Solution_Old = Solution;

  if (classical_rk4) Solution_New = Solution;

  /*--- Allocate and initialize solution for dual time strategy ---*/

  if (dual_time) {
    Solution_time_n = Solution;
    Solution_time_n1 = Solution;

    if (config->GetKind_DensityModel() != INC_DENSITYMODEL::CONSTANT) {
      Density_time_n.resize(nPoint) = su2double(0.0);
      Density_time_n1.resize(nPoint) = su2double(0.0);
    }
  }

  /*--- Initialize boolean flag for BC ---*/
  strongBC.resize(nPoint) = false;

}

bool CPBIncEulerVariable::SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) {

  bool physical = true;

  /*--- Set the value of the temperature ---*/

  SetTemperature(iPoint, FluidModel->GetTemperature(), TemperatureLimits);

  /*--- Set the value of the density ---*/
  
  const auto check_dens = SetDensity(iPoint, FluidModel->GetDensity());

  if (check_dens) physical = false;

  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/

  SetVelocity(iPoint);

  /*--- Set specific heats ---*/

  SetSpecificHeatCp(iPoint, FluidModel->GetCp());
  SetSpecificHeatCv(iPoint, FluidModel->GetCv());

  /*--- Set enthalpy ---*/

  SetEnthalpy(iPoint, FluidModel->GetEnthalpy());

  return physical;


}
