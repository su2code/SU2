/*!
 * \file CIncNSVariable.cpp
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

#include "../../include/variables/CIncNSVariable.hpp"
#include "../../include/fluid/CFluidModel.hpp"

CIncNSVariable::CIncNSVariable(su2double pressure, const su2double *velocity, su2double temperature,
                               unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config) :
                               CIncEulerVariable(pressure, velocity, temperature, npoint, ndim, nvar, config) {

  Vorticity.resize(nPoint,3);
  StrainMag.resize(nPoint);
  Tau_Wall.resize(nPoint) = su2double(-1.0);
  DES_LengthScale.resize(nPoint) = su2double(0.0);
  Max_Lambda_Visc.resize(nPoint);
  /*--- Allocate memory for the AuxVar and its gradient. See e.g. CIncEulerSolver::Source_Residual:
   * Axisymmetric: total-viscosity * y-vel / y-coord
   * Streamwise Periodic: eddy viscosity (mu_t) ---*/
  if (config->GetAxisymmetric() ||
      (config->GetStreamwise_Periodic_Temperature() && (config->GetKind_Turb_Model() != TURB_MODEL::NONE))) {
    nAuxVar = 1;
    AuxVar.resize(nPoint,nAuxVar) = su2double(0.0);
    Grad_AuxVar.resize(nPoint,nAuxVar,nDim);
  }
}

bool CIncNSVariable::SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel, const su2double *scalar) {

  bool physical = true;

  /*--- Set the value of the pressure ---*/

  SetPressure(iPoint);

  /*--- Set the value of the temperature directly ---*/

  su2double Temperature = Solution(iPoint, indices.Temperature());
  auto check_temp = SetTemperature(iPoint, Temperature);

  /*--- Use the fluid model to compute the new value of density.
  Note that the thermodynamic pressure is constant and decoupled
  from the dynamic pressure being iterated. ---*/

  /*--- Use the fluid model to compute the new value of density. ---*/

  FluidModel->SetTDState_T(Temperature, scalar);

  /*--- for FLAMELET: copy the LUT temperature into the solution ---*/
  Solution(iPoint,nDim+1) = FluidModel->GetTemperature();
  /*--- for FLAMELET: update the local temperature using LUT variables ---*/
  Temperature = Solution(iPoint,indices.Temperature());
  check_temp = SetTemperature(iPoint, Temperature);

  /*--- Set the value of the density ---*/

  const auto check_dens = SetDensity(iPoint, FluidModel->GetDensity());

  /*--- Non-physical solution found. Revert to old values. ---*/

  if (check_dens || check_temp) {

    /*--- Copy the old solution ---*/

    for (auto iVar = 0ul; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = Solution_Old(iPoint,iVar);

    /*--- Recompute the primitive variables ---*/

    Temperature = Solution(iPoint, indices.Temperature());
    SetTemperature(iPoint, Temperature);
    FluidModel->SetTDState_T(Temperature, scalar);
    SetDensity(iPoint, FluidModel->GetDensity());

    /*--- Flag this point as non-physical. ---*/

    physical = false;

  }

  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/

  SetVelocity(iPoint);

  /*--- Set laminar viscosity ---*/

  SetLaminarViscosity(iPoint, FluidModel->GetLaminarViscosity());

  /*--- Set eddy viscosity locally and in the fluid model. ---*/

  SetEddyViscosity(iPoint, eddy_visc);
  FluidModel->SetEddyViscosity(eddy_visc);

  /*--- Set thermal conductivity (effective value if RANS). ---*/

  SetThermalConductivity(iPoint, FluidModel->GetThermalConductivity());

  /*--- Set specific heats ---*/

  SetSpecificHeatCp(iPoint, FluidModel->GetCp());
  SetSpecificHeatCv(iPoint, FluidModel->GetCv());

  return physical;

}
