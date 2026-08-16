/*!
 * \file CIncNSVariable.cpp
 * \brief Definition of the variable classes for incompressible flow.
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

#include "../../include/variables/CIncNSVariable.hpp"
#include <type_traits>
#include "../../include/fluid/CFluidModel.hpp"

/*--- Template instantiations of PB and DB versions ---*/
template class CIncNSVariable<CPBIncEulerVariable>;
template class CIncNSVariable<CDBIncEulerVariable>;

template<class CIncEulerVariable>
CIncNSVariable<CIncEulerVariable>::CIncNSVariable(su2double density, su2double pressure, const su2double *velocity, su2double enthalpy,
                               unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config) :
                               CIncEulerVariable(density, pressure, velocity, enthalpy, npoint, ndim, nvar, config),
                               Energy(config->GetEnergy_Equation()) {

  this->Vorticity.resize(this->nPoint,3);
  this->StrainMag.resize(this->nPoint);
  Tau_Wall.resize(this->nPoint) = su2double(-1.0);
  DES_LengthScale.resize(this->nPoint) = su2double(0.0);
  lesMode.resize(this->nPoint) = su2double(0.0);
  this->Max_Lambda_Visc.resize(this->nPoint);
  /*--- Allocate memory for the AuxVar and its gradient. See e.g. CIncEulerSolver::Source_Residual:
   * Axisymmetric: total-viscosity * y-vel / y-coord
   * Streamwise Periodic: eddy viscosity (mu_t) ---*/
  if (config->GetAxisymmetric() ||
      (config->GetStreamwise_Periodic_Temperature() && (config->GetKind_Turb_Model() != TURB_MODEL::NONE))) {
    this->nAuxVar = 1;
    this->AuxVar.resize(this->nPoint,this->nAuxVar) = su2double(0.0);
    this->Grad_AuxVar.resize(this->nPoint,this->nAuxVar,this->nDim);
  }

  /*--- Check what kind of solver is instantiated ---*/

  pressure_based = std::is_same_v<CIncEulerVariable, CPBIncEulerVariable>;

}

template<class CIncEulerVariable>
bool CIncNSVariable<CIncEulerVariable>::SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel, const su2double *scalar) {

  bool physical = true;
  bool check_temp;
  su2double Enthalpy;

  /*--- Use the fluid model to compute the new value of density.
  Note that the thermodynamic pressure is constant and decoupled
  from the dynamic pressure being iterated. ---*/

  if (!pressure_based) {

    /*--- Set the value of the pressure ---*/

    this->SetPressure(iPoint);

    /*--- Set the value of the temperature ---*/

    Enthalpy = this->Solution(iPoint, this->nDim + 1);
    FluidModel->SetTDState_h(Enthalpy, scalar);
    su2double Temperature = FluidModel->GetTemperature();

    check_temp = this->SetTemperature(iPoint, Temperature, this->TemperatureLimits);

  }

  /*--- Set the value of the density ---*/

  const auto check_dens = this->SetDensity(iPoint, FluidModel->GetDensity());

  if (pressure_based) physical = !check_dens;
  else physical = !(check_dens || check_temp);

  /*--- Non-physical solution found. Revert to old values. ---*/

  if (!physical) {

    /*--- Copy the old solution ---*/

    for (auto iVar = 0ul; iVar < this->nVar; iVar++)
      this->Solution(iPoint,iVar) = this->Solution_Old(iPoint,iVar);

    /*--- Recompute the primitive variables ---*/

    if (!pressure_based) {
      Enthalpy = this->Solution(iPoint,this->nDim + 1);
      FluidModel->SetTDState_h(Enthalpy, scalar);
      this->SetTemperature(iPoint, FluidModel->GetTemperature(), this->TemperatureLimits);
    }
    this->SetDensity(iPoint, FluidModel->GetDensity());

  }

  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/

  this->SetVelocity(iPoint);

  /*--- Set laminar viscosity ---*/

  SetLaminarViscosity(iPoint, FluidModel->GetLaminarViscosity());

  /*--- Set eddy viscosity locally and in the fluid model. ---*/

  SetEddyViscosity(iPoint, eddy_visc);
  FluidModel->SetEddyViscosity(eddy_visc);

  /*--- Set thermal conductivity (effective value if RANS). ---*/

  SetThermalConductivity(iPoint, FluidModel->GetThermalConductivity());

  /*--- Set specific heats ---*/

  this->SetSpecificHeatCp(iPoint, FluidModel->GetCp());
  this->SetSpecificHeatCv(iPoint, FluidModel->GetCv());

  /*--- Set enthalpy ---*/

  this->SetEnthalpy(iPoint, FluidModel->GetEnthalpy());

  return physical;

}
