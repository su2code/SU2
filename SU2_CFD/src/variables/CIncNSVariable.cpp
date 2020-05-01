/*!
 * \file CIncNSVariable.cpp
 * \brief Definition of the variable classes for incompressible flow.
 * \author F. Palacios, T. Economon
 * \version 7.0.4 "Blackbird"
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


#include "../../include/variables/CIncNSVariable.hpp"


CIncNSVariable::CIncNSVariable(su2double pressure, const su2double *velocity, su2double temperature,
                               unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config) :
                               CIncEulerVariable(pressure, velocity, temperature, npoint, ndim, nvar, config) {
  Vorticity.resize(nPoint,3);
  StrainMag.resize(nPoint);
  DES_LengthScale.resize(nPoint) = su2double(0.0);
  Max_Lambda_Visc.resize(nPoint);
}

bool CIncNSVariable::SetVorticity_StrainMag() {

  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {

    /*--- Vorticity ---*/

    Vorticity(iPoint,0) = 0.0; Vorticity(iPoint,1) = 0.0;

    Vorticity(iPoint,2) = Gradient_Primitive(iPoint,2,0)-Gradient_Primitive(iPoint,1,1);

    if (nDim == 3) {
      Vorticity(iPoint,0) = Gradient_Primitive(iPoint,3,1)-Gradient_Primitive(iPoint,2,2);
      Vorticity(iPoint,1) = -(Gradient_Primitive(iPoint,3,0)-Gradient_Primitive(iPoint,1,2));
    }

    /*--- Strain Magnitude ---*/

    AD::StartPreacc();
    AD::SetPreaccIn(Gradient_Primitive[iPoint], nDim+1, nDim);

    su2double Div = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      Div += Gradient_Primitive(iPoint,iDim+1,iDim);

    StrainMag(iPoint) = 0.0;

    /*--- Add diagonal part ---*/

    for (unsigned long iDim = 0; iDim < nDim; iDim++) {
      StrainMag(iPoint) += pow(Gradient_Primitive(iPoint,iDim+1,iDim) - 1.0/3.0*Div, 2.0);
    }
    if (nDim == 2) {
      StrainMag(iPoint) += pow(1.0/3.0*Div, 2.0);
    }

    /*--- Add off diagonals ---*/

    StrainMag(iPoint) += 2.0*pow(0.5*(Gradient_Primitive(iPoint,1,1) + Gradient_Primitive(iPoint,2,0)), 2);

    if (nDim == 3) {
      StrainMag(iPoint) += 2.0*pow(0.5*(Gradient_Primitive(iPoint,1,2) + Gradient_Primitive(iPoint,3,0)), 2);
      StrainMag(iPoint) += 2.0*pow(0.5*(Gradient_Primitive(iPoint,2,2) + Gradient_Primitive(iPoint,3,1)), 2);
    }

    StrainMag(iPoint) = sqrt(2.0*StrainMag(iPoint));

    AD::SetPreaccOut(StrainMag(iPoint));
    AD::EndPreacc();
  }
  return false;
}


bool CIncNSVariable::SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) {

  unsigned short iVar;
  bool check_dens = false, check_temp = false, physical = true;

  /*--- Store the density from the previous iteration. ---*/

  Density_Old(iPoint) = GetDensity(iPoint);

  /*--- Set the value of the pressure ---*/

  SetPressure(iPoint);

  /*--- Set the value of the temperature directly ---*/

  su2double Temperature = Solution(iPoint,nDim+1);
  check_temp = SetTemperature(iPoint,Temperature);

  /*--- Use the fluid model to compute the new value of density.
  Note that the thermodynamic pressure is constant and decoupled
  from the dynamic pressure being iterated. ---*/

  /*--- Use the fluid model to compute the new value of density. ---*/

  FluidModel->SetTDState_T(Temperature);

  /*--- Set the value of the density ---*/

  check_dens = SetDensity(iPoint, FluidModel->GetDensity());

  /*--- Non-physical solution found. Revert to old values. ---*/

  if (check_dens || check_temp) {

    /*--- Copy the old solution ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = Solution_Old(iPoint,iVar);

    /*--- Recompute the primitive variables ---*/

    Temperature = Solution(iPoint,nDim+1);
    SetTemperature(iPoint, Temperature);
    FluidModel->SetTDState_T(Temperature);
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
