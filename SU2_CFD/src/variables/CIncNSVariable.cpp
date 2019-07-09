/*!
 * \file CIncNSVariable.cpp
 * \brief Definition of the variable classes for incompressible flow.
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

#include "../../include/variables/CIncNSVariable.hpp"


CIncNSVariable::CIncNSVariable(void) : CIncEulerVariable() { }

CIncNSVariable::CIncNSVariable(su2double val_pressure, su2double *val_velocity, su2double val_temperature,
                               unsigned short val_nDim, unsigned short val_nvar, CConfig *config) :
                               CIncEulerVariable(val_pressure, val_velocity, val_temperature, val_nDim, val_nvar, config) {
  DES_LengthScale = 0.0;
}

CIncNSVariable::CIncNSVariable(su2double *val_solution, unsigned short val_nDim, unsigned short val_nvar,
                               CConfig *config) : CIncEulerVariable(val_solution, val_nDim, val_nvar, config) {
  DES_LengthScale = 0.0;
}

CIncNSVariable::~CIncNSVariable(void) { }

bool CIncNSVariable::SetVorticity(void) {

  Vorticity[0] = 0.0; Vorticity[1] = 0.0;

  Vorticity[2] = Gradient_Primitive[2][0]-Gradient_Primitive[1][1];

  if (nDim == 3) {
    Vorticity[0] = Gradient_Primitive[3][1]-Gradient_Primitive[2][2];
    Vorticity[1] = -(Gradient_Primitive[3][0]-Gradient_Primitive[1][2]);
  }

  return false;

}

bool CIncNSVariable::SetStrainMag(void) {

  su2double Div;
  unsigned short iDim;

  AD::StartPreacc();
  AD::SetPreaccIn(Gradient_Primitive, nDim+1, nDim);

  Div = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Div += Gradient_Primitive[iDim+1][iDim];
  }

  StrainMag = 0.0;

  /*--- Add diagonal part ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    StrainMag += pow(Gradient_Primitive[iDim+1][iDim] - 1.0/3.0*Div, 2.0);
  }

  /*--- Add off diagonals ---*/

  StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][1] + Gradient_Primitive[2][0]), 2.0);

  if (nDim == 3) {
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][2] + Gradient_Primitive[3][0]), 2.0);
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[2][2] + Gradient_Primitive[3][1]), 2.0);
  }

  StrainMag = sqrt(2.0*StrainMag);

  AD::SetPreaccOut(StrainMag);
  AD::EndPreacc();

  return false;

}


bool CIncNSVariable::SetPrimVar(su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) {

  unsigned short iVar;
  bool check_dens = false, check_temp = false, physical = true;

  /*--- Store the density from the previous iteration. ---*/

  Density_Old = GetDensity();

  /*--- Set the value of the pressure ---*/

  SetPressure();

  /*--- Set the value of the temperature directly ---*/

  su2double Temperature = Solution[nDim+1];
  check_temp = SetTemperature(Temperature);

  /*--- Use the fluid model to compute the new value of density.
  Note that the thermodynamic pressure is constant and decoupled
  from the dynamic pressure being iterated. ---*/

  /*--- Use the fluid model to compute the new value of density. ---*/

  FluidModel->SetTDState_T(Temperature);

  /*--- Set the value of the density ---*/

  check_dens = SetDensity(FluidModel->GetDensity());

  /*--- Non-physical solution found. Revert to old values. ---*/

  if (check_dens || check_temp) {

    /*--- Copy the old solution ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];

    /*--- Recompute the primitive variables ---*/

    Temperature = Solution[nDim+1];
    SetTemperature(Temperature);
    FluidModel->SetTDState_T(Temperature);
    SetDensity(FluidModel->GetDensity());

    /*--- Flag this point as non-physical. ---*/

    physical = false;

  }

  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/

  SetVelocity();

  /*--- Set laminar viscosity ---*/

  SetLaminarViscosity(FluidModel->GetLaminarViscosity());

  /*--- Set eddy viscosity locally and in the fluid model. ---*/

  SetEddyViscosity(eddy_visc);
  FluidModel->SetEddyViscosity(eddy_visc);

  /*--- Set thermal conductivity (effective value if RANS). ---*/

  SetThermalConductivity(FluidModel->GetThermalConductivity());

  /*--- Set specific heats ---*/

  SetSpecificHeatCp(FluidModel->GetCp());
  SetSpecificHeatCv(FluidModel->GetCv());

  return physical;

}
