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

#include "../../include/variables/CPBIncNSVariable.hpp"


CPBIncNSVariable::CPBIncNSVariable(void) : CPBIncEulerVariable() { }

CPBIncNSVariable::CPBIncNSVariable(su2double val_pressure, su2double *val_velocity, 
                         unsigned short val_nDim, unsigned short val_nvar,
                         CConfig *config) : CPBIncEulerVariable(val_pressure, val_velocity, val_nDim, val_nvar, config) {
 
}

CPBIncNSVariable::CPBIncNSVariable(su2double *val_solution, unsigned short val_nDim,
                         unsigned short val_nvar, CConfig *config) : CPBIncEulerVariable(val_solution, val_nDim, val_nvar, config) {
  
}

CPBIncNSVariable::~CPBIncNSVariable(void) { }

bool CPBIncNSVariable::SetVorticity(void) {
  
  Vorticity[0] = 0.0; Vorticity[1] = 0.0;
  
  Vorticity[2] = Gradient_Primitive[2][0]-Gradient_Primitive[1][1];
  
  if (nDim == 3) {
    Vorticity[0] = Gradient_Primitive[3][1]-Gradient_Primitive[2][2];
    Vorticity[1] = -(Gradient_Primitive[3][0]-Gradient_Primitive[1][2]);
  }
  
  return false;
  
}

bool CPBIncNSVariable::SetStrainMag(void) {
  
  su2double Div;
  unsigned short iDim;
    
  Div = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Div += Gradient_Primitive[iDim+1][iDim];
  }
  
  StrainMag = 0.0;
  
  /*--- Add diagonal part ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    StrainMag += pow(Gradient_Primitive[iDim+1][iDim] - 1.0/3.0*Div, 2.0);
  }
  if (nDim == 2) {
    StrainMag += pow(1.0/3.0*Div, 2.0);
  }
  /*--- Add off diagonals ---*/

  StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][1] + Gradient_Primitive[2][0]), 2.0);

  if (nDim == 3) {
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][2] + Gradient_Primitive[3][0]), 2.0);
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[2][2] + Gradient_Primitive[3][1]), 2.0);
  }
  
  StrainMag = sqrt(2.0*StrainMag);

  return false;

}


bool CPBIncNSVariable::SetPrimVar(su2double Density_Inf, su2double Viscosity_Inf, su2double eddy_visc, su2double turb_ke, CConfig *config) {
      
  unsigned short iVar;
  bool physical = true;

  /*--- Set the value of the density ---*/
  
  SetDensity(Density_Inf);
  
  /*--- Set the value of the velocity and velocity^2 (requires density) ---*/
  
  SetVelocity();

  /*--- Set laminar viscosity ---*/
  
  SetLaminarViscosity(Viscosity_Inf);
  
  /*--- Set eddy viscosity locally and in the fluid model. ---*/
  SetEddyViscosity(eddy_visc);
  
  return physical; 
}


