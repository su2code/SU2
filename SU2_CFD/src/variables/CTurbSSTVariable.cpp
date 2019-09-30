/*!
 * \file CTurbSSTVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
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

#include "../../include/variables/CTurbSSTVariable.hpp"


CTurbSSTVariable::CTurbSSTVariable(void) : CTurbVariable() { }

CTurbSSTVariable::CTurbSSTVariable(su2double val_kine, su2double val_omega, su2double val_muT,
                                   unsigned short val_nDim, unsigned short val_nvar, su2double *constants,
                                   CConfig *config) : CTurbVariable(val_nDim, val_nvar, config) {

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  /*--- Initialization of variables ---*/

  Solution[0] = val_kine;     Solution_Old[0] = val_kine;
  Solution[1] = val_omega;  Solution_Old[1] = val_omega;

  sigma_om2 = constants[3];
  beta_star = constants[6];

  F1   = 1.0;
  F2   = 0.0;
  CDkw = 0.0;

  /*--- Initialization of eddy viscosity ---*/

  muT = val_muT;

  /*--- Allocate and initialize solution for the dual time strategy ---*/

  if (dual_time) {
    Solution_time_n[0]  = val_kine; Solution_time_n[1]  = val_omega;
    Solution_time_n1[0]  = val_kine; Solution_time_n1[1]  = val_omega;
  }

}

CTurbSSTVariable::~CTurbSSTVariable(void) {}

void CTurbSSTVariable::SetBlendingFunc(su2double val_viscosity, su2double val_dist, su2double val_density) {
  unsigned short iDim;
  su2double arg2, arg2A, arg2B, arg1;

  AD::StartPreacc();
  AD::SetPreaccIn(val_viscosity);  AD::SetPreaccIn(val_dist);
  AD::SetPreaccIn(val_density);
  AD::SetPreaccIn(Solution, nVar);
  AD::SetPreaccIn(Gradient, nVar, nDim);

  /*--- Cross diffusion ---*/

  CDkw = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    CDkw += Gradient[0][iDim]*Gradient[1][iDim];
  CDkw *= 2.0*val_density*sigma_om2/Solution[1];
  CDkw = max(CDkw, pow(10.0, -20.0));

  /*--- F1 ---*/

  arg2A = sqrt(Solution[0])/(beta_star*Solution[1]*val_dist+EPS*EPS);
  arg2B = 500.0*val_viscosity / (val_density*val_dist*val_dist*Solution[1]+EPS*EPS);
  arg2 = max(arg2A, arg2B);
  arg1 = min(arg2, 4.0*val_density*sigma_om2*Solution[0] / (CDkw*val_dist*val_dist+EPS*EPS));
  F1 = tanh(pow(arg1, 4.0));

  /*--- F2 ---*/

  arg2 = max(2.0*arg2A, arg2B);
  F2 = tanh(pow(arg2, 2.0));

  AD::SetPreaccOut(F1); AD::SetPreaccOut(F2); AD::SetPreaccOut(CDkw);
  AD::EndPreacc();

}
