/*!
 * \file CFEM_DielectricElastomer.cpp
 * \brief This file contains the routines for setting the tangent matrix and residual
 *        of a FEM nonlinear elastic structural problem.
 * \author R. Sanchez
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/elasticity/CFEM_DielectricElastomer.hpp"


CFEM_DielectricElastomer::CFEM_DielectricElastomer(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEANonlinearElasticity(val_nDim, val_nVar, config) {}


void CFEM_DielectricElastomer::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

  /*--- This reduces performance by now, but it is temporal ---*/

  if (nDim == 2){
    D_Mat[0][0] = 0.0;  D_Mat[0][1] = 0.0;  D_Mat[0][2] = 0.0;
    D_Mat[1][0] = 0.0;  D_Mat[1][1] = 0.0;  D_Mat[1][2] = 0.0;
    D_Mat[2][0] = 0.0;  D_Mat[2][1] = 0.0;  D_Mat[2][2] = 0.0;
  }
  else {
    D_Mat[0][0] = 0.0;  D_Mat[0][1] = 0.0;  D_Mat[0][2] = 0.0;  D_Mat[0][3] = 0.0;  D_Mat[0][4] = 0.0;  D_Mat[0][5] = 0.0;
    D_Mat[1][0] = 0.0;  D_Mat[1][1] = 0.0;  D_Mat[1][2] = 0.0;  D_Mat[1][3] = 0.0;  D_Mat[1][4] = 0.0;  D_Mat[1][5] = 0.0;
    D_Mat[2][0] = 0.0;  D_Mat[2][1] = 0.0;  D_Mat[2][2] = 0.0;  D_Mat[2][3] = 0.0;  D_Mat[2][4] = 0.0;  D_Mat[2][5] = 0.0;
    D_Mat[3][0] = 0.0;  D_Mat[3][1] = 0.0;  D_Mat[3][2] = 0.0;  D_Mat[3][3] = 0.0;  D_Mat[3][4] = 0.0;  D_Mat[3][5] = 0.0;
    D_Mat[4][0] = 0.0;  D_Mat[4][1] = 0.0;  D_Mat[4][2] = 0.0;  D_Mat[4][3] = 0.0;  D_Mat[4][4] = 0.0;  D_Mat[4][5] = 0.0;
    D_Mat[5][0] = 0.0;  D_Mat[5][1] = 0.0;  D_Mat[5][2] = 0.0;  D_Mat[5][3] = 0.0;  D_Mat[5][4] = 0.0;  D_Mat[5][5] = 0.0;
  }


}

void CFEM_DielectricElastomer::Compute_Stress_Tensor(CElement *element, CConfig *config) {

  unsigned short iDim, jDim;

  su2double E0 = 0.0, E1 = 0.0, E2 = 0.0;
  su2double E0_2 = 0.0, E1_2 = 0.0, E2_2 = 0.0;
  su2double E_2 = 0.0;

  Compute_FmT_Mat();

  for (iDim = 0; iDim < nDim; iDim++){
    EField_Curr_Unit[iDim] = 0.0;
    for (jDim = 0; jDim < nDim; jDim++){
      EField_Curr_Unit[iDim] += FmT_Mat[iDim][jDim] * EField_Ref_Unit[jDim];
    }
  }

  E0 = EFieldMod_Ref*EField_Curr_Unit[0];          E0_2 = pow(E0,2);
  E1 = EFieldMod_Ref*EField_Curr_Unit[1];          E1_2 = pow(E1,2);
  if (nDim == 3) {E2 = EFieldMod_Ref*EField_Curr_Unit[2];  E2_2 = pow(E2,2);}

  E_2 = E0_2+E1_2+E2_2;

  Stress_Tensor[0][0] = ke_DE*(E0_2-0.5*E_2);  Stress_Tensor[0][1] = ke_DE*E0*E1;      Stress_Tensor[0][2] = ke_DE*E0*E2;
  Stress_Tensor[1][0] = ke_DE*E1*E0;      Stress_Tensor[1][1] = ke_DE*(E1_2-0.5*E_2);  Stress_Tensor[1][2] = ke_DE*E1*E2;
  Stress_Tensor[2][0] = ke_DE*E2*E0;      Stress_Tensor[2][1] = ke_DE*E2*E1;      Stress_Tensor[2][2] = ke_DE*(E2_2-0.5*E_2);


}

