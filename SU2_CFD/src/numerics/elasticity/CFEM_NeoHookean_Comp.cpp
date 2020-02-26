/*!
 * \file CFEM_NeoHookean_Comp.cpp
 * \brief Definition of Neo-Hookean compressible material.
 * \author R. Sanchez
 * \version 7.0.1 "Blackbird"
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

#include "../../../include/numerics/elasticity/CFEM_NeoHookean_Comp.hpp"


CFEM_NeoHookean_Comp::CFEM_NeoHookean_Comp(unsigned short val_nDim,
                                           unsigned short val_nVar,
                                           CConfig *config) :
                                           CFEANonlinearElasticity(val_nDim, val_nVar, config) {
}

void CFEM_NeoHookean_Comp::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

  su2double j_red = 1.0;
  su2double fx = 0.0, fpx = 1.0;
  su2double xkm1 = 1.0, xk = 1.0;
  su2double cte = 0.0;

  unsigned short iNR, nNR;
  su2double NRTOL;

  // Maximum number of iterations and tolerance (relative)
  nNR = 10;
  NRTOL = 1E-25;

  // j_red: reduced jacobian, for the 2x2 submatrix of F
  j_red = F_Mat[0][0] * F_Mat[1][1] - F_Mat[1][0] * F_Mat[0][1];
  // cte: constant term in the NR method
  cte = Lambda*log(j_red) - Mu;

  // f(f33)  = mu*f33^2 + lambda*ln(f33) + (lambda*ln(j_red)-mu) = 0
  // f'(f33) = 2*mu*f33 + lambda/f33

  for (iNR = 0; iNR < nNR; iNR++) {
    fx  = Mu*pow(xk,2.0) + Lambda*log(xk) + cte;
    fpx = 2*Mu*xk + (Lambda / xk);
    xkm1 = xk - fx / fpx;
    if (((xkm1 - xk) / xk) < NRTOL) break;
    xk = xkm1;
  }

  f33 = xkm1;

}

void CFEM_NeoHookean_Comp::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

  su2double Mu_p = 0.0, Lambda_p = 0.0;

  /*--- This can be done in a better way ---*/
  if (J_F != 0.0) {
    Mu_p = (Mu - Lambda*log(J_F))/J_F;
    Lambda_p = Lambda/J_F;
  }

  /*--- Assuming plane strain ---*/

  su2double Lbd_2Mu = Lambda_p + 2.0 * Mu_p;

  if (nDim == 2) {
    D_Mat[0][0] = Lbd_2Mu;   D_Mat[0][1] = Lambda_p;  D_Mat[0][2] = 0.0;
    D_Mat[1][0] = Lambda_p;  D_Mat[1][1] = Lbd_2Mu;   D_Mat[1][2] = 0.0;
    D_Mat[2][0] = 0.0;       D_Mat[2][1] = 0.0;       D_Mat[2][2] = Mu_p;
  }
  else {
    D_Mat[0][0] = Lbd_2Mu;    D_Mat[0][1] = Lambda_p;   D_Mat[0][2] = Lambda_p;   D_Mat[0][3] = 0.0;   D_Mat[0][4] = 0.0;   D_Mat[0][5] = 0.0;
    D_Mat[1][0] = Lambda_p;   D_Mat[1][1] = Lbd_2Mu;    D_Mat[1][2] = Lambda_p;   D_Mat[1][3] = 0.0;   D_Mat[1][4] = 0.0;   D_Mat[1][5] = 0.0;
    D_Mat[2][0] = Lambda_p;   D_Mat[2][1] = Lambda_p;   D_Mat[2][2] = Lbd_2Mu;    D_Mat[2][3] = 0.0;   D_Mat[2][4] = 0.0;   D_Mat[2][5] = 0.0;
    D_Mat[3][0] = 0.0;        D_Mat[3][1] = 0.0;        D_Mat[3][2] = 0.0;        D_Mat[3][3] = Mu_p;  D_Mat[3][4] = 0.0;   D_Mat[3][5] = 0.0;
    D_Mat[4][0] = 0.0;        D_Mat[4][1] = 0.0;        D_Mat[4][2] = 0.0;        D_Mat[4][3] = 0.0;   D_Mat[4][4] = Mu_p;  D_Mat[4][5] = 0.0;
    D_Mat[5][0] = 0.0;        D_Mat[5][1] = 0.0;        D_Mat[5][2] = 0.0;        D_Mat[5][3] = 0.0;   D_Mat[5][4] = 0.0;   D_Mat[5][5] = Mu_p;
  }

}

void CFEM_NeoHookean_Comp::Compute_Stress_Tensor(CElement *element, CConfig *config) {

  unsigned short iVar,jVar;
  su2double Mu_J = 0.0, Lambda_J = 0.0;

  /*--- This can be done in a better way ---*/
  if (J_F != 0.0) {
    Mu_J = Mu/J_F;
    Lambda_J = Lambda/J_F;
  }

  for (iVar = 0; iVar < 3; iVar++) {
    for (jVar = 0; jVar < 3; jVar++) {
      su2double dij = deltaij(iVar,jVar);
      Stress_Tensor[iVar][jVar] = Mu_J * (b_Mat[iVar][jVar] - dij) + Lambda_J * log(J_F) * dij;
    }
  }

}

