/*!
 * \file nonlinear_models.cpp
 * \brief Definition of nonlinear constitutive models.
 * \author R. Sanchez
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

#include "../../../include/numerics/elasticity/nonlinear_models.hpp"


CFEM_NeoHookean_Comp::CFEM_NeoHookean_Comp(unsigned short val_nDim,
                                           unsigned short val_nVar,
                                           const CConfig *config) :
                                           CFEANonlinearElasticity(val_nDim, val_nVar, config) {
}

void CFEM_NeoHookean_Comp::Compute_Plane_Stress_Term(CElement *element, const CConfig *config) {

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

void CFEM_NeoHookean_Comp::Compute_Constitutive_Matrix(CElement *element, const CConfig *config) {

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

void CFEM_NeoHookean_Comp::Compute_Stress_Tensor(CElement *element, const CConfig *config) {

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

CFEM_Knowles_NearInc::CFEM_Knowles_NearInc(unsigned short val_nDim, unsigned short val_nVar,
                        const CConfig *config) : CFEANonlinearElasticity(val_nDim, val_nVar, config) {

  /* -- The formulation adopted for this material model has been described by:
   * --
   * -- Suchocki, C., A Finite Element Implementation of Knowles stored-energy function:
   * -- theory, coding and applications, Archive of Mechanical Engineering, Vol. 58, pp. 319-346 (2011).
   * --
   * -- DOI: 10.2478/v10180-011-0021-7
   */

  Bk = config->GetKnowles_B();
  Nk = config->GetKnowles_N();

  trbbar = 0.0;
  Ek     = 0.0;
  Pr     = 0.0;

}

void CFEM_Knowles_NearInc::Compute_Plane_Stress_Term(CElement *element, const CConfig *config) {

  SU2_MPI::Error("This material model cannot (yet) be used for plane stress.",CURRENT_FUNCTION);

}

void CFEM_Knowles_NearInc::Compute_Constitutive_Matrix(CElement *element, const CConfig *config) {

  /* -- Suchocki (2011) (full reference in class constructor). ---*/

  /*--- Computation of the tensor cijkl ---*/

  unsigned short iVar, jVar, kVar, lVar;

  for (iVar = 0; iVar < 3; iVar++){
    for (jVar = 0; jVar < 3; jVar++){
      for (kVar = 0; kVar < 3; kVar++){
        for (lVar = 0; lVar < 3; lVar++){
          cijkl[iVar][jVar][kVar][lVar] =
            term1 * ((1.0/2.0)*( deltaij(iVar,kVar)*b_Mat_Iso[jVar][lVar]
                      +deltaij(jVar,lVar)*b_Mat_Iso[iVar][kVar]
                      +deltaij(iVar,lVar)*b_Mat_Iso[jVar][kVar]
                      +deltaij(jVar,kVar)*b_Mat_Iso[iVar][lVar])
                 +(2.0/3.0)*(trbbar*deltaij(iVar,jVar)*deltaij(kVar,lVar)
                      -b_Mat_Iso[iVar][jVar]*deltaij(kVar,lVar)
                      -deltaij(iVar,jVar)*b_Mat_Iso[kVar][lVar]))
             +term2 * ( b_Mat_Iso[iVar][jVar]*b_Mat_Iso[kVar][lVar]
                - trbbar*(b_Mat_Iso[iVar][jVar]*deltaij(kVar,lVar)
                     +deltaij(iVar,jVar)*b_Mat_Iso[kVar][lVar])
                + pow(trbbar,2.0) * deltaij(iVar,jVar) * deltaij(kVar,lVar))
             +Kappa * (2.0 * J_F - 1.0) * deltaij(iVar,jVar) * deltaij(kVar,lVar);

        }
      }
    }
  }

  /*--- Reorganizing the tensor into the matrix D ---*/

  Assign_cijkl_D_Mat();


}

void CFEM_Knowles_NearInc::Compute_Stress_Tensor(CElement *element, const CConfig *config) {

  /* -- Suchocki (2011) (full reference in class constructor). ---*/

  unsigned short iVar, jVar;

  /*--- Compute the isochoric deformation gradient Fbar and left Cauchy-Green tensor bbar ---*/
  Compute_Isochoric_F_b();

  trbbar = (b_Mat_Iso[0][0] + b_Mat_Iso[1][1] + b_Mat_Iso[2][2]) / 3.0;
  term1 = (Mu / J_F) * pow((1 + (Bk / Nk) * (3.0 * trbbar - 3.0)), (Nk-1.0));
  term2 = 2.0 * (Mu / J_F) * (Bk * (Nk - 1.0) / Nk) *
      pow((1.0 + (Bk / Nk) * (3.0 * trbbar - 3.0)), (Nk-2.0));

  Ek = Kappa * (2.0 * J_F - 1.0);
  Pr = Kappa * (J_F - 1.0);

  for (iVar = 0; iVar < 3; iVar++){
    for (jVar = 0; jVar < 3; jVar++){
      Stress_Tensor[iVar][jVar] = term1 * (b_Mat_Iso[iVar][jVar] - deltaij(iVar,jVar)*trbbar) +
                                  deltaij(iVar,jVar) * Pr;
    }
  }

}

CFEM_DielectricElastomer::CFEM_DielectricElastomer(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config) :
                          CFEANonlinearElasticity(val_nDim, val_nVar, config) {

}

void CFEM_DielectricElastomer::Compute_Constitutive_Matrix(CElement *element, const CConfig *config) {

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

void CFEM_DielectricElastomer::Compute_Stress_Tensor(CElement *element, const CConfig *config) {

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

CFEM_IdealDE::CFEM_IdealDE(unsigned short val_nDim, unsigned short val_nVar,
                           const CConfig *config) : CFEANonlinearElasticity(val_nDim, val_nVar, config) {

  /* -- The formulation adopted for this material model has been described by:
   * --
   * -- Zhao, X. and Suo, Z., Method to analyze programmable deformation of
   * -- dielectric elastomer layers, Applied Physics Letters 93, 251902 (2008).
   * --
   * -- http://imechanica.org/node/4234
   */

  trbbar = 0.0;
  Eg     = 0.0;
  Eg23   = 0.0;
  Ek     = 0.0;
  Pr     = 0.0;

}

void CFEM_IdealDE::Compute_Plane_Stress_Term(CElement *element, const CConfig *config) {

  SU2_MPI::Error("This material model cannot (yet) be used for plane stress.", CURRENT_FUNCTION);

}

void CFEM_IdealDE::Compute_Constitutive_Matrix(CElement *element, const CConfig *config) {

  /* -- Zhao, X. and Suo, Z. (2008) (full reference in class constructor). ---*/

  if (nDim == 2){

    D_Mat[0][0] = Eg23*(b_Mat_Iso[0][0]+trbbar)+Ek;
    D_Mat[1][1] = Eg23*(b_Mat_Iso[1][1]+trbbar)+Ek;

    D_Mat[0][1] = -Eg23*(b_Mat_Iso[0][0]+b_Mat_Iso[1][1]-trbbar)+Ek;
    D_Mat[1][0] = -Eg23*(b_Mat_Iso[0][0]+b_Mat_Iso[1][1]-trbbar)+Ek;

    D_Mat[0][2] = Eg23*b_Mat_Iso[0][1]/2.0;
    D_Mat[2][0] = Eg23*b_Mat_Iso[0][1]/2.0;

    D_Mat[1][2] = Eg23*b_Mat_Iso[0][1]/2.0;
    D_Mat[2][1] = Eg23*b_Mat_Iso[0][1]/2.0;

    D_Mat[2][2] = Eg*(b_Mat_Iso[0][0]+b_Mat_Iso[1][1])/2.0;

  }
  else {
    SU2_MPI::Error("This material model cannot be used for 3D yet.", CURRENT_FUNCTION);
  }

}

void CFEM_IdealDE::Compute_Stress_Tensor(CElement *element, const CConfig *config) {

  /* -- Zhao, X. and Suo, Z. (2008) (full reference in class constructor). ---*/

  unsigned short iVar, jVar;
  su2double dij = 0.0;

  /*--- Compute the isochoric deformation gradient Fbar and left Cauchy-Green tensor bbar ---*/
  Compute_Isochoric_F_b();

  // Stress terms

  trbbar = (b_Mat_Iso[0][0] + b_Mat_Iso[1][1] + b_Mat_Iso[2][2]) / 3.0;
  Eg = Mu / J_F;
  Ek = Kappa * (2.0 * J_F - 1.0);
  Pr = Kappa * (J_F - 1.0);
  Eg23 = 2.0 * Eg / 3.0;

  // Stress tensor

  for (iVar = 0; iVar < 3; iVar++){
    for (jVar = 0; jVar < 3; jVar++){
      if (iVar == jVar) dij = 1.0;
      else if (iVar != jVar) dij = 0.0;
      Stress_Tensor[iVar][jVar] = Eg * ( b_Mat_Iso[iVar][jVar] - dij * trbbar) + dij * Pr ;
    }
  }

}
