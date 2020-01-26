/*!
 * \file CSourcePieceWise_TurbSST.cpp
 * \brief Implementation of numerics class CSourcePieceWise_TurbSST.
 * \author F. Palacios, T. Economon
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

#include "../../../include/numerics/turbulent/CSourcePieceWise_TurbSST.hpp"

CSourcePieceWise_TurbSST::CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar, const su2double *constants,
                                                   su2double val_kine_Inf, su2double val_omega_Inf, CConfig *config)
  : CNumerics(val_nDim, val_nVar, config) {

  incompressible   = (config->GetKind_Regime()     == INCOMPRESSIBLE);
  sustaining_terms = (config->GetKind_Turb_Model() == SST_SUST);

  /*--- Closure constants ---*/
  beta_star     = constants[6];
  sigma_omega_1 = constants[2];
  sigma_omega_2 = constants[3];
  beta_1        = constants[4];
  beta_2        = constants[5];
  alfa_1        = constants[8];
  alfa_2        = constants[9];
  a1            = constants[7];

  /*--- Set the ambient values of k and omega to the free stream values. ---*/
  kAmb     = val_kine_Inf;
  omegaAmb = val_omega_Inf;
}

CSourcePieceWise_TurbSST::~CSourcePieceWise_TurbSST(void) { }

void CSourcePieceWise_TurbSST::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(TurbVar_i, nVar);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);
  AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F2_i); AD::SetPreaccIn(CDkw_i);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);

  unsigned short iDim;
  su2double alfa_blended, beta_blended;
  su2double diverg, pk, pw, zeta;

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+6);

    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+4];
    Eddy_Viscosity_i = V_i[nDim+5];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7);

    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];
  }

  val_residual[0] = 0.0;        val_residual[1] = 0.0;
  val_Jacobian_i[0][0] = 0.0;    val_Jacobian_i[0][1] = 0.0;
  val_Jacobian_i[1][0] = 0.0;    val_Jacobian_i[1][1] = 0.0;

  /*--- Computation of blended constants for the source terms---*/

  alfa_blended = F1_i*alfa_1 + (1.0 - F1_i)*alfa_2;
  beta_blended = F1_i*beta_1 + (1.0 - F1_i)*beta_2;

  if (dist_i > 1e-10) {

   /*--- Production ---*/

   diverg = 0.0;
   for (iDim = 0; iDim < nDim; iDim++)
     diverg += PrimVar_Grad_i[iDim+1][iDim];

   /* if using UQ methodolgy, calculate production using perturbed Reynolds stress matrix */

   if (using_uq){
     SetReynoldsStressMatrix(TurbVar_i[0]);
     SetPerturbedRSM(TurbVar_i[0], config);
     SetPerturbedStrainMag(TurbVar_i[0]);
     pk = Eddy_Viscosity_i*PerturbedStrainMag*PerturbedStrainMag
     - 2.0/3.0*Density_i*TurbVar_i[0]*diverg;
   }
   else {
     pk = Eddy_Viscosity_i*StrainMag_i*StrainMag_i - 2.0/3.0*Density_i*TurbVar_i[0]*diverg;
   }


   pk = min(pk,20.0*beta_star*Density_i*TurbVar_i[1]*TurbVar_i[0]);
   pk = max(pk,0.0);

   zeta = max(TurbVar_i[1], StrainMag_i*F2_i/a1);

   /* if using UQ methodolgy, calculate production using perturbed Reynolds stress matrix */

   if (using_uq){
    pw = PerturbedStrainMag * PerturbedStrainMag - 2.0/3.0*zeta*diverg;
   }
   else {
     pw = StrainMag_i*StrainMag_i - 2.0/3.0*zeta*diverg;
   }
   pw = alfa_blended*Density_i*max(pw,0.0);

   /*--- Sustaining terms, if desired. Note that if the production terms are
         larger equal than the sustaining terms, the original formulation is
         obtained again. This is in contrast to the version in literature
         where the sustaining terms are simply added. This latter approach could
         lead to problems for very big values of the free-stream turbulence
         intensity. ---*/

   if ( sustaining_terms ) {
     const su2double sust_k = beta_star*Density_i*kAmb*omegaAmb;
     const su2double sust_w = beta_blended*Density_i*omegaAmb*omegaAmb;

     pk = max(pk, sust_k);
     pw = max(pw, sust_w);
   }

   /*--- Add the production terms to the residuals. ---*/

   val_residual[0] += pk*Volume;
   val_residual[1] += pw*Volume;

   /*--- Dissipation ---*/

   val_residual[0] -= beta_star*Density_i*TurbVar_i[1]*TurbVar_i[0]*Volume;
   val_residual[1] -= beta_blended*Density_i*TurbVar_i[1]*TurbVar_i[1]*Volume;

   /*--- Cross diffusion ---*/

   val_residual[1] += (1.0 - F1_i)*CDkw_i*Volume;

   /*--- Implicit part ---*/

   val_Jacobian_i[0][0] = -beta_star*TurbVar_i[1]*Volume;
   val_Jacobian_i[0][1] = -beta_star*TurbVar_i[0]*Volume;
   val_Jacobian_i[1][0] = 0.0;
   val_Jacobian_i[1][1] = -2.0*beta_blended*TurbVar_i[1]*Volume;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

void CSourcePieceWise_TurbSST::GetMeanRateOfStrainMatrix(su2double **S_ij)
{
    /* --- Calculate the rate of strain tensor, using mean velocity gradients --- */

  if (nDim == 3){
    S_ij[0][0] = PrimVar_Grad_i[1][0];
    S_ij[1][1] = PrimVar_Grad_i[2][1];
    S_ij[2][2] = PrimVar_Grad_i[3][2];
    S_ij[0][1] = 0.5 * (PrimVar_Grad_i[1][1] + PrimVar_Grad_i[2][0]);
    S_ij[0][2] = 0.5 * (PrimVar_Grad_i[1][2] + PrimVar_Grad_i[3][0]);
    S_ij[1][2] = 0.5 * (PrimVar_Grad_i[2][2] + PrimVar_Grad_i[3][1]);
    S_ij[1][0] = S_ij[0][1];
    S_ij[2][1] = S_ij[1][2];
    S_ij[2][0] = S_ij[0][2];
  }
  else {
    S_ij[0][0] = PrimVar_Grad_i[1][0];
    S_ij[1][1] = PrimVar_Grad_i[2][1];
    S_ij[2][2] = 0.0;
    S_ij[0][1] = 0.5 * (PrimVar_Grad_i[1][1] + PrimVar_Grad_i[2][0]);
    S_ij[0][2] = 0.0;
    S_ij[1][2] = 0.0;
    S_ij[1][0] = S_ij[0][1];
    S_ij[2][1] = S_ij[1][2];
    S_ij[2][0] = S_ij[0][2];

  }
}

void CSourcePieceWise_TurbSST::SetReynoldsStressMatrix(su2double turb_ke){
  unsigned short iDim, jDim;
  su2double **S_ij = new su2double* [3];
  su2double divVel = 0;
  su2double TWO3 = 2.0/3.0;



  for (iDim = 0; iDim < 3; iDim++){
    S_ij[iDim] = new su2double [3];
  }

  GetMeanRateOfStrainMatrix(S_ij);

    /* --- Using rate of strain matrix, calculate Reynolds stress tensor --- */

  for (iDim = 0; iDim < 3; iDim++){
    divVel += S_ij[iDim][iDim];
  }

  for (iDim = 0; iDim < 3; iDim++){
    for (jDim = 0; jDim < 3; jDim++){
      MeanReynoldsStress[iDim][jDim] = TWO3 * turb_ke * delta3[iDim][jDim]
      - Eddy_Viscosity_i / Density_i * (2 * S_ij[iDim][jDim] - TWO3 * divVel * delta3[iDim][jDim]);
    }
  }

  for (iDim = 0; iDim < 3; iDim++)
    delete [] S_ij[iDim];
  delete [] S_ij;
}

void CSourcePieceWise_TurbSST::SetPerturbedRSM(su2double turb_ke, CConfig *config){

  unsigned short iDim,jDim;

  /* --- Calculate anisotropic part of Reynolds Stress tensor --- */

  for (iDim = 0; iDim< 3; iDim++){
    for (jDim = 0; jDim < 3; jDim++){
      A_ij[iDim][jDim] = .5 * MeanReynoldsStress[iDim][jDim] / turb_ke - delta3[iDim][jDim] / 3.0;
      Eig_Vec[iDim][jDim] = A_ij[iDim][jDim];
    }
  }

  /* --- Get ordered eigenvectors and eigenvalues of A_ij --- */

  EigenDecomposition(A_ij, Eig_Vec, Eig_Val, 3);

  /* compute convex combination coefficients */
  su2double c1c = Eig_Val[2] - Eig_Val[1];
  su2double c2c = 2.0 * (Eig_Val[1] - Eig_Val[0]);
  su2double c3c = 3.0 * Eig_Val[0] + 1.0;

  /* define barycentric traingle corner points */
  Corners[0][0] = 1.0;
  Corners[0][1] = 0.0;
  Corners[1][0] = 0.0;
  Corners[1][1] = 0.0;
  Corners[2][0] = 0.5;
  Corners[2][1] = 0.866025;

  /* define barycentric coordinates */
  Barycentric_Coord[0] = Corners[0][0] * c1c + Corners[1][0] * c2c + Corners[2][0] * c3c;
  Barycentric_Coord[1] = Corners[0][1] * c1c + Corners[1][1] * c2c + Corners[2][1] * c3c;

  if (Eig_Val_Comp == 1) {
    /* 1C turbulence */
    New_Coord[0] = Corners[0][0];
    New_Coord[1] = Corners[0][1];
  }
  else if (Eig_Val_Comp == 2) {
    /* 2C turbulence */
    New_Coord[0] = Corners[1][0];
    New_Coord[1] = Corners[1][1];
  }
  else if (Eig_Val_Comp == 3) {
    /* 3C turbulence */
    New_Coord[0] = Corners[2][0];
    New_Coord[1] = Corners[2][1];
  }
  else {
    /* 2C turbulence */
    New_Coord[0] = Corners[1][0];
    New_Coord[1] = Corners[1][1];
  }
  /* calculate perturbed barycentric coordinates */

  Barycentric_Coord[0] = Barycentric_Coord[0] + (uq_delta_b) * (New_Coord[0] - Barycentric_Coord[0]);
  Barycentric_Coord[1] = Barycentric_Coord[1] + (uq_delta_b) * (New_Coord[1] - Barycentric_Coord[1]);

  /* rebuild c1c,c2c,c3c based on new barycentric coordinates */
  c3c = Barycentric_Coord[1] / Corners[2][1];
  c1c = Barycentric_Coord[0] - Corners[2][0] * c3c;
  c2c = 1 - c1c - c3c;

  /* build new anisotropy eigenvalues */
  Eig_Val[0] = (c3c - 1) / 3.0;
  Eig_Val[1] = 0.5 *c2c + Eig_Val[0];
  Eig_Val[2] = c1c + Eig_Val[1];

  /* permute eigenvectors if required */
  if (uq_permute) {
    for (iDim=0; iDim<3; iDim++) {
      for (jDim=0; jDim<3; jDim++) {
        New_Eig_Vec[iDim][jDim] = Eig_Vec[2-iDim][jDim];
      }
    }
  }

  else {
    for (iDim=0; iDim<3; iDim++) {
      for (jDim=0; jDim<3; jDim++) {
        New_Eig_Vec[iDim][jDim] = Eig_Vec[iDim][jDim];
      }
    }
  }

  EigenRecomposition(newA_ij, New_Eig_Vec, Eig_Val, 3);

  /* compute perturbed Reynolds stress matrix; use under-relaxation factor (urlx)*/
  for (iDim = 0; iDim< 3; iDim++){
    for (jDim = 0; jDim < 3; jDim++){
      MeanPerturbedRSM[iDim][jDim] = 2.0 * turb_ke * (newA_ij[iDim][jDim] + 1.0/3.0 * delta3[iDim][jDim]);
      MeanPerturbedRSM[iDim][jDim] = MeanReynoldsStress[iDim][jDim] +
      uq_urlx*(MeanPerturbedRSM[iDim][jDim] - MeanReynoldsStress[iDim][jDim]);
    }
  }

}

void CSourcePieceWise_TurbSST::SetPerturbedStrainMag(su2double turb_ke){
  unsigned short iDim, jDim;
  PerturbedStrainMag = 0;
  su2double **StrainRate = new su2double* [nDim];
  for (iDim= 0; iDim< nDim; iDim++){
    StrainRate[iDim] = new su2double [nDim];
  }

  /* compute perturbed strain rate tensor */

  for (iDim = 0; iDim < nDim; iDim++){
    for (jDim =0; jDim < nDim; jDim++){
      StrainRate[iDim][jDim] = MeanPerturbedRSM[iDim][jDim]
      - TWO3 * turb_ke * delta[iDim][jDim];
      StrainRate[iDim][jDim] = - StrainRate[iDim][jDim] * Density_i / (2 * Eddy_Viscosity_i);
    }
  }

  /*--- Add diagonal part ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    PerturbedStrainMag += pow(StrainRate[iDim][iDim], 2.0);
  }

  /*--- Add off diagonals ---*/

  PerturbedStrainMag += 2.0*pow(StrainRate[1][0], 2.0);

  if (nDim == 3) {
    PerturbedStrainMag += 2.0*pow(StrainRate[0][2], 2.0);
    PerturbedStrainMag += 2.0*pow(StrainRate[1][2], 2.0);
  }

  PerturbedStrainMag = sqrt(2.0*PerturbedStrainMag);

  for (iDim= 0; iDim< nDim; iDim++){
    delete [] StrainRate[iDim];
  }

  delete [] StrainRate;
}
