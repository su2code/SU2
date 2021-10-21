/*!
 * \file turb_sources.cpp
 * \brief Implementation of numerics classes for integration of
 *        turbulence source-terms.
 * \author F. Palacios, T. Economon
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/turbulent/turb_sources.hpp"

CSourceBase_TurbSA::CSourceBase_TurbSA(unsigned short val_nDim,
                                       unsigned short val_nVar,
                                       const CConfig* config) :
  CNumerics(val_nDim, val_nVar, config),
  incompressible(config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE),
  rotating_frame(config->GetRotating_Frame())
{
  /*--- Spalart-Allmaras closure constants ---*/

  cv1_3 = pow(7.1, 3.0);
  k2    = pow(0.41, 2.0);
  cb1   = 0.1355;
  cw2   = 0.3;
  ct3   = 1.2;
  ct4   = 0.5;
  cw3_6 = pow(2.0, 6.0);
  sigma = 2./3.;
  cb2   = 0.622;
  cb2_sigma = cb2/sigma;
  cw1 = cb1/k2+(1.0+cb2)/sigma;
  cr1 = 0.5;

  /*--- Setup the Jacobian pointer, we need to return su2double** but
   *    we know the Jacobian is 1x1 so we use this "trick" to avoid
   *    having to dynamically allocate. ---*/

  Jacobian_i = &Jacobian_Buffer;

}

CSourcePieceWise_TurbSST::CSourcePieceWise_TurbSST(unsigned short val_nDim,
                                                   unsigned short val_nVar,
                                                   const su2double *constants,
                                                   su2double val_kine_Inf,
                                                   su2double val_omega_Inf,
                                                   const CConfig* config) :
                          CNumerics(val_nDim, val_nVar, config) {

  incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
  sustaining_terms = (config->GetKind_Turb_Model() == SST_SUST);
  axisymmetric = config->GetAxisymmetric();

  /*--- Closure constants ---*/
  sigma_k_1     = constants[0];
  sigma_k_2     = constants[1];
  sigma_w_1     = constants[2];
  sigma_w_2     = constants[3];
  beta_1        = constants[4];
  beta_2        = constants[5];
  beta_star     = constants[6];
  a1            = constants[7];
  alfa_1        = constants[8];
  alfa_2        = constants[9];

  /*--- Set the ambient values of k and omega to the free stream values. ---*/
  kAmb     = val_kine_Inf;
  omegaAmb = val_omega_Inf;

  /*--- "Allocate" the Jacobian using the static buffer. ---*/
  Jacobian_i[0] = Jacobian_Buffer;
  Jacobian_i[1] = Jacobian_Buffer+2;

}

CNumerics::ResidualType<> CSourcePieceWise_TurbSST::ComputeResidual(const CConfig* config) {

  AD::StartPreacc();
  AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(ScalarVar_i, nVar);
  AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);
  AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F2_i); AD::SetPreaccIn(CDkw_i);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);
  AD::SetPreaccIn(Vorticity_i, 3);

  unsigned short iDim;
  su2double alfa_blended, beta_blended;
  su2double diverg, pk, pw, zeta;
  su2double VorticityMag = sqrt(Vorticity_i[0]*Vorticity_i[0] +
                                Vorticity_i[1]*Vorticity_i[1] +
                                Vorticity_i[2]*Vorticity_i[2]);

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

  Residual[0] = 0.0;       Residual[1] = 0.0;
  Jacobian_i[0][0] = 0.0;  Jacobian_i[0][1] = 0.0;
  Jacobian_i[1][0] = 0.0;  Jacobian_i[1][1] = 0.0;

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
     ComputePerturbedRSM(nDim, Eig_Val_Comp, uq_permute, uq_delta_b, uq_urlx,
                         PrimVar_Grad_i+1, Density_i, Eddy_Viscosity_i,
                         ScalarVar_i[0], MeanPerturbedRSM);
     SetPerturbedStrainMag(ScalarVar_i[0]);
     pk = Eddy_Viscosity_i*PerturbedStrainMag*PerturbedStrainMag
          - 2.0/3.0*Density_i*ScalarVar_i[0]*diverg;
   }
   else {
     pk = Eddy_Viscosity_i*StrainMag_i*StrainMag_i - 2.0/3.0*Density_i*ScalarVar_i[0]*diverg;
   }

   pk = min(pk,20.0*beta_star*Density_i*ScalarVar_i[1]*ScalarVar_i[0]);
   pk = max(pk,0.0);

   zeta = max(ScalarVar_i[1], VorticityMag*F2_i/a1);

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

   Residual[0] += pk*Volume;
   Residual[1] += pw*Volume;

   /*--- Dissipation ---*/

   Residual[0] -= beta_star*Density_i*ScalarVar_i[1]*ScalarVar_i[0]*Volume;
   Residual[1] -= beta_blended*Density_i*ScalarVar_i[1]*ScalarVar_i[1]*Volume;

   /*--- Cross diffusion ---*/

   Residual[1] += (1.0 - F1_i)*CDkw_i*Volume;

   /*--- Contribution due to 2D axisymmetric formulation ---*/

   if (axisymmetric) ResidualAxisymmetric(alfa_blended,zeta);

   /*--- Implicit part ---*/

   Jacobian_i[0][0] = -beta_star*ScalarVar_i[1]*Volume;
   Jacobian_i[0][1] = -beta_star*ScalarVar_i[0]*Volume;
   Jacobian_i[1][0] = 0.0;
   Jacobian_i[1][1] = -2.0*beta_blended*ScalarVar_i[1]*Volume;
  }

  AD::SetPreaccOut(Residual, nVar);
  AD::EndPreacc();

  return ResidualType<>(Residual, Jacobian_i, nullptr);

}

void CSourcePieceWise_TurbSST::SetPerturbedStrainMag(su2double turb_ke){

  /*--- Compute norm of perturbed strain rate tensor. ---*/

  PerturbedStrainMag = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++){
    for (unsigned short jDim = 0; jDim < nDim; jDim++){
      su2double StrainRate_ij = MeanPerturbedRSM[iDim][jDim] - TWO3 * turb_ke * delta[iDim][jDim];
      StrainRate_ij = - StrainRate_ij * Density_i / (2 * Eddy_Viscosity_i);

      PerturbedStrainMag += pow(StrainRate_ij, 2.0);
    }
  }
  PerturbedStrainMag = sqrt(2.0*PerturbedStrainMag);

}
