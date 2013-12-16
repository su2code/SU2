/*!
 * \file numerics_adjoint_mean.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.9
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
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

#include "../include/numerics_structure.hpp"
#include <limits>

CUpwRoe_AdjTNE2::CUpwRoe_AdjTNE2(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 unsigned short val_nPrimVar,
                                 unsigned short val_nPrimVarGrad,
                                 CConfig *config) : CNumerics(val_nDim,
                                                              val_nVar,
                                                              config) {
  
  /*--- Read configuration parameters ---*/
	implicit   = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  /*--- Define useful constants ---*/
  nVar         = val_nVar;
  nPrimVar     = val_nPrimVar;
  nPrimVarGrad = val_nPrimVarGrad;
  nDim         = val_nDim;
  nSpecies     = config->GetnSpecies();
  
  UnitNormal = new double[nDim];
  MeanU      = new double[nVar];
  MeanV      = new double[nPrimVar];
  MeandPdU   = new double[nVar];
  DiffPsi    = new double[nVar];
  Lambda     = new double[nVar];
  Ai     = new double* [nVar];
  Aj     = new double* [nVar];
  P      = new double* [nVar];
  invP   = new double* [nVar];
  PLPinv = new double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Ai[iVar]     = new double [nVar];
    Aj[iVar]     = new double [nVar];
    P[iVar]      = new double [nVar];
    invP[iVar]   = new double [nVar];
    PLPinv[iVar] = new double [nVar];
  }  
}

CUpwRoe_AdjTNE2::~CUpwRoe_AdjTNE2(void) {
  
  delete [] UnitNormal;
  delete [] MeanU;
  delete [] MeanV;
  delete [] MeandPdU;
  delete [] DiffPsi;
	delete [] Lambda;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] Ai[iVar];
    delete [] Aj[iVar];
    delete [] P[iVar];
    delete [] invP[iVar];
    delete [] PLPinv[iVar];
  }
  delete [] Ai;
  delete [] Aj;
  delete [] P;
  delete [] invP;
  delete [] PLPinv;
  
}

void CUpwRoe_AdjTNE2::ComputeResidual (double *val_residual_i,
                                       double *val_residual_j,
                                       double **val_Jacobian_ii,
                                       double **val_Jacobian_ij,
                                       double **val_Jacobian_ji,
                                       double **val_Jacobian_jj,
                                       CConfig *config) {
  
  unsigned short iDim, iVar, jVar, kVar;
  double Area, ProjVel;
  double MeanSoundSpeed;
  
  /*--- Roe flux: Fij = (Fi + Fj)/2 - 1/2*P|Lam|P^-1 * (Uj - Ui) ---*/
  // Notes:
  // 1) Non-conservative method, so for Fij -> Fi = A_i*Psi_i & Fj = A_i*Psi_j
  //    and Fji -> Fi = A_j*Psi_i & Fj = A_j*Psi_j
  // 2) Linear upwinding, so eigenvalue & eigenvector decomposition can be
  //    calculated using interface ij state (mean variables)
  
  /*--- Initialize the residuals (and Jacobians) ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual_i[iVar] = 0.0;
    val_residual_j[iVar] = 0.0;
  }
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_ii[iVar][jVar] = 0.0;
        val_Jacobian_ij[iVar][jVar] = 0.0;
        val_Jacobian_ji[iVar][jVar] = 0.0;
        val_Jacobian_jj[iVar][jVar] = 0.0;
      }
  }
  
  /*--- Calculate geometric quantities ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Calculate inviscid projected Jacobians ---*/
  GetInviscidProjJac(U_i, V_i, dPdU_i, Normal, 1.0, Ai);
  GetInviscidProjJac(U_j, V_j, dPdU_j, Normal, 1.0, Aj);
  
  /*--- Inviscid portion of the flux, A^T*Psi ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual_i[iVar] +=  0.5*Ai[jVar][iVar]*(Psi_i[jVar]+Psi_j[jVar]);
      val_residual_j[iVar] += -0.5*Aj[jVar][iVar]*(Psi_i[jVar]+Psi_j[jVar]);
    }
  }
  
  /*--- Calculate mean variables ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    MeanU[iVar] = 0.5*(U_i[iVar]+U_j[iVar]);
  for (iVar = 0; iVar < nPrimVar; iVar++)
    MeanV[iVar] = 0.5*(V_i[iVar]+V_j[iVar]);
  var->CalcdPdU(MeanV, config, MeandPdU);
  MeanSoundSpeed = sqrt((1.0+MeandPdU[nSpecies+nDim]) *
                        MeanV[P_INDEX]/MeanV[RHO_INDEX]);
  
  for (iVar = 0; iVar < nVar; iVar++)
    DiffPsi[iVar] = Psi_j[iVar] - Psi_i[iVar];
  
  ProjVel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ProjVel += MeanV[VEL_INDEX+iDim]*UnitNormal[iDim];
  
  /*--- Calculate eigenvalues of the interface state, ij ---*/
  for (iVar = 0; iVar < nSpecies+nDim-1; iVar++)
    Lambda[iVar] = ProjVel;
  Lambda[nSpecies+nDim-1] = ProjVel + MeanSoundSpeed;
  Lambda[nSpecies+nDim]   = ProjVel - MeanSoundSpeed;
  Lambda[nSpecies+nDim+1] = ProjVel;
  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = fabs(Lambda[iVar]);
  
  /*--- Calculate left and right eigenvector matrices ---*/
  CreateBasis(UnitNormal);
  GetPMatrix(MeanU, MeanV, MeandPdU, UnitNormal, l, m, P);
  GetPMatrix_inv(MeanU, MeanV, MeandPdU, UnitNormal, l, m, invP);
  
  /*--- Calculate eigenvalue/eigenvector decomposition ---*/
  // |PLPinv| = P x |Lambda| x inverse P
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      PLPinv[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        PLPinv[iVar][jVar] += P[iVar][kVar]*Lambda[kVar]*invP[kVar][jVar];
    }
  }
  
  /*--- Calculate the 'viscous' portion of the flux ---*/
  // 1/2*(P|Lam|P^-1)^T * (Uj - Ui)
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual_i[iVar] -= 0.5*PLPinv[jVar][iVar]*(Psi_i[jVar]-Psi_j[jVar])*Area;
      val_residual_j[iVar] += 0.5*PLPinv[jVar][iVar]*(Psi_i[jVar]-Psi_j[jVar])*Area;
    }
  }
  
  /*--- Populate Jacobian matrices ---*/
  // Note: Ai/j calculated using 'Normal', but PLPinv computed using UnitNormal.
  //       Only multiply PLP by area to properly account for integration.
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_ii[iVar][jVar] = 0.5*Ai[jVar][iVar] - 0.5*PLPinv[jVar][iVar]*Area;
        val_Jacobian_ij[iVar][jVar] = 0.5*Ai[jVar][iVar] + 0.5*PLPinv[jVar][iVar]*Area;
        val_Jacobian_ji[iVar][jVar] = -0.5*Aj[jVar][iVar] + 0.5*PLPinv[jVar][iVar]*Area;
        val_Jacobian_jj[iVar][jVar] = -0.5*Aj[jVar][iVar] - 0.5*PLPinv[jVar][iVar]*Area;
      }
    }
  }
}

CUpwSW_AdjTNE2::CUpwSW_AdjTNE2(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 unsigned short val_nPrimVar,
                                 unsigned short val_nPrimVarGrad,
                                 CConfig *config) : CNumerics(val_nDim,
                                                              val_nVar,
                                                              config) {
  
  /*--- Read configuration parameters ---*/
	implicit   = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  /*--- Define useful constants ---*/
  nVar         = val_nVar;
  nPrimVar     = val_nPrimVar;
  nPrimVarGrad = val_nPrimVarGrad;
  nDim         = val_nDim;
  nSpecies     = config->GetnSpecies();
  
  UnitNormal = new double[nDim];
  DiffPsi    = new double[nVar];
  Lambda_i   = new double[nVar];
  Lambda_j   = new double[nVar];
  Ai     = new double* [nVar];
  Aj     = new double* [nVar];
  P      = new double* [nVar];
  invP   = new double* [nVar];
  PLPinv = new double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Ai[iVar]     = new double [nVar];
    Aj[iVar]     = new double [nVar];
    P[iVar]      = new double [nVar];
    invP[iVar]   = new double [nVar];
    PLPinv[iVar] = new double [nVar];
  }
  
}

CUpwSW_AdjTNE2::~CUpwSW_AdjTNE2(void) {
  
  delete [] UnitNormal;
  delete [] DiffPsi;
	delete [] Lambda_i;
  delete [] Lambda_j;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] Ai[iVar];
    delete [] Aj[iVar];
    delete [] P[iVar];
    delete [] invP[iVar];
    delete [] PLPinv[iVar];
  }
  delete [] P;
  delete [] invP;
  delete [] PLPinv;
}

void CUpwSW_AdjTNE2::ComputeResidual (double *val_residual_i,
                                      double *val_residual_j,
                                      double **val_Jacobian_ii,
                                      double **val_Jacobian_ij,
                                      double **val_Jacobian_ji,
                                      double **val_Jacobian_jj,
                                      CConfig *config) {
  
  unsigned short iDim, iVar, jVar, kVar;
  double Area, ProjVel_i, ProjVel_j;
  
  /*--- Roe flux: Fij = (Fi + Fj)/2 - 1/2*P|Lam|P^-1 * (Uj - Ui) ---*/
  // Notes:
  // 1) Non-conservative method, so for Fij -> Fi = A_i*Psi_i & Fj = A_i*Psi_j
  //    and Fji -> Fi = A_j*Psi_i & Fj = A_j*Psi_j
  // 2) Linear upwinding, so eigenvalue & eigenvector decomposition can be
  //    calculated using interface ij state (mean variables)
  
  /*--- Initialize the residuals (and Jacobians) ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual_i[iVar] = 0.0;
    val_residual_j[iVar] = 0.0;
  }
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_ii[iVar][jVar] = 0.0;
        val_Jacobian_ij[iVar][jVar] = 0.0;
        val_Jacobian_ji[iVar][jVar] = 0.0;
        val_Jacobian_jj[iVar][jVar] = 0.0;
      }
  }
  
  /*--- Calculate geometric quantities ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  CreateBasis(UnitNormal);
  
  /*--- Calculate inviscid projected Jacobians ---*/
  GetInviscidProjJac(U_i, V_i, dPdU_i, Normal, 1.0, Ai);
  GetInviscidProjJac(U_j, V_j, dPdU_j, Normal, 1.0, Aj);
  
  /*--- Inviscid portion of the flux, A^T*Psi ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual_i[iVar] +=  0.5*Ai[jVar][iVar]*(Psi_i[jVar]+Psi_j[jVar]);
      val_residual_j[iVar] += -0.5*Aj[jVar][iVar]*(Psi_i[jVar]+Psi_j[jVar]);
    }
  }
  
  for (iVar = 0; iVar < nVar; iVar++)
    DiffPsi[iVar] = Psi_j[iVar] - Psi_i[iVar];
  
  ProjVel_i = 0.0;
  ProjVel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVel_i += V_i[VEL_INDEX+iDim]*UnitNormal[iDim];
    ProjVel_j += V_j[VEL_INDEX+iDim]*UnitNormal[iDim];
  }
  
  /*--- Calculate eigenvalues of the i & j states ---*/
  for (iVar = 0; iVar < nSpecies+nDim-1; iVar++) {
    Lambda_i[iVar] = ProjVel_i;
    Lambda_j[iVar] = ProjVel_j;
  }
  Lambda_i[nSpecies+nDim-1] = ProjVel_i + V_i[A_INDEX];
  Lambda_j[nSpecies+nDim-1] = ProjVel_j + V_j[A_INDEX];
  Lambda_i[nSpecies+nDim]   = ProjVel_i - V_i[A_INDEX];
  Lambda_j[nSpecies+nDim]   = ProjVel_j - V_j[A_INDEX];
  Lambda_i[nSpecies+nDim+1] = ProjVel_i;
  Lambda_j[nSpecies+nDim+1] = ProjVel_j;
  for (iVar = 0; iVar < nVar; iVar++) {
    Lambda_i[iVar] = fabs(Lambda_i[iVar]);
    Lambda_j[iVar] = fabs(Lambda_j[iVar]);
  }
  
  /*--- Calculate left and right eigenvector matrices ---*/
  GetPMatrix(U_i, V_i, dPdU_i, UnitNormal, l, m, P);
  GetPMatrix_inv(U_i, V_i, dPdU_i, UnitNormal, l, m, invP);
  
  /*--- Calculate eigenvalue/eigenvector decomposition for i ---*/
  // |PLPinv| = P x |Lambda| x inverse P
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      PLPinv[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        PLPinv[iVar][jVar] += P[iVar][kVar]*Lambda_i[kVar]*invP[kVar][jVar];
    }
  }
  
  /*--- Calculate the 'viscous' portion of the flux at i ---*/
  // 1/2*(P|Lam|P^-1)^T * (Uj - Ui)
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      val_residual_i[iVar] -= 0.5*PLPinv[jVar][iVar]*(Psi_i[jVar]-Psi_j[jVar])*Area;
  
  /*--- Populate Jacobian matrices ---*/
  // Note: Ai/j calculated using 'Normal', but PLPinv computed using UnitNormal.
  //       Only multiply PLP by area to properly account for integration.
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_ii[iVar][jVar] = 0.5*Ai[jVar][iVar] - 0.5*PLPinv[jVar][iVar]*Area;
        val_Jacobian_ij[iVar][jVar] = 0.5*Ai[jVar][iVar] + 0.5*PLPinv[jVar][iVar]*Area;
      }
    }
  }
  
  /*--- Calculate left and right eigenvector matrices for j ---*/
  GetPMatrix(U_j, V_j, dPdU_j, UnitNormal, l, m, P);
  GetPMatrix_inv(U_j, V_j, dPdU_j, UnitNormal, l, m, invP);
  
  /*--- Calculate eigenvalue/eigenvector decomposition for i ---*/
  // |PLPinv| = P x |Lambda| x inverse P
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      PLPinv[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        PLPinv[iVar][jVar] += P[iVar][kVar]*Lambda_j[kVar]*invP[kVar][jVar];
    }
  }
  
  /*--- Calculate the 'viscous' portion of the flux at j ---*/
  // 1/2*(P|Lam|P^-1)^T * (Uj - Ui)
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      val_residual_j[iVar] += 0.5*PLPinv[jVar][iVar]*(Psi_i[jVar]-Psi_j[jVar])*Area;
  
  /*--- Populate Jacobian matrices ---*/
  // Note: Ai/j calculated using 'Normal', but PLPinv computed using UnitNormal.
  //       Only multiply PLP by area to properly account for integration.
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_ji[iVar][jVar] = -0.5*Aj[jVar][iVar] + 0.5*PLPinv[jVar][iVar]*Area;
        val_Jacobian_jj[iVar][jVar] = -0.5*Aj[jVar][iVar] - 0.5*PLPinv[jVar][iVar]*Area;
      }
    }
  }
}


CCentJST_AdjTNE2::CCentJST_AdjTNE2(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	grid_movement = config->GetGrid_Movement();
	rotating_frame = config->GetRotating_Frame();
  
	Diff_Psi = new double [nVar]; Diff_Lapl = new double [nVar];
	Und_Lapl_i = new double [nVar]; Und_Lapl_j = new double [nVar];
	Velocity_i = new double [nDim]; Velocity_j = new double [nDim];
	MeanPhi = new double [nDim];
  
	Param_p = 0.3;
	Param_Kappa_2 = config->GetKappa_2nd_AdjTNE2();
	Param_Kappa_4 = config->GetKappa_4th_AdjTNE2();
	implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
}

CCentJST_AdjTNE2::~CCentJST_AdjTNE2(void) {
  
	delete [] Diff_Psi; delete [] Diff_Lapl;
	delete [] Und_Lapl_i; delete [] Und_Lapl_j;
	delete [] Velocity_i; delete [] Velocity_j;
	delete [] MeanPhi;
}

void CCentJST_AdjTNE2::ComputeResidual (double *val_resconv_i, double *val_resvisc_i,
                                        double *val_resconv_j, double *val_resvisc_j,
                                        double **val_Jacobian_ii, double **val_Jacobian_ij,
                                        double **val_Jacobian_ji, double **val_Jacobian_jj,
                                        CConfig *config) {
  
	/*--- Mean Values ---*/
	MeanPsiRho =  0.5*(Psi_i[0]+Psi_j[0]);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanPhi[iDim] =  0.5*(Psi_i[iDim+1]+Psi_j[iDim+1]);
	MeanPsiE =  0.5*(Psi_i[nVar-1]+Psi_j[nVar-1]);
  
	/*--- Point i convective residual evaluation ---*/
	ProjVelocity_i = 0; ProjPhi = 0; ProjPhi_Vel = 0; sq_vel = 0; Area = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1] / U_i[0];
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		ProjPhi += MeanPhi[iDim]*Normal[iDim];
		ProjPhi_Vel += MeanPhi[iDim]*Velocity_i[iDim];
		sq_vel += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
	Area = sqrt(Area);
	phis1 = ProjPhi + ProjVelocity_i*MeanPsiE;
	phis2 = MeanPsiRho + ProjPhi_Vel + Enthalpy_i*MeanPsiE;
  
	val_resconv_i[0] = ProjVelocity_i*MeanPsiRho - phis2*ProjVelocity_i + Gamma_Minus_One*phis1*sq_vel;
	for (iDim = 0; iDim < nDim; iDim++)
		val_resconv_i[iDim+1] = ProjVelocity_i*MeanPhi[iDim] + phis2*Normal[iDim] - Gamma_Minus_One*phis1*Velocity_i[iDim];
	val_resconv_i[nVar-1] = ProjVelocity_i*MeanPsiE + Gamma_Minus_One*phis1;
  
	/*--- Flux contributions due to grid movement at point i (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		val_resconv_i[0] -= ProjGridVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_i[iDim+1] -= ProjGridVel*MeanPhi[iDim];
		val_resconv_i[nVar-1] -= ProjGridVel*MeanPsiE;
	}
  
	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		val_Jacobian_ii[0][0] = 0.0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_ii[0][jDim+1] = -0.5*ProjVelocity_i*Velocity_i[jDim] + Gamma_Minus_One*sq_vel*0.5*Normal[jDim];
		val_Jacobian_ii[0][nVar-1] = 0.5*ProjVelocity_i*(Gamma_Minus_One*sq_vel - Enthalpy_i);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_Jacobian_ii[iDim+1][0] = 0.5*Normal[iDim];
			for (jDim = 0; jDim < nDim; jDim++)
				val_Jacobian_ii[iDim+1][jDim+1] = 0.5*Normal[iDim]*Velocity_i[jDim] - 0.5*Gamma_Minus_One*Velocity_i[iDim]*Normal[jDim];
			val_Jacobian_ii[iDim+1][iDim+1] += 0.5*ProjVelocity_i;
			val_Jacobian_ii[iDim+1][nVar-1] = 0.5*Enthalpy_i*Normal[iDim] - 0.5*Gamma_Minus_One*Velocity_i[iDim]*ProjVelocity_i;
		}
		val_Jacobian_ii[nVar-1][0] = 0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_ii[nVar-1][jDim+1] = 0.5*Gamma_Minus_One*Normal[jDim];
		val_Jacobian_ii[nVar-1][nVar-1] = 0.5*Gamma*ProjVelocity_i;
    
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ij[iVar][jVar] = val_Jacobian_ii[iVar][jVar];
    
		/*--- Jacobian contributions due to grid movement at point i (TDE) ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjGridVel;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjGridVel;
			}
		}
	}
  
  
	/*--- Point j convective residual evaluation ---*/
	ProjVelocity_j = 0; ProjPhi_Vel = 0; sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_j[iDim] = U_j[iDim+1] / U_j[0];
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		ProjPhi_Vel += MeanPhi[iDim]*Velocity_j[iDim];
		sq_vel += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}
  
	phis1 = ProjPhi + ProjVelocity_j*MeanPsiE;
	phis2 = MeanPsiRho + ProjPhi_Vel + Enthalpy_j*MeanPsiE;
  
	val_resconv_j[0] = -(ProjVelocity_j*MeanPsiRho - phis2*ProjVelocity_j + Gamma_Minus_One*phis1*sq_vel);
	for (iDim = 0; iDim < nDim; iDim++)
		val_resconv_j[iDim+1] = -(ProjVelocity_j*MeanPhi[iDim] + phis2*Normal[iDim] - Gamma_Minus_One*phis1*Velocity_j[iDim]);
	val_resconv_j[nVar-1] = -(ProjVelocity_j*MeanPsiE + Gamma_Minus_One*phis1);
  
	/*--- Flux contributions due to grid movement at point j (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		val_resconv_j[0] += ProjGridVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_j[iDim+1] += ProjGridVel*MeanPhi[iDim];
		val_resconv_j[nVar-1] += ProjGridVel*MeanPsiE;
	}
  
	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		val_Jacobian_jj[0][0] = 0.0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_jj[0][jDim+1] = 0.5*ProjVelocity_j*Velocity_j[jDim] - Gamma_Minus_One*sq_vel*0.5*Normal[jDim];
		val_Jacobian_jj[0][nVar-1] = -0.5*ProjVelocity_j*(Gamma_Minus_One*sq_vel - Enthalpy_j);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_Jacobian_jj[iDim+1][0] = -0.5*Normal[iDim];
			for (jDim = 0; jDim < nDim; jDim++)
				val_Jacobian_jj[iDim+1][jDim+1] = -0.5*Normal[iDim]*Velocity_j[jDim] + 0.5*Gamma_Minus_One*Velocity_j[iDim]*Normal[jDim];
			val_Jacobian_jj[iDim+1][iDim+1] -= 0.5*ProjVelocity_j;
			val_Jacobian_jj[iDim+1][nVar-1] = -0.5*Enthalpy_j*Normal[iDim] + 0.5*Gamma_Minus_One*Velocity_j[iDim]*ProjVelocity_j;
		}
		val_Jacobian_jj[nVar-1][0] = 0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_jj[nVar-1][jDim+1] = -0.5*Gamma_Minus_One*Normal[jDim];
		val_Jacobian_jj[nVar-1][nVar-1] = -0.5*Gamma*ProjVelocity_j;
    
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ji[iVar][jVar] = val_Jacobian_jj[iVar][jVar];
    
		/*--- Jacobian contributions due to grid movement at point j (TDE) ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjGridVel;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjGridVel;
			}
		}
	}
  
	/*--- Computes differences btw. variables and Laplacians ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
		Diff_Psi[iVar]  = Psi_i[iVar]-Psi_j[iVar];
	}
  
	/*--- Adjustment to projected velocity due to mesh motion (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel_i = 0.0; double ProjGridVel_j = 0.0; double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			ProjGridVel_i += GridVel_i[iDim]*Normal[iDim];
			ProjGridVel_j += GridVel_j[iDim]*Normal[iDim];
		}
		ProjVelocity_i -= ProjGridVel;
		ProjVelocity_j += ProjGridVel;
	}
  
	/*--- Compute the spectral radius and stretching factor ---*/
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);
  
	double sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	sc4 = sc2*sc2/4.0;
	Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
	Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;
  
	/*--- Compute viscous residual 1st- & 3rd-order dissipation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Residual = (Epsilon_2*Diff_Psi[iVar]-Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
		val_resvisc_i[iVar] = -Residual;
		val_resvisc_j[iVar] =  Residual;
		if (implicit) {
			val_Jacobian_ii[iVar][iVar] -= Epsilon_2 + double(Neighbor_i+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_ij[iVar][iVar] += Epsilon_2 + double(Neighbor_j+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_ji[iVar][iVar] += Epsilon_2 + double(Neighbor_i+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_jj[iVar][iVar] -= Epsilon_2 + double(Neighbor_j+1)*Epsilon_4*StretchingFactor*MeanLambda;
		}
	}
}


CCentLax_AdjTNE2::CCentLax_AdjTNE2(unsigned short val_nDim,
                                   unsigned short val_nVar,
                                   unsigned short val_nPrimVar,
                                   unsigned short val_nPrimVarGrad,
                                   CConfig *config) : CNumerics(val_nDim,
                                                                val_nVar,
                                                                config) {

	implicit   = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  nDim         = val_nDim;
  nSpecies     = config->GetnSpecies();
  nVar         = val_nVar;
  nPrimVar     = val_nPrimVar;
  nPrimVarGrad = val_nPrimVarGrad;
  
	DiffPsi   = new double [nVar];
  MeanPsi    = new double [nVar];
  
  Proj_Jac_Tensor_i = new double*[nVar];
  Proj_Jac_Tensor_j = new double*[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Proj_Jac_Tensor_i[iVar] = new double [nVar];
    Proj_Jac_Tensor_j[iVar] = new double [nVar];
  }
  
	Param_p = 0.3;
	Param_Kappa_0 = config->GetKappa_1st_AdjTNE2();
  
}

CCentLax_AdjTNE2::~CCentLax_AdjTNE2(void) {
  
	delete [] DiffPsi;
  delete [] MeanPsi;

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] Proj_Jac_Tensor_i[iVar];
    delete [] Proj_Jac_Tensor_j[iVar];
  }
  delete [] Proj_Jac_Tensor_i;
  delete [] Proj_Jac_Tensor_j;
}

void CCentLax_AdjTNE2::ComputeResidual (double *val_resconv_i,
                                        double *val_resvisc_i,
                                        double *val_resconv_j,
                                        double *val_resvisc_j,
                                        double **val_Jacobian_ii,
                                        double **val_Jacobian_ij,
                                        double **val_Jacobian_ji,
                                        double **val_Jacobian_jj,
                                        CConfig *config) {

  unsigned short iDim, iVar, jVar;
  double ProjVel_i, ProjVel_j;
  double Phi_i, Phi_j;
  double Local_Lambda_i, Local_Lambda_j, MeanLambda;
  double Residual;
  double StretchingFactor, sc2, Epsilon_0;
  
  /*--- Initialize the residuals ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_resconv_i[iVar] = 0.0;
    val_resconv_j[iVar] = 0.0;
    val_resvisc_i[iVar] = 0.0;
    val_resvisc_j[iVar] = 0.0;
  }
  
  /*--- Compute geometric parameters ---*/
	Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
	Area = sqrt(Area);
	for (iDim = 0; iDim < nDim; iDim++) {
		UnitNormal[iDim] = Normal[iDim]/Area;
    if (fabs(UnitNormal[iDim]) < EPS) UnitNormal[iDim] = EPS;
  }
  
  /*--- Calculate the mean & differences of the adjoint variables ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    MeanPsi[iVar] = 0.5 * (Psi_i[iVar]+Psi_j[iVar]);
    DiffPsi[iVar] = Psi_i[iVar]-Psi_j[iVar];
  }
  
  /*--- Calculate inviscid projected flux Jacobians ---*/
  GetInviscidProjJac(U_i, V_i, dPdU_i, Normal, 1.0, Proj_Jac_Tensor_i);
  GetInviscidProjJac(U_j, V_j, dPdU_j, Normal, 1.0, Proj_Jac_Tensor_j);
  
  /*--- Compute inviscid residual at point i, A^T*Psi ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_resconv_i[iVar] += Proj_Jac_Tensor_i[jVar][iVar]*MeanPsi[jVar];
      val_resconv_j[iVar] -= Proj_Jac_Tensor_j[jVar][iVar]*MeanPsi[jVar];
    }
  }
  
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_ii[iVar][jVar] = Proj_Jac_Tensor_i[jVar][iVar];
        val_Jacobian_ij[iVar][jVar] = Proj_Jac_Tensor_i[jVar][iVar];
        val_Jacobian_jj[iVar][jVar] = -Proj_Jac_Tensor_j[jVar][iVar];
        val_Jacobian_ji[iVar][jVar] = -Proj_Jac_Tensor_j[jVar][iVar];
      }
    }
  }
  
  /*--- Compute spectral radius ---*/
  ProjVel_i = 0.0;
  ProjVel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVel_i += V_i[VEL_INDEX+iDim]*Normal[iDim];
    ProjVel_j -= V_j[VEL_INDEX+iDim]*Normal[iDim];
  }
	Local_Lambda_i = (fabs(ProjVel_i)+V_i[A_INDEX]*Area);
	Local_Lambda_j = (fabs(ProjVel_j)+V_j[A_INDEX]*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
	/*--- Compute streching factor ---*/
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);
  
	sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	Epsilon_0 = Param_Kappa_0*sc2*double(nDim)/3.0;
  
	/*--- Artifical dissipation evaluation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Residual = Epsilon_0*StretchingFactor*MeanLambda*DiffPsi[iVar];
		val_resvisc_i[iVar] = -Residual;
		val_resvisc_j[iVar] =  Residual;
	}
  
	/*--- Contribution to implicit part ---*/
	if (implicit) {
		for (iVar = 0; iVar < nVar; iVar++) {
			val_Jacobian_ii[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_ij[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_ji[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_jj[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
		}
	}
}

CSource_AdjTNE2::CSource_AdjTNE2(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 unsigned short val_nPrimVar,
                                 unsigned short val_nPrimVarGrad,
                                 CConfig *config) : CNumerics(val_nDim,
                                                              val_nVar,
                                                              config) {
  
  unsigned short iDim, iVar;
  
  /*--- Assign booleans from CConfig ---*/
  implicit   = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  /*--- Set array sizes ---*/
  nDim         = val_nDim;
  nSpecies     = config->GetnSpecies();
  nVar         = val_nVar;
  nPrimVar     = val_nPrimVar;
  nPrimVarGrad = val_nPrimVarGrad;
  
  rhos = new double[nSpecies];
  vel = new double[nDim];
  
  GInvRho  = new double[nDim];
  GVeloRho = new double*[nDim];
  GPhiGInvRho = new double[nDim];
  GPsiEZetaTau = new double[nDim];
  tau      = new double*[nDim];
  eta      = new double*[nDim];
  pi       = new double*[nDim];
  zeta     = new double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    GVeloRho[iDim] = new double[nDim];
    tau[iDim]      = new double[nDim];
    eta[iDim]      = new double[nDim];
    pi[iDim]       = new double[nDim];
    zeta[iDim]     = new double[nDim];
  }
  Av2 = new double *[nVar];
  Av3 = new double *[nVar];
  Av4 = new double *[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Av2[iVar] = new double[nVar];
    Av3[iVar] = new double[nVar];
    Av4[iVar] = new double[nVar];
  }
  
  
//	unsigned short iDim;
//  
//	Velocity = new double [nVar];
//	GradDensity = new double [nDim];
//	GradInvDensity = new double [nDim];
//	dPoDensity2 = new double [nDim];
//	alpha = new double [nDim];
//	beta = new double [nDim];
//	Sigma_5_vec = new double [nDim];
//  
//	GradVel_o_Rho = new double* [nDim];
//	sigma = new double* [nDim];
//	Sigma_phi = new double* [nDim];
//	Sigma_5_Tensor = new double* [nDim];
//	Sigma = new double* [nDim];
//	for (iDim = 0; iDim < nDim; iDim++) {
//		GradVel_o_Rho[iDim] = new double [nDim];
//		sigma[iDim] = new double [nDim];
//		Sigma_phi[iDim] = new double [nDim];
//		Sigma_5_Tensor[iDim] = new double [nDim];
//		Sigma[iDim] = new double [nDim];
//	}
}

CSource_AdjTNE2::~CSource_AdjTNE2(void) {

  unsigned short iDim, iVar;
  
  delete [] rhos;
  delete [] vel;
  delete [] GInvRho;
  delete [] GPhiGInvRho;
  delete [] GPsiEZetaTau;
  for (iDim = 0; iDim < nDim; iDim++) {
    delete [] GVeloRho[iDim];
    delete [] tau[iDim];
    delete [] eta[iDim];
    delete [] pi[iDim];
    delete [] zeta[iDim];
  }
  delete [] GVeloRho;
  delete [] tau;
  delete [] eta;
  delete [] pi;
  delete [] zeta;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Av2[iVar];
    delete [] Av3[iVar];
    delete [] Av4[iVar];
  }
  delete [] Av2;
  delete [] Av3;
  delete [] Av4;
  
  
//	for (iDim = 0; iDim < nDim; iDim++) {
//		delete [] GradVel_o_Rho[iDim];
//		delete [] sigma[iDim];
//		delete [] Sigma_phi[iDim];
//		delete [] Sigma_5_Tensor[iDim];
//		delete [] Sigma[iDim];
//	}
//	delete [] GradVel_o_Rho;
//	delete [] sigma;
//	delete [] Sigma_phi;
//	delete [] Sigma_5_Tensor;
//	delete [] Sigma;
//  
//	delete [] Velocity;
//	delete [] GradDensity;
//	delete [] GradInvDensity;
//	delete [] dPoDensity2;
//	delete [] alpha;
//	delete [] beta;
//	delete [] Sigma_5_vec;
}

void CSource_AdjTNE2::ComputeSourceViscous (double *val_residual, CConfig *config) {
  
	unsigned short iDim, jDim, iSpecies, iVar, jVar;
  double rho, sqvel, rhoCvtr, Cvtrs, *Ms, *xi, Ru;
  double div_vel, div_velorho, velGInvRho, GPhiEta, GPsiEvelPi;
  double mu2, mu3, mu4;
  double **GPsi;
  
  /*--- Initialize arrays ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.0;
    for (jVar = 0; jVar < nVar; jVar++) {
      Av2[iVar][jVar] = 0.0;
      Av3[iVar][jVar] = 0.0;
      Av4[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Rename for convenience ---*/
  Ms      = config->GetMolar_Mass();
  xi      = config->GetRotationModes();
  Ru      = UNIVERSAL_GAS_CONSTANT;
  rhoCvtr = V_i[RHOCVTR_INDEX];
  GPsi    = PsiVar_Grad_i;
  
  /*--- Assign useful values ---*/
  // Mixture & species density
  rho = V_i[RHO_INDEX];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhos[iSpecies] = V_i[RHOS_INDEX+iSpecies];
  // Velocity
  sqvel = 0.0;
  div_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
		vel[iDim] = V_i[VEL_INDEX+iDim];
		sqvel += vel[iDim]*vel[iDim];
    div_vel += PrimVar_Grad_i[VEL_INDEX+iDim][iDim];
	}
  // Transport coefficients
  mu2 = Laminar_Viscosity_i;
  mu3 = Thermal_Conductivity_i;
  mu4 = Thermal_Conductivity_ve_i;

  /*--- Gradients of flow variables ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    GInvRho[iDim] = -1.0/(rho*rho)*PrimVar_Grad_i[RHO_INDEX][iDim];
    for (jDim = 0; jDim < nDim; jDim++)
      GVeloRho[iDim][jDim] = -vel[iDim]/(rho*rho)*(PrimVar_Grad_i[RHO_INDEX][jDim])
                           +  1.0/rho*(PrimVar_Grad_i[VEL_INDEX+iDim][jDim]);
  }
  div_velorho = 0.0;
  velGInvRho = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    div_velorho += GVeloRho[iDim][iDim];
    velGInvRho += vel[iDim]*GInvRho[iDim];
  }
  
  /*--- Supporting matrices ---*/
  // Tau: Stress tensor matrix
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      tau[iDim][jDim] = PrimVar_Grad_i[VEL_INDEX+iDim][jDim] +
                        PrimVar_Grad_i[VEL_INDEX+jDim][iDim];
    }
    tau[iDim][iDim] -= 2.0/3.0*div_vel;
  }
  // Eta
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      eta[iDim][jDim] = GVeloRho[iDim][jDim] + GVeloRho[jDim][iDim];
    }
    eta[iDim][iDim] -= 2.0/3.0*div_velorho;
  }
  // Pi
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      pi[iDim][jDim] = eta[iDim][jDim] - 1.0/rho*tau[iDim][jDim];
    }
  }
  // Zeta
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      zeta[iDim][jDim] = vel[jDim]*GInvRho[iDim] - 2.0/3.0*vel[iDim]*GInvRho[jDim];
    }
  }
  
  /*--- Calculate supporting quantities ---*/
  // GradPhi dot GradInvRho - e.g. (diPhix*di(1/rho), diPhiy*di(1/rho), ... )
  for (iDim = 0; iDim < nDim; iDim++) {
    GPhiGInvRho[iDim] = 0.0;
    for (jDim = 0; jDim < nDim; jDim++) {
      GPhiGInvRho[iDim] += GPsi[nSpecies+iDim][jDim]*GInvRho[jDim];
    }
  }
  // (di(PsiE)*(Zetai1 + 1/rho*Taui1), di(PsiE)*(Zetai2 + 1/rho*Taui2), ... )
  for (iDim = 0; iDim < nDim; iDim++) {
    GPsiEZetaTau[iDim] = 0.0;
    for (jDim = 0; jDim < nDim; jDim++) {
      GPsiEZetaTau[iDim] += GPsi[nSpecies+nDim][jDim] * (zeta[jDim][iDim] +
                                                         1/rho*tau[jDim][iDim]);
    }
  }
  // GPhi:Eta
  GPhiEta = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      GPhiEta += GPsi[nSpecies+iDim][jDim]*eta[jDim][iDim];
  // GPsiE_i*(u_j*pi_ij)
  GPsiEvelPi = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      GPsiEvelPi += GPsi[nSpecies+nDim][iDim]*vel[jDim]*pi[iDim][jDim];
  
  /*--- Contribution to viscous residual from Av2 ---*/
  // Mass conservation
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    val_residual[iSpecies] = (-GPhiEta + GPsiEvelPi)*mu2*Volume;
  }
  // X-momentum
  val_residual[nSpecies]   = ((GPhiGInvRho[0] + 1.0/3.0*GPsi[nSpecies][0]*GInvRho[0]) +
                              (GPsi[nSpecies+1][0]*GInvRho[1] - 2.0/3.0*GPsi[nSpecies+1][1]*GInvRho[0]) +
                              (GPsi[nSpecies+2][0]*GInvRho[2] - 2.0/3.0*GPsi[nSpecies+2][2]*GInvRho[0]) +
                              (GPsi[nSpecies+3][0]*velGInvRho + GPsiEZetaTau[0]))*mu2*Volume;
  // Y-momentum
  val_residual[nSpecies+1] = ((GPsi[nSpecies][1]*GInvRho[0]   - 2.0/3.0*GPsi[nSpecies][0]*GInvRho[1]) +
                              (GPhiGInvRho[1] + 1.0/3.0*GPsi[nSpecies+1][1]*GInvRho[1]) +
                              (GPsi[nSpecies+2][1]*GInvRho[2] - 2.0/3.0*GPsi[nSpecies+2][2]*GInvRho[1]) +
                              (GPsi[nSpecies+3][1]*velGInvRho + GPsiEZetaTau[1]))*mu2*Volume;
  // Z-momentum
  val_residual[nSpecies+2] = ((GPsi[nSpecies][2]*GInvRho[0]   - 2.0/3.0*GPsi[nSpecies][0]*GInvRho[2]) +
                              (GPsi[nSpecies+1][2]*GInvRho[1] - 2.0/3.0*GPsi[nSpecies+1][1]*GInvRho[2]) +
                              (GPhiGInvRho[2] + 1.0/3.0*GPsi[nSpecies+2][2]*GInvRho[2]) +
                              (GPsi[nSpecies+3][2]*velGInvRho + GPsiEZetaTau[2]))*mu2*Volume;
  
  
  /*--- Contribution to viscous residual from Av3 ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    // Mass conservation
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Cvtrs = (3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
      
      val_residual[iSpecies] += GPsi[nSpecies+3][iDim] * (1/rhoCvtr * (vel[0]*PrimVar_Grad_i[VEL_INDEX][iDim]   +
                                                                       vel[1]*PrimVar_Grad_i[VEL_INDEX+1][iDim] +
                                                                       vel[2]*PrimVar_Grad_i[VEL_INDEX+2][iDim] -
                                                                       Cvtrs*PrimVar_Grad_i[T_INDEX][iDim])
                                                          -dTdU_i[iSpecies]/rhoCvtr *
                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim])* mu3 * Volume;
    }
    // X-momentum
    val_residual[nSpecies]   += GPsi[nSpecies+3][iDim] * (-1/rhoCvtr*PrimVar_Grad_i[VEL_INDEX][iDim]
                                                          -dTdU_i[nSpecies]/rhoCvtr *
                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
    // Y-momentum
    val_residual[nSpecies+1] += GPsi[nSpecies+3][iDim] * (-1/rhoCvtr*PrimVar_Grad_i[VEL_INDEX+1][iDim]
                                                          -dTdU_i[nSpecies+1]/rhoCvtr *
                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
    // Z-momentum
    val_residual[nSpecies+2] += GPsi[nSpecies+3][iDim] * (-1/rhoCvtr*PrimVar_Grad_i[VEL_INDEX+2][iDim]
                                                          -dTdU_i[nSpecies+2]/rhoCvtr *
                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
    // Energy
    val_residual[nSpecies+3] += GPsi[nSpecies+3][iDim] * (-dTdU_i[nSpecies+3]/rhoCvtr *
                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
    // Vib.-el. energy
    val_residual[nSpecies+4] += GPsi[nSpecies+3][iDim] * (-dTdU_i[nSpecies+4]/rhoCvtr *
                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
  }
  
  
  
  
//  /*--- Momentum viscous Jacobian (Av2) ---*/
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (iDim = 0; iDim < nDim; iDim++) {
//      Av2[nSpecies][iSpecies]   -= eta[iDim][0]*UnitNormal[iDim];
//      Av2[nSpecies+1][iSpecies] -= eta[iDim][1]*UnitNormal[iDim];
//      Av2[nSpecies+2][iSpecies] -= eta[iDim][2]*UnitNormal[iDim];
//      for (jDim = 0; jDim < nDim; jDim++)
//        Av2[nSpecies+3][iSpecies] += vel[jDim]*pi[iDim][jDim]*UnitNormal[iDim];
//    }
//  }
//  // X-momentum
//  Av2[nSpecies][nSpecies]     = 4.0/3.0*GInvRho[0]*nx + GInvRho[1]*ny + GInvRho[2]*nz;
//  Av2[nSpecies][nSpecies+1]   = GInvRho[0]*ny - 2.0/3.0*GInvRho[1]*nx;
//  Av2[nSpecies][nSpecies+2]   = GInvRho[0]*nz - 2.0/3.0*GInvRho[2]*nx;
//  // Y-momentum
//  Av2[nSpecies+1][nSpecies]   = GInvRho[1]*nx - 2.0/3.0*GInvRho[0]*ny;
//  Av2[nSpecies+1][nSpecies+1] = GInvRho[0]*nx + 4.0/3.0*GInvRho[1]*ny + GInvRho[2]*nz;
//  Av2[nSpecies+1][nSpecies+2] = GInvRho[1]*nz - 2.0/3.0*GInvRho[2]*ny;
//  // Z-momentum
//  Av2[nSpecies+2][nSpecies]   = GInvRho[2]*nx - 2.0/3.0*GInvRho[0]*nz;
//  Av2[nSpecies+2][nSpecies+1] = GInvRho[2]*ny - 2.0/3.0*GInvRho[1]*nz;
//  Av2[nSpecies+2][nSpecies+2] = GInvRho[0]*nx + GInvRho[1]*ny + 4.0/3.0*GInvRho[2]*nz;
//  // Energy
//  for (iDim = 0; iDim < nDim; iDim++) {
//    Av2[nSpecies+3][nSpecies] += (zeta[iDim][0] + 1.0/rho*tau[iDim][0])*UnitNormal[iDim];
//    Av2[nSpecies+3][nSpecies+1] += (zeta[iDim][1] + 1.0/rho*tau[iDim][1])*UnitNormal[iDim];
//    Av2[nSpecies+3][nSpecies+2] += (zeta[iDim][2] + 1.0/rho*tau[iDim][2])*UnitNormal[iDim];
//  }
//  Av2[nSpecies+3][nSpecies]   += vel[0]*GInvRho[0]*nx;
//  Av2[nSpecies+3][nSpecies+1] += vel[1]*GInvRho[1]*ny;
//  Av2[nSpecies+3][nSpecies+2] += vel[2]*GInvRho[2]*nz;
//  
//  /*--- Total energy viscous Jacobian (Av3) ---*/
//  for (iDim = 0; iDim < nDim; iDim++) {
//    // Species
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      Cvtrs = (3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
//      Av3[nSpecies+3][iSpecies] += (1.0/rhoCvtr * (vel[0]*PrimVar_Grad_i[VEL_INDEX][iDim]  +
//                                                   vel[1]*PrimVar_Grad_i[VEL_INDEX+1][iDim] +
//                                                   vel[2]*PrimVar_Grad_i[VEL_INDEX+2][iDim] -
//                                                   Cvtrs*PrimVar_Grad_i[T_INDEX][iDim])
//                                    - dTdU_i[iSpecies]/rhoCvtr*PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * UnitNormal[iDim];
//    }
//    // X-momentum
//    Av3[nSpecies+3][nSpecies]   += (-1.0/rhoCvtr*PrimVar_Grad_i[VEL_INDEX][iDim] -
//                                    dTdU_i[nSpecies]/rhoCvtr * PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * UnitNormal[iDim];
//    // Y-momentum
//    Av3[nSpecies+3][nSpecies+1] += (-1.0/rhoCvtr*PrimVar_Grad_i[VEL_INDEX+1][iDim] -
//                                    dTdU_i[nSpecies+1]/rhoCvtr * PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * UnitNormal[iDim];
//    // Z-momentum
//    Av3[nSpecies+3][nSpecies+2] += (-1.0/rhoCvtr*PrimVar_Grad_i[VEL_INDEX+2][iDim] -
//                                    dTdU_i[nSpecies+2]/rhoCvtr * PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * UnitNormal[iDim];
//    // Total energy
//    Av3[nSpecies+3][nSpecies+3] += -dTdU_i[nSpecies+3]/rhoCvtr * PrimVar_Grad_i[RHOCVTR_INDEX][iDim] * UnitNormal[iDim];
//    // Vibrational-electronic energy
//    Av3[nSpecies+3][nSpecies+4] += -dTdU_i[nSpecies+4]/rhoCvtr * PrimVar_Grad_i[RHOCVTR_INDEX][iDim] * UnitNormal[iDim];
//  }
//  
//  /*--- Compute the residual, Grad(Psi^T)dot(muk * Avk) ---*/
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nVar; jVar++) {
//      val_residual[iVar] += Psi_i[jVar]*(mu2*Av2[jVar][iVar] +
//                                         mu3*Av3[jVar][iVar] +
//                                         mu4*Av4[jVar][iVar]  ) * Volume;
//    }
//  }

  
  
//  cout << "Normal: " << Normal[0] << " " << Normal[1] << " " << Normal[2] << endl;
//  cout << "UnitNormal: " << UnitNormal[0] << " " << UnitNormal[1] << " " << UnitNormal[2] << endl;
//  cout << "Area; " << Area << endl;
//  cout << "Av2: " << endl;
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nVar; jVar++) {
//      cout << Av2[iVar][jVar] << "\t";
//    }
//    cout << endl;
//  }
//  cout << endl << endl;
//  cout << "Av3: " << endl;
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nVar; jVar++) {
//      cout << Av2[iVar][jVar] << "\t";
//    }
//    cout << endl;
//  }
//  cout << endl << endl;
//  cout << "Av4: " << endl;
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nVar; jVar++) {
//      cout << Av2[iVar][jVar] << "\t";
//    }
//    cout << endl;
//  }
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nDim; jVar++) {
//      cout << GPsi[iVar][jVar] << "\t";
//    }
//    cout << endl;
//  }
//  
//  cout << "Residual" << endl;
//  for (iVar = 0; iVar < nVar; iVar++)
//    cout << val_residual[iVar] << endl;
//  cin.get();
  
  
//	/*--- Required gradients of the flow variables, point j ---*/
//	for (iDim = 0; iDim < nDim; iDim++) {
//		/*--- grad density ---*/
//		GradDensity[iDim] = PrimVar_Grad_i[nDim+2][iDim];
//		/*--- grad (1/rho) ---*/
//		GradInvDensity[iDim] = -GradDensity[iDim]*invDensitysq;
//		/*--- Computation of the derivatives of P/(Density^2) ---*/
//		dPoDensity2[iDim] = (PrimVar_Grad_i[nVar-1][iDim]*Density - 2.0*GradDensity[iDim]*Pressure)*invDensitycube;
//		/*--- Abbreviations: alpha, beta, sigma_5_vec ---*/
//		alpha[iDim] = Gamma*mu_tot_2*GradInvDensity[iDim];
//		beta[iDim] = Gamma/Gamma_Minus_One*mu_tot_2*dPoDensity2[iDim];
//		Sigma_5_vec[iDim] = Gamma*mu_tot_2*PsiVar_Grad_i[nVar-1][iDim];
//	}
  
//	/*--- Definition of tensors and derivatives of velocity over density ---*/
//	double div_vel = 0, div_phi = 0, vel_gradpsi5 = 0;
//	for (iDim = 0; iDim < nDim; iDim++) {
//		div_vel += PrimVar_Grad_i[iDim+1][iDim];
//		div_phi += PsiVar_Grad_i[iDim+1][iDim];
//		vel_gradpsi5 += Velocity[iDim]*PsiVar_Grad_i[nVar-1][iDim];
//		for (jDim = 0; jDim < nDim; jDim++) {
//			sigma[iDim][jDim] = mu_tot_1*(PrimVar_Grad_i[iDim+1][jDim]+PrimVar_Grad_i[jDim+1][iDim]);
//			Sigma_phi[iDim][jDim] = mu_tot_1*(PsiVar_Grad_i[iDim+1][jDim]+PsiVar_Grad_i[jDim+1][iDim]);
//			Sigma_5_Tensor[iDim][jDim] = mu_tot_1*(Velocity[jDim]*PsiVar_Grad_i[nVar-1][iDim]+Velocity[iDim]*PsiVar_Grad_i[nVar-1][jDim]);
//			GradVel_o_Rho[iDim][jDim] = (PrimVar_Grad_i[iDim+1][jDim]*Density - Velocity[iDim]*GradDensity[jDim])*invDensitysq;
//		}
//	}
//	for (iDim = 0; iDim < nDim; iDim++) {
//		sigma[iDim][iDim] -= TWO3*mu_tot_1*div_vel;
//		Sigma_phi[iDim][iDim] -= TWO3*mu_tot_1*div_phi;
//		Sigma_5_Tensor[iDim][iDim] -= TWO3*mu_tot_1*vel_gradpsi5;
//	}
//	for (iDim = 0; iDim < nDim; iDim++)
//		for (jDim = 0; jDim < nDim; jDim++)
//			Sigma[iDim][jDim] = Sigma_phi[iDim][jDim] + Sigma_5_Tensor[iDim][jDim];
//  
//	/*--- Vector-Tensors products ---*/
//	double gradT_gradpsi5 = 0, sigma_gradpsi = 0, vel_sigma_gradpsi5 = 0;
//	for (iDim = 0; iDim < nDim; iDim++) {
//		gradT_gradpsi5 += PrimVar_Grad_i[0][iDim]*PsiVar_Grad_i[nVar-1][iDim];
//		for (jDim = 0; jDim < nDim; jDim++) {
//			sigma_gradpsi += sigma[iDim][jDim]*PsiVar_Grad_i[jDim+1][iDim];
//			vel_sigma_gradpsi5 += Velocity[iDim]*sigma[iDim][jDim]*PsiVar_Grad_i[nVar-1][jDim];
//		}
//	}
//  
//	/*--- Residuals ---*/
//	double alpha_gradpsi5 = 0, beta_gradpsi5 = 0, Sigma_gradvel_o_rho = 0, Sigma5_vel_gradvel = 0;
//	for (iDim = 0; iDim < nDim; iDim++) {
//		alpha_gradpsi5 += alpha[iDim]*PsiVar_Grad_i[nVar-1][iDim];
//		beta_gradpsi5 += beta[iDim]*PsiVar_Grad_i[nVar-1][iDim];
//		for (jDim = 0; jDim < nDim; jDim++) {
//			Sigma_gradvel_o_rho += Sigma[iDim][jDim]*GradVel_o_Rho[iDim][jDim];
//			Sigma5_vel_gradvel += Sigma_5_vec[iDim]*(Velocity[jDim]*PrimVar_Grad_i[jDim+1][iDim]);
//		}
//	}
//	val_residual[0] = (-vel_sigma_gradpsi5/Density - Sigma_gradvel_o_rho + 0.5*sq_vel*alpha_gradpsi5 -
//                     beta_gradpsi5 + Sigma5_vel_gradvel/Density) * Volume;
//	for (iDim = 0; iDim < nDim; iDim++)
//		for (jDim = 0; jDim < nDim; jDim++)
//			val_residual[iDim+1] = (sigma[iDim][jDim]*PsiVar_Grad_i[nVar-1][jDim]/Density +
//                              Sigma[iDim][jDim]*GradInvDensity[jDim] - Velocity[iDim]*alpha_gradpsi5 -
//                              Sigma_5_vec[jDim]*PrimVar_Grad_i[iDim+1][jDim]/Density) * Volume;
//	val_residual[nVar-1] = alpha_gradpsi5 * Volume;
  
}

//CSourceConservative_AdjFlow::CSourceConservative_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
//  
//	Gamma = config->GetGamma();
//	Gamma_Minus_One = Gamma - 1.0;
//  
//	Velocity = new double [nDim];
//	Residual_i = new double [nVar];
//	Residual_j = new double [nVar];
//	Mean_Residual = new double [nVar];
//  
//	Mean_PrimVar_Grad = new double* [nVar];
//	for (unsigned short iVar = 0; iVar < nVar; iVar++)
//		Mean_PrimVar_Grad[iVar] = new double [nDim];
//}
//
//CSourceConservative_AdjFlow::~CSourceConservative_AdjFlow(void) {
//	delete [] Mean_Residual;
//	delete [] Residual_j;
//	delete [] Residual_i;
//	delete [] Velocity;
//  
//	for (unsigned short iVar = 0; iVar < nVar; iVar++)
//		delete [] Mean_PrimVar_Grad[iVar];
//	delete [] Mean_PrimVar_Grad;
//}

void CSource_AdjTNE2::ComputeSourceConservative (double *val_residual,
                                                 CConfig *config) {
//	unsigned short iDim, jDim, iVar;
//	double rho, nu, Ji, fv1, fv2, Omega, Shat, dist_sq, Ji_2, Ji_3, one_o_oneplusJifv1;
//	double r, g, g_6, glim, dfw_g, dg_r, dr_nuhat, dr_Shat, Ms_coeff, invOmega;
//  
//	/*--- CLOUSURE CONSTANTS ---*/
//	double cv1_3 = 7.1*7.1*7.1;
//	double k2 = 0.41*0.41;
//	double cb1 = 0.1355;
//	double cw2 = 0.3;
//	double cw3_6 = pow(2.0,6.0);
//	double sigma = 2./3.;
//	double cb2 = 0.622;
//	double cw1 = cb1/k2+(1+cb2)/sigma;
//  
//	for (iVar = 0; iVar < nVar; iVar++) {
//		Residual_i[iVar] = 0.0;
//		Residual_j[iVar] = 0.0;
//	}
//  
//	/*--- iPoint ---*/
//  
//	/*--- Density and velocities ---*/
//	rho = U_i[0];
//	for (iDim = 0; iDim < nDim; iDim++)
//		Velocity[iDim] = U_i[iDim+1]/rho;
//  
//	/*--- Vorticity ---*/
//	Omega = (PrimVar_Grad_i[1][1]-PrimVar_Grad_i[2][0])*(PrimVar_Grad_i[1][1]-PrimVar_Grad_i[2][0]);
//	if (nDim == 3) Omega += (PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) +
//    (PrimVar_Grad_i[2][2]-PrimVar_Grad_i[3][1])*(PrimVar_Grad_i[2][2]-PrimVar_Grad_i[3][1]);
//	Omega = sqrt(Omega);
//	invOmega = 1.0/(Omega + TURB_EPS);
//	//	invOmega = min(1.0/Omega, max_invOmega);
//  
//	/*--- Compute Ms_coeff -> coming from partial derivatives ---*/
//	Ms_coeff = 0.0;
//	if (dist_i > 0) {
//		dist_sq = dist_i*dist_i;
//		nu = Laminar_Viscosity_i/rho;
//		Ji = TurbVar_i[0]/nu;
//		Ji_2 = Ji*Ji;
//		Ji_3 = Ji_2*Ji;
//		fv1 = Ji_3/(Ji_3+cv1_3);
//		one_o_oneplusJifv1 = 1.0/(1.0+Ji*fv1);
//		fv2 = 1.0 - Ji*one_o_oneplusJifv1;
//		Shat = max(Omega + TurbVar_i[0]*fv2/(k2*dist_sq),TURB_EPS);
//    
//		r = min(TurbVar_i[0]/(Shat*k2*dist_sq),10.);
//		g = r + cw2*(pow(r,6.)-r);
//		g_6 = pow(g,6.);
//		glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);
//    
//		dfw_g  = glim*cw3_6/(g_6+cw3_6);
//		dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0);
//		dr_nuhat = 1.0/(Shat*k2*dist_sq);
//		dr_Shat = -dr_nuhat*TurbVar_i[0]/Shat;
//    
//		Ms_coeff = (cb1*TurbVar_i[0]-cw1*TurbVar_i[0]*TurbVar_i[0]/dist_sq*dfw_g*dg_r*dr_Shat);
//	}
//	Ms_coeff *= TurbPsi_i[0]*invOmega/rho;
//  
//	/*--- Compute residual of iPoint ---*/
//	for (iDim = 0; iDim < nDim; iDim++) {
//		for (jDim = 0; jDim < nDim; jDim++) {
//			Residual_i[0] -= Ms_coeff*(Velocity[jDim]*PrimVar_Grad_i[jDim+1][iDim]*Normal[iDim] -
//                                 Velocity[jDim]*PrimVar_Grad_i[iDim+1][jDim]*Normal[iDim]);
//			Residual_i[iDim+1] += Ms_coeff*(PrimVar_Grad_i[iDim+1][jDim]*Normal[jDim] -
//                                      PrimVar_Grad_i[jDim+1][iDim]*Normal[jDim]);
//		}
//	}
//  
//	/*--- jPoint ---*/
//  
//	/*--- Density and velocities ---*/
//	rho = U_j[0];
//	for (iDim = 0; iDim < nDim; iDim++)
//		Velocity[iDim] = U_j[iDim+1]/rho;
//  
//	/*--- Vorticity ---*/
//	Omega = (PrimVar_Grad_j[1][1]-PrimVar_Grad_j[2][0])*(PrimVar_Grad_j[1][1]-PrimVar_Grad_j[2][0]);
//	if (nDim == 3) Omega += (PrimVar_Grad_j[1][2]-PrimVar_Grad_j[3][0])*(PrimVar_Grad_j[1][2]-PrimVar_Grad_j[3][0]) +
//    (PrimVar_Grad_j[2][2]-PrimVar_Grad_j[3][1])*(PrimVar_Grad_j[2][2]-PrimVar_Grad_j[3][1]);
//	Omega = sqrt(Omega);
//	invOmega = 1.0/(Omega + TURB_EPS);
//	//	invOmega = min(1.0/Omega, max_invOmega);
//  
//	/*--- Compute Ms_coeff -> coming from partial derivatives ---*/
//	Ms_coeff = 0.0;
//	if (dist_j > 0) {
//		dist_sq = dist_j*dist_j;
//		nu = Laminar_Viscosity_j/rho;
//		Ji = TurbVar_j[0]/nu;
//		Ji_2 = Ji*Ji;
//		Ji_3 = Ji_2*Ji;
//		fv1 = Ji_3/(Ji_3+cv1_3);
//		one_o_oneplusJifv1 = 1.0/(1.0+Ji*fv1);
//		fv2 = 1.0 - Ji*one_o_oneplusJifv1;
//		Shat = max(Omega + TurbVar_j[0]*fv2/(k2*dist_sq),TURB_EPS);
//    
//		r = min(TurbVar_j[0]/(Shat*k2*dist_sq),10.);
//		g = r + cw2*(pow(r,6.)-r);
//		g_6 = pow(g,6.);
//		glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);
//    
//		dfw_g  = glim*cw3_6/(g_6+cw3_6);
//		dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0);
//		dr_nuhat = 1.0/(Shat*k2*dist_sq);
//		dr_Shat = -dr_nuhat*TurbVar_j[0]/Shat;
//    
//		Ms_coeff = (cb1*TurbVar_j[0]-cw1*TurbVar_j[0]*TurbVar_j[0]/dist_sq*dfw_g*dg_r*dr_Shat);
//	}
//	Ms_coeff *= TurbPsi_j[0]*invOmega/rho;
//  
//	/*--- Compute residual of jPoint ---*/
//	for (iDim = 0; iDim < nDim; iDim++) {
//		for (jDim = 0; jDim < nDim; jDim++) {
//			Residual_j[0] -= Ms_coeff*(Velocity[jDim]*PrimVar_Grad_j[jDim+1][iDim]*Normal[iDim] -
//                                 Velocity[jDim]*PrimVar_Grad_j[iDim+1][jDim]*Normal[iDim]);
//			Residual_j[iDim+1] += Ms_coeff*(PrimVar_Grad_j[iDim+1][jDim]*Normal[jDim] -
//                                      PrimVar_Grad_j[jDim+1][iDim]*Normal[jDim]);
//		}
//	}
//  
//	/*--- MEAN RESIDUAL ---*/
//	for (iVar = 0; iVar < nVar; iVar++)
//		val_residual[iVar] = 0.5*(Residual_i[iVar] + Residual_j[iVar]);
  
}
