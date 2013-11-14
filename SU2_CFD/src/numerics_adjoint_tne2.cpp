/*!
 * \file numerics_adjoint_mean.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.8
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
  P      = new double* [nVar];
  invP   = new double* [nVar];
  PLPinv = new double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
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
  
  unsigned short iDim, iSpecies, iVar, jVar, kVar;
  double Area, ProjVel_i, ProjVel_j;
  
  /*--- SW flux: Fij = (F^+i + F^-j) ---*/
  
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
  
  /*--- Calculate projected velocties ---*/
  ProjVel_i = 0.0; ProjVel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVel_i += V_i[VEL_INDEX+iDim]*UnitNormal[iDim];
    ProjVel_j += V_j[VEL_INDEX+iDim]*UnitNormal[iDim];
  }
  
  /*--- Flow eigenvalues at i (Lambda+) --- */
  for (iSpecies = 0; iSpecies < nSpecies+nDim-1; iSpecies++)
    Lambda_i[iSpecies]      = 0.5*(ProjVel_i + fabs(ProjVel_i));
  Lambda_i[nSpecies+nDim-1] = 0.5*(     ProjVel_i + V_i[A_INDEX] +
                                   fabs(ProjVel_i + V_i[A_INDEX])  );
  Lambda_i[nSpecies+nDim]   = 0.5*(     ProjVel_i - V_i[A_INDEX] +
                                   fabs(ProjVel_i - V_i[A_INDEX])  );
  Lambda_i[nSpecies+nDim+1] = 0.5*(ProjVel_i + fabs(ProjVel_i));
  
  /*--- Compute projected P, invP, and Lambda ---*/
  GetPMatrix    (U_i, V_i, dPdU_i, UnitNormal, l, m, P   );
  GetPMatrix_inv(U_i, V_i, dPdU_i, UnitNormal, l, m, invP);
  
  /*--- Projected flux (f+) at i ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      PLPinv[iVar][jVar] = 0.0;
      /*--- Compute Proj_ModJac_Tensor = P x Lambda+ x inverse P ---*/
      for (kVar = 0; kVar < nVar; kVar++)
        PLPinv[iVar][jVar] += P[iVar][kVar]*Lambda_i[kVar]*invP[kVar][jVar];
    }
  }
  
  /*--- Calculate i+ fluxes ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual_i[iVar] += PLPinv[jVar][iVar]*Psi_i[jVar]*Area;
      val_Jacobian_ii[iVar][jVar] = PLPinv[jVar][iVar]*Area;
    }
  }
  
  /*--- Flow eigenvalues at i (Lambda-) --- */
  for (iSpecies = 0; iSpecies < nSpecies+nDim-1; iSpecies++)
    Lambda_i[iSpecies]      = 0.5*(ProjVel_i - fabs(ProjVel_i));
  Lambda_i[nSpecies+nDim-1] = 0.5*(     ProjVel_i + V_i[A_INDEX] -
                                   fabs(ProjVel_i + V_i[A_INDEX])  );
  Lambda_i[nSpecies+nDim]   = 0.5*(     ProjVel_i - V_i[A_INDEX] -
                                   fabs(ProjVel_i - V_i[A_INDEX])  );
  Lambda_i[nSpecies+nDim+1] = 0.5*(ProjVel_i - fabs(ProjVel_i));
  
  /*--- Projected flux (f-) at i ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      PLPinv[iVar][jVar] = 0.0;
      /*--- Compute Proj_ModJac_Tensor = P x Lambda+ x inverse P ---*/
      for (kVar = 0; kVar < nVar; kVar++)
        PLPinv[iVar][jVar] += P[iVar][kVar]*Lambda_i[kVar]*invP[kVar][jVar];
    }
  }
  
  /*--- Calculate i- fluxes ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual_i[iVar] += PLPinv[jVar][iVar]*Psi_j[jVar]*Area;
      val_Jacobian_ij[iVar][jVar] = PLPinv[jVar][iVar]*Area;
    }
  }
  
  //////// j -> i ////////
  /*--- Flow eigenvalues at j (Lambda+) --- */
  for (iVar = 0; iVar < nSpecies+nDim-1; iVar++)
    Lambda_j[iVar]          = -0.5*(ProjVel_j + fabs(ProjVel_j));
  Lambda_j[nSpecies+nDim-1] = -0.5*(     ProjVel_j + V_j[A_INDEX] +
                                    fabs(ProjVel_j + V_j[A_INDEX])  );
  Lambda_j[nSpecies+nDim]   = -0.5*(     ProjVel_j - V_j[A_INDEX] +
                                    fabs(ProjVel_j - V_j[A_INDEX])  );
  Lambda_j[nSpecies+nDim+1] = -0.5*(ProjVel_i + fabs(ProjVel_i));
  
  /*--- Compute projected P, invP, and Lambda ---*/
  GetPMatrix(U_j, V_j, dPdU_j, UnitNormal, l, m, P);
  GetPMatrix_inv(U_j, V_j, dPdU_j, UnitNormal, l, m, invP);
  
	/*--- Projected flux (f-) ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      PLPinv[iVar][jVar] = 0.0;
      /*--- Compute Proj_ModJac_Tensor = P x Lambda- x inverse P ---*/
      for (kVar = 0; kVar < nVar; kVar++)
        PLPinv[iVar][jVar] += P[iVar][kVar]*Lambda_j[kVar]*invP[kVar][jVar];
    }
  }
  
  /*--- Calculate j+ fluxes ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual_j[iVar] += PLPinv[jVar][iVar]*Psi_j[jVar]*Area;
      val_Jacobian_jj[iVar][jVar] = PLPinv[jVar][iVar]*Area;
    }
  }
  
  /*--- Flow eigenvalues at j (Lambda-) --- */
  for (iVar = 0; iVar < nSpecies+nDim-1; iVar++)
    Lambda_j[iVar]          = -0.5*(ProjVel_j - fabs(ProjVel_j));
  Lambda_j[nSpecies+nDim-1] = -0.5*(     ProjVel_j + V_j[A_INDEX] -
                                    fabs(ProjVel_j + V_j[A_INDEX])  );
  Lambda_j[nSpecies+nDim]   = -0.5*(     ProjVel_j - V_j[A_INDEX] -
                                    fabs(ProjVel_j - V_j[A_INDEX])  );
  Lambda_j[nSpecies+nDim+1] = -0.5*(ProjVel_j - fabs(ProjVel_j));
  
  /*--- Compute projected P, invP, and Lambda ---*/
  GetPMatrix(U_j, V_j, dPdU_j, UnitNormal, l, m, P);
  GetPMatrix_inv(U_j, V_j, dPdU_j, UnitNormal, l, m, invP);
  
	/*--- Projected flux (f-) ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      PLPinv[iVar][jVar] = 0.0;
      /*--- Compute Proj_ModJac_Tensor = P x Lambda- x inverse P ---*/
      for (kVar = 0; kVar < nVar; kVar++)
        PLPinv[iVar][jVar] += P[iVar][kVar]*Lambda_j[kVar]*invP[kVar][jVar];
    }
  }
  
  /*--- Calculate j- fluxes ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual_j[iVar] += PLPinv[jVar][iVar]*Psi_i[jVar]*Area;
      val_Jacobian_ji[iVar][jVar] = PLPinv[jVar][iVar]*Area;
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
