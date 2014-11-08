/*!
 * \file numerics_adjoint_mean.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.0.0 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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
  MeandTdU   = new double[nVar];
  MeandTvedU = new double[nVar];
  MeanEve    = new double[nSpecies];
  MeanCvve   = new double[nSpecies];
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
  delete [] MeandTdU;
  delete [] MeandTvedU;
  delete [] MeanEve;
  delete [] MeanCvve;
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
  
  double soundspeed_i, soundspeed_j;
  
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
  
  CreateBasis(UnitNormal);
  for (iVar = 0; iVar < nVar; iVar++)
    DiffPsi[iVar] = Psi_j[iVar] - Psi_i[iVar];
  
  
  /*--- Calculate eigenvalues of state i ---*/
  
  soundspeed_i = V_i[A_INDEX];
  ProjVel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ProjVel += V_i[VEL_INDEX+iDim];
  for (iVar = 0; iVar < nSpecies+nDim-1; iVar++)
    Lambda[iVar] = ProjVel;
  Lambda[nSpecies+nDim-1] = ProjVel + soundspeed_i;
  Lambda[nSpecies+nDim]   = ProjVel - soundspeed_i;
  Lambda[nSpecies+nDim+1] = ProjVel;
  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = fabs(Lambda[iVar]);
  
  // Calculate left and right eigenvector matrices
  GetPMatrix(U_i, V_i, dPdU_i, UnitNormal, l, m, P);
  GetPMatrix_inv(U_i, V_i, dPdU_i, UnitNormal, l, m, invP);
  
  // Calculate eigenvalue/eigenvector decomposition
  // |PLPinv| = P x |Lambda| x inverse P
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      PLPinv[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        PLPinv[iVar][jVar] += P[iVar][kVar]*Lambda[kVar]*invP[kVar][jVar];
    }
  }
  
  // Calculate the 'viscous' portion of the flux
  // 1/2*(P|Lam|P^-1)^T * (Uj - Ui)
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      val_residual_i[iVar] -= 0.5*PLPinv[jVar][iVar]*(Psi_i[jVar]-Psi_j[jVar])*Area;
  
  
  // Populate Jacobian matrices
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
  
  
  
  /*--- Calculate eigenvalues of state j ---*/
  
  soundspeed_j = V_j[A_INDEX];
  ProjVel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ProjVel += V_j[VEL_INDEX+iDim];
  for (iVar = 0; iVar < nSpecies+nDim-1; iVar++)
    Lambda[iVar] = ProjVel;
  Lambda[nSpecies+nDim-1] = ProjVel + soundspeed_j;
  Lambda[nSpecies+nDim]   = ProjVel - soundspeed_j;
  Lambda[nSpecies+nDim+1] = ProjVel;
  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = fabs(Lambda[iVar]);
  
  // Calculate left and right eigenvector matrices
  GetPMatrix(U_j, V_j, dPdU_j, UnitNormal, l, m, P);
  GetPMatrix_inv(U_j, V_j, dPdU_j, UnitNormal, l, m, invP);
  
  
  // Calculate eigenvalue/eigenvector decomposition
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
  
//  ////////// OLD
//  
//  /*--- Calculate mean variables ---*/
//  for (iVar = 0; iVar < nVar; iVar++)
//    MeanU[iVar] = 0.5*(U_i[iVar]+U_j[iVar]);
//  var->Cons2PrimVar(config, MeanU, MeanV, MeandPdU, MeandTdU,
//                    MeandTvedU, MeanEve, MeanCvve);
//  MeanSoundSpeed = MeanV[A_INDEX];
//  
//  for (iVar = 0; iVar < nVar; iVar++)
//    DiffPsi[iVar] = Psi_j[iVar] - Psi_i[iVar];
//  
//  ProjVel = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++)
//    ProjVel += MeanV[VEL_INDEX+iDim]*UnitNormal[iDim];
//  
//  /*--- Calculate eigenvalues of the interface state, ij ---*/
//  for (iVar = 0; iVar < nSpecies+nDim-1; iVar++)
//    Lambda[iVar] = ProjVel;
//  Lambda[nSpecies+nDim-1] = ProjVel + MeanSoundSpeed;
//  Lambda[nSpecies+nDim]   = ProjVel - MeanSoundSpeed;
//  Lambda[nSpecies+nDim+1] = ProjVel;
//  for (iVar = 0; iVar < nVar; iVar++)
//    Lambda[iVar] = fabs(Lambda[iVar]);
//  
//  /*--- Calculate left and right eigenvector matrices ---*/
//  CreateBasis(UnitNormal);
//  GetPMatrix(MeanU, MeanV, MeandPdU, UnitNormal, l, m, P);
//  GetPMatrix_inv(MeanU, MeanV, MeandPdU, UnitNormal, l, m, invP);
//  
//  /*--- Calculate eigenvalue/eigenvector decomposition ---*/
//  // |PLPinv| = P x |Lambda| x inverse P
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nVar; jVar++) {
//      PLPinv[iVar][jVar] = 0.0;
//      for (kVar = 0; kVar < nVar; kVar++)
//        PLPinv[iVar][jVar] += P[iVar][kVar]*Lambda[kVar]*invP[kVar][jVar];
//    }
//  }
//  
//  /*--- Calculate the 'viscous' portion of the flux ---*/
//  // 1/2*(P|Lam|P^-1)^T * (Uj - Ui)
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nVar; jVar++) {
//      val_residual_i[iVar] -= 0.5*PLPinv[jVar][iVar]*(Psi_i[jVar]-Psi_j[jVar])*Area;
//      val_residual_j[iVar] += 0.5*PLPinv[jVar][iVar]*(Psi_i[jVar]-Psi_j[jVar])*Area;
//    }
//  }
//  
//  
//  /*--- Populate Jacobian matrices ---*/
//  // Note: Ai/j calculated using 'Normal', but PLPinv computed using UnitNormal.
//  //       Only multiply PLP by area to properly account for integration.
//  if (implicit) {
//    for (iVar = 0; iVar < nVar; iVar++) {
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_ii[iVar][jVar] = 0.5*Ai[jVar][iVar] - 0.5*PLPinv[jVar][iVar]*Area;
//        val_Jacobian_ij[iVar][jVar] = 0.5*Ai[jVar][iVar] + 0.5*PLPinv[jVar][iVar]*Area;
//        val_Jacobian_ji[iVar][jVar] = -0.5*Aj[jVar][iVar] + 0.5*PLPinv[jVar][iVar]*Area;
//        val_Jacobian_jj[iVar][jVar] = -0.5*Aj[jVar][iVar] - 0.5*PLPinv[jVar][iVar]*Area;
//      }
//    }
//  }
//  ////////// OLD
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
  Pi      = new double* [nVar];
  Pj      = new double* [nVar];
  invPi   = new double* [nVar];
  invPj   = new double* [nVar];
  PLPinvi = new double* [nVar];
  PLPinvj = new double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Pi[iVar]      = new double [nVar];
    Pj[iVar]      = new double [nVar];
    invPi[iVar]   = new double [nVar];
    invPj[iVar]   = new double [nVar];
    PLPinvi[iVar] = new double [nVar];
    PLPinvj[iVar] = new double [nVar];
  }
  
}

CUpwSW_AdjTNE2::~CUpwSW_AdjTNE2(void) {
  
  delete [] UnitNormal;
  delete [] DiffPsi;
	delete [] Lambda_i;
  delete [] Lambda_j;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] Pi[iVar];
    delete [] Pj[iVar];
    delete [] invPi[iVar];
    delete [] invPj[iVar];
    delete [] PLPinvi[iVar];
    delete [] PLPinvj[iVar];
  }
  delete [] Pi;
  delete [] Pj;
  delete [] invPi;
  delete [] invPj;
  delete [] PLPinvi;
  delete [] PLPinvj;
}

void CUpwSW_AdjTNE2::ComputeResidual (double *val_residual_i,
                                      double *val_residual_j,
                                      double **val_Jacobian_ii,
                                      double **val_Jacobian_ij,
                                      double **val_Jacobian_ji,
                                      double **val_Jacobian_jj,
                                      CConfig *config) {

  unsigned short iDim, iVar, jVar, kVar;
  double uni, unj;
  
  /*--- Steger-Warming Flux Vector Splitting Method ---*/
  //
  // Phi^c = Psi^T (A^c \cdot n)
  // \hat{Phi}^c = Psi_i^T (Ai^c \cdot n)^+ + Psi_j^T (Aj^c \cdot n)
  //
  // Where A+/- = P (Lambda+/-) Pinv

  
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

  
  /*---+++ Calculate the i->j flux +++---*/
  
  
  /*--- Calculate auxiliary quantities ---*/
  uni = 0; unj = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    uni += V_i[VEL_INDEX+iDim]*UnitNormal[iDim];
    unj += V_j[VEL_INDEX+iDim]*UnitNormal[iDim];
  }
  
  /*--- Calculate the eigenvalues ---*/
  for (iVar = 0; iVar < nSpecies+nDim-1; iVar++) {
    Lambda_i[iVar] = uni;
    Lambda_j[iVar] = unj;
  }
  Lambda_i[nSpecies+nDim-1] = uni + V_i[A_INDEX];
  Lambda_j[nSpecies+nDim-1] = unj + V_j[A_INDEX];
  Lambda_i[nSpecies+nDim]   = uni - V_i[A_INDEX];
  Lambda_j[nSpecies+nDim]   = unj - V_j[A_INDEX];
  Lambda_i[nSpecies+nDim+1] = uni;
  Lambda_j[nSpecies+nDim+1] = unj;
  
  /*--- Calculate left and right eigenvector matrices ---*/
  GetPMatrix    (U_i, V_i, dPdU_i, UnitNormal, l, m, Pi   );
  GetPMatrix_inv(U_i, V_i, dPdU_i, UnitNormal, l, m, invPi);
  GetPMatrix    (U_j, V_j, dPdU_j, UnitNormal, l, m, Pj   );
  GetPMatrix_inv(U_j, V_j, dPdU_j, UnitNormal, l, m, invPj);

  /*--- Calculate Ai+ and Aj- ---*/
  // |PLPinv| = P x |Lambda| x inverse P
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      PLPinvi[iVar][jVar] = 0.0;
      PLPinvj[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++) {
        PLPinvi[iVar][jVar] += Pi[iVar][kVar]
                             * (Lambda_i[kVar] + fabs(Lambda_i[kVar]))
                             * invPi[kVar][jVar];
        PLPinvj[iVar][jVar] += Pj[iVar][kVar]
                             * (Lambda_j[kVar] - fabs(Lambda_j[kVar]))
                             * invPj[kVar][jVar];
      }
    }
  }
  
  /*--- Populate the flux residual: Ai+^T Psi_i + Aj-^T Psi_j ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual_i[iVar] = 0;
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual_i[iVar] += (PLPinvi[jVar][iVar]*Psi_i[jVar] +
                               PLPinvj[jVar][iVar]*Psi_j[jVar]  )*Area;
    }
  }
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_ii[iVar][jVar] = PLPinvi[jVar][iVar]*Area;
        val_Jacobian_ij[iVar][jVar] = PLPinvj[jVar][iVar]*Area;
      }
  }
  
  
  /*---+++ Calculate the j->i flux +++---*/
  
  
  /*--- Reverse the direction of the normal vector ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    UnitNormal[iDim] = -UnitNormal[iDim];
    uni = -uni;
    unj = -unj;
  }
  CreateBasis(UnitNormal);
  
  /*--- Calculate the eigenvalues ---*/
  for (iVar = 0; iVar < nSpecies+nDim-1; iVar++) {
    Lambda_i[iVar] = uni;
    Lambda_j[iVar] = unj;
  }
  Lambda_i[nSpecies+nDim-1] = uni + V_i[A_INDEX];
  Lambda_j[nSpecies+nDim-1] = unj + V_j[A_INDEX];
  Lambda_i[nSpecies+nDim]   = uni - V_i[A_INDEX];
  Lambda_j[nSpecies+nDim]   = unj - V_j[A_INDEX];
  Lambda_i[nSpecies+nDim+1] = uni;
  Lambda_j[nSpecies+nDim+1] = unj;
  
  /*--- Calculate left and right eigenvector matrices ---*/
  GetPMatrix    (U_i, V_i, dPdU_i, UnitNormal, l, m, Pi   );
  GetPMatrix_inv(U_i, V_i, dPdU_i, UnitNormal, l, m, invPi);
  GetPMatrix    (U_j, V_j, dPdU_j, UnitNormal, l, m, Pj   );
  GetPMatrix_inv(U_j, V_j, dPdU_j, UnitNormal, l, m, invPj);
  
  /*--- Calculate Aj+ and Ai- ---*/
  // |PLPinv| = P x |Lambda| x inverse P
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      PLPinvi[iVar][jVar] = 0.0;
      PLPinvj[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++) {
        PLPinvi[iVar][jVar] += Pi[iVar][kVar]
                             * (Lambda_i[kVar] - fabs(Lambda_i[kVar]))
                             * invPi[kVar][jVar];
        PLPinvj[iVar][jVar] += Pj[iVar][kVar]
                             * (Lambda_j[kVar] + fabs(Lambda_j[kVar]))
                             * invPj[kVar][jVar];
      }
    }
  }
  
  /*--- Populate the flux residual: Aj+^T Psi_j + Ai-^T Psi_i ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual_j[iVar] = 0;
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual_j[iVar] += (PLPinvj[jVar][iVar]*Psi_j[jVar] +
                               PLPinvi[jVar][iVar]*Psi_i[jVar]  )*Area;
    }
  }
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_ji[iVar][jVar] = PLPinvi[jVar][iVar]*Area;
        val_Jacobian_jj[iVar][jVar] = PLPinvj[jVar][iVar]*Area;
      }
  }
  
  
  ////////////////// OLD //////////////////
//  unsigned short iDim, iVar, jVar, kVar;
//  double Area, ProjVel_i, ProjVel_j;
//
//  /*--- Initialize the residuals (and Jacobians) ---*/
//  for (iVar = 0; iVar < nVar; iVar++) {
//    val_residual_i[iVar] = 0.0;
//    val_residual_j[iVar] = 0.0;
//  }
//  if (implicit) {
//    for (iVar = 0; iVar < nVar; iVar++)
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_ii[iVar][jVar] = 0.0;
//        val_Jacobian_ij[iVar][jVar] = 0.0;
//        val_Jacobian_ji[iVar][jVar] = 0.0;
//        val_Jacobian_jj[iVar][jVar] = 0.0;
//      }
//  }
//  
//  /*--- Calculate geometric quantities ---*/
//  Area = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++)
//    Area += Normal[iDim]*Normal[iDim];
//  Area = sqrt(Area);
//  for (iDim = 0; iDim < nDim; iDim++)
//    UnitNormal[iDim] = Normal[iDim]/Area;
//  CreateBasis(UnitNormal);
//  
//  /*--- Calculate inviscid projected Jacobians ---*/
//  GetInviscidProjJac(U_i, V_i, dPdU_i, Normal, 1.0, Ai);
//  GetInviscidProjJac(U_j, V_j, dPdU_j, Normal, 1.0, Aj);
//  
//  /*--- Inviscid portion of the flux, A^T*Psi ---*/
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nVar; jVar++) {
//      val_residual_i[iVar] +=  0.5*Ai[jVar][iVar]*(Psi_i[jVar]+Psi_j[jVar]);
//      val_residual_j[iVar] += -0.5*Aj[jVar][iVar]*(Psi_i[jVar]+Psi_j[jVar]);
//    }
//  }
//  
//  for (iVar = 0; iVar < nVar; iVar++)
//    DiffPsi[iVar] = Psi_j[iVar] - Psi_i[iVar];
//  
//  ProjVel_i = 0.0;
//  ProjVel_j = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++) {
//    ProjVel_i += V_i[VEL_INDEX+iDim]*UnitNormal[iDim];
//    ProjVel_j += V_j[VEL_INDEX+iDim]*UnitNormal[iDim];
//  }
//  
//  /*--- Calculate eigenvalues of the i & j states ---*/
//  for (iVar = 0; iVar < nSpecies+nDim-1; iVar++) {
//    Lambda_i[iVar] = ProjVel_i;
//    Lambda_j[iVar] = ProjVel_j;
//  }
//  Lambda_i[nSpecies+nDim-1] = ProjVel_i + V_i[A_INDEX];
//  Lambda_j[nSpecies+nDim-1] = ProjVel_j + V_j[A_INDEX];
//  Lambda_i[nSpecies+nDim]   = ProjVel_i - V_i[A_INDEX];
//  Lambda_j[nSpecies+nDim]   = ProjVel_j - V_j[A_INDEX];
//  Lambda_i[nSpecies+nDim+1] = ProjVel_i;
//  Lambda_j[nSpecies+nDim+1] = ProjVel_j;
//  for (iVar = 0; iVar < nVar; iVar++) {
//    Lambda_i[iVar] = fabs(Lambda_i[iVar]);
//    Lambda_j[iVar] = fabs(Lambda_j[iVar]);
//  }
//  
//  /*--- Calculate left and right eigenvector matrices ---*/
//  GetPMatrix(U_i, V_i, dPdU_i, UnitNormal, l, m, P);
//  GetPMatrix_inv(U_i, V_i, dPdU_i, UnitNormal, l, m, invP);
//  
//  /*--- Calculate eigenvalue/eigenvector decomposition for i ---*/
//  // |PLPinv| = P x |Lambda| x inverse P
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nVar; jVar++) {
//      PLPinv[iVar][jVar] = 0.0;
//      for (kVar = 0; kVar < nVar; kVar++)
//        PLPinv[iVar][jVar] += P[iVar][kVar]*Lambda_i[kVar]*invP[kVar][jVar];
//    }
//  }
//  
//  /*--- Calculate the 'viscous' portion of the flux at i ---*/
//  // 1/2*(P|Lam|P^-1)^T * (Uj - Ui)
//  for (iVar = 0; iVar < nVar; iVar++)
//    for (jVar = 0; jVar < nVar; jVar++)
//      val_residual_i[iVar] -= 0.5*PLPinv[jVar][iVar]*(Psi_i[jVar]-Psi_j[jVar])*Area;
//  
//  /*--- Populate Jacobian matrices ---*/
//  // Note: Ai/j calculated using 'Normal', but PLPinv computed using UnitNormal.
//  //       Only multiply PLP by area to properly account for integration.
//  if (implicit) {
//    for (iVar = 0; iVar < nVar; iVar++) {
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_ii[iVar][jVar] = 0.5*Ai[jVar][iVar] - 0.5*PLPinv[jVar][iVar]*Area;
//        val_Jacobian_ij[iVar][jVar] = 0.5*Ai[jVar][iVar] + 0.5*PLPinv[jVar][iVar]*Area;
//      }
//    }
//  }
//  
//  /*--- Calculate left and right eigenvector matrices for j ---*/
//  GetPMatrix(U_j, V_j, dPdU_j, UnitNormal, l, m, P);
//  GetPMatrix_inv(U_j, V_j, dPdU_j, UnitNormal, l, m, invP);
//  
//  /*--- Calculate eigenvalue/eigenvector decomposition for i ---*/
//  // |PLPinv| = P x |Lambda| x inverse P
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nVar; jVar++) {
//      PLPinv[iVar][jVar] = 0.0;
//      for (kVar = 0; kVar < nVar; kVar++)
//        PLPinv[iVar][jVar] += P[iVar][kVar]*Lambda_j[kVar]*invP[kVar][jVar];
//    }
//  }
//  
//  /*--- Calculate the 'viscous' portion of the flux at j ---*/
//  // 1/2*(P|Lam|P^-1)^T * (Uj - Ui)
//  for (iVar = 0; iVar < nVar; iVar++)
//    for (jVar = 0; jVar < nVar; jVar++)
//      val_residual_j[iVar] += 0.5*PLPinv[jVar][iVar]*(Psi_i[jVar]-Psi_j[jVar])*Area;
//  
//  /*--- Populate Jacobian matrices ---*/
//  // Note: Ai/j calculated using 'Normal', but PLPinv computed using UnitNormal.
//  //       Only multiply PLP by area to properly account for integration.
//  if (implicit) {
//    for (iVar = 0; iVar < nVar; iVar++) {
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_ji[iVar][jVar] = -0.5*Aj[jVar][iVar] + 0.5*PLPinv[jVar][iVar]*Area;
//        val_Jacobian_jj[iVar][jVar] = -0.5*Aj[jVar][iVar] - 0.5*PLPinv[jVar][iVar]*Area;
//      }
//    }
//  }
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
  
	sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j)) /
            (double(Neighbor_i)*double(Neighbor_j));
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


CAvgGrad_AdjTNE2::CAvgGrad_AdjTNE2(unsigned short val_nDim,
                                   unsigned short val_nVar,
                                   CConfig *config) : CNumerics(val_nDim,
                                                                val_nVar,
                                                                config) {
	unsigned short iDim, iSpecies;
  
  implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  
  nDim         = val_nDim;
  nSpecies     = config->GetnSpecies();
  nVar         = val_nVar;
  
	vel   = new double[nDim];
  vel_i = new double[nDim];
  vel_j = new double[nDim];
  hs_i  = new double[nSpecies];
  hs_j  = new double[nSpecies];
  
  DdYk   = new double[nSpecies];
  dYdrs  = new double*[nSpecies];
  dJddrs = new double*[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    dYdrs[iSpecies]  = new double[nSpecies];
    dJddrs[iSpecies] = new double[nSpecies];
  }
  
  GPsirsdotn = new double [nSpecies];
  Mean_GradPsirs = new double *[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Mean_GradPsirs[iSpecies] = new double[nDim];
  
	Mean_GradPhi = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Mean_GradPhi[iDim] = new double [nDim];
  
	Mean_GradPsiE = new double [nDim];
  Mean_GradPsiEve = new double [nDim];
	
  Edge_Vector = new double [nDim];
  
  SigmaVel = new double*[nDim];
  SigmaPhi  = new double*[nDim];
  SigmaPsiE = new double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    SigmaVel[iDim] = new double[nDim];
    SigmaPhi[iDim]  = new double[nDim];
    SigmaPsiE[iDim] = new double[nDim];
  }
}

CAvgGrad_AdjTNE2::~CAvgGrad_AdjTNE2(void) {
  unsigned short iDim, iSpecies;
  
  delete [] vel;
  delete [] vel_i;
  delete [] vel_j;
  delete [] hs_i;
  delete [] hs_j;
  
	delete [] Edge_Vector;
  
  delete [] DdYk;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    delete [] dJddrs[iSpecies];
    delete [] dYdrs[iSpecies];
  }
  delete [] dJddrs;
  delete [] dYdrs;
  
  delete [] GPsirsdotn;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    delete [] Mean_GradPsirs[iSpecies];
  delete [] Mean_GradPsirs;
  
	delete [] Mean_GradPsiE;
  delete [] Mean_GradPsiEve;
  
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Mean_GradPhi[iDim];
  delete [] Mean_GradPhi;
  
  for (iDim = 0; iDim < nDim; iDim++) {
    delete [] SigmaVel[iDim];
    delete [] SigmaPhi[iDim];
    delete [] SigmaPsiE[iDim];
  }
  delete [] SigmaVel;
  delete [] SigmaPhi;
  delete [] SigmaPsiE;
  
}

void CAvgGrad_AdjTNE2::ComputeResidual(double *val_residual_i,
                                       double *val_residual_j,
                                       double **val_Jacobian_ii,
                                       double **val_Jacobian_ij,
                                       double **val_Jacobian_ji,
                                       double **val_Jacobian_jj,
                                       CConfig *config) {

  
  unsigned short iDim, jDim, iSpecies, jSpecies, iVar, jVar;
  double mu_i, mu_j, ktr_i, ktr_j, kve_i, kve_j, *D_i, *D_j;
  double rho, rho_i, rho_j, un_i, un_j, u2_i, u2_j, Ys;
  double GdotPhi, GPsiEdotVel, GPsiEdotn, GPsiEvedotn, SigVelGPhi;
  double dij, theta, thetax, thetay, thetaz, etax, etay, etaz;
  
  /*--- Initialize residuals ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual_i[iVar] = 0.0;
    val_residual_j[iVar] = 0.0;
    for (jVar = 0; jVar < nVar; jVar++) {
      val_Jacobian_ii[iVar][jVar] = 0.0;
      val_Jacobian_ij[iVar][jVar] = 0.0;
      val_Jacobian_ji[iVar][jVar] = 0.0;
      val_Jacobian_jj[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Calculate geometric quantities ---*/
  Area = 0.0;
  dij = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Area += Normal[iDim]*Normal[iDim];
    dij  += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
  }
  Area = sqrt(Area);
  dij  = sqrt(dij);
  theta = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    UnitNormal[iDim] = Normal[iDim]/Area;
    theta += UnitNormal[iDim]*UnitNormal[iDim];
  }
  if (nDim == 2) {
    thetax = theta + (UnitNormal[0]*UnitNormal[0])/3.0;
    thetay = theta + (UnitNormal[1]*UnitNormal[1])/3.0;
    etaz   = UnitNormal[0]*UnitNormal[1]/3.0;
  } else {
    thetax = theta + (UnitNormal[0]*UnitNormal[0])/3.0;
    thetay = theta + (UnitNormal[1]*UnitNormal[1])/3.0;
    thetaz = theta + (UnitNormal[2]*UnitNormal[2])/3.0;
    etax   = UnitNormal[1]*UnitNormal[2]/3.0;
    etay   = UnitNormal[0]*UnitNormal[2]/3.0;
    etaz   = UnitNormal[0]*UnitNormal[1]/3.0;
  }
  
  /*--- Get flow state (Rename for convenience) ---*/
  D_i   = Diffusion_Coeff_i;
  D_j   = Diffusion_Coeff_j;
  mu_i  = Laminar_Viscosity_i;
  mu_j  = Laminar_Viscosity_j;
  ktr_i = Thermal_Conductivity_i;
  ktr_j = Thermal_Conductivity_j;
  kve_i = Thermal_Conductivity_ve_i;
  kve_j = Thermal_Conductivity_ve_j;
  rho_i = V_i[RHO_INDEX];
  rho_j = V_j[RHO_INDEX];
  rho   = 0.5*(rho_i+rho_j);
  un_i = 0.0; un_j = 0.0;
  u2_i = 0.0; u2_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    vel_i[iDim] = V_i[VEL_INDEX+iDim];
    vel_j[iDim] = V_j[VEL_INDEX+iDim];
    vel[iDim] = 0.5*(vel_i[iDim]+vel_j[iDim]);
    un_i += vel_i[iDim]*UnitNormal[iDim];
    un_j += vel_j[iDim]*UnitNormal[iDim];
    u2_i += vel_i[iDim]*vel_i[iDim];
    u2_j += vel_j[iDim]*vel_j[iDim];
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    hs_i[iSpecies]  = var->CalcHs (config, V_i[T_INDEX], eve_i[iSpecies], iSpecies);
    hs_j[iSpecies]  = var->CalcHs (config, V_j[T_INDEX], eve_j[iSpecies], iSpecies);
  }
  
  /*--- Calculate mean gradients ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      Mean_GradPsirs[iSpecies][iDim] = 0.5*(PsiVar_Grad_i[iSpecies][iDim] +
                                            PsiVar_Grad_j[iSpecies][iDim]  );
    
    Mean_GradPsiE[iDim]   =  0.5*(PsiVar_Grad_i[nSpecies+nDim][iDim] +
                                  PsiVar_Grad_j[nSpecies+nDim][iDim]  );
    Mean_GradPsiEve[iDim] = 0.5*(PsiVar_Grad_i[nSpecies+nDim+1][iDim] +
                                 PsiVar_Grad_j[nSpecies+nDim+1][iDim]  );
		for (jDim = 0; jDim < nDim; jDim++)
      Mean_GradPhi[iDim][jDim] =  0.5*(PsiVar_Grad_i[nSpecies+iDim][jDim] +
                                       PsiVar_Grad_j[nSpecies+iDim][jDim]  );
  }
  
  /*--- Calculate auxiliary quantities for SigmaPhi ---*/
  GdotPhi     = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    GdotPhi     += Mean_GradPhi[iDim][iDim];
  
  /*--- Project mean gradient of PsiRs, PsiE & PsiEve into normal ---*/
  GPsiEdotn   = 0.0;
  GPsiEvedotn = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    GPsirsdotn[iSpecies] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      GPsirsdotn[iSpecies] +=Mean_GradPsirs[iSpecies][iDim]*Normal[iDim];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    GPsiEdotn   += Mean_GradPsiE[iDim]*Normal[iDim];
    GPsiEvedotn += Mean_GradPsiEve[iDim]*Normal[iDim];
  }
  
  /*--- Initialize SigmaPhi ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      SigmaPhi[iDim][jDim] = 0.0;
  
  /*--- Calculate SigmaPhi ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      SigmaPhi[iDim][jDim] += Mean_GradPhi[iDim][jDim] +
                              Mean_GradPhi[jDim][iDim];
    }
    SigmaPhi[iDim][iDim]  -= 2.0/3.0*GdotPhi;
  }
  
  /*---+++ Residual at node i +++---*/
  
  // k = 1
  /*--- Calculate auxiliary quantities ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dYdrs[iSpecies][jSpecies]  = 0.0;
      dJddrs[iSpecies][jSpecies] = 0.0;
    }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dYdrs[iSpecies][jSpecies] += 1.0/rho_i*(-V_i[RHOS_INDEX+iSpecies]/rho_i);
    }
    dYdrs[iSpecies][iSpecies] += 1.0/rho_i;
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    DdYk[iSpecies] = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      DdYk[iSpecies] += rho_i*D_i[jSpecies]*dYdrs[jSpecies][iSpecies];
    }
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Ys = V_i[RHOS_INDEX+iSpecies]/rho_i;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dJddrs[iSpecies][jSpecies] = -rho_i*D_i[iSpecies]*dYdrs[iSpecies][jSpecies]
                                 + Ys*DdYk[jSpecies];
    }
  }
  
//  //////////////////// DEBUG ////////////////////
//  double UnitNormal[3], tmp, tmp2;
//  double *Fv_old, *Fv_new, d;
//  double **GU, **GV, **Dxx;
//  UnitNormal[0] = 0.0;
//  UnitNormal[1] = 1.0;
//  UnitNormal[2] = 0.0;
//  Fv_new = new double[nVar];
//  Fv_old = new double[nVar];
//  
//  Dxx = new double*[nVar];
//  for (iVar = 0; iVar < nVar; iVar++)
//    Dxx[iVar] = new double[nVar];
//  
//  GU = ConsVar_Grad_i;
//  GV = PrimVar_Grad_i;
//  
//  for (iVar = 0; iVar < nVar; iVar++)
//    for (jVar = 0; jVar < nVar; jVar++)
//      Dxx[iVar][jVar] = 0.0;
//  
//  /*--- Calculate the k=1 residual at i ---*/
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//      Dxx[iSpecies][jSpecies] +=(-dJddrs[iSpecies][jSpecies]);
//      Dxx[nSpecies+nDim][iSpecies] +=(-dJddrs[jSpecies][iSpecies])*hs_i[jSpecies];
//      Dxx[nSpecies+nDim+1][iSpecies] +=(-dJddrs[jSpecies][iSpecies])*eve_i[jSpecies];
//    }
//  }
//  
//  cout << "C: " << endl;
//  for (iSpecies = 0; iSpecies < nVar; iSpecies++) {
//    for (jSpecies = 0; jSpecies < nVar; jSpecies++)
//      cout << Dxx[jSpecies][iSpecies] << "\t";
//    cout << endl;
//  }
//  
//  
//  cout << endl << "FD: " << endl;
//  // finite difference gradient
//  iDim = 1;
//  for (iVar = 0; iVar < nVar; iVar++) {
//    // set displacement value
//    
//    d = 0.0001*GU[iVar][iDim];
////    cout << "d: " << d << endl;
////    cin.get();
//    
//    // calculate viscous flux
//    GetViscousProjFlux(V_i, GV, UnitNormal, D_i, mu_i, ktr_i, kve_i, config);
//    
//    // copy solution
//    for (jVar = 0; jVar < nVar; jVar++)
//      Fv_old[jVar] = Proj_Flux_Tensor[jVar];
//    
//    // perturb solution
//    GU[iVar][iDim] += d;
//    var->GradCons2GradPrimVar(config, U_i, V_i, GU, GV);
//    
//    // calculate viscous flux
//    GetViscousProjFlux(V_i, GV, UnitNormal, D_i, mu_i, ktr_i, kve_i, config);
//    
//    // copy solution
//    for (jVar = 0; jVar < nVar; jVar++)
//      Fv_new[jVar] = Proj_Flux_Tensor[jVar];
//    
//    // return solution to original value
//    GU[iVar][iDim] -= d;
//    var->GradCons2GradPrimVar(config, U_i, V_i, GU, GV);
//    
//    // display FD gradient
//    for (jVar = 0; jVar < nVar; jVar++)
//      cout << (Fv_new[jVar]-Fv_old[jVar])/d << "\t";
//    cout << endl;
//  }
//  
//  cin.get();
//  
//  //////////////////// DEBUG ////////////////////
  
  
  
  /*--- Calculate the k=1 residual at i ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      val_residual_i[iSpecies] +=
          GPsirsdotn[jSpecies]*(-dJddrs[jSpecies][iSpecies])
        + GPsiEdotn           *(-dJddrs[jSpecies][iSpecies])*hs_i[jSpecies]
        + GPsiEvedotn         *(-dJddrs[jSpecies][iSpecies])*eve_i[jSpecies];
    }
  }
  
  // k = 2
  /*--- Calculate auxiliary quantities for SigmaPsiE ---*/
  GPsiEdotVel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    GPsiEdotVel += Mean_GradPsiE[iDim]*vel_i[iDim];
  
  /*--- Initialize SigmaPsiE & SigmaVel ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++) {
      SigmaPsiE[iDim][jDim] = 0.0;
      SigmaVel[iDim][jDim] = 0.0;
    }
  
  /*--- Calculate SigmaPsiE & SigmaVel ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      SigmaPsiE[iDim][jDim] += Mean_GradPsiE[iDim]*vel_i[jDim] +
                               Mean_GradPsiE[jDim]*vel_i[iDim];
      SigmaVel[iDim][jDim] += vel_i[iDim]*UnitNormal[jDim] + vel_i[jDim]*UnitNormal[iDim];
    }
    SigmaPsiE[iDim][iDim] -= 2.0/3.0*GPsiEdotVel;
    SigmaVel[iDim][iDim] -= 2.0/3.0*un_i;
  }
  
  /*--- Calculate the k=2 residual at i (SigmaPhi + SigmaPsiE) dot n ---*/
  SigVelGPhi = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      SigVelGPhi += SigmaVel[iDim][jDim]*Mean_GradPhi[jDim][iDim];

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    val_residual_i[iSpecies] += -mu_i/rho_i*(SigVelGPhi*Area +
                                             u2_i*GPsiEdotn +
                                             un_i*GPsiEdotVel*Area/3.0);
  
  
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      val_residual_i[nSpecies+iDim] += mu_i/rho_i*(SigmaPhi[iDim][jDim] +
                                                   SigmaPsiE[iDim][jDim]  )
                                     * Normal[jDim];
    }
  }
  
  // k = 3
  /*--- Calculate the k=3 residual at i dT/dU * (GradPsiE dot n) ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual_i[iVar] += ktr_i*dTdU_i[iVar]*GPsiEdotn;
  
  // k = 4
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual_i[iVar] += kve_i*dTvedU_i[iVar]*(GPsiEvedotn+GPsiEdotn);
  
  /*--- Calculate Jacobians for implicit time-stepping ---*/
  if (implicit) {
    
    
    /*--- Jacobian from k = 1 viscous flux ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        // species density
        val_Jacobian_ij[iSpecies][jSpecies] += -theta/dij*dJddrs[jSpecies][iSpecies]*Area;
        // total energy
        val_Jacobian_ij[iSpecies][nSpecies+nDim] += -theta/dij*dJddrs[jSpecies][iSpecies]*hs_i[jSpecies]*Area;
        // vib.-el. energy
        val_Jacobian_ij[iSpecies][nSpecies+nDim] += -theta/dij*dJddrs[jSpecies][iSpecies]*eve_i[jSpecies]*Area;
      }
    }
    
    /*--- Jacobian from k = 2 viscous flux ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        for (jDim = 0; jDim < nDim; jDim++) {
          val_Jacobian_ij[iSpecies][nSpecies+iDim] += -mu_i/rho_i*SigmaVel[jDim][iDim]*UnitNormal[jDim]/dij*Area;
        }
      }
      val_Jacobian_ij[iSpecies][nSpecies+nDim] += -mu_i/rho_i*(u2_i*theta/dij - un_i*un_i/(3*dij))*Area;
    }
    if (nDim == 2) {
      // x-momentum
      val_Jacobian_ij[nSpecies][nSpecies]     += mu_i/(rho_i*dij) * thetax * Area;
      val_Jacobian_ij[nSpecies][nSpecies+1]   += mu_i/(rho_i*dij) * etaz   * Area;
      val_Jacobian_ij[nSpecies][nSpecies+2]   += mu_i/(rho_i*dij) *
                                                 (vel_i[0]*theta +
                                                  un_i*UnitNormal[0]/3.0)*Area;
      // y-momentum
      val_Jacobian_ij[nSpecies+1][nSpecies]   += mu_i/(rho_i*dij) * etaz   * Area;
      val_Jacobian_ij[nSpecies+1][nSpecies+1] += mu_i/(rho_i*dij) * thetay * Area;
      val_Jacobian_ij[nSpecies+1][nSpecies+2] += mu_i/(rho_i*dij) *
                                                 (vel_i[1]*theta +
                                                  un_i*UnitNormal[1]/3.0)*Area;
    } else {
      // x-momentum
      val_Jacobian_ij[nSpecies][nSpecies]     += mu_i/(rho_i*dij) * thetax * Area;
      val_Jacobian_ij[nSpecies][nSpecies+1]   += mu_i/(rho_i*dij) * etaz   * Area;
      val_Jacobian_ij[nSpecies][nSpecies+2]   += mu_i/(rho_i*dij) * etay   * Area;
      val_Jacobian_ij[nSpecies][nSpecies+3]   += mu_i/(rho_i*dij) *
                                                 (vel_i[0]*theta +
                                                  un_i*UnitNormal[0]/3.0)*Area;
      // y-momentum
      val_Jacobian_ij[nSpecies+1][nSpecies]   += mu_i/(rho_i*dij) * etaz   * Area;
      val_Jacobian_ij[nSpecies+1][nSpecies+1] += mu_i/(rho_i*dij) * thetay * Area;
      val_Jacobian_ij[nSpecies+1][nSpecies+2] += mu_i/(rho_i*dij) * etax   * Area;
      val_Jacobian_ij[nSpecies+1][nSpecies+3] += mu_i/(rho_i*dij) *
                                                 (vel_i[1]*theta +
                                                  un_i*UnitNormal[1]/3.0)*Area;
      // z-momentum
      val_Jacobian_ij[nSpecies+2][nSpecies]   += mu_i/(rho_i*dij) * etay   * Area;
      val_Jacobian_ij[nSpecies+2][nSpecies+1] += mu_i/(rho_i*dij) * etax   * Area;
      val_Jacobian_ij[nSpecies+2][nSpecies+2] += mu_i/(rho_i*dij) * thetaz * Area;
      val_Jacobian_ij[nSpecies+2][nSpecies+3] += mu_i/(rho_i*dij) *
                                                 (vel_i[2]*theta +
                                                  un_i*UnitNormal[2]/3.0)*Area;
    }
    
    /*--- Jacobian from k = 3 viscous flux ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      val_Jacobian_ij[iVar][nSpecies+nDim] += ktr_i*dTdU_i[iVar]*(theta/dij)*Area;
    
    /*--- Jacobian from k = 4 viscous flux ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      val_Jacobian_ij[iVar][nSpecies+nDim]   += kve_i*dTvedU_i[iVar]*(theta/dij)*Area;
      val_Jacobian_ij[iVar][nSpecies+nDim+1] += kve_i*dTvedU_i[iVar]*(theta/dij)*Area;
    }
    

    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_ii[iVar][jVar] = -val_Jacobian_ij[iVar][jVar];
  }
  
  /*---+++ Residual at node j +++---*/
  // k = 1
  
  /*--- Calculate auxiliary quantities ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dYdrs[iSpecies][jSpecies] = 0.0;
      dJddrs[iSpecies][jSpecies] = 0.0;
    }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dYdrs[iSpecies][jSpecies] += 1.0/rho_j*(-V_j[RHOS_INDEX+iSpecies]/rho_j);
    }
    dYdrs[iSpecies][iSpecies] += 1.0/rho_j;
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    DdYk[iSpecies] = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      DdYk[iSpecies] += rho_j*D_j[jSpecies]*dYdrs[jSpecies][iSpecies];
    }
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Ys = V_j[RHOS_INDEX+iSpecies]/rho_j;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dJddrs[iSpecies][jSpecies] = -rho_j*D_j[iSpecies]*dYdrs[iSpecies][jSpecies]
                                 + Ys*DdYk[jSpecies];
    }
  }
  
  
  /*--- Calculate the k=1 residual at j ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      val_residual_j[iSpecies] += GPsirsdotn[jSpecies]*(-dJddrs[jSpecies][iSpecies])
                                + GPsiEdotn           *(-dJddrs[jSpecies][iSpecies])*hs_j[jSpecies]
                                + GPsiEvedotn         *(-dJddrs[jSpecies][iSpecies])*eve_j[jSpecies];
    }
  }
  
  
  //k = 2
  /*--- Calculate auxiliary quantities for SigmaPsiE ---*/
  GPsiEdotVel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    GPsiEdotVel += Mean_GradPsiE[iDim]*vel_j[iDim];
  
  /*--- Initialize SigmaPsiE ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++) {
      SigmaPsiE[iDim][jDim] = 0.0;
      SigmaVel[iDim][jDim] = 0.0;
    }
  
  /*--- Calculate SigmaPsiE & SigmaVel ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      SigmaPsiE[iDim][jDim] += Mean_GradPsiE[iDim]*vel_j[jDim] +
                               Mean_GradPsiE[jDim]*vel_j[iDim];
      SigmaVel[iDim][jDim] += vel_j[iDim]*UnitNormal[jDim]+vel_j[jDim]*UnitNormal[iDim];
    }
    SigmaPsiE[iDim][iDim] -= 2.0/3.0*GPsiEdotVel;
    SigmaVel[iDim][iDim]  -= 2.0/3.0*un_j;
  }
  
  /*--- Calculate the k=2 residual at j (SigmaPhi + SigmaPsiE) dot n ---*/
  SigVelGPhi = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      SigVelGPhi += SigmaVel[iDim][jDim]*Mean_GradPhi[jDim][iDim];
  
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    val_residual_j[iSpecies] += -mu_j/rho_j*(SigVelGPhi*Area +
                                             u2_j*GPsiEdotn  +
                                             un_j*GPsiEdotVel*Area/3.0);
  
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      val_residual_j[nSpecies+iDim] += mu_j/rho_j*(SigmaPhi[iDim][jDim] +
                                                   SigmaPsiE[iDim][jDim]  )
                                     * Normal[jDim];
    }
  }
  
  // k = 3
  /*--- Calculate the k=3 residual at j dT/dU * (GradPsiE dot n) ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual_j[iVar] += ktr_j*dTdU_j[iVar]*GPsiEdotn;
  
  
  // k = 4
  /*--- Calculate the k=4 residual at j ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual_j[iVar] += kve_j*dTvedU_j[iVar]*(GPsiEvedotn+GPsiEdotn);
  
  
  /*--- Calculate Jacobians for implicit time-stepping ---*/
  if (implicit) {
    
    /*--- Jacobian from k = 1 viscous flux ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        // species density
        val_Jacobian_jj[iSpecies][jSpecies] += -theta/dij*dJddrs[jSpecies][iSpecies]*Area;
        // total energy
        val_Jacobian_jj[iSpecies][nSpecies+nDim] += -theta/dij*dJddrs[jSpecies][iSpecies]*hs_j[jSpecies]*Area;
        // vib.-el. energy
        val_Jacobian_jj[iSpecies][nSpecies+nDim] += -theta/dij*dJddrs[jSpecies][iSpecies]*eve_j[jSpecies]*Area;
      }
    }
    
    /*--- Jacobian from k = 2 viscous flux ---*/
    if (nDim == 2) {
      // x-momentum
      val_Jacobian_jj[nSpecies][nSpecies]     += mu_j/(rho_j*dij) * thetax * Area;
      val_Jacobian_jj[nSpecies][nSpecies+1]   += mu_j/(rho_j*dij) * etaz   * Area;
      val_Jacobian_jj[nSpecies][nSpecies+2]   += mu_j/(rho_j*dij) *
                                                 (vel_j[0]*theta +
                                                  un_j*UnitNormal[0]/3.0)*Area;
      // y-momentum
      val_Jacobian_jj[nSpecies+1][nSpecies]   += mu_j/(rho_j*dij) * etaz   * Area;
      val_Jacobian_jj[nSpecies+1][nSpecies+1] += mu_j/(rho_j*dij) * thetay * Area;
      val_Jacobian_jj[nSpecies+1][nSpecies+2] += mu_j/(rho_j*dij) *
                                               (vel_j[1]*theta +
                                                un_j*UnitNormal[1]/3.0)*Area;
    } else {
      // x-momentum
      val_Jacobian_jj[nSpecies][nSpecies]     += mu_j/(rho_j*dij) * thetax * Area;
      val_Jacobian_jj[nSpecies][nSpecies+1]   += mu_j/(rho_j*dij) * etaz   * Area;
      val_Jacobian_jj[nSpecies][nSpecies+2]   += mu_j/(rho_j*dij) * etay   * Area;
      val_Jacobian_jj[nSpecies][nSpecies+3]   += mu_j/(rho_j*dij) *
                                                 (vel_j[0]*theta +
                                                 un_j*UnitNormal[0]/3.0)*Area;
      // y-momentum
      val_Jacobian_jj[nSpecies+1][nSpecies]   += mu_j/(rho_j*dij) * etaz   * Area;
      val_Jacobian_jj[nSpecies+1][nSpecies+1] += mu_j/(rho_j*dij) * thetay * Area;
      val_Jacobian_jj[nSpecies+1][nSpecies+2] += mu_j/(rho_j*dij) * etax   * Area;
      val_Jacobian_jj[nSpecies+1][nSpecies+3] += mu_j/(rho_j*dij) *
                                                 (vel_j[1]*theta +
                                                  un_j*UnitNormal[1]/3.0)*Area;
      // z-momentum
      val_Jacobian_jj[nSpecies+2][nSpecies]   += mu_j/(rho_j*dij) * etay   * Area;
      val_Jacobian_jj[nSpecies+2][nSpecies+1] += mu_j/(rho_j*dij) * etax   * Area;
      val_Jacobian_jj[nSpecies+2][nSpecies+2] += mu_j/(rho_j*dij) * thetaz * Area;
      val_Jacobian_jj[nSpecies+2][nSpecies+3] += mu_j/(rho_j*dij) *
                                               (vel_j[2]*theta +
                                                un_j*UnitNormal[2]/3.0)*Area;
    }
    
    /*--- Jacobian from k = 3 viscous flux ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      val_Jacobian_jj[iVar][nSpecies+nDim] += ktr_j*dTdU_j[iVar]*(theta/dij)*Area;
    
    /*--- Jacobian from k = 4 viscous flux ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      val_Jacobian_jj[iVar][nSpecies+nDim]   += kve_j*dTvedU_j[iVar]*(theta/dij)*Area;
      val_Jacobian_jj[iVar][nSpecies+nDim+1] += kve_j*dTvedU_j[iVar]*(theta/dij)*Area;
    }
    
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_ji[iVar][jVar] = -val_Jacobian_jj[iVar][jVar];
  }
}

CAvgGradCorrected_AdjTNE2::CAvgGradCorrected_AdjTNE2(unsigned short val_nDim,
                                                     unsigned short val_nVar,
                                                     CConfig *config) : CNumerics(val_nDim,
                                                                                  val_nVar,
                                                                                  config) {
	unsigned short iDim, iSpecies, iVar;
  
  implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  
  nDim         = val_nDim;
  nSpecies     = config->GetnSpecies();
  nVar         = val_nVar;
  
	vel   = new double[nDim];
  vel_i = new double[nDim];
  vel_j = new double[nDim];
  hs_i  = new double[nSpecies];
  hs_j  = new double[nSpecies];
  
  DdYk   = new double[nSpecies];
  dYdrs  = new double*[nSpecies];
  dJddrs = new double*[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    dYdrs[iSpecies]  = new double[nSpecies];
    dJddrs[iSpecies] = new double[nSpecies];
  }
  
  Proj_Mean_GPsi = new double[nVar];
  Mean_GPsi = new double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GPsi[iVar] = new double[nDim];
  
  GPsirsdotn = new double [nSpecies];
  Mean_GradPsirs = new double *[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Mean_GradPsirs[iSpecies] = new double[nDim];
  
	Mean_GradPhi = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Mean_GradPhi[iDim] = new double [nDim];
  
	Mean_GradPsiE = new double [nDim];
  Mean_GradPsiEve = new double [nDim];
	
  Edge_Vector = new double [nDim];
  
  SigmaPhi  = new double*[nDim];
  SigmaPsiE = new double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    SigmaPhi[iDim]  = new double[nDim];
    SigmaPsiE[iDim] = new double[nDim];
  }
  
}

CAvgGradCorrected_AdjTNE2::~CAvgGradCorrected_AdjTNE2(void) {
  unsigned short iDim, iSpecies, iVar;
  
  delete [] vel;
  delete [] vel_i;
  delete [] vel_j;
  delete [] hs_i;
  delete [] hs_j;
  
	delete [] Edge_Vector;
  
  delete [] DdYk;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    delete [] dJddrs[iSpecies];
    delete [] dYdrs[iSpecies];
  }
  delete [] dJddrs;
  delete [] dYdrs;
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GPsi[iVar];
  delete [] Mean_GPsi;
  delete [] Proj_Mean_GPsi;
  
  delete [] GPsirsdotn;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    delete [] Mean_GradPsirs[iSpecies];
  delete [] Mean_GradPsirs;
  
	delete [] Mean_GradPsiE;
  delete [] Mean_GradPsiEve;
  
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Mean_GradPhi[iDim];
  delete [] Mean_GradPhi;
  
  for (iDim = 0; iDim < nDim; iDim++) {
    delete [] SigmaPhi[iDim];
    delete [] SigmaPsiE[iDim];
  }
  delete [] SigmaPhi;
  delete [] SigmaPsiE;
  
}

void CAvgGradCorrected_AdjTNE2::ComputeResidual(double *val_residual_i,
                                                double *val_residual_j,
                                                double **val_Jacobian_ii,
                                                double **val_Jacobian_ij,
                                                double **val_Jacobian_ji,
                                                double **val_Jacobian_jj,
                                                CConfig *config) {
  
  
  unsigned short iDim, jDim, iSpecies, jSpecies, iVar, jVar;
  double mu_i, mu_j, ktr_i, ktr_j, kve_i, kve_j, *D_i, *D_j;
  double rho, rho_i, rho_j, un, Ys;
  double GdotPhi, GPsiEdotVel, GPsiEdotn, GPsiEvedotn;
  double dij, theta, thetax, thetay, thetaz, etax, etay, etaz;
  
  /*--- Initialize residuals ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual_i[iVar] = 0.0;
    val_residual_j[iVar] = 0.0;
    for (jVar = 0; jVar < nVar; jVar++) {
      val_Jacobian_ii[iVar][jVar] = 0.0;
      val_Jacobian_ij[iVar][jVar] = 0.0;
      val_Jacobian_ji[iVar][jVar] = 0.0;
      val_Jacobian_jj[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Calculate geometric quantities ---*/
  Area = 0.0;
  dij = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Area += Normal[iDim]*Normal[iDim];
    dij  += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
  }
  Area = sqrt(Area);
  dij  = sqrt(dij);
  theta = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    UnitNormal[iDim] = Normal[iDim]/Area;
    theta += UnitNormal[iDim]*UnitNormal[iDim];
  }
  if (nDim == 2) {
    thetax = theta + (UnitNormal[0]*UnitNormal[0])/3.0;
    thetay = theta + (UnitNormal[1]*UnitNormal[1])/3.0;
    etaz   = UnitNormal[0]*UnitNormal[1]/3.0;
  } else {
    thetax = theta + (UnitNormal[0]*UnitNormal[0])/3.0;
    thetay = theta + (UnitNormal[1]*UnitNormal[1])/3.0;
    thetaz = theta + (UnitNormal[2]*UnitNormal[2])/3.0;
    etax   = UnitNormal[1]*UnitNormal[2]/3.0;
    etay   = UnitNormal[0]*UnitNormal[2]/3.0;
    etaz   = UnitNormal[0]*UnitNormal[1]/3.0;
  }
  
  /*--- Get flow state (Rename for convenience) ---*/
  D_i   = Diffusion_Coeff_i;
  D_j   = Diffusion_Coeff_j;
  mu_i  = Laminar_Viscosity_i;
  mu_j  = Laminar_Viscosity_j;
  ktr_i = Thermal_Conductivity_i;
  ktr_j = Thermal_Conductivity_j;
  kve_i = Thermal_Conductivity_ve_i;
  kve_j = Thermal_Conductivity_ve_j;
  rho_i = V_i[RHO_INDEX];
  rho_j = V_j[RHO_INDEX];
  rho   = 0.5*(rho_i+rho_j);
  for (iDim = 0; iDim < nDim; iDim++) {
    vel_i[iDim] = V_i[VEL_INDEX+iDim];
    vel_j[iDim] = V_j[VEL_INDEX+iDim];
    vel[iDim] = 0.5*(vel_i[iDim]+vel_j[iDim]);
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    hs_i[iSpecies]  = var->CalcHs (config, V_i[T_INDEX], eve_i[iSpecies], iSpecies);
    hs_j[iSpecies]  = var->CalcHs (config, V_j[T_INDEX], eve_j[iSpecies], iSpecies);
  }
  
  /*--- Compute vector going from iPoint to jPoint ---*/
	dij = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dij += Edge_Vector[iDim]*Edge_Vector[iDim];
	}
  dij = sqrt(dij);
  
	/*--- Mean gradient approximation ---*/
  // Note: Projection of the mean gradient in the direction of the edge, weiss correction
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Mean_GPsi[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GPsi[iVar][iDim] = 0.5*(PsiVar_Grad_i[iVar][iDim] +
                                   PsiVar_Grad_j[iVar][iDim]  );
			Proj_Mean_GPsi[iVar] += Mean_GPsi[iVar][iDim]*Edge_Vector[iDim];
		}
		for (iDim = 0; iDim < nDim; iDim++)
			Mean_GPsi[iVar][iDim] -= (Proj_Mean_GPsi[iVar]     -
                                (Psi_j[iVar]-Psi_i[iVar]) )
                             * Edge_Vector[iDim]/(dij*dij);
	}
  
  /*--- Assign mean gradients ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      Mean_GradPsirs[iSpecies][iDim] = Mean_GPsi[iSpecies][iDim];
    Mean_GradPsiE[iDim]              = Mean_GPsi[nSpecies+nDim][iDim];
    Mean_GradPsiEve[iDim]            = Mean_GPsi[nSpecies+nDim+1][iDim];
		for (jDim = 0; jDim < nDim; jDim++)
      Mean_GradPhi[iDim][jDim]       =  Mean_GPsi[nSpecies+iDim][jDim];
  }
  
  /*--- Calculate auxiliary quantities for SigmaPhi ---*/
  GdotPhi = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    GdotPhi += Mean_GradPhi[iDim][iDim];
  
  /*--- Project mean gradient of PsiRs, PsiE & PsiEve into normal ---*/
  GPsiEdotn   = 0.0;
  GPsiEvedotn = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    GPsirsdotn[iSpecies] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      GPsirsdotn[iSpecies] +=Mean_GradPsirs[iSpecies][iDim]*Normal[iDim];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    GPsiEdotn   += Mean_GradPsiE[iDim]*Normal[iDim];
    GPsiEvedotn += Mean_GradPsiEve[iDim]*Normal[iDim];
  }
  
  /*--- Initialize SigmaPhi ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      SigmaPhi[iDim][jDim] = 0.0;
  
  /*--- Calculate SigmaPhi ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      SigmaPhi[iDim][jDim] += Mean_GradPhi[iDim][jDim] +
                              Mean_GradPhi[jDim][iDim];
    }
    SigmaPhi[iDim][iDim]  -= 2.0/3.0*GdotPhi;
  }
  
  
  /*---+++ Residual at node i +++---*/
  
  // k = 1
  /*--- Calculate auxiliary quantities ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dYdrs[iSpecies][jSpecies]  = 0.0;
      dJddrs[iSpecies][jSpecies] = 0.0;
    }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dYdrs[iSpecies][jSpecies] += 1.0/rho_i*(-V_i[RHOS_INDEX+iSpecies]/rho_i);
    }
    dYdrs[iSpecies][iSpecies] += 1.0/rho_i;
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    DdYk[iSpecies] = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      DdYk[iSpecies] += rho_i*D_i[jSpecies]*dYdrs[jSpecies][iSpecies];
    }
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Ys = V_i[RHOS_INDEX+iSpecies]/rho_i;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dJddrs[iSpecies][jSpecies] = -rho_i*D_i[iSpecies]*dYdrs[iSpecies][jSpecies]
                                 + Ys*DdYk[jSpecies];
    }
  }
  
  /*--- Calculate the k=1 residual at i ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      val_residual_i[iSpecies] +=
          GPsirsdotn[jSpecies]*(-dJddrs[jSpecies][iSpecies])
        + GPsiEdotn           *(-dJddrs[jSpecies][iSpecies])*hs_i[jSpecies]
        + GPsiEvedotn         *(-dJddrs[jSpecies][iSpecies])*eve_i[jSpecies];
    }
  }
  
  // k = 2
  /*--- Calculate auxiliary quantities for SigmaPsiE ---*/
  GPsiEdotVel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    GPsiEdotVel += Mean_GradPsiE[iDim]*vel_i[iDim];
  
  /*--- Initialize SigmaPsiE ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      SigmaPsiE[iDim][jDim] = 0.0;
  
  /*--- Calculate SigmaPsiE ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      SigmaPsiE[iDim][jDim] += Mean_GradPsiE[iDim]*vel_i[jDim] +
                               Mean_GradPsiE[jDim]*vel_i[iDim];
    }
    SigmaPsiE[iDim][iDim] -= 2.0/3.0*GPsiEdotVel;
  }
  
  /*--- Calculate the k=2 residual at i (SigmaPhi + SigmaPsiE) dot n ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      val_residual_i[nSpecies+iDim] += mu_i/rho_i*(SigmaPhi[iDim][jDim] +
                                                   SigmaPsiE[iDim][jDim]  )
                                     * Normal[jDim];
    }
  }
  
  // k = 3
  /*--- Calculate the k=3 residual at i dT/dU * (GradPsiE dot n) ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual_i[iVar] += ktr_i*dTdU_i[iVar]*GPsiEdotn;
  
  // k = 4
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual_i[iVar] += kve_i*dTvedU_i[iVar]*(GPsiEvedotn+GPsiEdotn);
  
  /*--- Calculate Jacobians for implicit time-stepping ---*/
  if (implicit) {
    
    /*--- Calculate projected velocity at node i ---*/
    un = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      un += vel_i[iDim]*UnitNormal[iDim];
    
    /*--- Jacobian from k = 1 viscous flux ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        // species density
        val_Jacobian_ij[iSpecies][jSpecies] += -theta/dij*dJddrs[jSpecies][iSpecies]*Area;
        // total energy
        val_Jacobian_ij[iSpecies][nSpecies+nDim] += -theta/dij*dJddrs[jSpecies][iSpecies]*hs_i[jSpecies]*Area;
        // vib.-el. energy
        val_Jacobian_ij[iSpecies][nSpecies+nDim] += -theta/dij*dJddrs[jSpecies][iSpecies]*eve_i[jSpecies]*Area;
      }
    }
    
    /*--- Jacobian from k = 2 viscous flux ---*/
    
    if (nDim == 2) {
      // x-momentum
      val_Jacobian_ij[nSpecies][nSpecies]     += mu_i/(rho_i*dij) * thetax * Area;
      val_Jacobian_ij[nSpecies][nSpecies+1]   += mu_i/(rho_i*dij) * etaz   * Area;
      val_Jacobian_ij[nSpecies][nSpecies+2]   += mu_i/(rho_i*dij) *
                                                 (vel_i[0]*theta +
                                                  un*UnitNormal[0]/3.0)*Area;
      // y-momentum
      val_Jacobian_ij[nSpecies+1][nSpecies]   += mu_i/(rho_i*dij) * etaz   * Area;
      val_Jacobian_ij[nSpecies+1][nSpecies+1] += mu_i/(rho_i*dij) * thetay * Area;
      val_Jacobian_ij[nSpecies+1][nSpecies+2] += mu_i/(rho_i*dij) *
                                                 (vel_i[1]*theta +
                                                  un*UnitNormal[1]/3.0)*Area;
      
    } else {
      // x-momentum
      val_Jacobian_ij[nSpecies][nSpecies]     += mu_i/(rho_i*dij) * thetax * Area;
      val_Jacobian_ij[nSpecies][nSpecies+1]   += mu_i/(rho_i*dij) * etaz   * Area;
      val_Jacobian_ij[nSpecies][nSpecies+2]   += mu_i/(rho_i*dij) * etay   * Area;
      val_Jacobian_ij[nSpecies][nSpecies+3]   += mu_i/(rho_i*dij) *
                                                 (vel_i[0]*theta +
                                                  un*UnitNormal[0]/3.0)*Area;
      // y-momentum
      val_Jacobian_ij[nSpecies+1][nSpecies]   += mu_i/(rho_i*dij) * etaz   * Area;
      val_Jacobian_ij[nSpecies+1][nSpecies+1] += mu_i/(rho_i*dij) * thetay * Area;
      val_Jacobian_ij[nSpecies+1][nSpecies+2] += mu_i/(rho_i*dij) * etax   * Area;
      val_Jacobian_ij[nSpecies+1][nSpecies+3] += mu_i/(rho_i*dij) *
                                                 (vel_i[1]*theta +
                                                  un*UnitNormal[1]/3.0)*Area;
      // z-momentum
      val_Jacobian_ij[nSpecies+2][nSpecies]   += mu_i/(rho_i*dij) * etay   * Area;
      val_Jacobian_ij[nSpecies+2][nSpecies+1] += mu_i/(rho_i*dij) * etax   * Area;
      val_Jacobian_ij[nSpecies+2][nSpecies+2] += mu_i/(rho_i*dij) * thetaz * Area;
      val_Jacobian_ij[nSpecies+2][nSpecies+3] += mu_i/(rho_i*dij) *
                                                 (vel_i[2]*theta +
                                                  un*UnitNormal[2]/3.0)*Area;
    }
    
    /*--- Jacobian from k = 3 viscous flux ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      val_Jacobian_ij[iVar][nSpecies+nDim] += ktr_i*(theta/dij)*dTdU_i[iVar]*Area;
    
    /*--- Jacobian from k = 4 viscous flux ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      val_Jacobian_ij[iVar][nSpecies+nDim]   += kve_i*(theta/dij)*dTvedU_i[iVar]*Area;
      val_Jacobian_ij[iVar][nSpecies+nDim+1] += kve_i*(theta/dij)*dTvedU_i[iVar]*Area;
    }
    
    
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_ii[iVar][jVar] = -val_Jacobian_ij[iVar][jVar];
  }
  
  /*---+++ Residual at node j +++---*/
  // k = 1
  
  /*--- Calculate auxiliary quantities ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dYdrs[iSpecies][jSpecies] = 0.0;
      dJddrs[iSpecies][jSpecies] = 0.0;
    }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dYdrs[iSpecies][jSpecies] += 1.0/rho_j*(-V_j[RHOS_INDEX+iSpecies]/rho_j);
    }
    dYdrs[iSpecies][iSpecies] += 1.0/rho_j;
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    DdYk[iSpecies] = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      DdYk[iSpecies] += rho_j*D_j[jSpecies]*dYdrs[jSpecies][iSpecies];
    }
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Ys = V_j[RHOS_INDEX+iSpecies]/rho_j;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dJddrs[iSpecies][jSpecies] = -rho_j*D_j[iSpecies]*dYdrs[iSpecies][jSpecies]
      + Ys*DdYk[jSpecies];
    }
  }
  
  
  /*--- Calculate the k=1 residual at j ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      val_residual_j[iSpecies] += GPsirsdotn[jSpecies]*(-dJddrs[jSpecies][iSpecies])
      + GPsiEdotn           *(-dJddrs[jSpecies][iSpecies])*hs_j[jSpecies]
      + GPsiEvedotn         *(-dJddrs[jSpecies][iSpecies])*eve_j[jSpecies];
    }
  }
  
  
  //k = 2
  /*--- Calculate auxiliary quantities for SigmaPsiE ---*/
  GPsiEdotVel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    GPsiEdotVel += Mean_GradPsiE[iDim]*vel_j[iDim];
  
  /*--- Initialize SigmaPsiE ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      SigmaPsiE[iDim][jDim] = 0.0;
  
  /*--- Calculate SigmaPsiE ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      SigmaPsiE[iDim][jDim] += Mean_GradPsiE[iDim]*vel_j[jDim] +
      Mean_GradPsiE[jDim]*vel_j[iDim];
    }
    SigmaPsiE[iDim][iDim] -= 2.0/3.0*GPsiEdotVel;
  }
  
  /*--- Calculate the residual at j (SigmaPhi + SigmaPsiE) dot n ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      val_residual_j[nSpecies+iDim] += mu_j/rho_j*(SigmaPhi[iDim][jDim] +
                                                   SigmaPsiE[iDim][jDim]  )
                                     * Normal[jDim];
    }
  }
  
  // k = 3
  /*--- Calculate the k=3 residual at j dT/dU * (GradPsiE dot n) ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual_j[iVar] += ktr_j*dTdU_j[iVar]*GPsiEdotn;
  
  
  // k = 4
  /*--- Calculate the k=4 residual at j ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual_j[iVar] += kve_j*dTvedU_j[iVar]*(GPsiEvedotn+GPsiEdotn);
  
  
  /*--- Calculate Jacobians for implicit time-stepping ---*/
  if (implicit) {
    
    /*--- Calculate projected velocity at node i ---*/
    un = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      un += vel_j[iDim]*UnitNormal[iDim];
    
    /*--- Jacobian from k = 1 viscous flux ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        // species density
        val_Jacobian_jj[iSpecies][jSpecies] += -theta/dij*dJddrs[jSpecies][iSpecies]*Area;
        // total energy
        val_Jacobian_jj[iSpecies][nSpecies+nDim] += -theta/dij*dJddrs[jSpecies][iSpecies]*hs_j[jSpecies]*Area;
        // vib.-el. energy
        val_Jacobian_jj[iSpecies][nSpecies+nDim] += -theta/dij*dJddrs[jSpecies][iSpecies]*eve_j[jSpecies]*Area;
      }
    }
    
    
    /*--- Jacobian from k = 2 viscous flux ---*/
    if (nDim == 2) {
      // x-momentum
      val_Jacobian_jj[nSpecies][nSpecies]     += mu_j/(rho_j*dij) * thetax * Area;
      val_Jacobian_jj[nSpecies][nSpecies+1]   += mu_j/(rho_j*dij) * etaz   * Area;
      val_Jacobian_jj[nSpecies][nSpecies+2]   += mu_j/(rho_j*dij) *
                                                 (vel_j[0]*theta +
                                                  un*UnitNormal[0]/3.0)*Area;
      // y-momentum
      val_Jacobian_jj[nSpecies+1][nSpecies]   += mu_j/(rho_j*dij) * etaz   * Area;
      val_Jacobian_jj[nSpecies+1][nSpecies+1] += mu_j/(rho_j*dij) * thetay * Area;
      val_Jacobian_jj[nSpecies+1][nSpecies+2] += mu_j/(rho_j*dij) *
                                                 (vel_j[1]*theta +
                                                  un*UnitNormal[1]/3.0)*Area;
    } else {
      // x-momentum
      val_Jacobian_jj[nSpecies][nSpecies]     += mu_j/(rho_j*dij) * thetax * Area;
      val_Jacobian_jj[nSpecies][nSpecies+1]   += mu_j/(rho_j*dij) * etaz   * Area;
      val_Jacobian_jj[nSpecies][nSpecies+2]   += mu_j/(rho_j*dij) * etay   * Area;
      val_Jacobian_jj[nSpecies][nSpecies+3]   += mu_j/(rho_j*dij) *
                                                 (vel_j[0]*theta +
                                                  un*UnitNormal[0]/3.0)*Area;
      // y-momentum
      val_Jacobian_jj[nSpecies+1][nSpecies]   += mu_j/(rho_j*dij) * etaz   * Area;
      val_Jacobian_jj[nSpecies+1][nSpecies+1] += mu_j/(rho_j*dij) * thetay * Area;
      val_Jacobian_jj[nSpecies+1][nSpecies+2] += mu_j/(rho_j*dij) * etax   * Area;
      val_Jacobian_jj[nSpecies+1][nSpecies+3] += mu_j/(rho_j*dij) *
                                                 (vel_j[1]*theta +
                                                  un*UnitNormal[1]/3.0)*Area;
      // z-momentum
      val_Jacobian_jj[nSpecies+2][nSpecies]   += mu_j/(rho_j*dij) * etay   * Area;
      val_Jacobian_jj[nSpecies+2][nSpecies+1] += mu_j/(rho_j*dij) * etax   * Area;
      val_Jacobian_jj[nSpecies+2][nSpecies+2] += mu_j/(rho_j*dij) * thetaz * Area;
      val_Jacobian_jj[nSpecies+2][nSpecies+3] += mu_j/(rho_j*dij) *
                                                 (vel_j[2]*theta +
                                                  un*UnitNormal[2]/3.0)*Area;
    }
    
    /*--- Jacobian from k = 3 viscous flux ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      val_Jacobian_jj[iVar][nSpecies+nDim] += ktr_j*dTdU_j[iVar]*(theta/dij)*Area;
    
    /*--- Jacobian from k = 4 viscous flux ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      val_Jacobian_jj[iVar][nSpecies+nDim]   += kve_j*dTvedU_j[iVar]*(theta/dij)*Area;
      val_Jacobian_jj[iVar][nSpecies+nDim+1] += kve_j*dTvedU_j[iVar]*(theta/dij)*Area;
    }
    
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_ji[iVar][jVar] = -val_Jacobian_jj[iVar][jVar];
  }
}


CSource_AdjTNE2::CSource_AdjTNE2(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 unsigned short val_nPrimVar,
                                 unsigned short val_nPrimVarGrad,
                                 CConfig *config) : CNumerics(val_nDim,
                                                              val_nVar,
                                                              config) {
  
  unsigned short iDim, iSpecies, jSpecies, iVar;
  
  /*--- Assign booleans from CConfig ---*/
  implicit   = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
  
  /*--- Set array sizes ---*/
  nDim         = val_nDim;
  nSpecies     = config->GetnSpecies();
  nVar         = val_nVar;
  nPrimVar     = val_nPrimVar;
  nPrimVarGrad = val_nPrimVarGrad;
  
  rhos  = new double[nSpecies];
  vel   = new double[nDim];
  V     = new double[nPrimVar];
  GV = new double*[nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    GV[iVar] = new double[nDim];
  
  
  GInvRho      = new double[nDim];
  GVeloRho     = new double*[nDim];
  GPhiGInvRho  = new double[nDim];
  GPsiEZetaTau = new double[nDim];
  tau          = new double*[nDim];
  eta          = new double*[nDim];
  pi           = new double*[nDim];
  zeta         = new double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    GVeloRho[iDim] = new double[nDim];
    tau[iDim]      = new double[nDim];
    eta[iDim]      = new double[nDim];
    pi[iDim]       = new double[nDim];
    zeta[iDim]     = new double[nDim];
  }

  Y     = new double[nSpecies];
  hs    = new double[nSpecies];
  Cvtrs = new double[nSpecies];
  SIk   = new double[nDim];
  GY    = new double*[nSpecies];
  SdIdr = new double*[nSpecies];
  Js    = new double*[nSpecies];
  dIdr  = new double**[nSpecies];
  dJdr  = new double**[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    GY[iSpecies]    = new double[nDim];
    SdIdr[iSpecies] = new double[nDim];
    Js[iSpecies]    = new double[nDim];
    dIdr[iSpecies]  = new double*[nSpecies];
    dJdr[iSpecies]  = new double*[nSpecies];
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      dIdr[iSpecies][jSpecies] = new double[nDim];
      dJdr[iSpecies][jSpecies] = new double[nDim];
    }
  }
}

CSource_AdjTNE2::~CSource_AdjTNE2(void) {

  unsigned short iDim, iSpecies, jSpecies, iVar;
  
  delete [] rhos;
  delete [] vel;
  delete [] Y;
  delete [] hs;
  delete [] Cvtrs;
  delete [] V;
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    delete [] GV[iVar];
  delete [] GV;
  

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
  
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      delete [] dIdr[iSpecies][jSpecies];
      delete [] dJdr[iSpecies][jSpecies];
    }
    delete [] dIdr[iSpecies];
    delete [] dJdr[iSpecies];
    delete [] GY[iSpecies];
    delete [] Js[iSpecies];
    delete [] SdIdr[iSpecies];
  }
  delete [] GY;
  delete [] SIk;
  delete [] dIdr;
  delete [] dJdr;
  delete [] SdIdr;
  
}

void CSource_AdjTNE2::ComputeSourceViscous (double *val_residual, CConfig *config) {
  
  /*--- This subroutine computes the viscous source term, GPsi \cdot Av ---*/
  // Note: The viscous fluxes are split into four different flux vectors,
  // one for each component (diffusion, momentum transport, thermal transport,
  // and vib-el. thermal transport).
  
  unsigned short iDim, jDim, iSpecies, iVar;
  double rho, Tve, rCvtr, rCvve;
  double mu, ktr, kve;
  double *Ms, *xi, Ru;
  double div_vel, div_velorho, velGInvRho;
  double **GPsi, GPsiEGvel[3], GPhiEta, GPsiEvelPi;
  
  
  /*--- Initialize residual vector ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = 0.0;
  
  /*--- Rename variables for convenience ---*/
  rho   = V_i[RHO_INDEX];
  Tve   = V_i[TVE_INDEX];
  rCvtr = V_i[RHOCVTR_INDEX];
  rCvve = V_i[RHOCVVE_INDEX];
  GPsi  = PsiVar_Grad_i;
  Ru    = UNIVERSAL_GAS_CONSTANT;
  Ms    = config->GetMolar_Mass();
  xi    = config->GetRotationModes();
  mu    = Laminar_Viscosity_i;
  ktr   = Thermal_Conductivity_i;
  kve   = Thermal_Conductivity_ve_i;
  div_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    vel[iDim] = V_i[VEL_INDEX+iDim];
    div_vel  += PrimVar_Grad_i[VEL_INDEX+iDim][iDim];
  }
  
  /*--- Calculate specific heat ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Cvtrs[iSpecies] = (3.0+xi[iSpecies])/2.0*Ru/Ms[iSpecies];
  }
  
  /*--- Convert from species density to mass-fraction, Y ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    V[RHOS_INDEX+iSpecies] = V_i[RHOS_INDEX+iSpecies]/rho;
    for (iDim = 0; iDim < nDim; iDim++)
      GV[RHOS_INDEX+iSpecies][iDim] = 1.0/rho*(PrimVar_Grad_i[RHOS_INDEX+iSpecies][iDim] -
                                               V[RHOS_INDEX+iSpecies]*
                                               PrimVar_Grad_i[RHO_INDEX][iDim]);
  }
  
  /*--- Copy remaining quantities to V & GV ---*/
  for (iVar = nSpecies; iVar < nPrimVar; iVar++)
    V[iVar] = V_i[iVar];
  for (iVar = nSpecies; iVar < nPrimVarGrad; iVar++)
    for (iDim = 0; iDim < nDim; iDim++)
      GV[iVar][iDim] = PrimVar_Grad_i[iVar][iDim];
  
  /*--- Calculate supporting quantities ---*/
  // GPsiEGvel = GPsiE * G{u,v,w} = {GPsiEdotGu, GPsiEdotGv, GPsiEdotGw}
  for (iDim = 0; iDim < nDim; iDim++) {
    GPsiEGvel[iDim] = 0.0;
    for (jDim = 0; jDim < nDim; jDim++) {
      GPsiEGvel[iDim] += GPsi[nSpecies+nDim][jDim]*GV[VEL_INDEX+iDim][jDim];
    }
  }
  // G(1/rho) and G(u/rho)
  for (iDim = 0; iDim < nDim; iDim++) {
    GInvRho[iDim] = -1.0/(rho*rho)*PrimVar_Grad_i[RHO_INDEX][iDim];
    for (jDim = 0; jDim < nDim; jDim++)
      GVeloRho[iDim][jDim] = -vel[iDim]/(rho*rho)*(PrimVar_Grad_i[RHO_INDEX][jDim])
                           +  1.0/rho*(PrimVar_Grad_i[VEL_INDEX+iDim][jDim]);
  }
  div_velorho = 0.0;
  velGInvRho  = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    div_velorho += GVeloRho[iDim][iDim];
    velGInvRho  += vel[iDim]*GInvRho[iDim];
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
  
  /*--- Contribution ot viscous source from Av2 (momentum transport) ---*/
  // Fv2 = [0, ..., 0, tx, ty, tz, utx+vty+wtz, 0]
  // Mass conservation
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    val_residual[iSpecies] += (-GPhiEta + GPsiEvelPi)*mu*Volume;
  }
  if (nDim == 2) {
    // X-momentum
    val_residual[nSpecies]   += ((GPhiGInvRho[0] + 1.0/3.0*GPsi[nSpecies][0]*GInvRho[0]) +
                                 (GPsi[nSpecies+1][0]*GInvRho[1] - 2.0/3.0*GPsi[nSpecies+1][1]*GInvRho[0]) +
                                 (GPsi[nSpecies+nDim][0]*velGInvRho + GPsiEZetaTau[0]))*mu*Volume;
    // Y-momentum
    val_residual[nSpecies+1] += ((GPsi[nSpecies][1]*GInvRho[0]   - 2.0/3.0*GPsi[nSpecies][0]*GInvRho[1]) +
                                 (GPhiGInvRho[1] + 1.0/3.0*GPsi[nSpecies+1][1]*GInvRho[1]) +
                                 (GPsi[nSpecies+nDim][1]*velGInvRho + GPsiEZetaTau[1]))*mu*Volume;
  } else {
    // X-momentum
    val_residual[nSpecies]   += ((GPhiGInvRho[0] + 1.0/3.0*GPsi[nSpecies][0]*GInvRho[0]) +
                                 (GPsi[nSpecies+1][0]*GInvRho[1] - 2.0/3.0*GPsi[nSpecies+1][1]*GInvRho[0]) +
                                 (GPsi[nSpecies+2][0]*GInvRho[2] - 2.0/3.0*GPsi[nSpecies+2][2]*GInvRho[0]) +
                                 (GPsi[nSpecies+3][0]*velGInvRho + GPsiEZetaTau[0]))*mu*Volume;
    // Y-momentum
    val_residual[nSpecies+1] += ((GPsi[nSpecies][1]*GInvRho[0]   - 2.0/3.0*GPsi[nSpecies][0]*GInvRho[1]) +
                                 (GPhiGInvRho[1] + 1.0/3.0*GPsi[nSpecies+1][1]*GInvRho[1]) +
                                 (GPsi[nSpecies+2][1]*GInvRho[2] - 2.0/3.0*GPsi[nSpecies+2][2]*GInvRho[1]) +
                                 (GPsi[nSpecies+3][1]*velGInvRho + GPsiEZetaTau[1]))*mu*Volume;
    // Z-momentum
    val_residual[nSpecies+2] += ((GPsi[nSpecies][2]*GInvRho[0]   - 2.0/3.0*GPsi[nSpecies][0]*GInvRho[2]) +
                                 (GPsi[nSpecies+1][2]*GInvRho[1] - 2.0/3.0*GPsi[nSpecies+1][1]*GInvRho[2]) +
                                 (GPhiGInvRho[2] + 1.0/3.0*GPsi[nSpecies+2][2]*GInvRho[2]) +
                                 (GPsi[nSpecies+3][2]*velGInvRho + GPsiEZetaTau[2]))*mu*Volume;
  }
  
  
  /*--- Contribution to viscous source from Av3 (thermal transport) ---*/
  // Fv3 = [0, ..., 0, 0, 0, 0, GT, 0]
  // Av3 = [0, ..., 0, 0, 0, 0, dGT/dU, 0]
  // By fundamental theorem of calculus, d/dU(GT) = G(dT/dU)
  // Residual contribution is GPsi^T * ktr*Av3 * Volume
  
  // Species continuity
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    for (iDim = 0; iDim < nDim; iDim++) {
      val_residual[iSpecies] += GPsiEGvel[iDim]*V[VEL_INDEX+iDim]/rCvtr*ktr*Volume;
      val_residual[iSpecies] += -(Cvtrs[iSpecies]/rCvtr)*GPsi[nSpecies+nDim][iDim]*GV[T_INDEX][iDim]*ktr*Volume;
      val_residual[iSpecies] += -(dTdU_i[iSpecies]/rCvtr)*GPsi[nSpecies+nDim][iDim]*GV[RHOCVTR_INDEX][iDim]*ktr*Volume;
    }
  // Momentum
  for (iDim = 0; iDim < nDim; iDim++) {
    val_residual[nSpecies+iDim] += -1.0/rCvtr*GPsiEGvel[iDim]*ktr*Volume;
    for (jDim = 0; jDim < nDim; jDim++) {
      val_residual[nSpecies+iDim] += -dTdU_i[nSpecies+iDim]/rCvtr*GPsi[nSpecies+nDim][jDim]*GV[RHOCVTR_INDEX][jDim]*ktr*Volume;
    }
  }
  // Energy
  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[nSpecies+nDim] += -dTdU_i[nSpecies+nDim]/rCvtr*GPsi[nSpecies+nDim][iDim]*GV[RHOCVTR_INDEX][iDim]*ktr*Volume;
  // Vib - el. Energy
  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[nSpecies+nDim+1] += -dTdU_i[nSpecies+nDim+1]/rCvtr*GPsi[nSpecies+nDim][iDim]*GV[RHOCVTR_INDEX][iDim]*ktr*Volume;
  
  
  /*--- Contribution to viscous source term from Av4 (vib.-el. transport) ---*/
  // Fv4 = [0, ..., 0, 0, 0, 0, GTve, GTve]
  // Av4 = [0, ..., 0, 0, 0, 0, dGTve/dU, dGTve/dU]
  // By fundamental theorem of calculus, d/dU(GTve) = G(dTve/dU)
  // Residual contribution is GPsi^T * kve*Av4 * Volume
  
  // Species continuity
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      val_residual[iSpecies] += -Cvve_i[iSpecies]/rCvve*(GPsi[nSpecies+nDim][iDim]+
                                                         GPsi[nSpecies+nDim+1][iDim])*GV[TVE_INDEX][iDim]*kve*Volume;
      val_residual[iSpecies] += -dTvedU_i[iSpecies]/rCvve*(GPsi[nSpecies+nDim][iDim]+
                                                           GPsi[nSpecies+nDim+1][iDim])*GV[RHOCVVE_INDEX][iDim]*kve*Volume;
    }
  }
  // Vib-el energy
  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[nSpecies+nDim+1] += -dTvedU_i[nSpecies+nDim+1]/rCvve*(GPsi[nSpecies+nDim][iDim]+
                                                                       GPsi[nSpecies+nDim+1][iDim])*GV[RHOCVVE_INDEX][iDim]*kve*Volume;
  
  
  
//	unsigned short iDim, jDim, iSpecies, jSpecies, iVar, jVar;
//  double rho, sqvel, rhoCvtr, rhoCvve, YDbar, *Ms, *xi, Ru;
//  double div_vel, div_velorho, velGInvRho, GPhiEta, GPsiEvelPi;
//  double GPsiEGr, GPsiEveGr, DGPsiEGrs, DGPsiEveGrs;
//  double *Ds, mu2, mu3, mu4;
//  double **GPsi, **GV, **GU;
//  double *PrimVar;
//  
//  PrimVar = new double[nPrimVar];
//  
//  /*--- Initialize arrays ---*/
//  for (iVar = 0; iVar < nVar; iVar++) {
//    val_residual[iVar] = 0.0;
//  }
//  
//  /*--- Rename for convenience ---*/
//  Ms      = config->GetMolar_Mass();
//  xi      = config->GetRotationModes();
//  Ru      = UNIVERSAL_GAS_CONSTANT;
//  rhoCvtr = V_i[RHOCVTR_INDEX];
//  rhoCvve = V_i[RHOCVVE_INDEX];
//  GPsi    = PsiVar_Grad_i;
//  GV      = PrimVar_Grad_i;
//  GU      = ConsVar_Grad_i;
//  
//  /*--- Assign useful values ---*/
//  // Mixture & species density, mass fraction
//  rho = V_i[RHO_INDEX];
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    rhos[iSpecies] = V_i[RHOS_INDEX+iSpecies];
//    Y[iSpecies]    = rhos[iSpecies]/rho;
//  }
//  // Velocity
//  sqvel = 0.0;
//  div_vel = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++) {
//		vel[iDim] = V_i[VEL_INDEX+iDim];
//		sqvel += vel[iDim]*vel[iDim];
//    div_vel += PrimVar_Grad_i[VEL_INDEX+iDim][iDim];
//	}
//  // Species enthalpy and vib.-el. energy
//  for (iSpecies     = 0; iSpecies < nSpecies; iSpecies++) {
//    eves[iSpecies]  = var->CalcEve(config, V_i[TVE_INDEX], iSpecies);
//    hs[iSpecies]    = var->CalcHs(config, V_i[T_INDEX], eves[iSpecies], iSpecies);
//    Cvves[iSpecies] = var->CalcCvve(V_i[TVE_INDEX], config, iSpecies);
//    Cvtrs[iSpecies] = (3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
//  }
//  
//  // Transport coefficients
//  mu2 = Laminar_Viscosity_i;
//  mu3 = Thermal_Conductivity_i;
//  mu4 = Thermal_Conductivity_ve_i;
//  Ds  = Diffusion_Coeff_i;
//
//  /*--- Gradients of flow variables ---*/
//  for (iDim = 0; iDim < nDim; iDim++) {
//    GInvRho[iDim] = -1.0/(rho*rho)*PrimVar_Grad_i[RHO_INDEX][iDim];
//    for (jDim = 0; jDim < nDim; jDim++)
//      GVeloRho[iDim][jDim] = -vel[iDim]/(rho*rho)*(PrimVar_Grad_i[RHO_INDEX][jDim])
//                           +  1.0/rho*(PrimVar_Grad_i[VEL_INDEX+iDim][jDim]);
//  }
//  div_velorho = 0.0;
//  velGInvRho = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++) {
//    div_velorho += GVeloRho[iDim][iDim];
//    velGInvRho += vel[iDim]*GInvRho[iDim];
//  }
//  
//  /*--- Supporting matrices ---*/
//  // Tau: Stress tensor matrix
//  for (iDim = 0; iDim < nDim; iDim++) {
//    for (jDim = 0; jDim < nDim; jDim++) {
//      tau[iDim][jDim] = PrimVar_Grad_i[VEL_INDEX+iDim][jDim] +
//                        PrimVar_Grad_i[VEL_INDEX+jDim][iDim];
//    }
//    tau[iDim][iDim] -= 2.0/3.0*div_vel;
//  }
//  // Eta
//  for (iDim = 0; iDim < nDim; iDim++) {
//    for (jDim = 0; jDim < nDim; jDim++) {
//      eta[iDim][jDim] = GVeloRho[iDim][jDim] + GVeloRho[jDim][iDim];
//    }
//    eta[iDim][iDim] -= 2.0/3.0*div_velorho;
//  }
//  // Pi
//  for (iDim = 0; iDim < nDim; iDim++) {
//    for (jDim = 0; jDim < nDim; jDim++) {
//      pi[iDim][jDim] = eta[iDim][jDim] - 1.0/rho*tau[iDim][jDim];
//    }
//  }
//  // Zeta
//  for (iDim = 0; iDim < nDim; iDim++) {
//    for (jDim = 0; jDim < nDim; jDim++) {
//      zeta[iDim][jDim] = vel[jDim]*GInvRho[iDim] - 2.0/3.0*vel[iDim]*GInvRho[jDim];
//    }
//  }
//  
//  /*--- Calculate supporting quantities ---*/
//  // GradPhi dot GradInvRho - e.g. (diPhix*di(1/rho), diPhiy*di(1/rho), ... )
//  for (iDim = 0; iDim < nDim; iDim++) {
//    GPhiGInvRho[iDim] = 0.0;
//    for (jDim = 0; jDim < nDim; jDim++) {
//      GPhiGInvRho[iDim] += GPsi[nSpecies+iDim][jDim]*GInvRho[jDim];
//    }
//  }
//  // (di(PsiE)*(Zetai1 + 1/rho*Taui1), di(PsiE)*(Zetai2 + 1/rho*Taui2), ... )
//  for (iDim = 0; iDim < nDim; iDim++) {
//    GPsiEZetaTau[iDim] = 0.0;
//    for (jDim = 0; jDim < nDim; jDim++) {
//      GPsiEZetaTau[iDim] += GPsi[nSpecies+nDim][jDim] * (zeta[jDim][iDim] +
//                                                         1/rho*tau[jDim][iDim]);
//    }
//  }
//  // GPhi:Eta
//  GPhiEta = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++)
//    for (jDim = 0; jDim < nDim; jDim++)
//      GPhiEta += GPsi[nSpecies+iDim][jDim]*eta[jDim][iDim];
//  // GPsiE_i*(u_j*pi_ij)
//  GPsiEvelPi = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++)
//    for (jDim = 0; jDim < nDim; jDim++)
//      GPsiEvelPi += GPsi[nSpecies+nDim][iDim]*vel[jDim]*pi[iDim][jDim];
//  
//  /*--- Calculate diffusion-related useful values ---*/
//  // Grad(Y)
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (iDim = 0; iDim < nDim; iDim++)
//      GY[iSpecies][iDim] = 1/rho*(GV[RHOS_INDEX+iSpecies][iDim] -
//                                  Y[iSpecies]*GV[RHO_INDEX][iDim]);
//  }
//  
//  // Sum_k (Ik)
//  for (iDim = 0; iDim < nDim; iDim++) {
//    SIk[iDim] = 0.0;
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//      SIk[iDim] += -rho*Ds[iSpecies]*GY[iSpecies][iDim];
//  }
//  
//  // Diffusion flux, J
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (iDim = 0; iDim < nDim; iDim++)
//      Js[iSpecies][iDim] = -rho*Ds[iSpecies]*GY[iSpecies][iDim] - Y[iSpecies]*SIk[iDim];
//    cout << "********" << endl;
//    cout << "Ds: " << Ds[iSpecies] << endl;
//    cout << "rho: " << rho << endl;
//    cout << "GY: " << GY[iSpecies][0] << "\t" << GY[iSpecies][1] << "\t" << GY[iSpecies][2] << endl;
//    cout << "Y: " << Y[iSpecies] << endl;
//    cout << "SIk: " << SIk[0] << "\t" << SIk[1] << "\t" << SIk[2] << endl;
//    cout << "J: " << Js[iSpecies][0] << "\t" << Js[iSpecies][1] << "\t" << Js[iSpecies][2] << endl;
//  }
//  cin.get();
//  
//  // Initialize arrays
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
//      for (iDim = 0; iDim < nDim; iDim++) {
//        dIdr[iSpecies][jSpecies][iDim] = 0.0;
//        dJdr[iSpecies][jSpecies][iDim] = 0.0;
//      }
//  
//  // dI/dr
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        dIdr[jSpecies][iSpecies][iDim] += -Ds[jSpecies]/rho*Y[jSpecies]*GV[RHO_INDEX][iDim];
//      }
//    }
//    for (iDim = 0; iDim < nDim; iDim++)
//      dIdr[iSpecies][iSpecies][iDim] += Ds[iSpecies]/rho*GV[RHO_INDEX][iDim];
//  }
//  
//  // Sum_k (dIk/dr)
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (iDim = 0; iDim < nDim; iDim++) {
//      SdIdr[iSpecies][iDim] = 0.0;
//      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//        SdIdr[iSpecies][iDim] += dIdr[jSpecies][iSpecies][iDim];
//      }
//    }
//  }
//  
//  // dJ/dr
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        dJdr[jSpecies][iSpecies][iDim] += (dIdr[jSpecies][iSpecies][iDim]
//                                           +1/rho*Y[jSpecies]*SIk[iDim]
//                                           -Y[jSpecies]*SdIdr[iSpecies][iDim]);
//      }
//    }
//    for (iDim = 0; iDim < nDim; iDim++)
//      dJdr[iSpecies][iSpecies][iDim] += -1/rho*SIk[iDim];
//  }
//  
//
//  
//
//  /*--- Contribution to viscous residual from Av1 ---*/
//  // Species mass
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (jSpecies =0; jSpecies < nSpecies; jSpecies++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        val_residual[iSpecies] +=
//            -GPsi[jSpecies][iDim]       *(dJdr[jSpecies][iSpecies][iDim])*Volume
//            -GPsi[nSpecies+nDim][iDim]  *(dJdr[jSpecies][iSpecies][iDim]*hs[jSpecies] +
//                                          Js[jSpecies][iDim]*
//                                          ((Ru/Ms[jSpecies]+Cvtrs[jSpecies])*dTdU_i[iSpecies]+
//                                           Cvves[jSpecies]*dTvedU_i[iSpecies]))*Volume
//            -GPsi[nSpecies+nDim+1][iDim]*(dJdr[jSpecies][iSpecies][iDim]*eves[jSpecies] +
//                                          Js[jSpecies][iDim]*
//                                          (Cvves[jSpecies]*dTvedU_i[iSpecies]))*Volume;
//      }
//    }
//  }
//  // Remaining terms
//  for (iVar = nSpecies; iVar < nVar; iVar++) {
//    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        val_residual[iVar] +=
//            -GPsi[nSpecies+nDim][iDim]*(Js[jSpecies][iDim]*
//                                        ((Ru/Ms[jSpecies]+Cvtrs[jSpecies])*dTdU_i[iVar]+
//                                         Cvves[jSpecies]*dTvedU_i[iVar]))*Volume
//            -GPsi[nSpecies+nDim+1][iDim]*(Js[jSpecies][iDim]*
//                                          (Cvves[jSpecies]*dTvedU_i[iVar]))*Volume;
//      }
//    }
//  }
//  
//  /*--- Contribution to viscous residual from Av2 ---*/
//  // Mass conservation
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    val_residual[iSpecies] += (-GPhiEta + GPsiEvelPi)*mu2*Volume;
//  }
//  if (iDim == 2) {
//    // X-momentum
//    val_residual[nSpecies]   += ((GPhiGInvRho[0] + 1.0/3.0*GPsi[nSpecies][0]*GInvRho[0]) +
//                                 (GPsi[nSpecies+1][0]*GInvRho[1] - 2.0/3.0*GPsi[nSpecies+1][1]*GInvRho[0]) +
//                                 (GPsi[nSpecies+nDim][0]*velGInvRho + GPsiEZetaTau[0]))*mu2*Volume;
//    // Y-momentum
//    val_residual[nSpecies+1] += ((GPsi[nSpecies][1]*GInvRho[0]   - 2.0/3.0*GPsi[nSpecies][0]*GInvRho[1]) +
//                                 (GPhiGInvRho[1] + 1.0/3.0*GPsi[nSpecies+1][1]*GInvRho[1]) +
//                                 (GPsi[nSpecies+nDim][1]*velGInvRho + GPsiEZetaTau[1]))*mu2*Volume;
//  } else {
//    // X-momentum
//    val_residual[nSpecies]   += ((GPhiGInvRho[0] + 1.0/3.0*GPsi[nSpecies][0]*GInvRho[0]) +
//                                 (GPsi[nSpecies+1][0]*GInvRho[1] - 2.0/3.0*GPsi[nSpecies+1][1]*GInvRho[0]) +
//                                 (GPsi[nSpecies+2][0]*GInvRho[2] - 2.0/3.0*GPsi[nSpecies+2][2]*GInvRho[0]) +
//                                 (GPsi[nSpecies+3][0]*velGInvRho + GPsiEZetaTau[0]))*mu2*Volume;
//    // Y-momentum
//    val_residual[nSpecies+1] += ((GPsi[nSpecies][1]*GInvRho[0]   - 2.0/3.0*GPsi[nSpecies][0]*GInvRho[1]) +
//                                 (GPhiGInvRho[1] + 1.0/3.0*GPsi[nSpecies+1][1]*GInvRho[1]) +
//                                 (GPsi[nSpecies+2][1]*GInvRho[2] - 2.0/3.0*GPsi[nSpecies+2][2]*GInvRho[1]) +
//                                 (GPsi[nSpecies+3][1]*velGInvRho + GPsiEZetaTau[1]))*mu2*Volume;
//    // Z-momentum
//    val_residual[nSpecies+2] += ((GPsi[nSpecies][2]*GInvRho[0]   - 2.0/3.0*GPsi[nSpecies][0]*GInvRho[2]) +
//                                 (GPsi[nSpecies+1][2]*GInvRho[1] - 2.0/3.0*GPsi[nSpecies+1][1]*GInvRho[2]) +
//                                 (GPhiGInvRho[2] + 1.0/3.0*GPsi[nSpecies+2][2]*GInvRho[2]) +
//                                 (GPsi[nSpecies+3][2]*velGInvRho + GPsiEZetaTau[2]))*mu2*Volume;
//  }
//  
//  for (iDim = 0; iDim < nDim; iDim++) {
//    // Mass conservation
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      
//      for (jDim = 0; jDim < nDim; jDim++)
//        val_residual[iSpecies] += GPsi[nSpecies+nDim][iDim] * (1/rhoCvtr * (vel[jDim]*PrimVar_Grad_i[VEL_INDEX+jDim][iDim])) * mu3 * Volume;
//      val_residual[iSpecies] +=
//          GPsi[nSpecies+nDim][iDim] * (1/rhoCvtr * (Cvtrs[iSpecies]*PrimVar_Grad_i[T_INDEX][iDim])
//                                       -dTdU_i[iSpecies]/rhoCvtr *
//                                       PrimVar_Grad_i[RHOCVTR_INDEX][iDim])* mu3 * Volume;
//      
////      val_residual[iSpecies] += GPsi[nSpecies+3][iDim] * (1/rhoCvtr * (vel[0]*PrimVar_Grad_i[VEL_INDEX][iDim]   +
////                                                                       vel[1]*PrimVar_Grad_i[VEL_INDEX+1][iDim] +
////                                                                       vel[2]*PrimVar_Grad_i[VEL_INDEX+2][iDim] -
////                                                                       Cvtrs[iSpecies]*PrimVar_Grad_i[T_INDEX][iDim])
////                                                          -dTdU_i[iSpecies]/rhoCvtr *
////                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim])* mu3 * Volume;
//    }
//    
//    // Momentum eqns.
//    for (jDim = 0; jDim < nDim; jDim++) {
//      val_residual[nSpecies+jDim] +=
//          GPsi[nSpecies+nDim][iDim] * (-1/rhoCvtr*PrimVar_Grad_i[VEL_INDEX+jDim][iDim]
//                                       -dTdU_i[nSpecies+jDim]/rhoCvtr *
//                                       PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
//    }
////    // X-momentum
////    val_residual[nSpecies]   += GPsi[nSpecies+3][iDim] * (-1/rhoCvtr*PrimVar_Grad_i[VEL_INDEX][iDim]
////                                                          -dTdU_i[nSpecies]/rhoCvtr *
////                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
////    // Y-momentum
////    val_residual[nSpecies+1] += GPsi[nSpecies+3][iDim] * (-1/rhoCvtr*PrimVar_Grad_i[VEL_INDEX+1][iDim]
////                                                          -dTdU_i[nSpecies+1]/rhoCvtr *
////                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
////    // Z-momentum
////    val_residual[nSpecies+2] += GPsi[nSpecies+3][iDim] * (-1/rhoCvtr*PrimVar_Grad_i[VEL_INDEX+2][iDim]
////                                                          -dTdU_i[nSpecies+2]/rhoCvtr *
////                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
//    // Energy
//    val_residual[nSpecies+nDim] +=
//        GPsi[nSpecies+nDim][iDim] * (-dTdU_i[nSpecies+nDim]/rhoCvtr *
//                                     PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
//    // Vib.-el. energy
//    val_residual[nSpecies+nDim+1] +=
//        GPsi[nSpecies+nDim][iDim] * (-dTdU_i[nSpecies+nDim+1]/rhoCvtr *
//                                     PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
//  }
//  
//  /*--- Contribution to viscous residual from Av4 ---*/
//  for (iDim = 0; iDim < nDim; iDim++) {
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      val_residual[iSpecies] += ((GPsi[nSpecies+nDim][iDim]+GPsi[nSpecies+nDim+1][iDim]) *
//                                 (-Cvves[iSpecies]/rhoCvve*PrimVar_Grad_i[TVE_INDEX][iDim]
//                                  -dTvedU_i[iSpecies]/rhoCvve*PrimVar_Grad_i[RHOCVVE_INDEX][iDim])) *
//                                mu4 * Volume;
//    }
//    val_residual[nSpecies+nDim+1] += ((GPsi[nSpecies+nDim][iDim]+GPsi[nSpecies+nDim+1][iDim]) *
//                                      (-dTvedU_i[nSpecies+nDim+1]/rhoCvve * PrimVar_Grad_i[RHOCVVE_INDEX][iDim])) *
//                                     mu4 * Volume;
//    
//  }
}

//void CSource_AdjTNE2::ComputeSourceViscous (double *val_residual, CConfig *config) {
//  
//	unsigned short iDim, jDim, iSpecies, jSpecies, iVar, jVar;
//  double rho, sqvel, rhoCvtr, rhoCvve, YDbar, *Ms, *xi, Ru;
//  double div_vel, div_velorho, velGInvRho, GPhiEta, GPsiEvelPi;
//  double GPsiEGr, GPsiEveGr, DGPsiEGrs, DGPsiEveGrs;
//  double *Ds, mu2, mu3, mu4;
//  double **GPsi, **GV;
//  
//  /*--- Initialize arrays ---*/
//  for (iVar = 0; iVar < nVar; iVar++) {
//    val_residual[iVar] = 0.0;
//  }
//  
//  /*--- Rename for convenience ---*/
//  Ms      = config->GetMolar_Mass();
//  xi      = config->GetRotationModes();
//  Ru      = UNIVERSAL_GAS_CONSTANT;
//  rhoCvtr = V_i[RHOCVTR_INDEX];
//  rhoCvve = V_i[RHOCVVE_INDEX];
//  GPsi    = PsiVar_Grad_i;
//  GV      = PrimVar_Grad_i;
//  
//  /*--- Assign useful values ---*/
//  // Mixture & species density, mass fraction
//  rho = V_i[RHO_INDEX];
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    rhos[iSpecies] = V_i[RHOS_INDEX+iSpecies];
//    Y[iSpecies]    = rhos[iSpecies]/rho;
//  }
//  // Velocity
//  sqvel = 0.0;
//  div_vel = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++) {
//		vel[iDim] = V_i[VEL_INDEX+iDim];
//		sqvel += vel[iDim]*vel[iDim];
//    div_vel += PrimVar_Grad_i[VEL_INDEX+iDim][iDim];
//	}
//  // Species enthalpy and vib.-el. energy
//  for (iSpecies     = 0; iSpecies < nSpecies; iSpecies++) {
//    eves[iSpecies]  = var->CalcEve(config, V_i[TVE_INDEX], iSpecies);
//    hs[iSpecies]    = var->CalcHs(config, V_i[T_INDEX], eves[iSpecies], iSpecies);
//    Cvves[iSpecies] = var->CalcCvve(V_i[TVE_INDEX], config, iSpecies);
//    Cvtrs[iSpecies] = (3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
//  }
//  
//  // Transport coefficients
//  mu2 = Laminar_Viscosity_i;
//  mu3 = Thermal_Conductivity_i;
//  mu4 = Thermal_Conductivity_ve_i;
//  Ds  = Diffusion_Coeff_i;
//  
//  /*--- Gradients of flow variables ---*/
//  for (iDim = 0; iDim < nDim; iDim++) {
//    GInvRho[iDim] = -1.0/(rho*rho)*PrimVar_Grad_i[RHO_INDEX][iDim];
//    for (jDim = 0; jDim < nDim; jDim++)
//      GVeloRho[iDim][jDim] = -vel[iDim]/(rho*rho)*(PrimVar_Grad_i[RHO_INDEX][jDim])
//      +  1.0/rho*(PrimVar_Grad_i[VEL_INDEX+iDim][jDim]);
//  }
//  div_velorho = 0.0;
//  velGInvRho = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++) {
//    div_velorho += GVeloRho[iDim][iDim];
//    velGInvRho += vel[iDim]*GInvRho[iDim];
//  }
//  
//  /*--- Supporting matrices ---*/
//  // Tau: Stress tensor matrix
//  for (iDim = 0; iDim < nDim; iDim++) {
//    for (jDim = 0; jDim < nDim; jDim++) {
//      tau[iDim][jDim] = PrimVar_Grad_i[VEL_INDEX+iDim][jDim] +
//      PrimVar_Grad_i[VEL_INDEX+jDim][iDim];
//    }
//    tau[iDim][iDim] -= 2.0/3.0*div_vel;
//  }
//  // Eta
//  for (iDim = 0; iDim < nDim; iDim++) {
//    for (jDim = 0; jDim < nDim; jDim++) {
//      eta[iDim][jDim] = GVeloRho[iDim][jDim] + GVeloRho[jDim][iDim];
//    }
//    eta[iDim][iDim] -= 2.0/3.0*div_velorho;
//  }
//  // Pi
//  for (iDim = 0; iDim < nDim; iDim++) {
//    for (jDim = 0; jDim < nDim; jDim++) {
//      pi[iDim][jDim] = eta[iDim][jDim] - 1.0/rho*tau[iDim][jDim];
//    }
//  }
//  // Zeta
//  for (iDim = 0; iDim < nDim; iDim++) {
//    for (jDim = 0; jDim < nDim; jDim++) {
//      zeta[iDim][jDim] = vel[jDim]*GInvRho[iDim] - 2.0/3.0*vel[iDim]*GInvRho[jDim];
//    }
//  }
//  
//  /*--- Calculate supporting quantities ---*/
//  // GradPhi dot GradInvRho - e.g. (diPhix*di(1/rho), diPhiy*di(1/rho), ... )
//  for (iDim = 0; iDim < nDim; iDim++) {
//    GPhiGInvRho[iDim] = 0.0;
//    for (jDim = 0; jDim < nDim; jDim++) {
//      GPhiGInvRho[iDim] += GPsi[nSpecies+iDim][jDim]*GInvRho[jDim];
//    }
//  }
//  // (di(PsiE)*(Zetai1 + 1/rho*Taui1), di(PsiE)*(Zetai2 + 1/rho*Taui2), ... )
//  for (iDim = 0; iDim < nDim; iDim++) {
//    GPsiEZetaTau[iDim] = 0.0;
//    for (jDim = 0; jDim < nDim; jDim++) {
//      GPsiEZetaTau[iDim] += GPsi[nSpecies+nDim][jDim] * (zeta[jDim][iDim] +
//                                                         1/rho*tau[jDim][iDim]);
//    }
//  }
//  // GPhi:Eta
//  GPhiEta = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++)
//    for (jDim = 0; jDim < nDim; jDim++)
//      GPhiEta += GPsi[nSpecies+iDim][jDim]*eta[jDim][iDim];
//  // GPsiE_i*(u_j*pi_ij)
//  GPsiEvelPi = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++)
//    for (jDim = 0; jDim < nDim; jDim++)
//      GPsiEvelPi += GPsi[nSpecies+nDim][iDim]*vel[jDim]*pi[iDim][jDim];
//  
//  /*--- Calculate diffusion-related useful values ---*/
//  // Grad(Y)
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (iDim = 0; iDim < nDim; iDim++)
//      GY[iSpecies][iDim] = 1/rho*(GV[RHOS_INDEX+iSpecies][iDim] -
//                                  Y[iSpecies]*GV[RHO_INDEX][iDim]);
//  }
//  
//  // Sum_k (Ik)
//  for (iDim = 0; iDim < nDim; iDim++) {
//    SIk[iDim] = 0.0;
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//      SIk[iDim] += -rho*Ds[iSpecies]*GY[iSpecies][iDim];
//  }
//  
//  // Diffusion flux, J
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    for (iDim = 0; iDim < nDim; iDim++)
//      Js[iSpecies][iDim] = -rho*Ds[iSpecies]*GY[iSpecies][iDim] - Y[iSpecies]*SIk[iDim];
//  
//  // Initialize arrays
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
//      for (iDim = 0; iDim < nDim; iDim++) {
//        dIdr[iSpecies][jSpecies][iDim] = 0.0;
//        dJdr[iSpecies][jSpecies][iDim] = 0.0;
//      }
//  
//  // dI/dr
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        dIdr[jSpecies][iSpecies][iDim] += -Ds[jSpecies]/rho*Y[jSpecies]*GV[RHO_INDEX][iDim];
//      }
//    }
//    for (iDim = 0; iDim < nDim; iDim++)
//      dIdr[iSpecies][iSpecies][iDim] += Ds[iSpecies]/rho*GV[RHO_INDEX][iDim];
//  }
//  
//  // Sum_k (dIk/dr)
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (iDim = 0; iDim < nDim; iDim++) {
//      SdIdr[iSpecies][iDim] = 0.0;
//      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//        SdIdr[iSpecies][iDim] += dIdr[jSpecies][iSpecies][iDim];
//      }
//    }
//  }
//  
//  // dJ/dr
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        dJdr[jSpecies][iSpecies][iDim] += (dIdr[jSpecies][iSpecies][iDim]
//                                           +1/rho*Y[jSpecies]*SIk[iDim]
//                                           -Y[jSpecies]*SdIdr[iSpecies][iDim]);
//      }
//    }
//    for (iDim = 0; iDim < nDim; iDim++)
//      dJdr[iSpecies][iSpecies][iDim] += -1/rho*SIk[iDim];
//  }
//  
//  //////////////////// DEBUG ////////////////////
//  double UnitNormal[3], tmp, tmp2;
//  double *Fv_old, *Fv_new, d;
//  UnitNormal[0] = 1.0;
//  UnitNormal[1] = 0.0;
//  UnitNormal[2] = 0.0;
//  Fv_new = new double[nVar];
//  Fv_old = new double[nVar];
//
//  /*--- Contribution to viscous residual from Av1 ---*/
//  double **Av1;
//  Av1 = new double*[nVar];
//  for (iVar = 0; iVar < nVar; iVar++)
//    Av1[iVar] = new double[nVar];
//
//  for (iVar = 0; iVar < nVar; iVar++)
//    for (jVar = 0; jVar < nVar; jVar++)
//      Av1[iVar][jVar] = 0.0;
//
//  //
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (jSpecies =0; jSpecies < nSpecies; jSpecies++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        Av1[iSpecies][jSpecies] +=
//            -UnitNormal[iDim]*(dJdr[iSpecies][jSpecies][iDim]);
//        Av1[nSpecies+nDim][iSpecies] +=
//            -UnitNormal[iDim]*(dJdr[jSpecies][iSpecies][iDim]*hs[jSpecies] +
//                               Js[jSpecies][iDim]*
//                               ((Ru/Ms[jSpecies]+Cvtrs[jSpecies])*dTdU_i[iSpecies]+
//                                Cvves[jSpecies]*dTvedU_i[iSpecies]));
//        Av1[nSpecies+nDim+1][iSpecies] +=
//            -UnitNormal[iDim]*(dJdr[jSpecies][iSpecies][iDim]*eves[jSpecies] +
//                               Js[jSpecies][iDim]*
//                               (Cvves[jSpecies]*dTvedU_i[iSpecies]));
//      }
//    }
//  }
//  // Remaining terms
//  for (iVar = nSpecies; iVar < nVar; iVar++) {
//    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        Av1[nSpecies+nDim][iVar] +=
//        -UnitNormal[iDim]*(Js[jSpecies][iDim]*
//                           ((Ru/Ms[jSpecies]+Cvtrs[jSpecies])*dTdU_i[iVar]+
//                            Cvves[jSpecies]*dTvedU_i[iVar]));
//        Av1[nSpecies+nDim+1][iVar] +=
//        -UnitNormal[iDim]*(Js[jSpecies][iDim]*
//                           (Cvves[jSpecies]*dTvedU_i[iVar]));
//      }
//    }
//  }
//
//  cout << endl << "Av1: " << endl;
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nVar; jVar++) {
//      cout << Av1[jVar][iVar] << "\t";
//    }
//    cout << endl;
//  }
//
//  cout << endl << "FD: " << endl;
//  // finite difference gradient
//  for (iVar = 0; iVar < nVar; iVar++) {
//    // set displacement value
//    d = 0.0001*U_i[iVar];
//
//    // calculate viscous flux
//    GetViscousProjFlux(V_i, GV, UnitNormal, Ds, mu2, mu3, mu4, config);
//
//    // copy solution
//    for (jVar = 0; jVar < nVar; jVar++)
//      Fv_old[jVar] = Proj_Flux_Tensor[jVar];
//
//    // perturb solution
//    U_i[iVar] += d;
//    var->Cons2PrimVar(config, U_i, V_i, dPdU_i, dTdU_i, dTvedU_i);
//
//    // calculate viscous flux
//    GetViscousProjFlux(V_i, GV, UnitNormal, Ds, mu2, mu3, mu4, config);
//
//    // copy solution
//    for (jVar = 0; jVar < nVar; jVar++)
//      Fv_new[jVar] = Proj_Flux_Tensor[jVar];
//
//    // return solution to original value
//    U_i[iVar] -= d;
//    var->Cons2PrimVar(config, U_i, V_i, dPdU_i, dTdU_i, dTvedU_i);
//
//    // display FD gradient
//    for (jVar = 0; jVar < nVar; jVar++)
//      cout << (Fv_new[jVar]-Fv_old[jVar])/d << "\t";
//    cout << endl;
//  }
//
//  cin.get();
//  
//  //////////////////// DEBUG ////////////////////
//  
//  
//  /*--- Contribution to viscous residual from Av1 ---*/
//  // Species mass
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (jSpecies =0; jSpecies < nSpecies; jSpecies++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        val_residual[iSpecies] +=
//        -GPsi[jSpecies][iDim]       *(dJdr[jSpecies][iSpecies][iDim])*Volume
//        -GPsi[nSpecies+nDim][iDim]  *(dJdr[jSpecies][iSpecies][iDim]*hs[jSpecies] +
//                                      Js[jSpecies][iDim]*
//                                      ((Ru/Ms[jSpecies]+Cvtrs[jSpecies])*dTdU_i[iSpecies]+
//                                       Cvves[jSpecies]*dTvedU_i[iSpecies]))*Volume
//        -GPsi[nSpecies+nDim+1][iDim]*(dJdr[jSpecies][iSpecies][iDim]*eves[jSpecies] +
//                                      Js[jSpecies][iDim]*
//                                      (Cvves[jSpecies]*dTvedU_i[iSpecies]))*Volume;
//      }
//    }
//  }
//  // Remaining terms
//  for (iVar = nSpecies; iVar < nVar; iVar++) {
//    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//      for (iDim = 0; iDim < nDim; iDim++) {
//        val_residual[iVar] +=
//        -GPsi[nSpecies+nDim][iDim]*(Js[jSpecies][iDim]*
//                                    ((Ru/Ms[jSpecies]+Cvtrs[jSpecies])*dTdU_i[iVar]+
//                                     Cvves[jSpecies]*dTvedU_i[iVar]))*Volume
//        -GPsi[nSpecies+nDim+1][iDim]*(Js[jSpecies][iDim]*
//                                      (Cvves[jSpecies]*dTvedU_i[iVar]))*Volume;
//      }
//    }
//  }
//  
//  /*--- Contribution to viscous residual from Av2 ---*/
//  // Mass conservation
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    val_residual[iSpecies] += (-GPhiEta + GPsiEvelPi)*mu2*Volume;
//  }
//  if (iDim == 2) {
//    // X-momentum
//    val_residual[nSpecies]   += ((GPhiGInvRho[0] + 1.0/3.0*GPsi[nSpecies][0]*GInvRho[0]) +
//                                 (GPsi[nSpecies+1][0]*GInvRho[1] - 2.0/3.0*GPsi[nSpecies+1][1]*GInvRho[0]) +
//                                 (GPsi[nSpecies+nDim][0]*velGInvRho + GPsiEZetaTau[0]))*mu2*Volume;
//    // Y-momentum
//    val_residual[nSpecies+1] += ((GPsi[nSpecies][1]*GInvRho[0]   - 2.0/3.0*GPsi[nSpecies][0]*GInvRho[1]) +
//                                 (GPhiGInvRho[1] + 1.0/3.0*GPsi[nSpecies+1][1]*GInvRho[1]) +
//                                 (GPsi[nSpecies+nDim][1]*velGInvRho + GPsiEZetaTau[1]))*mu2*Volume;
//  } else {
//    // X-momentum
//    val_residual[nSpecies]   += ((GPhiGInvRho[0] + 1.0/3.0*GPsi[nSpecies][0]*GInvRho[0]) +
//                                 (GPsi[nSpecies+1][0]*GInvRho[1] - 2.0/3.0*GPsi[nSpecies+1][1]*GInvRho[0]) +
//                                 (GPsi[nSpecies+2][0]*GInvRho[2] - 2.0/3.0*GPsi[nSpecies+2][2]*GInvRho[0]) +
//                                 (GPsi[nSpecies+3][0]*velGInvRho + GPsiEZetaTau[0]))*mu2*Volume;
//    // Y-momentum
//    val_residual[nSpecies+1] += ((GPsi[nSpecies][1]*GInvRho[0]   - 2.0/3.0*GPsi[nSpecies][0]*GInvRho[1]) +
//                                 (GPhiGInvRho[1] + 1.0/3.0*GPsi[nSpecies+1][1]*GInvRho[1]) +
//                                 (GPsi[nSpecies+2][1]*GInvRho[2] - 2.0/3.0*GPsi[nSpecies+2][2]*GInvRho[1]) +
//                                 (GPsi[nSpecies+3][1]*velGInvRho + GPsiEZetaTau[1]))*mu2*Volume;
//    // Z-momentum
//    val_residual[nSpecies+2] += ((GPsi[nSpecies][2]*GInvRho[0]   - 2.0/3.0*GPsi[nSpecies][0]*GInvRho[2]) +
//                                 (GPsi[nSpecies+1][2]*GInvRho[1] - 2.0/3.0*GPsi[nSpecies+1][1]*GInvRho[2]) +
//                                 (GPhiGInvRho[2] + 1.0/3.0*GPsi[nSpecies+2][2]*GInvRho[2]) +
//                                 (GPsi[nSpecies+3][2]*velGInvRho + GPsiEZetaTau[2]))*mu2*Volume;
//  }
//  
//  
//  /*--- Contribution to viscous residual from Av3 ---*/
//  for (iDim = 0; iDim < nDim; iDim++) {
//    // Mass conservation
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      
//      for (jDim = 0; jDim < nDim; jDim++)
//        val_residual[iSpecies] += GPsi[nSpecies+nDim][iDim] * (1/rhoCvtr * (vel[jDim]*PrimVar_Grad_i[VEL_INDEX+jDim][iDim])) * mu3 * Volume;
//      val_residual[iSpecies] +=
//      GPsi[nSpecies+nDim][iDim] * (1/rhoCvtr * (Cvtrs[iSpecies]*PrimVar_Grad_i[T_INDEX][iDim])
//                                   -dTdU_i[iSpecies]/rhoCvtr *
//                                   PrimVar_Grad_i[RHOCVTR_INDEX][iDim])* mu3 * Volume;
//      
//      //      val_residual[iSpecies] += GPsi[nSpecies+3][iDim] * (1/rhoCvtr * (vel[0]*PrimVar_Grad_i[VEL_INDEX][iDim]   +
//      //                                                                       vel[1]*PrimVar_Grad_i[VEL_INDEX+1][iDim] +
//      //                                                                       vel[2]*PrimVar_Grad_i[VEL_INDEX+2][iDim] -
//      //                                                                       Cvtrs[iSpecies]*PrimVar_Grad_i[T_INDEX][iDim])
//      //                                                          -dTdU_i[iSpecies]/rhoCvtr *
//      //                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim])* mu3 * Volume;
//    }
//    
//    // Momentum eqns.
//    for (jDim = 0; jDim < nDim; jDim++) {
//      val_residual[nSpecies+jDim] +=
//      GPsi[nSpecies+nDim][iDim] * (-1/rhoCvtr*PrimVar_Grad_i[VEL_INDEX+jDim][iDim]
//                                   -dTdU_i[nSpecies+jDim]/rhoCvtr *
//                                   PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
//    }
//    //    // X-momentum
//    //    val_residual[nSpecies]   += GPsi[nSpecies+3][iDim] * (-1/rhoCvtr*PrimVar_Grad_i[VEL_INDEX][iDim]
//    //                                                          -dTdU_i[nSpecies]/rhoCvtr *
//    //                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
//    //    // Y-momentum
//    //    val_residual[nSpecies+1] += GPsi[nSpecies+3][iDim] * (-1/rhoCvtr*PrimVar_Grad_i[VEL_INDEX+1][iDim]
//    //                                                          -dTdU_i[nSpecies+1]/rhoCvtr *
//    //                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
//    //    // Z-momentum
//    //    val_residual[nSpecies+2] += GPsi[nSpecies+3][iDim] * (-1/rhoCvtr*PrimVar_Grad_i[VEL_INDEX+2][iDim]
//    //                                                          -dTdU_i[nSpecies+2]/rhoCvtr *
//    //                                                          PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
//    // Energy
//    val_residual[nSpecies+nDim] +=
//    GPsi[nSpecies+nDim][iDim] * (-dTdU_i[nSpecies+nDim]/rhoCvtr *
//                                 PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
//    // Vib.-el. energy
//    val_residual[nSpecies+nDim+1] +=
//    GPsi[nSpecies+nDim][iDim] * (-dTdU_i[nSpecies+nDim+1]/rhoCvtr *
//                                 PrimVar_Grad_i[RHOCVTR_INDEX][iDim]) * mu3 * Volume;
//  }
//  
//  /*--- Contribution to viscous residual from Av4 ---*/
//  for (iDim = 0; iDim < nDim; iDim++) {
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      val_residual[iSpecies] += ((GPsi[nSpecies+nDim][iDim]+GPsi[nSpecies+nDim+1][iDim]) *
//                                 (-Cvves[iSpecies]/rhoCvve*PrimVar_Grad_i[TVE_INDEX][iDim]
//                                  -dTvedU_i[iSpecies]/rhoCvve*PrimVar_Grad_i[RHOCVVE_INDEX][iDim])) *
//      mu4 * Volume;
//    }
//    val_residual[nSpecies+nDim+1] += ((GPsi[nSpecies+nDim][iDim]+GPsi[nSpecies+nDim+1][iDim]) *
//                                      (-dTvedU_i[nSpecies+nDim+1]/rhoCvve * PrimVar_Grad_i[RHOCVVE_INDEX][iDim])) *
//    mu4 * Volume;
//    
//  }
//}


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
