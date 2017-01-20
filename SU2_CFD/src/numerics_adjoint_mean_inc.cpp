/*!
 * \file numerics_adjoint_mean_inc.cpp
 * \brief This file contains the numerical methods for adjoint incompressible flow.
 * \author F. Palacios, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

CUpwRoeArtComp_AdjFlow::CUpwRoeArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  
  MeanVelocity = new su2double [nDim];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  Lambda = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  Proj_Jac_Tensor_i = new su2double*[nVar];
  Proj_Jac_Tensor_j = new su2double*[nVar];
  Proj_ModJac_Tensor = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
    Proj_Jac_Tensor_i[iVar] = new su2double[nVar];
    Proj_Jac_Tensor_j[iVar] = new su2double[nVar];
    Proj_ModJac_Tensor[iVar] = new su2double[nVar];
  }
  
}

CUpwRoeArtComp_AdjFlow::~CUpwRoeArtComp_AdjFlow(void) {
  
  delete [] MeanVelocity;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] Lambda;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
    delete [] Proj_Jac_Tensor_i[iVar];
    delete [] Proj_Jac_Tensor_j[iVar];
    delete [] Proj_ModJac_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  delete [] Proj_Jac_Tensor_i;
  delete [] Proj_Jac_Tensor_j;
  delete [] Proj_ModJac_Tensor;
  
}

void CUpwRoeArtComp_AdjFlow::ComputeResidual (su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii,
                                          su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) {
  
  /*--- Compute face area ---*/
  
  Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*--- Compute and unitary normal vector ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    UnitNormal[iDim] = Normal[iDim]/Area;
    if (fabs(UnitNormal[iDim]) < EPS) UnitNormal[iDim] = EPS;
  }
  
  /*--- Set the variables at point i, and j ---*/
  
  Pressure_i = V_i[0];          Pressure_j = V_j[0];
  DensityInc_i =  V_i[nDim+1];  DensityInc_j = V_j[nDim+1];
  BetaInc2_i = V_i[nDim+2];     BetaInc2_j = V_j[nDim+2];

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
  }
  
  /*--- Jacobians of the inviscid flux, scaled by 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
  
  GetInviscidArtCompProjJac(&DensityInc_i, Velocity_i, &BetaInc2_i, Normal, 0.5, Proj_Jac_Tensor_i);
  GetInviscidArtCompProjJac(&DensityInc_j, Velocity_j, &BetaInc2_j, Normal, 0.5, Proj_Jac_Tensor_j);
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual_i[iVar] = 0.0; val_residual_j[iVar] = 0.0;
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual_i[iVar] += Proj_Jac_Tensor_i[jVar][iVar]*(Psi_i[jVar] + Psi_j[jVar]);
      val_residual_j[iVar] -= Proj_Jac_Tensor_j[jVar][iVar]*(Psi_i[jVar] + Psi_j[jVar]);
    }
  }
  
  /*--- Mean variables at points iPoint and jPoint ---*/
  
  MeanDensity = 0.5*(DensityInc_i + DensityInc_j);
  MeanPressure = 0.5*(Pressure_i + Pressure_j);
  MeanBetaInc2 = 0.5*(BetaInc2_i + BetaInc2_j);
  
  ProjVelocity = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim] + Velocity_j[iDim]);
    ProjVelocity += MeanVelocity[iDim]*Normal[iDim];
  }
  
  MeanSoundSpeed = sqrt(ProjVelocity*ProjVelocity + (MeanBetaInc2/MeanDensity) * Area * Area);
  
  /*--- Compute P, inverse P, and store eigenvalues ---*/
  
  GetPArtCompMatrix_inv(&MeanDensity, MeanVelocity, &MeanBetaInc2, UnitNormal, invP_Tensor);
  GetPArtCompMatrix(&MeanDensity, MeanVelocity, &MeanBetaInc2, UnitNormal, P_Tensor);
  
  /*--- Flow eigenvalues ---*/
  
  if (nDim == 2) {
    Lambda[0] = ProjVelocity;
    Lambda[1] = ProjVelocity + MeanSoundSpeed;
    Lambda[2] = ProjVelocity - MeanSoundSpeed;
  }
  if (nDim == 3) {
    Lambda[0] = ProjVelocity;
    Lambda[1] = ProjVelocity;
    Lambda[2] = ProjVelocity + MeanSoundSpeed;
    Lambda[3] = ProjVelocity - MeanSoundSpeed;
  }
  
  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = fabs(Lambda[iVar]);
  
  
  /*--- Flux approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_ij = 0.0;
      
      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
      Proj_ModJac_Tensor[iVar][jVar] = 0.5*Proj_ModJac_Tensor_ij;
    }
  }
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual_i[iVar] -= Proj_ModJac_Tensor[jVar][iVar]*(Psi_i[jVar] - Psi_j[jVar]);
      val_residual_j[iVar] += Proj_ModJac_Tensor[jVar][iVar]*(Psi_i[jVar] - Psi_j[jVar]);
    }
  
  /*--- Implicit contributions, Transpose the matrices and store the Jacobians. Note the negative
   sign for the ji and jj Jacobians bc the normal direction is flipped. ---*/
  
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_ii[jVar][iVar] = Proj_Jac_Tensor_i[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar];
        val_Jacobian_ij[jVar][iVar] = Proj_Jac_Tensor_i[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar];
        val_Jacobian_ji[jVar][iVar] = -(Proj_Jac_Tensor_j[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar]);
        val_Jacobian_jj[jVar][iVar] = -(Proj_Jac_Tensor_j[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar]);
      }
    }
  }
}

CCentJSTArtComp_AdjFlow::CCentJSTArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Diff_Psi = new su2double [nVar]; Diff_Lapl = new su2double [nVar];
  Und_Lapl_i = new su2double [nVar]; Und_Lapl_j = new su2double [nVar];
  Velocity_i = new su2double [nDim]; Velocity_j = new su2double [nDim];
  Proj_Jac_Tensor_i = new su2double*[nVar];
  Proj_Jac_Tensor_j = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Jac_Tensor_i[iVar] = new su2double[nVar];
    Proj_Jac_Tensor_j[iVar] = new su2double[nVar];
  }
  
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_AdjFlow();
  Param_Kappa_4 = config->GetKappa_4th_AdjFlow();
  implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  
}

CCentJSTArtComp_AdjFlow::~CCentJSTArtComp_AdjFlow(void) {
  
  delete [] Diff_Psi; delete [] Diff_Lapl;
  delete [] Und_Lapl_i; delete [] Und_Lapl_j;
  delete [] Velocity_i; delete [] Velocity_j;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Proj_Jac_Tensor_i[iVar];
    delete [] Proj_Jac_Tensor_j[iVar];
  }
  delete [] Proj_Jac_Tensor_i;
  delete [] Proj_Jac_Tensor_j;
  
}

void CCentJSTArtComp_AdjFlow::ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                                           su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                                           CConfig *config) {
  
  /*--- Set the variables at point i ---*/
  Pressure_i = U_i[0]; Pressure_j = U_j[0];
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i;
    Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
  }
  
  /*--- Jacobians of the inviscid flux, scaled by 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
  GetInviscidArtCompProjJac(&DensityInc_i, Velocity_i, &BetaInc2_i, Normal, 0.5, Proj_Jac_Tensor_i);
  GetInviscidArtCompProjJac(&DensityInc_j, Velocity_j, &BetaInc2_j, Normal, 0.5, Proj_Jac_Tensor_j);
  
  for (iVar = 0; iVar < nDim; iVar++) {
    val_resconv_i[iVar] = 0.0; val_resconv_j[iVar] = 0.0;
    for (jVar = 0; jVar < nVar; jVar++) {
      val_resconv_i[iVar] += Proj_Jac_Tensor_i[jVar][iVar]*(Psi_i[jVar]+Psi_j[jVar]);
      val_resconv_j[iVar] -= Proj_Jac_Tensor_j[jVar][iVar]*(Psi_i[jVar]+Psi_j[jVar]);
    }
  }
  
  /*--- Jacobians of the inviscid flux ---*/
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_ii[iVar][jVar] = Proj_Jac_Tensor_i[jVar][iVar];
        val_Jacobian_ij[iVar][jVar] = Proj_Jac_Tensor_i[jVar][iVar];
        val_Jacobian_ji[iVar][jVar] = -Proj_Jac_Tensor_j[jVar][iVar];
        val_Jacobian_jj[iVar][jVar] = -Proj_Jac_Tensor_j[jVar][iVar];
      }
  }
  
  /*--- Computes differences btw. variables and Laplacians ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
    Diff_Psi[iVar]  = Psi_i[iVar]-Psi_j[iVar];
  }
  
  /*--- Compute the local espectral radius and the stretching factor ---*/
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  SoundSpeed_i = sqrt(ProjVelocity_i*ProjVelocity_i + (BetaInc2_i/DensityInc_i)*Area*Area);
  SoundSpeed_j = sqrt(ProjVelocity_j*ProjVelocity_j + (BetaInc2_j/DensityInc_j)*Area*Area);
  
  /*--- Compute the spectral radius and stretching factor ---*/
  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  /*--- Compute streching factor ---*/
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc2 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  sc4 = sc2*sc2/4.0;
  
  Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
  Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;
  
  /*--- Compute viscous residual 1st- & 3rd-order dissipation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Residual = (Epsilon_2*Diff_Psi[iVar]-Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
    val_resvisc_i[iVar] = -Residual;
    val_resvisc_j[iVar] = Residual;
    
    if (implicit) {
      val_Jacobian_ii[iVar][iVar] -= Epsilon_2 + Epsilon_4*su2double(Neighbor_i+1)*StretchingFactor*MeanLambda;
      val_Jacobian_ij[iVar][iVar] += Epsilon_2 + Epsilon_4*su2double(Neighbor_j+1)*StretchingFactor*MeanLambda;
      val_Jacobian_ji[iVar][iVar] += Epsilon_2 + Epsilon_4*su2double(Neighbor_i+1)*StretchingFactor*MeanLambda;
      val_Jacobian_jj[iVar][iVar] -= Epsilon_2 + Epsilon_4*su2double(Neighbor_j+1)*StretchingFactor*MeanLambda;
    }
  }
  
}

CCentLaxArtComp_AdjFlow::CCentLaxArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Diff_Psi = new su2double [nVar];   MeanPhi = new su2double [nDim];
  Velocity_i = new su2double [nDim]; Velocity_j = new su2double [nDim];
  Proj_Jac_Tensor_i = new su2double*[nVar];
  Proj_Jac_Tensor_j = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Jac_Tensor_i[iVar] = new su2double[nVar];
    Proj_Jac_Tensor_j[iVar] = new su2double[nVar];
  }
  
  implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_AdjFlow();
  
}

CCentLaxArtComp_AdjFlow::~CCentLaxArtComp_AdjFlow(void) {
  
  delete [] Diff_Psi; delete [] MeanPhi;
  delete [] Velocity_i; delete [] Velocity_j;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Proj_Jac_Tensor_i[iVar];
    delete [] Proj_Jac_Tensor_j[iVar];
  }
  delete [] Proj_Jac_Tensor_i;
  delete [] Proj_Jac_Tensor_j;
  
}

void CCentLaxArtComp_AdjFlow::ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                                           su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
                                           CConfig *config) {
  
  /*--- Set the variables at point i ---*/
  Pressure_i = U_i[0]; Pressure_j = U_j[0];
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i;
    Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
  }
  
  /*--- Jacobians of the inviscid flux, scaled by 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
  GetInviscidArtCompProjJac(&DensityInc_i, Velocity_i, &BetaInc2_i, Normal, 0.5, Proj_Jac_Tensor_i);
  GetInviscidArtCompProjJac(&DensityInc_j, Velocity_j, &BetaInc2_j, Normal, 0.5, Proj_Jac_Tensor_j);
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_resconv_i[iVar]  = 0.0; val_resconv_j[iVar]  = 0.0;
    for (jVar = 0; jVar < nVar; jVar++) {
      val_resconv_i[iVar] += Proj_Jac_Tensor_i[jVar][iVar]*0.5*(Psi_i[jVar]+Psi_j[jVar]);
      val_resconv_j[iVar] -= Proj_Jac_Tensor_j[jVar][iVar]*0.5*(Psi_i[jVar]+Psi_j[jVar]);
    }
  }
  
  /*--- Jacobians of the inviscid flux ---*/
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_ii[iVar][jVar] = Proj_Jac_Tensor_i[jVar][iVar];
        val_Jacobian_ij[iVar][jVar] = Proj_Jac_Tensor_i[jVar][iVar];
        val_Jacobian_ji[iVar][jVar] = -Proj_Jac_Tensor_j[jVar][iVar];
        val_Jacobian_jj[iVar][jVar] = -Proj_Jac_Tensor_j[jVar][iVar];
      }
  }
  
  /*--- Computes differences btw. variables ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_Psi[iVar] = Psi_i[iVar]-Psi_j[iVar];
  
  /*--- Compute the local espectral radius and the stretching factor ---*/
  ProjVelocity_i = 0; ProjVelocity_j = 0; Area = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  SoundSpeed_i = sqrt(ProjVelocity_i*ProjVelocity_i + (BetaInc2_i/DensityInc_i)*Area*Area);
  SoundSpeed_j = sqrt(ProjVelocity_j*ProjVelocity_j + (BetaInc2_j/DensityInc_j)*Area*Area);
  
  /*--- Compute spectral radius ---*/
  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  /*--- Compute streching factor ---*/
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc2 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc2*su2double(nDim)/3.0;
  
  /*--- Artifical dissipation evaluation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Residual = Epsilon_0*StretchingFactor*MeanLambda*Diff_Psi[iVar];
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

CAvgGradArtComp_AdjFlow::CAvgGradArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Incompressible flow, primitive variables nDim+1, (P, vx, vy, vz) ---*/
  Mean_GradPsiVar = new su2double* [nVar];
  
  /*--- Incompressible flow, gradient primitive variables nDim+1, (P, vx, vy, vz) ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradPsiVar[iVar] = new su2double [nDim];
  
}

CAvgGradArtComp_AdjFlow::~CAvgGradArtComp_AdjFlow(void) {
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradPsiVar[iVar];
  delete [] Mean_GradPsiVar;
  
}

void CAvgGradArtComp_AdjFlow::ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                                          su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                                          su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) {
  
  /*--- Normalized normal vector ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
  Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  
  /*--- Mean Viscosities ---*/
  
  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  
  /*--- Mean gradient approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (iDim = 0; iDim < nDim; iDim++)
      Mean_GradPsiVar[iVar][iDim] = 0.5*(PsiVar_Grad_i[iVar][iDim] + PsiVar_Grad_j[iVar][iDim]);
  
  /*--- Get projected flux tensor ---*/
  
  GetViscousArtCompProjFlux(Mean_GradPsiVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  
  /*--- Update viscous residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual_i[iVar] = Proj_Flux_Tensor[iVar];
    val_residual_j[iVar] = Proj_Flux_Tensor[iVar];
  }
  
  /*--- Implicit part ---*/
  
  if (implicit) {
    
    dist_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
    dist_ij = sqrt(dist_ij);
    
    if (dist_ij == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_ii[iVar][jVar] = 0.0;
          val_Jacobian_jj[iVar][jVar] = 0.0;
          val_Jacobian_ij[iVar][jVar] = 0.0;
          val_Jacobian_ji[iVar][jVar] = 0.0;
        }
      }
    }
    else {
      GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, dist_ij, UnitNormal,
                                Area, val_Jacobian_ii, val_Jacobian_jj);
      GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, dist_ij, UnitNormal,
                                Area, val_Jacobian_ji, val_Jacobian_ij);
    }
    
  }
  
}

CAvgGradCorrectedArtComp_AdjFlow::CAvgGradCorrectedArtComp_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  PsiVar_i = new su2double [nVar];
  PsiVar_j = new su2double [nVar];
  Proj_Mean_GradPsiVar_Edge = new su2double [nVar];
  Edge_Vector = new su2double [nDim];
  
  Mean_GradPsiVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradPsiVar[iVar] = new su2double [nDim];
  
}

CAvgGradCorrectedArtComp_AdjFlow::~CAvgGradCorrectedArtComp_AdjFlow(void) {
  
  delete [] PsiVar_i;
  delete [] PsiVar_j;
  delete [] Proj_Mean_GradPsiVar_Edge;
  delete [] Edge_Vector;
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradPsiVar[iVar];
  delete [] Mean_GradPsiVar;
  
}

void CAvgGradCorrectedArtComp_AdjFlow::ComputeResidual(su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) {
  
  /*--- Normalized normal vector ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Conversion to Primitive Variables ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    PsiVar_i[iVar] = Psi_i[iVar];
    PsiVar_j[iVar] = Psi_j[iVar];
  }
  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
  Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  
  /*--- Mean Viscosities ---*/
  
  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }
  
  /*--- Projection of the mean gradient in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradPsiVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPsiVar[iVar][iDim] = 0.5*(PsiVar_Grad_i[iVar][iDim] + PsiVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradPsiVar_Edge[iVar] += Mean_GradPsiVar[iVar][iDim]*Edge_Vector[iDim];
    }
    if (dist_ij_2 != 0.0) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Mean_GradPsiVar[iVar][iDim] -= (Proj_Mean_GradPsiVar_Edge[iVar] -
                                        (PsiVar_j[iVar]-PsiVar_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
      }
    }
  }
  
  /*--- Get projected flux tensor ---*/
  
  GetViscousArtCompProjFlux(Mean_GradPsiVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  
  /*--- Update viscous residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual_i[iVar] = Proj_Flux_Tensor[iVar];
    val_residual_j[iVar] = Proj_Flux_Tensor[iVar];
  }
  
  /*--- Implicit part ---*/
  
  if (implicit) {
    
    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_ii[iVar][jVar] = 0.0;
          val_Jacobian_jj[iVar][jVar] = 0.0;
          val_Jacobian_ij[iVar][jVar] = 0.0;
          val_Jacobian_ji[iVar][jVar] = 0.0;
        }
      }
    }
    else {
      GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, sqrt(dist_ij_2), UnitNormal,
                                Area, val_Jacobian_ii, val_Jacobian_jj);
      GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, sqrt(dist_ij_2), UnitNormal,
                                Area, val_Jacobian_ji, val_Jacobian_ij);
    }
    
  }
  
}
