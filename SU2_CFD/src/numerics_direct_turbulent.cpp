/*!
 * \file numerics_direct_turbulent.cpp
 * \brief This file contains all the convective term discretization.
 * \author F. Palacios, A. Bueno
 * \version 3.2.5 "eagle"
 *
 * Copyright (C) 2012-2014 SU2 <https://github.com/su2code>.
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
#include "../include/numerics_machine_learning_turbulent.hpp"
#include <limits>

CUpwSca_TurbSA::CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                               CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();
  
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  
}

CUpwSca_TurbSA::~CUpwSca_TurbSA(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  
}

void CUpwSca_TurbSA::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  q_ij = 0.0;
  
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  } else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  val_residual[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_j[0][0] = a1;
  }
  
}

CAvgGrad_TurbSA::CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  
  Edge_Vector = new double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new double [nVar];
  Proj_Mean_GradTurbVar_Edge = new double [nVar];
  Mean_GradTurbVar = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new double [nDim];
  
}

CAvgGrad_TurbSA::~CAvgGrad_TurbSA(void) {
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGrad_TurbSA::ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {
  
  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+nu_e*proj_vector_ij)/sigma;
  }
  
}

CAvgGradCorrected_TurbSA::CAvgGradCorrected_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iVar;
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  
  Edge_Vector = new double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new double [nVar];
  Proj_Mean_GradTurbVar_Edge = new double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new double [nVar];
  Mean_GradTurbVar = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new double [nDim];
  
}

CAvgGradCorrected_TurbSA::~CAvgGradCorrected_TurbSA(void) {
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGradCorrected_TurbSA::ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {
  
  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
  }
  
}

CSourcePieceWise_TurbSA::CSourcePieceWise_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  rotating_frame = config->GetRotating_Frame();
  
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
  
  /*--- Create values for interfacing with the functions ---*/
  
  SAInputs = new SpalartAllmarasInputs(nDim);
  SAConstants = new SpalartAllmarasConstants;
  
  nResidual = 4; nJacobian = 1;
  testResidual = new double[nResidual];
  testJacobian = new double[nJacobian];
  DUiDXj = new double*[nDim];
  for(int i=0; i < nDim; i++)
    DUiDXj[i] = new double[nDim];
  DNuhatDXj = new double[nDim];
  
}

CSourcePieceWise_TurbSA::~CSourcePieceWise_TurbSA(void) {
  
  delete SAInputs;
  delete SAConstants;
  delete testResidual;
  delete testJacobian;
  for (int i=0; i < nDim; i++)
    delete DUiDXj[i];
  delete DUiDXj;
  delete DNuhatDXj;
  
}

void CSourcePieceWise_TurbSA::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  if (incompressible) {
    Density_i = V_i[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
  }
  
  val_residual[0] = 0.0;
  Production      = 0.0;
  Destruction     = 0.0;
  CrossProduction = 0.0;
  val_Jacobian_i[0][0] = 0.0;
  
  /*--- Computation of vorticity ---*/
  
  Vorticity = (PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1])*(PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1]);
  if (nDim == 3) Vorticity += ( (PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2])*(PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2])
                               + (PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) );
  Omega = sqrt(Vorticity);
  
  /*--- Rotational correction term ---*/
  
  if (rotating_frame) {
    div = PrimVar_Grad_i[1][0] + PrimVar_Grad_i[2][1];
    if (nDim == 3) div += PrimVar_Grad_i[3][2];
    StrainMag = 0.0;
    StrainMag += pow(PrimVar_Grad_i[1][0] - 1.0/3.0*div,2.0);
    StrainMag += pow(PrimVar_Grad_i[2][1] - 1.0/3.0*div,2.0);
    if (nDim == 3) StrainMag += pow(PrimVar_Grad_i[3][2] - 1.0/3.0*div,2.0);
    StrainMag += 2.0*pow(0.5*(PrimVar_Grad_i[1][1]+PrimVar_Grad_i[2][0]),2.0);
    if (nDim == 3) {
      StrainMag += 2.0*pow(0.5*(PrimVar_Grad_i[1][2]+PrimVar_Grad_i[3][0]),2.0);
      StrainMag += 2.0*pow(0.5*(PrimVar_Grad_i[2][2]+PrimVar_Grad_i[3][1]),2.0);
    }
    StrainMag = sqrt(2.0*StrainMag);
    Omega += 2.0*min(0.0, StrainMag-Omega);
  }
  
  if (dist_i > 1e-10) {
    
    /*--- Production term ---*/
    
    dist_i_2 = dist_i*dist_i;
    nu = Laminar_Viscosity_i/Density_i;
    Ji = TurbVar_i[0]/nu;
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    ft2 = ct3*exp(-ct4*Ji_2);
    S = Omega;
    inv_k2_d2 = 1.0/(k2*dist_i_2);
    
    Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
    Shat = max(Shat, 1.0e-10);
    inv_Shat = 1.0/Shat;
    
    /*--- Production term ---*/;

//    Original SA model
//    Production = cb1*(1.0-ft2)*Shat*TurbVar_i[0]*Volume;
    
    Production = cb1*Shat*TurbVar_i[0]*Volume;

    /*--- Destruction term ---*/
    
    r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
    g = r + cw2*(pow(r,6.0)-r);
    g_6 =	pow(g,6.0);
    glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
    fw = g*glim;
    
//    Original SA model
//    Destruction = (cw1*fw-cb1*ft2/k2)*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
    
    Destruction = cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;

    /*--- Diffusion term ---*/
    
    norm2_Grad = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];
    CrossProduction = cb2_sigma*norm2_Grad*Volume;
    
    val_residual[0] = Production - Destruction + CrossProduction;
    
    /*--- Implicit part, production term ---*/
    
    dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
    dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
    if ( Shat <= 1.0e-10 ) dShat = 0.0;
    else dShat = (fv2+TurbVar_i[0]*dfv2)*inv_k2_d2;
    val_Jacobian_i[0][0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;
    
    /*--- mplicit part, destruction term ---*/
    
    dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
    if (r == 10.0) dr = 0.0;
    dg = dr*(1.+cw2*(6.0*pow(r,5.0)-1.0));
    dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
    val_Jacobian_i[0][0] -= cw1*(dfw*TurbVar_i[0] +	2.0*fw)*TurbVar_i[0]/dist_i_2*Volume;
  }
  
}

CUpwSca_TurbSST::CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();
  
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  
}

CUpwSca_TurbSST::~CUpwSca_TurbSST(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  
}

void CUpwSca_TurbSST::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  if (incompressible) {
    Density_i = V_i[nDim+1];
    Density_j = V_j[nDim+1];
  }
  else {
    Density_i = V_i[nDim+2];
    Density_j = V_j[nDim+2];
  }
  
  q_ij = 0.0;
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  
  val_residual[0] = a0*Density_i*TurbVar_i[0]+a1*Density_j*TurbVar_j[0];
  val_residual[1] = a0*Density_i*TurbVar_i[1]+a1*Density_j*TurbVar_j[1];
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;		val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[1][0] = 0.0;		val_Jacobian_i[1][1] = a0;
    
    val_Jacobian_j[0][0] = a1;		val_Jacobian_j[0][1] = 0.0;
    val_Jacobian_j[1][0] = 0.0;		val_Jacobian_j[1][1] = a1;
  }
  
}

CAvgGrad_TurbSST::CAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar, double *constants, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma_k1  = constants[0];
  sigma_om1 = constants[2];
  sigma_k2  = constants[1];
  sigma_om2 = constants[3];
  
  Edge_Vector = new double [nDim];
  Proj_Mean_GradTurbVar_Normal = new double [nVar];
  Proj_Mean_GradTurbVar_Edge = new double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new double [nVar];
  Mean_GradTurbVar = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new double [nDim];
  
}

CAvgGrad_TurbSST::~CAvgGrad_TurbSST(void) {
  
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Normal;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGrad_TurbSST::ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {
  
  double sigma_kine_i, sigma_kine_j, sigma_omega_i, sigma_omega_j;
  double diff_i_kine, diff_i_omega, diff_j_kine, diff_j_omega;
  
  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute the blended constant for the viscous terms ---*/
  sigma_kine_i  = F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2;
  sigma_kine_j  = F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2;
  sigma_omega_i = F1_i*sigma_om1 + (1.0 - F1_i)*sigma_om2;
  sigma_omega_j = F1_j*sigma_om1 + (1.0 - F1_j)*sigma_om2;
  
  /*--- Compute mean effective viscosity ---*/
  diff_i_kine  = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine  = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_omega = Laminar_Viscosity_i + sigma_omega_i*Eddy_Viscosity_i;
  diff_j_omega = Laminar_Viscosity_j + sigma_omega_j*Eddy_Viscosity_j;
  
  diff_kine  = 0.5*(diff_i_kine + diff_j_kine);    // Could instead use weighted average!
  diff_omega = 0.5*(diff_i_omega + diff_j_omega);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
  }
  
  val_residual[0] = diff_kine*Proj_Mean_GradTurbVar_Corrected[0];
  val_residual[1] = diff_omega*Proj_Mean_GradTurbVar_Corrected[1];
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {
    Jacobian_i[0][0] = -diff_kine*proj_vector_ij/Density_i;		Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;									    Jacobian_i[1][1] = -diff_omega*proj_vector_ij/Density_i;
    
    Jacobian_j[0][0] = diff_kine*proj_vector_ij/Density_j; 		Jacobian_j[0][1] = 0.0;
    Jacobian_j[1][0] = 0.0;									    Jacobian_j[1][1] = diff_omega*proj_vector_ij/Density_j;
  }
  
}


CAvgGradCorrected_TurbSST::CAvgGradCorrected_TurbSST(unsigned short val_nDim, unsigned short val_nVar, double *constants, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma_k1  = constants[0];
  sigma_om1 = constants[2];
  sigma_k2  = constants[1];
  sigma_om2 = constants[3];
  
  Edge_Vector = new double [nDim];
  Proj_Mean_GradTurbVar_Normal = new double [nVar];
  Proj_Mean_GradTurbVar_Edge = new double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new double [nVar];
  Mean_GradTurbVar = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new double [nDim];
  
}

CAvgGradCorrected_TurbSST::~CAvgGradCorrected_TurbSST(void) {
  
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Normal;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGradCorrected_TurbSST::ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {
  
  double sigma_kine_i, sigma_kine_j, sigma_omega_i, sigma_omega_j;
  double diff_i_kine, diff_i_omega, diff_j_kine, diff_j_omega;
  
  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute the blended constant for the viscous terms ---*/
  sigma_kine_i  = F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2;
  sigma_kine_j  = F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2;
  sigma_omega_i = F1_i*sigma_om1 + (1.0 - F1_i)*sigma_om2;
  sigma_omega_j = F1_j*sigma_om1 + (1.0 - F1_j)*sigma_om2;
  
  /*--- Compute mean effective viscosity ---*/
  diff_i_kine  = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine  = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_omega = Laminar_Viscosity_i + sigma_omega_i*Eddy_Viscosity_i;
  diff_j_omega = Laminar_Viscosity_j + sigma_omega_j*Eddy_Viscosity_j;
  
  diff_kine  = 0.5*(diff_i_kine + diff_j_kine);    // Could instead use weighted average!
  diff_omega = 0.5*(diff_i_omega + diff_j_omega);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  dist_ij_2 = 0.0; proj_vector_ij = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  }
  
  val_residual[0] = diff_kine*Proj_Mean_GradTurbVar_Corrected[0];
  val_residual[1] = diff_omega*Proj_Mean_GradTurbVar_Corrected[1];
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {
    Jacobian_i[0][0] = -diff_kine*proj_vector_ij/Density_i;		Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;									    Jacobian_i[1][1] = -diff_omega*proj_vector_ij/Density_i;
    
    Jacobian_j[0][0] = diff_kine*proj_vector_ij/Density_j; 		Jacobian_j[0][1] = 0.0;
    Jacobian_j[1][0] = 0.0;									    Jacobian_j[1][1] = diff_omega*proj_vector_ij/Density_j;
  }
  
}

CSourcePieceWise_TurbSST::CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar, double *constants,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  /*--- Closure constants ---*/
  beta_star     = constants[6];
  sigma_omega_1 = constants[2];
  sigma_omega_2 = constants[3];
  beta_1        = constants[4];
  beta_2        = constants[5];
  alfa_1        = constants[8];
  alfa_2        = constants[9];
  a1            = constants[7];
}

CSourcePieceWise_TurbSST::~CSourcePieceWise_TurbSST(void) { }

void CSourcePieceWise_TurbSST::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  unsigned short iDim;
  double alfa_blended, beta_blended;
  double diverg, pk, pw, zeta;
  
  if (incompressible) {
    Density_i = V_i[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];
  }
  
  val_residual[0] = 0.0;        val_residual[1] = 0.0;
  val_Jacobian_i[0][0] = 0.0;		val_Jacobian_i[0][1] = 0.0;
  val_Jacobian_i[1][0] = 0.0;		val_Jacobian_i[1][1] = 0.0;
  
  /*--- Computation of blended constants for the source terms---*/
  
  alfa_blended = F1_i*alfa_1 + (1.0 - F1_i)*alfa_2;
  beta_blended = F1_i*beta_1 + (1.0 - F1_i)*beta_2;
  
  if (dist_i > 1e-10) {
    
    /*--- Production ---*/
    
    diverg = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      diverg += PrimVar_Grad_i[iDim+1][iDim];
    
    pk = Eddy_Viscosity_i*StrainMag*StrainMag - 2.0/3.0*Density_i*TurbVar_i[0]*diverg;
    pk = min(pk,20.0*beta_star*Density_i*TurbVar_i[1]*TurbVar_i[0]);
    pk = max(pk,0.0);
    
    zeta = max(TurbVar_i[1],StrainMag*F2_i/a1);
    pw = StrainMag*StrainMag - 2.0/3.0*zeta*diverg;
    pw = max(pw,0.0);
    
    val_residual[0] += pk*Volume;
    val_residual[1] += alfa_blended*Density_i*pw*Volume;
    
    /*--- Dissipation ---*/
    
    val_residual[0] -= beta_star*Density_i*TurbVar_i[1]*TurbVar_i[0]*Volume;
    val_residual[1] -= beta_blended*Density_i*TurbVar_i[1]*TurbVar_i[1]*Volume;
    
    /*--- Cross diffusion ---*/
    
    val_residual[1] += (1.0 - F1_i)*CDkw*Volume;
    
    /*--- Implicit part ---*/
    
    val_Jacobian_i[0][0] = -beta_star*TurbVar_i[1]*Volume;		val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[1][0] = 0.0;                               val_Jacobian_i[1][1] = -2.0*beta_blended*TurbVar_i[1]*Volume;
  }
  
}

CUpwSca_TurbML::CUpwSca_TurbML(unsigned short val_nDim, unsigned short val_nVar,
                               CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();
  
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  
}

CUpwSca_TurbML::~CUpwSca_TurbML(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  
}

void CUpwSca_TurbML::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  
  q_ij = 0.0;
  
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  } else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  val_residual[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_j[0][0] = a1;
  }
  
  
}

CAvgGrad_TurbML::CAvgGrad_TurbML(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  
  Edge_Vector = new double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new double [nVar];
  Proj_Mean_GradTurbVar_Edge = new double [nVar];
  Mean_GradTurbVar = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new double [nDim];
  
}

CAvgGrad_TurbML::~CAvgGrad_TurbML(void) {
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGrad_TurbML::ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {
  
  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+nu_e*proj_vector_ij)/sigma;
  }
  
}

CAvgGradCorrected_TurbML::CAvgGradCorrected_TurbML(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iVar;
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  
  Edge_Vector = new double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new double [nVar];
  Proj_Mean_GradTurbVar_Edge = new double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new double [nVar];
  Mean_GradTurbVar = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new double [nDim];
  
}

CAvgGradCorrected_TurbML::~CAvgGradCorrected_TurbML(void) {
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGradCorrected_TurbML::ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {
  
  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
  }
  
}

CSourcePieceWise_TurbML::CSourcePieceWise_TurbML(unsigned short val_nDim, unsigned short val_nVar,
                                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  double *uinf = config->GetVelocity_FreeStreamND();
  for (unsigned short i = 0; i < nDim; i++){
    uInfinity += uinf[i] * uinf[i];
  }
  uInfinity = sqrt(uInfinity);
  
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  //transition     = (config->GetKind_Trans_Model() == LM);
  transition = false; // Debugging, -AA
  rotating_frame = config->GetRotating_Frame();
  
  /* Create values for interfacing with the functions */
  SAInputs = new SpalartAllmarasInputs(nDim);
  SAConstants = new SpalartAllmarasConstants;
  SAOtherOutputs = new SpalartAllmarasOtherOutputs;
  
  SANondimInputs = new CSANondimInputs(val_nDim);
  
  nResidual = 4;
  nJacobian = 1;
  
  SAResidual = new double[nResidual];
  SANondimResidual = new double[nResidual];
  Residual = new double[nResidual];
  NondimResidual = new double[nResidual];
  ResidualDiff = new double[nResidual];
  NondimResidualDiff = new double[nResidual];
  
  SAJacobian = new double[nJacobian];
  
  //testResidual = new double[nResidual];
  //testJacobian = new double[nJacobian];
  DUiDXj = new double*[nDim];
  for(int i=0; i < nDim; i++){
    DUiDXj[i] = new double[nDim];
  }
  DNuhatDXj = new double[nDim];
  
  // Construct the nnet
  string readFile = config->GetML_Turb_Model_File();
//  string checkFile = config->GetML_Turb_Model_Check_File();
//  cout << "Loading ML file from " << readFile << endl;
//  CNeurNet* Net = new CNeurNet(readFile, checkFile);
  
  this->featureset = config->GetML_Turb_Model_FeatureSet();
  
  cout << "in constructor, featureset is " << featureset << endl;
  
  CScalePredictor* Pred = new CScalePredictor(readFile);
  this->MLModel = Pred;
  cout << "ML File successfully read " << endl;
}

CSourcePieceWise_TurbML::~CSourcePieceWise_TurbML(void) {
  
  delete MLModel;
  delete SAInputs;
  delete SAConstants;
  
  delete SAResidual;
  delete SANondimResidual;
  delete Residual;
  delete NondimResidual;
  delete ResidualDiff;
  delete NondimResidualDiff;
  delete SAJacobian;
//  delete testResidual;
//  delete testJacobian;
  for (int i=0; i < nDim; i++){
    delete DUiDXj[i];
  }
  delete DUiDXj;
  delete DNuhatDXj;
  
  delete SANondimInputs;
}

void CSourcePieceWise_TurbML::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  if (incompressible) {
    Density_i = V_i[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
  }
  
  /* Intialize */
  // Note that the Production, destruction, etc. are all volume independent
  
  for (int i= 0; i < nResidual; i++){
    SAResidual[i] = 0;
    SANondimResidual[i] = 0;
    Residual[i] = 0;
    NondimResidual[i] = 0;
    ResidualDiff[i] = 0;
    NondimResidualDiff[i] = 0;
  }
  
  val_residual[0] = 0.0;
  val_Jacobian_i[0][0] = 0.0;
  
  NuhatGradNorm = 0;
  for (int i =0; i < nDim; i++){
    for (int j=0; j < nDim; j++){
      DUiDXj[i][j] = PrimVar_Grad_i[i+1][j];
    }
    DNuhatDXj[i] = TurbVar_Grad_i[0][i];
    NuhatGradNorm += TurbVar_Grad_i[0][i] * TurbVar_Grad_i[0][i];
  }
  
  /* Call Spalart-Allmaras (for comparison) */
  SAInputs->Set(DUiDXj, DNuhatDXj, rotating_frame, transition, dist_i, Laminar_Viscosity_i, Density_i, TurbVar_i[0], intermittency);
  
  SpalartAllmarasSourceTerm(SAInputs, SAConstants,SAResidual, SAJacobian, SAOtherOutputs);
  this->SANondimInputs -> Set(SAInputs);

  for (int i=0; i < nResidual; i++){
    SANondimResidual[i] = SAResidual[i];
  }
  SANondimInputs->NondimensionalizeSource(nResidual, SANondimResidual);
  
  
  fw = SAOtherOutputs->fw;
  
  // Need the individual terms of the NuHat Norm
  double dNuHatDXBar = DNuhatDXj[0] / sqrt(SANondimInputs->SourceNondim);
  double dNuHatDYBar = DNuhatDXj[1] / sqrt(SANondimInputs->SourceNondim);
  double dUDXBar = DUiDXj[0][0] / SANondimInputs->OmegaNondim;
  double dVDXBar = DUiDXj[1][0] / SANondimInputs->OmegaNondim;
  double dUDYBar = DUiDXj[0][1] / SANondimInputs->OmegaNondim;
  double dVDYBar = DUiDXj[1][1] / SANondimInputs->OmegaNondim;
  double Turbulent_Kinematic_Viscosity = TurbVar_i[0];
  
  double nuRef = 1.0;
  double nu = Laminar_Viscosity_i / Density_i;
  double nuscale = nu / nuRef;
  double distalt = dist_i;
  double nuhatalt = Turbulent_Kinematic_Viscosity / nuscale;
  double omega = SANondimInputs->OmegaBar * SANondimInputs->OmegaNondim;
  double omegaalt = omega / nuscale;
  double nuhatgradmagalt = SANondimInputs->NuHatGradNorm / (nuscale * nuscale);
  double omeganondimeralt = (1/distalt) * (nuhatalt / distalt);
  double sourcenondimeralt = (nuhatalt / distalt) * (nuhatalt / distalt);
  double nondim_nuhatgradmagalt = nuhatgradmagalt / sourcenondimeralt;
  double nondimOmegaAlt = omegaalt / omeganondimeralt;
  
//  double Laminar_Kinematic_Viscosity = Laminar_Viscosity_i / Density_i;
  
  int nInputMLVariables = 0;
  int nOutputMLVariables = 0;
  double* netInput = NULL;
  double* netOutput = NULL;
  
  if (featureset.compare("SA") == 0){
    // Set the output equal to the spalart allmaras output.
    for (int i = 0; i < nResidual; i++){
      Residual[i] = SAResidual[i];
      NondimResidual[i] = SANondimResidual[i];
    }
  }else if (featureset.compare("nondim_production")==0){
    nInputMLVariables = 2;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    netInput[0] = SANondimInputs->Chi;
    netInput[1] = SANondimInputs->OmegaBar;
    
    // Predict using Nnet
    MLModel->Predict(netInput, netOutput);

    // Gather all the appropriate variables
    NondimResidual[0] = netOutput[0];
    NondimResidual[1] = SANondimResidual[1];
    NondimResidual[2] = SANondimResidual[2];
    NondimResidual[3] = NondimResidual[0] - NondimResidual[1] + NondimResidual[2];
    
    for (int i=0; i < nResidual; i++){
      Residual[i] = NondimResidual[i];
      //cout << "NondimResidual " << i <<" "<< NondimResidual[i] << endl;
    }
    SANondimInputs->DimensionalizeSource(nResidual, Residual);
    /*
    for (int i=0; i < nResidual; i++){
      cout << "DimResidual " << i << " " << Residual[i] << endl;
    }
     */
  }else if(featureset.compare("nondim_production_log") == 0){
    nInputMLVariables = 2;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    netInput[0] = log10(SANondimInputs->Chi);
    netInput[1] = log10(SANondimInputs->OmegaBar);
    
    // Predict using Nnet
    MLModel->Predict(netInput, netOutput);
    
    // Gather all the appropriate variables
    NondimResidual[0] = netOutput[0];
    NondimResidual[1] = SANondimResidual[1];
    NondimResidual[2] = SANondimResidual[2];
    NondimResidual[3] = NondimResidual[0] - NondimResidual[1] + NondimResidual[2];
    
    for (int i=0; i < nResidual; i++){
      Residual[i] = NondimResidual[i];
//      cout << "NondimResidual " << i << NondimResidual[i] << endl;
    }
    
    SANondimInputs->DimensionalizeSource(nResidual, Residual);
  /*
    for (int i=0; i < nResidual; i++){
      cout << "DimResidual " << i << Residual[i] << endl;
    }
   */
  }else if(featureset.compare("nondim_production_logchi") == 0){
    nInputMLVariables = 2;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    netInput[0] = log10(SANondimInputs->Chi);
    netInput[1] = SANondimInputs->OmegaBar;
    
    // Predict using Nnet
    MLModel->Predict(netInput, netOutput);
    
    // Gather all the appropriate variables
    NondimResidual[0] = netOutput[0];
    NondimResidual[1] = SANondimResidual[1];
    NondimResidual[2] = SANondimResidual[2];
    NondimResidual[3] = NondimResidual[0] - NondimResidual[1] + NondimResidual[2];
    
    for (int i=0; i < nResidual; i++){
      Residual[i] = NondimResidual[i];
    }
    SANondimInputs->DimensionalizeSource(nResidual, Residual);
  }else if(featureset.compare("production")==0){
//    cout <<"In production" << endl;
    nInputMLVariables = 3;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    netInput[0] = SANondimInputs->SourceNondim;
    netInput[1] = SANondimInputs->Chi;
    netInput[2] = SANondimInputs->OmegaBar;
    
//    cout << "Net inputs ";
//    for (int i = 0; i < 3; i++){
//      cout << "\t" << netInput[i];
//    }
//    cout << endl;
    
    // Predict using Nnet
    MLModel->Predict(netInput, netOutput);
    
    // Gather the appropriate values
    Residual[0] = netOutput[0];
    Residual[1] = SAResidual[1];
    Residual[2] = SAResidual[2];
    Residual[3] = Residual[0] - Residual[1] + Residual[2];

//    cout << "ML Production " << Residual[0] << endl;
//    cout << "SA Production " << SAResidual[0] << endl;
    
    for (int i=0; i < nResidual; i++){
      NondimResidual[i] = Residual[i];
    }
    SANondimInputs->NondimensionalizeSource(nResidual, NondimResidual);
    
  }else if (featureset.compare("nondim_destruction")==0){
    nInputMLVariables = 2;
    nOutputMLVariables = 1;
    
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    netInput[0] = SANondimInputs->Chi;
    netInput[1] = SANondimInputs->OmegaBar;
    
    MLModel->Predict(netInput, netOutput);
    
    NondimResidual[0] = SANondimResidual[0];
    NondimResidual[1] = netOutput[0];
    NondimResidual[2] = SANondimResidual[2];
    NondimResidual[3] = NondimResidual[0] - NondimResidual[1] + NondimResidual[2];
    
    for (int i=0; i < nResidual; i++){
      Residual[i] = NondimResidual[i];
    }
    SANondimInputs->DimensionalizeSource(nResidual, Residual);
  }else if(featureset.compare("destruction")==0){
      nInputMLVariables = 3;
      nOutputMLVariables = 1;
      netInput = new double[nInputMLVariables];
      netOutput = new double[nOutputMLVariables];
    
      netInput[0] = SANondimInputs->SourceNondim;
      netInput[1] = SANondimInputs->Chi;
      netInput[2] = SANondimInputs->OmegaBar;
      
      // Predict using Nnet
      MLModel->Predict(netInput, netOutput);
      
      // Gather the appropriate values
      Residual[0] = SAResidual[0];
      Residual[1] = netOutput[0];
      Residual[2] = SAResidual[2];
      Residual[3] = Residual[0] - Residual[1] + Residual[2];
      
      for (int i=0; i < nResidual; i++){
        NondimResidual[i] = Residual[i];
      }
      SANondimInputs->NondimensionalizeSource(nResidual, NondimResidual);
  }else if (featureset.compare("nondim_crossproduction")==0){
    nInputMLVariables = 2;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    netInput[0] = SANondimInputs->Chi;
    netInput[1] = SANondimInputs->NuHatGradNormBar;
    
    //cout << "IN nondim_crossproduction " << endl;
    //  cout << "chi " << netInput[0] <<  endl;
    //    cout << "grad norm bar "<< netInput[1] << endl;
    
    MLModel->Predict(netInput, netOutput);
    
    NondimResidual[0] = SANondimResidual[0];
    NondimResidual[1] = SANondimResidual[1];
    NondimResidual[2] = netOutput[0];
    NondimResidual[3] = NondimResidual[0] - NondimResidual[1] + NondimResidual[2];
    
    for (int i=0; i < nResidual; i++){
      Residual[i] = NondimResidual[i];
    }
    SANondimInputs->DimensionalizeSource(nResidual, Residual);
  }else if(featureset.compare("cross_production")==0){
    nInputMLVariables = 3;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    netInput[0] = SANondimInputs->SourceNondim;
    netInput[1] = SANondimInputs->Chi;
    netInput[2] = SANondimInputs->NuHatGradNormBar;
    
    // Predict using Nnet
    MLModel->Predict(netInput, netOutput);
    
    // Gather the appropriate values
    Residual[0] = SAResidual[0];
    Residual[1] = SAResidual[1];
    Residual[2] = netOutput[0];
    Residual[3] = Residual[0] - Residual[1] + Residual[2];
    
    for (int i=0; i < nResidual; i++){
      NondimResidual[i] = Residual[i];
    }
    SANondimInputs->NondimensionalizeSource(nResidual, NondimResidual);
  }else if (featureset.compare("nondim_source")==0){
    nInputMLVariables = 3;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    netInput[0] = SANondimInputs->Chi;
    netInput[1] = SANondimInputs->OmegaBar;
    netInput[2] = SANondimInputs->NuHatGradNormBar;
  
    // Predict using Nnet
    MLModel->Predict(netInput, netOutput);
    
    NondimResidual[0] = 0;
    NondimResidual[1] = 0;
    NondimResidual[2] = 0;
    NondimResidual[3] = netOutput[0];
    
    for (int i=0; i < nResidual; i++){
      Residual[i] = NondimResidual[i];
    }
    SANondimInputs->DimensionalizeSource(nResidual, Residual);
  }else if(featureset.compare("source")==0){
    nInputMLVariables =4;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    netInput[0] = SANondimInputs->SourceNondim;
    netInput[1] = SANondimInputs->Chi;
    netInput[2] = SANondimInputs->OmegaBar;
    netInput[3] = SANondimInputs->NuHatGradNormBar;
    
    
    // Predict using Nnet
    MLModel->Predict(netInput, netOutput);
    
    // Gather the appropriate values
    Residual[0] = 0;
    Residual[1] = 0;
    Residual[2] = 0;
    Residual[3] = netOutput[0];
    
    for (int i=0; i < nResidual; i++){
      NondimResidual[i] = Residual[i];
    }
    SANondimInputs->NondimensionalizeSource(nResidual, NondimResidual);
  }else if(featureset.compare("source_alt")==0){
    nInputMLVariables =4;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    netInput[0] = sourcenondimeralt;
    netInput[1] = SANondimInputs->Chi;
    netInput[2] = nondimOmegaAlt;
    netInput[3] = nondim_nuhatgradmagalt;
    
    
    // Predict using Nnet
    MLModel->Predict(netInput, netOutput);

    netOutput[0] *= nuscale * nuscale;
    
    // Gather the appropriate values
    Residual[0] = 0;
    Residual[1] = 0;
    Residual[2] = 0;
    Residual[3] = netOutput[0];
    
    for (int i=0; i < nResidual; i++){
      NondimResidual[i] = Residual[i];
    }
    SANondimInputs->NondimensionalizeSource(nResidual, NondimResidual);
  }else if(featureset.compare("source_dim_alt") == 0){
    nInputMLVariables = 4;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    
    
    netInput[0] = nuhatalt;
    netInput[1] = omegaalt;
    netInput[2] = nuhatgradmagalt;
    netInput[3] = dist_i;
    
    // Predict using Nnet
    MLModel->Predict(netInput, netOutput);
    
    // Need to scale the output back
    netOutput[0] *= nuscale * nuscale;
    
    // Gather the appropriate values
    Residual[0] = 0;
    Residual[1] = 0;
    Residual[2] = 0;
    Residual[3] = netOutput[0];
    
    for (int i=0; i < nResidual; i++){
      NondimResidual[i] = Residual[i];
    }
    SANondimInputs->NondimensionalizeSource(nResidual, NondimResidual);
    
  }else if(featureset.compare("source_all")==0){
    nInputMLVariables = 8;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    netInput[0] = SANondimInputs->SourceNondim;
    netInput[1] = SANondimInputs->Chi;
    netInput[2] = dNuHatDXBar;
    netInput[3] = dNuHatDYBar;
    netInput[4] = dUDXBar;
    netInput[5] = dUDYBar;
    netInput[6] = dVDXBar;
    netInput[7] = dVDYBar;
    
    // Predict using Nnet
    MLModel->Predict(netInput, netOutput);
    
    // Gather the appropriate values
    Residual[0] = 0;
    Residual[1] = 0;
    Residual[2] = 0;
    Residual[3] = netOutput[0];
    
    for (int i=0; i < nResidual; i++){
      NondimResidual[i] = Residual[i];
    }
    SANondimInputs->NondimensionalizeSource(nResidual, NondimResidual);
    
  }else if (featureset.compare("fw_hifi")==0){
    throw("doesn't work");
    nInputMLVariables = 2;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    double chi = SANondimInputs->Chi;
    double omegaBar = SANondimInputs->OmegaBar;
    // Karthik nondimensionalizes by d / vhat whereas I do by /(v + vhat)
    omegaBar *= 1 + 1/chi;
    
    netInput[0] = chi;
    netInput[1] = omegaBar;
    
    MLModel->Predict(netInput, netOutput);
    
    double safw = SAOtherOutputs->fw;
    double newfw = netOutput[0];
    if (newfw < -1){
      newfw = 1;
    }
    if (newfw > 6){
      newfw = 6;
    }
    // The output is the value of fw. Need to replace the destruction term with the new computation
    double turbKinVisc = SAInputs->Turbulent_Kinematic_Viscosity;
    double dist2 = SAInputs->dist * SAInputs->dist;
    double newdestruction = SAConstants->cw1 * (newfw +safw) * turbKinVisc * turbKinVisc / dist2;
    
    for (int i= 0; i < nResidual; i++){
      Residual[i] = SAResidual[i];
    }
    Residual[1] = newdestruction;
    Residual[3] = Residual[0] - Residual[1] + Residual[2];
    
    for (int i= 0; i < nResidual; i++){
      NondimResidual[i] = Residual[i];
    }
    SANondimInputs->NondimensionalizeSource(nResidual, NondimResidual);
    
  }else if (featureset.compare("fw_hifi_2")==0){
    nInputMLVariables = 2;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    
    double chi = SANondimInputs->Chi;
    double omegaBar = SANondimInputs->OmegaBar;
    // Karthik nondimensionalizes by d / vhat whereas I do by /(v + vhat)
    omegaBar *= 1 + 1/chi;
    
    netInput[0] = chi;
    netInput[1] = omegaBar;
    
    MLModel->Predict(netInput, netOutput);
    
    double safw = SAOtherOutputs->fw;
    double newfw = netOutput[0];
    // The output is the value of fw. Need to replace the destruction term with the new computation
    double turbKinVisc = SAInputs->Turbulent_Kinematic_Viscosity;
    double dist2 = SAInputs->dist * SAInputs->dist;
    double newdestruction = SAConstants->cw1 * (newfw +safw) * turbKinVisc * turbKinVisc / dist2;
    
    for (int i= 0; i < nResidual; i++){
      Residual[i] = SAResidual[i];
    }
    Residual[1] = newdestruction;
    Residual[3] = Residual[0] - Residual[1] + Residual[2];
    
    for (int i= 0; i < nResidual; i++){
      NondimResidual[i] = Residual[i];
    }
    SANondimInputs->NondimensionalizeSource(nResidual, NondimResidual);
    
  }else if(featureset.compare("fw") == 0){
    nInputMLVariables = 2;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    double chi = SANondimInputs->Chi;
    double omegaBar = SANondimInputs->OmegaBar;
    netInput[0] = chi;
    netInput[1] = omegaBar;
    MLModel->Predict(netInput, netOutput);
    
    // The output is fw. Replicate the destruction term.
    double fw_ml = netOutput[0];
    double mul_dest = SAConstants->cw1 * fw_ml;
    Residual[0] = SAResidual[0];
    Residual[1] = mul_dest * Turbulent_Kinematic_Viscosity * Turbulent_Kinematic_Viscosity / (dist_i * dist_i);
    Residual[2] = SAResidual[2];
    Residual[3] = Residual[0] - Residual[1] + Residual[2];
    for (int i= 0; i < nResidual; i++){
      NondimResidual[i] = Residual[i];
    }
    SANondimInputs->NondimensionalizeSource(nResidual, NondimResidual);
  }else if(featureset.compare("mul_destruction") == 0){
    nInputMLVariables = 2;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    double chi = SANondimInputs->Chi;
    double omegaBar = SANondimInputs->OmegaBar;
    netInput[0] = chi;
    netInput[1] = omegaBar;
    MLModel->Predict(netInput, netOutput);
    
    // The output is a multiplier to the destruction term. Replicate the
    // destruction term
    double mul_dest = netOutput[0];
    Residual[0] = SAResidual[0];
    Residual[1] = mul_dest * Turbulent_Kinematic_Viscosity * Turbulent_Kinematic_Viscosity / (dist_i * dist_i);
    Residual[2] = SAResidual[2];
    Residual[3] = Residual[0] - Residual[1] + Residual[2];
    for (int i= 0; i < nResidual; i++){
      NondimResidual[i] = Residual[i];
    }
    SANondimInputs->NondimensionalizeSource(nResidual, NondimResidual);
  }else if(featureset.compare("mul_production")==0){
    nInputMLVariables = 2;
    nOutputMLVariables = 1;
    netInput = new double[nInputMLVariables];
    netOutput = new double[nOutputMLVariables];
    double chi = SANondimInputs->Chi;
    double omegaBar = SANondimInputs->OmegaBar;
    netInput[0] = chi;
    netInput[1] = omegaBar;
    MLModel->Predict(netInput, netOutput);
    // The output is a multiplier to the destruction term. Replicate the
    // production term
    double mul_prod = netOutput[0];
    Residual[0] = mul_prod * Turbulent_Kinematic_Viscosity * SAOtherOutputs->Omega;
    Residual[1] = SAResidual[1];
    Residual[2] = SAResidual[2];
    Residual[3] = Residual[0] - Residual[1] + Residual[2];
    for (int i= 0; i < nResidual; i++){
      NondimResidual[i] = Residual[i];
    }
    SANondimInputs->NondimensionalizeSource(nResidual, NondimResidual);
    
  }else{
    cout << "None of the conditions met" << endl;
    cout << "featureset is " << featureset << endl;
    throw "ML_Turb_Model_Nondimensionalization not recognized";
  }
  delete [] netInput;
  delete [] netOutput;
  
  // Hack if the wall distance is too low
  if (dist_i < 1e-6){
    for (int i= 0; i < nResidual; i++){
      Residual[i] = 0;
      NondimResidual[i] = 0;
      SAResidual[i] = 0;
      SANondimResidual[i] = 0;
    }
  }
  
  // Compute Shivaji Medida's BL vs. Wake equation
  double strainRateMag = 0;
  for (int i= 0; i < nDim; i++){
    for (int j = 0; j < nDim; j++){
      double sij = 0.5 * (DUiDXj[i][j] + DUiDXj[j][i]);
      strainRateMag += 2 * (sij * sij);
    }
  }
  
  strainRateMag = sqrt(strainRateMag);
  double ReS = Density_i * strainRateMag * dist_i * dist_i / (0.09 * Laminar_Viscosity_i);
  fWake = exp(- (1e-10 * ReS * ReS));
  double magU = 0;
  for (unsigned short i = 0; i < nDim; i++){
    magU += V_i[1+i] * V_i[1+i];
  }
  magU = sqrt(magU);
  isInBL = fWake > 0.5 && (magU < uInfinity * 0.99);
  
  
  // Now that we have found the ML Residual and the SA residual, see if there are
  // any special hacks that we should use
  
  unsigned short nStrings = config->GetNumML_Turb_Model_Extra();
  string *extraString = config->GetML_Turb_Model_Extra();
  
  bool hasBlOnly = false;
  for (int i= 0; i < nStrings; i++){
    if (extraString[i].compare("BlOnly") == 0){
      hasBlOnly = true;
      break;
    }
  }
  
  if (nStrings > 0){
    if (extraString[0].compare("FlatplateBlOnlyCutoff") == 0){
        // Only use ML in the boundary layer and have a sharp cutoff
      if ((Coord_i[0] < 0) || (Coord_i[1]) > 0.06 ){
        // Not in the BL, so just use the SA residual
        for (int i = 0; i < nResidual; i++){
          Residual[i] = SAResidual[i];
          NondimResidual[i] = SANondimResidual[i];
        }
      }
    }
    if (hasBlOnly){
      // Only use ML in the boundary layer (where isInBL == true)
      if (!isInBL){
        // Then use SA
        for (int i = 0; i < nResidual; i++){
          Residual[i] = SAResidual[i];
          NondimResidual[i] = SANondimResidual[i];
        }
      }
    }
  }
  
  // Compute the differences
  for (int i = 0; i < nResidual; i++){
    ResidualDiff[i] = Residual[i] - SAResidual[i];
    NondimResidualDiff[i] = NondimResidual[i] - SANondimResidual[i];
  }
  
  // Store The residual for the outer structure
  val_residual[0] = Residual[3] * Volume;
  val_Jacobian_i[0][0] = SAJacobian[0] * Volume;
  
}

int CSourcePieceWise_TurbML::NumResidual(){
  return this->nResidual;
}
