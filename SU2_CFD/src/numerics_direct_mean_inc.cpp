/*!
 * \file numerics_direct_mean_inc.cpp
 * \brief This file contains the numerical methods for incompressible flow.
 * \author F. Palacios, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
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

CUpwArtComp_Flow::CUpwArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  variable_density = (config->GetKind_DensityModel() == VARIABLE);
  energy           = config->GetEnergy_Equation();
  grid_movement    = config->GetGrid_Movement();

  Diff_V       = new su2double[nVar];
  Velocity_i   = new su2double[nDim];
  Velocity_j   = new su2double[nDim];
  MeanVelocity = new su2double[nDim];
  ProjFlux_i   = new su2double[nVar];
  ProjFlux_j   = new su2double[nVar];
  Lambda       = new su2double[nVar];
  Epsilon      = new su2double[nVar];
  Precon       = new su2double*[nVar];
  invPrecon_A  = new su2double*[nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Precon[iVar]      = new su2double[nVar];
    invPrecon_A[iVar] = new su2double[nVar];
  }
  
}

CUpwArtComp_Flow::~CUpwArtComp_Flow(void) {
  
  delete [] Diff_V;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Precon[iVar];
    delete [] invPrecon_A[iVar];
  }
  delete [] Precon;
  delete [] invPrecon_A;
  
}

void CUpwArtComp_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+8); AD::SetPreaccIn(V_j, nDim+8); AD::SetPreaccIn(Normal, nDim);

  /*--- Face area (norm or the normal vector) ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*--- Compute and unitary normal vector ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    UnitNormal[iDim] = Normal[iDim]/Area;
    if (fabs(UnitNormal[iDim]) < EPS) UnitNormal[iDim] = EPS;
  }
  
  /*--- Set primitive variables at points iPoint and jPoint ---*/
    
  Pressure_i    = V_i[0];             Pressure_j    = V_j[0];
  Temperature_i = V_i[nDim+1];        Temperature_j = V_j[nDim+1];
  DensityInc_i  = V_i[nDim+2];        DensityInc_j  = V_j[nDim+2];
  BetaInc2_i    = V_i[nDim+3];        BetaInc2_j    = V_j[nDim+3];
  Cp_i          = V_i[nDim+7];        Cp_j          = V_j[nDim+7];
  Enthalpy_i    = Cp_i*Temperature_i; Enthalpy_j    = Cp_j*Temperature_j;

  ProjVelocity = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim]    = V_i[iDim+1];
    Velocity_j[iDim]    = V_j[iDim+1];
    MeanVelocity[iDim]  = 0.5*(Velocity_i[iDim] + Velocity_j[iDim]);
    ProjVelocity       += MeanVelocity[iDim]*Normal[iDim];
  }
  
  /*--- Mean variables at points iPoint and jPoint ---*/
  
  MeanDensity     = 0.5*(DensityInc_i  + DensityInc_j);
  MeanPressure    = 0.5*(Pressure_i    + Pressure_j);
  MeanBetaInc2    = 0.5*(BetaInc2_i    + BetaInc2_j);
  MeanEnthalpy    = 0.5*(Enthalpy_i    + Enthalpy_j);
  MeanCp          = 0.5*(Cp_i          + Cp_j);
  MeanTemperature = 0.5*(Temperature_i + Temperature_j);

  /*--- Artificial sound speed based on eigs of preconditioned system ---*/

  MeanSoundSpeed = sqrt(MeanBetaInc2*Area*Area);

  /*--- We need the derivative of the equation of state to build the
   preconditioning matrix. For now, the only option is the ideal gas
   law, but in the future, dRhodT should be in the fluid model. ---*/

  MeandRhodT = 0.0; dRhodT_i = 0.0; dRhodT_j = 0.0;
  if (variable_density) {
    MeandRhodT = -MeanDensity/MeanTemperature;
    dRhodT_i   = -Density_i/Temperature_i;
    dRhodT_j   = -Density_j/Temperature_j;
  }

  /*--- Compute ProjFlux_i ---*/

  GetInviscidArtCompProjFlux(&DensityInc_i, Velocity_i, &Pressure_i, &BetaInc2_i, &Enthalpy_i, Normal, ProjFlux_i);
  
  /*--- Compute ProjFlux_j ---*/

  GetInviscidArtCompProjFlux(&DensityInc_j, Velocity_j, &Pressure_j, &BetaInc2_j, &Enthalpy_j, Normal, ProjFlux_j);

  /*--- Eigenvalues of the preconditioned system ---*/
  
  if (nDim == 2) {
    Lambda[0] = ProjVelocity;
    Lambda[1] = ProjVelocity;
    Lambda[2] = ProjVelocity - MeanSoundSpeed;
    Lambda[3] = ProjVelocity + MeanSoundSpeed;
  }
  if (nDim == 3) {
    Lambda[0] = ProjVelocity;
    Lambda[1] = ProjVelocity;
    Lambda[2] = ProjVelocity;
    Lambda[3] = ProjVelocity - MeanSoundSpeed;
    Lambda[4] = ProjVelocity + MeanSoundSpeed;
  }
  
  /*--- Absolute value of the eigenvalues ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = fabs(Lambda[iVar]);

  /*--- Build the preconditioning matrix using mean values ---*/

  GetPreconditioner(&MeanDensity, MeanVelocity, &MeanBetaInc2, &MeanCp, &MeanTemperature, &MeandRhodT, Precon);

  /*--- Build the absolute value of the preconditioned Jacobian, i.e.,
   |A_precon| = P x |Lambda| x inv(P), where P diagonalizes the matrix
   inv(Precon) x dF/dV and Lambda is the diag. matrix of its eigenvalues. ---*/

  GetPreconditionedProjJac(&MeanDensity, Lambda, &MeanBetaInc2, UnitNormal, invPrecon_A);

  /*--- Difference of primitive variables at iPoint and jPoint ---*/

  Diff_V[0] = Pressure_j - Pressure_i;
  for (iDim = 0; iDim < nDim; iDim++)
    Diff_V[iDim+1] = Velocity_j[iDim] - Velocity_i[iDim];
  Diff_V[nDim+1] = Temperature_j - Temperature_i;

  /*--- Build the inviscid Jacobian w.r.t. the primitive variables ---*/

  if (implicit) {
    GetInviscidArtCompProjJac(&DensityInc_i, Velocity_i, &BetaInc2_i, &Cp_i, &Temperature_i, &dRhodT_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidArtCompProjJac(&DensityInc_j, Velocity_j, &BetaInc2_j, &Cp_j, &Temperature_j, &dRhodT_j, Normal, 0.5, val_Jacobian_j);
  }

  /*--- Compute dissipation as Precon x |A_precon| x dV. If implicit,
   store Precon x |A_precon| from dissipation term. ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += Precon[iVar][kVar]*invPrecon_A[kVar][jVar];
      val_residual[iVar] -= 0.5*Proj_ModJac_Tensor_ij*Diff_V[jVar];
      if (implicit) {
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij;
      }
    }
  }

  if (!energy) {
    val_residual[nDim+1] = 0.0;
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        val_Jacobian_i[iVar][nDim+1] = 0.0;
        val_Jacobian_j[iVar][nDim+1] = 0.0;

        val_Jacobian_i[nDim+1][iVar] = 0.0;
        val_Jacobian_j[nDim+1][iVar] = 0.0;
      }
    }
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}

CCentJSTArtComp_Flow::CCentJSTArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  variable_density = (config->GetKind_DensityModel() == VARIABLE);
  energy           = config->GetEnergy_Equation();
  grid_movement    = config->GetGrid_Movement();

  /*--- Artifical dissipation part ---*/

  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_Flow();
  Param_Kappa_4 = config->GetKappa_4th_Flow();
  
  /*--- Allocate some structures ---*/

  Diff_V       = new su2double [nVar];
  Diff_Lapl    = new su2double [nVar];
  Velocity_i   = new su2double [nDim];
  Velocity_j   = new su2double [nDim];
  MeanVelocity = new su2double [nDim];
  ProjFlux     = new su2double [nVar];
  Precon       = new su2double*[nVar];

  for (iVar = 0; iVar < nVar; iVar++)
    Precon[iVar] = new su2double[nVar];

}

CCentJSTArtComp_Flow::~CCentJSTArtComp_Flow(void) {
  
  delete [] Diff_V;
  delete [] Diff_Lapl;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;

  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Precon[iVar];
  delete [] Precon;

}

void CCentJSTArtComp_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  /*--- Primitive variables at point i and j ---*/
  
  Pressure_i    = V_i[0];             Pressure_j    = V_j[0];
  Temperature_i = V_i[nDim+1];        Temperature_j = V_j[nDim+1];
  DensityInc_i  = V_i[nDim+2];        DensityInc_j  = V_j[nDim+2];
  BetaInc2_i    = V_i[nDim+3];        BetaInc2_j    = V_j[nDim+3];
  Cp_i          = V_i[nDim+7];        Cp_j          = V_j[nDim+7];
  Enthalpy_i    = Cp_i*Temperature_i; Enthalpy_j    = Cp_j*Temperature_j;

  Area = 0.0;
  sq_vel_i = 0.0; sq_vel_j = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim]    = V_i[iDim+1];
    Velocity_j[iDim]    = V_j[iDim+1];
    MeanVelocity[iDim]  =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
    sq_vel_i           += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j           += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
    ProjVelocity_i     += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j     += Velocity_j[iDim]*Normal[iDim];
    Area               += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  /*--- Compute mean values of the variables ---*/
  
  MeanDensity     = 0.5*(DensityInc_i  + DensityInc_j);
  MeanPressure    = 0.5*(Pressure_i    + Pressure_j);
  MeanBetaInc2    = 0.5*(BetaInc2_i    + BetaInc2_j);
  MeanEnthalpy    = 0.5*(Enthalpy_i    + Enthalpy_j);
  MeanCp          = 0.5*(Cp_i          + Cp_j);
  MeanTemperature = 0.5*(Temperature_i + Temperature_j);

  /*--- We need the derivative of the equation of state to build the
   preconditioning matrix. For now, the only option is the ideal gas
   law, but in the future, dRhodT should be in the fluid model. ---*/

  MeandRhodT = 0.0;
  if (variable_density) {
    MeandRhodT = -MeanDensity/MeanTemperature;
  }

  /*--- Get projected flux tensor ---*/

  GetInviscidArtCompProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanBetaInc2, &MeanEnthalpy, Normal, ProjFlux);
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = ProjFlux[iVar];
  }

  /*--- Jacobians of the inviscid flux ---*/
  
  if (implicit) {
    GetInviscidArtCompProjJac(&MeanDensity, MeanVelocity, &MeanBetaInc2, &MeanCp, &MeanTemperature, &MeandRhodT, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
      }
    }
  }
  
  /*--- Computes differences between Laplacians and conservative variables ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
    Diff_V[iVar]    = V_i[iVar]-V_j[iVar];
  }

  /*--- Build the preconditioning matrix using mean values ---*/

  GetPreconditioner(&MeanDensity, MeanVelocity, &MeanBetaInc2, &MeanCp, &MeanTemperature, &MeandRhodT, Precon);

  /*--- Compute the local spectral radius of the preconditioned system
   and the stretching factor. ---*/

  SoundSpeed_i = sqrt(BetaInc2_i*Area*Area);
  SoundSpeed_j = sqrt(BetaInc2_j*Area*Area);
  
  Local_Lambda_i = fabs(ProjVelocity_i)+SoundSpeed_i;
  Local_Lambda_j = fabs(ProjVelocity_j)+SoundSpeed_j;

  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);

  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc2 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  sc4 = sc2*sc2/4.0;
  
  Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
  Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;
  
  /*--- Compute viscous part of the residual ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual[iVar] += Precon[iVar][jVar]*(Epsilon_2*Diff_V[jVar] - Epsilon_4*Diff_Lapl[jVar])*StretchingFactor*MeanLambda;
      if (implicit) {
        val_Jacobian_i[iVar][jVar] += Precon[iVar][jVar]*(Epsilon_2 + Epsilon_4*su2double(Neighbor_i+1))*StretchingFactor*MeanLambda;
        val_Jacobian_j[iVar][jVar] -= Precon[iVar][jVar]*(Epsilon_2 + Epsilon_4*su2double(Neighbor_j+1))*StretchingFactor*MeanLambda;
      }
    }
  }

  /*--- Remove energy contributions if not solving the energy equation. ---*/

  if (!energy) {
    val_residual[nDim+1] = 0.0;
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        val_Jacobian_i[iVar][nDim+1] = 0.0;
        val_Jacobian_j[iVar][nDim+1] = 0.0;

        val_Jacobian_i[nDim+1][iVar] = 0.0;
        val_Jacobian_j[nDim+1][iVar] = 0.0;
      }
    }
  }
  
}

CCentLaxArtComp_Flow::CCentLaxArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  variable_density = (config->GetKind_DensityModel() == VARIABLE);
  grid_movement    = config->GetGrid_Movement();
  energy           = config->GetEnergy_Equation();

  /*--- Artificial dissipation part ---*/

  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_Flow();
  
  /*--- Allocate some structures ---*/

  Diff_V       = new su2double[nVar];
  Velocity_i   = new su2double[nDim];
  Velocity_j   = new su2double[nDim];
  MeanVelocity = new su2double[nDim];
  ProjFlux     = new su2double[nVar];
  Precon       = new su2double*[nVar];

  for (iVar = 0; iVar < nVar; iVar++)
    Precon[iVar] = new su2double[nVar];
  
}

CCentLaxArtComp_Flow::~CCentLaxArtComp_Flow(void) {
  
  delete [] Diff_V;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;

  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Precon[iVar];
  delete [] Precon;

}

void CCentLaxArtComp_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  /*--- Primitive variables at point i and j ---*/
  
  Pressure_i    = V_i[0];             Pressure_j    = V_j[0];
  Temperature_i = V_i[nDim+1];        Temperature_j = V_j[nDim+1];
  DensityInc_i  = V_i[nDim+2];        DensityInc_j  = V_j[nDim+2];
  BetaInc2_i    = V_i[nDim+3];        BetaInc2_j    = V_j[nDim+3];
  Cp_i          = V_i[nDim+7];        Cp_j          = V_j[nDim+7];
  Enthalpy_i    = Cp_i*Temperature_i; Enthalpy_j    = Cp_j*Temperature_j;

  Area = 0.0;
  sq_vel_i = 0.0; sq_vel_j = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim]    = V_i[iDim+1];
    Velocity_j[iDim]    = V_j[iDim+1];
    MeanVelocity[iDim]  = 0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
    sq_vel_i           += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j           += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
    ProjVelocity_i     += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j     += Velocity_j[iDim]*Normal[iDim];
    Area               += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);

  /*--- Compute mean values of the variables ---*/

  MeanDensity     = 0.5*(DensityInc_i  + DensityInc_j);
  MeanPressure    = 0.5*(Pressure_i    + Pressure_j);
  MeanBetaInc2    = 0.5*(BetaInc2_i    + BetaInc2_j);
  MeanEnthalpy    = 0.5*(Enthalpy_i    + Enthalpy_j);
  MeanCp          = 0.5*(Cp_i          + Cp_j);
  MeanTemperature = 0.5*(Temperature_i + Temperature_j);

  /*--- We need the derivative of the equation of state to build the
   preconditioning matrix. For now, the only option is the ideal gas
   law, but in the future, dRhodT should be in the fluid model. ---*/

  MeandRhodT = 0.0;
  if (variable_density) {
    MeandRhodT = -MeanDensity/MeanTemperature;
  }

  /*--- Get projected flux tensor ---*/

  GetInviscidArtCompProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanBetaInc2, &MeanEnthalpy, Normal, ProjFlux);

  /*--- Compute inviscid residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = ProjFlux[iVar];
  }

  /*--- Jacobians of the inviscid flux ---*/

  if (implicit) {
    GetInviscidArtCompProjJac(&MeanDensity, MeanVelocity, &MeanBetaInc2, &MeanCp, &MeanTemperature, &MeandRhodT, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
      }
    }
  }
  
  /*--- Computes differences btw. conservative variables ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_V[iVar] = V_i[iVar]-V_j[iVar];

  /*--- Build the preconditioning matrix using mean values ---*/

  GetPreconditioner(&MeanDensity, MeanVelocity, &MeanBetaInc2, &MeanCp, &MeanTemperature, &MeandRhodT, Precon);

  /*--- Compute the local espectral radius of the preconditioned system
   and the stretching factor. ---*/

  SoundSpeed_i = sqrt(BetaInc2_i*Area*Area);
  SoundSpeed_j = sqrt(BetaInc2_j*Area*Area);

  Local_Lambda_i = fabs(ProjVelocity_i)+SoundSpeed_i;
  Local_Lambda_j = fabs(ProjVelocity_j)+SoundSpeed_j;

  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);

  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc0 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc0*su2double(nDim)/3.0;
  
  /*--- Compute viscous part of the residual ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual[iVar] += Precon[iVar][jVar]*Epsilon_0*Diff_V[jVar]*StretchingFactor*MeanLambda;
      if (implicit) {
        val_Jacobian_i[iVar][jVar] += Precon[iVar][jVar]*Epsilon_0*StretchingFactor*MeanLambda;
        val_Jacobian_j[iVar][jVar] -= Precon[iVar][jVar]*Epsilon_0*StretchingFactor*MeanLambda;
      }
    }
  }
  
  /*--- Remove energy contributions if we aren't solving the energy equation. ---*/

  if (!energy) {
    val_residual[nDim+1] = 0.0;
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        val_Jacobian_i[iVar][nDim+1] = 0.0;
        val_Jacobian_j[iVar][nDim+1] = 0.0;

        val_Jacobian_i[nDim+1][iVar] = 0.0;
        val_Jacobian_j[nDim+1][iVar] = 0.0;
      }
    }
  }

}

CAvgGradArtComp_Flow::CAvgGradArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  energy   = config->GetEnergy_Equation();

  /*--- Incompressible flow, primitive variables nDim+2, (P, vx, vy, vz, T) ---*/ 
  
  Mean_GradPrimVar = new su2double*[nVar];
  
  /*--- Incompressible flow, gradient primitive variables nDim+2, (P, vx, vy, vz, T) ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradPrimVar[iVar] = new su2double[nDim];
  
}

CAvgGradArtComp_Flow::~CAvgGradArtComp_Flow(void) {
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
  
}

void CAvgGradArtComp_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  /*--- Normalized normal vector ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Density and transport properties ---*/
  
  Laminar_Viscosity_i    = V_i[nDim+4];  Laminar_Viscosity_j    = V_j[nDim+4];
  Eddy_Viscosity_i       = V_i[nDim+5];  Eddy_Viscosity_j       = V_j[nDim+5];
  Thermal_Conductivity_i = V_i[nDim+6];  Thermal_Conductivity_j = V_j[nDim+6];

  /*--- Mean transport properties ---*/
  
  Mean_Laminar_Viscosity    = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity       = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_Thermal_Conductivity = 0.5*(Thermal_Conductivity_i + Thermal_Conductivity_j);

  /*--- Mean gradient approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (iDim = 0; iDim < nDim; iDim++)
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
  
  /*--- Get projected flux tensor ---*/
  
  GetViscousArtCompProjFlux(Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, Mean_Thermal_Conductivity);
  
  /*--- Update viscous residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
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
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    else {
      GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, dist_ij, UnitNormal,
                                Area, val_Jacobian_i, val_Jacobian_j);
    }
    
  }

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  if (implicit) {
    su2double Edge_Vector[3];
    su2double dist_ij_2 = 0.0, proj_vector_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
      dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
      proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
    }
    if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
    else proj_vector_ij = proj_vector_ij/dist_ij_2;

    val_Jacobian_i[nDim+1][nDim+1] = -Mean_Thermal_Conductivity*proj_vector_ij;
    val_Jacobian_j[nDim+1][nDim+1] =  Mean_Thermal_Conductivity*proj_vector_ij;
  }

  if (!energy) {
    val_residual[nDim+1] = 0.0;
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        val_Jacobian_i[iVar][nDim+1] = 0.0;
        val_Jacobian_j[iVar][nDim+1] = 0.0;

        val_Jacobian_i[nDim+1][iVar] = 0.0;
        val_Jacobian_j[nDim+1][iVar] = 0.0;
      }
    }
  }

}

CAvgGradCorrectedArtComp_Flow::CAvgGradCorrectedArtComp_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  energy   = config->GetEnergy_Equation();

  PrimVar_i = new su2double [nVar];
  PrimVar_j = new su2double [nVar];
  Proj_Mean_GradPrimVar_Edge = new su2double [nVar];
  Edge_Vector = new su2double [nDim];
  
  Mean_GradPrimVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradPrimVar[iVar] = new su2double [nDim];
  
}

CAvgGradCorrectedArtComp_Flow::~CAvgGradCorrectedArtComp_Flow(void) {
  
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Proj_Mean_GradPrimVar_Edge;
  delete [] Edge_Vector;
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
  
}

void CAvgGradCorrectedArtComp_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+8);   AD::SetPreaccIn(V_j, nDim+8);
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(PrimVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(PrimVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(Normal, nDim);
  
  /*--- Normalized normal vector ---*/

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Conversion to Primitive Variables (P, u, v, w, T) ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
  }
  
  /*--- Density and transport properties ---*/
  
  DensityInc_i           = V_i[nDim+2];  DensityInc_j           = V_j[nDim+2];
  Laminar_Viscosity_i    = V_i[nDim+4];  Laminar_Viscosity_j    = V_j[nDim+4];
  Eddy_Viscosity_i       = V_i[nDim+5];  Eddy_Viscosity_j       = V_j[nDim+5];
  Thermal_Conductivity_i = V_i[nDim+6];  Thermal_Conductivity_j = V_j[nDim+6];
  Cp_i                   = V_i[nDim+7];  Cp_j                   = V_j[nDim+7];

  /*--- Mean transport properties ---*/
  
  Mean_Laminar_Viscosity    = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity       = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_Thermal_Conductivity = 0.5*(Thermal_Conductivity_i + Thermal_Conductivity_j);
  Mean_Cp                   = 0.5*(Cp_i + Cp_j);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }
  
  /*--- Projection of the mean gradient in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
    }
    if (dist_ij_2 != 0.0) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                         (PrimVar_j[iVar]-PrimVar_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
      }
    }
  }
  
  /*--- Get projected flux tensor ---*/
  
  GetViscousArtCompProjFlux(Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, Mean_Thermal_Conductivity);
  
  /*--- Update viscous residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  }
  
  /*--- Implicit part for conservative portion ---*/
  
  if (implicit) {
    
    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    else {
      GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, sqrt(dist_ij_2), UnitNormal,
                                Area, val_Jacobian_i, val_Jacobian_j);
    }
  }

    /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
    
    if (implicit) {
      su2double proj_vector_ij = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
      if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
      else proj_vector_ij = proj_vector_ij/dist_ij_2;

      val_Jacobian_i[nDim+1][nDim+1] = -Mean_Thermal_Conductivity*proj_vector_ij;
      val_Jacobian_j[nDim+1][nDim+1] =  Mean_Thermal_Conductivity*proj_vector_ij;
    }

  if (!energy) {
    val_residual[nDim+1] = 0.0;
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        val_Jacobian_i[iVar][nDim+1] = 0.0;
        val_Jacobian_j[iVar][nDim+1] = 0.0;

        val_Jacobian_i[nDim+1][iVar] = 0.0;
        val_Jacobian_j[nDim+1][iVar] = 0.0;
      }
    }
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}

CSourceIncBodyForce::CSourceIncBodyForce(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  /*--- Store the pointer to the constant body force vector. ---*/

  Body_Force_Vector = new su2double[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Body_Force_Vector[iDim] = config->GetBody_Force_Vector()[iDim];

}

CSourceIncBodyForce::~CSourceIncBodyForce(void) {

  if (Body_Force_Vector != NULL) delete [] Body_Force_Vector;

}

void CSourceIncBodyForce::ComputeResidual(su2double *val_residual, CConfig *config) {

  unsigned short iDim;
  su2double DensityInc_0 = config->GetDensity_FreeStreamND();
  su2double Force_Ref    = config->GetForce_Ref();

    /*--- Zero the continuity contribution ---*/

    val_residual[0] = 0.0;

  /*--- Momentum contribution. Note that this form assumes we have
   subtracted the operating density * gravity, i.e., removed the
   hydrostatic pressure component (important for pressure BCs). ---*/

    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = -Volume * (DensityInc_i - DensityInc_0) * Body_Force_Vector[iDim] / Force_Ref;

    /*--- Zero the temperature contribution ---*/

    val_residual[nDim+1] = 0.0;

}

CSourceBoussinesq::CSourceBoussinesq(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  /*--- Store the pointer to the constant body force vector. ---*/

  Gravity_Vector = new su2double[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Gravity_Vector[iDim] = 0.0;

  /*--- Gravity is downward in y-dir for 2D and downward z-dir for 3D. ---*/

  Gravity_Vector[nDim-1] = -STANDARD_GRAVITY;

}

CSourceBoussinesq::~CSourceBoussinesq(void) {

  if (Gravity_Vector != NULL) delete [] Gravity_Vector;

}

void CSourceBoussinesq::ComputeResidual(su2double *val_residual, CConfig *config) {

  unsigned short iDim;
  su2double Force_Ref = config->GetForce_Ref();
  su2double T0        = config->GetTemperature_FreeStreamND();
  su2double Beta      = config->GetThermal_Expansion_CoeffND();

  /*--- Zero the continuity contribution ---*/

  val_residual[0] = 0.0;

  /*--- Momentum contribution. Note that this form assumes we have
   subtracted the operating density * gravity, i.e., removed the
   hydrostatic pressure component (important for pressure BCs). ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[iDim+1] = Volume * DensityInc_i * ( Beta * (U_i[nDim+1] - T0)) * Gravity_Vector[iDim] / Force_Ref; 

  /*--- Zero the energy contribution ---*/

  val_residual[nDim+1] = 0.0;

}
