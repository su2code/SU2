/*!
 * \file numerics_direct_mean.cpp
 * \brief This file contains the numerical methods for compressible flow.
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

CCentJST_Flow::CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  grid_movement = config->GetGrid_Movement();
  
  /*--- Artifical dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_Flow();
  Param_Kappa_4 = config->GetKappa_4th_Flow();
  
  /*--- Allocate some structures ---*/
  Diff_U = new su2double [nVar];
  Diff_Lapl = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  MeanVelocity = new su2double [nDim];
  ProjFlux = new su2double [nVar];
  
}

CCentJST_Flow::~CCentJST_Flow(void) {
  delete [] Diff_U;
  delete [] Diff_Lapl;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
}

void CCentJST_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                                    CConfig *config) {
  
  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);
  AD::SetPreaccIn(Sensor_i);    AD::SetPreaccIn(Sensor_j);
  AD::SetPreaccIn(Lambda_i);    AD::SetPreaccIn(Lambda_j);
  AD::SetPreaccIn(Und_Lapl_i, nVar); AD::SetPreaccIn(Und_Lapl_j, nVar);
  if (grid_movement) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }

  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/
  
  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i = V_i[nDim+2];                        Density_j = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  SoundSpeed_i = V_i[nDim+4];                     SoundSpeed_j = V_j[nDim+4];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;
  
  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }

  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Compute mean values of the variables ---*/
  
  MeanDensity = 0.5*(Density_i+Density_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  MeanEnergy = 0.5*(Energy_i+Energy_j);
  
  /*--- Get projected flux tensor ---*/
  
  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, ProjFlux);
  
  /*--- Residual of the inviscid flux ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];
  
  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/

  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }

  /*--- Adjustment due to grid motion ---*/
  
  if (grid_movement) {
    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar] + U_j[iVar]);
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }
  
  /*--- Computes differences btw. Laplacians and conservative variables,
   with a correction for the enthalpy ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  }
  Diff_U[nVar-1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;
  
  /*--- Compute the local spectral radius and the stretching factor ---*/
  
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  /*--- Adjustment due to mesh motion ---*/
  
  if (grid_movement) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }
  
  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc2 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  sc4 = sc2*sc2/4.0;
  
  Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
  Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;
  
  /*--- Compute viscous part of the residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += (Epsilon_2*Diff_U[iVar] - Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
  
  /*--- Jacobian computation ---*/

  if (implicit) {

    cte_0 = (Epsilon_2 + Epsilon_4*su2double(Neighbor_i+1))*StretchingFactor*MeanLambda;
    cte_1 = (Epsilon_2 + Epsilon_4*su2double(Neighbor_j+1))*StretchingFactor*MeanLambda;
    
    for (iVar = 0; iVar < (nVar-1); iVar++) {
      val_Jacobian_i[iVar][iVar] += cte_0;
      val_Jacobian_j[iVar][iVar] -= cte_1;
    }
    
    /*--- Last row of Jacobian_i ---*/
    
    val_Jacobian_i[nVar-1][0] += cte_0*Gamma_Minus_One*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_i[nVar-1][iDim+1] -= cte_0*Gamma_Minus_One*Velocity_i[iDim];
    val_Jacobian_i[nVar-1][nVar-1] += cte_0*Gamma;
    
    /*--- Last row of Jacobian_j ---*/
    
    val_Jacobian_j[nVar-1][0] -= cte_1*Gamma_Minus_One*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_j[nVar-1][iDim+1] += cte_1*Gamma_Minus_One*Velocity_j[iDim];
    val_Jacobian_j[nVar-1][nVar-1] -= cte_1*Gamma;
    
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}

CCentJST_KE_Flow::CCentJST_KE_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  grid_movement = config->GetGrid_Movement();

  /*--- Artifical dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_Flow();
  Param_Kappa_4 = config->GetKappa_4th_Flow();

  /*--- Allocate some structures ---*/
  Diff_U = new su2double [nVar];
  Diff_Lapl = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  MeanVelocity = new su2double [nDim];
  ProjFlux = new su2double [nVar];

}

CCentJST_KE_Flow::~CCentJST_KE_Flow(void) {
  delete [] Diff_U;
  delete [] Diff_Lapl;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
}

void CCentJST_KE_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                                    CConfig *config) {

  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);
  AD::SetPreaccIn(Sensor_i);    AD::SetPreaccIn(Sensor_j);
  AD::SetPreaccIn(Lambda_i);    AD::SetPreaccIn(Lambda_j);
  AD::SetPreaccIn(Und_Lapl_i, nVar); AD::SetPreaccIn(Und_Lapl_j, nVar);

  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/

  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i = V_i[nDim+2];                        Density_j = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  SoundSpeed_i = V_i[nDim+4];                     SoundSpeed_j = V_j[nDim+4];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;

  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }

  /*--- Recompute conservative variables ---*/

  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;

  /*--- Compute mean values of the variables ---*/

  MeanDensity = 0.5*(Density_i+Density_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  MeanEnergy = 0.5*(Energy_i+Energy_j);

  /*--- Get projected flux tensor ---*/

  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, ProjFlux);

  /*--- Residual of the inviscid flux ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];

  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/

  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }

  /*--- Adjustment due to grid motion ---*/

  if (grid_movement) {
    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }

  /*--- Computes differences btw. Laplacians and conservative variables,
   with a correction for the enthalpy ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  }
  Diff_U[nVar-1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;

  /*--- Compute the local spectral radius and the stretching factor ---*/

  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);

  /*--- Adjustment due to mesh motion ---*/

  if (grid_movement) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }

  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);

  sc2 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  sc4 = sc2*sc2/4.0;

  Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;

  /*--- Compute viscous part of the residual ---*/

  for (iVar = 0; iVar < nVar; iVar++)
      val_residual[iVar] += Epsilon_2*(Diff_U[iVar])*StretchingFactor*MeanLambda;

  /*--- Jacobian computation ---*/

  if (implicit) {

    cte_0 = Epsilon_2*StretchingFactor*MeanLambda;

    for (iVar = 0; iVar < (nVar-1); iVar++) {
      val_Jacobian_i[iVar][iVar] += cte_0;
      val_Jacobian_j[iVar][iVar] -= cte_0;
    }

    /*--- Last row of Jacobian_i ---*/

    val_Jacobian_i[nVar-1][0] += cte_0*Gamma_Minus_One*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_i[nVar-1][iDim+1] -= cte_0*Gamma_Minus_One*Velocity_i[iDim];
    val_Jacobian_i[nVar-1][nVar-1] += cte_0*Gamma;

    /*--- Last row of Jacobian_j ---*/

    val_Jacobian_j[nVar-1][0] -= cte_1*Gamma_Minus_One*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_j[nVar-1][iDim+1] += cte_1*Gamma_Minus_One*Velocity_j[iDim];
    val_Jacobian_j[nVar-1][nVar-1] -= cte_1*Gamma;

  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}


CCentLax_Flow::CCentLax_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  
  /*--- Artifical dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_Flow();
  
  /*--- Allocate some structures ---*/
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  MeanVelocity = new su2double [nDim];
  ProjFlux = new su2double [nVar];
  
}

CCentLax_Flow::~CCentLax_Flow(void) {
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
  
}

void CCentLax_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                                    CConfig *config) {
  
  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/
  
  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i = V_i[nDim+2];                        Density_j = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  SoundSpeed_i = V_i[nDim+4];                     SoundSpeed_j = V_j[nDim+4];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;
  
  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }
  
  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Compute mean values of the variables ---*/
  
  MeanDensity = 0.5*(Density_i+Density_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  MeanEnergy = 0.5*(Energy_i+Energy_j);
  
  /*--- Get projected flux tensor ---*/
  
  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, ProjFlux);
  
  /*--- Residual of the inviscid flux ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];
  
  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
  
  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }
  
  /*--- Adjustment due to grid motion ---*/
  
  if (grid_movement) {
    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }
  
  /*--- Computes differences btw. conservative variables,
   with a correction for the enthalpy ---*/
  
  for (iVar = 0; iVar < nDim+1; iVar++)
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  Diff_U[nDim+1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;
  
  /*--- Compute the local spectral radius and the stretching factor ---*/
  
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  /*--- Adjustment due to grid motion ---*/
  if (grid_movement) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }
  
  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc0 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc0*su2double(nDim)/3.0;
  
  /*--- Compute viscous part of the residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;
  
  /*--- Jacobian computation ---*/

  if (implicit) {
    cte = Epsilon_0*StretchingFactor*MeanLambda;
    for (iVar = 0; iVar < (nVar-1); iVar++) {
      val_Jacobian_i[iVar][iVar] += cte;
      val_Jacobian_j[iVar][iVar] -= cte;
    }
    
    /*--- Last row of Jacobian_i ---*/
    
    val_Jacobian_i[nVar-1][0] += cte*Gamma_Minus_One*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_i[nVar-1][iDim+1] -= cte*Gamma_Minus_One*Velocity_i[iDim];
    val_Jacobian_i[nVar-1][nVar-1] += cte*Gamma;
    
    /*--- Last row of Jacobian_j ---*/
    
    val_Jacobian_j[nVar-1][0] -= cte*Gamma_Minus_One*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_j[nVar-1][iDim+1] += cte*Gamma_Minus_One*Velocity_j[iDim];
    val_Jacobian_j[nVar-1][nVar-1] -= cte*Gamma;
    
  }
  
}

CUpwCUSP_Flow::CUpwCUSP_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  grid_movement = config->GetGrid_Movement();
  
  /*--- Allocate some structures ---*/
  Diff_U = new su2double [nVar];
  Diff_Flux = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  MeanVelocity = new su2double [nDim];
  ProjFlux = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Jacobian = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian[iVar] = new su2double [nVar];
  }
}

CUpwCUSP_Flow::~CUpwCUSP_Flow(void) {
  delete [] Diff_U;
  delete [] Diff_Flux;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian[iVar];
  }
}

void CUpwCUSP_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                                     CConfig *config) {
  
  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/
  
  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i  = V_i[nDim+2];                       Density_j  = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  Energy_i   = Enthalpy_i - Pressure_i/Density_i; Energy_j   = Enthalpy_j - Pressure_j/Density_j;

  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += Velocity_j[iDim]*Velocity_j[iDim];
  }
  
  SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel_i));
  SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel_j));

  /*-- Face area ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit normal ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Compute mean values of the variables ---*/
  
  MeanDensity = 0.5*(Density_i+Density_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
  ProjVelocity = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
    ProjVelocity +=  MeanVelocity[iDim]*UnitNormal[iDim];
  }
  MeanSoundSpeed = 0.5*(SoundSpeed_i+SoundSpeed_j);
  MeanEnergy = 0.5*(Energy_i+Energy_j);
  
  /*--- Get projected flux tensor ---*/
  
  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, ProjFlux);
  
  /*--- Residual of the inviscid flux ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];
  
  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/

  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }
  
  /*--- Computes differences conservative variables,
   with a correction for the enthalpy ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  Diff_U[nVar-1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;
  
  /*--- Computes differences projected fluxes,
   with a correction for the enthalpy ---*/
  
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, UnitNormal, ProjFlux_i);
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, UnitNormal, ProjFlux_j);
  
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_Flux[iVar] = ProjFlux_i[iVar]-ProjFlux_j[iVar];
  
  /*--- Compute dissipation parameters ---*/
  
  Mach = ProjVelocity / MeanSoundSpeed;
  
  LamdaNeg = ProjVelocity - MeanSoundSpeed;
  LamdaPos = ProjVelocity + MeanSoundSpeed;
  
  if ((0.0 <= Mach) && (Mach < 1.0)) Beta = + max(0.0, (ProjVelocity + LamdaNeg)/(ProjVelocity - LamdaNeg));
  if ((-1.0 <= Mach) && (Mach < 0.0)) Beta = - max(0.0, (ProjVelocity + LamdaPos)/(ProjVelocity - LamdaPos));
  if (fabs(Mach) >= 1.0) Beta = Mach/fabs(Mach);
  
  if (Beta == 0.0) Nu_c = fabs(ProjVelocity);
  if ((Beta > 0.0) && ((0.0 < Mach) && (Mach < 1.0))) Nu_c = - (1.0+Beta)*LamdaNeg;
  if ((Beta < 0.0) && ((-1.0 < Mach) && (Mach < 0.0))) Nu_c = (1.0-Beta)*LamdaPos;
  if (fabs(Mach) >= 1) Nu_c = 0.0;
  
  /*--- Compute the residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += (0.5*Nu_c*Diff_U[iVar] + 0.5*Beta*Diff_Flux[iVar])*Area;

  /*--- Jacobian computation ---*/

  if (implicit) {
    
    cte_0 = 0.5*Nu_c*Area;
    cte_1 = 0.5*Beta*Area;
    
    for (iVar = 0; iVar < (nVar-1); iVar++) {
      val_Jacobian_i[iVar][iVar] += cte_0;
      val_Jacobian_j[iVar][iVar] -= cte_0;
    }
    
    /*--- Last row of Jacobian_i (solution difference contribution) ---*/
    
    val_Jacobian_i[nVar-1][0] += cte_0*Gamma_Minus_One*0.5*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_i[nVar-1][iDim+1] -= cte_0*Gamma_Minus_One*Velocity_i[iDim];
    val_Jacobian_i[nVar-1][nVar-1] += cte_0*Gamma;
    
    /*--- Last row of Jacobian_j (solution difference contribution) ---*/
    
    val_Jacobian_j[nVar-1][0] -= cte_0*Gamma_Minus_One*0.5*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_j[nVar-1][iDim+1] += cte_0*Gamma_Minus_One*Velocity_j[iDim];
    val_Jacobian_j[nVar-1][nVar-1] -= cte_0*Gamma;
    
    /*--- Flux difference contribution ---*/
    
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 1.0, Jacobian);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[iVar][jVar] += cte_1*Jacobian[iVar][jVar];
    
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 1.0, Jacobian);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] -= cte_1*Jacobian[iVar][jVar];
    
  }
  
}

CUpwAUSM_Flow::CUpwAUSM_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  delta_vel  = new su2double [nDim];
  delta_wave = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwAUSM_Flow::~CUpwAUSM_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwAUSM_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  /*--- Face area (norm or the normal vector) ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel)));
  
  /*--- Primitive variables at point j ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
  }
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel)));
  
  /*--- Projected velocities ---*/
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  mL  = ProjVelocity_i/SoundSpeed_i;
  mR  = ProjVelocity_j/SoundSpeed_j;
  
  if (fabs(mL) <= 1.0) mLP = 0.25*(mL+1.0)*(mL+1.0);
  else mLP = 0.5*(mL+fabs(mL));
  
  if (fabs(mR) <= 1.0) mRM = -0.25*(mR-1.0)*(mR-1.0);
  else mRM = 0.5*(mR-fabs(mR));
  
  mF = mLP + mRM;
  
  if (fabs(mL) <= 1.0) pLP = 0.25*Pressure_i*(mL+1.0)*(mL+1.0)*(2.0-mL);
  else pLP = 0.5*Pressure_i*(mL+fabs(mL))/mL;
  
  if (fabs(mR) <= 1.0) pRM = 0.25*Pressure_j*(mR-1.0)*(mR-1.0)*(2.0+mR);
  else pRM = 0.5*Pressure_j*(mR-fabs(mR))/mR;
  
  pF = pLP + pRM;
  Phi = fabs(mF);
  
  val_residual[0] = 0.5*(mF*((Density_i*SoundSpeed_i)+(Density_j*SoundSpeed_j))-Phi*((Density_j*SoundSpeed_j)-(Density_i*SoundSpeed_i)));
  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[iDim+1] = 0.5*(mF*((Density_i*SoundSpeed_i*Velocity_i[iDim])+(Density_j*SoundSpeed_j*Velocity_j[iDim]))
                                -Phi*((Density_j*SoundSpeed_j*Velocity_j[iDim])-(Density_i*SoundSpeed_i*Velocity_i[iDim])))+UnitNormal[iDim]*pF;
  val_residual[nVar-1] = 0.5*(mF*((Density_i*SoundSpeed_i*Enthalpy_i)+(Density_j*SoundSpeed_j*Enthalpy_j))-Phi*((Density_j*SoundSpeed_j*Enthalpy_j)-(Density_i*SoundSpeed_i*Enthalpy_i)));

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] *= Area;
  
  /*--- Roe's Jacobian for AUSM (this must be fixed) ---*/
  if (implicit) {
    
    /*--- Mean Roe variables iPoint and jPoint ---*/
    R = sqrt(fabs(Density_j/Density_i));
    RoeDensity = R*Density_i;
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
      sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
    }
    RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
    RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
    
    /*--- Compute P and Lambda (do it with the Normal) ---*/
    GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
    
    ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
      ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
    }
    
    /*--- Flow eigenvalues and Entropy correctors ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Lambda[iDim] = ProjVelocity;
    Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
    Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
    
    /*--- Compute inverse P ---*/
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobias of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
    
    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }
}

CUpwAUSMPLUSUP_Flow::CUpwAUSMPLUSUP_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  Minf = config->GetMach();
  Kp = 0.25; Ku = 0.75; sigma = 1.0;

  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  delta_vel  = new su2double [nDim];
  delta_wave = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwAUSMPLUSUP_Flow::~CUpwAUSMPLUSUP_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwAUSMPLUSUP_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  /*--- Face area (norm or the normal vector) ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel)));
  
  /*--- Primitive variables at point j ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
  }
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel)));
  
  /*--- Projected velocities ---*/
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  astarL = sqrt(2.0*(Gamma-1.0)/(Gamma+1.0)*Enthalpy_i); astarR = sqrt(2.0*(Gamma-1.0)/(Gamma+1.0)*Enthalpy_j);

 // ahatL = astarL*astarL/max(astarL,fabs(ProjVelocity_i));
 // ahatR = astarR*astarR/max(astarR,fabs(ProjVelocity_j));

  ahatL = astarL*astarL/max(astarL,ProjVelocity_i);
  ahatR = astarR*astarR/max(astarR,-ProjVelocity_j);
 
  /*--- Interface speed of sound ---*/ 
  aF = min(ahatL,ahatR);

  mL  = ProjVelocity_i/aF;
  mR  = ProjVelocity_j/aF;
  
  rhoF = 0.5*(Density_i+Density_j);
  MFsq = 0.5*(mL*mL+mR*mR);

  param1 = max(MFsq, Minf*Minf);
  Mrefsq = (min(1.0, param1));
  fa = 2.0*sqrt(Mrefsq)-Mrefsq;

  alpha = 3.0/16.0*(-4.0+5.0*fa*fa);
  beta = 1.0/8.0;

  /*--- Pressure diffusion term ---*/
  Mp = -(Kp/fa)*max((1.0-sigma*MFsq),0.0)*(Pressure_j-Pressure_i)/(rhoF*aF*aF);

  if (fabs(mL) <= 1.0) mLP = 0.25*(mL+1.0)*(mL+1.0)+beta*(mL*mL-1.0)*(mL*mL-1.0);
  else mLP = 0.5*(mL+fabs(mL));
  
  if (fabs(mR) <= 1.0) mRM = -0.25*(mR-1.0)*(mR-1.0)-beta*(mR*mR-1.0)*(mR*mR-1.0);
  else mRM = 0.5*(mR-fabs(mR));
  
  mF = mLP + mRM + Mp;
  
  if (fabs(mL) <= 1.0) pLP = Pressure_i*(0.25*(mL+1.0)*(mL+1.0)*(2.0-mL)+alpha*mL*(mL*mL-1.0)*(mL*mL-1.0));
  else pLP = 0.5*Pressure_i*(mL+fabs(mL))/mL;
  
  if (fabs(mR) <= 1.0) pRM = Pressure_j*(0.25*(mR-1.0)*(mR-1.0)*(2.0+mR)-alpha*mR*(mR*mR-1.0)*(mR*mR-1.0));  
  else pRM = 0.5*Pressure_j*(mR-fabs(mR))/mR;
  
  /*--- Velocity diffusion term ---*/
  Pu = -Ku*fa*(pLP/Pressure_i)*(pRM/Pressure_j)*(Density_i+Density_j)*aF*(ProjVelocity_j-ProjVelocity_i);

  pF = pLP+pRM+Pu;
  Phi = fabs(mF);
 
  mfP=0.5*(mF+Phi);
  mfM=0.5*(mF-Phi);
 
  val_residual[0] = aF*(mfP*Density_i+mfM*Density_j);  
  
  for (iDim = 0; iDim < nDim; iDim++)
  	val_residual[iDim+1] = aF*(mfP*Density_i*Velocity_i[iDim]+mfM*Density_j*Velocity_j[iDim])+UnitNormal[iDim]*pF;
  val_residual[nVar-1] = aF*(mfP*Density_i*Enthalpy_i+mfM*Density_j*Enthalpy_j);
  

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] *= Area;
  
  /*--- Roe's Jacobian for AUSM (this must be fixed) ---*/
  if (implicit) {
    
    /*--- Mean Roe variables iPoint and jPoint ---*/
    R = sqrt(fabs(Density_j/Density_i));
    RoeDensity = R*Density_i;
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
      sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
    }
    RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
    RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
    
    /*--- Compute P and Lambda (do it with the Normal) ---*/
    GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
    
    ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
      ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
    }
    
    /*--- Flow eigenvalues and Entropy correctors ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Lambda[iDim] = ProjVelocity;
    Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
    Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
    
    /*--- Compute inverse P ---*/
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobias of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
    
    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }
}

CUpwAUSMPLUSUP2_Flow::CUpwAUSMPLUSUP2_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  Minf = config->GetMach();
  Kp = 0.25; sigma = 1.0;

  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  delta_vel  = new su2double [nDim];
  delta_wave = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwAUSMPLUSUP2_Flow::~CUpwAUSMPLUSUP2_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwAUSMPLUSUP2_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  /*--- Face area (norm or the normal vector) ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
  }
  sq_veli=sq_vel; 
  
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel)));
  
  /*--- Primitive variables at point j ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
  }
  sq_velj=sq_vel; 

  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel)));
  
  /*--- Projected velocities ---*/
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  astarL = sqrt(2.0*(Gamma-1.0)/(Gamma+1.0)*Enthalpy_i); astarR = sqrt(2.0*(Gamma-1.0)/(Gamma+1.0)*Enthalpy_j);

 // ahatL = astarL*astarL/max(astarL,fabs(ProjVelocity_i));
 // ahatR = astarR*astarR/max(astarR,fabs(ProjVelocity_j));

  ahatL = astarL*astarL/max(astarL,ProjVelocity_i);
  ahatR = astarR*astarR/max(astarR,-ProjVelocity_j);
 
  /*--- Interface speed of sound ---*/ 
  aF = min(ahatL,ahatR);

  mL  = ProjVelocity_i/aF;
  mR  = ProjVelocity_j/aF;
  
  rhoF = 0.5*(Density_i+Density_j);
  MFsq = 0.5*(mL*mL+mR*mR);

  param1 = max(MFsq, Minf*Minf);
  Mrefsq = (min(1.0, param1));
  fa = 2.0*sqrt(Mrefsq)-Mrefsq;

  alpha = 3.0/16.0*(-4.0+5.0*fa*fa);
  beta = 1.0/8.0;

  /*--- Pressure diffusion term ---*/
  Mp = -(Kp/fa)*max((1.0-sigma*MFsq),0.0)*(Pressure_j-Pressure_i)/(rhoF*aF*aF);

  if (fabs(mL) <= 1.0) mLP = 0.25*(mL+1.0)*(mL+1.0)+beta*(mL*mL-1.0)*(mL*mL-1.0);
  else mLP = 0.5*(mL+fabs(mL));
  
  if (fabs(mR) <= 1.0) mRM = -0.25*(mR-1.0)*(mR-1.0)-beta*(mR*mR-1.0)*(mR*mR-1.0);
  else mRM = 0.5*(mR-fabs(mR));
  
  mF = mLP + mRM + Mp;
  
  if (fabs(mL) <= 1.0) pLP = (0.25*(mL+1.0)*(mL+1.0)*(2.0-mL)+alpha*mL*(mL*mL-1.0)*(mL*mL-1.0));
  else pLP = 0.5*(mL+fabs(mL))/mL;
  
  if (fabs(mR) <= 1.0) pRM = (0.25*(mR-1.0)*(mR-1.0)*(2.0+mR)-alpha*mR*(mR*mR-1.0)*(mR*mR-1.0));  
  else pRM = 0.5*(mR-fabs(mR))/mR;
  
  /*... Modified pressure flux ...*/
  //Use this definition
  pFi=sqrt(0.5*(sq_veli+sq_velj))*(pLP+pRM-1.0)*0.5*(Density_j+Density_i)*aF;
  //pFi=sqrt(0.5*(sq_veli+sq_velj))*(pLP+pRM-1.0)*0.5*(Pressure_i+Pressure_j)/aF;  
  pF=0.5*(Pressure_j+Pressure_i)+0.5*(pLP-pRM)*(Pressure_i-Pressure_j)+pFi;  

  Phi = fabs(mF);
 
  mfP=0.5*(mF+Phi);
  mfM=0.5*(mF-Phi);
 
  val_residual[0] = aF*(mfP*Density_i+mfM*Density_j);  
  
  for (iDim = 0; iDim < nDim; iDim++)
  	val_residual[iDim+1] = aF*(mfP*Density_i*Velocity_i[iDim]+mfM*Density_j*Velocity_j[iDim])+UnitNormal[iDim]*pF;
  val_residual[nVar-1] = aF*(mfP*Density_i*Enthalpy_i+mfM*Density_j*Enthalpy_j);
  

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] *= Area;
  
  /*--- Roe's Jacobian for AUSM (this must be fixed) ---*/
  if (implicit) {
    
    /*--- Mean Roe variables iPoint and jPoint ---*/
    R = sqrt(fabs(Density_j/Density_i));
    RoeDensity = R*Density_i;
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
      sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
    }
    RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
    RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
    
    /*--- Compute P and Lambda (do it with the Normal) ---*/
    GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
    
    ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
      ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
    }
    
    /*--- Flow eigenvalues and Entropy correctors ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Lambda[iDim] = ProjVelocity;
    Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
    Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
    
    /*--- Compute inverse P ---*/
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobias of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
    
    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }
}

CUpwSLAU_Flow::CUpwSLAU_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config, bool val_low_dissipation) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  slau_low_diss = val_low_dissipation;
  
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  delta_vel  = new su2double [nDim];
  delta_wave = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwSLAU_Flow::~CUpwSLAU_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwSLAU_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
   
  /*--- Face area (norm or the normal vector) ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel)));
  
  /*--- Primitive variables at point j ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
  }
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel)));
  
  /*--- Projected velocities ---*/
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
    
  aF = 0.5 * (SoundSpeed_i + SoundSpeed_j);
  
  mL  = ProjVelocity_i/aF;
  mR  = ProjVelocity_j/aF;
 
  /*--- Smooth function of the local Mach number---*/
  aux_slau = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    aux_slau += Velocity_i[iDim]*Velocity_i[iDim];
    aux_slau += Velocity_j[iDim]*Velocity_j[iDim];
  }
  Mach_tilde = min(1.0, (1.0/aF) * sqrt(aux_slau/2.0));  
  Chi = pow((1.0 - Mach_tilde),2.0);
  f_rho = -max(min(mL,0.0),-1.0) * min(max(mR,0.0),1.0);

  /*--- Mean normal velocity with density weighting ---*/
  Vn_Mag = (Density_i*fabs(ProjVelocity_i) + Density_j*fabs(ProjVelocity_j)) / (Density_i + Density_j);
  Vn_MagL= (1.0 - f_rho)*Vn_Mag + f_rho*fabs(ProjVelocity_i);
  Vn_MagR= (1.0 - f_rho)*Vn_Mag + f_rho*fabs(ProjVelocity_j);  
  
  /*--- Mass flux function ---*/

  mF = 0.5 * (Density_i * (ProjVelocity_i + Vn_MagL) + Density_j * (ProjVelocity_j - Vn_MagR) - (Chi/aF)*(Pressure_j-Pressure_i));
  
  /*--- Pressure function ---*/

  if (fabs(mL) < 1.0) BetaL = 0.25*(2.0-mL)*pow((mL+1.0),2.0);
  else {
      if (mL >= 0) BetaL = 1.0;
      else BetaL = 0.0;
  }
    
  if (fabs(mR) < 1.0) BetaR = 0.25*(2.0+mR)*pow((mR-1.0),2.0);
  else {
      if (mR >= 0 ) BetaR = 0.0;
      else BetaR = 1.0;
  }
    
  if (slau_low_diss){
    SetRoe_Dissipation(Dissipation_i, Dissipation_j, Sensor_i, Sensor_j, Dissipation_ij, config);
  }
  
  pF = 0.5 * (Pressure_i + Pressure_j) + 0.5 * (BetaL - BetaR) * (Pressure_i - Pressure_j) + Dissipation_ij*(1.0 - Chi) * (BetaL + BetaR - 1.0) *  0.5 * (Pressure_i + Pressure_j);
  
  val_residual[0] = 0.5*(mF+fabs(mF)) + 0.5*(mF-fabs(mF));
  
  for (iDim = 0; iDim < nDim; iDim++) {
    val_residual[iDim+1] = 0.5*(mF+fabs(mF)) * Velocity_i[iDim];
    val_residual[iDim+1]+= 0.5*(mF-fabs(mF)) * Velocity_j[iDim] ;
    val_residual[iDim+1]+= pF*UnitNormal[iDim];
  }

  val_residual[nVar-1] = 0.5*(mF+fabs(mF))*(Enthalpy_i) + 0.5*(mF-fabs(mF))*(Enthalpy_j);

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] *= Area;
  
  /*--- Roe's Jacobian for AUSM (this must be fixed) ---*/
  if (implicit) {
    
    /*--- Mean Roe variables iPoint and jPoint ---*/
    R = sqrt(fabs(Density_j/Density_i));
    RoeDensity = R*Density_i;
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
      sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
    }
    RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
    RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
    
    /*--- Compute P and Lambda (do it with the Normal) ---*/
    GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
    
    ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
      ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
    }
    
    /*--- Flow eigenvalues and Entropy correctors ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Lambda[iDim] = ProjVelocity;
    Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
    Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
    
    /*--- Compute inverse P ---*/
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobias of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
    
    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }
}

CUpwSLAU2_Flow::CUpwSLAU2_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config, bool val_low_dissipation) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  slau_low_dissipation = val_low_dissipation;
  
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  delta_vel  = new su2double [nDim];
  delta_wave = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwSLAU2_Flow::~CUpwSLAU2_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwSLAU2_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
   
  /*--- Face area (norm or the normal vector) ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel)));
  
  /*--- Primitive variables at point j ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
  }
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel)));
  
  /*--- Projected velocities ---*/
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
    
  aF = 0.5 * (SoundSpeed_i + SoundSpeed_j);
  
  mL  = ProjVelocity_i/aF;
  mR  = ProjVelocity_j/aF;
 
  /*--- Smooth function of the local Mach number---*/
  aux_slau = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    aux_slau += Velocity_i[iDim]*Velocity_i[iDim];
    aux_slau += Velocity_j[iDim]*Velocity_j[iDim];
  }
  Mach_tilde = min(1.0, (1.0/aF) * sqrt(aux_slau/2.0));  
  Chi = pow((1.0 - Mach_tilde),2.0);
  f_rho = -max(min(mL,0.0),-1.0) * min(max(mR,0.0),1.0);

  /*--- Mean normal velocity with density weighting ---*/
  Vn_Mag = (Density_i*fabs(ProjVelocity_i) + Density_j*fabs(ProjVelocity_j)) / (Density_i + Density_j);
  Vn_MagL= (1.0 - f_rho)*Vn_Mag + f_rho*fabs(ProjVelocity_i);
  Vn_MagR= (1.0 - f_rho)*Vn_Mag + f_rho*fabs(ProjVelocity_j);  
  
  /*--- Mass flux function ---*/
  mF = 0.5 * (Density_i * (ProjVelocity_i + Vn_MagL) + Density_j * (ProjVelocity_j - Vn_MagR) - (Chi/aF)*(Pressure_j-Pressure_i));
  
  /*--- Pressure function ---*/
  if (fabs(mL) < 1.0) BetaL = 0.25*(2.0-mL)*pow((mL+1.0),2.0);
  else {
    if (mL >= 0) BetaL = 1.0;
    else BetaL = 0.0;
  }
    
  if (fabs(mR) < 1.0) BetaR = 0.25*(2.0+mR)*pow((mR-1.0),2.0);
  else {
    if (mR >= 0 ) BetaR = 0.0;
    else BetaR = 1.0;
  }
  
  if (slau_low_dissipation){
    SetRoe_Dissipation(Dissipation_i, Dissipation_j, Sensor_i, Sensor_j, Dissipation_ij, config);
  }
  
  /*--- Pressure Flux ---*/
  
  pF = 0.5 * (Pressure_i + Pressure_j) + 0.5 * (BetaL - BetaR) * (Pressure_i - Pressure_j) + Dissipation_ij * sqrt(aux_slau/2.0) * (BetaL + BetaR - 1.0) * aF * 0.5 * (Density_i + Density_j);
  
  val_residual[0] = 0.5*(mF+fabs(mF)) + 0.5*(mF-fabs(mF));
  
  for (iDim = 0; iDim < nDim; iDim++) {
    val_residual[iDim+1] = 0.5*(mF+fabs(mF)) * Velocity_i[iDim];
    val_residual[iDim+1]+= 0.5*(mF-fabs(mF)) * Velocity_j[iDim] ;
    val_residual[iDim+1]+= pF*UnitNormal[iDim];
  }

  val_residual[nVar-1] = 0.5*(mF+fabs(mF))*(Enthalpy_i) + 0.5*(mF-fabs(mF))*(Enthalpy_j);

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] *= Area;
  
  /*--- Roe's Jacobian for AUSM (this must be fixed) ---*/
  if (implicit) {
    
    /*--- Mean Roe variables iPoint and jPoint ---*/
    R = sqrt(fabs(Density_j/Density_i));
    RoeDensity = R*Density_i;
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
      sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
    }
    RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
    RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
    
    /*--- Compute P and Lambda (do it with the Normal) ---*/
    GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
    
    ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
      ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
    }
    
    /*--- Flow eigenvalues and Entropy correctors ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Lambda[iDim] = ProjVelocity;
    Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
    Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
    
    /*--- Compute inverse P ---*/
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobias of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
    
    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }
}

CUpwHLLC_Flow::CUpwHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  kappa = config->GetRoe_Kappa();
  grid_movement = config->GetGrid_Movement();
  
  Gamma = config->GetGamma();

  Gamma_Minus_One = Gamma - 1.0;
  
  IntermediateState = new su2double [nVar];
  dSm_dU            = new su2double [nVar];
  dPI_dU            = new su2double [nVar];
  drhoStar_dU       = new su2double [nVar];
  dpStar_dU         = new su2double [nVar];
  dEStar_dU         = new su2double [nVar];

  Velocity_i        = new su2double [nDim];
  Velocity_j        = new su2double [nDim];
  RoeVelocity       = new su2double [nDim];  
  
}

CUpwHLLC_Flow::~CUpwHLLC_Flow(void) {
  
  delete [] IntermediateState;
  delete [] dSm_dU;
  delete [] dPI_dU;
  delete [] drhoStar_dU;
  delete [] dpStar_dU;
  delete [] dEStar_dU;

  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  
}

void CUpwHLLC_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  /*--- Face area (norm or the normal vector) ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim] * Normal[iDim];

  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim] / Area;

  /*-- Fluid velocity at node i,j ---*/  

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim]  = V_i[iDim+1];
    Velocity_j[iDim]  = V_j[iDim+1];
  }

  /*--- Primitive variables at point i ---*/

  Pressure_i = V_i[nDim+1];
  Density_i  = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];

  /*--- Primitive variables at point j ---*/
  
  Pressure_j = V_j[nDim+1];
  Density_j  = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];


  sq_vel_i = 0.0;
  sq_vel_j = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    sq_vel_i += Velocity_i[iDim] * Velocity_i[iDim];
    sq_vel_j += Velocity_j[iDim] * Velocity_j[iDim];
  }

  Energy_i = Enthalpy_i - Pressure_i / Density_i;
  Energy_j = Enthalpy_j - Pressure_j / Density_j;

  SoundSpeed_i = sqrt( (Enthalpy_i - 0.5 * sq_vel_i) * Gamma_Minus_One );
  SoundSpeed_j = sqrt( (Enthalpy_j - 0.5 * sq_vel_j) * Gamma_Minus_One );
   
  /*--- Projected velocities ---*/
  
  ProjVelocity_i = 0; 
  ProjVelocity_j = 0;

  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim] * UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim] * UnitNormal[iDim];
  }

  /*--- Projected Grid Velocity ---*/

  ProjInterfaceVel = 0;

  if (grid_movement) {

  for (iDim = 0; iDim < nDim; iDim++)
    ProjInterfaceVel += 0.5 * ( GridVel_i[iDim] + GridVel_j[iDim] )*UnitNormal[iDim];

  SoundSpeed_i -= ProjInterfaceVel;
    SoundSpeed_j += ProjInterfaceVel;

        ProjVelocity_i -= ProjInterfaceVel; 
        ProjVelocity_j -= ProjInterfaceVel;
  }  
  
  /*--- Roe's averaging ---*/

  Rrho = ( sqrt(Density_i) + sqrt(Density_j) );

  sq_velRoe        = 0.0;
  RoeProjVelocity  = - ProjInterfaceVel;

  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = ( Velocity_i[iDim] * sqrt(Density_i) + Velocity_j[iDim] * sqrt(Density_j) ) / Rrho;
    sq_velRoe        +=  RoeVelocity[iDim] * RoeVelocity[iDim];
    RoeProjVelocity  +=  RoeVelocity[iDim] * UnitNormal[iDim];
  } 

  /*--- Mean Roe variables iPoint and jPoint ---*/
    
  RoeDensity = sqrt( Density_i * Density_j );
  RoeEnthalpy = ( sqrt(Density_j) * Enthalpy_j + sqrt(Density_i) * Enthalpy_i) / Rrho;

  /*--- Roe-averaged speed of sound ---*/

  //RoeSoundSpeed2 = Gamma_Minus_One * ( RoeEnthalpy - 0.5 * sq_velRoe );
  RoeSoundSpeed  = sqrt( Gamma_Minus_One * ( RoeEnthalpy - 0.5 * sq_velRoe  ) ) - ProjInterfaceVel;


  /*--- Speed of sound at L and R ---*/
  
  sL = min( RoeProjVelocity - RoeSoundSpeed, ProjVelocity_i - SoundSpeed_i);
  sR = max( RoeProjVelocity + RoeSoundSpeed, ProjVelocity_j + SoundSpeed_j);
  
  /*--- speed of contact surface ---*/

  RHO = Density_j * (sR - ProjVelocity_j) - Density_i * (sL - ProjVelocity_i);
  sM = ( Pressure_i - Pressure_j - Density_i * ProjVelocity_i * ( sL - ProjVelocity_i ) + Density_j * ProjVelocity_j * ( sR - ProjVelocity_j ) ) / RHO;
  
  /*--- Pressure at right and left (Pressure_j=Pressure_i) side of contact surface ---*/
  
  pStar = Density_j * ( ProjVelocity_j - sR ) * ( ProjVelocity_j - sM ) + Pressure_j;


if (sM > 0.0) {

  if (sL > 0.0) {

    /*--- Compute Left Flux ---*/

    val_residual[0] = Density_i * ProjVelocity_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = Density_i * Velocity_i[iDim] * ProjVelocity_i + Pressure_i * UnitNormal[iDim];
    val_residual[nVar-1] = Enthalpy_i * Density_i * ProjVelocity_i;
  }
  else {

    /*--- Compute Flux Left Star from Left Star State ---*/

                rhoSL = ( sL - ProjVelocity_i ) / ( sL - sM );

    IntermediateState[0] = rhoSL * Density_i;
    for (iDim = 0; iDim < nDim; iDim++)
      IntermediateState[iDim+1] = rhoSL * ( Density_i * Velocity_i[iDim] + ( pStar - Pressure_i ) / ( sL - ProjVelocity_i ) * UnitNormal[iDim] ) ;
    IntermediateState[nVar-1] = rhoSL * ( Density_i * Energy_i - ( Pressure_i * ProjVelocity_i - pStar * sM) / ( sL - ProjVelocity_i ) );


    val_residual[0] = sM * IntermediateState[0];
    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = sM * IntermediateState[iDim+1] + pStar * UnitNormal[iDim];
    val_residual[nVar-1] = sM * ( IntermediateState[nVar-1] + pStar ) + pStar * ProjInterfaceVel;
  }
  }
  else {

  if (sR < 0.0) {

    /*--- Compute Right Flux ---*/

    val_residual[0] = Density_j * ProjVelocity_j;  
    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = Density_j * Velocity_j[iDim] * ProjVelocity_j + Pressure_j * UnitNormal[iDim];
    val_residual[nVar-1] = Enthalpy_j * Density_j * ProjVelocity_j;
  }
  else {

    /*--- Compute Flux Right Star from Right Star State ---*/

                rhoSR = ( sR - ProjVelocity_j ) / ( sR - sM );

    IntermediateState[0] = rhoSR * Density_j;
    for (iDim = 0; iDim < nDim; iDim++)
      IntermediateState[iDim+1] = rhoSR * ( Density_j * Velocity_j[iDim] + ( pStar - Pressure_j ) / ( sR - ProjVelocity_j ) * UnitNormal[iDim] ) ;
    IntermediateState[nVar-1] = rhoSR * ( Density_j * Energy_j - ( Pressure_j * ProjVelocity_j - pStar * sM ) / ( sR - ProjVelocity_j ) );


    val_residual[0] = sM * IntermediateState[0];
    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = sM * IntermediateState[iDim+1] + pStar * UnitNormal[iDim];
    val_residual[nVar-1] = sM * (IntermediateState[nVar-1] + pStar ) + pStar * ProjInterfaceVel;
  }
  }


  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] *= Area;


  if (implicit) {

  if (sM > 0.0) {

    if (sL > 0.0) {

      /*--- Compute Jacobian based on Left State ---*/
  
      for (iVar = 0; iVar < nVar; iVar++) 
        for (jVar = 0; jVar < nVar; jVar++) 
          val_Jacobian_j[iVar][jVar] = 0;

      GetInviscidProjJac(Velocity_i, &Energy_i, UnitNormal, 1.0, val_Jacobian_i);

    }
    else {
      /*--- Compute Jacobian based on Left Star State ---*/

      EStar = IntermediateState[nVar-1];
      Omega = 1/(sL-sM);
      OmegaSM = Omega * sM;


      /*--------- Left Jacobian ---------*/


      /*--- Computing pressure derivatives d/dU_L (PI) ---*/

      dPI_dU[0] = 0.5 * Gamma_Minus_One * sq_vel_i;
      for (iDim = 0; iDim < nDim; iDim++)      
        dPI_dU[iDim+1] = - Gamma_Minus_One * Velocity_i[iDim];
      dPI_dU[nVar-1] = Gamma_Minus_One;
      

      /*--- Computing d/dU_L (Sm) ---*/

      dSm_dU[0] = ( - ProjVelocity_i * ProjVelocity_i + sM * sL + dPI_dU[0] ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = ( UnitNormal[iDim] * ( 2 * ProjVelocity_i - sL - sM ) + dPI_dU[iDim+1] ) / RHO;
      dSm_dU[nVar-1] = dPI_dU[nVar-1] / RHO;

      
      /*--- Computing d/dU_L (rhoStar) ---*/

      drhoStar_dU[0] = Omega * ( sL + IntermediateState[0] * dSm_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        drhoStar_dU[iDim+1] = Omega * ( - UnitNormal[iDim] + IntermediateState[0] * dSm_dU[iDim+1] );
      drhoStar_dU[nVar-1] = Omega * IntermediateState[0] * dSm_dU[nVar-1];
      

      /*--- Computing d/dU_L (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_i * (sR - ProjVelocity_j) * dSm_dU[iVar];


      /*--- Computing d/dU_L (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );
      
      dEStar_dU[0] += Omega * ProjVelocity_i * ( Enthalpy_i - dPI_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        dEStar_dU[iDim+1] += Omega * ( - UnitNormal[iDim] * Enthalpy_i - ProjVelocity_i * dPI_dU[iDim+1] );
      dEStar_dU[nVar-1] += Omega * ( sL - ProjVelocity_i - ProjVelocity_i * dPI_dU[nVar-1] );



      /*--- Jacobian First Row ---*/
            
      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_i[0][iVar] = sM * drhoStar_dU[iVar] + IntermediateState[0] * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (jDim = 0; jDim < nDim; jDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_i[jDim+1][iVar] = ( OmegaSM + 1 ) * ( UnitNormal[jDim] * dpStar_dU[iVar] + IntermediateState[jDim+1] * dSm_dU[iVar] );

        val_Jacobian_i[jDim+1][0] += OmegaSM * Velocity_i[jDim] * ProjVelocity_i;

        val_Jacobian_i[jDim+1][jDim+1] += OmegaSM * (sL - ProjVelocity_i);
        
        for (iDim = 0; iDim < nDim; iDim++)
          val_Jacobian_i[jDim+1][iDim+1] -= OmegaSM * Velocity_i[jDim] * UnitNormal[iDim];

        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_i[jDim+1][iVar] -= OmegaSM * dPI_dU[iVar] * UnitNormal[jDim];
      }

      /*--- Jacobian Last Row ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_i[nVar-1][iVar] = sM * ( dEStar_dU[iVar] + dpStar_dU[iVar] ) + ( EStar + pStar ) * dSm_dU[iVar];




      /*--------- Right Jacobian ---------*/


      /*--- Computing d/dU_R (Sm) ---*/
      
      dSm_dU[0] = ( ProjVelocity_j * ProjVelocity_j - sM * sR - 0.5 * Gamma_Minus_One * sq_vel_j ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = - ( UnitNormal[iDim] * ( 2 * ProjVelocity_j - sR - sM) - Gamma_Minus_One * Velocity_j[iDim] ) / RHO;
      dSm_dU[nVar-1]  = - Gamma_Minus_One / RHO;
      
      
      /*--- Computing d/dU_R (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_j * (sL - ProjVelocity_i) * dSm_dU[iVar];


      /*--- Computing d/dU_R (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );
      


      /*--- Jacobian First Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_j[0][iVar] = IntermediateState[0] * ( OmegaSM + 1 ) * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_j[iDim+1][iVar] = ( OmegaSM + 1 ) * ( IntermediateState[iDim+1] * dSm_dU[iVar] + UnitNormal[iDim] * dpStar_dU[iVar] );
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_j[nVar-1][iVar] = sM * (dEStar_dU[iVar] + dpStar_dU[iVar]) + (EStar + pStar) * dSm_dU[iVar];
    }
  }
  else {
    if (sR < 0.0) {

      /*--- Compute Jacobian based on Right State ---*/
  
      for (iVar = 0; iVar < nVar; iVar++) 
        for (jVar = 0; jVar < nVar; jVar++) 
          val_Jacobian_i[iVar][jVar] = 0;

      GetInviscidProjJac(Velocity_j, &Energy_j, UnitNormal, 1.0, val_Jacobian_j);
    
    }
    else {
      /*--- Compute Jacobian based on Right Star State ---*/

      EStar = IntermediateState[nVar-1];
      Omega = 1/(sR-sM);
      OmegaSM = Omega * sM;


      /*--------- Left Jacobian ---------*/


      /*--- Computing d/dU_L (Sm) ---*/

      dSm_dU[0] = ( - ProjVelocity_i * ProjVelocity_i + sM * sL + 0.5 * Gamma_Minus_One * sq_vel_i ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = ( UnitNormal[iDim] * ( 2 * ProjVelocity_i - sL - sM ) - Gamma_Minus_One * Velocity_i[iDim] ) / RHO;
      dSm_dU[nVar-1] = Gamma_Minus_One / RHO;
      

      /*--- Computing d/dU_L (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_i * (sR - ProjVelocity_j) * dSm_dU[iVar];


      /*--- Computing d/dU_L (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );
      


      /*--- Jacobian First Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_i[0][iVar] = IntermediateState[0] * ( OmegaSM + 1 ) * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_i[iDim+1][iVar] = (OmegaSM + 1) * ( IntermediateState[iDim+1] * dSm_dU[iVar] + UnitNormal[iDim] * dpStar_dU[iVar] );
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_i[nVar-1][iVar] = sM * (dEStar_dU[iVar] + dpStar_dU[iVar]) + (EStar + pStar) * dSm_dU[iVar];



      /*--------- Right Jacobian ---------*/
      
      
      /*--- Computing pressure derivatives d/dU_R (PI) ---*/

      dPI_dU[0] = 0.5 * Gamma_Minus_One * sq_vel_j;
      for (iDim = 0; iDim < nDim; iDim++)      
        dPI_dU[iDim+1] = - Gamma_Minus_One * Velocity_j[iDim];
      dPI_dU[nVar-1] = Gamma_Minus_One;



      /*--- Computing d/dU_R (Sm) ---*/
      
      dSm_dU[0] = - ( - ProjVelocity_j * ProjVelocity_j + sM * sR + dPI_dU[0] ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = - ( UnitNormal[iDim] * ( 2 * ProjVelocity_j - sR - sM) + dPI_dU[iDim+1] ) / RHO;
      dSm_dU[nVar-1]  = - dPI_dU[nVar-1] / RHO;
      

      /*--- Computing d/dU_R (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_j * (sL - ProjVelocity_i) * dSm_dU[iVar];
      

      /*--- Computing d/dU_R (rhoStar) ---*/

      drhoStar_dU[0] = Omega * ( sR + IntermediateState[0] * dSm_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        drhoStar_dU[iDim+1] = Omega * ( - UnitNormal[iDim] + IntermediateState[0] * dSm_dU[iDim+1] );
      drhoStar_dU[nVar-1] = Omega * IntermediateState[0] * dSm_dU[nVar-1];


      /*--- Computing d/dU_R (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );
      
      dEStar_dU[0] += Omega * ProjVelocity_j * ( Enthalpy_j - dPI_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        dEStar_dU[iDim+1] += Omega * ( - UnitNormal[iDim] * Enthalpy_j - ProjVelocity_j * dPI_dU[iDim+1] );
      dEStar_dU[nVar-1] += Omega * ( sR - ProjVelocity_j - ProjVelocity_j * dPI_dU[nVar-1] );



      /*--- Jacobian First Row ---*/
            
      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_j[0][iVar] = sM * drhoStar_dU[iVar] + IntermediateState[0] * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (jDim = 0; jDim < nDim; jDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_j[jDim+1][iVar] = ( OmegaSM + 1 ) * ( UnitNormal[jDim] * dpStar_dU[iVar] + IntermediateState[jDim+1] * dSm_dU[iVar] );

        val_Jacobian_j[jDim+1][0] += OmegaSM * Velocity_j[jDim] * ProjVelocity_j;

        val_Jacobian_j[jDim+1][jDim+1] += OmegaSM * (sR - ProjVelocity_j);
        
        for (iDim = 0; iDim < nDim; iDim++)
          val_Jacobian_j[jDim+1][iDim+1] -= OmegaSM * Velocity_j[jDim] * UnitNormal[iDim];

        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_j[jDim+1][iVar] -= OmegaSM * dPI_dU[iVar] * UnitNormal[jDim];
      }

      /*--- Jacobian Last Row ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_j[nVar-1][iVar] = sM * ( dEStar_dU[iVar] + dpStar_dU[iVar] ) + ( EStar + pStar ) * dSm_dU[iVar];
  
    }
  }


  /*--- Jacobians of the inviscid flux, scale = k because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
  
  Area *= kappa;  

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_Jacobian_i[iVar][jVar] *=   Area;
      val_Jacobian_j[iVar][jVar] *=   Area;
    }
  }
}

}

CUpwGeneralHLLC_Flow::CUpwGeneralHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  kappa = config->GetRoe_Kappa();
  grid_movement = config->GetGrid_Movement();
  
  Gamma = config->GetGamma();
  
  IntermediateState = new su2double [nVar];
  dSm_dU            = new su2double [nVar];
  dPI_dU            = new su2double [nVar];
  drhoStar_dU       = new su2double [nVar];
  dpStar_dU         = new su2double [nVar];
  dEStar_dU         = new su2double [nVar];

  Velocity_i        = new su2double [nDim];
  Velocity_j        = new su2double [nDim];
  RoeVelocity       = new su2double [nDim];

}

CUpwGeneralHLLC_Flow::~CUpwGeneralHLLC_Flow(void) {
  
  delete [] IntermediateState;
  delete [] dSm_dU;
  delete [] dPI_dU;
  delete [] drhoStar_dU;
  delete [] dpStar_dU;
  delete [] dEStar_dU;

  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;

}

void CUpwGeneralHLLC_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  /*--- Face area (norm or the normal vector) ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim] * Normal[iDim];

  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim] / Area;

  /*-- Fluid velocity at node i,j ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim]  = V_i[iDim+1];
    Velocity_j[iDim]  = V_j[iDim+1];
  }

  /*--- Primitive variables at point i ---*/

  Pressure_i = V_i[nDim+1];
  Density_i  = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];

  /*--- Primitive variables at point j ---*/
  
  Pressure_j = V_j[nDim+1];
  Density_j  = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];


  sq_vel_i = 0.0;
  sq_vel_j = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    sq_vel_i += Velocity_i[iDim] * Velocity_i[iDim];
    sq_vel_j += Velocity_j[iDim] * Velocity_j[iDim];
  }

  Energy_i         = Enthalpy_i - Pressure_i / Density_i;
  StaticEnthalpy_i = Enthalpy_i - 0.5 * sq_vel_i;
  StaticEnergy_i   = Energy_i - 0.5 * sq_vel_i;
  
  Kappa_i = S_i[1] / Density_i;
  Chi_i   = S_i[0] - Kappa_i * StaticEnergy_i;
  SoundSpeed_i = sqrt(Chi_i + StaticEnthalpy_i * Kappa_i);
  

  Energy_j         = Enthalpy_j - Pressure_j / Density_j;
  StaticEnthalpy_j = Enthalpy_j - 0.5 * sq_vel_j;
  StaticEnergy_j   = Energy_j - 0.5 * sq_vel_j;

  Kappa_j = S_j[1] / Density_j;
  Chi_j   = S_j[0] - Kappa_j * StaticEnergy_j;
  SoundSpeed_j = sqrt(Chi_j + StaticEnthalpy_j * Kappa_j);
   
  /*--- Projected velocities ---*/
  
  ProjVelocity_i = 0.0; 
  ProjVelocity_j = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim] * UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim] * UnitNormal[iDim];
  }
  

  /*--- Projected Grid Velocity ---*/

  ProjInterfaceVel = 0;

  if (grid_movement) {

  for (iDim = 0; iDim < nDim; iDim++)
    ProjInterfaceVel += 0.5 * ( GridVel_i[iDim] + GridVel_j[iDim] )*UnitNormal[iDim];

  SoundSpeed_i -= ProjInterfaceVel;
    SoundSpeed_j += ProjInterfaceVel;

        ProjVelocity_i -= ProjInterfaceVel; 
        ProjVelocity_j -= ProjInterfaceVel;
  }  

  /*--- Roe's averaging ---*/

  Rrho = ( sqrt(Density_i) + sqrt(Density_j) );

  sq_velRoe        = 0.0;
  RoeProjVelocity  = - ProjInterfaceVel;

  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = ( Velocity_i[iDim] * sqrt(Density_i) + Velocity_j[iDim] * sqrt(Density_j) ) / Rrho;
    sq_velRoe        +=  RoeVelocity[iDim] * RoeVelocity[iDim];
    RoeProjVelocity  +=  RoeVelocity[iDim] * UnitNormal[iDim];
  } 

  /*--- Mean Roe variables iPoint and jPoint ---*/
    
  RoeKappa = 0.5 * ( Kappa_i + Kappa_j );
  RoeChi   = 0.5 * ( Chi_i + Chi_j );
  RoeDensity = sqrt( Density_i * Density_j );
  RoeEnthalpy = ( sqrt(Density_j) * Enthalpy_j + sqrt(Density_i) * Enthalpy_i) / Rrho;

  VinokurMontagne();

  /*--- Roe-averaged speed of sound ---*/

  //RoeSoundSpeed2 = RoeChi + RoeKappa * ( RoeEnthalpy - 0.5 * sq_velRoe );
  RoeSoundSpeed  = sqrt( RoeChi + RoeKappa * ( RoeEnthalpy - 0.5 * sq_velRoe ) ) - ProjInterfaceVel;

  /*--- Speed of sound at L and R ---*/
  
  sL = min( RoeProjVelocity - RoeSoundSpeed, ProjVelocity_i - SoundSpeed_i );
  sR = max( RoeProjVelocity + RoeSoundSpeed, ProjVelocity_j + SoundSpeed_j );
  
  /*--- speed of contact surface ---*/

  RHO = Density_j * (sR - ProjVelocity_j) - Density_i * (sL - ProjVelocity_i);
  sM = ( Pressure_i - Pressure_j - Density_i * ProjVelocity_i * ( sL - ProjVelocity_i ) + Density_j * ProjVelocity_j * ( sR - ProjVelocity_j ) ) / RHO;
  
  /*--- Pressure at right and left (Pressure_j=Pressure_i) side of contact surface ---*/
  
  pStar = Density_j * ( ProjVelocity_j - sR ) * ( ProjVelocity_j - sM ) + Pressure_j;


if (sM > 0.0) {

  if (sL > 0.0) {

    /*--- Compute Left Flux ---*/

    val_residual[0] = Density_i * ProjVelocity_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = Density_i * Velocity_i[iDim] * ProjVelocity_i + Pressure_i * UnitNormal[iDim];
    val_residual[nVar-1] = Enthalpy_i * Density_i * ProjVelocity_i;
  }
  else {

    /*--- Compute Flux Left Star from Left Star State ---*/

                rhoSL = ( sL - ProjVelocity_i ) / ( sL - sM );

    IntermediateState[0] = rhoSL * Density_i;
    for (iDim = 0; iDim < nDim; iDim++)
      IntermediateState[iDim+1] = rhoSL * ( Density_i * Velocity_i[iDim] + ( pStar - Pressure_i ) / ( sL - ProjVelocity_i ) * UnitNormal[iDim] ) ;
    IntermediateState[nVar-1] = rhoSL * ( Density_i * Energy_i - ( Pressure_i * ProjVelocity_i - pStar * sM) / ( sL - ProjVelocity_i ) );


    val_residual[0] = sM * IntermediateState[0];
    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = sM * IntermediateState[iDim+1] + pStar * UnitNormal[iDim];
    val_residual[nVar-1] = sM * ( IntermediateState[nVar-1] + pStar )  + pStar * ProjInterfaceVel;
  }
  }
  else {

  if (sR < 0.0) {

    /*--- Compute Right Flux ---*/

    val_residual[0] = Density_j * ProjVelocity_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = Density_j * Velocity_j[iDim] * ProjVelocity_j + Pressure_j * UnitNormal[iDim];
    val_residual[nVar-1] = Enthalpy_j * Density_j * ProjVelocity_j;
  }
  else {

    /*--- Compute Flux Right Star from Right Star State ---*/

                rhoSR = ( sR - ProjVelocity_j ) / ( sR - sM );

    IntermediateState[0] = rhoSR * Density_j;
    for (iDim = 0; iDim < nDim; iDim++)
      IntermediateState[iDim+1] = rhoSR * ( Density_j * Velocity_j[iDim] + ( pStar - Pressure_j ) / ( sR - ProjVelocity_j ) * UnitNormal[iDim] ) ;
    IntermediateState[nVar-1] = rhoSR * ( Density_j * Energy_j - ( Pressure_j * ProjVelocity_j - pStar * sM ) / ( sR - ProjVelocity_j ) );


    val_residual[0] = sM * IntermediateState[0];
    for (iDim = 0; iDim < nDim; iDim++)
      val_residual[iDim+1] = sM * IntermediateState[iDim+1] + pStar * UnitNormal[iDim];
    val_residual[nVar-1] = sM * (IntermediateState[nVar-1] + pStar )  + pStar * ProjInterfaceVel;
  }
  }

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] *= Area;


  if (implicit) {

  if (sM > 0.0) {

    if (sL > 0.0) {

      /*--- Compute Jacobian based on Left State ---*/
  
      for (iVar = 0; iVar < nVar; iVar++) 
        for (jVar = 0; jVar < nVar; jVar++) 
          val_Jacobian_j[iVar][jVar] = 0;


      GetInviscidProjJac(Velocity_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, 1.0, val_Jacobian_i);

    }
    else {
      /*--- Compute Jacobian based on Left Star State ---*/

      EStar = IntermediateState[nVar-1];
      Omega = 1/(sL-sM);
      OmegaSM = Omega * sM;


      /*--------- Left Jacobian ---------*/


      /*--- Computing pressure derivatives d/dU_L (PI) ---*/

      dPI_dU[0] = Chi_i - 0.5 * Kappa_i * sq_vel_i;
      for (iDim = 0; iDim < nDim; iDim++)      
        dPI_dU[iDim+1] = - Kappa_i * Velocity_i[iDim];
      dPI_dU[nVar-1] = Kappa_i;

      
      /*--- Computing d/dU_L (Sm) ---*/

      dSm_dU[0] = ( - ProjVelocity_i * ProjVelocity_i + sM * sL + dPI_dU[0] ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = ( UnitNormal[iDim] * ( 2 * ProjVelocity_i - sL - sM ) + dPI_dU[iDim+1] ) / RHO;
      dSm_dU[nVar-1] = dPI_dU[nVar-1] / RHO;


      /*--- Computing d/dU_L (rhoStar) ---*/

      drhoStar_dU[0] = Omega * ( sL + IntermediateState[0] * dSm_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        drhoStar_dU[iDim+1] = Omega * ( - UnitNormal[iDim] + IntermediateState[0] * dSm_dU[iDim+1] );
      drhoStar_dU[nVar-1] = Omega * IntermediateState[0] * dSm_dU[nVar-1];

      
      /*--- Computing d/dU_L (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_i * (sR - ProjVelocity_j) * dSm_dU[iVar];


      /*--- Computing d/dU_L (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );
      
      dEStar_dU[0] += Omega * ProjVelocity_i * ( Enthalpy_i - dPI_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        dEStar_dU[iDim+1] += Omega * ( - UnitNormal[iDim] * Enthalpy_i - ProjVelocity_i * dPI_dU[iDim+1] );
      dEStar_dU[nVar-1] += Omega * ( sL - ProjVelocity_i - ProjVelocity_i * dPI_dU[nVar-1] );



      /*--- Jacobian First Row ---*/
            
      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_i[0][iVar] = sM * drhoStar_dU[iVar] + IntermediateState[0] * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (jDim = 0; jDim < nDim; jDim++) {

        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_i[jDim+1][iVar] = ( OmegaSM + 1 ) * ( UnitNormal[jDim] * dpStar_dU[iVar] + IntermediateState[jDim+1] * dSm_dU[iVar] );

        val_Jacobian_i[jDim+1][0] += OmegaSM * Velocity_i[jDim] * ProjVelocity_i;

        val_Jacobian_i[jDim+1][jDim+1] += OmegaSM * (sL - ProjVelocity_i);
        
        for (iDim = 0; iDim < nDim; iDim++)
          val_Jacobian_i[jDim+1][iDim+1] -= OmegaSM * Velocity_i[jDim] * UnitNormal[iDim];

        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_i[jDim+1][iVar] -= OmegaSM * dPI_dU[iVar] * UnitNormal[jDim];
      }

      /*--- Jacobian Last Row ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_i[nVar-1][iVar] = sM * ( dEStar_dU[iVar] + dpStar_dU[iVar] ) + ( EStar + pStar ) * dSm_dU[iVar];




      /*--------- Right Jacobian ---------*/

      
      /*--- Computing pressure derivatives d/dU_R (PI) ---*/

      dPI_dU[0] = Chi_j - 0.5 * Kappa_j * sq_vel_j;
      for (iDim = 0; iDim < nDim; iDim++)      
        dPI_dU[iDim+1] = - Kappa_j * Velocity_j[iDim];
      dPI_dU[nVar-1] = Kappa_j;


      /*--- Computing d/dU_R (Sm) ---*/
      
      dSm_dU[0] = ( ProjVelocity_j * ProjVelocity_j - sM * sR - dPI_dU[0] ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = - ( UnitNormal[iDim] * ( 2 * ProjVelocity_j - sR - sM) + dPI_dU[iDim+1] ) / RHO;
      dSm_dU[nVar-1]  = - dPI_dU[nVar-1] / RHO;
      

      /*--- Computing d/dU_R (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_j * (sL - ProjVelocity_i) * dSm_dU[iVar];


      /*--- Computing d/dU_R (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );
      


      /*--- Jacobian First Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_j[0][iVar] = IntermediateState[0] * ( OmegaSM + 1 ) * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_j[iDim+1][iVar] = ( OmegaSM + 1 ) * ( IntermediateState[iDim+1] * dSm_dU[iVar] + UnitNormal[iDim] * dpStar_dU[iVar] );
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_j[nVar-1][iVar] = sM * (dEStar_dU[iVar] + dpStar_dU[iVar]) + (EStar + pStar) * dSm_dU[iVar];
    }
  }
  else {
    if (sR < 0.0) {

      /*--- Compute Jacobian based on Right State ---*/
  
      for (iVar = 0; iVar < nVar; iVar++) 
        for (jVar = 0; jVar < nVar; jVar++) 
          val_Jacobian_i[iVar][jVar] = 0;

      GetInviscidProjJac(Velocity_j, &Enthalpy_j, &Chi_j, &Kappa_j, UnitNormal, 1.0, val_Jacobian_j);
    
    }
    else {
      /*--- Compute Jacobian based on Right Star State ---*/

      EStar = IntermediateState[nVar-1];
      Omega = 1/(sR-sM);
      OmegaSM = Omega * sM;


      /*--------- Left Jacobian ---------*/


      /*--- Computing pressure derivatives d/dU_L (PI) ---*/

      dPI_dU[0] = Chi_i - 0.5 * Kappa_i * sq_vel_i;
      for (iDim = 0; iDim < nDim; iDim++)      
        dPI_dU[iDim+1] = - Kappa_i * Velocity_i[iDim];
      dPI_dU[nVar-1] = Kappa_i;


      /*--- Computing d/dU_L (Sm) ---*/

      dSm_dU[0] = ( - ProjVelocity_i * ProjVelocity_i + sM * sL + dPI_dU[0] ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = ( UnitNormal[iDim] * ( 2 * ProjVelocity_i - sL - sM ) + dPI_dU[iDim+1] ) / RHO;
      dSm_dU[nVar-1] = dPI_dU[nVar-1] / RHO;
      
      
      /*--- Computing d/dU_L (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_i * (sR - ProjVelocity_j) * dSm_dU[iVar];


      /*--- Computing d/dU_L (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );
      


      /*--- Jacobian First Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_i[0][iVar] = IntermediateState[0] * ( OmegaSM + 1 ) * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (iDim = 0; iDim < nDim; iDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_i[iDim+1][iVar] = (OmegaSM + 1) * ( IntermediateState[iDim+1] * dSm_dU[iVar] + UnitNormal[iDim] * dpStar_dU[iVar] );
      }

      /*--- Jacobian Last Row ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_i[nVar-1][iVar] = sM * (dEStar_dU[iVar] + dpStar_dU[iVar]) + (EStar + pStar) * dSm_dU[iVar];



      /*--------- Right Jacobian ---------*/

      
      /*--- Computing pressure derivatives d/dU_R (PI) ---*/

      dPI_dU[0] = Chi_j - 0.5 * Kappa_j * sq_vel_j;
      for (iDim = 0; iDim < nDim; iDim++)      
        dPI_dU[iDim+1] = - Kappa_j * Velocity_j[iDim];
      dPI_dU[nVar-1] = Kappa_j;


      /*--- Computing d/dU_R (Sm) ---*/
      
      dSm_dU[0] = - ( - ProjVelocity_j * ProjVelocity_j + sM * sR + dPI_dU[0] ) / RHO;
      for (iDim = 0; iDim < nDim; iDim++)
        dSm_dU[iDim+1] = - ( UnitNormal[iDim] * ( 2 * ProjVelocity_j - sR - sM) + dPI_dU[iDim+1] ) / RHO;
      dSm_dU[nVar-1]  = - dPI_dU[nVar-1] / RHO;


      /*--- Computing d/dU_R (pStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dpStar_dU[iVar] = Density_j * (sL - ProjVelocity_i) * dSm_dU[iVar];

      
      /*--- Computing d/dU_R (rhoStar) ---*/

      drhoStar_dU[0] = Omega * ( sR + IntermediateState[0] * dSm_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        drhoStar_dU[iDim+1] = Omega * ( - UnitNormal[iDim] + IntermediateState[0] * dSm_dU[iDim+1] );
      drhoStar_dU[nVar-1] = Omega * IntermediateState[0] * dSm_dU[nVar-1];


      /*--- Computing d/dU_R (EStar) ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        dEStar_dU[iVar] = Omega * ( sM * dpStar_dU[iVar] + ( EStar + pStar ) * dSm_dU[iVar] );
      
      dEStar_dU[0] += Omega * ProjVelocity_j * ( Enthalpy_j - dPI_dU[0] );
      for (iDim = 0; iDim < nDim; iDim++)
        dEStar_dU[iDim+1] += Omega * ( - UnitNormal[iDim] * Enthalpy_j - ProjVelocity_j * dPI_dU[iDim+1] );
      dEStar_dU[nVar-1] += Omega * ( sR - ProjVelocity_j - ProjVelocity_j * dPI_dU[nVar-1] );



      /*--- Jacobian First Row ---*/
            
      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_j[0][iVar] = sM * drhoStar_dU[iVar] + IntermediateState[0] * dSm_dU[iVar];

      /*--- Jacobian Middle Rows ---*/

      for (jDim = 0; jDim < nDim; jDim++) {
        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_j[jDim+1][iVar] = ( OmegaSM + 1 ) * ( UnitNormal[jDim] * dpStar_dU[iVar] + IntermediateState[jDim+1] * dSm_dU[iVar] );

        val_Jacobian_j[jDim+1][0] += OmegaSM * Velocity_j[jDim] * ProjVelocity_j;

        val_Jacobian_j[jDim+1][jDim+1] += OmegaSM * (sR - ProjVelocity_j);
        
        for (iDim = 0; iDim < nDim; iDim++)
          val_Jacobian_j[jDim+1][iDim+1] -= OmegaSM * Velocity_j[jDim] * UnitNormal[iDim];

        for (iVar = 0; iVar < nVar; iVar++)
          val_Jacobian_j[jDim+1][iVar] -= OmegaSM * dPI_dU[iVar] * UnitNormal[jDim];
      }
      
      /*--- Jacobian Last Row ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        val_Jacobian_j[nVar-1][iVar] = sM * ( dEStar_dU[iVar] + dpStar_dU[iVar] ) + ( EStar + pStar ) * dSm_dU[iVar];  
    }
  }


  /*--- Jacobians of the inviscid flux, scale = kappa because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/

  Area *= kappa;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_Jacobian_i[iVar][jVar] *= Area;
      val_Jacobian_j[iVar][jVar] *= Area;
    }
  }

  }

}

void CUpwGeneralHLLC_Flow::VinokurMontagne() {

  su2double delta_rhoStaticEnergy, delta_rho, delta_p, err_P, s, D;

  delta_rho = Density_j - Density_i;
  delta_p   = Pressure_j - Pressure_i;

  RoeKappaStaticEnthalpy = 0.5 * ( StaticEnthalpy_i * Kappa_i + StaticEnthalpy_j * Kappa_j );

  s = RoeChi + RoeKappaStaticEnthalpy;

  D = s*s * delta_rho * delta_rho + delta_p * delta_p;

  delta_rhoStaticEnergy = Density_j * StaticEnergy_j - Density_i * StaticEnergy_i;

  err_P = delta_p - RoeChi * delta_rho - RoeKappa * delta_rhoStaticEnergy;

  if (abs((D - delta_p*err_P)/Density_i) > 1e-3 && abs(delta_rho/Density_i) > 1e-3 && s/Density_i > 1e-3) {

    RoeKappa = ( D * RoeKappa ) / ( D - delta_p * err_P );
    RoeChi   = ( D * RoeChi+ s*s * delta_rho * err_P ) / ( D - delta_p * err_P );

  }
}

#ifdef CHECK

int UgpWithCvCompFlow::calcEulerFluxMatrices_HLLC(su2double (*val_Jacobian_i)[5], su2double (*val_Jacobian_j)[5], su2double (*val_Jacobian_i_Scal)[6], su2double (*val_Jacobian_j_Scal)[6],
                                                  const su2double Density_i, const su2double *uL, const su2double pL, const su2double TL, const su2double h0, const su2double RL, const su2double gammaL, const su2double *scalL, const su2double kL,
                                                  const su2double Density_j, const su2double *uR, const su2double pR, const su2double TR, const su2double h1, const su2double RR, const su2double gammaR, const su2double *scalR, const su2double kR,
                                                  const su2double area, const su2double *nVec, const int nScal, const su2double surfVeloc)
{

  su2double unL  = vecDotVec3d(uL, nVec);
  su2double uLuL = vecDotVec3d(uL, uL);
  su2double cL   = sqrt(gammaL*pL/Density_i);
  su2double hL   = gammaL/(gammaL-1.0)*pL/Density_i + 0.5*uLuL + kL;
  //  su2double hL   = h0 + 0.5*uLuL + kL;
  su2double eL   = hL*Density_i-pL;
  
  su2double unR  = vecDotVec3d(uR, nVec);
  su2double uRuR = vecDotVec3d(uR, uR);
  su2double cR   = sqrt(gammaR*pR/Density_j);
  su2double hR   = gammaR/(gammaR-1.0)*pR/Density_j + 0.5*uRuR + kR;
  //  su2double hR   = h1 + 0.5*uRuR + kR;
  su2double eR   = hR*Density_j-pR;
  
  
  // Roe's aveaging
  su2double Rrho = sqrt(Density_j/Density_i);
  su2double tmp = 1.0/(1.0+Rrho);
  su2double velRoe[3];
  for (int i=0; i<3; i++)
    velRoe[i] = tmp*(uL[i] + uR[i]*Rrho);
  su2double uRoe  = vecDotVec3d(velRoe, nVec);
  su2double hRoe = tmp*(hL + hR*Rrho);
  
  //  su2double cRoe  = sqrt((gammaL-1.0)*(hRoe- 0.5*vecDotVec3d(velRoe, velRoe)));
  su2double gamPdivRho = tmp*((gammaL*pL/Density_i+0.5*(gammaL-1.0)*uLuL) + (gammaR*pR/Density_j+0.5*(gammaR-1.0)*uRuR)*Rrho);
  su2double cRoe  = sqrt(gamPdivRho - ((gammaL+gammaR)*0.5-1.0)*0.5*vecDotVec3d(velRoe, velRoe));
  
  // speed of sound at L and R
  su2double sL = min(uRoe-cRoe, unL-cL);
  su2double sR = max(uRoe+cRoe, unR+cR);
  
  // speed of contact surface
  su2double sM = (pL-pR-Density_i*unL*(sL-unL)+Density_j*unR*(sR-unR))/(Density_j*(sR-unR)-Density_i*(sL-unL));
  
  // pressure at right and left (pR=pL) side of contact surface
  su2double pStar = Density_j*(unR-sR)*(unR-sM)+pR;
  
  if (sM >= 0.0) {
    
    if (sL > 0.0) {
      
      su2double nVecArea[3];
      for (int i=0; i<3; i++) nVecArea[i] = nVec[i]*area;
      
      calcJacobianA(val_Jacobian_i, uL, pL, Density_i, nVecArea, 0.5*(gammaL+gammaR), 0.0);
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_j[i][j] = 0.0;
      
    }
    else {
      
      su2double invSLmSs = 1.0/(sL-sM);
      su2double sLmuL = sL-unL;
      su2double rhoSL = Density_i*sLmuL*invSLmSs;
      su2double rhouSL[3];
      
      for (int i=0; i<3; i++)
        rhouSL[i] = (Density_i*uL[i]*sLmuL+(pStar-pL)*nVec[i])*invSLmSs;
      
      su2double eSL = (sLmuL*eL-pL*unL+pStar*sM)*invSLmSs;
      su2double gammaLM1 = (gammaL-1.0);
      su2double gammaRM1 = (gammaR-1.0);
      su2double invrhotld = 1.0/(Density_j*(sR-unR)-Density_i*(sL-unL));
      
      su2double dSMdUL[5], dSMdUR[5];
      su2double dpsdUL[5], dpsdUR[5];
      
      dSMdUL[0] = -unL*unL + uLuL*gammaLM1/2.0 + sM*sL;
      dSMdUL[1] =  nVec[0]*(2.0*unL-sL-sM) - gammaLM1*uL[0];
      dSMdUL[2] =  nVec[1]*(2.0*unL-sL-sM) - gammaLM1*uL[1];
      dSMdUL[3] =  nVec[2]*(2.0*unL-sL-sM) - gammaLM1*uL[2];
      dSMdUL[4] =  gammaLM1;
      
      for (iVar = 0; iVar < nVar; iVar++)
      {
        dSMdUL[i] *= invrhotld;
        dpsdUL[i] = Density_j*(sR-unR)*dSMdUL[i];
      }
      
      dSMdUR[0] =  unR*unR - uRuR*gammaRM1/2.0 - sM*sR;
      dSMdUR[1] = -nVec[0]*(2.0*unR-sR-sM) + gammaRM1*uR[0];
      dSMdUR[2] = -nVec[1]*(2.0*unR-sR-sM) + gammaRM1*uR[1];
      dSMdUR[3] = -nVec[2]*(2.0*unR-sR-sM) + gammaRM1*uR[2];
      dSMdUR[4] = -gammaRM1;
      
      for (iVar = 0; iVar < nVar; iVar++)
      {
        dSMdUR[i] *= invrhotld;
        dpsdUR[i] = Density_i*(sL-unL)*dSMdUR[i];
      }
      
      calcSubSonicJacobeanHLLC(val_Jacobian_i, val_Jacobian_j,
                               Density_i, uL, pL, eL, unL, uLuL, sL,
                               rhoSL, rhouSL, eSL, dSMdUL,
                               dSMdUR, dpsdUL, dpsdUR, sM, pStar, 0.5*(gammaL+gammaR), nVec);
      
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[0][i] =  val_Jacobian_i[0][i]*sM + dSMdUL[i]*rhoSL;
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[1][i] =  val_Jacobian_i[1][i]*sM + dSMdUL[i]*rhouSL[0] + dpsdUL[i]*nVec[0];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[2][i] =  val_Jacobian_i[2][i]*sM + dSMdUL[i]*rhouSL[1] + dpsdUL[i]*nVec[1];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[3][i] =  val_Jacobian_i[3][i]*sM + dSMdUL[i]*rhouSL[2] + dpsdUL[i]*nVec[2];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[4][i] = (val_Jacobian_i[4][i]+dpsdUL[i])*sM + (eSL+pStar)*dSMdUL[i];
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_i[i][j] *= area;
      
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[0][i] =  val_Jacobian_j[0][i]*sM + dSMdUR[i]*rhoSL;
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[1][i] =  val_Jacobian_j[1][i]*sM + dSMdUR[i]*rhouSL[0] + dpsdUR[i]*nVec[0];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[2][i] =  val_Jacobian_j[2][i]*sM + dSMdUR[i]*rhouSL[1] + dpsdUR[i]*nVec[1];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[3][i] =  val_Jacobian_j[3][i]*sM + dSMdUR[i]*rhouSL[2] + dpsdUR[i]*nVec[2];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[4][i] = (val_Jacobian_j[4][i]+dpsdUR[i])*sM + (eSL+pStar)*dSMdUR[i];
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_j[i][j] *= area;
      
    }
  }
  
  else {
    
    if (sR >= 0.0) {
      
      su2double invSRmSs = 1.0/(sR-sM);
      su2double sRmuR = sR-unR;
      su2double rhoSR = Density_j*sRmuR*invSRmSs;
      su2double rhouSR[3];
      for (int i=0; i<3; i++)
        rhouSR[i] = (Density_j*uR[i]*sRmuR+(pStar-pR)*nVec[i])*invSRmSs;
      su2double eSR = (sRmuR*eR-pR*unR+pStar*sM)*invSRmSs;
      su2double gammaLM1 = (gammaL-1.0);
      su2double gammaRM1 = (gammaR-1.0);
      su2double invrhotld = 1.0/(Density_j*(sR-unR)-Density_i*(sL-unL));
      
      su2double dSMdUL[5], dSMdUR[5];
      su2double dpsdUL[5], dpsdUR[5];
      
      dSMdUL[0] = -unL*unL + uLuL*gammaLM1/2.0 + sM*sL;
      dSMdUL[1] =  nVec[0]*(2.0*unL-sL-sM) - gammaLM1*uL[0];
      dSMdUL[2] =  nVec[1]*(2.0*unL-sL-sM) - gammaLM1*uL[1];
      dSMdUL[3] =  nVec[2]*(2.0*unL-sL-sM) - gammaLM1*uL[2];
      dSMdUL[4] =  gammaLM1;
      
      for (iVar = 0; iVar < nVar; iVar++) {
        dSMdUL[i] *= invrhotld;
        dpsdUL[i] = Density_j*(sR-unR)*dSMdUL[i];
      }
      
      dSMdUR[0] =  unR*unR - uRuR*gammaRM1/2.0 - sM*sR;
      dSMdUR[1] = -nVec[0]*(2.0*unR-sR-sM) + gammaRM1*uR[0];
      dSMdUR[2] = -nVec[1]*(2.0*unR-sR-sM) + gammaRM1*uR[1];
      dSMdUR[3] = -nVec[2]*(2.0*unR-sR-sM) + gammaRM1*uR[2];
      dSMdUR[4] = -gammaRM1;
      
      for (iVar = 0; iVar < nVar; iVar++) {
        dSMdUR[i] *= invrhotld;
        dpsdUR[i] = Density_i*(sL-unL)*dSMdUR[i];
      }
      
      calcSubSonicJacobeanHLLC(val_Jacobian_j, val_Jacobian_i,
                               Density_j, uR, pR, eR, unR, uRuR, sR,
                               rhoSR, rhouSR, eSR,
                               dSMdUR, dSMdUL, dpsdUR, dpsdUL, sM, pStar, 0.5*(gammaL+gammaR), nVec);
      
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[0][i] =  val_Jacobian_i[0][i]*sM + dSMdUL[i]*rhoSR;
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[1][i] =  val_Jacobian_i[1][i]*sM + dSMdUL[i]*rhouSR[0] + dpsdUL[i]*nVec[0];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[2][i] =  val_Jacobian_i[2][i]*sM + dSMdUL[i]*rhouSR[1] + dpsdUL[i]*nVec[1];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[3][i] =  val_Jacobian_i[3][i]*sM + dSMdUL[i]*rhouSR[2] + dpsdUL[i]*nVec[2];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_i[4][i] = (val_Jacobian_i[4][i]+dpsdUL[i])*sM + (eSR+pStar)*dSMdUL[i];
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_i[i][j] *= area;
      
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[0][i] =  val_Jacobian_j[0][i]*sM + dSMdUR[i]*rhoSR;
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[1][i] =  val_Jacobian_j[1][i]*sM + dSMdUR[i]*rhouSR[0] + dpsdUR[i]*nVec[0];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[2][i] =  val_Jacobian_j[2][i]*sM + dSMdUR[i]*rhouSR[1] + dpsdUR[i]*nVec[1];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[3][i] =  val_Jacobian_j[3][i]*sM + dSMdUR[i]*rhouSR[2] + dpsdUR[i]*nVec[2];
      for (iVar = 0; iVar < nVar; iVar++)  val_Jacobian_j[4][i] = (val_Jacobian_j[4][i]+dpsdUR[i])*sM + (eSR+pStar)*dSMdUR[i];
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_j[i][j] *= area;
      
    }
    
    else {
      
      su2double nVecArea[3];
      for (int i=0; i<3; i++)        nVecArea[i] = nVec[i]*area;
      calcJacobianA(val_Jacobian_j, uR, pR, Density_j, nVecArea, 0.5*(gammaL+gammaR), 0.0);
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          val_Jacobian_i[i][j] = 0.0;
      
    }
    
  }
  
}

void UgpWithCvCompFlow::calcSubSonicJacobeanHLLC(su2double (*AL)[5], su2double (*AR)[5],
                                                 su2double Density_i, const su2double *uL, su2double pL, su2double eL, su2double qL, su2double psiL, su2double SL,
                                                 su2double rhoSL, su2double *rhouSL, su2double eSL,
                                                 su2double *dSMdUL, su2double *dSMdUR, su2double *dpsdUL, su2double *dpsdUR, su2double SM, su2double pS,
                                                 su2double gamma, const su2double *nV) // nV is not normalized
{
  
  su2double gammaMinus1 = (gamma-1.0);
  su2double omL = 1.0/(SL-SM);
  
  AL[0][0] =  SL    + rhoSL*dSMdUL[0];
  AL[0][1] = -nV[0] + rhoSL*dSMdUL[1];
  AL[0][2] = -nV[1] + rhoSL*dSMdUL[2];
  AL[0][3] = -nV[2] + rhoSL*dSMdUL[3];
  AL[0][4] =        + rhoSL*dSMdUL[4];
  
  AL[1][0] =    qL*uL[0]       - nV[0]*psiL*gammaMinus1/2.0   + nV[0]*dpsdUL[0] + rhouSL[0]*dSMdUL[0];
  AL[1][1] =  SL - qL          + nV[0]*(gamma-2.0)*uL[0]      + nV[0]*dpsdUL[1] + rhouSL[0]*dSMdUL[1];
  AL[1][2] =     - uL[0]*nV[1] + nV[0]*gammaMinus1*uL[1]      + nV[0]*dpsdUL[2] + rhouSL[0]*dSMdUL[2];
  AL[1][3] =     - uL[0]*nV[2] + nV[0]*gammaMinus1*uL[2]      + nV[0]*dpsdUL[3] + rhouSL[0]*dSMdUL[3];
  AL[1][4] = -gammaMinus1*nV[0]                               + nV[0]*dpsdUL[4] + rhouSL[0]*dSMdUL[4];
  
  AL[2][0] =    qL*uL[1]       - nV[1]*psiL*gammaMinus1/2.0   + nV[1]*dpsdUL[0] + rhouSL[1]*dSMdUL[0];
  AL[2][1] =     - uL[1]*nV[0] + nV[1]*gammaMinus1*uL[0]      + nV[1]*dpsdUL[1] + rhouSL[1]*dSMdUL[1];
  AL[2][2] =  SL - qL          + nV[1]*(gamma-2.0)*uL[1]      + nV[1]*dpsdUL[2] + rhouSL[1]*dSMdUL[2];
  AL[2][3] =     - uL[1]*nV[2] + nV[1]*gammaMinus1*uL[2]      + nV[1]*dpsdUL[3] + rhouSL[1]*dSMdUL[3];
  AL[2][4] = -gammaMinus1*nV[1]                               + nV[1]*dpsdUL[4] + rhouSL[1]*dSMdUL[4];
  
  AL[3][0] =    qL*uL[2]       - nV[2]*psiL*gammaMinus1/2.0   + nV[2]*dpsdUL[0] + rhouSL[2]*dSMdUL[0];
  AL[3][1] =     - uL[2]*nV[0] + nV[2]*gammaMinus1*uL[0]      + nV[2]*dpsdUL[1] + rhouSL[2]*dSMdUL[1];
  AL[3][2] =     - uL[2]*nV[1] + nV[2]*gammaMinus1*uL[1]      + nV[2]*dpsdUL[2] + rhouSL[2]*dSMdUL[2];
  AL[3][3] =  SL - qL          + nV[2]*(gamma-2.0)*uL[2]      + nV[2]*dpsdUL[3] + rhouSL[2]*dSMdUL[3];
  AL[3][4] = -gammaMinus1*nV[2]                               + nV[2]*dpsdUL[4] + rhouSL[2]*dSMdUL[4];
  
  AL[4][0] =      qL*(eL+pL)/Density_i - qL*psiL*(gamma-1.0)/2.0   + SM*dpsdUL[0] + (pS+eSL)*dSMdUL[0];
  AL[4][1] = - nV[0]*(eL+pL)/Density_i + gammaMinus1*uL[0]*qL      + SM*dpsdUL[1] + (pS+eSL)*dSMdUL[1];
  AL[4][2] = - nV[1]*(eL+pL)/Density_i + gammaMinus1*uL[1]*qL      + SM*dpsdUL[2] + (pS+eSL)*dSMdUL[2];
  AL[4][3] = - nV[2]*(eL+pL)/Density_i + gammaMinus1*uL[2]*qL      + SM*dpsdUL[3] + (pS+eSL)*dSMdUL[3];
  AL[4][4] =   SL-qL*gamma                                    + SM*dpsdUL[4] + (pS+eSL)*dSMdUL[4];
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      AL[i][j] *= omL;
  
  for (iVar = 0; iVar < nVar; iVar++)    AR[0][i] = omL*rhoSL*dSMdUR[i];
  for (iVar = 0; iVar < nVar; iVar++)    AR[1][i] = omL*(nV[0]*dpsdUR[i]+rhouSL[0]*dSMdUR[i]);
  for (iVar = 0; iVar < nVar; iVar++)    AR[2][i] = omL*(nV[1]*dpsdUR[i]+rhouSL[1]*dSMdUR[i]);
  for (iVar = 0; iVar < nVar; iVar++)    AR[3][i] = omL*(nV[2]*dpsdUR[i]+rhouSL[2]*dSMdUR[i]);
  for (iVar = 0; iVar < nVar; iVar++)    AR[4][i] = omL*(dpsdUR[i]*SM+(pS+eSL)*dSMdUR[i]);
  
}

void UgpWithCvCompFlow::calcJacobianA(su2double (*A)[5], const su2double *vel, su2double pp, su2double rrho, const su2double *nV, su2double gamma, su2double surfVeloc) // nV is not normalized
{
 
  su2double kapm1 = (gamma - 1.0);
  
  su2double nVel[3];
  nVel[0] = vel[0]*nV[0];
  nVel[1] = vel[1]*nV[1];
  nVel[2] = vel[2]*nV[2];
  su2double U_k = nVel[0]+nVel[1]+nVel[2];
  su2double vSquHlf = 0.5*vecDotVec3d(vel, vel);
  su2double c = sqrt(gamma*pp/rrho);
  su2double inv_kap_m1 = 1.0/kapm1;
  
  A[0][0] =-surfVeloc;
  A[0][1] = nV[0];
  A[0][2] = nV[1];
  A[0][3] = nV[2];
  A[0][4] = 0.0;
  
  A[1][0] = -vel[0]*(nVel[1]+nVel[2])+nV[0]*(kapm1*vSquHlf-vel[0]*vel[0]);
  A[1][1] = (2.-gamma)*nVel[0]+U_k-surfVeloc;
  A[1][2] = vel[0]*nV[1]-kapm1*vel[1]*nV[0];
  A[1][3] = vel[0]*nV[2]-kapm1*vel[2]*nV[0];
  A[1][4] = kapm1*nV[0];
  
  A[2][0] = -vel[1]*(nVel[0]+nVel[2])+nV[1]*(kapm1*vSquHlf-vel[1]*vel[1]);
  A[2][1] = -kapm1*vel[0]*nV[1]+ vel[1]*nV[0];
  A[2][2] = (2.-gamma)*nVel[1]+U_k-surfVeloc;
  A[2][3] = vel[1]*nV[2]-kapm1*vel[2]*nV[1];
  A[2][4] = kapm1*nV[1];
  
  A[3][0] = -vel[2]*(nVel[0]+nVel[1])+nV[2]*(kapm1*vSquHlf-vel[2]*vel[2]);
  A[3][1] = -kapm1*vel[0]*nV[2]+vel[2]*nV[0];
  A[3][2] = -kapm1*vel[1]*nV[2]+vel[2]*nV[1];
  A[3][3] = (2.-gamma)*nVel[2]+U_k-surfVeloc;
  A[3][4] = kapm1*nV[2];
  
  A[4][0] = U_k*((gamma-2.)*vSquHlf-c*c*inv_kap_m1);
  A[4][1] = c*c*inv_kap_m1*nV[0]-kapm1*vel[0]*(nVel[1]+nVel[2])-(kapm1*vel[0]*vel[0]-vSquHlf)*nV[0];
  A[4][2] = c*c*inv_kap_m1*nV[1]-kapm1*vel[1]*(nVel[0]+nVel[2])-(kapm1*vel[1]*vel[1]-vSquHlf)*nV[1];
  A[4][3] = c*c*inv_kap_m1*nV[2]-kapm1*vel[2]*(nVel[0]+nVel[1])-(kapm1*vel[2]*vel[2]-vSquHlf)*nV[2];
  A[4][4] = gamma*U_k-surfVeloc;
  
}


#endif


CUpwRoe_Flow::CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config, bool val_low_dissipation) :
                             CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  kappa = config->GetRoe_Kappa(); // 1 is unstable

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  roe_low_dissipation = val_low_dissipation;
  
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  delta_vel  = new su2double [nDim];
  delta_wave = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwRoe_Flow::~CUpwRoe_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwRoe_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};
  su2double ProjGridVel = 0.0;

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+4); AD::SetPreaccIn(V_j, nDim+4); AD::SetPreaccIn(Normal, nDim);
  if (grid_movement) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  if (roe_low_dissipation){
    AD::SetPreaccIn(Sensor_i); AD::SetPreaccIn(Sensor_j);
    AD::SetPreaccIn(Dissipation_i); AD::SetPreaccIn(Dissipation_j);
    AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  }
  
  /*--- Face area (norm or the normal vector) ---*/

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_i[iDim] = V_i[iDim+1];
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Pressure_i*Gamma/Density_i));
 
  /*--- Primitive variables at point j ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_j[iDim] = V_j[iDim+1];
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Pressure_j*Gamma/Density_j));

  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Roe-averaged variables at interface between i & j ---*/
  
  R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  
  RoeSoundSpeed2 = (Gamma-1)*(RoeEnthalpy-0.5*sq_vel);
  
  /*--- Negative RoeSoundSpeed2, the jump
   variables is too large, exit the subrotuine
   without computing the fluxes ---*/
  
  if (RoeSoundSpeed2 <= 0.0) {
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] = 0.0;
      if (implicit){
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    AD::SetPreaccOut(val_residual, nVar);
    AD::EndPreacc();
    return;
  }
  
  RoeSoundSpeed = sqrt(RoeSoundSpeed2);
  
  /*--- Compute ProjFlux_i ---*/
  
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i);
  
  /*--- Compute ProjFlux_j ---*/
  
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j);
  
  /*--- Compute P and Lambda (do it with the Normal) ---*/
  
  GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
  
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- Projected velocity adjustment due to mesh motion ---*/
  
  if (grid_movement) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    }
    ProjVelocity   -= ProjGridVel;
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }
  
  /*--- Flow eigenvalues and entropy correctors ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;
  
  Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
  Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
  
  /*--- Compute absolute value with Mavriplis' entropy correction ---*/
  
  MaxLambda = fabs(ProjVelocity) + RoeSoundSpeed;
  Delta = config->GetEntropyFix_Coeff();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Lambda[iVar] = max(fabs(Lambda[iVar]), Delta*MaxLambda);
  }
  
  /*--- Compute inverse P ---*/
  
  GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
  
  /*--- Jacobians of the inviscid flux, scaled by
   kappa because val_resconv ~ kappa*(fc_i+fc_j)*Normal ---*/
  if (implicit) {
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, kappa, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, kappa, val_Jacobian_j);
  }
  
  /*--- Diference variables iPoint and jPoint ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = U_j[iVar]-U_i[iVar];
  
  if (roe_low_dissipation)
    SetRoe_Dissipation(Dissipation_i, Dissipation_j, Sensor_i, Sensor_j, Dissipation_ij, config);
  
  /*--- Roe's Flux approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    
    val_residual[iVar] = kappa*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_ij = 0.0;
      
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];

        val_residual[iVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area*Dissipation_ij;
        if(implicit){
          val_Jacobian_i[iVar][jVar] += (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
          val_Jacobian_j[iVar][jVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
        }
    }
  }
  
  /*--- Jacobian contributions due to grid motion ---*/
  
  if (grid_movement) {
    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      
      /*--- Implicit terms ---*/
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}

CUpwGeneralRoe_Flow::CUpwGeneralRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  kappa = config->GetRoe_Kappa(); // 1 is unstable


  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  delta_vel  = new su2double [nDim];
  delta_wave = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];

  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwGeneralRoe_Flow::~CUpwGeneralRoe_Flow(void) {

  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;

}

void CUpwGeneralRoe_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+4); AD::SetPreaccIn(V_j, nDim+4); AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(S_i, 2); AD::SetPreaccIn(S_j, 2);
  if (grid_movement) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  /*--- Face area (norm or the normal vector) ---*/

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
  Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  /*-- Unit Normal ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Primitive variables at point i ---*/

  Velocity2_i = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
  }

  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  StaticEnthalpy_i = Enthalpy_i - 0.5*Velocity2_i;
  StaticEnergy_i = StaticEnthalpy_i - Pressure_i/Density_i;

  Kappa_i = S_i[1]/Density_i;
  Chi_i = S_i[0] - Kappa_i*StaticEnergy_i;
  SoundSpeed_i = sqrt(Chi_i + StaticEnthalpy_i*Kappa_i);

  /*--- Primitive variables at point j ---*/


  Velocity2_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    Velocity2_j += Velocity_j[iDim]*Velocity_j[iDim];
  }

  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;

  StaticEnthalpy_j = Enthalpy_j - 0.5*Velocity2_j;
  StaticEnergy_j = StaticEnthalpy_j - Pressure_j/Density_j;

  Kappa_j = S_j[1]/Density_j;
  Chi_j = S_j[0] - Kappa_j*StaticEnergy_j;
  SoundSpeed_j = sqrt(Chi_j + StaticEnthalpy_j*Kappa_j);

  /*--- Recompute conservative variables ---*/

  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;

//  /*--- Roe-averaged variables at interface between i & j ---*/

    ComputeRoeAverage();

    if (RoeSoundSpeed2 <= 0.0) {
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] = 0.0;
      for (jVar = 0; jVar < nVar; jVar++) {
      val_Jacobian_i[iVar][iVar] = 0.0;
      val_Jacobian_j[iVar][iVar] = 0.0;
      }
    }
      return;
    }

    RoeSoundSpeed = sqrt(RoeSoundSpeed2);

  /*--- Compute ProjFlux_i ---*/
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i);

  /*--- Compute ProjFlux_j ---*/
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j);

  /*--- Compute P and Lambda (do it with the Normal) ---*/

  GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, &RoeEnthalpy, &RoeChi, &RoeKappa, UnitNormal, P_Tensor);

  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }

  /*--- Projected velocity adjustment due to mesh motion ---*/
  if (grid_movement) {
    su2double ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    }
    ProjVelocity   -= ProjGridVel;
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }

  /*--- Flow eigenvalues and entropy correctors ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;

  Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
  Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;

  /*--- Compute absolute value with Mavriplis' entropy correction ---*/

  MaxLambda = fabs(ProjVelocity) + RoeSoundSpeed;
  Delta = config->GetEntropyFix_Coeff();

  for (iVar = 0; iVar < nVar; iVar++) {
    Lambda[iVar] = max(fabs(Lambda[iVar]), Delta*MaxLambda);
   }

//  /*--- Harten and Hyman (1983) entropy correction ---*/
//  for (iDim = 0; iDim < nDim; iDim++)
//    Epsilon[iDim] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i, ProjVelocity_j-Lambda[iDim]));
//
//  Epsilon[nVar-2] = 4.0*max(0.0, max(Lambda[nVar-2]-(ProjVelocity_i+SoundSpeed_i),(ProjVelocity_j+SoundSpeed_j)-Lambda[nVar-2]));
//  Epsilon[nVar-1] = 4.0*max(0.0, max(Lambda[nVar-1]-(ProjVelocity_i-SoundSpeed_i),(ProjVelocity_j-SoundSpeed_j)-Lambda[nVar-1]));
//
//  for (iVar = 0; iVar < nVar; iVar++)
//    if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
//      Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
//    else
//      Lambda[iVar] = fabs(Lambda[iVar]);

//  for (iVar = 0; iVar < nVar; iVar++)
//    Lambda[iVar] = fabs(Lambda[iVar]);

  if (!implicit) {

    /*--- Compute wave amplitudes (characteristics) ---*/
    proj_delta_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
      proj_delta_vel += delta_vel[iDim]*Normal[iDim];
    }
    delta_p = Pressure_j - Pressure_i;
    delta_rho = Density_j - Density_i;
    proj_delta_vel = proj_delta_vel/Area;

    if (nDim == 2) {
      delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
      delta_wave[1] = UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1];
      delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
      delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    } else {
      delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
      delta_wave[1] = UnitNormal[0]*delta_vel[2]-UnitNormal[2]*delta_vel[0];
      delta_wave[2] = UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1];
      delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
      delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    }

    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
      for (jVar = 0; jVar < nVar; jVar++)
        val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
    }

    /*--- Flux contribution due to grid motion ---*/
    if (grid_movement) {
      ProjVelocity = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      for (iVar = 0; iVar < nVar; iVar++) {
        val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      }
    }
  }
  else {

    /*--- Compute inverse P ---*/

    GetPMatrix_inv(invP_Tensor, &RoeDensity, RoeVelocity, &RoeSoundSpeed, &RoeChi , &RoeKappa, UnitNormal);

     /*--- Jacobians of the inviscid flux, scaled by
      kappa because val_resconv ~ kappa*(fc_i+fc_j)*Normal ---*/

    GetInviscidProjJac(Velocity_i, &Enthalpy_i, &Chi_i, &Kappa_i, Normal, kappa, val_Jacobian_i);

    GetInviscidProjJac(Velocity_j, &Enthalpy_j, &Chi_j, &Kappa_j, Normal, kappa, val_Jacobian_j);


    /*--- Diference variables iPoint and jPoint ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Diff_U[iVar] = U_j[iVar]-U_i[iVar];

    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] = kappa*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;

        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/

        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];

        val_residual[iVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area;
        val_Jacobian_i[iVar][jVar] += (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
      }
    }

    /*--- Jacobian contributions due to grid motion ---*/
    if (grid_movement) {
      ProjVelocity = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      for (iVar = 0; iVar < nVar; iVar++) {
        val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
        /*--- Implicit terms ---*/
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }

  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}


void CUpwGeneralRoe_Flow::ComputeRoeAverage() {

  //su2double delta_rhoStaticEnergy, err_P, s, D;
  // su2double tol = 10-6;

  R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  sq_vel = 0;  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }

  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  delta_rho = Density_j - Density_i;
  delta_p = Pressure_j - Pressure_i;
  RoeKappa = 0.5*(Kappa_i + Kappa_j);
  RoeKappa = (Kappa_i + Kappa_j + 4*RoeKappa)/6;
  RoeChi = 0.5*(Chi_i + Chi_j);
  RoeChi = (Chi_i + Chi_j + 4*RoeChi)/6;


//  RoeKappaStaticEnthalpy = 0.5*(StaticEnthalpy_i*Kappa_i + StaticEnthalpy_j*Kappa_j);
//  RoeKappaStaticEnthalpy = (StaticEnthalpy_i*Kappa_i + StaticEnthalpy_j*Kappa_j + 4*RoeKappaStaticEnthalpy)/6;
//  s = RoeChi + RoeKappaStaticEnthalpy;
//  D = s*s*delta_rho*delta_rho + delta_p*delta_p;
//  delta_rhoStaticEnergy = Density_j*StaticEnergy_j - Density_i*StaticEnergy_i;
//  err_P = delta_p - RoeChi*delta_rho - RoeKappa*delta_rhoStaticEnergy;
//
//
//  if (abs((D - delta_p*err_P)/Density_i)>1e-3 && abs(delta_rho/Density_i)>1e-3 && s/Density_i > 1e-3) {
//
//    RoeKappa = (D*RoeKappa)/(D - delta_p*err_P);
//    RoeChi = (D*RoeChi+ s*s*delta_rho*err_P)/(D - delta_p*err_P);
//
//  }

  RoeSoundSpeed2 = RoeChi + RoeKappa*(RoeEnthalpy-0.5*sq_vel);

}

CUpwL2Roe_Flow::CUpwL2Roe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
    
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  kappa = config->GetRoe_Kappa(); // 1 is unstable
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  delta_vel  = new su2double [nDim];
  delta_wave = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
      P_Tensor[iVar] = new su2double [nVar];
      invP_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwL2Roe_Flow::~CUpwL2Roe_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
      delete [] P_Tensor[iVar];
      delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwL2Roe_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
    
  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};
  su2double ProjGridVel = 0.0;
  su2double zeta,Mach_i,Mach_j;
  //unsigned short ssw_L=0, ssw_R=0;
  
  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+4); AD::SetPreaccIn(V_j, nDim+4); AD::SetPreaccIn(Normal, nDim);
  
  /*--- Face area (norm or the normal vector) ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  
  Mach_i = 0.0;
  for (iDim = 0; iDim < nDim; iDim++){
    Velocity_i[iDim] = V_i[iDim+1];
    Mach_i += pow(Velocity_i[iDim],2.0);}
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Pressure_i*Gamma/Density_i));
  Mach_i = sqrt(Mach_i)/SoundSpeed_i;
  
  /*--- Primitive variables at point j ---*/
  
  Mach_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++){
    Velocity_j[iDim] = V_j[iDim+1];
    Mach_j += pow(Velocity_j[iDim],2.0);}
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Pressure_j*Gamma/Density_j));
  Mach_j = sqrt(Mach_j)/SoundSpeed_j;
  
  /*--- L2Roe: a low dissipation version of Roe's approximate Riemann solver for low Mach numbers. IJNMF 2015 ---*/
  
  zeta = min(1.0,max(Mach_i,Mach_j));
  zeta = max(zeta,0.05);
  
  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Roe-averaged variables at interface between i & j ---*/
  
  R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  
  RoeSoundSpeed2 = (Gamma-1)*(RoeEnthalpy-0.5*sq_vel);
  
  /*--- Negative RoeSoundSpeed2, the jump
   variables is too large, exit the subrotuine
   without computing the fluxes ---*/
  
  if (RoeSoundSpeed2 <= 0.0) {
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] = 0.0;
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_i[iVar][iVar] = 0.0;
        val_Jacobian_j[iVar][iVar] = 0.0;
      }
    }
    AD::SetPreaccOut(val_residual, nVar);
    AD::EndPreacc();
    return;
  }
  
  RoeSoundSpeed = sqrt(RoeSoundSpeed2);
  
  /*--- Compute ProjFlux_i ---*/
  
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i);
  
  /*--- Compute ProjFlux_j ---*/
  
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j);
  
  /*--- Compute P and Lambda (do it with the Normal) ---*/
  
  GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
  
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- Projected velocity adjustment due to mesh motion ---*/
  
  if (grid_movement) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    }
    ProjVelocity   -= ProjGridVel;
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }
  
  /*--- Flow eigenvalues and entropy correctors ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;
  
  Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
  Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
  
  /*--- Compute absolute value with Mavriplis' entropy correction ---*/
  
  MaxLambda = fabs(ProjVelocity) + RoeSoundSpeed;
  Delta = config->GetEntropyFix_Coeff();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Lambda[iVar] = max(fabs(Lambda[iVar]), Delta*MaxLambda);
  }

  /*--- Compute wave amplitudes (characteristics) ---*/
  
  proj_delta_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
    proj_delta_vel += delta_vel[iDim]*Normal[iDim];
  }
  delta_p = Pressure_j - Pressure_i;
  delta_rho = Density_j - Density_i;
  proj_delta_vel = (proj_delta_vel/Area)*zeta;
  
  if (nDim == 2) {
    delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
    delta_wave[1] = (UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1])*zeta;
    delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
  } else {
    delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
    delta_wave[1] = (UnitNormal[0]*delta_vel[2]-UnitNormal[2]*delta_vel[0])*zeta;
    delta_wave[2] = (UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1])*zeta;
    delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
  }
  
  if (!implicit) {
    
    /*--- Roe's Flux approximation ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
      for (jVar = 0; jVar < nVar; jVar++)
        val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
    }
    
    /*--- Flux contribution due to grid motion ---*/
    
    if (grid_movement) {
      ProjVelocity = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      for (iVar = 0; iVar < nVar; iVar++) {
        val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      }
    }
  }
  
  else {
      
    /*--- Compute inverse P ---*/
    
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobians of the inviscid flux, scaled by
     kappa because val_resconv ~ kappa*(fc_i+fc_j)*Normal ---*/
    
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, kappa, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, kappa, val_Jacobian_j);
    
    /*--- Diference variables iPoint and jPoint ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      Diff_U[iVar] = U_j[iVar]-U_i[iVar];
    
    /*--- Roe's Flux approximation ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      
      val_residual[iVar] = kappa*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
        
        //val_residual[iVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area;
        val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
        val_Jacobian_i[iVar][jVar] += (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
      }
      
    }

    /*--- Jacobian contributions due to grid motion ---*/
    
    if (grid_movement) {
      ProjVelocity = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      for (iVar = 0; iVar < nVar; iVar++) {
        val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
        
        /*--- Implicit terms ---*/
        
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
        
      }
    }
  
  }  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
    
}

CUpwLMRoe_Flow::CUpwLMRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
    
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  kappa = config->GetRoe_Kappa(); // 1 is unstable
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  delta_vel  = new su2double [nDim];
  delta_wave = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwLMRoe_Flow::~CUpwLMRoe_Flow(void) {
    
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
      delete [] P_Tensor[iVar];
      delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
    
}

void CUpwLMRoe_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
    
  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};
  su2double ProjGridVel = 0.0;
  su2double zeta,Mach_i,Mach_j;
  
  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+4); AD::SetPreaccIn(V_j, nDim+4); AD::SetPreaccIn(Normal, nDim);
  
  /*--- Face area (norm or the normal vector) ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  
  Mach_i = 0.0;
  for (iDim = 0; iDim < nDim; iDim++){
    Velocity_i[iDim] = V_i[iDim+1];
    Mach_i += pow(Velocity_i[iDim],2.0);}
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Pressure_i*Gamma/Density_i));
  Mach_i = sqrt(Mach_i)/SoundSpeed_i;
  
  /*--- Primitive variables at point j ---*/
  
  Mach_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++){
    Velocity_j[iDim] = V_j[iDim+1];
    Mach_j += pow(Velocity_j[iDim],2.0);}
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Pressure_j*Gamma/Density_j));
  Mach_j = sqrt(Mach_j)/SoundSpeed_j;
  
  /*--- Rieper, A low-Mach number fix for Roe's approximate Riemman Solver, JCP 2011 ---*/
  
  zeta = min(1.0,max(Mach_i,Mach_j));
  zeta = max(0.05,zeta);

  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Roe-averaged variables at interface between i & j ---*/
  
  R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  
  RoeSoundSpeed2 = (Gamma-1)*(RoeEnthalpy-0.5*sq_vel);
  
  /*--- Negative RoeSoundSpeed2, the jump
   variables is too large, exit the subrotuine
   without computing the fluxes ---*/
  
  if (RoeSoundSpeed2 <= 0.0) {
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] = 0.0;
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_i[iVar][iVar] = 0.0;
        val_Jacobian_j[iVar][iVar] = 0.0;
      }
    }
    AD::SetPreaccOut(val_residual, nVar);
    AD::EndPreacc();
    return;
  }
  
  RoeSoundSpeed = sqrt(RoeSoundSpeed2);
  
  /*--- Compute ProjFlux_i ---*/
  
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i);
  
  /*--- Compute ProjFlux_j ---*/
  
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j);
  
  /*--- Compute P and Lambda (do it with the Normal) ---*/
  
  GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
  
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- Projected velocity adjustment due to mesh motion ---*/
  
  if (grid_movement) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
        ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    }
    ProjVelocity   -= ProjGridVel;
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }
  
  /*--- Flow eigenvalues and entropy correctors ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;
  
  Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
  Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
  
  /*--- Compute absolute value with Mavriplis' entropy correction ---*/
  
  MaxLambda = fabs(ProjVelocity) + RoeSoundSpeed;
  Delta = config->GetEntropyFix_Coeff();
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Lambda[iVar] = max(fabs(Lambda[iVar]), Delta*MaxLambda);
  }
  
  /*--- Compute wave amplitudes (characteristics) ---*/
  
  proj_delta_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
    proj_delta_vel += delta_vel[iDim]*Normal[iDim];
  }
  delta_p = Pressure_j - Pressure_i;
  delta_rho = Density_j - Density_i;
  proj_delta_vel = (proj_delta_vel/Area)*zeta;
  
  if (nDim == 2) {
    delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
    delta_wave[1] = (UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1]);
    delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
  } else {
    delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
    delta_wave[1] = (UnitNormal[0]*delta_vel[2]-UnitNormal[2]*delta_vel[0]);
    delta_wave[2] = (UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1]);
    delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
  }
  
  if (!implicit) {
    
    /*--- Roe's Flux approximation ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
      for (jVar = 0; jVar < nVar; jVar++)
        val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
    }
    
    /*--- Flux contribution due to grid motion ---*/
    
    if (grid_movement) {
      ProjVelocity = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      for (iVar = 0; iVar < nVar; iVar++) {
        val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      }
    }
  }
  
  else {
    
    /*--- Compute inverse P ---*/
    
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobians of the inviscid flux, scaled by
     kappa because val_resconv ~ kappa*(fc_i+fc_j)*Normal ---*/
    
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, kappa, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, kappa, val_Jacobian_j);
    
    /*--- Diference variables iPoint and jPoint ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      Diff_U[iVar] = U_j[iVar]-U_i[iVar];
    
    /*--- Roe's Flux approximation ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
        
      val_residual[iVar] = kappa*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        
        for (kVar = 0; kVar < nVar; kVar++)
            Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
        
        //val_residual[iVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area;
        val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
        val_Jacobian_i[iVar][jVar] += (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
        
      }
      
    }
    
    /*--- Jacobian contributions due to grid motion ---*/
    
    if (grid_movement) {
      ProjVelocity = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      for (iVar = 0; iVar < nVar; iVar++) {
        val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
        
        /*--- Implicit terms ---*/
        
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
    
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}


CUpwMSW_Flow::CUpwMSW_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  /*--- Set booleans from CConfig settings ---*/
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Allocate arrays ---*/
  Diff_U   = new su2double [nVar];
  Fc_i     = new su2double [nVar];
  Fc_j     = new su2double [nVar];
  Lambda_i = new su2double [nVar];
  Lambda_j = new su2double [nVar];
  
  u_i       = new su2double [nDim];
  u_j       = new su2double [nDim];
  ust_i    = new su2double [nDim];
  ust_j    = new su2double [nDim];
  Vst_i    = new su2double [nPrimVar];
  Vst_j    = new su2double [nPrimVar];
  Ust_i    = new su2double [nVar];
  Ust_j    = new su2double [nVar];
  
  Velst_i    = new su2double [nDim];
  Velst_j    = new su2double [nDim];
  
  P_Tensor    = new su2double* [nVar];
  invP_Tensor  = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar]    = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }
  
}

CUpwMSW_Flow::~CUpwMSW_Flow(void) {
  
  delete [] Diff_U;
  delete [] Fc_i;
  delete [] Fc_j;
  delete [] Lambda_i;
  delete [] Lambda_j;
  
  delete [] u_i;
  delete [] u_j;
  delete [] ust_i;
  delete [] ust_j;
  delete [] Ust_i;
  delete [] Vst_i;
  delete [] Ust_j;
  delete [] Vst_j;
  delete [] Velst_i;
  delete [] Velst_j;
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;

}

void CUpwMSW_Flow::ComputeResidual(su2double *val_residual,
                                   su2double **val_Jacobian_i,
                                   su2double **val_Jacobian_j, CConfig *config) {
  
  unsigned short iDim, iVar, jVar, kVar;
  su2double P_i, P_j;
  su2double ProjVel_i, ProjVel_j, ProjVelst_i, ProjVelst_j;
  su2double sqvel_i, sqvel_j;
  su2double alpha, w, dp, onemw;
  su2double Proj_ModJac_Tensor_i, Proj_ModJac_Tensor_j;
  
  /*--- Set parameters in the numerical method ---*/
  alpha = 6.0;
  
  /*--- Calculate supporting geometry parameters ---*/
  
  Area = 0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Initialize flux & Jacobian vectors ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Fc_i[iVar] = 0.0;
    Fc_j[iVar] = 0.0;
  }
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_i[iVar][jVar] = 0.0;
        val_Jacobian_j[iVar][jVar] = 0.0;
      }
    }
  }
  
  /*--- Load variables from nodes i & j ---*/
  
  rhos_i = V_i[0];
  rhos_j = V_j[0];
  for (iDim = 0; iDim < nDim; iDim++) {
    u_i[iDim] = V_i[iDim+1];
    u_j[iDim] = V_j[iDim+1];
  }
  P_i = V_i[nDim+1];
  P_j = V_j[nDim+1];
  
  /*--- Calculate supporting quantities ---*/
  
  sqvel_i   = 0.0; sqvel_j   = 0.0;
  ProjVel_i = 0.0; ProjVel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    sqvel_i   += u_i[iDim]*u_i[iDim];
    sqvel_j   += u_j[iDim]*u_j[iDim];
    ProjVel_i += u_i[iDim]*UnitNormal[iDim];
    ProjVel_j += u_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- Calculate the state weighting function ---*/
  
  dp = fabs(P_j-P_i) / min(P_j, P_i);
  w = 0.5 * (1.0/(pow(alpha*dp,2.0) +1.0));
  onemw = 1.0 - w;
  
  /*--- Calculate weighted state vector (*) for i & j ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Ust_i[iVar] = onemw*U_i[iVar] + w*U_j[iVar];
    Ust_j[iVar] = onemw*U_j[iVar] + w*U_i[iVar];
  }
  for (iVar = 0; iVar < nDim+5; iVar++) {
    Vst_i[iVar] = onemw*V_i[iVar] + w*V_j[iVar];
    Vst_j[iVar] = onemw*V_j[iVar] + w*V_i[iVar];
  }
  ProjVelst_i = onemw*ProjVel_i + w*ProjVel_j;
  ProjVelst_j = onemw*ProjVel_j + w*ProjVel_i;
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Velst_i[iDim] = Vst_i[iDim+1];
    Velst_j[iDim] = Vst_j[iDim+1];
  }
  
  /*--- Flow eigenvalues at i (Lambda+) ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
  Lambda_i[iDim]      = 0.5*(ProjVelst_i + fabs(ProjVelst_i));
  }

  Lambda_i[nDim] = 0.5*( ProjVelst_i + Vst_i[nDim+4] + fabs(ProjVelst_i + Vst_i[nDim+4])  );
  Lambda_i[nDim+1]   = 0.5*( ProjVelst_i - Vst_i[nDim+4] + fabs(ProjVelst_i - Vst_i[nDim+4])  );
  
  /*--- Compute projected P, invP, and Lambda ---*/
  
  GetPMatrix(&Vst_i[nDim+2], Velst_i, &Vst_i[nDim+4], UnitNormal, P_Tensor);
  GetPMatrix_inv(&Vst_i[nDim+2], Velst_i, &Vst_i[nDim+4], UnitNormal, invP_Tensor);
  
  /*--- Projected flux (f+) at i ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_i = 0.0;
      
      /*--- Compute Proj_ModJac_Tensor = P x Lambda+ x inverse P ---*/
      
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_i += P_Tensor[iVar][kVar]*Lambda_i[kVar]*invP_Tensor[kVar][jVar];
      Fc_i[iVar] += Proj_ModJac_Tensor_i*U_i[jVar]*Area;
      if (implicit)
        val_Jacobian_i[iVar][jVar] += Proj_ModJac_Tensor_i*Area;
    }
  }
  
  /*--- Flow eigenvalues at j (Lambda-) ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Lambda_j[iDim]          = 0.5*(ProjVelst_j - fabs(ProjVelst_j));
  }
  Lambda_j[nDim] = 0.5*(     ProjVelst_j + Vst_j[nDim+4] -
                                   fabs(ProjVelst_j + Vst_j[nDim+4])  );
  Lambda_j[nDim+1]   = 0.5*(     ProjVelst_j - Vst_j[nDim+4] -
                                   fabs(ProjVelst_j - Vst_j[nDim+4])  );
  
  /*--- Compute projected P, invP, and Lambda ---*/
  
  GetPMatrix(&Vst_j[nDim+2], Velst_j, &Vst_j[nDim+4], UnitNormal, P_Tensor);
  GetPMatrix_inv(&Vst_j[nDim+2], Velst_j, &Vst_j[nDim+4], UnitNormal, invP_Tensor);
  
  /*--- Projected flux (f-) ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_j = 0.0;
      /*--- Compute Proj_ModJac_Tensor = P x Lambda- x inverse P ---*/
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_j += P_Tensor[iVar][kVar]*Lambda_j[kVar]*invP_Tensor[kVar][jVar];
      Fc_j[iVar] += Proj_ModJac_Tensor_j*U_j[jVar]*Area;
      if (implicit)
        val_Jacobian_j[iVar][jVar] += Proj_ModJac_Tensor_j*Area;
    }
  }
  
  /*--- Flux splitting ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = Fc_i[iVar]+Fc_j[iVar];
  }
  
}

CUpwTurkel_Flow::CUpwTurkel_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Beta_min = config->GetminTurkelBeta();
  Beta_max = config->GetmaxTurkelBeta();
  
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  absPeJac = new su2double* [nVar];
  invRinvPe = new su2double* [nVar];
  R_Tensor  = new su2double* [nVar];
  Matrix    = new su2double* [nVar];
  Art_Visc  = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    absPeJac[iVar] = new su2double [nVar];
    invRinvPe[iVar] = new su2double [nVar];
    Matrix[iVar] = new su2double [nVar];
    Art_Visc[iVar] = new su2double [nVar];
    R_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwTurkel_Flow::~CUpwTurkel_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] absPeJac[iVar];
    delete [] invRinvPe[iVar];
    delete [] Matrix[iVar];
    delete [] Art_Visc[iVar];
    delete [] R_Tensor[iVar];
  }
  delete [] Matrix;
  delete [] Art_Visc;
  delete [] absPeJac;
  delete [] invRinvPe;
  delete [] R_Tensor;
  
}

void CUpwTurkel_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  /*--- Face area (norm or the normal vector) ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_i[iDim] = V_i[iDim+1];
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Pressure_i*Gamma/Density_i));

  /*--- Primitive variables at point j ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_j[iDim] = V_j[iDim+1];
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Pressure_j*Gamma/Density_j));

  /*--- Recompute conservative variables ---*/
  
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Roe-averaged variables at interface between i & j ---*/
  
  R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
  RoePressure = RoeDensity/Gamma*RoeSoundSpeed*RoeSoundSpeed;
  
  /*--- Compute ProjFlux_i ---*/
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i);
  
  /*--- Compute ProjFlux_j ---*/
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j);
  
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- Projected velocity adjustment due to mesh motion ---*/
  if (grid_movement) {
    su2double ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    }
    ProjVelocity   -= ProjGridVel;
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }
  
  /*--- First few flow eigenvalues of A.Normal with the normal---*/
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;
  
  local_Mach = sqrt(sq_vel)/RoeSoundSpeed;
  Beta      = max(Beta_min, min(local_Mach, Beta_max));
  Beta2      = Beta*Beta;
  
  one_m_Betasqr        = 1.0 - Beta2;  // 1-Beta*Beta
  one_p_Betasqr        = 1.0 + Beta2;  // 1+Beta*Beta
  sqr_one_m_Betasqr_Lam1 = pow((one_m_Betasqr*Lambda[0]),2); // [(1-Beta^2)*Lambda[0]]^2
  sqr_two_Beta_c_Area    = pow(2.0*Beta*RoeSoundSpeed*Area,2); // [2*Beta*c*Area]^2
  
  /*--- The rest of the flow eigenvalues of preconditioned matrix---*/
  Lambda[nVar-2] = 0.5 * ( one_p_Betasqr*Lambda[0] + sqrt( sqr_one_m_Betasqr_Lam1 + sqr_two_Beta_c_Area));
  Lambda[nVar-1] = 0.5 * ( one_p_Betasqr*Lambda[0] - sqrt( sqr_one_m_Betasqr_Lam1 + sqr_two_Beta_c_Area));
  
  s_hat = 1.0/Area * (Lambda[nVar-1] - Lambda[0]*Beta2);
  r_hat = 1.0/Area * (Lambda[nVar-2] - Lambda[0]*Beta2);
  t_hat = 0.5/Area * (Lambda[nVar-1] - Lambda[nVar-2]);
  rhoB2a2 = RoeDensity*Beta2*RoeSoundSpeed*RoeSoundSpeed;
  
  /*--- Diference variables iPoint and jPoint and absolute value of the eigen values---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_U[iVar] = U_j[iVar]-U_i[iVar];
    Lambda[iVar] = fabs(Lambda[iVar]);
  }
  
  /*--- Compute the absolute Preconditioned Jacobian in entropic Variables (do it with the Unitary Normal) ---*/
  GetPrecondJacobian(Beta2, r_hat, s_hat, t_hat, rhoB2a2, Lambda, UnitNormal, absPeJac);
  
  /*--- Compute the matrix from entropic variables to conserved variables ---*/
  GetinvRinvPe(Beta2, RoeEnthalpy, RoeSoundSpeed, RoeDensity, RoeVelocity, invRinvPe);
  
  /*--- Compute the matrix from entropic variables to conserved variables ---*/
  GetRMatrix(RoePressure, RoeSoundSpeed, RoeDensity, RoeVelocity, R_Tensor);
  
  if (implicit) {
    /*--- Jacobians of the inviscid flux, scaled by
     0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
  }
  
  for (iVar = 0; iVar < nVar; iVar ++) {
    for (jVar = 0; jVar < nVar; jVar ++) {
      Matrix[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Matrix[iVar][jVar]  += absPeJac[iVar][kVar]*R_Tensor[kVar][jVar];
    }
  }
  
  for (iVar = 0; iVar < nVar; iVar ++) {
    for (jVar = 0; jVar < nVar; jVar ++) {
      Art_Visc[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Art_Visc[iVar][jVar]  += invRinvPe[iVar][kVar]*Matrix[kVar][jVar];
    }
  }
  
  /*--- Roe's Flux approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual[iVar] -= 0.5*Art_Visc[iVar][jVar]*Diff_U[jVar];
      if (implicit) {
        val_Jacobian_i[iVar][jVar] += 0.5*Art_Visc[iVar][jVar];
        val_Jacobian_j[iVar][jVar] -= 0.5*Art_Visc[iVar][jVar];
      }
    }
  }
  
  /*--- Contributions due to mesh motion---*/
  if (grid_movement) {
    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      /*--- Implicit terms ---*/
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }
  
}

CAvgGrad_Base::CAvgGrad_Base(unsigned short val_nDim,
                             unsigned short val_nVar,
                             unsigned short val_nPrimVar,
                             bool val_correct_grad,
                             CConfig *config)
    : CNumerics(val_nDim, val_nVar, config),
      nPrimVar(val_nPrimVar),
      correct_gradient(val_correct_grad) {

  unsigned short iVar, iDim;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  TauWall_i = 0; TauWall_j = 0;

  PrimVar_i = new su2double [nPrimVar];
  PrimVar_j = new su2double [nPrimVar];
  Mean_PrimVar = new su2double [nPrimVar];

  Mean_GradPrimVar = new su2double* [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++)
    Mean_GradPrimVar[iVar] = new su2double [nDim];

  Edge_Vector = new su2double[nDim];

  if (correct_gradient) {
    Proj_Mean_GradPrimVar_Edge = new su2double[val_nPrimVar];
  } else {
    Proj_Mean_GradPrimVar_Edge = NULL;
  }

  tau_jacobian_i = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    tau_jacobian_i[iDim] = new su2double [nVar];
  }
  heat_flux_vector = new su2double[nDim];
  heat_flux_jac_i = new su2double[nVar];

}

CAvgGrad_Base::~CAvgGrad_Base() {

  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Mean_PrimVar;
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;

  if (tau_jacobian_i != NULL) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      delete [] tau_jacobian_i[iDim];
    }
    delete [] tau_jacobian_i;
  }
  if (heat_flux_vector != NULL) {
    delete [] heat_flux_vector;
  }
  if (heat_flux_jac_i != NULL) {
    delete [] heat_flux_jac_i;
  }
  
  delete [] Edge_Vector;
  if (Proj_Mean_GradPrimVar_Edge != NULL)
    delete [] Proj_Mean_GradPrimVar_Edge;
}

void CAvgGrad_Base::CorrectGradient(su2double** GradPrimVar,
                                    const su2double* val_PrimVar_i,
                                    const su2double* val_PrimVar_j,
                                    const su2double* val_edge_vector,
                                    const su2double val_dist_ij_2,
                                    const unsigned short val_nPrimVar) {
  for (unsigned short iVar = 0; iVar < val_nPrimVar; iVar++) {
    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Proj_Mean_GradPrimVar_Edge[iVar] += GradPrimVar[iVar][iDim]*val_edge_vector[iDim];
    }
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                 (val_PrimVar_j[iVar]-val_PrimVar_i[iVar]))*val_edge_vector[iDim] / val_dist_ij_2;
    }
  }
}

void CAvgGrad_Base::SetStressTensor(const su2double *val_primvar,
                           const su2double* const *val_gradprimvar,
                           const su2double val_turb_ke,
                           const su2double val_laminar_viscosity,
                           const su2double val_eddy_viscosity) {

  unsigned short iDim, jDim;
  const su2double Density = val_primvar[nDim+2];
  const su2double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;

  su2double div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    div_vel += val_gradprimvar[iDim+1][iDim];

  /* --- If UQ methodology is used, calculate tau using the perturbed reynolds stress tensor --- */

  if (using_uq){
    for (iDim = 0 ; iDim < nDim; iDim++)
      for (jDim = 0 ; jDim < nDim; jDim++)
        tau[iDim][jDim] = val_laminar_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
        - TWO3*val_laminar_viscosity*div_vel*delta[iDim][jDim] - Density * MeanPerturbedRSM[iDim][jDim];

  } else {

    for (iDim = 0 ; iDim < nDim; iDim++)
      for (jDim = 0 ; jDim < nDim; jDim++)
        tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
                          - TWO3*total_viscosity*div_vel*delta[iDim][jDim]
                          - TWO3*Density*val_turb_ke*delta[iDim][jDim];
  }
}

void CAvgGrad_Base::AddQCR(const su2double* const *val_gradprimvar) {

  su2double den_aux, c_cr1= 0.3, O_ik, O_jk;
  unsigned short iDim, jDim, kDim;

  /*--- Denominator Antisymmetric normalized rotation tensor ---*/

  den_aux = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    for (jDim = 0 ; jDim < nDim; jDim++)
      den_aux += val_gradprimvar[iDim+1][jDim] * val_gradprimvar[iDim+1][jDim];
  den_aux = sqrt(max(den_aux,1E-10));

  /*--- Adding the QCR contribution ---*/

  for (iDim = 0 ; iDim < nDim; iDim++){
    for (jDim = 0 ; jDim < nDim; jDim++){
      for (kDim = 0 ; kDim < nDim; kDim++){
        O_ik = (val_gradprimvar[iDim+1][kDim] - val_gradprimvar[kDim+1][iDim])/ den_aux;
        O_jk = (val_gradprimvar[jDim+1][kDim] - val_gradprimvar[kDim+1][jDim])/ den_aux;
        tau[iDim][jDim] -= c_cr1 * ((O_ik * tau[jDim][kDim]) + (O_jk * tau[iDim][kDim]));
      }
    }
  }
}

void CAvgGrad_Base::AddTauWall(const su2double *val_normal,
                               const su2double val_tau_wall) {

  unsigned short iDim, jDim;
  su2double TauNormal, TauElem[3], TauTangent[3], WallShearStress, Area, UnitNormal[3];

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += val_normal[iDim]*val_normal[iDim];
  Area = sqrt(Area);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = val_normal[iDim]/Area;

  /*--- First, compute wall shear stress as the magnitude of the wall-tangential
   component of the shear stress tensor---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    TauElem[iDim] = 0.0;
    for (jDim = 0; jDim < nDim; jDim++)
      TauElem[iDim] += tau[iDim][jDim]*UnitNormal[jDim];
  }

  TauNormal = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    TauNormal += TauElem[iDim] * UnitNormal[iDim];

  for (iDim = 0; iDim < nDim; iDim++)
    TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];

  WallShearStress = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    WallShearStress += TauTangent[iDim]*TauTangent[iDim];
  WallShearStress = sqrt(WallShearStress);

  /*--- Scale the stress tensor by the ratio of the wall shear stress
   to the computed representation of the shear stress ---*/

  for (iDim = 0 ; iDim < nDim; iDim++)
    for (jDim = 0 ; jDim < nDim; jDim++)
      tau[iDim][jDim] = tau[iDim][jDim]*(val_tau_wall/WallShearStress);
}

void CAvgGrad_Base::GetMeanRateOfStrainMatrix(su2double **S_ij) const
{
  /* --- Calculate the rate of strain tensor, using mean velocity gradients --- */

  if (nDim == 3){
    S_ij[0][0] = Mean_GradPrimVar[1][0];
    S_ij[1][1] = Mean_GradPrimVar[2][1];
    S_ij[2][2] = Mean_GradPrimVar[3][2];
    S_ij[0][1] = 0.5 * (Mean_GradPrimVar[1][1] + Mean_GradPrimVar[2][0]);
    S_ij[0][2] = 0.5 * (Mean_GradPrimVar[1][2] + Mean_GradPrimVar[3][0]);
    S_ij[1][2] = 0.5 * (Mean_GradPrimVar[2][2] + Mean_GradPrimVar[3][1]);
    S_ij[1][0] = S_ij[0][1];
    S_ij[2][1] = S_ij[1][2];
    S_ij[2][0] = S_ij[0][2];
  }
  else {
    S_ij[0][0] = Mean_GradPrimVar[1][0];
    S_ij[1][1] = Mean_GradPrimVar[2][1];
    S_ij[2][2] = 0.0;
    S_ij[0][1] = 0.5 * (Mean_GradPrimVar[1][1] + Mean_GradPrimVar[2][0]);
    S_ij[0][2] = 0.0;
    S_ij[1][2] = 0.0;
    S_ij[1][0] = S_ij[0][1];
    S_ij[2][1] = S_ij[1][2];
    S_ij[2][0] = S_ij[0][2];

  }
}

void CAvgGrad_Base::SetReynoldsStressMatrix(su2double turb_ke){
  unsigned short iDim, jDim;
  su2double **S_ij = new su2double* [3];
  su2double muT = Mean_Eddy_Viscosity;
  su2double divVel = 0;
  su2double density;
  su2double TWO3 = 2.0/3.0;
  density = Mean_PrimVar[nDim+2];

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
      - muT / density * (2 * S_ij[iDim][jDim] - TWO3 * divVel * delta3[iDim][jDim]);
    }
  }

  for (iDim = 0; iDim < 3; iDim++)
    delete [] S_ij[iDim];
  delete [] S_ij;
}

void CAvgGrad_Base::SetPerturbedRSM(su2double turb_ke, CConfig *config){

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
  else if (Eig_Val_Comp== 2) {
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

  /* rebuild c1c,c2c,c3c based on perturbed barycentric coordinates */
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

  /* compute perturbed Reynolds stress matrix; use under-relaxation factor (uq_urlx)*/
  for (iDim = 0; iDim< 3; iDim++){
    for (jDim = 0; jDim < 3; jDim++){
      MeanPerturbedRSM[iDim][jDim] = 2.0 * turb_ke * (newA_ij[iDim][jDim] + 1.0/3.0 * delta3[iDim][jDim]);
      MeanPerturbedRSM[iDim][jDim] = MeanReynoldsStress[iDim][jDim] +
      uq_urlx*(MeanPerturbedRSM[iDim][jDim] - MeanReynoldsStress[iDim][jDim]);
    }
  }

}


void CAvgGrad_Base::SetTauJacobian(const su2double *val_Mean_PrimVar,
                                   const su2double val_laminar_viscosity,
                                   const su2double val_eddy_viscosity,
                                   const su2double val_dist_ij,
                                   const su2double *val_normal) {

  /*--- QCR and wall functions are **not** accounted for here ---*/

  const su2double Density = val_Mean_PrimVar[nDim+2];
  const su2double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
  const su2double xi = total_viscosity/(Density*val_dist_ij);

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      // Jacobian w.r.t. momentum
      tau_jacobian_i[iDim][jDim+1] = -xi*(delta[iDim][jDim] + val_normal[iDim]*val_normal[jDim]/3.0);
    }
    // Jacobian w.r.t. density
    tau_jacobian_i[iDim][0] = 0;
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
       tau_jacobian_i[iDim][0] -= tau_jacobian_i[iDim][jDim+1]*val_Mean_PrimVar[jDim+1];
    }
    // Jacobian w.r.t. energy
    tau_jacobian_i[iDim][nDim+1] = 0;
  }
}

void CAvgGrad_Base::SetIncTauJacobian(const su2double val_laminar_viscosity,
                                      const su2double val_eddy_viscosity,
                                      const su2double val_dist_ij,
                                      const su2double *val_normal) {

  const su2double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
  const su2double xi = total_viscosity/val_dist_ij;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    tau_jacobian_i[iDim][0] = 0;
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      tau_jacobian_i[iDim][jDim+1] = -xi*(delta[iDim][jDim] + val_normal[iDim]*val_normal[jDim]/3.0);
    }
    tau_jacobian_i[iDim][nDim+1] = 0;
  }
}

void CAvgGrad_Base::GetViscousProjFlux(const su2double *val_primvar,
                                       const su2double *val_normal) {

  /*--- Primitive variables -> [Temp vel_x vel_y vel_z Pressure] ---*/

  if (nDim == 2) {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2]+
        heat_flux_vector[0];
    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2]+
        heat_flux_vector[1];
  } else {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][2];
    Flux_Tensor[4][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2] + tau[0][2]*val_primvar[3] +
        heat_flux_vector[0];
    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][2];
    Flux_Tensor[4][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2] + tau[1][2]*val_primvar[3] +
        heat_flux_vector[1];
    Flux_Tensor[0][2] = 0.0;
    Flux_Tensor[1][2] = tau[2][0];
    Flux_Tensor[2][2] = tau[2][1];
    Flux_Tensor[3][2] = tau[2][2];
    Flux_Tensor[4][2] = tau[2][0]*val_primvar[1] + tau[2][1]*val_primvar[2] + tau[2][2]*val_primvar[3] +
        heat_flux_vector[2];
  }
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Proj_Flux_Tensor[iVar] = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim] * val_normal[iDim];
  }
  
}

void CAvgGrad_Base::GetViscousProjJacs(const su2double *val_Mean_PrimVar,
                                       const su2double val_dS,
                                       const su2double *val_Proj_Visc_Flux,
                                       su2double **val_Proj_Jac_Tensor_i,
                                       su2double **val_Proj_Jac_Tensor_j) {

  const su2double Density = val_Mean_PrimVar[nDim+2];
  const su2double factor = 0.5/Density;

  if (nDim == 2) {

    val_Proj_Jac_Tensor_i[0][0] = 0.0;
    val_Proj_Jac_Tensor_i[0][1] = 0.0;
    val_Proj_Jac_Tensor_i[0][2] = 0.0;
    val_Proj_Jac_Tensor_i[0][3] = 0.0;
    val_Proj_Jac_Tensor_i[1][0] = val_dS*tau_jacobian_i[0][0];
    val_Proj_Jac_Tensor_i[1][1] = val_dS*tau_jacobian_i[0][1];
    val_Proj_Jac_Tensor_i[1][2] = val_dS*tau_jacobian_i[0][2];
    val_Proj_Jac_Tensor_i[1][3] = val_dS*tau_jacobian_i[0][3];
    val_Proj_Jac_Tensor_i[2][0] = val_dS*tau_jacobian_i[1][0];
    val_Proj_Jac_Tensor_i[2][1] = val_dS*tau_jacobian_i[1][1];
    val_Proj_Jac_Tensor_i[2][2] = val_dS*tau_jacobian_i[1][2];
    val_Proj_Jac_Tensor_i[2][3] = val_dS*tau_jacobian_i[1][3];
    const su2double contraction = tau_jacobian_i[0][0]*val_Mean_PrimVar[1] +
                                  tau_jacobian_i[1][0]*val_Mean_PrimVar[2];
    val_Proj_Jac_Tensor_i[3][0] = val_dS*(contraction - heat_flux_jac_i[0]);
    val_Proj_Jac_Tensor_i[3][1] = -val_dS*(tau_jacobian_i[0][0] + heat_flux_jac_i[1]); 
    val_Proj_Jac_Tensor_i[3][2] = -val_dS*(tau_jacobian_i[1][0] + heat_flux_jac_i[2]); 
    val_Proj_Jac_Tensor_i[3][3] = -val_dS*heat_flux_jac_i[3];
    
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      for (unsigned short jVar = 0; jVar < nVar; jVar++)
        val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];

    const su2double proj_viscousflux_vel= val_Proj_Visc_Flux[1]*val_Mean_PrimVar[1] + 
                                          val_Proj_Visc_Flux[2]*val_Mean_PrimVar[2];
    val_Proj_Jac_Tensor_i[3][0] -= factor*proj_viscousflux_vel;
    val_Proj_Jac_Tensor_j[3][0] -= factor*proj_viscousflux_vel;
    val_Proj_Jac_Tensor_i[3][1] += factor*val_Proj_Visc_Flux[1];
    val_Proj_Jac_Tensor_j[3][1] += factor*val_Proj_Visc_Flux[1];
    val_Proj_Jac_Tensor_i[3][2] += factor*val_Proj_Visc_Flux[2];
    val_Proj_Jac_Tensor_j[3][2] += factor*val_Proj_Visc_Flux[2];
    
    
  } else {

    val_Proj_Jac_Tensor_i[0][0] = 0.0;
    val_Proj_Jac_Tensor_i[0][1] = 0.0;
    val_Proj_Jac_Tensor_i[0][2] = 0.0;
    val_Proj_Jac_Tensor_i[0][3] = 0.0;
    val_Proj_Jac_Tensor_i[0][4] = 0.0;
    val_Proj_Jac_Tensor_i[1][0] = val_dS*tau_jacobian_i[0][0];
    val_Proj_Jac_Tensor_i[1][1] = val_dS*tau_jacobian_i[0][1];
    val_Proj_Jac_Tensor_i[1][2] = val_dS*tau_jacobian_i[0][2];
    val_Proj_Jac_Tensor_i[1][3] = val_dS*tau_jacobian_i[0][3];
    val_Proj_Jac_Tensor_i[1][4] = val_dS*tau_jacobian_i[0][4];
    val_Proj_Jac_Tensor_i[2][0] = val_dS*tau_jacobian_i[1][0];
    val_Proj_Jac_Tensor_i[2][1] = val_dS*tau_jacobian_i[1][1];
    val_Proj_Jac_Tensor_i[2][2] = val_dS*tau_jacobian_i[1][2];
    val_Proj_Jac_Tensor_i[2][3] = val_dS*tau_jacobian_i[1][3];
    val_Proj_Jac_Tensor_i[2][4] = val_dS*tau_jacobian_i[1][4];
    val_Proj_Jac_Tensor_i[3][0] = val_dS*tau_jacobian_i[2][0];
    val_Proj_Jac_Tensor_i[3][1] = val_dS*tau_jacobian_i[2][1];
    val_Proj_Jac_Tensor_i[3][2] = val_dS*tau_jacobian_i[2][2];
    val_Proj_Jac_Tensor_i[3][3] = val_dS*tau_jacobian_i[2][3];
    val_Proj_Jac_Tensor_i[3][4] = val_dS*tau_jacobian_i[2][4];
    const su2double contraction = tau_jacobian_i[0][0]*val_Mean_PrimVar[1] +
                                  tau_jacobian_i[1][0]*val_Mean_PrimVar[2] +
                                  tau_jacobian_i[2][0]*val_Mean_PrimVar[3];
    val_Proj_Jac_Tensor_i[4][0] = val_dS*(contraction - heat_flux_jac_i[0]);
    val_Proj_Jac_Tensor_i[4][1] = -val_dS*(tau_jacobian_i[0][0] + heat_flux_jac_i[1]); 
    val_Proj_Jac_Tensor_i[4][2] = -val_dS*(tau_jacobian_i[1][0] + heat_flux_jac_i[2]); 
    val_Proj_Jac_Tensor_i[4][3] = -val_dS*(tau_jacobian_i[2][0] + heat_flux_jac_i[3]); 
    val_Proj_Jac_Tensor_i[4][4] = -val_dS*heat_flux_jac_i[4];

    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      for (unsigned short jVar = 0; jVar < nVar; jVar++)
        val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];
    
    const su2double proj_viscousflux_vel= val_Proj_Visc_Flux[1]*val_Mean_PrimVar[1] + 
                                          val_Proj_Visc_Flux[2]*val_Mean_PrimVar[2] +
                                          val_Proj_Visc_Flux[3]*val_Mean_PrimVar[3];
    val_Proj_Jac_Tensor_i[4][0] -= factor*proj_viscousflux_vel;
    val_Proj_Jac_Tensor_j[4][0] -= factor*proj_viscousflux_vel;
    val_Proj_Jac_Tensor_i[4][1] += factor*val_Proj_Visc_Flux[1];
    val_Proj_Jac_Tensor_j[4][1] += factor*val_Proj_Visc_Flux[1];
    val_Proj_Jac_Tensor_i[4][2] += factor*val_Proj_Visc_Flux[2];
    val_Proj_Jac_Tensor_j[4][2] += factor*val_Proj_Visc_Flux[2];
    val_Proj_Jac_Tensor_i[4][3] += factor*val_Proj_Visc_Flux[3];
    val_Proj_Jac_Tensor_j[4][3] += factor*val_Proj_Visc_Flux[3];

  }

}

CAvgGrad_Flow::CAvgGrad_Flow(unsigned short val_nDim,
                             unsigned short val_nVar,
                             bool val_correct_grad,
                             CConfig *config)
    : CAvgGrad_Base(val_nDim, val_nVar, val_nDim+3, val_correct_grad, config) {

}

CAvgGrad_Flow::~CAvgGrad_Flow(void) {

}

void CAvgGrad_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+9);   AD::SetPreaccIn(V_j, nDim+9);
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);
  AD::SetPreaccIn(PrimVar_Grad_j, nDim+1, nDim);
  AD::SetPreaccIn(turb_ke_i); AD::SetPreaccIn(turb_ke_j);
  AD::SetPreaccIn(Normal, nDim);

  unsigned short iVar, jVar, iDim;

  /*--- Normalized normal vector ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  for (iVar = 0; iVar < nPrimVar; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }
  

  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }

  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+5]; Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6]; Eddy_Viscosity_j = V_j[nDim+6];

  /*--- Mean Viscosities and turbulent kinetic energy---*/
  
  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke = 0.5*(turb_ke_i + turb_ke_j);
  
  /*--- Mean gradient approximation ---*/

  for (iVar = 0; iVar < nDim+1; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
    }
  }

  /*--- Projection of the mean gradient in the direction of the edge ---*/

  if (correct_gradient && dist_ij_2 != 0.0) {
    CorrectGradient(Mean_GradPrimVar, PrimVar_i, PrimVar_j, Edge_Vector,
                    dist_ij_2, nDim+1);
  }
  
  /*--- Wall shear stress values (wall functions) ---*/
  
  if (TauWall_i > 0.0 && TauWall_j > 0.0) Mean_TauWall = 0.5*(TauWall_i + TauWall_j);
  else if (TauWall_i > 0.0) Mean_TauWall = TauWall_i;
  else if (TauWall_j > 0.0) Mean_TauWall = TauWall_j;
  else Mean_TauWall = -1.0;

  /* --- If using UQ methodology, set Reynolds Stress tensor and perform perturbation--- */

  if (using_uq){
    SetReynoldsStressMatrix(Mean_turb_ke);
    SetPerturbedRSM(Mean_turb_ke, config);
  }

  /*--- Get projected flux tensor ---*/

  SetStressTensor(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke,
         Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  if (config->GetQCR()) AddQCR(Mean_GradPrimVar);
  if (Mean_TauWall > 0) AddTauWall(Normal, Mean_TauWall);

  SetHeatFluxVector(Mean_GradPrimVar, Mean_Laminar_Viscosity,
                    Mean_Eddy_Viscosity);

  GetViscousProjFlux(Mean_PrimVar, Normal);

  /*--- Update viscous residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
  /*--- Compute the implicit part ---*/
  
  if (implicit) {
    
    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    } else {
      const su2double dist_ij = sqrt(dist_ij_2);
      SetTauJacobian(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
                     dist_ij, UnitNormal);
      SetHeatFluxJacobian(Mean_PrimVar, Mean_Laminar_Viscosity,
                          Mean_Eddy_Viscosity, dist_ij, UnitNormal);
      GetViscousProjJacs(Mean_PrimVar, Area, Proj_Flux_Tensor,
                         val_Jacobian_i, val_Jacobian_j);
    }
    
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}



void CAvgGrad_Flow::SetHeatFluxVector(const su2double* const *val_gradprimvar,
                                      const su2double val_laminar_viscosity,
                                      const su2double val_eddy_viscosity) {

  const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  const su2double heat_flux_factor = Cp * (val_laminar_viscosity/Prandtl_Lam + val_eddy_viscosity/Prandtl_Turb);

  /*--- Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure] ---*/

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    heat_flux_vector[iDim] = heat_flux_factor*val_gradprimvar[0][iDim];
  }
}

void CAvgGrad_Flow::SetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                                        const su2double val_laminar_viscosity,
                                        const su2double val_eddy_viscosity,
                                        const su2double val_dist_ij,
                                        const su2double *val_normal) {
  su2double sqvel = 0.0;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
  }

  const su2double Density = val_Mean_PrimVar[nDim+2];
  const su2double Pressure = val_Mean_PrimVar[nDim+1];
  const su2double phi = Gamma_Minus_One/Density;

  /*--- R times partial derivatives of temp. ---*/

  const su2double R_dTdu0 = -Pressure/(Density*Density) + 0.5*sqvel*phi;
  const su2double R_dTdu1 = -phi*val_Mean_PrimVar[1];
  const su2double R_dTdu2 = -phi*val_Mean_PrimVar[2];

  const su2double heat_flux_factor = val_laminar_viscosity/Prandtl_Lam + val_eddy_viscosity/Prandtl_Turb;
  const su2double cpoR = Gamma/Gamma_Minus_One; // cp over R
  const su2double conductivity_over_Rd = cpoR*heat_flux_factor/val_dist_ij;

  heat_flux_jac_i[0] = conductivity_over_Rd * R_dTdu0;
  heat_flux_jac_i[1] = conductivity_over_Rd * R_dTdu1;
  heat_flux_jac_i[2] = conductivity_over_Rd * R_dTdu2;

  if (nDim == 2) {

    const su2double R_dTdu3 = phi;
    heat_flux_jac_i[3] = conductivity_over_Rd * R_dTdu3;

  } else {

    const su2double R_dTdu3 = -phi*val_Mean_PrimVar[3];
    const su2double R_dTdu4 = phi;
    heat_flux_jac_i[3] = conductivity_over_Rd * R_dTdu3;
    heat_flux_jac_i[4] = conductivity_over_Rd * R_dTdu4;

  }
}


CGeneralAvgGrad_Flow::CGeneralAvgGrad_Flow(unsigned short val_nDim,
                                           unsigned short val_nVar,
                                           bool val_correct_grad,
                                           CConfig *config)
    : CAvgGrad_Base(val_nDim, val_nVar, val_nDim+4, val_correct_grad, config) {
  
  Mean_SecVar = new su2double [2];
  
}

CGeneralAvgGrad_Flow::~CGeneralAvgGrad_Flow(void) {
  
  delete [] Mean_SecVar;
  
}

void CGeneralAvgGrad_Flow::SetHeatFluxVector(const su2double* const *val_gradprimvar,
                                             const su2double val_laminar_viscosity,
                                             const su2double val_eddy_viscosity,
                                             const su2double val_thermal_conductivity,
                                             const su2double val_heat_capacity_cp) {

  const su2double heat_flux_factor = val_thermal_conductivity + val_heat_capacity_cp*val_eddy_viscosity/Prandtl_Turb;

  /*--- Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure] ---*/
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    heat_flux_vector[iDim] = heat_flux_factor*val_gradprimvar[0][iDim];
  }
}

void CGeneralAvgGrad_Flow::SetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                                               const su2double *val_Mean_SecVar,
                                               const su2double val_eddy_viscosity,
                                               const su2double val_thermal_conductivity,
                                               const su2double val_heat_capacity_cp,
                                               const su2double val_dist_ij) {
  /* Viscous flux Jacobians for arbitrary equations of state */

  //order of val_mean_primitives: T, vx, vy, vz, P, rho, ht
  //order of secondary:dTdrho_e, dTde_rho

  su2double sqvel = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
  }

  su2double rho = val_Mean_PrimVar[nDim+2];
  su2double P= val_Mean_PrimVar[nDim+1];
  su2double h= val_Mean_PrimVar[nDim+3];
  su2double dTdrho_e= val_Mean_SecVar[0];
  su2double dTde_rho= val_Mean_SecVar[1];

  su2double dTdu0= dTdrho_e + dTde_rho*(-(h-P/rho) + sqvel)*(1/rho);
  su2double dTdu1= dTde_rho*(-val_Mean_PrimVar[1])*(1/rho);
  su2double dTdu2= dTde_rho*(-val_Mean_PrimVar[2])*(1/rho);

  su2double total_conductivity = val_thermal_conductivity + val_heat_capacity_cp*val_eddy_viscosity/Prandtl_Turb;
  su2double factor2 = total_conductivity/val_dist_ij;

  heat_flux_jac_i[0] = factor2*dTdu0;
  heat_flux_jac_i[1] = factor2*dTdu1;
  heat_flux_jac_i[2] = factor2*dTdu2;

  if (nDim == 2) {

    su2double dTdu3= dTde_rho*(1/rho);
    heat_flux_jac_i[3] = factor2*dTdu3;

  } else {

    su2double dTdu3= dTde_rho*(-val_Mean_PrimVar[3])*(1/rho);
    su2double dTdu4= dTde_rho*(1/rho);
    heat_flux_jac_i[3] = factor2*dTdu3;
    heat_flux_jac_i[4] = factor2*dTdu4;

  }

}

void CGeneralAvgGrad_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+9);   AD::SetPreaccIn(V_j, nDim+9);
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(S_i, 4); AD::SetPreaccIn(S_j, 4);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);
  AD::SetPreaccIn(PrimVar_Grad_j, nDim+1, nDim);
  AD::SetPreaccIn(turb_ke_i); AD::SetPreaccIn(turb_ke_j);
  AD::SetPreaccIn(Normal, nDim);

  unsigned short iVar, jVar, iDim;

  /*--- Normalized normal vector ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Mean primitive variables ---*/
  
  for (iVar = 0; iVar < nPrimVar; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }

  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }
  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+5];    Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6];       Eddy_Viscosity_j = V_j[nDim+6];
  Thermal_Conductivity_i = V_i[nDim+7]; Thermal_Conductivity_j = V_j[nDim+7];
  Cp_i = V_i[nDim+8]; Cp_j = V_j[nDim+8];

  /*--- Mean secondary variables ---*/
  
  for (iVar = 0; iVar < 2; iVar++) {
    Mean_SecVar[iVar] = 0.5*(S_i[iVar+2]+S_j[iVar+2]);
  }
  
  /*--- Mean Viscosities and turbulent kinetic energy---*/
  
  Mean_Laminar_Viscosity    = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity       = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke              = 0.5*(turb_ke_i + turb_ke_j);
  Mean_Thermal_Conductivity = 0.5*(Thermal_Conductivity_i + Thermal_Conductivity_j);
  Mean_Cp                   = 0.5*(Cp_i + Cp_j);

  /*--- Mean gradient approximation ---*/
  
  for (iVar = 0; iVar < nDim+1; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
    }
  }

  /*--- Projection of the mean gradient in the direction of the edge ---*/

  if (correct_gradient && dist_ij_2 != 0.0) {
    CorrectGradient(Mean_GradPrimVar, PrimVar_i, PrimVar_j, Edge_Vector,
                    dist_ij_2, nDim+1);
  }
  
  /* --- If using UQ methodology, set Reynolds Stress tensor and perform perturbation--- */

  if (using_uq){
    SetReynoldsStressMatrix(Mean_turb_ke);
    SetPerturbedRSM(Mean_turb_ke, config);
  }

  /*--- Get projected flux tensor ---*/
  
  SetStressTensor(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke,
         Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);

  SetHeatFluxVector(Mean_GradPrimVar, Mean_Laminar_Viscosity,
                    Mean_Eddy_Viscosity, Mean_Thermal_Conductivity, Mean_Cp);

  GetViscousProjFlux(Mean_PrimVar, Normal);
  
  /*--- Update viscous residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
  /*--- Compute the implicit part ---*/
  
  if (implicit) {

    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    } else {
      const su2double dist_ij = sqrt(dist_ij_2);
      SetTauJacobian(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
                     dist_ij, UnitNormal);
      SetHeatFluxJacobian(Mean_PrimVar, Mean_SecVar, Mean_Eddy_Viscosity,
                          Mean_Thermal_Conductivity, Mean_Cp, dist_ij);
      GetViscousProjJacs(Mean_PrimVar, Area, Proj_Flux_Tensor,
                         val_Jacobian_i, val_Jacobian_j);
    }
    
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}

CSourceGravity::CSourceGravity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
}

CSourceGravity::~CSourceGravity(void) { }

void CSourceGravity::ComputeResidual(su2double *val_residual, CConfig *config) {
  unsigned short iVar;
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = 0.0;
  
  /*--- Evaluate the source term  ---*/
  val_residual[nDim] = Volume * U_i[0] * STANDARD_GRAVITY;
  
}

CSourceBodyForce::CSourceBodyForce(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  /*--- Store the pointer to the constant body force vector. ---*/

  Body_Force_Vector = new su2double[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Body_Force_Vector[iDim] = config->GetBody_Force_Vector()[iDim];

}

CSourceBodyForce::~CSourceBodyForce(void) {

  if (Body_Force_Vector != NULL) delete [] Body_Force_Vector;

}

void CSourceBodyForce::ComputeResidual(su2double *val_residual, CConfig *config) {
  
  unsigned short iDim;
  su2double Force_Ref = config->GetForce_Ref();
  
  /*--- Zero the continuity contribution ---*/
  
  val_residual[0] = 0.0;
  
  /*--- Momentum contribution ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[iDim+1] = -Volume * U_i[0] * Body_Force_Vector[iDim] / Force_Ref;
  
  /*--- Energy contribution ---*/
  
  val_residual[nDim+1] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[nDim+1] += -Volume * U_i[iDim+1] * Body_Force_Vector[iDim] / Force_Ref;
  
}

CSourceRotatingFrame_Flow::CSourceRotatingFrame_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
}

CSourceRotatingFrame_Flow::~CSourceRotatingFrame_Flow(void) { }

void CSourceRotatingFrame_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {
  
  unsigned short iDim, iVar, jVar;
  su2double Omega[3] = {0,0,0}, Momentum[3] = {0,0,0};
  
  bool implicit     = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Retrieve the angular velocity vector from config. ---*/
  
  Omega[0] = config->GetRotation_Rate_X(config->GetiZone())/config->GetOmega_Ref();
  Omega[1] = config->GetRotation_Rate_Y(config->GetiZone())/config->GetOmega_Ref();
  Omega[2] = config->GetRotation_Rate_Z(config->GetiZone())/config->GetOmega_Ref();
  
  /*--- Get the momentum vector at the current node. ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Momentum[iDim] = U_i[iDim+1];
  
  /*--- Calculate rotating frame source term as ( Omega X Rho-U ) ---*/
  
  if (nDim == 2) {
    val_residual[0] = 0.0;
    val_residual[1] = (Omega[1]*Momentum[2] - Omega[2]*Momentum[1])*Volume;
    val_residual[2] = (Omega[2]*Momentum[0] - Omega[0]*Momentum[2])*Volume;
    val_residual[3] = 0.0;
  } else {
    val_residual[0] = 0.0;
    val_residual[1] = (Omega[1]*Momentum[2] - Omega[2]*Momentum[1])*Volume;
    val_residual[2] = (Omega[2]*Momentum[0] - Omega[0]*Momentum[2])*Volume;
    val_residual[3] = (Omega[0]*Momentum[1] - Omega[1]*Momentum[0])*Volume;
    val_residual[4] = 0.0;
  }
  
  /*--- Calculate the source term Jacobian ---*/
  
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[iVar][jVar] = 0.0;
    if (nDim == 2) {
      val_Jacobian_i[1][2] = -Omega[2]*Volume;
      val_Jacobian_i[2][1] =  Omega[2]*Volume;
    } else {
      val_Jacobian_i[1][2] = -Omega[2]*Volume;
      val_Jacobian_i[1][3] =  Omega[1]*Volume;
      val_Jacobian_i[2][1] =  Omega[2]*Volume;
      val_Jacobian_i[2][3] = -Omega[0]*Volume;
      val_Jacobian_i[3][1] = -Omega[1]*Volume;
      val_Jacobian_i[3][2] =  Omega[0]*Volume;
    }
  }
  
}

CSourceAxisymmetric_Flow::CSourceAxisymmetric_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
}

CSourceAxisymmetric_Flow::~CSourceAxisymmetric_Flow(void) { }

void CSourceAxisymmetric_Flow::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config) {
  
  su2double yinv, Pressure_i, Enthalpy_i, Velocity_i, sq_vel;
  unsigned short iDim, iVar, jVar;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  if (Coord_i[1] > EPS) {
    
    yinv = 1.0/Coord_i[1];
    
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i = U_i[iDim+1] / U_i[0];
      sq_vel += Velocity_i *Velocity_i;
    }
    
    Pressure_i = (Gamma-1.0)*U_i[0]*(U_i[nDim+1]/U_i[0]-0.5*sq_vel);
    Enthalpy_i = (U_i[nDim+1] + Pressure_i) / U_i[0];
    
    val_residual[0] = yinv*Volume*U_i[2];
    val_residual[1] = yinv*Volume*U_i[1]*U_i[2]/U_i[0];
    val_residual[2] = yinv*Volume*(U_i[2]*U_i[2]/U_i[0]);
    val_residual[3] = yinv*Volume*Enthalpy_i*U_i[2];
    
    if (implicit) {
      Jacobian_i[0][0] = 0.0;
      Jacobian_i[0][1] = 0.0;
      Jacobian_i[0][2] = 1.0;
      Jacobian_i[0][3] = 0.0;
      
      Jacobian_i[1][0] = -U_i[1]*U_i[2]/(U_i[0]*U_i[0]);
      Jacobian_i[1][1] = U_i[2]/U_i[0];
      Jacobian_i[1][2] = U_i[1]/U_i[0];
      Jacobian_i[1][3] = 0.0;
      
      Jacobian_i[2][0] = -U_i[2]*U_i[2]/(U_i[0]*U_i[0]);
      Jacobian_i[2][1] = 0.0;
      Jacobian_i[2][2] = 2*U_i[2]/U_i[0];
      Jacobian_i[2][3] = 0.0;
      
      Jacobian_i[3][0] = -Gamma*U_i[2]*U_i[3]/(U_i[0]*U_i[0]) + (Gamma-1)*U_i[2]*(U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]*U_i[0]);
      Jacobian_i[3][1] = -(Gamma-1)*U_i[2]*U_i[1]/(U_i[0]*U_i[0]);
      Jacobian_i[3][2] = Gamma*U_i[3]/U_i[0] - 1/2*(Gamma-1)*( (U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]) + 2*U_i[2]*U_i[2]/(U_i[0]*U_i[0]) );
      Jacobian_i[3][3] = Gamma*U_i[2]/U_i[0];
      
      for (iVar=0; iVar < nVar; iVar++)
        for (jVar=0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] *= yinv*Volume;
      
    }
    
  }
  
  else {
    
    for (iVar=0; iVar < nVar; iVar++)
      val_residual[iVar] = 0.0;
    
    if (implicit) {
      for (iVar=0; iVar < nVar; iVar++) {
        for (jVar=0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;
      }
    }
    
  }
  
}

CSourceWindGust::CSourceWindGust(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
}

CSourceWindGust::~CSourceWindGust(void) { }

void CSourceWindGust::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {
  
  su2double u_gust, v_gust, du_gust_dx, du_gust_dy, du_gust_dt, dv_gust_dx, dv_gust_dy, dv_gust_dt, smx, smy, se, rho, u, v, p;
  unsigned short GustDir = config->GetGust_Dir(); //Gust direction
  
  u_gust = WindGust_i[0];
  v_gust = WindGust_i[1];
  
  if (GustDir == X_DIR) {
    du_gust_dx = WindGustDer_i[0];
    du_gust_dy = WindGustDer_i[1];
    du_gust_dt = WindGustDer_i[2];
    dv_gust_dx = 0.0;
    dv_gust_dy = 0.0;
    dv_gust_dt = 0.0;
  } else {
    du_gust_dx = 0.0;
    du_gust_dy = 0.0;
    du_gust_dt = 0.0;
    dv_gust_dx = WindGustDer_i[0];
    dv_gust_dy = WindGustDer_i[1];
    dv_gust_dt = WindGustDer_i[2];
    
  }
  
  /*--- Primitive variables at point i ---*/
  u = V_i[1];
  v = V_i[2];
  p = V_i[nDim+1];
  rho = V_i[nDim+2];
  
  /*--- Source terms ---*/
  smx = rho*(du_gust_dt + (u+u_gust)*du_gust_dx + (v+v_gust)*du_gust_dy);
  smy = rho*(dv_gust_dt + (u+u_gust)*dv_gust_dx + (v+v_gust)*dv_gust_dy);
  se = u*smx + v*smy + p*(du_gust_dx + dv_gust_dy);
  
  if (nDim == 2) {
    val_residual[0] = 0.0;
    val_residual[1] = smx*Volume;
    val_residual[2] = smy*Volume;
    val_residual[3] = se*Volume;
  } else {
    SU2_MPI::Error("You should only be in the gust source term in two dimensions", CURRENT_FUNCTION);
  }
  
  /*--- For now the source term Jacobian is just set to zero ---*/
  
  unsigned short iVar, jVar;
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Calculate the source term Jacobian ---*/
  
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[iVar][jVar] = 0.0;
  }
  
}
