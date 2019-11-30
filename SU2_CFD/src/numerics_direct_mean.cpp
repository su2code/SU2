/*!
 * \file numerics_direct_mean.cpp
 * \brief This file contains the numerical methods for compressible flow.
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


#include "../include/numerics_structure.hpp"
#include <limits>

CCentBase_Flow::CCentBase_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
  fix_factor = config->GetCent_Jac_Fix_Factor();

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Allocate required structures ---*/
  Diff_U = new su2double [nVar];
  Diff_Lapl = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  MeanVelocity = new su2double [nDim];
  ProjFlux = new su2double [nVar];
}

CCentBase_Flow::~CCentBase_Flow(void) {
  delete [] Diff_U;
  delete [] Diff_Lapl;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
}

void CCentBase_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                                    CConfig *config) {

  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  bool preacc = SetPreaccInVars();

  if (preacc) {
    AD::SetPreaccIn(Normal, nDim);
    AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);
    AD::SetPreaccIn(Lambda_i);    AD::SetPreaccIn(Lambda_j);
    if (dynamic_grid) {
      AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
    }
  }

  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/
  
  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i  = V_i[nDim+2];                       Density_j  = V_j[nDim+2];
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
  
  if (dynamic_grid) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjGridVel * 0.5*(U_i[iVar] + U_j[iVar]);
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjGridVel;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjGridVel;
      }
    }
  }
  
  /*--- Compute the local spectral radius and the stretching factor ---*/
  
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  /*--- Adjustment due to mesh motion ---*/
  
  if (dynamic_grid) {
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }
  
  /*--- Dissipation term ---*/
  
  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  /*--- Compute differences btw. conservative variables, with a correction for enthalpy ---*/
  
  for (iVar = 0; iVar < nVar-1; iVar++) {
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  }
  Diff_U[nVar-1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;
  
  DissipationTerm(val_residual, val_Jacobian_i, val_Jacobian_j);

  if (preacc) {
    AD::SetPreaccOut(val_residual, nVar);
    AD::EndPreacc();
  }
}

void CCentBase_Flow::ScalarDissipationJacobian(su2double **val_Jacobian_i, su2double **val_Jacobian_j) {

  /*--- n-1 diagonal entries ---*/

  for (iVar = 0; iVar < (nVar-1); iVar++) {
    val_Jacobian_i[iVar][iVar] += fix_factor*cte_0;
    val_Jacobian_j[iVar][iVar] -= fix_factor*cte_1;
  }
  
  /*--- Last row of Jacobian_i ---*/
  
  val_Jacobian_i[nVar-1][0] += fix_factor*cte_0*Gamma_Minus_One*sq_vel_i;
  for (iDim = 0; iDim < nDim; iDim++)
    val_Jacobian_i[nVar-1][iDim+1] -= fix_factor*cte_0*Gamma_Minus_One*Velocity_i[iDim];
  val_Jacobian_i[nVar-1][nVar-1] += fix_factor*cte_0*Gamma;
  
  /*--- Last row of Jacobian_j ---*/
  
  val_Jacobian_j[nVar-1][0] -= fix_factor*cte_1*Gamma_Minus_One*sq_vel_j;
  for (iDim = 0; iDim < nDim; iDim++)
    val_Jacobian_j[nVar-1][iDim+1] += fix_factor*cte_1*Gamma_Minus_One*Velocity_j[iDim];
  val_Jacobian_j[nVar-1][nVar-1] -= fix_factor*cte_1*Gamma;

}

CCentJST_Flow::CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
               CCentBase_Flow(val_nDim, val_nVar, config) {

  /*--- Artifical dissipation parameters ---*/
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_Flow();
  Param_Kappa_4 = config->GetKappa_4th_Flow();

}

CCentJST_Flow::~CCentJST_Flow(void) {

}

void CCentJST_Flow::DissipationTerm(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j) {

  /*--- Compute differences btw. Laplacians ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
  }

  /*--- Compute dissipation coefficients ---*/

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
    
    ScalarDissipationJacobian(val_Jacobian_i, val_Jacobian_j);
  }
}

bool CCentJST_Flow::SetPreaccInVars(void) {
  AD::StartPreacc();
  AD::SetPreaccIn(Sensor_i);  AD::SetPreaccIn(Und_Lapl_i, nVar);
  AD::SetPreaccIn(Sensor_j);  AD::SetPreaccIn(Und_Lapl_j, nVar);
  return true;
}

CCentJST_KE_Flow::CCentJST_KE_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                  CCentBase_Flow(val_nDim, val_nVar, config) {

  /*--- Artifical dissipation parameters ---*/
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_Flow();

}

CCentJST_KE_Flow::~CCentJST_KE_Flow(void) {

}

void CCentJST_KE_Flow::DissipationTerm(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j) {

  /*--- Compute dissipation coefficient ---*/

  sc2 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;

  /*--- Compute viscous part of the residual ---*/

  for (iVar = 0; iVar < nVar; iVar++)
      val_residual[iVar] += Epsilon_2*(Diff_U[iVar])*StretchingFactor*MeanLambda;

  /*--- Jacobian computation ---*/

  if (implicit) {

    cte_0 = Epsilon_2*StretchingFactor*MeanLambda;
    cte_1 = cte_0;

    ScalarDissipationJacobian(val_Jacobian_i, val_Jacobian_j);
  }
}

bool CCentJST_KE_Flow::SetPreaccInVars(void) {
  AD::StartPreacc();
  AD::SetPreaccIn(Sensor_i);  AD::SetPreaccIn(Sensor_j);
  return true;
}

CCentLax_Flow::CCentLax_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
               CCentBase_Flow(val_nDim, val_nVar, config) {

  /*--- Artifical dissipation parameters ---*/
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_Flow();

}

CCentLax_Flow::~CCentLax_Flow(void) {

}

void CCentLax_Flow::DissipationTerm(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j) {

  /*--- Compute dissipation coefficient ---*/

  sc0 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc0*su2double(nDim)/3.0;
  
  /*--- Compute viscous part of the residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;
  
  /*--- Jacobian computation ---*/

  if (implicit) {
    
    cte_0 = Epsilon_0*StretchingFactor*MeanLambda;
    cte_1 = cte_0;
    
    ScalarDissipationJacobian(val_Jacobian_i, val_Jacobian_j);
  }
}

bool CCentLax_Flow::SetPreaccInVars(void) {
  AD::StartPreacc();
  return true;
}

CUpwCUSP_Flow::CUpwCUSP_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  if (config->GetDynamic_Grid() && (SU2_MPI::GetRank() == MASTER_NODE))
    cout << "WARNING: Grid velocities are NOT yet considered by the CUSP scheme." << endl;
  
  /*--- Allocate some structures ---*/
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
}

CUpwCUSP_Flow::~CUpwCUSP_Flow(void) {
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
}

void CUpwCUSP_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                                     CConfig *config) {
  
  unsigned short iDim, iVar;
  su2double Diff_U[5] = {0.0,0.0,0.0,0.0,0.0};
  
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+4);
  AD::SetPreaccIn(V_j, nDim+4);
  
  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/
  
  Pressure_i = V_i[nDim+1];     Pressure_j = V_j[nDim+1];
  Density_i  = V_i[nDim+2];     Density_j  = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];     Enthalpy_j = V_j[nDim+3];
  su2double Energy_i = Enthalpy_i - Pressure_i/Density_i;
  su2double Energy_j = Enthalpy_j - Pressure_j/Density_j;

  su2double sq_vel_i = 0.0, sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += Velocity_j[iDim]*Velocity_j[iDim];
  }

  /*-- Face area and unit normal ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Computes differences of conservative variables, with a correction for the enthalpy ---*/
  
  Diff_U[0] = Density_i - Density_j;
  for (iDim = 0; iDim < nDim; iDim++)
    Diff_U[iDim+1] = Density_i*Velocity_i[iDim] - Density_j*Velocity_j[iDim];
  Diff_U[nVar-1] = Density_i*Enthalpy_i - Density_j*Enthalpy_j;
  
  /*--- Get left and right fluxes ---*/
  
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, UnitNormal, ProjFlux_i);
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, UnitNormal, ProjFlux_j);
  
  /*--- Compute dissipation parameters based on Roe-averaged values ---*/
  
  su2double Beta, Nu_c;
  
  su2double R = sqrt(Density_j/Density_i), ProjVelocity = 0.0, sq_vel = 0.0;
  
  for (iDim = 0; iDim < nDim; iDim++) {
    su2double MeanVel = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1.0);
    ProjVelocity += MeanVel*UnitNormal[iDim];
    sq_vel += MeanVel*MeanVel;
  }
  su2double MeanEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1.0);
  su2double MeanSoundSpeed = sqrt(Gamma_Minus_One*fabs(MeanEnthalpy-0.5*sq_vel));
  
  su2double Mach = ProjVelocity / MeanSoundSpeed;
  
  su2double tmp1 = 0.5*(Gamma+1.0)/Gamma*ProjVelocity;
  su2double tmp2 = sqrt(pow(tmp1-ProjVelocity/Gamma, 2.0) + pow(MeanSoundSpeed,2.0)/Gamma);
  su2double LamdaNeg = tmp1 - tmp2, LamdaPos = tmp1 + tmp2;
  
  if (fabs(Mach) >= 1.0) Beta = Mach/fabs(Mach);
  else if (Mach  >= 0.0) Beta = max(0.0, (ProjVelocity + LamdaNeg)/(ProjVelocity - LamdaNeg));
  else                   Beta =-max(0.0, (ProjVelocity + LamdaPos)/(ProjVelocity - LamdaPos));
  
  if (fabs(Mach) >= 1.0) Nu_c = 0.0;
  else {
    if      (Beta > 0.0) Nu_c =-(1.0+Beta)*LamdaNeg;
    else if (Beta < 0.0) Nu_c = (1.0-Beta)*LamdaPos;
    /*--- Limit the minimum scalar dissipation ---*/
    else Nu_c = max(fabs(ProjVelocity), config->GetEntropyFix_Coeff()*MeanSoundSpeed);
  }  
  
  /*--- Compute the residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = 0.5*((1.0+Beta)*ProjFlux_i[iVar] + (1.0-Beta)*ProjFlux_j[iVar] + Nu_c*Diff_U[iVar])*Area;

  /*--- Jacobian computation ---*/

  if (implicit) {
    
    /*--- Flux average and difference contributions ---*/
    
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5*(1.0+Beta), val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5*(1.0-Beta), val_Jacobian_j);
    
    /*--- Solution difference (scalar dissipation) contribution ---*/
    
    su2double cte_0 = 0.5*Nu_c*Area*config->GetCent_Jac_Fix_Factor();
    
    /*--- n-1 diagonal entries ---*/
    
    for (iVar = 0; iVar < (nVar-1); iVar++) {
      val_Jacobian_i[iVar][iVar] += cte_0;
      val_Jacobian_j[iVar][iVar] -= cte_0;
    }
    
    /*--- Last rows ---*/
    
    val_Jacobian_i[nVar-1][0] += cte_0*Gamma_Minus_One*0.5*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_i[nVar-1][iDim+1] -= cte_0*Gamma_Minus_One*Velocity_i[iDim];
    val_Jacobian_i[nVar-1][nVar-1] += cte_0*Gamma;
    
    val_Jacobian_j[nVar-1][0] -= cte_0*Gamma_Minus_One*0.5*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_j[nVar-1][iDim+1] += cte_0*Gamma_Minus_One*Velocity_j[iDim];
    val_Jacobian_j[nVar-1][nVar-1] -= cte_0*Gamma;
    
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}

CUpwAUSM_Flow::CUpwAUSM_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  if (config->GetDynamic_Grid() && (SU2_MPI::GetRank() == MASTER_NODE))
    cout << "WARNING: Grid velocities are NOT yet considered in AUSM-type schemes." << endl;
  
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
  
  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+4);
  AD::SetPreaccIn(V_j, nDim+4);
  
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
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
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

CUpwAUSMPLUS_SLAU_Base_Flow::CUpwAUSMPLUS_SLAU_Base_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                             CNumerics(val_nDim, val_nVar, config) {
  
  if (config->GetDynamic_Grid() && (SU2_MPI::GetRank() == MASTER_NODE))
    cout << "WARNING: Grid velocities are NOT yet considered in AUSM-type schemes." << endl;
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  UseAccurateJacobian = config->GetUse_Accurate_Jacobians();
  HasAnalyticalDerivatives = false;
  FinDiffStep = 1e-4;
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  psi_i = new su2double [nVar];
  psi_j = new su2double [nVar];
  
  RoeVelocity = new su2double [nDim];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwAUSMPLUS_SLAU_Base_Flow::~CUpwAUSMPLUS_SLAU_Base_Flow(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] psi_i;
  delete [] psi_j;
  
  delete [] RoeVelocity;
  delete [] Lambda;
  delete [] Epsilon;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwAUSMPLUS_SLAU_Base_Flow::ComputeMassAndPressureFluxes(CConfig *config, su2double &mdot, su2double &pressure)
{
  /*--- For schemes that fit in the general form of AUSM+up and SLAU schemes you can inherit from this class
   and implement only the specifics, which should be the face mass flux (per unit area) and the face pressure.
     For implicit solution methods this class can either approximate the flux Jacobians (using those of the Roe
   scheme) or compute accurate ones. This is done either numerically, differentiating "mdot" and "pressure"
   using 1st order finite differences, or analytically if you use this function to set the values of
   "dmdot_dVi/j", "dpres_dVi/j" and set "HasAnalyticalDerivatives" to true in the ctor of the derived class.
     For accurate numerical differentiation "mdot" and "pressure" can be functions of, at most, the velocities,
   pressures, densities, and enthalpies at nodes i/j. This is also the order expected for the partial derivatives
   of "mdot" and "pressure" in "d?_dVi/j" (in case they are known analytically, see the AUSM+up implementation).
  ---*/
}

void CUpwAUSMPLUS_SLAU_Base_Flow::ApproximateJacobian(su2double **val_Jacobian_i, su2double **val_Jacobian_j) {

  unsigned short iDim, iVar, jVar, kVar;
  su2double R, RoeDensity, RoeEnthalpy, RoeSoundSpeed, ProjVelocity, sq_vel, Energy_i, Energy_j;
  
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  Energy_j = Enthalpy_j - Pressure_j/Density_j;

  /*--- Mean Roe variables iPoint and jPoint ---*/
  
  R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  ProjVelocity = 0.0;
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    ProjVelocity += RoeVelocity[iDim]*UnitNormal[iDim];
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
  
  /*--- Compute P and Lambda (do it with the Normal) ---*/

  GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);

  /*--- Flow eigenvalues and Entropy correctors ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;
  Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
  Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
  
  /*--- Compute inverse P ---*/
  GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
  
  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
  GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
  GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
  
  /*--- Roe's Flux approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      su2double Proj_ModJac_Tensor_ij = 0.0;
      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
      val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
      val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
    }
  }

}

void CUpwAUSMPLUS_SLAU_Base_Flow::AccurateJacobian(CConfig *config, su2double **val_Jacobian_i, su2double **val_Jacobian_j) {

  /*--- Compute Jacobians using a mixed (numerical/analytical) formulation ---*/
  
  unsigned short iDim, iVar, jVar;
  
  /*--- If not computed analytically, numerically differentiate the fluxes wrt primitives ---*/
  
  if (!HasAnalyticalDerivatives) {
    
    /*--- Create arrays of pointers to the primitive variables so
     we can loop through and perturb them in a general way. ---*/
    
    su2double *primitives_i[6], *primitives_j[6];
    
    for (iDim = 0; iDim < nDim; ++iDim) {
      primitives_i[iDim] = &Velocity_i[iDim];
      primitives_j[iDim] = &Velocity_j[iDim];
    }
    primitives_i[ nDim ] = &Pressure_i;  primitives_j[ nDim ] = &Pressure_j;
    primitives_i[nDim+1] = &Density_i;   primitives_j[nDim+1] = &Density_j;
    primitives_i[nDim+2] = &Enthalpy_i;  primitives_j[nDim+2] = &Enthalpy_j;
    
    /*--- Initialize the gradient arrays with the negative of the quantity,
     then for forward finite differences we add to it and divide. ---*/
    
    for (iVar = 0; iVar < 6; ++iVar) {
      dmdot_dVi[iVar] = -MassFlux;  dpres_dVi[iVar] = -Pressure;
      dmdot_dVj[iVar] = -MassFlux;  dpres_dVj[iVar] = -Pressure;
    }
    
    for (iVar = 0; iVar < nDim+3; ++iVar) {
      /*--- Perturb side i ---*/
      su2double epsilon = FinDiffStep * max(1.0, fabs(*primitives_i[iVar]));
      *primitives_i[iVar] += epsilon;
      ComputeMassAndPressureFluxes(config, MassFlux, Pressure);
      dmdot_dVi[iVar] += MassFlux;  dpres_dVi[iVar] += Pressure;
      dmdot_dVi[iVar] /= epsilon;   dpres_dVi[iVar] /= epsilon;
      *primitives_i[iVar] -= epsilon;
      
      /*--- Perturb side j ---*/
      epsilon = FinDiffStep * max(1.0, fabs(*primitives_j[iVar]));
      *primitives_j[iVar] += epsilon;
      ComputeMassAndPressureFluxes(config, MassFlux, Pressure);
      dmdot_dVj[iVar] += MassFlux;  dpres_dVj[iVar] += Pressure;
      dmdot_dVj[iVar] /= epsilon;   dpres_dVj[iVar] /= epsilon;
      *primitives_j[iVar] -= epsilon;
    }
  }

  /*--- Differentiation of fluxes wrt conservatives assuming ideal gas ---*/
  
  su2double dmdot_dUi[5], dmdot_dUj[5], dpres_dUi[5], dpres_dUj[5];
  su2double sq_veli = 0.0, sq_velj = 0.0, dHi_drhoi = 0.0, dHj_drhoj = 0.0;
  su2double oneOnRhoi = 1.0/Density_i, oneOnRhoj = 1.0/Density_j;
  
  for (jVar = 0; jVar < nVar; ++jVar) {

    /*--- Partial derivatives of the primitives wrt conservative "jVar" ---*/
    su2double dVi_dUi[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    su2double dVj_dUj[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    if (jVar == 0) { // Density
      for (iDim = 0; iDim < nDim; ++iDim) {
        // -u,v,w / rho
        dVi_dUi[iDim] = -Velocity_i[iDim] * oneOnRhoi;
        dVj_dUj[iDim] = -Velocity_j[iDim] * oneOnRhoj;
        // ||V||^2
        sq_veli += Velocity_i[iDim] * Velocity_i[iDim];
        sq_velj += Velocity_j[iDim] * Velocity_j[iDim];
      }
      dVi_dUi[nDim] = 0.5*Gamma_Minus_One*sq_veli;
      dVj_dUj[nDim] = 0.5*Gamma_Minus_One*sq_velj;
      
      dVi_dUi[nDim+1] = dVj_dUj[nDim+1] = 1.0;
      
      dHi_drhoi = 0.5*(Gamma-2.0)*sq_veli - Gamma*Pressure_i/((Gamma-1.0)*Density_i);
      dHj_drhoj = 0.5*(Gamma-2.0)*sq_velj - Gamma*Pressure_j/((Gamma-1.0)*Density_j);
      dVi_dUi[nDim+2] = dHi_drhoi * oneOnRhoi;
      dVj_dUj[nDim+2] = dHj_drhoj * oneOnRhoj;
    }
    else if (jVar == nVar-1) { // rho*Energy
      dVi_dUi[nDim] = dVj_dUj[nDim] = Gamma_Minus_One;
      dVi_dUi[nDim+2] = Gamma * oneOnRhoi;
      dVj_dUj[nDim+2] = Gamma * oneOnRhoj;
    }
    else { // Momentum
      dVi_dUi[jVar-1] = oneOnRhoi;
      dVj_dUj[jVar-1] = oneOnRhoj;
      
      dVi_dUi[nDim] = -Gamma_Minus_One*Velocity_i[jVar-1];
      dVj_dUj[nDim] = -Gamma_Minus_One*Velocity_j[jVar-1];
      
      dVi_dUi[nDim+2] = dVi_dUi[nDim] * oneOnRhoi;
      dVj_dUj[nDim+2] = dVj_dUj[nDim] * oneOnRhoj;
    }
    
    /*--- Dot product to complete chain rule ---*/
    dmdot_dUi[jVar] = 0.0;  dpres_dUi[jVar] = 0.0;
    dmdot_dUj[jVar] = 0.0;  dpres_dUj[jVar] = 0.0;
    for (iVar = 0; iVar < 6; ++iVar) {
      dmdot_dUi[jVar] += dmdot_dVi[iVar]*dVi_dUi[iVar];
      dpres_dUi[jVar] += dpres_dVi[iVar]*dVi_dUi[iVar];
      dmdot_dUj[jVar] += dmdot_dVj[iVar]*dVj_dUj[iVar];
      dpres_dUj[jVar] += dpres_dVj[iVar]*dVj_dUj[iVar];
    }
  }
  
  /*--- Assemble final Jacobians (assuming phi = |mdot|) ---*/
  
  su2double mdot_hat, psi_hat[5];
  
  if (MassFlux > 0.0) {
    mdot_hat = Area*MassFlux*oneOnRhoi;
    for (iVar = 0; iVar < nVar; ++iVar) psi_hat[iVar] = Area*psi_i[iVar];
  }
  else {
    mdot_hat = Area*MassFlux*oneOnRhoj;
    for (iVar = 0; iVar < nVar; ++iVar) psi_hat[iVar] = Area*psi_j[iVar];
  }
  
  /*--- Contribution from the mass flux derivatives ---*/
  for (iVar = 0; iVar < nVar; ++iVar) {
    for (jVar = 0; jVar < nVar; ++jVar) {
      val_Jacobian_i[iVar][jVar] = psi_hat[iVar] * dmdot_dUi[jVar];
      val_Jacobian_j[iVar][jVar] = psi_hat[iVar] * dmdot_dUj[jVar];
    }
  }
  
  /*--- Contribution from the pressure derivatives ---*/
  for (iDim = 0; iDim < nDim; ++iDim) {
    for (jVar = 0; jVar < nVar; ++jVar) {
      val_Jacobian_i[iDim+1][jVar] += Normal[iDim] * dpres_dUi[jVar];
      val_Jacobian_j[iDim+1][jVar] += Normal[iDim] * dpres_dUj[jVar];
    }
  }
  
  /*--- Contributions from the derivatives of PSI wrt the conservatives ---*/
  if (MassFlux > 0.0) {
    /*--- Velocity terms ---*/
    for (iDim = 0; iDim < nDim; ++iDim) {
      val_Jacobian_i[iDim+1][0]      -= mdot_hat*Velocity_i[iDim];
      val_Jacobian_i[iDim+1][iDim+1] += mdot_hat;
      val_Jacobian_i[nVar-1][iDim+1] -= mdot_hat*Gamma_Minus_One*Velocity_i[iDim];
    }
    /*--- Energy terms ---*/
    val_Jacobian_i[nVar-1][0]      += mdot_hat*dHi_drhoi;
    val_Jacobian_i[nVar-1][nVar-1] += mdot_hat*Gamma;
  }
  else {
    /*--- Velocity terms ---*/
    for (iDim = 0; iDim < nDim; ++iDim) {
      val_Jacobian_j[iDim+1][0]      -= mdot_hat*Velocity_j[iDim];
      val_Jacobian_j[iDim+1][iDim+1] += mdot_hat;
      val_Jacobian_j[nVar-1][iDim+1] -= mdot_hat*Gamma_Minus_One*Velocity_j[iDim];
    }
    /*--- Energy terms ---*/
    val_Jacobian_j[nVar-1][0]      += mdot_hat*dHj_drhoj;
    val_Jacobian_j[nVar-1][nVar-1] += mdot_hat*Gamma;
  }

}

void CUpwAUSMPLUS_SLAU_Base_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  unsigned short iDim, iVar;

  /*--- Space to start preaccumulation ---*/

  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+4);
  AD::SetPreaccIn(V_j, nDim+4);

  /*--- Variables for the general form and primitives for mass flux and pressure calculation.  ---*/
  /*--- F_{1/2} = ||A|| ( 0.5 * mdot * (psi_i+psi_j) - 0.5 * |mdot| * (psi_i-psi_j) + N * pf ) ---*/

  psi_i[0] = 1.0;  psi_j[0] = 1.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    /*--- Velocities ---*/
    Velocity_i[iDim] = psi_i[iDim+1] = V_i[iDim+1];
    Velocity_j[iDim] = psi_j[iDim+1] = V_j[iDim+1];
  }

  /*--- Pressure and Density ---*/
  Pressure_i = V_i[nDim+1];  Pressure_j = V_j[nDim+1];
  Density_i  = V_i[nDim+2];  Density_j  = V_j[nDim+2];

  /*--- Enthalpy ---*/
  Enthalpy_i = psi_i[nVar-1] = V_i[nDim+3];
  Enthalpy_j = psi_j[nVar-1] = V_j[nDim+3];

  /*--- Face area (norm or the normal vector) ---*/

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Mass and pressure fluxes defined by derived classes ---*/

  ComputeMassAndPressureFluxes(config, MassFlux, Pressure);
  DissFlux = fabs(MassFlux);

  val_residual[0] = MassFlux;

  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[iDim+1] = 0.5*MassFlux*(psi_i[iDim+1]+psi_j[iDim+1]) +
                           0.5*DissFlux*(psi_i[iDim+1]-psi_j[iDim+1]) +
                           UnitNormal[iDim]*Pressure;

  val_residual[nVar-1] = 0.5*MassFlux*(psi_i[nVar-1]+psi_j[nVar-1]) +
                         0.5*DissFlux*(psi_i[nVar-1]-psi_j[nVar-1]);

  for (iVar = 0; iVar < nVar; iVar++) val_residual[iVar] *= Area;
  
  /*--- Space to end preaccumulation ---*/
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
  /*--- If required, compute Jacobians, either approximately (Roe) or numerically ---*/
  
  if (!implicit) return;
  
  if (UseAccurateJacobian)
    AccurateJacobian(config, val_Jacobian_i, val_Jacobian_j);
  else
    ApproximateJacobian(val_Jacobian_i, val_Jacobian_j);

}


CUpwAUSMPLUSUP_Flow::CUpwAUSMPLUSUP_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                     CUpwAUSMPLUS_SLAU_Base_Flow(val_nDim, val_nVar, config) {

  HasAnalyticalDerivatives = true;
  Minf = config->GetMach();
  Kp = 0.25;
  Ku = 0.75;
  sigma = 1.0;

  if (Minf < EPS)
    SU2_MPI::Error("AUSM+Up requires a reference Mach number (\"MACH_NUMBER\") greater than 0.", CURRENT_FUNCTION);
}

CUpwAUSMPLUSUP_Flow::~CUpwAUSMPLUSUP_Flow(void) {

}

void CUpwAUSMPLUSUP_Flow::ComputeMassAndPressureFluxes(CConfig *config, su2double &mdot, su2double &pressure) {

  /*--- Projected velocities ---*/

  su2double ProjVelocity_i = 0.0, ProjVelocity_j = 0.0;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }

  /*--- Compute interface speed of sound (aF) ---*/ 

  su2double astarL = sqrt(2.0*(Gamma-1.0)/(Gamma+1.0)*Enthalpy_i);
  su2double astarR = sqrt(2.0*(Gamma-1.0)/(Gamma+1.0)*Enthalpy_j);

  su2double ahatL = astarL*astarL/max(astarL, ProjVelocity_i);
  su2double ahatR = astarR*astarR/max(astarR,-ProjVelocity_j);

  su2double aF = min(ahatL,ahatR);

  /*--- Left and right pressures and Mach numbers ---*/
  
  su2double mLP, betaLP, mRM, betaRM;

  su2double mL = ProjVelocity_i/aF;
  su2double mR = ProjVelocity_j/aF;

  su2double MFsq = 0.5*(mL*mL+mR*mR);
  su2double param1 = max(MFsq, Minf*Minf);
  su2double Mrefsq = min(1.0, param1);

  su2double fa = 2.0*sqrt(Mrefsq)-Mrefsq;

  su2double alpha = 3.0/16.0*(-4.0+5.0*fa*fa);
  su2double beta = 1.0/8.0;

  if (fabs(mL) <= 1.0) {
    su2double p1 = 0.25*(mL+1.0)*(mL+1.0);
    su2double p2 = (mL*mL-1.0)*(mL*mL-1.0);

    mLP = p1 + beta*p2;
    betaLP = p1*(2.0-mL) + alpha*mL*p2;
  }
  else {
    mLP = 0.5*(mL+fabs(mL));
    betaLP = mLP/mL;
  }

  if (fabs(mR) <= 1.0) {
    su2double p1 = 0.25*(mR-1.0)*(mR-1.0);
    su2double p2 = (mR*mR-1.0)*(mR*mR-1.0);

    mRM = -p1 - beta*p2;
    betaRM = p1*(2.0+mR) - alpha*mR*p2;
  }
  else {
    mRM = 0.5*(mR-fabs(mR));
    betaRM = mRM/mR;
  }

  /*--- Pressure and velocity diffusion terms ---*/

  su2double rhoF = 0.5*(Density_i+Density_j);
  su2double Mp = -(Kp/fa)*max((1.0-sigma*MFsq),0.0)*(Pressure_j-Pressure_i)/(rhoF*aF*aF);

  su2double Pu = -Ku*fa*betaLP*betaRM*2.0*rhoF*aF*(ProjVelocity_j-ProjVelocity_i);

  /*--- Finally the fluxes ---*/

  su2double mF = mLP + mRM + Mp;
  mdot = aF * (max(mF,0.0)*Density_i + min(mF,0.0)*Density_j);

  pressure = betaLP*Pressure_i + betaRM*Pressure_j + Pu;
  
  if (!implicit || !UseAccurateJacobian) return;
  
  /*--- Analytical differentiation of the face mass flux and
   pressure (in reverse mode, "?_b" denotes dmot_d?). ---*/
  
  /*--- limited mean Mach number (used in division...) ---*/
  su2double MF = max(numeric_limits<passivedouble>::epsilon(),sqrt(MFsq));
  
  for (int outVar=0; outVar<2; ++outVar) {
    
    su2double aF_b    = 0.0, mF_b    = 0.0, MF_b  = 0.0, rhoF_b = 0.0, fa_b   = 0.0, alpha_b = 0.0,
              rho_i_b = 0.0, rho_j_b = 0.0, p_i_b = 0.0, p_j_b  = 0.0, Vn_i_b = 0.0, Vn_j_b  = 0.0,
              mR_b    = 0.0, mL_b    = 0.0, betaLP_b = 0.0,  betaRM_b = 0.0,  tmp = 0.0;
    
    if (outVar==0) {
      /*--- mdot = ... ---*/
      if (mF > 0.0) {
        aF_b += mF*Density_i;
        mF_b += aF*Density_i;
        rho_i_b += mF*aF;
      }
      else {
        aF_b += mF*Density_j;
        mF_b += aF*Density_j;
        rho_j_b += mF*aF;
      }
      
      /*--- Mp = ... ---*/
      if (sigma*MFsq < 1.0) {
        rhoF_b -= Mp/rhoF * mF_b;
        fa_b -= Mp/fa * mF_b;
        aF_b -= 2.0*Mp/aF * mF_b;
        MF_b += 2.0*sigma*MF*(Kp/fa)*(Pressure_j-Pressure_i)/(rhoF*aF*aF) * mF_b;
        tmp = -(Kp/fa)*(1.0-sigma*MFsq)/(rhoF*aF*aF);
        p_i_b -= tmp * mF_b;
        p_j_b += tmp * mF_b;
      }
      
      /*--- rhoF = ... ---*/
      rho_i_b += 0.5*rhoF_b;  rho_j_b += 0.5*rhoF_b;
      
      /*--- mRM = ... ---*/
      if (fabs(mR) < 1.0) mR_b += (1.0-mR)*(0.5+4.0*beta*mR*(mR+1.0)) * mF_b;
      else if (mR <=-1.0) mR_b += mF_b;
      
      /*--- mLP = ... ---*/
      if (fabs(mL) < 1.0) mL_b += (1.0+mL)*(0.5+4.0*beta*mL*(mL-1.0)) * mF_b;
      else if (mL >= 1.0) mL_b += mF_b;
    }
    else {
      /*--- pressure = ... ---*/
      p_i_b += betaLP;  betaLP_b += Pressure_i;
      p_j_b += betaRM;  betaRM_b += Pressure_j;
      
      /*--- Pu = ... ---*/
      rhoF_b += Pu/rhoF;
      fa_b += Pu/fa;
      aF_b += Pu/aF;
      tmp = -Ku*fa*2.0*rhoF*aF*(ProjVelocity_j-ProjVelocity_i);
      betaLP_b += tmp*betaRM;
      betaRM_b += tmp*betaLP;
      tmp = -Ku*fa*betaLP*betaRM*2.0*rhoF*aF;
      Vn_i_b -= tmp;
      Vn_j_b += tmp;
      
      /*--- rhoF = ... ---*/
      rho_i_b += 0.5*rhoF_b;  rho_j_b += 0.5*rhoF_b;
      
      /*--- betaRM = ... ---*/
      if (fabs(mR) < 1.0) {
        tmp = mR*mR-1.0;
        mR_b += tmp*(0.75-alpha*(5.0*tmp+4.0)) * betaRM_b;
        alpha_b -= mR*tmp*tmp * betaRM_b;
      }
      
      /*--- betaLP = ... ---*/
      if (fabs(mL) < 1.0) {
        tmp = mL*mL-1.0;
        mL_b -= tmp*(0.75-alpha*(5.0*tmp+4.0)) * betaLP_b;
        alpha_b += mL*tmp*tmp * betaLP_b;
      }
      
      /*--- alpha = ... ---*/
      fa_b += 1.875*fa * alpha_b;
    }
    
    /*--- steps shared by both ---*/
    /*--- fa = ... ---*/
    su2double Mref_b = 2.0*(1.0-sqrt(Mrefsq)) * fa_b;
    
    /*--- Mrefsq = ... ---*/
    if (MF < 1.0 && MF > Minf) MF_b += Mref_b;
    
    /*--- MFsq = ... ---*/
    mL_b += 0.5*mL/MF * MF_b;  mR_b += 0.5*mR/MF * MF_b;
    
    /*--- mL/R = ... ---*/
    Vn_i_b += mL_b/aF;  Vn_j_b += mR_b/aF;
    aF_b -= (mL*mL_b+mR*mR_b)/aF;
    
    /*--- aF,ahat,astar = f(H_i,H_j) ---*/
    su2double astar_b = aF_b, H_i_b, H_j_b;
    
    if (ahatL < ahatR) {
      if (astarL <= ProjVelocity_i) {
        tmp = astarL/ProjVelocity_i;
        astar_b *= 2.0*tmp;
        Vn_i_b -= tmp*tmp * aF_b;
      }
      H_i_b = sqrt(0.5*(Gamma-1.0)/((Gamma+1.0)*Enthalpy_i)) * astar_b;
      H_j_b = 0.0;
    }
    else {
      if (astarR <= -ProjVelocity_j) {
        tmp = -astarR/ProjVelocity_j;
        astar_b *= 2.0*tmp;
        Vn_j_b += tmp*tmp * aF_b;
      }
      H_j_b = sqrt(0.5*(Gamma-1.0)/((Gamma+1.0)*Enthalpy_j)) * astar_b;
      H_i_b = 0.0;
    }
    
    /*--- store derivatives ---*/
    su2double *target_i = (outVar==0 ? dmdot_dVi : dpres_dVi),
              *target_j = (outVar==0 ? dmdot_dVj : dpres_dVj);
    target_i[5] = target_j[5] = 0.0;
    
    /*--- ProjVelocity = ... ---*/
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      target_i[iDim] = UnitNormal[iDim] * Vn_i_b;
      target_j[iDim] = UnitNormal[iDim] * Vn_j_b;
    }
    target_i[ nDim ] = p_i_b;   target_j[ nDim ] = p_j_b;
    target_i[nDim+1] = rho_i_b; target_j[nDim+1] = rho_j_b;
    target_i[nDim+2] = H_i_b;   target_j[nDim+2] = H_j_b;
  }
}

CUpwAUSMPLUSUP2_Flow::CUpwAUSMPLUSUP2_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                      CUpwAUSMPLUS_SLAU_Base_Flow(val_nDim, val_nVar, config) {
  
  Minf = config->GetMach();
  Kp = 0.25;
  sigma = 1.0;

  if (Minf < EPS)
    SU2_MPI::Error("AUSM+Up2 requires a reference Mach number (\"MACH_NUMBER\") greater than 0.", CURRENT_FUNCTION);
}

CUpwAUSMPLUSUP2_Flow::~CUpwAUSMPLUSUP2_Flow(void) {

}

void CUpwAUSMPLUSUP2_Flow::ComputeMassAndPressureFluxes(CConfig *config, su2double &mdot, su2double &pressure) {

  /*--- Projected velocities and squared magnitude ---*/

  su2double ProjVelocity_i = 0.0, ProjVelocity_j = 0.0, sq_vel = 0.0;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];

    sq_vel += 0.5*(Velocity_i[iDim]*Velocity_i[iDim] + Velocity_j[iDim]*Velocity_j[iDim]);
  }

  /*--- Compute interface speed of sound (aF) ---*/

  su2double astarL = sqrt(2.0*(Gamma-1.0)/(Gamma+1.0)*Enthalpy_i);
  su2double astarR = sqrt(2.0*(Gamma-1.0)/(Gamma+1.0)*Enthalpy_j);

  su2double ahatL = astarL*astarL/max(astarL, ProjVelocity_i);
  su2double ahatR = astarR*astarR/max(astarR,-ProjVelocity_j);

  su2double aF = min(ahatL,ahatR);

  /*--- Left and right pressure functions and Mach numbers ---*/

  su2double mLP, pLP, mRM, pRM;

  su2double mL = ProjVelocity_i/aF;
  su2double mR = ProjVelocity_j/aF;

  su2double MFsq = 0.5*(mL*mL+mR*mR);
  su2double param1 = max(MFsq, Minf*Minf);
  su2double Mrefsq = min(1.0, param1);

  su2double fa = 2.0*sqrt(Mrefsq)-Mrefsq;

  su2double alpha = 3.0/16.0*(-4.0+5.0*fa*fa);
  su2double beta = 1.0/8.0;

  if (fabs(mL) <= 1.0) {
    su2double p1 = 0.25*(mL+1.0)*(mL+1.0);
    su2double p2 = (mL*mL-1.0)*(mL*mL-1.0);

    mLP = p1 + beta*p2;
    pLP = p1*(2.0-mL) + alpha*mL*p2;
  }
  else {
    mLP = 0.5*(mL+fabs(mL));
    pLP = mLP/mL;
  }

  if (fabs(mR) <= 1.0) {
    su2double p1 = 0.25*(mR-1.0)*(mR-1.0);
    su2double p2 = (mR*mR-1.0)*(mR*mR-1.0);

    mRM = -p1 - beta*p2;
    pRM =  p1*(2.0+mR) - alpha*mR*p2;
  }
  else {
    mRM = 0.5*(mR-fabs(mR));
    pRM = mRM/mR;
  }

  /*--- Mass flux with pressure diffusion term ---*/

  su2double rhoF = 0.5*(Density_i+Density_j);
  su2double Mp = -(Kp/fa)*max((1.0-sigma*MFsq),0.0)*(Pressure_j-Pressure_i)/(rhoF*aF*aF);

  su2double mF = mLP + mRM + Mp;
  mdot = aF * (max(mF,0.0)*Density_i + min(mF,0.0)*Density_j);

  /*--- Modified pressure flux ---*/

  pressure = 0.5*(Pressure_j+Pressure_i) + 0.5*(pLP-pRM)*(Pressure_i-Pressure_j) + sqrt(sq_vel)*(pLP+pRM-1.0)*rhoF*aF;

}

CUpwSLAU_Flow::CUpwSLAU_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config, bool val_low_dissipation) :
               CUpwAUSMPLUS_SLAU_Base_Flow(val_nDim, val_nVar, config) {

  slau_low_diss = val_low_dissipation;
  slau2 = false;
}

CUpwSLAU_Flow::~CUpwSLAU_Flow(void) {

}

void CUpwSLAU_Flow::ComputeMassAndPressureFluxes(CConfig *config, su2double &mdot, su2double &pressure) {

  /*--- Project velocities and speed of sound ---*/

  su2double ProjVelocity_i = 0.0, ProjVelocity_j = 0.0, sq_veli = 0.0, sq_velj = 0.0;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];

    sq_veli += Velocity_i[iDim]*Velocity_i[iDim];
    sq_velj += Velocity_j[iDim]*Velocity_j[iDim];
  }

  su2double Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_veli)));

  su2double Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_velj)));

  /*--- Compute interface speed of sound (aF), and left/right Mach number ---*/

  su2double aF = 0.5 * (SoundSpeed_i + SoundSpeed_j);
  su2double mL = ProjVelocity_i/aF;
  su2double mR = ProjVelocity_j/aF;

  /*--- Smooth function of the local Mach number---*/

  su2double Mach_tilde = min(1.0, (1.0/aF) * sqrt(0.5*(sq_veli+sq_velj)));  
  su2double Chi = pow((1.0 - Mach_tilde),2.0);
  su2double f_rho = -max(min(mL,0.0),-1.0) * min(max(mR,0.0),1.0);

  /*--- Mean normal velocity with density weighting ---*/

  su2double Vn_Mag = (Density_i*fabs(ProjVelocity_i) + Density_j*fabs(ProjVelocity_j)) / (Density_i + Density_j);
  su2double Vn_MagL= (1.0 - f_rho)*Vn_Mag + f_rho*fabs(ProjVelocity_i);
  su2double Vn_MagR= (1.0 - f_rho)*Vn_Mag + f_rho*fabs(ProjVelocity_j);  

  /*--- Mass flux function ---*/

  mdot = 0.5 * (Density_i*(ProjVelocity_i+Vn_MagL) + Density_j*(ProjVelocity_j-Vn_MagR) - (Chi/aF)*(Pressure_j-Pressure_i));

  /*--- Pressure function ---*/

  su2double BetaL, BetaR, Dissipation_ij;

  if (fabs(mL) < 1.0) BetaL = 0.25*(2.0-mL)*pow((mL+1.0),2.0);
  else if (mL >= 0)   BetaL = 1.0;
  else                BetaL = 0.0;

  if (fabs(mR) < 1.0) BetaR = 0.25*(2.0+mR)*pow((mR-1.0),2.0);
  else if (mR >= 0)   BetaR = 0.0;
  else                BetaR = 1.0;

  if (slau_low_diss)
    SetRoe_Dissipation(Dissipation_i, Dissipation_j, Sensor_i, Sensor_j, Dissipation_ij, config);
  else
    Dissipation_ij = 1.0;

  pressure = 0.5*(Pressure_i+Pressure_j) + 0.5*(BetaL-BetaR)*(Pressure_i-Pressure_j);

  if (!slau2) pressure += Dissipation_ij*(1.0-Chi)*(BetaL+BetaR-1.0)*0.5*(Pressure_i+Pressure_j);
  else        pressure += Dissipation_ij*sqrt(0.5*(sq_veli+sq_velj))*(BetaL+BetaR-1.0)*aF*0.5*(Density_i+Density_j);

}

CUpwSLAU2_Flow::CUpwSLAU2_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config, bool val_low_dissipation) :
                CUpwSLAU_Flow(val_nDim, val_nVar, config, val_low_dissipation) {

  /*--- The difference between SLAU and SLAU2 is minimal, so we derive from SLAU and set this flag
   so that the ComputeMassAndPressureFluxes function modifies the pressure according to SLAU2.
   This is safe since this constructor is guaranteed to execute after SLAU's one. ---*/
  slau2 = true;
}

CUpwSLAU2_Flow::~CUpwSLAU2_Flow(void) {

}

CUpwHLLC_Flow::CUpwHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  kappa = config->GetRoe_Kappa();
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
  
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

  if (dynamic_grid) {

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
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
  
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

  if (dynamic_grid) {

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


CUpwRoeBase_Flow::CUpwRoeBase_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config,
                                   bool val_low_dissipation) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
  kappa = config->GetRoe_Kappa(); // 1 is unstable

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  roe_low_dissipation = val_low_dissipation;
  
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Conservatives_i = new su2double [nVar];
  Conservatives_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwRoeBase_Flow::~CUpwRoeBase_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Conservatives_i;
  delete [] Conservatives_j;
  delete [] Lambda;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwRoeBase_Flow::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
/*---
 CUpwRoeBase_Flow::ComputeResidual initializes the residual (flux) and its Jacobians with the standard Roe averaging
 fc_{1/2} = kappa*(fc_i+fc_j)*Normal. It then calls this method, which derived classes specialize, to account for
 the dissipation part.
---*/
}

void CUpwRoeBase_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  unsigned short iVar, jVar, iDim;
  su2double ProjGridVel = 0.0, Energy_i, Energy_j;

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+4); AD::SetPreaccIn(V_j, nDim+4); AD::SetPreaccIn(Normal, nDim);
  if (dynamic_grid) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  if (roe_low_dissipation){
    AD::SetPreaccIn(Sensor_i); AD::SetPreaccIn(Sensor_j);
    AD::SetPreaccIn(Dissipation_i); AD::SetPreaccIn(Dissipation_j);
  }
  
  /*--- Face area (norm or the normal vector) and unit normal ---*/

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_i[iDim] = V_i[iDim+1];
  Pressure_i = V_i[nDim+1];
  Density_i  = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
 
  /*--- Primitive variables at point j ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_j[iDim] = V_j[iDim+1];
  Pressure_j = V_j[nDim+1];
  Density_j  = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  
  /*--- Compute variables that are common to the derived schemes ---*/
  
  /*--- Roe-averaged variables at interface between i & j ---*/
  
  su2double R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  su2double sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  RoeSoundSpeed2 = (Gamma-1)*(RoeEnthalpy-0.5*sq_vel);
  
  /*--- Negative RoeSoundSpeed^2, the jump variables is too large, clear fluxes and exit. ---*/
  
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
  
  /*--- P tensor ---*/
  
  GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
  
  /*--- Projected velocity adjusted for mesh motion ---*/
  
  ProjVelocity = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ProjVelocity += RoeVelocity[iDim]*UnitNormal[iDim];
  
  if (dynamic_grid) {
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    ProjVelocity -= ProjGridVel;
  }
  
  /*--- Flow eigenvalues ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;
  
  Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
  Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
  
  /*--- Apply Mavriplis' entropy correction to eigenvalues ---*/
  
  su2double MaxLambda = fabs(ProjVelocity) + RoeSoundSpeed;
  
  for (iVar = 0; iVar < nVar; iVar++) 
    Lambda[iVar] = max(fabs(Lambda[iVar]), config->GetEntropyFix_Coeff()*MaxLambda);
  
  /*--- Reconstruct conservative variables ---*/
  
  Conservatives_i[0] = Density_i;
  Conservatives_j[0] = Density_j;
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Conservatives_i[iDim+1] = Density_i*Velocity_i[iDim];
    Conservatives_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  Conservatives_i[nDim+1] = Density_i*Energy_i;
  Conservatives_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Compute left and right fluxes ---*/
  
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i);
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j);
  
  /*--- Initialize residual (flux) and Jacobians ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = kappa*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
  
  if (implicit) {
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, kappa, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, kappa, val_Jacobian_j);
  }
  
  /*--- Finalize in children class ---*/
  
  FinalizeResidual(val_residual, val_Jacobian_i, val_Jacobian_j, config);
  
  /*--- Correct for grid motion ---*/
  
  if (dynamic_grid) {
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjGridVel*Area * 0.5*(Conservatives_i[iVar]+Conservatives_j[iVar]);
      
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjGridVel*Area;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjGridVel*Area;
      }
    }
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}

CUpwRoe_Flow::CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config,
              bool val_low_dissipation) : CUpwRoeBase_Flow(val_nDim, val_nVar, config, val_low_dissipation) {}

CUpwRoe_Flow::~CUpwRoe_Flow() {}

void CUpwRoe_Flow::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                    su2double **val_Jacobian_j, CConfig *config) {

  unsigned short iVar, jVar, kVar;
  
  /*--- Compute inverse P tensor ---*/
  GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
  
  /*--- Diference between conservative variables at jPoint and iPoint ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = Conservatives_j[iVar]-Conservatives_i[iVar];
  
  /*--- Low dissipation formulation ---*/
  if (roe_low_dissipation)
    SetRoe_Dissipation(Dissipation_i, Dissipation_j, Sensor_i, Sensor_j, Dissipation_ij, config);
  else
    Dissipation_ij = 1.0;
  
  /*--- Standard Roe "dissipation" ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      su2double Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];

      /*--- Update residual and Jacobians ---*/
      val_residual[iVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area*Dissipation_ij;
      
      if(implicit){
        val_Jacobian_i[iVar][jVar] += (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }

}

CUpwL2Roe_Flow::CUpwL2Roe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                CUpwRoeBase_Flow(val_nDim, val_nVar, config, false) {}

CUpwL2Roe_Flow::~CUpwL2Roe_Flow() {}

void CUpwL2Roe_Flow::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                      su2double **val_Jacobian_j, CConfig *config) {

  /*--- L2Roe: a low dissipation version of Roe's approximate Riemann solver for low Mach numbers. IJNMF 2015 ---*/
  
  unsigned short iVar, jVar, kVar, iDim;
  
  /*--- Clamped Mach number ---*/
  
  su2double M_i = 0.0, M_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    M_i += Velocity_i[iDim]*Velocity_i[iDim];
    M_j += Velocity_j[iDim]*Velocity_j[iDim];
  }
  M_i = sqrt(M_i / fabs(Pressure_i*Gamma/Density_i));
  M_j = sqrt(M_j / fabs(Pressure_j*Gamma/Density_j));
  
  su2double zeta = max(0.05,min(max(M_i,M_j),1.0));
  
  /*--- Compute wave amplitudes (characteristics) ---*/
  
  su2double proj_delta_vel = 0.0, delta_vel[3];
  for (iDim = 0; iDim < nDim; iDim++) {
    delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
    proj_delta_vel += delta_vel[iDim]*UnitNormal[iDim];
  }
  proj_delta_vel *= zeta;
  su2double delta_p = Pressure_j - Pressure_i;
  su2double delta_rho = Density_j - Density_i;
  
  su2double delta_wave[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  if (nDim == 2) {
    delta_wave[0] = delta_rho - delta_p/RoeSoundSpeed2;
    delta_wave[1] = (UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1])*zeta;
    delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
  } else {
    delta_wave[0] = delta_rho - delta_p/RoeSoundSpeed2;
    delta_wave[1] = (UnitNormal[0]*delta_vel[2]-UnitNormal[2]*delta_vel[0])*zeta;
    delta_wave[2] = (UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1])*zeta;
    delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
  }
  
  /*--- Update residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (kVar = 0; kVar < nVar; kVar++)
      val_residual[iVar] -= (1.0-kappa)*Lambda[kVar]*delta_wave[kVar]*P_Tensor[iVar][kVar]*Area;
  
  if (!implicit) return;
  
  /*--- If implicit use the Jacobians of the standard Roe scheme as an approximation ---*/
  
  GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      su2double Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];

      val_Jacobian_i[iVar][jVar] += (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
      val_Jacobian_j[iVar][jVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
    }
  }
  
}

CUpwLMRoe_Flow::CUpwLMRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                CUpwRoeBase_Flow(val_nDim, val_nVar, config, false) {}

CUpwLMRoe_Flow::~CUpwLMRoe_Flow() {}

void CUpwLMRoe_Flow::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                      su2double **val_Jacobian_j, CConfig *config) {

  /*--- Rieper, A low-Mach number fix for Roe's approximate Riemman Solver, JCP 2011 ---*/
  
  unsigned short iVar, jVar, kVar, iDim;
  
  /*--- Clamped Mach number ---*/
  
  su2double M_i = 0.0, M_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    M_i += Velocity_i[iDim]*Velocity_i[iDim];
    M_j += Velocity_j[iDim]*Velocity_j[iDim];
  }
  M_i = sqrt(M_i / fabs(Pressure_i*Gamma/Density_i));
  M_j = sqrt(M_j / fabs(Pressure_j*Gamma/Density_j));
  
  su2double zeta = max(0.05,min(max(M_i,M_j),1.0));
  
  /*--- Compute wave amplitudes (characteristics) ---*/
  
  su2double proj_delta_vel = 0.0, delta_vel[3];
  for (iDim = 0; iDim < nDim; iDim++) {
    delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
    proj_delta_vel += delta_vel[iDim]*UnitNormal[iDim];
  }
  proj_delta_vel *= zeta;
  su2double delta_p = Pressure_j - Pressure_i;
  su2double delta_rho = Density_j - Density_i;
  
  su2double delta_wave[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  if (nDim == 2) {
    delta_wave[0] = delta_rho - delta_p/RoeSoundSpeed2;
    delta_wave[1] = (UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1]);
    delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
  } else {
    delta_wave[0] = delta_rho - delta_p/RoeSoundSpeed2;
    delta_wave[1] = (UnitNormal[0]*delta_vel[2]-UnitNormal[2]*delta_vel[0]);
    delta_wave[2] = (UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1]);
    delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
  }
  
  /*--- Update residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (kVar = 0; kVar < nVar; kVar++)
      val_residual[iVar] -= (1.0-kappa)*Lambda[kVar]*delta_wave[kVar]*P_Tensor[iVar][kVar]*Area;
  
  if (!implicit) return;
  
  /*--- If implicit use the Jacobians of the standard Roe scheme as an approximation ---*/
  
  GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      su2double Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];

      val_Jacobian_i[iVar][jVar] += (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
      val_Jacobian_j[iVar][jVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
    }
  }

}

CUpwGeneralRoe_Flow::CUpwGeneralRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
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
  if (dynamic_grid) {
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
  if (dynamic_grid) {
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
    if (dynamic_grid) {
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
    if (dynamic_grid) {
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

CUpwMSW_Flow::CUpwMSW_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  if (config->GetDynamic_Grid() && (SU2_MPI::GetRank() == MASTER_NODE))
    cout << "WARNING: Grid velocities are NOT yet considered in the MSW scheme." << endl;
  
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
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
  
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
  if (dynamic_grid) {
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
  if (dynamic_grid) {
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
                        - TWO3*total_viscosity*div_vel*delta[iDim][jDim];
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
  
  for (iDim = 0; iDim < 3; iDim++){
    Omega[iDim] = config->GetRotation_Rate(iDim)/config->GetOmega_Ref();
  }
  
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
