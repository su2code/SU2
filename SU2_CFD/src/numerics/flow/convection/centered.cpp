/*!
 * \file centered.cpp
 * \brief Implementations of centered schemes.
 * \author F. Palacios, T. Economon
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/flow/convection/centered.hpp"

CCentBase_Flow::CCentBase_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
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
  ProjFlux = new su2double [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CCentBase_Flow::~CCentBase_Flow(void) {
  delete [] Diff_U;
  delete [] Diff_Lapl;
  delete [] ProjFlux;

  if (Jacobian_i != nullptr) {
    for (iVar = 0; iVar < nVar; iVar++) {
      delete [] Jacobian_i[iVar];
      delete [] Jacobian_j[iVar];
    }
    delete [] Jacobian_i;
    delete [] Jacobian_j;
  }
}

CNumerics::ResidualType<> CCentBase_Flow::ComputeResidual(const CConfig* config) {

  su2double U_i[5] = {0.0}, U_j[5] = {0.0};

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

  /*--- Residual of the inviscid flux ---*/

  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, ProjFlux);

  /*--- Jacobians of the inviscid flux, scale = 0.5 because ProjFlux ~ 0.5*(fc_i+fc_j)*Normal ---*/

  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        Jacobian_j[iVar][jVar] = Jacobian_i[iVar][jVar];
  }

  /*--- Adjustment due to grid motion ---*/

  if (dynamic_grid) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

    for (iVar = 0; iVar < nVar; iVar++) {
      ProjFlux[iVar] -= ProjGridVel * 0.5*(U_i[iVar] + U_j[iVar]);
      if (implicit) {
        Jacobian_i[iVar][iVar] -= 0.5*ProjGridVel;
        Jacobian_j[iVar][iVar] -= 0.5*ProjGridVel;
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

  DissipationTerm(ProjFlux, Jacobian_i, Jacobian_j);

  if (preacc) {
    AD::SetPreaccOut(ProjFlux, nVar);
    AD::EndPreacc();
  }

  return ResidualType<>(ProjFlux, Jacobian_i, Jacobian_j);

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

CCentLax_Flow::CCentLax_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
               CCentBase_Flow(val_nDim, val_nVar, config) {

  /*--- Artifical dissipation parameters ---*/
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_Flow();

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

CCentJST_KE_Flow::CCentJST_KE_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
                  CCentBase_Flow(val_nDim, val_nVar, config) {

  /*--- Artifical dissipation parameters ---*/
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_Flow();

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

CCentJST_Flow::CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
               CCentBase_Flow(val_nDim, val_nVar, config) {

  /*--- Artifical dissipation parameters ---*/
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_Flow();
  Param_Kappa_4 = config->GetKappa_4th_Flow();

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

CCentLaxInc_Flow::CCentLaxInc_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) : CNumerics(val_nDim, val_nVar, config) {

  implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  variable_density = (config->GetKind_DensityModel() == VARIABLE);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
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
  Jacobian_i   = new su2double* [nVar];
  Jacobian_j   = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Precon[iVar] = new su2double[nVar];
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }

}

CCentLaxInc_Flow::~CCentLaxInc_Flow(void) {

  delete [] Diff_V;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;

  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Precon[iVar];
  delete [] Precon;

  if (Jacobian_i != nullptr) {
    for (iVar = 0; iVar < nVar; iVar++) {
      delete [] Jacobian_i[iVar];
      delete [] Jacobian_j[iVar];
    }
    delete [] Jacobian_i;
    delete [] Jacobian_j;
  }

}

CNumerics::ResidualType<> CCentLaxInc_Flow::ComputeResidual(const CConfig* config) {

  su2double U_i[5] = {0.0}, U_j[5] = {0.0};
  su2double ProjGridVel = 0.0, ProjVelocity = 0.0;

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

  GetInviscidIncProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanBetaInc2, &MeanEnthalpy, Normal, ProjFlux);

  /*--- Jacobians of the inviscid flux ---*/

  if (implicit) {
    GetInviscidIncProjJac(&MeanDensity, MeanVelocity, &MeanBetaInc2, &MeanCp, &MeanTemperature, &MeandRhodT, Normal, 0.5, Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_j[iVar][jVar] = Jacobian_i[iVar][jVar];
      }
    }
  }

  /*--- Corrections due to grid motion ---*/
  if (dynamic_grid) {

    /*--- Recompute conservative variables ---*/

    U_i[0] = DensityInc_i; U_j[0] = DensityInc_j;
    for (iDim = 0; iDim < nDim; iDim++) {
      U_i[iDim+1] = DensityInc_i*Velocity_i[iDim]; U_j[iDim+1] = DensityInc_j*Velocity_j[iDim];
    }
    U_i[nDim+1] = DensityInc_i*Enthalpy_i; U_j[nDim+1] = DensityInc_j*Enthalpy_j;

    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

    /*--- Residual contributions ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      ProjFlux[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);

      /*--- Jacobian contributions ---*/
      /*--- Implicit terms ---*/
      if (implicit) {
        for (iDim = 0; iDim < nDim; iDim++){
          Jacobian_i[iDim+1][iDim+1] -= 0.5*ProjVelocity*DensityInc_i;
          Jacobian_j[iDim+1][iDim+1] -= 0.5*ProjVelocity*DensityInc_j;
        }
        Jacobian_i[nDim+1][nDim+1] -= 0.5*ProjVelocity*DensityInc_i*Cp_i;
        Jacobian_j[nDim+1][nDim+1] -= 0.5*ProjVelocity*DensityInc_j*Cp_j;
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

  /*--- Projected velocity adjustment due to mesh motion ---*/

  if (dynamic_grid) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    }
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }

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
      ProjFlux[iVar] += Precon[iVar][jVar]*Epsilon_0*Diff_V[jVar]*StretchingFactor*MeanLambda;
      if (implicit) {
        Jacobian_i[iVar][jVar] += Precon[iVar][jVar]*Epsilon_0*StretchingFactor*MeanLambda;
        Jacobian_j[iVar][jVar] -= Precon[iVar][jVar]*Epsilon_0*StretchingFactor*MeanLambda;
      }
    }
  }

  /*--- Remove energy contributions if we aren't solving the energy equation. ---*/

  if (!energy) {
    ProjFlux[nDim+1] = 0.0;
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[iVar][nDim+1] = 0.0;
        Jacobian_j[iVar][nDim+1] = 0.0;

        Jacobian_i[nDim+1][iVar] = 0.0;
        Jacobian_j[nDim+1][iVar] = 0.0;
      }
    }
  }

  return ResidualType<>(ProjFlux, Jacobian_i, Jacobian_j);

}


CCentJSTInc_Flow::CCentJSTInc_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) : CNumerics(val_nDim, val_nVar, config) {

  implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  variable_density = (config->GetKind_DensityModel() == VARIABLE);
  energy           = config->GetEnergy_Equation();
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

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
  Jacobian_i   = new su2double* [nVar];
  Jacobian_j   = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Precon[iVar] = new su2double[nVar];
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }

}

CCentJSTInc_Flow::~CCentJSTInc_Flow(void) {

  delete [] Diff_V;
  delete [] Diff_Lapl;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;

  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Precon[iVar];
  delete [] Precon;

  if (Jacobian_i != nullptr) {
    for (iVar = 0; iVar < nVar; iVar++) {
      delete [] Jacobian_i[iVar];
      delete [] Jacobian_j[iVar];
    }
    delete [] Jacobian_i;
    delete [] Jacobian_j;
  }

}

CNumerics::ResidualType<> CCentJSTInc_Flow::ComputeResidual(const CConfig* config) {

  su2double U_i[5] = {0.0}, U_j[5] = {0.0};
  su2double ProjGridVel = 0.0;

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

  GetInviscidIncProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanBetaInc2, &MeanEnthalpy, Normal, ProjFlux);

  /*--- Jacobians of the inviscid flux ---*/

  if (implicit) {
    GetInviscidIncProjJac(&MeanDensity, MeanVelocity, &MeanBetaInc2, &MeanCp, &MeanTemperature, &MeandRhodT, Normal, 0.5, Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_j[iVar][jVar] = Jacobian_i[iVar][jVar];
      }
    }
  }

  /*--- Corrections due to grid motion ---*/
  if (dynamic_grid) {

    /*--- Recompute conservative variables ---*/

    U_i[0] = DensityInc_i; U_j[0] = DensityInc_j;
    for (iDim = 0; iDim < nDim; iDim++) {
      U_i[iDim+1] = DensityInc_i*Velocity_i[iDim]; U_j[iDim+1] = DensityInc_j*Velocity_j[iDim];
    }
    U_i[nDim+1] = DensityInc_i*Enthalpy_i; U_j[nDim+1] = DensityInc_j*Enthalpy_j;

    su2double ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

    /*--- Residual contributions ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      ProjFlux[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);

      /*--- Jacobian contributions ---*/
      /*--- Implicit terms ---*/
      if (implicit) {
        for (iDim = 0; iDim < nDim; iDim++){
          Jacobian_i[iDim+1][iDim+1] -= 0.5*ProjVelocity*DensityInc_i;
          Jacobian_j[iDim+1][iDim+1] -= 0.5*ProjVelocity*DensityInc_j;
        }
        Jacobian_i[nDim+1][nDim+1] -= 0.5*ProjVelocity*DensityInc_i*Cp_i;
        Jacobian_j[nDim+1][nDim+1] -= 0.5*ProjVelocity*DensityInc_j*Cp_j;
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

  /*--- Projected velocity adjustment due to mesh motion ---*/

  if (dynamic_grid) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    }
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }

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
      ProjFlux[iVar] += Precon[iVar][jVar]*(Epsilon_2*Diff_V[jVar] - Epsilon_4*Diff_Lapl[jVar])*StretchingFactor*MeanLambda;
      if (implicit) {
        Jacobian_i[iVar][jVar] += Precon[iVar][jVar]*(Epsilon_2 + Epsilon_4*su2double(Neighbor_i+1))*StretchingFactor*MeanLambda;
        Jacobian_j[iVar][jVar] -= Precon[iVar][jVar]*(Epsilon_2 + Epsilon_4*su2double(Neighbor_j+1))*StretchingFactor*MeanLambda;
      }
    }
  }

  /*--- Remove energy contributions if not solving the energy equation. ---*/

  if (!energy) {
    ProjFlux[nDim+1] = 0.0;
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[iVar][nDim+1] = 0.0;
        Jacobian_j[iVar][nDim+1] = 0.0;

        Jacobian_i[nDim+1][iVar] = 0.0;
        Jacobian_j[nDim+1][iVar] = 0.0;
      }
    }
  }

  return ResidualType<>(ProjFlux, Jacobian_i, Jacobian_j);

}
