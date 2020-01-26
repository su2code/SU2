/*!
 * \file CUpwAUSMPLUS_SLAU_Base_Flow.cpp
 * \brief Implementation of numerics class CUpwAUSMPLUS_SLAU_Base_Flow.
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

#include "../../../../include/numerics/flow/convection_upwind/CUpwAUSMPLUS_SLAU_Base_Flow.hpp"

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

