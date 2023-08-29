/*!
 * \file ausm_slau.cpp
 * \brief Implementations of the AUSM-family of schemes.
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/flow/convection/ausm_slau.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwAUSMPLUS_SLAU_Base_Flow::CUpwAUSMPLUS_SLAU_Base_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
                             CNumerics(val_nDim, val_nVar, config) {

  if (config->GetDynamic_Grid() && (SU2_MPI::GetRank() == MASTER_NODE))
    cout << "WARNING: Grid velocities are NOT yet considered in AUSM-type schemes." << endl;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  UseAccurateJacobian = config->GetUse_Accurate_Jacobians();
  HasAnalyticalDerivatives = false;
  FinDiffStep = 1e-4;

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  psi_i = new su2double [nVar];
  psi_j = new su2double [nVar];

  Flux = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CUpwAUSMPLUS_SLAU_Base_Flow::~CUpwAUSMPLUS_SLAU_Base_Flow() {

  delete [] psi_i;
  delete [] psi_j;

  delete [] Flux;
  delete [] Lambda;
  delete [] Epsilon;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  delete [] Jacobian_i;
  delete [] Jacobian_j;

}

void CUpwAUSMPLUS_SLAU_Base_Flow::ComputeMassAndPressureFluxes(const CConfig* config, su2double &mdot, su2double &pressure)
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

void CUpwAUSMPLUS_SLAU_Base_Flow::AccurateJacobian(const CConfig* config, su2double **val_Jacobian_i, su2double **val_Jacobian_j) {

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

CNumerics::ResidualType<> CUpwAUSMPLUS_SLAU_Base_Flow::ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
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

  Area = GeometryToolbox::Norm(nDim, Normal);

  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Mass and pressure fluxes defined by derived classes ---*/

  ComputeMassAndPressureFluxes(config, MassFlux, Pressure);
  DissFlux = fabs(MassFlux);

  Flux[0] = MassFlux;

  for (iDim = 0; iDim < nDim; iDim++)
    Flux[iDim+1] = 0.5*MassFlux*(psi_i[iDim+1]+psi_j[iDim+1]) +
                   0.5*DissFlux*(psi_i[iDim+1]-psi_j[iDim+1]) +
                   UnitNormal[iDim]*Pressure;

  Flux[nVar-1] = 0.5*MassFlux*(psi_i[nVar-1]+psi_j[nVar-1]) +
                 0.5*DissFlux*(psi_i[nVar-1]-psi_j[nVar-1]);

  for (iVar = 0; iVar < nVar; iVar++) Flux[iVar] *= Area;

  /*--- Space to end preaccumulation ---*/

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  /*--- If required, compute Jacobians, either approximately (Roe) or numerically ---*/

  if (implicit) {
    if (UseAccurateJacobian)
      AccurateJacobian(config, Jacobian_i, Jacobian_j);
    else
      ApproximateJacobian(Jacobian_i, Jacobian_j);
  }

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

CUpwAUSMPLUSUP_Flow::CUpwAUSMPLUSUP_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
                     CUpwAUSMPLUS_SLAU_Base_Flow(val_nDim, val_nVar, config) {

  HasAnalyticalDerivatives = true;
  Minf = config->GetMach();
  Kp = 0.25;
  Ku = 0.75;
  sigma = 1.0;

  if (Minf < EPS)
    SU2_MPI::Error("AUSM+Up requires a reference Mach number (\"MACH_NUMBER\") greater than 0.", CURRENT_FUNCTION);
}

void CUpwAUSMPLUSUP_Flow::ComputeMassAndPressureFluxes(const CConfig* config, su2double &mdot, su2double &pressure) {

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

CUpwAUSMPLUSUP2_Flow::CUpwAUSMPLUSUP2_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
                      CUpwAUSMPLUS_SLAU_Base_Flow(val_nDim, val_nVar, config) {

  Minf = config->GetMach();
  Kp = 0.25;
  sigma = 1.0;

  if (Minf < EPS)
    SU2_MPI::Error("AUSM+Up2 requires a reference Mach number (\"MACH_NUMBER\") greater than 0.", CURRENT_FUNCTION);
}

void CUpwAUSMPLUSUP2_Flow::ComputeMassAndPressureFluxes(const CConfig* config, su2double &mdot, su2double &pressure) {

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

CUpwSLAU_Flow::CUpwSLAU_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config, bool val_low_dissipation) :
               CUpwAUSMPLUS_SLAU_Base_Flow(val_nDim, val_nVar, config) {

  slau_low_diss = val_low_dissipation;
  slau2 = false;
}

CUpwSLAU2_Flow::CUpwSLAU2_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config, bool val_low_dissipation) :
                CUpwSLAU_Flow(val_nDim, val_nVar, config, val_low_dissipation) {

  /*--- The difference between SLAU and SLAU2 is minimal, so we derive from SLAU and set this flag
   so that the ComputeMassAndPressureFluxes function modifies the pressure according to SLAU2.
   This is safe since this constructor is guaranteed to execute after SLAU's one. ---*/
  slau2 = true;
}

void CUpwSLAU_Flow::ComputeMassAndPressureFluxes(const CConfig* config, su2double &mdot, su2double &pressure) {

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
    Dissipation_ij = GetRoe_Dissipation(Dissipation_i, Dissipation_j, Sensor_i, Sensor_j, config);
  else
    Dissipation_ij = 1.0;

  pressure = 0.5*(Pressure_i+Pressure_j) + 0.5*(BetaL-BetaR)*(Pressure_i-Pressure_j);

  if (!slau2) pressure += Dissipation_ij*(1.0-Chi)*(BetaL+BetaR-1.0)*0.5*(Pressure_i+Pressure_j);
  else        pressure += Dissipation_ij*sqrt(0.5*(sq_veli+sq_velj))*(BetaL+BetaR-1.0)*aF*0.5*(Density_i+Density_j);

}

CUpwAUSM_Flow::CUpwAUSM_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) : CNumerics(val_nDim, val_nVar, config) {

  if (config->GetDynamic_Grid() && (SU2_MPI::GetRank() == MASTER_NODE))
    cout << "WARNING: Grid velocities are NOT yet considered in AUSM-type schemes." << endl;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  Flux = new su2double [nVar];
  Diff_U = new su2double [nVar];
  delta_wave = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CUpwAUSM_Flow::~CUpwAUSM_Flow() {

  delete [] Flux;
  delete [] Diff_U;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  delete [] Jacobian_i;
  delete [] Jacobian_j;

}

CNumerics::ResidualType<> CUpwAUSM_Flow::ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+4);
  AD::SetPreaccIn(V_j, nDim+4);

  /*--- Face area (norm or the normal vector) ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

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

  Flux[0] = 0.5*(mF*((Density_i*SoundSpeed_i)+(Density_j*SoundSpeed_j))-Phi*((Density_j*SoundSpeed_j)-(Density_i*SoundSpeed_i)));
  for (iDim = 0; iDim < nDim; iDim++)
    Flux[iDim+1] = 0.5*(mF*((Density_i*SoundSpeed_i*Velocity_i[iDim])+(Density_j*SoundSpeed_j*Velocity_j[iDim]))
                      -Phi*((Density_j*SoundSpeed_j*Velocity_j[iDim])-(Density_i*SoundSpeed_i*Velocity_i[iDim])))+UnitNormal[iDim]*pF;
  Flux[nVar-1] = 0.5*(mF*((Density_i*SoundSpeed_i*Enthalpy_i)+(Density_j*SoundSpeed_j*Enthalpy_j))-Phi*((Density_j*SoundSpeed_j*Enthalpy_j)-(Density_i*SoundSpeed_i*Enthalpy_i)));

  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] *= Area;

  AD::SetPreaccOut(Flux, nVar);
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
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, Jacobian_j);

    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
        Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}
