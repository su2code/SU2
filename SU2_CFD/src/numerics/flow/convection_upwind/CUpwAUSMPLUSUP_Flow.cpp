/*!
 * \file CUpwAUSMPLUSUP_Flow.cpp
 * \brief Implementation of numerics class CUpwAUSMPLUSUP_Flow.
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

#include "../../../../include/numerics/flow/convection_upwind/CUpwAUSMPLUSUP_Flow.hpp"

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
