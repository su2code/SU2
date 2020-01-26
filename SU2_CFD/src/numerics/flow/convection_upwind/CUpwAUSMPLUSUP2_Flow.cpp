/*!
 * \file CUpwAUSMPLUSUP2_Flow.cpp
 * \brief Implementation of numerics class CUpwAUSMPLUSUP2_Flow.
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

#include "../../../../include/numerics/flow/convection_upwind/CUpwAUSMPLUSUP2_Flow.hpp"

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
