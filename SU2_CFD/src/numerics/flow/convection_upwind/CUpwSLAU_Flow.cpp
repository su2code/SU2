/*!
 * \file CUpwSLAU_Flow.cpp
 * \brief Implementation of numerics class CUpwSLAU_Flow.
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

#include "../../../../include/numerics/flow/convection_upwind/CUpwSLAU_Flow.hpp"

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
