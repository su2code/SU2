/*!
 * \file slau.cpp
 * \brief Implementations of the SLAU-family of schemes in NEMO.
 * \author W. Maier
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/NEMO/convection/slau.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwSLAU_NEMO::CUpwSLAU_NEMO(unsigned short val_nDim, unsigned short val_nVar,
                                                         unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
                                                         CConfig *config) : CNEMONumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad,
                                                         config) {

  rhos_i = new su2double [nSpecies];
  rhos_j = new su2double [nSpecies];
  mF_s   = new su2double [nSpecies];
  u_i    = new su2double [nDim];
  u_j    = new su2double [nDim];
  Flux   = new su2double [nVar];

}

CUpwSLAU_NEMO::~CUpwSLAU_NEMO(void) {
  
  delete [] rhos_i;
  delete [] rhos_j;
  delete [] mF_s;
  delete [] u_i;
  delete [] u_j;
  delete [] Flux;
}

CNumerics::ResidualType<> CUpwSLAU_NEMO::ComputeResidual(const CConfig *config) {

  unsigned short iDim, iVar, iSpecies;
  su2double rho_i, rho_j,
      e_ve_i, e_ve_j, mL, mR, mF, pF;

  /*--- Compute geometric quantities ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Pull stored primitive variables ---*/
  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, a, c]
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    rhos_i[iSpecies] = V_i[RHOS_INDEX+iSpecies];
    rhos_j[iSpecies] = V_j[RHOS_INDEX+iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    u_i[iDim] = V_i[VEL_INDEX+iDim];
    u_j[iDim] = V_j[VEL_INDEX+iDim];
  }

  P_i       = V_i[P_INDEX];   P_j       = V_j[P_INDEX];
  h_i       = V_i[H_INDEX];   h_j       = V_j[H_INDEX];
  a_i       = V_i[A_INDEX];   a_j       = V_j[A_INDEX];
  rho_i     = V_i[RHO_INDEX]; rho_j     = V_j[RHO_INDEX];

  e_ve_i  = 0; e_ve_j  = 0;
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    e_ve_i += (V_i[RHOS_INDEX+iSpecies]*eve_i[iSpecies])/rho_i;
    e_ve_j += (V_j[RHOS_INDEX+iSpecies]*eve_j[iSpecies])/rho_j;
  }

  /*--- Projected and squared velocities ---*/
  ProjVel_i = GeometryToolbox::DotProduct(nDim, u_i, UnitNormal);
  ProjVel_j = GeometryToolbox::DotProduct(nDim, u_j, UnitNormal);
  
  su2double sq_vel_i = GeometryToolbox::SquaredNorm(nDim,u_i);
  su2double sq_vel_j = GeometryToolbox::SquaredNorm(nDim,u_j);

  /*--- Calculate interface soundspeed and L/R Mach numbers ---*/
  aF = 0.5*(a_i + a_j);
  mL = ProjVel_i/a_i;
  mR = ProjVel_j/a_j;

  /*--- Smooth function of the local Mach number---*/

  Mach_tilde = min(1.0, (1.0/aF) * sqrt(0.5*(sq_vel_i+sq_vel_j)));
  Chi = pow((1.0 - Mach_tilde),2.0);
  f_rho = -max(min(mL,0.0),-1.0) * min(max(mR,0.0),1.0);

  /*--- Mean normal velocity with density weighting ---*/
  Vn_Mag = (rho_i*fabs(ProjVel_i) + rho_j*fabs(ProjVel_j)) / (rho_i + rho_j);
  Vn_MagL= (1.0 - f_rho)*Vn_Mag + f_rho*fabs(ProjVel_i);
  Vn_MagR= (1.0 - f_rho)*Vn_Mag + f_rho*fabs(ProjVel_j);

  mF = 0.5 * (rho_i * (ProjVel_i + Vn_MagL) + rho_j * (ProjVel_j - Vn_MagR) - (Chi/aF)*(P_j-P_i));

  /*--- Species extension mass flux function ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++ ){
    su2double Ys_i = rhos_i[iSpecies]/rho_i; 
    su2double Ys_j = rhos_j[iSpecies]/rho_j;
    mF_s[iSpecies]  = 0.5*mF*(Ys_i+Ys_j);
  }

  /*--- Pressure function ---*/

  if (fabs(mL) < 1.0) BetaL = 0.25*(2.0-mL)*pow((mL+1.0),2.0);
  else if (mL >= 0)   BetaL = 1.0;
  else                BetaL = 0.0;
  
  if (fabs(mR) < 1.0) BetaR = 0.25*(2.0+mR)*pow((mR-1.0),2.0);
  else if (mR >= 0 )  BetaR = 0.0;
  else                BetaR = 1.0;
  
  //Roe Dissipation not implemented
  Dissipation_ij = 1.0;
  pF = 0.5 * (P_i + P_j) + 0.5 * (BetaL - BetaR) * (P_i - P_j)
      + Dissipation_ij*(1.0 - Chi) * (BetaL + BetaR - 1.0) *  0.5 * (P_i + P_j);

  //TODO this could be dumb.....should just be mF_s???
  for (iSpecies=0;iSpecies<nSpecies;iSpecies++){
    Flux[iSpecies] = 0.5*(mF_s[iSpecies]+fabs(mF_s[iSpecies])) +
                     0.5*(mF_s[iSpecies]-fabs(mF_s[iSpecies]));
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    Flux[nSpecies+iDim] = 0.5*(mF+fabs(mF)) * u_i[iDim];
    Flux[nSpecies+iDim]+= 0.5*(mF-fabs(mF)) * u_j[iDim] ;
    Flux[nSpecies+iDim]+= pF*UnitNormal[iDim];
  }
  Flux[nSpecies+nDim]   = 0.5*(mF+fabs(mF))*(h_i) + 0.5*(mF-fabs(mF))*(h_j);
  Flux[nSpecies+nDim+1] = 0.5*(mF+fabs(mF))*(e_ve_i) + 0.5*(mF-fabs(mF))*(e_ve_j);

  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] *= Area;

  if (implicit){
    SU2_MPI::Error("NEMO SLAU: Impicit not operational.", CURRENT_FUNCTION);
  }

  return ResidualType<>(Flux, nullptr, nullptr);
}
