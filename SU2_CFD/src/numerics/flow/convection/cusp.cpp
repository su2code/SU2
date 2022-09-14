/*!
 * \file cusp.cpp
 * \brief Implementation of the CUSP scheme.
 * \author F. Palacios, T. Economon
 * \version 7.3.1 "Blackbird"
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

#include "../../../../include/numerics/flow/convection/cusp.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwCUSP_Flow::CUpwCUSP_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  if (config->GetDynamic_Grid() && (SU2_MPI::GetRank() == MASTER_NODE))
    cout << "WARNING: Grid velocities are NOT yet considered by the CUSP scheme." << endl;

  /*--- Allocate some structures ---*/
  Flux = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CUpwCUSP_Flow::~CUpwCUSP_Flow(void) {
  delete [] Flux;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Jacobian_i;
  delete [] Jacobian_j;
}

CNumerics::ResidualType<> CUpwCUSP_Flow::ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  unsigned short iDim, iVar;
  su2double Diff_U[5] = {0.0};

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

  Area = GeometryToolbox::Norm(nDim, Normal);

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
    Flux[iVar] = 0.5*((1.0+Beta)*ProjFlux_i[iVar] + (1.0-Beta)*ProjFlux_j[iVar] + Nu_c*Diff_U[iVar])*Area;

  /*--- Jacobian computation ---*/

  if (implicit) {

    /*--- Flux average and difference contributions ---*/

    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5*(1.0+Beta), Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5*(1.0-Beta), Jacobian_j);

    /*--- Solution difference (scalar dissipation) contribution ---*/

    su2double cte_0 = 0.5*Nu_c*Area*config->GetCent_Jac_Fix_Factor();

    /*--- n-1 diagonal entries ---*/

    for (iVar = 0; iVar < (nVar-1); iVar++) {
      Jacobian_i[iVar][iVar] += cte_0;
      Jacobian_j[iVar][iVar] -= cte_0;
    }

    /*--- Last rows ---*/

    Jacobian_i[nVar-1][0] += cte_0*Gamma_Minus_One*0.5*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      Jacobian_i[nVar-1][iDim+1] -= cte_0*Gamma_Minus_One*Velocity_i[iDim];
    Jacobian_i[nVar-1][nVar-1] += cte_0*Gamma;

    Jacobian_j[nVar-1][0] -= cte_0*Gamma_Minus_One*0.5*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      Jacobian_j[nVar-1][iDim+1] += cte_0*Gamma_Minus_One*Velocity_j[iDim];
    Jacobian_j[nVar-1][nVar-1] -= cte_0*Gamma;

  }

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}
