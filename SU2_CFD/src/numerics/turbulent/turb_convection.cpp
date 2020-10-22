/*!
 * \file turb_convection.cpp
 * \brief Implementation of numerics classes to compute convective
 *        fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.0.3 "Blackbird"
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

#include "../../../include/numerics/turbulent/turb_convection.hpp"

CUpwScalar::CUpwScalar(unsigned short val_nDim,
                       unsigned short val_nVar,
                       const CConfig* config,
                       bool val_muscl) :
  CNumerics(val_nDim, val_nVar, config),
  incompressible(config->GetKind_Regime() == INCOMPRESSIBLE),
  dynamic_grid(config->GetDynamic_Grid())
{
  Flux = new su2double [nVar] ();
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (auto iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar] ();
    Jacobian_j[iVar] = new su2double [nVar] ();
  }
}

CUpwScalar::~CUpwScalar(void) {

  delete [] Flux;
  for (auto iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Jacobian_i;
  delete [] Jacobian_j;
}

CNumerics::ResidualType<> CUpwScalar::ComputeResidual(const CConfig* config) {

  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_i, nVar);  AD::SetPreaccIn(TurbVar_j, nVar);
  if (dynamic_grid) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }

  ExtraADPreaccIn();

  Area = 0.0;
  for (auto iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  for (auto iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Primitive variables ---*/

  for (auto iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
  }

  Pressure_i = V_i[nDim+1];
  Pressure_j = V_j[nDim+1];

  Density_i  = V_i[nDim+2];
  Density_j  = V_j[nDim+2];

  turb_ke_i = tkeNeeded ? TurbVar_i[0] : 0.0;
  turb_ke_j = tkeNeeded ? TurbVar_j[0] : 0.0;

  Energy_i = Pressure_i/(Gamma_Minus_One*Density_i)+turb_ke_i;
  Energy_j = Pressure_j/(Gamma_Minus_One*Density_j)+turb_ke_j;
  for (auto iDim = 0; iDim < nDim; iDim++) {
    Energy_i += 0.5*pow(Velocity_i[iDim],2);
    Energy_j += 0.5*pow(Velocity_j[iDim],2);
  }

  Enthalpy_i = Energy_i + Pressure_i/Density_i;
  Enthalpy_j = Energy_j + Pressure_j/Density_j;

  SoundSpeed_i = sqrt(fabs(Pressure_i*Gamma/Density_i));
  SoundSpeed_j = sqrt(fabs(Pressure_j*Gamma/Density_j));
    
  /*--- Compute variables that are common to the derived schemes ---*/

  /*--- Roe-averaged variables at interface between i & j ---*/

  R = sqrt(fabs(Density_j/Density_i));
  inv_R_Plus_One = 1./(R+1.);

  Lambda[0] = 0.0;
  ProjVel_i = 0.0;
  ProjVel_j = 0.0;
  RoeSqVel  = 0.0;
  if (dynamic_grid) {
    for (auto iDim = 0; iDim < nDim; iDim++) {
      su2double Velocity_i = V_i[iDim+1] - GridVel_i[iDim];
      su2double Velocity_j = V_j[iDim+1] - GridVel_j[iDim];
      Lambda[0] += 0.5*(Velocity_i+Velocity_j)*UnitNormal[iDim];
    }
  }
  else {
    for (auto iDim = 0; iDim < nDim; iDim++) {
      const su2double RoeVelocity = (R*Velocity_j[iDim]+Velocity_i[iDim])*inv_R_Plus_One;
      Lambda[0] += RoeVelocity*UnitNormal[iDim];
      ProjVel_i += Velocity_i[iDim]*UnitNormal[iDim];
      ProjVel_j += Velocity_j[iDim]*UnitNormal[iDim];
      RoeSqVel  += RoeVelocity*RoeVelocity;
    }
  }

  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)*inv_R_Plus_One;
  RoeTke = (R*turb_ke_j+turb_ke_i)*inv_R_Plus_One;
  RoeSoundSpeed2 = Gamma_Minus_One*(RoeEnthalpy-0.5*RoeSqVel-RoeTke);

  /*--- Negative RoeSoundSpeed^2, the jump variables are too large, clear fluxes and exit ---*/

  if (RoeSoundSpeed2 <= 0.0) {
    for (auto iVar = 0; iVar < nVar; iVar++) {
      Flux[iVar] = 0.0;
      for (auto jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] = 0.0;
        Jacobian_j[iVar][jVar] = 0.0;
      }
    }

    AD::SetPreaccOut(Flux, nVar);
    AD::EndPreacc();

    return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
  }

  /*--- Positive RoeSoundSpeed^2, compute fluxes and Jacobian ---*/

  FinishResidualCalc(config);
  
  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

CUpwSca_TurbSA::CUpwSca_TurbSA(unsigned short val_nDim,
                               unsigned short val_nVar,
                               const CConfig* config,
                               bool val_muscl) :
                CUpwScalar(val_nDim, val_nVar, config, val_muscl) { }

void CUpwSca_TurbSA::ExtraADPreaccIn() {
  AD::SetPreaccIn(V_i, nDim+1);
  AD::SetPreaccIn(V_j, nDim+1);
}

void CUpwSca_TurbSA::FinishResidualCalc(const CConfig* config) {

  Flux[0] = 0.5*(ProjVel_i*TurbVar_i[0]+ProjVel_j*TurbVar_j[0])*Area;

  Jacobian_i[0][0] = 0.5*ProjVel_i*Area;
  Jacobian_j[0][0] = 0.5*ProjVel_j*Area;
}

CUpwSca_TurbSST::CUpwSca_TurbSST(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 const CConfig* config,
                                 bool val_muscl) :
                 CUpwScalar(val_nDim, val_nVar, config, val_muscl) { }

void CUpwSca_TurbSST::ExtraADPreaccIn() {
  AD::SetPreaccIn(V_i, nDim+3);
  AD::SetPreaccIn(V_j, nDim+3);
}

void CUpwSca_TurbSST::FinishResidualCalc(const CConfig* config) {

  RoeOmega = (R*TurbVar_j[1]+TurbVar_i[1])*inv_R_Plus_One;
  RoeSoundSpeed = sqrt(RoeSoundSpeed2);

  const su2double dir = 1.0 - 2.0*(Lambda[0] < 0);

  Lambda[0] = Lambda[0];
  Lambda[1] = Lambda[0] + RoeSoundSpeed;
  Lambda[2] = Lambda[0] - RoeSoundSpeed;

  /*--- Harten and Hyman (1983) entropy correction ---*/

  Epsilon[0] = 4.0*max((Lambda[0]-ProjVel_i), 
                       (ProjVel_j-Lambda[0]));
  Epsilon[1] = 4.0*max((Lambda[1]-(ProjVel_i+SoundSpeed_i)),
                       ((ProjVel_j+SoundSpeed_j)-Lambda[1]));
  Epsilon[2] = 4.0*max((Lambda[2]-(ProjVel_i-SoundSpeed_i)),
                       ((ProjVel_j-SoundSpeed_j)-Lambda[2]));

  for (auto iVar = 0; iVar < 3; iVar++) {
    Epsilon[iVar] = max(Epsilon[iVar], 0.0);
    Lambda[iVar] = (fabs(Lambda[iVar]) < Epsilon[iVar]) ? su2double(0.5*(Lambda[iVar]*Lambda[iVar]/Epsilon[iVar] + Epsilon[iVar]))
                                                        : su2double(fabs(Lambda[iVar]));
  }

  /*--- Fluxes ---*/

  const su2double rkv_i = Density_i*TurbVar_i[0]*ProjVel_i;
  const su2double rkv_j = Density_j*TurbVar_j[0]*ProjVel_j;

  const su2double rov_i = Density_i*TurbVar_i[1]*ProjVel_i;
  const su2double rov_j = Density_j*TurbVar_j[1]*ProjVel_j;

  /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/

  const su2double Diff_rk = Density_j*TurbVar_j[0]-Density_i*TurbVar_i[0];
  const su2double Diff_ro = Density_j*TurbVar_j[1]-Density_i*TurbVar_i[1];

  // const su2double Lambda_G = (Lambda[0]-0.5*Lambda[1]-0.5*Lambda[2])*(Gamma - FIVE3)/RoeSoundSpeed2;
  const su2double Lambda_G = (Lambda[0]-0.5*Lambda[1]-0.5*Lambda[2])*Gamma_Minus_One/RoeSoundSpeed2; 

  const su2double Diss_rk_rk = Lambda[0]+RoeTke*Lambda_G;
  const su2double Diss_ro_ro = Lambda[0];
  const su2double Diss_ro_rk = RoeOmega*Lambda_G;

  Flux[0] = 0.5*(rkv_i+rkv_j-Diss_rk_rk*Diff_rk)*Area;
  Flux[1] = 0.5*(rov_i+rov_j-Diss_ro_rk*Diff_rk
                            -Diss_ro_ro*Diff_ro)*Area;

  Jacobian_i[0][0] = 0.5*(ProjVel_i+Diss_rk_rk)*Area;
  Jacobian_j[0][0] = 0.5*(ProjVel_j-Diss_rk_rk)*Area;

  Jacobian_i[1][0] =  0.5*Diss_ro_rk*Area;
  Jacobian_j[1][0] = -0.5*Diss_ro_rk*Area;

  Jacobian_i[1][1] = 0.5*(ProjVel_i+Diss_ro_ro)*Area;
  Jacobian_j[1][1] = 0.5*(ProjVel_j-Diss_ro_ro)*Area;
}
