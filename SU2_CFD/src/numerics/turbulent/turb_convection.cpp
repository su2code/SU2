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
  Flux = new su2double [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }

  muscl_kappa = 0.5*config->GetMUSCL_Kappa();
  muscl = val_muscl;
}

CUpwScalar::~CUpwScalar(void) {

  delete [] Flux;
  if (Jacobian_i != nullptr) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete [] Jacobian_i[iVar];
      delete [] Jacobian_j[iVar];
    }
    delete [] Jacobian_i;
    delete [] Jacobian_j;
  }
}

void CUpwScalar::GetMUSCLJac(const su2double val_kappa, su2double **val_Jacobian,
                             const su2double *lim_i, const su2double *lim_j,
                             const su2double *val_density, const su2double *val_density_n) {
  const bool wasActive = AD::BeginPassive();

  const su2double weight = (*val_density)/(*val_density_n);
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    val_Jacobian[iVar][iVar] *= weight*(1.0+val_kappa*(lim_j[iVar]-lim_i[iVar]));


  AD::EndPassive(wasActive);
}

CNumerics::ResidualType<> CUpwScalar::ComputeResidual(const CConfig* config) {

  unsigned short iDim;

  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_i, nVar);  AD::SetPreaccIn(TurbVar_j, nVar);
  if (dynamic_grid) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }

  ExtraADPreaccIn();

  Density_i = V_i[nDim+2];
  Density_j = V_j[nDim+2];

  const su2double R = sqrt(fabs(Density_j/Density_i));
  su2double sq_vel = 0.0;

  q_ij = 0.0;
  a0   = 0.0;
  a1   = 0.0;
  Area = 0.0;
  if (dynamic_grid) {
    for (iDim = 0; iDim < nDim; iDim++) {
      su2double Velocity_i = V_i[iDim+1] - GridVel_i[iDim];
      su2double Velocity_j = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i+Velocity_j)*Normal[iDim];
    }
  }
  else {
    for (iDim = 0; iDim < nDim; iDim++) {
      const su2double RoeVelocity = (R*V_j[iDim+1]+V_i[iDim+1])/(R+1.);
      q_ij += RoeVelocity*Normal[iDim];
      a0   += V_i[iDim+1]*Normal[iDim];
      a1   += V_j[iDim+1]*Normal[iDim];

      sq_vel += RoeVelocity*RoeVelocity;
      Area   += Normal[iDim]*Normal[iDim];
    }
  }

  Area = sqrt(Area);

  const su2double RoeEnthalpy = (R*V_j[nDim+3]+V_i[nDim+3])/(R+1.);
  const su2double RoeTke = (R*TurbVar_j[0]+TurbVar_i[0])/(R+1.);
  const su2double RoeSoundSpeed2 = Gamma_Minus_One*(RoeEnthalpy-0.5*sq_vel-RoeTke);

  const su2double RoeSoundSpeed = sqrt(RoeSoundSpeed2);
  const su2double MaxLambda = config->GetEntropyFix_Coeff()*(fabs(q_ij) + RoeSoundSpeed*Area);

  q_ij = (fabs(q_ij) >= MaxLambda) ? su2double(0.5*fabs(q_ij)) 
                                   : su2double(0.25*(q_ij*q_ij/MaxLambda+MaxLambda));
  a0  *= 0.5;
  a1  *= 0.5;

  FinishResidualCalc(config);

  // if (muscl) {

  //   /*--- Extract nodal values ---*/

  //   const su2double Density_n_i = Vn_i[nDim+2];
  //   const su2double Density_n_j = Vn_j[nDim+2];

  //   /*--- Compute Jacobian wrt extrapolation ---*/

  //   GetMUSCLJac(muscl_kappa, Jacobian_i, Limiter_i, Limiter_j, &Density_i, &Density_n_i);
  //   GetMUSCLJac(muscl_kappa, Jacobian_j, Limiter_j, Limiter_i, &Density_j, &Density_n_j);
  // }
  
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

  Flux[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];

  Jacobian_i[0][0] = a0;
  Jacobian_j[0][0] = a1;
}

CUpwSca_TurbSST::CUpwSca_TurbSST(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 const CConfig* config,
                                 bool val_muscl) :
                 CUpwScalar(val_nDim, val_nVar, config, val_muscl) { }

void CUpwSca_TurbSST::ExtraADPreaccIn() {
  AD::SetPreaccIn(V_i, nDim+4);
  AD::SetPreaccIn(V_j, nDim+4);
}

void CUpwSca_TurbSST::FinishResidualCalc(const CConfig* config) {

  Flux[0] = a0*Density_i*TurbVar_i[0]+a1*Density_j*TurbVar_j[0]-q_ij*(Density_j*TurbVar_j[0]-Density_i*TurbVar_i[0]);
  Flux[1] = a0*Density_i*TurbVar_i[1]+a1*Density_j*TurbVar_j[1]-q_ij*(Density_j*TurbVar_j[1]-Density_i*TurbVar_i[1]);

  Jacobian_i[0][0] = a0+q_ij;  Jacobian_i[0][1] = 0.0;
  Jacobian_i[1][0] = 0.0; Jacobian_i[1][1] = a0+q_ij;

  Jacobian_j[0][0] = a1-q_ij;  Jacobian_j[0][1] = 0.0;
  Jacobian_j[1][0] = 0.0; Jacobian_j[1][1] = a1-q_ij;
}
