/*!
 * \file turb_convection.cpp
 * \brief Implementation of numerics classes to compute convective
 *        fluxes in turbulence problems.
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

#include "../../../include/numerics/turbulent/turb_convection.hpp"

CUpwScalar::CUpwScalar(unsigned short val_nDim,
                       unsigned short val_nVar,
                       const CConfig* config) :
  CNumerics(val_nDim, val_nVar, config),
  implicit(config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT),
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

  q_ij = 0.0;
  if (dynamic_grid) {
    for (iDim = 0; iDim < nDim; iDim++) {
      su2double Velocity_i = V_i[iDim+1] - GridVel_i[iDim];
      su2double Velocity_j = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i+Velocity_j)*Normal[iDim];
    }
  }
  else {
    for (iDim = 0; iDim < nDim; iDim++) {
      q_ij += 0.5*(V_i[iDim+1]+V_j[iDim+1])*Normal[iDim];
    }
  }

  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));

  FinishResidualCalc(config);

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

CUpwSca_TurbSA::CUpwSca_TurbSA(unsigned short val_nDim,
                               unsigned short val_nVar,
                               const CConfig* config) :
                CUpwScalar(val_nDim, val_nVar, config) { }

void CUpwSca_TurbSA::ExtraADPreaccIn() {
  AD::SetPreaccIn(V_i, nDim+1);
  AD::SetPreaccIn(V_j, nDim+1);
}

void CUpwSca_TurbSA::FinishResidualCalc(const CConfig* config) {

  Flux[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];

  if (implicit) {
    Jacobian_i[0][0] = a0;
    Jacobian_j[0][0] = a1;
  }
}

CUpwSca_TurbSST::CUpwSca_TurbSST(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 const CConfig* config) :
                 CUpwScalar(val_nDim, val_nVar, config) { }

void CUpwSca_TurbSST::ExtraADPreaccIn() {
  AD::SetPreaccIn(V_i, nDim+3);
  AD::SetPreaccIn(V_j, nDim+3);
}

void CUpwSca_TurbSST::FinishResidualCalc(const CConfig* config) {

  Flux[0] = a0*Density_i*TurbVar_i[0]+a1*Density_j*TurbVar_j[0];
  Flux[1] = a0*Density_i*TurbVar_i[1]+a1*Density_j*TurbVar_j[1];

  if (implicit) {
    Jacobian_i[0][0] = a0;    Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;   Jacobian_i[1][1] = a0;

    Jacobian_j[0][0] = a1;    Jacobian_j[0][1] = 0.0;
    Jacobian_j[1][0] = 0.0;   Jacobian_j[1][1] = a1;
  }
}
