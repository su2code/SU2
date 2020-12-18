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
                       const CConfig* config) :
  CNumerics(val_nDim, val_nVar, config),
  incompressible(config->GetKind_Regime() == INCOMPRESSIBLE),
  dynamic_grid(config->GetDynamic_Grid())
{

  nPrimVarTot = nVar + (nDim+1);

  Flux = new su2double [nVar] ();
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (auto iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nPrimVarTot] ();
    Jacobian_j[iVar] = new su2double [nPrimVarTot] ();
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

  /*--- Primitive variables ---*/

  for (auto iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
  }

  Density_i = V_i[nDim+2];
  Density_j = V_j[nDim+2];

  a_ij = 0.0;
  if (dynamic_grid) {
    for (auto iDim = 0; iDim < nDim; iDim++) {
      su2double Velocity_i = V_i[iDim+1] - GridVel_i[iDim];
      su2double Velocity_j = V_j[iDim+1] - GridVel_j[iDim];
      a_ij += 0.5*(Velocity_i+Velocity_j)*Normal[iDim];
    }
  }
  else {
    for (auto iDim = 0; iDim < nDim; iDim++)
      a_ij += 0.5*(V_i[iDim+1]+V_j[iDim+1])*Normal[iDim];
  }

  a_i = 0.5*(a_ij+fabs(a_ij));
  a_j = 0.5*(a_ij-fabs(a_ij));

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

  Flux[0] = a_i*TurbVar_i[0]+a_j*TurbVar_j[0];

  Jacobian_i[0][0] = a_i;
  Jacobian_j[0][0] = a_j;
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

  Flux[0] = a_i*Density_i*TurbVar_i[0]+a_j*Density_j*TurbVar_j[0];
  Flux[1] = a_i*Density_i*TurbVar_i[1]+a_j*Density_j*TurbVar_j[1];

  Jacobian_i[0][0] = Jacobian_i[1][1] = a_i;
  Jacobian_j[0][0] = Jacobian_j[1][1] = a_j;

  /*--- Jacobian wrt flow variables ---*/
  for (auto iVar = 0; iVar < nVar; iVar++) {
    const su2double dai_daij = 0.5*(1.0+a_ij/fabs(a_ij));
    const su2double daj_daij = 0.5*(1.0-a_ij/fabs(a_ij));
    const su2double dF_daij = Density_i*TurbVar_i[iVar]*dai_daij+Density_j*TurbVar_j[iVar]*daj_daij;
    for (auto iDim = 0; iDim < nDim; iDim++) {
      Jacobian_i[iVar][iDim+2] = 0.5*dF_daij/Density_i*Normal[iDim];
      Jacobian_j[iVar][iDim+2] = 0.5*dF_daij/Density_j*Normal[iDim];
    }
    Jacobian_i[iVar][nDim+2] = -a_i*TurbVar_i[iVar];
    Jacobian_j[iVar][nDim+2] = -a_j*TurbVar_j[iVar];
  }
}
