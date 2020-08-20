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

void CUpwScalar::GetMUSCLJac(su2double **jac_i, su2double **jac_j,
                             const su2double *lim_i, const su2double *lim_j,
                             const su2double *r_i, const su2double *r_j,
                             const su2double *r_n_i, const su2double *r_n_j) {
  const bool wasActive = AD::BeginPassive();

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    const su2double dFidUi = jac_i[iVar][iVar]*(*r_i)*(1.0-muscl_kappa*lim_i[iVar])/(*r_n_i);
    const su2double dFidUj = jac_i[iVar][iVar]*(*r_i)*muscl_kappa*lim_i[iVar]/(*r_n_j);
    const su2double dFjdUi = jac_j[iVar][iVar]*(*r_j)*muscl_kappa*lim_j[iVar]/(*r_n_i);
    const su2double dFjdUj = jac_j[iVar][iVar]*(*r_j)*(1.0-muscl_kappa*lim_j[iVar])/(*r_n_j);
    
    jac_i[iVar][iVar] = dFidUi + dFjdUi;
    jac_j[iVar][iVar] = dFidUj + dFjdUj;
  }


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

  Lambda = 0.0;
  ProjVel_i   = 0.0;
  ProjVel_j   = 0.0;
  if (dynamic_grid) {
    for (iDim = 0; iDim < nDim; iDim++) {
      su2double Velocity_i = V_i[iDim+1] - GridVel_i[iDim];
      su2double Velocity_j = V_j[iDim+1] - GridVel_j[iDim];
      Lambda += 0.5*(Velocity_i+Velocity_j)*Normal[iDim];
    }
  }
  else {
    for (iDim = 0; iDim < nDim; iDim++) {
      const su2double RoeVelocity = (R*V_j[iDim+1]+V_i[iDim+1])/(R+1.);
      Lambda    += 0.5*RoeVelocity*Normal[iDim];
      ProjVel_i += 0.5*V_i[iDim+1]*Normal[iDim];
      ProjVel_j += 0.5*V_j[iDim+1]*Normal[iDim];
    }
  }

  /*--- Harten and Hyman (1983) entropy correction ---*/

  Epsilon = 4.0*max(0.0, max(Lambda-ProjVel_i, ProjVel_j-Lambda));
  Lambda  = (fabs(Lambda) < Epsilon) ? su2double(0.5*(Lambda*Lambda/Epsilon + Epsilon))
                                     : su2double(fabs(Lambda));

  FinishResidualCalc(config);

  if (muscl) {

    /*--- Extract nodal values ---*/

    const su2double Density_n_i = Vn_i[nDim+2];
    const su2double Density_n_j = Vn_j[nDim+2];

    /*--- Compute Jacobian wrt extrapolation ---*/

    GetMUSCLJac(Jacobian_i, Jacobian_j, Limiter_i, Limiter_j, 
                &Density_i, &Density_j, &Density_n_i, &Density_n_j);
  }
  
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

  Flux[0] = ProjVel_i*TurbVar_i[0]+ProjVel_j*TurbVar_j[0];

  Jacobian_i[0][0] = ProjVel_i;
  Jacobian_j[0][0] = ProjVel_j;
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

  Flux[0] = ProjVel_i*Density_i*TurbVar_i[0]+ProjVel_j*Density_j*TurbVar_j[0]-Lambda*(Density_j*TurbVar_j[0]-Density_i*TurbVar_i[0]);
  Flux[1] = ProjVel_i*Density_i*TurbVar_i[1]+ProjVel_j*Density_j*TurbVar_j[1]-Lambda*(Density_j*TurbVar_j[1]-Density_i*TurbVar_i[1]);

  Jacobian_i[0][0] = ProjVel_i+Lambda;  Jacobian_i[0][1] = 0.0;
  Jacobian_i[1][0] = 0.0; Jacobian_i[1][1] = ProjVel_i+Lambda;

  Jacobian_j[0][0] = ProjVel_j-Lambda;  Jacobian_j[0][1] = 0.0;
  Jacobian_j[1][0] = 0.0; Jacobian_j[1][1] = ProjVel_j-Lambda;
}
