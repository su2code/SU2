/*!
 * \file lax.cpp
 * \brief Implementations of Lax centered scheme.
 * \author F. Palacios, S.R. Copeland, W. Maier, C. Garbacz
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/NEMO/convection/lax.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CCentLax_NEMO::CCentLax_NEMO(unsigned short val_nDim,
                             unsigned short val_nVar,
                             unsigned short val_nPrimVar,
                             unsigned short val_nPrimVarGrad,
                             CConfig *config) : CNEMONumerics(val_nDim, val_nVar,
                                                              val_nPrimVar,
                                                              val_nPrimVarGrad,
                                                              config) {

  /*--- Artifical dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_Flow();

  /*--- Allocate some structures ---*/
  Diff_U   = new su2double[nVar];
  MeanU    = new su2double[nVar];
  MeanV    = new su2double[nPrimVar];
  MeandPdU = new su2double[nVar];
  ProjFlux = new su2double[nVar];
  Flux     = new su2double[nVar];

  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CCentLax_NEMO::~CCentLax_NEMO() {

  delete [] Diff_U;
  delete [] MeanU;
  delete [] MeanV;
  delete [] MeandPdU;
  delete [] ProjFlux;
  delete [] Flux;
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Jacobian_i;
  delete [] Jacobian_j;
}

CNumerics::ResidualType<> CCentLax_NEMO::ComputeResidual(const CConfig *config) {

  /*--- Calculate geometrical quantities ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  for (auto iDim = 0ul; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Rename for convenience ---*/
  Density_i = V_i[RHO_INDEX]; Density_j = V_j[RHO_INDEX];
  Enthalpy_i = V_i[H_INDEX]; Enthalpy_j = V_j[H_INDEX];
  SoundSpeed_i = V_i[A_INDEX]; SoundSpeed_j = V_j[A_INDEX];

  /*--- Compute mean quantities for the variables ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++)
    MeanU[iVar] = 0.5*(U_i[iVar]+U_j[iVar]);
  for (auto iVar = 0ul; iVar < nPrimVar; iVar++)
    MeanV[iVar] = 0.5*(V_i[iVar]+V_j[iVar]);

  /*--- Compute NonEq specific variables ---*/
  const auto& mean_eves = fluidmodel->ComputeSpeciesEve(MeanV[TVE_INDEX]);
  fluidmodel->ComputedPdU(MeanV, mean_eves, MeandPdU);

  /*--- Get projected flux tensor ---*/
  GetInviscidProjFlux(MeanU, MeanV, Normal, ProjFlux);

  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_res ~ 0.5*(fc_i+fc_j)*Normal ---*/
  if (implicit) {
    GetInviscidProjJac(MeanU, MeanV, MeandPdU, Normal, 0.5, Jacobian_i);

    for (auto iVar = 0ul; iVar < nVar; iVar++)
      for (auto jVar = 0ul; jVar < nVar; jVar++)
        Jacobian_j[iVar][jVar] = Jacobian_i[iVar][jVar];
  }

  /*--- Compute the local spectral radius and the stretching factor ---*/
  ProjVelocity_i = GeometryToolbox::DotProduct(nDim, &V_i[VEL_INDEX],Normal);
  ProjVelocity_j = GeometryToolbox::DotProduct(nDim, &V_j[VEL_INDEX],Normal);

  /*--- Dissipation --*/
  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

  Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);

  /*--- Computes differences btw. conservative variables ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++)
    Diff_U[iVar] = U_i[iVar] - U_j[iVar];
  Diff_U[nSpecies+nDim] = Density_i*Enthalpy_i - Density_j*Enthalpy_j;

  /*--- Compute dissipation coefficient ---*/
  sc0 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc0*su2double(nDim)/3.0;

  /*--- Compute viscous part of the residual ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    Flux[iVar] = ProjFlux[iVar]+Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;
  }

  if (implicit) {

    cte = Epsilon_0*StretchingFactor*MeanLambda;

    for (auto iVar = 0u; iVar < nSpecies+nDim; iVar++) {
      Jacobian_i[iVar][iVar] += cte;
      Jacobian_j[iVar][iVar] -= cte;
    }

    /*--- Last rows: CAREFUL!! You have differences of \rho_Enthalpy, not differences of \rho_Energy ---*/
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      Jacobian_i[nSpecies+nDim][iSpecies] += cte*dPdU_i[iSpecies];
    for (auto iDim = 0ul; iDim < nDim; iDim++)
      Jacobian_i[nSpecies+nDim][nSpecies+iDim]   += cte*dPdU_i[nSpecies+iDim];
    Jacobian_i[nSpecies+nDim][nSpecies+nDim]     += cte*(1+dPdU_i[nSpecies+nDim]);
    Jacobian_i[nSpecies+nDim][nSpecies+nDim+1]   += cte*dPdU_i[nSpecies+nDim+1];
    Jacobian_i[nSpecies+nDim+1][nSpecies+nDim+1] += cte;

    /*--- Last row of Jacobian_j ---*/
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      Jacobian_j[nSpecies+nDim][iSpecies] -= cte*dPdU_j[iSpecies];
    for (auto iDim = 0ul; iDim < nDim; iDim++)
      Jacobian_j[nSpecies+nDim][nSpecies+iDim]   -= cte*dPdU_j[nSpecies+nDim];
    Jacobian_j[nSpecies+nDim][nSpecies+nDim]     -= cte*(1+dPdU_j[nSpecies+nDim]);
    Jacobian_j[nSpecies+nDim][nSpecies+nDim+1]   -= cte*dPdU_j[nSpecies+nDim+1];
    Jacobian_j[nSpecies+nDim+1][nSpecies+nDim+1] -= cte;
  }
  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}
