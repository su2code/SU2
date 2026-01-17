/*!
 * \file fvs.cpp
 * \brief Implementations of Flux-Vector-Splitting schemes.
 * \author F. Palacios, T. Economon
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/flow/convection/fvs.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwMSW_Flow::CUpwMSW_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
  CNumerics(val_nDim, val_nVar, config), alpha(config->GetMSW_Alpha()) {

  if (config->GetDynamic_Grid() && SU2_MPI::GetRank() == MASTER_NODE) {
    std::cout << "WARNING: Grid velocities are NOT yet considered in the MSW scheme." << std::endl;
  }

  /*--- Allocate arrays ---*/
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = &buf_Jacobian_i[iVar * MAXNVAR];
    Jacobian_j[iVar] = &buf_Jacobian_j[iVar * MAXNVAR];
  }
}

CNumerics::ResidualType<> CUpwMSW_Flow::ComputeResidual(const CConfig* config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim + 4);
  AD::SetPreaccIn(V_j, nDim + 4);
  AD::SetPreaccIn(Normal, nDim);

  /*--- Calculate supporting geometry parameters ---*/

  const su2double Area = GeometryToolbox::Norm(nDim, Normal);

  for (auto iDim = 0u; iDim < nDim; ++iDim) {
    UnitNormal[iDim] = Normal[iDim] / Area;
  }

  /*--- Load variables from nodes i & j ---*/

  const su2double P_i = V_i[nDim + 1];
  const su2double P_j = V_j[nDim + 1];
  const su2double rho_i = V_i[nDim + 2];
  const su2double rho_j = V_j[nDim + 2];
  const su2double H_i = V_i[nDim + 3];
  const su2double H_j = V_j[nDim + 3];
  /*--- Recompute the speed of sound because it is not MUSCL-reconstructed. ---*/
  const su2double sqvel_i = GeometryToolbox::SquaredNorm(nDim, V_i + 1);
  const su2double sqvel_j = GeometryToolbox::SquaredNorm(nDim, V_j + 1);
  const su2double c_i = sqrt(fmax((Gamma - 1) * (H_i - 0.5 * sqvel_i), EPS));
  const su2double c_j = sqrt(fmax((Gamma - 1) * (H_j - 0.5 * sqvel_j), EPS));

  /*--- Recompute conservatives ---*/

  su2double U_i[MAXNVAR] = {}, U_j[MAXNVAR] = {};
  U_i[0] = rho_i;
  U_j[0] = rho_j;
  for (auto iDim = 0u; iDim < nDim; ++iDim) {
    U_i[iDim + 1] = rho_i * V_i[iDim + 1];
    U_j[iDim + 1] = rho_j * V_j[iDim + 1];
  }
  U_i[nDim + 1] = rho_i * H_i - P_i;
  U_j[nDim + 1] = rho_j * H_j - P_j;

  /*--- Calculate the state weighting function ---*/

  const su2double dp = fabs(P_j - P_i) / fmin(P_j, P_i);
  const su2double w = 0.5 * (1 / (pow(alpha * dp, 2) + 1));
  const su2double onemw = 1 - w;

  /*--- Calculate weighted state vector (*) for i & j ---*/

  su2double Vst_i[MAXNDIM + 5] = {}, Vst_j[MAXNDIM + 5] = {};
  for (auto iVar = 0; iVar < nDim + 4; ++iVar) {
    Vst_i[iVar] = onemw * V_i[iVar] + w * V_j[iVar];
    Vst_j[iVar] = onemw * V_j[iVar] + w * V_i[iVar];
  }
  Vst_i[nDim + 4] = onemw * c_i + w * c_j;
  Vst_j[nDim + 4] = onemw * c_j + w * c_i;

  su2double Velst_i[MAXNDIM] = {}, Velst_j[MAXNDIM] = {};
  for (auto iDim = 0u; iDim < nDim; ++iDim) {
    Velst_i[iDim] = Vst_i[iDim + 1];
    Velst_j[iDim] = Vst_j[iDim + 1];
  }
  const su2double ProjVelst_i = GeometryToolbox::DotProduct(MAXNDIM, Velst_i, UnitNormal);
  const su2double ProjVelst_j = GeometryToolbox::DotProduct(MAXNDIM, Velst_j, UnitNormal);

  /*--- Flow eigenvalues at i (Lambda+) ---*/

  su2double Lambda[MAXNVAR] = {};
  for (auto iDim = 0u; iDim < nDim; ++iDim) {
    Lambda[iDim] = fmax(ProjVelst_i, 0);
  }
  Lambda[nDim] = fmax(ProjVelst_i + Vst_i[nDim + 4], 0);
  Lambda[nDim + 1] = fmax(ProjVelst_i - Vst_i[nDim + 4], 0);

  /*--- Compute projected P, invP, and Lambda ---*/

  su2double P_Tensor[MAXNVAR][MAXNVAR], invP_Tensor[MAXNVAR][MAXNVAR];
  GetPMatrix(Vst_i[nDim + 2], Velst_i, Vst_i[nDim + 4], UnitNormal, P_Tensor);
  GetPMatrix_inv(Vst_i[nDim + 2], Velst_i, Vst_i[nDim + 4], UnitNormal, invP_Tensor);

  /*--- Projected flux (f+) at i ---*/

  for (auto iVar = 0u; iVar < MAXNVAR; ++iVar) {
    Fc[iVar] = 0;
  }
  auto UpdateFlux = [&](const auto* U, auto* Jacobian) {
    for (auto iVar = 0u; iVar < nVar; ++iVar) {
      for (auto jVar = 0u; jVar < nVar; ++jVar) {
        su2double Proj_ModJac_Tensor_i = 0.0;

        /*--- Compute Proj_ModJac_Tensor = P x Lambda+ x inverse P ---*/

        for (auto kVar = 0u; kVar < nVar; ++kVar) {
          Proj_ModJac_Tensor_i += P_Tensor[iVar][kVar] * Lambda[kVar] * invP_Tensor[kVar][jVar];
        }
        Fc[iVar] += Proj_ModJac_Tensor_i * Area * U[jVar];
        Jacobian[iVar][jVar] = Proj_ModJac_Tensor_i * Area;
      }
    }
  };
  UpdateFlux(U_i, Jacobian_i);

  /*--- Flow eigenvalues at j (Lambda-) ---*/

  for (auto iDim = 0u; iDim < nDim; ++iDim) {
    Lambda[iDim] = fmin(ProjVelst_j, 0);
  }
  Lambda[nDim] = fmin(ProjVelst_j + Vst_j[nDim + 4], 0);
  Lambda[nDim + 1] = fmin(ProjVelst_j - Vst_j[nDim + 4], 0);

  /*--- Compute projected P, invP, and Lambda ---*/

  GetPMatrix(Vst_j[nDim + 2], Velst_j, Vst_j[nDim + 4], UnitNormal, P_Tensor);
  GetPMatrix_inv(Vst_j[nDim + 2], Velst_j, Vst_j[nDim + 4], UnitNormal, invP_Tensor);

  /*--- Projected flux (f-) ---*/

  UpdateFlux(U_j, Jacobian_j);

  AD::SetPreaccOut(Fc, nVar);
  AD::EndPreacc();

  return ResidualType<>(Fc, Jacobian_i, Jacobian_j);

}
