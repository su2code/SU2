/*!
 * \file msw.cpp
 * \brief Implementations of the modified Steger-Warming scheme.
 * \author ADL Stanford, S.R. Copeland, W. Maier, C. Garbacz
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

#include "../../../../include/numerics/NEMO/convection/msw.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwMSW_NEMO::CUpwMSW_NEMO(unsigned short val_nDim, unsigned short val_nVar,
                           unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
                           CConfig *config) : CNEMONumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad,
                                                          config) {

  /*--- Allocate arrays ---*/
  Diff_U   = new su2double [nVar];
  Fc_i     = new su2double [nVar];
  Fc_j     = new su2double [nVar];
  Lambda_i = new su2double [nVar];
  Lambda_j = new su2double [nVar];

  rhos_i   = new su2double [nSpecies];
  rhos_j   = new su2double [nSpecies];
  rhosst_i = new su2double [nSpecies];
  rhosst_j = new su2double [nSpecies];
  ust_i    = new su2double [nDim];
  ust_j    = new su2double [nDim];
  Vst_i    = new su2double [nPrimVar];
  Vst_j    = new su2double [nPrimVar];
  Ust_i    = new su2double [nVar];
  Ust_j    = new su2double [nVar];
  dPdUst_i = new su2double [nVar];
  dPdUst_j = new su2double [nVar];

  P_Tensor    = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    P_Tensor[iVar]    = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }

  Flux   = new su2double[nVar];
  Flux   = new su2double[nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CUpwMSW_NEMO::~CUpwMSW_NEMO() {

  delete [] Diff_U;
  delete [] Fc_i;
  delete [] Fc_j;
  delete [] Lambda_i;
  delete [] Lambda_j;

  delete [] rhos_i;
  delete [] rhos_j;
  delete [] rhosst_i;
  delete [] rhosst_j;
  delete [] ust_i;
  delete [] ust_j;
  delete [] Ust_i;
  delete [] Vst_i;
  delete [] Ust_j;
  delete [] Vst_j;
  delete [] dPdUst_i;
  delete [] dPdUst_j;

  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  delete [] Flux;
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Jacobian_i;
  delete [] Jacobian_j;
}

CNumerics::ResidualType<> CUpwMSW_NEMO::ComputeResidual(const CConfig *config) {

  /*--- Set parameters in the numerical method ---*/
  const su2double alpha   = 5.0;
  const su2double epsilon = 0.0;

  /*--- Calculate supporting geometry parameters ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  for (auto iDim = 0ul; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Initialize flux & Jacobian vectors ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    Fc_i[iVar] = 0.0;
    Fc_j[iVar] = 0.0;
  }
  if (implicit) {
    for (auto iVar = 0ul; iVar < nVar; iVar++) {
      for (auto jVar = 0ul; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] = 0.0;
        Jacobian_j[iVar][jVar] = 0.0;
      }
    }
  }

  /*--- Load variables from nodes i & j ---*/
  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
    rhos_i[iSpecies] = V_i[RHOS_INDEX+iSpecies];
    rhos_j[iSpecies] = V_j[RHOS_INDEX+iSpecies];
  }
  for (auto iDim = 0ul; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[VEL_INDEX+iDim];
    Velocity_j[iDim] = V_j[VEL_INDEX+iDim];
  }
  Pressure_i = V_i[P_INDEX];
  Pressure_j = V_j[P_INDEX];

  /*--- Calculate velocity quantities ---*/
  ProjVelocity_i = GeometryToolbox::DotProduct(nDim, Velocity_i, UnitNormal);
  ProjVelocity_j = GeometryToolbox::DotProduct(nDim, Velocity_j, UnitNormal);

  /*--- Calculate the state weighting function ---*/
  const su2double dp = fabs(Pressure_j-Pressure_i) / min(Pressure_j, Pressure_i);
  const su2double w = 0.5 * (1.0/(pow(alpha*dp,2.0) +1.0));
  const su2double onemw = 1.0 - w;

  /*--- Calculate weighted state vector (*) for i & j ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    Ust_i[iVar] = onemw*U_i[iVar] + w*U_j[iVar];
    Ust_j[iVar] = onemw*U_j[iVar] + w*U_i[iVar];
  }
  for (auto iVar = 0ul; iVar < nPrimVar; iVar++) {
    Vst_i[iVar] = onemw*V_i[iVar] + w*V_j[iVar];
    Vst_j[iVar] = onemw*V_j[iVar] + w*V_i[iVar];
  }
  const su2double ProjVelst_i = onemw*ProjVelocity_i + w*ProjVelocity_j;
  const su2double ProjVelst_j = onemw*ProjVelocity_j + w*ProjVelocity_i;

  const auto& eves_st_i = fluidmodel->ComputeSpeciesEve(Vst_i[TVE_INDEX]);
  const auto& eves_st_j = fluidmodel->ComputeSpeciesEve(Vst_j[TVE_INDEX]);

  fluidmodel->ComputedPdU(Vst_i, eves_st_i, dPdUst_i);
  fluidmodel->ComputedPdU(Vst_j, eves_st_j, dPdUst_j);

  /*--- Flow eigenvalues at i (Lambda+) ---*/
  for (auto iVar = 0; iVar < nSpecies+nDim-1; iVar++)
    Lambda_i[iVar]          = 0.5*(ProjVelst_i + sqrt(ProjVelst_i*ProjVelst_i +
                                                      epsilon*epsilon));
  Lambda_i[nSpecies+nDim-1] = 0.5*(ProjVelst_i + Vst_i[A_INDEX] +
                             sqrt((ProjVelst_i + Vst_i[A_INDEX])*
                                  (ProjVelst_i + Vst_i[A_INDEX])+
                                               epsilon*epsilon));
  Lambda_i[nSpecies+nDim]   = 0.5*(ProjVelst_i - Vst_i[A_INDEX] +
                             sqrt((ProjVelst_i - Vst_i[A_INDEX])*
                                  (ProjVelst_i - Vst_i[A_INDEX])+
                                               epsilon*epsilon));
  Lambda_i[nSpecies+nDim+1] = 0.5*(ProjVelst_i + sqrt(ProjVelst_i*ProjVelst_i +
                                                      epsilon*epsilon));

  /*--- Compute projected P, invP, and Lambda ---*/
  su2double l[MAXNDIM], m[MAXNDIM];
  CreateBasis(UnitNormal,l,m);
  GetPMatrix(Ust_i, Vst_i, dPdUst_i, UnitNormal, l, m, P_Tensor);
  GetPMatrix_inv(Ust_i, Vst_i, dPdUst_i, UnitNormal, l, m, invP_Tensor);

  /*--- Projected flux (f+) at i ---*/
  su2double Proj_ModJac_Tensor_i;
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    for (auto jVar = 0ul; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_i = 0.0;

      /*--- Compute Proj_ModJac_Tensor = P x Lambda+ x inverse P ---*/
      for (auto kVar = 0ul; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_i += P_Tensor[iVar][kVar]*Lambda_i[kVar]*invP_Tensor[kVar][jVar];
      Fc_i[iVar] += Proj_ModJac_Tensor_i*U_i[jVar]*Area;
      if (implicit)
        Jacobian_i[iVar][jVar] += Proj_ModJac_Tensor_i*Area;
    }
  }

  /*--- Flow eigenvalues at j (Lambda-) ---*/
  for (auto iVar = 0; iVar < nSpecies+nDim-1; iVar++)
    Lambda_j[iVar]          = 0.5*(ProjVelst_j - sqrt(ProjVelst_j*ProjVelst_j +
                                                      epsilon*epsilon));
  Lambda_j[nSpecies+nDim-1] = 0.5*(ProjVelst_j + Vst_j[A_INDEX] -
                             sqrt((ProjVelst_j + Vst_j[A_INDEX])*
                                  (ProjVelst_j + Vst_j[A_INDEX])+
                                               epsilon*epsilon));
  Lambda_j[nSpecies+nDim]   = 0.5*(ProjVelst_j - Vst_j[A_INDEX] -
                             sqrt((ProjVelst_j - Vst_j[A_INDEX])*
                                  (ProjVelst_j - Vst_j[A_INDEX])+
                                                epsilon*epsilon)                 );
  Lambda_j[nSpecies+nDim+1] = 0.5*(ProjVelst_j - sqrt(ProjVelst_j*ProjVelst_j+
                                                      epsilon*epsilon));

  /*--- Compute projected P, invP, and Lambda ---*/
  CreateBasis(UnitNormal,l,m);
  GetPMatrix(Ust_j, Vst_j, dPdUst_j, UnitNormal, l, m, P_Tensor);
  GetPMatrix_inv(Ust_j, Vst_j, dPdUst_j, UnitNormal, l, m, invP_Tensor);

  /*--- Projected flux (f-) ---*/
  su2double Proj_ModJac_Tensor_j;
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    for (auto jVar = 0ul; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_j = 0.0;

      /*--- Compute Proj_ModJac_Tensor = P x Lambda- x inverse P ---*/
      for (auto kVar = 0ul; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_j += P_Tensor[iVar][kVar]*Lambda_j[kVar]*invP_Tensor[kVar][jVar];
      Fc_j[iVar] += Proj_ModJac_Tensor_j*U_j[jVar]*Area;
      if (implicit)
        Jacobian_j[iVar][jVar] += Proj_ModJac_Tensor_j*Area;
    }
  }

  /*--- Flux splitting ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    Flux[iVar] = Fc_i[iVar]+Fc_j[iVar];
  }

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}
