/*!
 * \file fvs.cpp
 * \brief Implementations of Flux-Vector-Splitting schemes.
 * \author F. Palacios, T. Economon
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

#include "../../../../include/numerics/flow/convection/fvs.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwMSW_Flow::CUpwMSW_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) : CNumerics(val_nDim, val_nVar, config) {

  if (config->GetDynamic_Grid() && (SU2_MPI::GetRank() == MASTER_NODE))
    cout << "WARNING: Grid velocities are NOT yet considered in the MSW scheme." << endl;

  /*--- Set booleans from CConfig settings ---*/
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Allocate arrays ---*/
  Fc_i     = new su2double [nVar];
  Fc_j     = new su2double [nVar];
  Lambda_i = new su2double [nVar];
  Lambda_j = new su2double [nVar];

  u_i      = new su2double [nDim];
  u_j      = new su2double [nDim];
  ust_i    = new su2double [nDim];
  ust_j    = new su2double [nDim];
  Vst_i    = new su2double [nDim+5];
  Vst_j    = new su2double [nDim+5];

  Velst_i    = new su2double [nDim];
  Velst_j    = new su2double [nDim];

  P_Tensor   = new su2double* [nVar];
  invP_Tensor= new su2double* [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }

}

CUpwMSW_Flow::~CUpwMSW_Flow() {

  delete [] Fc_i;
  delete [] Fc_j;
  delete [] Lambda_i;
  delete [] Lambda_j;

  delete [] u_i;
  delete [] u_j;
  delete [] ust_i;
  delete [] ust_j;
  delete [] Vst_i;
  delete [] Vst_j;
  delete [] Velst_i;
  delete [] Velst_j;

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  delete [] Jacobian_i;
  delete [] Jacobian_j;

}

CNumerics::ResidualType<> CUpwMSW_Flow::ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  unsigned short iDim, iVar, jVar, kVar;
  su2double rho_i, rho_j, P_i, P_j, H_i, H_j;
  su2double ProjVel_i, ProjVel_j, ProjVelst_i, ProjVelst_j;
  su2double sqvel_i, sqvel_j;
  su2double alpha, w, dp, onemw;
  su2double Proj_ModJac_Tensor_i, Proj_ModJac_Tensor_j;

  /*--- Set parameters in the numerical method ---*/
  alpha = 6.0;

  /*--- Calculate supporting geometry parameters ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Initialize flux & Jacobian vectors ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Fc_i[iVar] = 0.0;
    Fc_j[iVar] = 0.0;
  }
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] = 0.0;
        Jacobian_j[iVar][jVar] = 0.0;
      }
    }
  }

  /*--- Load variables from nodes i & j ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    u_i[iDim] = V_i[iDim+1];
    u_j[iDim] = V_j[iDim+1];
  }
  P_i = V_i[nDim+1];
  P_j = V_j[nDim+1];
  rho_i = V_i[nDim+2];
  rho_j = V_j[nDim+2];
  H_i = V_i[nDim+3];
  H_j = V_j[nDim+3];

  /*--- Recompute conservatives ---*/

  su2double U_i[5] = {0.0}, U_j[5] = {0.0};

  U_i[0] = rho_i; U_j[0] = rho_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = rho_i*u_i[iDim];
    U_j[iDim+1] = rho_j*u_j[iDim];
  }
  U_i[nDim+1] = rho_i*H_i - P_i;
  U_j[nDim+1] = rho_j*H_j - P_j;

  /*--- Calculate supporting quantities ---*/

  sqvel_i   = 0.0; sqvel_j   = 0.0;
  ProjVel_i = 0.0; ProjVel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    sqvel_i   += u_i[iDim]*u_i[iDim];
    sqvel_j   += u_j[iDim]*u_j[iDim];
    ProjVel_i += u_i[iDim]*UnitNormal[iDim];
    ProjVel_j += u_j[iDim]*UnitNormal[iDim];
  }

  /*--- Calculate the state weighting function ---*/

  dp = fabs(P_j-P_i) / min(P_j, P_i);
  w = 0.5 * (1.0/(pow(alpha*dp,2.0) +1.0));
  onemw = 1.0 - w;

  /*--- Calculate weighted state vector (*) for i & j ---*/

  for (iVar = 0; iVar < nDim+5; iVar++) {
    Vst_i[iVar] = onemw*V_i[iVar] + w*V_j[iVar];
    Vst_j[iVar] = onemw*V_j[iVar] + w*V_i[iVar];
  }
  ProjVelst_i = onemw*ProjVel_i + w*ProjVel_j;
  ProjVelst_j = onemw*ProjVel_j + w*ProjVel_i;

  for (iDim = 0; iDim < nDim; iDim++) {
    Velst_i[iDim] = Vst_i[iDim+1];
    Velst_j[iDim] = Vst_j[iDim+1];
  }

  /*--- Flow eigenvalues at i (Lambda+) ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
  Lambda_i[iDim] = 0.5*(ProjVelst_i + fabs(ProjVelst_i));
  }

  Lambda_i[nDim]   = 0.5*( ProjVelst_i + Vst_i[nDim+4] + fabs(ProjVelst_i + Vst_i[nDim+4]) );
  Lambda_i[nDim+1] = 0.5*( ProjVelst_i - Vst_i[nDim+4] + fabs(ProjVelst_i - Vst_i[nDim+4]) );

  /*--- Compute projected P, invP, and Lambda ---*/

  GetPMatrix(&Vst_i[nDim+2], Velst_i, &Vst_i[nDim+4], UnitNormal, P_Tensor);
  GetPMatrix_inv(&Vst_i[nDim+2], Velst_i, &Vst_i[nDim+4], UnitNormal, invP_Tensor);

  /*--- Projected flux (f+) at i ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_i = 0.0;

      /*--- Compute Proj_ModJac_Tensor = P x Lambda+ x inverse P ---*/

      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_i += P_Tensor[iVar][kVar]*Lambda_i[kVar]*invP_Tensor[kVar][jVar];
      Fc_i[iVar] += Proj_ModJac_Tensor_i*U_i[jVar]*Area;
      if (implicit)
        Jacobian_i[iVar][jVar] += Proj_ModJac_Tensor_i*Area;
    }
  }

  /*--- Flow eigenvalues at j (Lambda-) ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    Lambda_j[iDim] = 0.5*(ProjVelst_j - fabs(ProjVelst_j));
  }
  Lambda_j[nDim]   = 0.5*( ProjVelst_j + Vst_j[nDim+4] - fabs(ProjVelst_j + Vst_j[nDim+4]) );
  Lambda_j[nDim+1] = 0.5*( ProjVelst_j - Vst_j[nDim+4] - fabs(ProjVelst_j - Vst_j[nDim+4]) );

  /*--- Compute projected P, invP, and Lambda ---*/

  GetPMatrix(&Vst_j[nDim+2], Velst_j, &Vst_j[nDim+4], UnitNormal, P_Tensor);
  GetPMatrix_inv(&Vst_j[nDim+2], Velst_j, &Vst_j[nDim+4], UnitNormal, invP_Tensor);

  /*--- Projected flux (f-) ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_j = 0.0;
      /*--- Compute Proj_ModJac_Tensor = P x Lambda- x inverse P ---*/
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_j += P_Tensor[iVar][kVar]*Lambda_j[kVar]*invP_Tensor[kVar][jVar];
      Fc_j[iVar] += Proj_ModJac_Tensor_j*U_j[jVar]*Area;
      if (implicit)
        Jacobian_j[iVar][jVar] += Proj_ModJac_Tensor_j*Area;
    }
  }

  /*--- Flux splitting, use the i flux as final output. ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Fc_i[iVar] += Fc_j[iVar];
  }

  return ResidualType<>(Fc_i, Jacobian_i, Jacobian_j);

}
