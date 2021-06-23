/*!
 * \file roe.cpp
 * \brief Implementations of Roe-type schemes in NEMO.
 * \author S. R. Copeland, W. Maier, C. Garbacz
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/NEMO/convection/roe.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwRoe_NEMO::CUpwRoe_NEMO(unsigned short val_nDim, unsigned short val_nVar,
                           unsigned short val_nPrimVar,
                           unsigned short val_nPrimVarGrad,
                           CConfig *config) : CNEMONumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad,
                                                          config) {

  unsigned short iVar;

  /*--- Allocate arrays ---*/
  Diff_U      = new su2double  [nVar];
  RoeU        = new su2double  [nVar];
  RoeV        = new su2double  [nPrimVar];
  RoedPdU     = new su2double  [nVar];
  Lambda      = new su2double  [nVar];
  Epsilon     = new su2double  [nVar];
  P_Tensor    = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];

  roe_eves.resize(nSpecies,0.0);

  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar]    = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }

  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];

  Flux   = new su2double[nVar];

}

CUpwRoe_NEMO::~CUpwRoe_NEMO(void) {

  unsigned short iVar;

  delete [] Diff_U;
  delete [] RoeU;
  delete [] RoeV;
  delete [] RoedPdU;
  delete [] Lambda;
  delete [] Epsilon;

  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }

  delete [] P_Tensor;
  delete [] invP_Tensor;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Flux;
}

CNumerics::ResidualType<> CUpwRoe_NEMO::ComputeResidual(const CConfig *config) {

  unsigned short iDim, iSpecies, iVar, jVar, kVar;

  /*--- Face area (norm or the normal vector) ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  /*--- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Calculate Roe averaged variables ---*/
  R = sqrt(abs(V_j[RHO_INDEX]/V_i[RHO_INDEX]));

  for (iVar = 0; iVar < nVar; iVar++)
    RoeU[iVar] = (R*U_j[iVar] + U_i[iVar])/(R+1);

  for (iVar = 0; iVar < nPrimVar; iVar++)
    RoeV[iVar] = (R*V_j[iVar] + V_i[iVar])/(R+1);

  auto& roe_eves = fluidmodel->ComputeSpeciesEve(RoeV[TVE_INDEX]);

  /*--- Calculate derivatives of pressure ---*/
  fluidmodel->ComputedPdU(RoeV, roe_eves, RoedPdU);

  /*--- Calculate dual grid tangent vectors for P & invP ---*/
  CreateBasis(UnitNormal);

  /*--- Compute the inviscid projected fluxes ---*/
  GetInviscidProjFlux(U_i, V_i, Normal, ProjFlux_i);
  GetInviscidProjFlux(U_j, V_j, Normal, ProjFlux_j);

  /*--- Compute projected P, invP, and Lambda ---*/
  GetPMatrix    (RoeU, RoeV, RoedPdU, UnitNormal, l, m, P_Tensor);
  GetPMatrix_inv(RoeU, RoeV, RoedPdU, UnitNormal, l, m, invP_Tensor);

  /*--- Compute projected velocities ---*/
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeV[VEL_INDEX+iDim] * UnitNormal[iDim];
    ProjVelocity_i += V_i[VEL_INDEX+iDim]  * UnitNormal[iDim];
    ProjVelocity_j += V_j[VEL_INDEX+iDim]  * UnitNormal[iDim];
  }

  RoeSoundSpeed = sqrt((1.0+RoedPdU[nSpecies+nDim])*
                            RoeV[P_INDEX]/RoeV[RHO_INDEX]);

  /*--- Calculate eigenvalues ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Lambda[iSpecies] = ProjVelocity;

  for (iDim = 0; iDim < nDim-1; iDim++)
    Lambda[nSpecies+iDim] = ProjVelocity;

  Lambda[nSpecies+nDim-1] = ProjVelocity + RoeSoundSpeed;
  Lambda[nSpecies+nDim]   = ProjVelocity - RoeSoundSpeed;
  Lambda[nSpecies+nDim+1] = ProjVelocity;

  /*--- Harten and Hyman (1983) entropy correction ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Epsilon[iSpecies] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,
                                         ProjVelocity_j-Lambda[iDim] ));
  for (iDim = 0; iDim < nDim-1; iDim++)
    Epsilon[nSpecies+iDim] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,
                                              ProjVelocity_j-Lambda[iDim] ));
  Epsilon[nSpecies+nDim-1] = 4.0*max(0.0, max(Lambda[nSpecies+nDim-1]-(ProjVelocity_i+V_i[A_INDEX]),
                                     (ProjVelocity_j+V_j[A_INDEX])-Lambda[nSpecies+nDim-1]));
  Epsilon[nSpecies+nDim]   = 4.0*max(0.0, max(Lambda[nSpecies+nDim]-(ProjVelocity_i-V_i[A_INDEX]),
                                     (ProjVelocity_j-V_j[A_INDEX])-Lambda[nSpecies+nDim]));
  Epsilon[nSpecies+nDim+1] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,
                                              ProjVelocity_j-Lambda[iDim] ));
  for (iVar = 0; iVar < nVar; iVar++)
    if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
      Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
    else
      Lambda[iVar] = fabs(Lambda[iVar]);

  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = fabs(Lambda[iVar]);

  /*--- Calculate inviscid projected Jacobians ---*/
  // Note: Scaling value is 0.5 because inviscid flux is based on 0.5*(Fc_i+Fc_j)
  //if (implicit){
  //  GetInviscidProjJac(U_i, V_i, dPdU_i, Normal, 0.5, val_Jacobian_i);
  //  GetInviscidProjJac(U_j, V_j, dPdU_j, Normal, 0.5, val_Jacobian_j);
  //}

  /*--- Difference of conserved variables at iPoint and jPoint ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = U_j[iVar]-U_i[iVar];

  /*--- Roe's Flux approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Flux[iVar] = 0.5 * (ProjFlux_i[iVar] + ProjFlux_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {

      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];

      Flux[iVar] -= 0.5*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area;
      //if (implicit){
      //  val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
      //  val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      //}
    }
  }

  return ResidualType<>(Flux, nullptr, nullptr);
}
