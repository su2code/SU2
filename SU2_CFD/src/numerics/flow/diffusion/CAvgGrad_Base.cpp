/*!
 * \file CAvgGrad_Base.cpp
 * \brief Implementation of numerics class CAvgGrad_Base.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/flow/diffusion/CAvgGrad_Base.hpp"

CAvgGrad_Base::CAvgGrad_Base(unsigned short val_nDim,
                             unsigned short val_nVar,
                             unsigned short val_nPrimVar,
                             bool val_correct_grad,
                             CConfig *config)
    : CNumerics(val_nDim, val_nVar, config),
      nPrimVar(val_nPrimVar),
      correct_gradient(val_correct_grad) {

  unsigned short iVar, iDim;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  TauWall_i = 0; TauWall_j = 0;

  PrimVar_i = new su2double [nPrimVar];
  PrimVar_j = new su2double [nPrimVar];
  Mean_PrimVar = new su2double [nPrimVar];

  Mean_GradPrimVar = new su2double* [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++)
    Mean_GradPrimVar[iVar] = new su2double [nDim];

  Edge_Vector = new su2double[nDim];

  if (correct_gradient) {
    Proj_Mean_GradPrimVar_Edge = new su2double[val_nPrimVar];
  } else {
    Proj_Mean_GradPrimVar_Edge = NULL;
  }

  tau_jacobian_i = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    tau_jacobian_i[iDim] = new su2double [nVar];
  }
  heat_flux_vector = new su2double[nDim];
  heat_flux_jac_i = new su2double[nVar];

}

CAvgGrad_Base::~CAvgGrad_Base() {

  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Mean_PrimVar;
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;

  if (tau_jacobian_i != NULL) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      delete [] tau_jacobian_i[iDim];
    }
    delete [] tau_jacobian_i;
  }
  if (heat_flux_vector != NULL) {
    delete [] heat_flux_vector;
  }
  if (heat_flux_jac_i != NULL) {
    delete [] heat_flux_jac_i;
  }
  
  delete [] Edge_Vector;
  if (Proj_Mean_GradPrimVar_Edge != NULL)
    delete [] Proj_Mean_GradPrimVar_Edge;
}

void CAvgGrad_Base::CorrectGradient(su2double** GradPrimVar,
                                    const su2double* val_PrimVar_i,
                                    const su2double* val_PrimVar_j,
                                    const su2double* val_edge_vector,
                                    const su2double val_dist_ij_2,
                                    const unsigned short val_nPrimVar) {
  for (unsigned short iVar = 0; iVar < val_nPrimVar; iVar++) {
    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Proj_Mean_GradPrimVar_Edge[iVar] += GradPrimVar[iVar][iDim]*val_edge_vector[iDim];
    }
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                 (val_PrimVar_j[iVar]-val_PrimVar_i[iVar]))*val_edge_vector[iDim] / val_dist_ij_2;
    }
  }
}

void CAvgGrad_Base::SetStressTensor(const su2double *val_primvar,
                           const su2double* const *val_gradprimvar,
                           const su2double val_turb_ke,
                           const su2double val_laminar_viscosity,
                           const su2double val_eddy_viscosity) {

  unsigned short iDim, jDim;
  const su2double Density = val_primvar[nDim+2];
  const su2double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;

  su2double div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    div_vel += val_gradprimvar[iDim+1][iDim];

  /* --- If UQ methodology is used, calculate tau using the perturbed reynolds stress tensor --- */

  if (using_uq){
    for (iDim = 0 ; iDim < nDim; iDim++)
      for (jDim = 0 ; jDim < nDim; jDim++)
        tau[iDim][jDim] = val_laminar_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
        - TWO3*val_laminar_viscosity*div_vel*delta[iDim][jDim] - Density * MeanPerturbedRSM[iDim][jDim];

  } else {

    for (iDim = 0 ; iDim < nDim; iDim++)
      for (jDim = 0 ; jDim < nDim; jDim++)
        tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
                        - TWO3*total_viscosity*div_vel*delta[iDim][jDim];
  }
}

void CAvgGrad_Base::AddQCR(const su2double* const *val_gradprimvar) {

  su2double den_aux, c_cr1= 0.3, O_ik, O_jk;
  unsigned short iDim, jDim, kDim;

  /*--- Denominator Antisymmetric normalized rotation tensor ---*/

  den_aux = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    for (jDim = 0 ; jDim < nDim; jDim++)
      den_aux += val_gradprimvar[iDim+1][jDim] * val_gradprimvar[iDim+1][jDim];
  den_aux = sqrt(max(den_aux,1E-10));

  /*--- Adding the QCR contribution ---*/

  for (iDim = 0 ; iDim < nDim; iDim++){
    for (jDim = 0 ; jDim < nDim; jDim++){
      for (kDim = 0 ; kDim < nDim; kDim++){
        O_ik = (val_gradprimvar[iDim+1][kDim] - val_gradprimvar[kDim+1][iDim])/ den_aux;
        O_jk = (val_gradprimvar[jDim+1][kDim] - val_gradprimvar[kDim+1][jDim])/ den_aux;
        tau[iDim][jDim] -= c_cr1 * ((O_ik * tau[jDim][kDim]) + (O_jk * tau[iDim][kDim]));
      }
    }
  }
}

void CAvgGrad_Base::AddTauWall(const su2double *val_normal,
                               const su2double val_tau_wall) {

  unsigned short iDim, jDim;
  su2double TauNormal, TauElem[3], TauTangent[3], WallShearStress, Area, UnitNormal[3];

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += val_normal[iDim]*val_normal[iDim];
  Area = sqrt(Area);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = val_normal[iDim]/Area;

  /*--- First, compute wall shear stress as the magnitude of the wall-tangential
   component of the shear stress tensor---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    TauElem[iDim] = 0.0;
    for (jDim = 0; jDim < nDim; jDim++)
      TauElem[iDim] += tau[iDim][jDim]*UnitNormal[jDim];
  }

  TauNormal = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    TauNormal += TauElem[iDim] * UnitNormal[iDim];

  for (iDim = 0; iDim < nDim; iDim++)
    TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];

  WallShearStress = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    WallShearStress += TauTangent[iDim]*TauTangent[iDim];
  WallShearStress = sqrt(WallShearStress);

  /*--- Scale the stress tensor by the ratio of the wall shear stress
   to the computed representation of the shear stress ---*/

  for (iDim = 0 ; iDim < nDim; iDim++)
    for (jDim = 0 ; jDim < nDim; jDim++)
      tau[iDim][jDim] = tau[iDim][jDim]*(val_tau_wall/WallShearStress);
}

void CAvgGrad_Base::GetMeanRateOfStrainMatrix(su2double **S_ij) const
{
  /* --- Calculate the rate of strain tensor, using mean velocity gradients --- */

  if (nDim == 3){
    S_ij[0][0] = Mean_GradPrimVar[1][0];
    S_ij[1][1] = Mean_GradPrimVar[2][1];
    S_ij[2][2] = Mean_GradPrimVar[3][2];
    S_ij[0][1] = 0.5 * (Mean_GradPrimVar[1][1] + Mean_GradPrimVar[2][0]);
    S_ij[0][2] = 0.5 * (Mean_GradPrimVar[1][2] + Mean_GradPrimVar[3][0]);
    S_ij[1][2] = 0.5 * (Mean_GradPrimVar[2][2] + Mean_GradPrimVar[3][1]);
    S_ij[1][0] = S_ij[0][1];
    S_ij[2][1] = S_ij[1][2];
    S_ij[2][0] = S_ij[0][2];
  }
  else {
    S_ij[0][0] = Mean_GradPrimVar[1][0];
    S_ij[1][1] = Mean_GradPrimVar[2][1];
    S_ij[2][2] = 0.0;
    S_ij[0][1] = 0.5 * (Mean_GradPrimVar[1][1] + Mean_GradPrimVar[2][0]);
    S_ij[0][2] = 0.0;
    S_ij[1][2] = 0.0;
    S_ij[1][0] = S_ij[0][1];
    S_ij[2][1] = S_ij[1][2];
    S_ij[2][0] = S_ij[0][2];

  }
}

void CAvgGrad_Base::SetReynoldsStressMatrix(su2double turb_ke){
  unsigned short iDim, jDim;
  su2double **S_ij = new su2double* [3];
  su2double muT = Mean_Eddy_Viscosity;
  su2double divVel = 0;
  su2double density;
  su2double TWO3 = 2.0/3.0;
  density = Mean_PrimVar[nDim+2];

  for (iDim = 0; iDim < 3; iDim++){
    S_ij[iDim] = new su2double [3];
  }


  GetMeanRateOfStrainMatrix(S_ij);

  /* --- Using rate of strain matrix, calculate Reynolds stress tensor --- */

  for (iDim = 0; iDim < 3; iDim++){
    divVel += S_ij[iDim][iDim];
  }

  for (iDim = 0; iDim < 3; iDim++){
    for (jDim = 0; jDim < 3; jDim++){
      MeanReynoldsStress[iDim][jDim] = TWO3 * turb_ke * delta3[iDim][jDim]
      - muT / density * (2 * S_ij[iDim][jDim] - TWO3 * divVel * delta3[iDim][jDim]);
    }
  }

  for (iDim = 0; iDim < 3; iDim++)
    delete [] S_ij[iDim];
  delete [] S_ij;
}

void CAvgGrad_Base::SetPerturbedRSM(su2double turb_ke, CConfig *config){

  unsigned short iDim,jDim;

  /* --- Calculate anisotropic part of Reynolds Stress tensor --- */

  for (iDim = 0; iDim< 3; iDim++){
    for (jDim = 0; jDim < 3; jDim++){
      A_ij[iDim][jDim] = .5 * MeanReynoldsStress[iDim][jDim] / turb_ke - delta3[iDim][jDim] / 3.0;
      Eig_Vec[iDim][jDim] = A_ij[iDim][jDim];
    }
  }

  /* --- Get ordered eigenvectors and eigenvalues of A_ij --- */

  EigenDecomposition(A_ij, Eig_Vec, Eig_Val, 3);

  /* compute convex combination coefficients */
  su2double c1c = Eig_Val[2] - Eig_Val[1];
  su2double c2c = 2.0 * (Eig_Val[1] - Eig_Val[0]);
  su2double c3c = 3.0 * Eig_Val[0] + 1.0;

  /* define barycentric traingle corner points */
  Corners[0][0] = 1.0;
  Corners[0][1] = 0.0;
  Corners[1][0] = 0.0;
  Corners[1][1] = 0.0;
  Corners[2][0] = 0.5;
  Corners[2][1] = 0.866025;

  /* define barycentric coordinates */
  Barycentric_Coord[0] = Corners[0][0] * c1c + Corners[1][0] * c2c + Corners[2][0] * c3c;
  Barycentric_Coord[1] = Corners[0][1] * c1c + Corners[1][1] * c2c + Corners[2][1] * c3c;

  if (Eig_Val_Comp == 1) {
    /* 1C turbulence */
    New_Coord[0] = Corners[0][0];
    New_Coord[1] = Corners[0][1];
  }
  else if (Eig_Val_Comp== 2) {
    /* 2C turbulence */
    New_Coord[0] = Corners[1][0];
    New_Coord[1] = Corners[1][1];
  }
  else if (Eig_Val_Comp == 3) {
    /* 3C turbulence */
    New_Coord[0] = Corners[2][0];
    New_Coord[1] = Corners[2][1];
  }
  else {
    /* 2C turbulence */
    New_Coord[0] = Corners[1][0];
    New_Coord[1] = Corners[1][1];
  }

  /* calculate perturbed barycentric coordinates */
  Barycentric_Coord[0] = Barycentric_Coord[0] + (uq_delta_b) * (New_Coord[0] - Barycentric_Coord[0]);
  Barycentric_Coord[1] = Barycentric_Coord[1] + (uq_delta_b) * (New_Coord[1] - Barycentric_Coord[1]);

  /* rebuild c1c,c2c,c3c based on perturbed barycentric coordinates */
  c3c = Barycentric_Coord[1] / Corners[2][1];
  c1c = Barycentric_Coord[0] - Corners[2][0] * c3c;
  c2c = 1 - c1c - c3c;

  /* build new anisotropy eigenvalues */
  Eig_Val[0] = (c3c - 1) / 3.0;
  Eig_Val[1] = 0.5 *c2c + Eig_Val[0];
  Eig_Val[2] = c1c + Eig_Val[1];

  /* permute eigenvectors if required */
  if (uq_permute) {
    for (iDim=0; iDim<3; iDim++) {
      for (jDim=0; jDim<3; jDim++) {
        New_Eig_Vec[iDim][jDim] = Eig_Vec[2-iDim][jDim];
      }
    }
  }

  else {
    for (iDim=0; iDim<3; iDim++) {
      for (jDim=0; jDim<3; jDim++) {
        New_Eig_Vec[iDim][jDim] = Eig_Vec[iDim][jDim];
      }
    }
  }

  EigenRecomposition(newA_ij, New_Eig_Vec, Eig_Val, 3);

  /* compute perturbed Reynolds stress matrix; use under-relaxation factor (uq_urlx)*/
  for (iDim = 0; iDim< 3; iDim++){
    for (jDim = 0; jDim < 3; jDim++){
      MeanPerturbedRSM[iDim][jDim] = 2.0 * turb_ke * (newA_ij[iDim][jDim] + 1.0/3.0 * delta3[iDim][jDim]);
      MeanPerturbedRSM[iDim][jDim] = MeanReynoldsStress[iDim][jDim] +
      uq_urlx*(MeanPerturbedRSM[iDim][jDim] - MeanReynoldsStress[iDim][jDim]);
    }
  }

}


void CAvgGrad_Base::SetTauJacobian(const su2double *val_Mean_PrimVar,
                                   const su2double val_laminar_viscosity,
                                   const su2double val_eddy_viscosity,
                                   const su2double val_dist_ij,
                                   const su2double *val_normal) {

  /*--- QCR and wall functions are **not** accounted for here ---*/

  const su2double Density = val_Mean_PrimVar[nDim+2];
  const su2double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
  const su2double xi = total_viscosity/(Density*val_dist_ij);

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      // Jacobian w.r.t. momentum
      tau_jacobian_i[iDim][jDim+1] = -xi*(delta[iDim][jDim] + val_normal[iDim]*val_normal[jDim]/3.0);
    }
    // Jacobian w.r.t. density
    tau_jacobian_i[iDim][0] = 0;
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
       tau_jacobian_i[iDim][0] -= tau_jacobian_i[iDim][jDim+1]*val_Mean_PrimVar[jDim+1];
    }
    // Jacobian w.r.t. energy
    tau_jacobian_i[iDim][nDim+1] = 0;
  }
}

void CAvgGrad_Base::SetIncTauJacobian(const su2double val_laminar_viscosity,
                                      const su2double val_eddy_viscosity,
                                      const su2double val_dist_ij,
                                      const su2double *val_normal) {

  const su2double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
  const su2double xi = total_viscosity/val_dist_ij;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    tau_jacobian_i[iDim][0] = 0;
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      tau_jacobian_i[iDim][jDim+1] = -xi*(delta[iDim][jDim] + val_normal[iDim]*val_normal[jDim]/3.0);
    }
    tau_jacobian_i[iDim][nDim+1] = 0;
  }
}

void CAvgGrad_Base::GetViscousProjFlux(const su2double *val_primvar,
                                       const su2double *val_normal) {

  /*--- Primitive variables -> [Temp vel_x vel_y vel_z Pressure] ---*/

  if (nDim == 2) {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2]+
        heat_flux_vector[0];
    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2]+
        heat_flux_vector[1];
  } else {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][2];
    Flux_Tensor[4][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2] + tau[0][2]*val_primvar[3] +
        heat_flux_vector[0];
    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][2];
    Flux_Tensor[4][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2] + tau[1][2]*val_primvar[3] +
        heat_flux_vector[1];
    Flux_Tensor[0][2] = 0.0;
    Flux_Tensor[1][2] = tau[2][0];
    Flux_Tensor[2][2] = tau[2][1];
    Flux_Tensor[3][2] = tau[2][2];
    Flux_Tensor[4][2] = tau[2][0]*val_primvar[1] + tau[2][1]*val_primvar[2] + tau[2][2]*val_primvar[3] +
        heat_flux_vector[2];
  }
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Proj_Flux_Tensor[iVar] = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim] * val_normal[iDim];
  }
  
}

void CAvgGrad_Base::GetViscousProjJacs(const su2double *val_Mean_PrimVar,
                                       const su2double val_dS,
                                       const su2double *val_Proj_Visc_Flux,
                                       su2double **val_Proj_Jac_Tensor_i,
                                       su2double **val_Proj_Jac_Tensor_j) {

  const su2double Density = val_Mean_PrimVar[nDim+2];
  const su2double factor = 0.5/Density;

  if (nDim == 2) {

    val_Proj_Jac_Tensor_i[0][0] = 0.0;
    val_Proj_Jac_Tensor_i[0][1] = 0.0;
    val_Proj_Jac_Tensor_i[0][2] = 0.0;
    val_Proj_Jac_Tensor_i[0][3] = 0.0;
    val_Proj_Jac_Tensor_i[1][0] = val_dS*tau_jacobian_i[0][0];
    val_Proj_Jac_Tensor_i[1][1] = val_dS*tau_jacobian_i[0][1];
    val_Proj_Jac_Tensor_i[1][2] = val_dS*tau_jacobian_i[0][2];
    val_Proj_Jac_Tensor_i[1][3] = val_dS*tau_jacobian_i[0][3];
    val_Proj_Jac_Tensor_i[2][0] = val_dS*tau_jacobian_i[1][0];
    val_Proj_Jac_Tensor_i[2][1] = val_dS*tau_jacobian_i[1][1];
    val_Proj_Jac_Tensor_i[2][2] = val_dS*tau_jacobian_i[1][2];
    val_Proj_Jac_Tensor_i[2][3] = val_dS*tau_jacobian_i[1][3];
    const su2double contraction = tau_jacobian_i[0][0]*val_Mean_PrimVar[1] +
                                  tau_jacobian_i[1][0]*val_Mean_PrimVar[2];
    val_Proj_Jac_Tensor_i[3][0] = val_dS*(contraction - heat_flux_jac_i[0]);
    val_Proj_Jac_Tensor_i[3][1] = -val_dS*(tau_jacobian_i[0][0] + heat_flux_jac_i[1]); 
    val_Proj_Jac_Tensor_i[3][2] = -val_dS*(tau_jacobian_i[1][0] + heat_flux_jac_i[2]); 
    val_Proj_Jac_Tensor_i[3][3] = -val_dS*heat_flux_jac_i[3];
    
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      for (unsigned short jVar = 0; jVar < nVar; jVar++)
        val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];

    const su2double proj_viscousflux_vel= val_Proj_Visc_Flux[1]*val_Mean_PrimVar[1] + 
                                          val_Proj_Visc_Flux[2]*val_Mean_PrimVar[2];
    val_Proj_Jac_Tensor_i[3][0] -= factor*proj_viscousflux_vel;
    val_Proj_Jac_Tensor_j[3][0] -= factor*proj_viscousflux_vel;
    val_Proj_Jac_Tensor_i[3][1] += factor*val_Proj_Visc_Flux[1];
    val_Proj_Jac_Tensor_j[3][1] += factor*val_Proj_Visc_Flux[1];
    val_Proj_Jac_Tensor_i[3][2] += factor*val_Proj_Visc_Flux[2];
    val_Proj_Jac_Tensor_j[3][2] += factor*val_Proj_Visc_Flux[2];
    
    
  } else {

    val_Proj_Jac_Tensor_i[0][0] = 0.0;
    val_Proj_Jac_Tensor_i[0][1] = 0.0;
    val_Proj_Jac_Tensor_i[0][2] = 0.0;
    val_Proj_Jac_Tensor_i[0][3] = 0.0;
    val_Proj_Jac_Tensor_i[0][4] = 0.0;
    val_Proj_Jac_Tensor_i[1][0] = val_dS*tau_jacobian_i[0][0];
    val_Proj_Jac_Tensor_i[1][1] = val_dS*tau_jacobian_i[0][1];
    val_Proj_Jac_Tensor_i[1][2] = val_dS*tau_jacobian_i[0][2];
    val_Proj_Jac_Tensor_i[1][3] = val_dS*tau_jacobian_i[0][3];
    val_Proj_Jac_Tensor_i[1][4] = val_dS*tau_jacobian_i[0][4];
    val_Proj_Jac_Tensor_i[2][0] = val_dS*tau_jacobian_i[1][0];
    val_Proj_Jac_Tensor_i[2][1] = val_dS*tau_jacobian_i[1][1];
    val_Proj_Jac_Tensor_i[2][2] = val_dS*tau_jacobian_i[1][2];
    val_Proj_Jac_Tensor_i[2][3] = val_dS*tau_jacobian_i[1][3];
    val_Proj_Jac_Tensor_i[2][4] = val_dS*tau_jacobian_i[1][4];
    val_Proj_Jac_Tensor_i[3][0] = val_dS*tau_jacobian_i[2][0];
    val_Proj_Jac_Tensor_i[3][1] = val_dS*tau_jacobian_i[2][1];
    val_Proj_Jac_Tensor_i[3][2] = val_dS*tau_jacobian_i[2][2];
    val_Proj_Jac_Tensor_i[3][3] = val_dS*tau_jacobian_i[2][3];
    val_Proj_Jac_Tensor_i[3][4] = val_dS*tau_jacobian_i[2][4];
    const su2double contraction = tau_jacobian_i[0][0]*val_Mean_PrimVar[1] +
                                  tau_jacobian_i[1][0]*val_Mean_PrimVar[2] +
                                  tau_jacobian_i[2][0]*val_Mean_PrimVar[3];
    val_Proj_Jac_Tensor_i[4][0] = val_dS*(contraction - heat_flux_jac_i[0]);
    val_Proj_Jac_Tensor_i[4][1] = -val_dS*(tau_jacobian_i[0][0] + heat_flux_jac_i[1]); 
    val_Proj_Jac_Tensor_i[4][2] = -val_dS*(tau_jacobian_i[1][0] + heat_flux_jac_i[2]); 
    val_Proj_Jac_Tensor_i[4][3] = -val_dS*(tau_jacobian_i[2][0] + heat_flux_jac_i[3]); 
    val_Proj_Jac_Tensor_i[4][4] = -val_dS*heat_flux_jac_i[4];

    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      for (unsigned short jVar = 0; jVar < nVar; jVar++)
        val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];
    
    const su2double proj_viscousflux_vel= val_Proj_Visc_Flux[1]*val_Mean_PrimVar[1] + 
                                          val_Proj_Visc_Flux[2]*val_Mean_PrimVar[2] +
                                          val_Proj_Visc_Flux[3]*val_Mean_PrimVar[3];
    val_Proj_Jac_Tensor_i[4][0] -= factor*proj_viscousflux_vel;
    val_Proj_Jac_Tensor_j[4][0] -= factor*proj_viscousflux_vel;
    val_Proj_Jac_Tensor_i[4][1] += factor*val_Proj_Visc_Flux[1];
    val_Proj_Jac_Tensor_j[4][1] += factor*val_Proj_Visc_Flux[1];
    val_Proj_Jac_Tensor_i[4][2] += factor*val_Proj_Visc_Flux[2];
    val_Proj_Jac_Tensor_j[4][2] += factor*val_Proj_Visc_Flux[2];
    val_Proj_Jac_Tensor_i[4][3] += factor*val_Proj_Visc_Flux[3];
    val_Proj_Jac_Tensor_j[4][3] += factor*val_Proj_Visc_Flux[3];

  }

}
