/*!
 * \file flow_diffusion.cpp
 * \brief Implementation of numerics classes for discretization
 *        of viscous fluxes in fluid flow problems.
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

#include "../../../include/numerics/flow/flow_diffusion.hpp"

CAvgGrad_Base::CAvgGrad_Base(unsigned short val_nDim,
                             unsigned short val_nVar,
                             unsigned short val_nPrimVar,
                             bool val_correct_grad,
                             const CConfig* config)
    : CNumerics(val_nDim, val_nVar, config),
      nPrimVar(val_nPrimVar),
      correct_gradient(val_correct_grad) {

  unsigned short iVar, iDim;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  TauWall_i = -1.0; TauWall_j = -1.0;

  Mean_PrimVar = new su2double [nPrimVar];

  Mean_GradPrimVar = new su2double* [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++)
    Mean_GradPrimVar[iVar] = new su2double [nDim];
        
  Mean_GradTurbVar = new su2double[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Mean_GradTurbVar[iDim] = 0.;
  }

  Proj_Mean_GradPrimVar_Edge = new su2double[val_nPrimVar];
  Proj_Mean_GradTurbVar_Edge = 0.0;

  tau_jacobian_i = new su2double* [nDim];
  tau_jacobian_j = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    tau_jacobian_i[iDim] = new su2double [nVar];
    tau_jacobian_j[iDim] = new su2double [nVar];
  }

  heat_flux_jac_i = new su2double[nVar];
  heat_flux_jac_j = new su2double[nVar];

  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  Jacobian_ic = new su2double** [nDim];
  Jacobian_jc = new su2double** [nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Jacobian_ic[iDim] = new su2double* [nVar];
    Jacobian_jc[iDim] = new su2double* [nVar];
  }
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Jacobian_ic[iDim][iVar] = new su2double [nVar];
      Jacobian_jc[iDim][iVar] = new su2double [nVar];
    }
  }

}

CAvgGrad_Base::~CAvgGrad_Base() {

  delete [] Mean_PrimVar;

  if (Mean_GradPrimVar != nullptr) {
    for (unsigned short iVar = 0; iVar < nPrimVar; iVar++)
      delete [] Mean_GradPrimVar[iVar];
    delete [] Mean_GradPrimVar;
  }
  
  if (Mean_GradTurbVar != nullptr) delete [] Mean_GradTurbVar;

  delete [] Proj_Mean_GradPrimVar_Edge;

  if (tau_jacobian_i != nullptr) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      delete [] tau_jacobian_i[iDim];
    }
    delete [] tau_jacobian_i;
  }
  if (tau_jacobian_j != nullptr) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      delete [] tau_jacobian_j[iDim];
    }
    delete [] tau_jacobian_j;
  }

  delete [] heat_flux_jac_i;
  delete [] heat_flux_jac_j;

  if (Jacobian_i != nullptr) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete [] Jacobian_i[iVar];
      delete [] Jacobian_j[iVar];
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        delete [] Jacobian_ic[iDim][iVar];
        delete [] Jacobian_jc[iDim][iVar];
      }
    }
    delete [] Jacobian_i;
    delete [] Jacobian_j;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      delete [] Jacobian_ic[iDim];
      delete [] Jacobian_jc[iDim];
    }
    delete [] Jacobian_ic;
    delete [] Jacobian_jc;
  }

}

void CAvgGrad_Base::CorrectGradient(su2double** GradPrimVar,
                                    const su2double* val_PrimVar_i,
                                    const su2double* val_PrimVar_j,
                                    const su2double* val_edge_vector,
                                    const su2double val_proj_vector,
                                    const unsigned short val_nPrimVar) {
  for (unsigned short iVar = 0; iVar < val_nPrimVar; iVar++) {
    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Proj_Mean_GradPrimVar_Edge[iVar] += GradPrimVar[iVar][iDim]*val_edge_vector[iDim];
    }
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                 (val_PrimVar_j[iVar]-val_PrimVar_i[iVar]))*val_edge_vector[iDim] / val_proj_vector;
    }
  }

  // if (TurbVar_Grad_i != NULL) {
  //   Proj_Mean_GradTurbVar_Edge = 0.0;
  //   for (unsigned short iDim = 0; iDim < nDim; iDim++) {
  //     Proj_Mean_GradTurbVar_Edge += Mean_GradTurbVar[iDim]*val_edge_vector[iDim];
  //   }
  //   for (unsigned short iDim = 0; iDim < nDim; iDim++) {
  //     Mean_GradTurbVar[iDim] -= (Proj_Mean_GradTurbVar_Edge -
  //                                (turb_ke_j-turb_ke_i))*val_edge_vector[iDim] / val_proj_vector;
  //   }
  // }
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
                        - TWO3*total_viscosity*div_vel*delta[iDim][jDim] - TWO3*Density*val_turb_ke*delta[iDim][jDim];
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
  su2double TauNormal, TauElem[3], TauTangent[3], Area, UnitNormal[3];

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

void CAvgGrad_Base::SetPerturbedRSM(su2double turb_ke, const CConfig* config){

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
                                   const su2double val_proj_vector,
                                   const su2double val_area,
                                   const su2double *val_normal) {

  /*--- QCR and wall functions are **not** accounted for here ---*/
  /*--- BCM: account for wall functions ---*/
  
  const su2double WF_Factor = (Mean_TauWall > 0) ? Mean_TauWall/WallShearStress : su2double(1.0);
  const su2double Density_i = V_i[nDim+2];
  const su2double Density_j = V_j[nDim+2];
  const su2double Viscosity = val_laminar_viscosity + val_eddy_viscosity;
  const su2double xi_i = WF_Factor*Viscosity/(Density_i*dist_ij_2);
  const su2double xi_j = WF_Factor*Viscosity/(Density_j*dist_ij_2);

  /*--- BCM: contribution of TKE viscous stress in energy equation ---*/

  su2double xi_k_i = 0.0, xi_k_j = 0.0;

  // if (TurbVar_Grad_i != NULL) {
  //   const su2double sigma_k1 = 0.85;
  //   const su2double sigma_k2 = 1.0;

  //   const su2double sigma_k_i = F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2;
  //   const su2double sigma_k_j = F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2;

  //   const su2double Visc_k_i = Laminar_Viscosity_i + sigma_k_i*Eddy_Viscosity_i;
  //   const su2double Visc_k_j = Laminar_Viscosity_j + sigma_k_j*Eddy_Viscosity_j;

  //   const su2double Mean_Visc_k = 0.5*(Visc_k_i+Visc_k_j);

  //   xi_k_i = Mean_Visc_k/(Density_i*dist_ij_2);
  //   xi_k_j = Mean_Visc_k/(Density_j*dist_ij_2);
  // }

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      // Jacobian w.r.t. momentum
      tau_jacobian_i[iDim][jDim+1] = -xi_i*(Edge_Vector[iDim]*Normal[jDim] 
                                          - 2./3.*Edge_Vector[jDim]*Normal[iDim] 
                                          + delta[iDim][jDim]*val_proj_vector);
      tau_jacobian_j[iDim][jDim+1] =  xi_j*(Edge_Vector[iDim]*Normal[jDim] 
                                          - 2./3.*Edge_Vector[jDim]*Normal[iDim] 
                                          + delta[iDim][jDim]*val_proj_vector);
    }
    // Jacobian w.r.t. density
    tau_jacobian_i[iDim][0] = 0;
    tau_jacobian_j[iDim][0] = 0;
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      tau_jacobian_i[iDim][0] -= tau_jacobian_i[iDim][jDim+1]*V_i[jDim+1];
      tau_jacobian_j[iDim][0] -= tau_jacobian_j[iDim][jDim+1]*V_j[jDim+1];

      // tau_jacobian_i[iDim][jDim+1] -= xi_k_i*(Edge_Vector[iDim]*Normal[jDim] 
      //                                 - 2./3.*Edge_Vector[jDim]*Normal[iDim] 
      //                                     + delta[iDim][jDim]*val_proj_vector)*turb_ke_i;
      // tau_jacobian_j[iDim][jDim+1] += xi_k_j*(Edge_Vector[iDim]*Normal[jDim] 
      //                                     - 2./3.*Edge_Vector[jDim]*Normal[iDim] 
      //                                     + delta[iDim][jDim]*val_proj_vector)*turb_ke_j;
    }
    // Jacobian w.r.t. energy
    tau_jacobian_i[iDim][nDim+1] = 0;
    tau_jacobian_j[iDim][nDim+1] = 0;
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

  // if (TurbVar_Grad_i != NULL) {
  //   const su2double sigma_k1 = 0.85;
  //   const su2double sigma_k2 = 1.0;

  //   const su2double sigma_k_i = F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2;
  //   const su2double sigma_k_j = F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2;

  //   const su2double Visc_k_i = Laminar_Viscosity_i + sigma_k_i*Eddy_Viscosity_i;
  //   const su2double Visc_k_j = Laminar_Viscosity_j + sigma_k_j*Eddy_Viscosity_j;

  //   const su2double Mean_Visc_k = 0.5*(Visc_k_i+Visc_k_j);

  //   for (unsigned short iDim = 0; iDim < nDim; iDim++) {
  //     Flux_Tensor[nVar-1][iDim] += Mean_Visc_k*Mean_GradTurbVar[iDim];
  //   }
  // }

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Proj_Flux_Tensor[iVar] = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim] * val_normal[iDim];
  }

}

void CAvgGrad_Base::GetViscousProjJacs(const su2double *val_Mean_PrimVar,
                                       su2double val_proj_vector,
                                       su2double val_dS,
                                       const su2double *val_Proj_Visc_Flux,
                                       su2double **val_Proj_Jac_Tensor_i,
                                       su2double **val_Proj_Jac_Tensor_j,
                                       const CConfig *config) {
  
  const su2double Density = val_Mean_PrimVar[nDim+2];
  const su2double factor_i = 0.5/V_i[nDim+2];
  const su2double factor_j = 0.5/V_j[nDim+2];

  if (nDim == 2) {

    Jacobian_i[1][0] = tau_jacobian_i[0][0];
    Jacobian_i[1][1] = tau_jacobian_i[0][1];
    Jacobian_i[1][2] = tau_jacobian_i[0][2];
    Jacobian_i[2][0] = tau_jacobian_i[1][0];
    Jacobian_i[2][1] = tau_jacobian_i[1][1];
    Jacobian_i[2][2] = tau_jacobian_i[1][2];

    Jacobian_j[1][0] = tau_jacobian_j[0][0];
    Jacobian_j[1][1] = tau_jacobian_j[0][1];
    Jacobian_j[1][2] = tau_jacobian_j[0][2];
    Jacobian_j[2][0] = tau_jacobian_j[1][0];
    Jacobian_j[2][1] = tau_jacobian_j[1][1];
    Jacobian_j[2][2] = tau_jacobian_j[1][2];
    
    const su2double contraction_i = tau_jacobian_i[0][0]*val_Mean_PrimVar[1] +
                                    tau_jacobian_i[1][0]*val_Mean_PrimVar[2];
    
    Jacobian_i[3][0] = contraction_i + heat_flux_jac_i[0];
    Jacobian_i[3][1] = tau_jacobian_i[0][1]*val_Mean_PrimVar[1]
                     + tau_jacobian_i[0][2]*val_Mean_PrimVar[2]
                     + heat_flux_jac_i[1];
    Jacobian_i[3][2] = tau_jacobian_i[1][1]*val_Mean_PrimVar[1]
                     + tau_jacobian_i[1][2]*val_Mean_PrimVar[2]
                     + heat_flux_jac_i[2];
    Jacobian_i[3][3] = heat_flux_jac_i[3];
    
    const su2double contraction_j = tau_jacobian_j[0][0]*val_Mean_PrimVar[1] +
                                    tau_jacobian_j[1][0]*val_Mean_PrimVar[2];
    
    Jacobian_j[3][0] = contraction_j + heat_flux_jac_j[0];
    Jacobian_j[3][1] = tau_jacobian_j[0][1]*val_Mean_PrimVar[1]
                     + tau_jacobian_j[0][2]*val_Mean_PrimVar[2]
                     + heat_flux_jac_j[1];
    Jacobian_j[3][2] = tau_jacobian_j[1][1]*val_Mean_PrimVar[1]
                     + tau_jacobian_j[1][2]*val_Mean_PrimVar[2]
                     + heat_flux_jac_j[2];
    Jacobian_j[3][3] = heat_flux_jac_j[3];

    // CorrectJacobian(val_proj_vector, val_dS, val_Proj_Jac_Tensor_i, val_Proj_Jac_Tensor_j, Density, config);

    const su2double proj_flux_vel_i = val_Proj_Visc_Flux[1]*V_i[1] +
                                      val_Proj_Visc_Flux[2]*V_i[2];

    Jacobian_i[3][0] -= factor_i*proj_flux_vel_i;
    Jacobian_i[3][1] += factor_i*val_Proj_Visc_Flux[1];
    Jacobian_i[3][2] += factor_i*val_Proj_Visc_Flux[2];

    const su2double proj_flux_vel_j = val_Proj_Visc_Flux[1]*V_j[1] +
                                      val_Proj_Visc_Flux[2]*V_j[2];

    Jacobian_j[3][0] -= factor_j*proj_flux_vel_j;
    Jacobian_j[3][1] += factor_j*val_Proj_Visc_Flux[1];
    Jacobian_j[3][2] += factor_j*val_Proj_Visc_Flux[2];

  } else {

    Jacobian_i[1][0] = tau_jacobian_i[0][0];
    Jacobian_i[1][1] = tau_jacobian_i[0][1];
    Jacobian_i[1][2] = tau_jacobian_i[0][2];
    Jacobian_i[1][3] = tau_jacobian_i[0][3];
    Jacobian_i[2][0] = tau_jacobian_i[1][0];
    Jacobian_i[2][1] = tau_jacobian_i[1][1];
    Jacobian_i[2][2] = tau_jacobian_i[1][2];
    Jacobian_i[2][3] = tau_jacobian_i[1][3];
    Jacobian_i[3][0] = tau_jacobian_i[2][0];
    Jacobian_i[3][1] = tau_jacobian_i[2][1];
    Jacobian_i[3][2] = tau_jacobian_i[2][2];
    Jacobian_i[3][3] = tau_jacobian_i[2][3];

    Jacobian_j[1][0] = tau_jacobian_j[0][0];
    Jacobian_j[1][1] = tau_jacobian_j[0][1];
    Jacobian_j[1][2] = tau_jacobian_j[0][2];
    Jacobian_j[1][3] = tau_jacobian_j[0][3];
    Jacobian_j[2][0] = tau_jacobian_j[1][0];
    Jacobian_j[2][1] = tau_jacobian_j[1][1];
    Jacobian_j[2][2] = tau_jacobian_j[1][2];
    Jacobian_j[2][3] = tau_jacobian_j[1][3];
    Jacobian_j[3][0] = tau_jacobian_j[2][0];
    Jacobian_j[3][1] = tau_jacobian_j[2][1];
    Jacobian_j[3][2] = tau_jacobian_j[2][2];
    Jacobian_j[3][3] = tau_jacobian_j[2][3];
    
    const su2double contraction_i = tau_jacobian_i[0][0]*val_Mean_PrimVar[1] +
                                    tau_jacobian_i[1][0]*val_Mean_PrimVar[2] +
                                    tau_jacobian_i[2][0]*val_Mean_PrimVar[3];

    Jacobian_i[4][0] = contraction_i + heat_flux_jac_i[0];
    Jacobian_i[4][1] = tau_jacobian_i[0][1]*val_Mean_PrimVar[1]
                     + tau_jacobian_i[0][2]*val_Mean_PrimVar[2]
                     + tau_jacobian_i[0][3]*val_Mean_PrimVar[3]
                     + heat_flux_jac_i[1];
    Jacobian_i[4][2] = tau_jacobian_i[1][1]*val_Mean_PrimVar[1]
                     + tau_jacobian_i[1][2]*val_Mean_PrimVar[2]
                     + tau_jacobian_i[1][3]*val_Mean_PrimVar[3]
                     + heat_flux_jac_i[2];
    Jacobian_i[4][3] = tau_jacobian_i[2][1]*val_Mean_PrimVar[1]
                     + tau_jacobian_i[2][2]*val_Mean_PrimVar[2]
                     + tau_jacobian_i[2][3]*val_Mean_PrimVar[3]
                     + heat_flux_jac_i[3];
    Jacobian_i[4][4] = heat_flux_jac_i[4];
    
    const su2double contraction_j = tau_jacobian_j[0][0]*val_Mean_PrimVar[1] +
                                    tau_jacobian_j[1][0]*val_Mean_PrimVar[2] +
                                    tau_jacobian_j[2][0]*val_Mean_PrimVar[3];
    
    Jacobian_j[4][0] = contraction_j + heat_flux_jac_j[0];
    Jacobian_j[4][1] = tau_jacobian_j[0][1]*val_Mean_PrimVar[1]
                     + tau_jacobian_j[0][2]*val_Mean_PrimVar[2]
                     + tau_jacobian_j[0][3]*val_Mean_PrimVar[3]
                     + heat_flux_jac_j[1];
    Jacobian_j[4][2] = tau_jacobian_j[1][1]*val_Mean_PrimVar[1]
                     + tau_jacobian_j[1][2]*val_Mean_PrimVar[2]
                     + tau_jacobian_j[1][3]*val_Mean_PrimVar[3]
                     + heat_flux_jac_j[2];
    Jacobian_j[4][3] = tau_jacobian_j[2][1]*val_Mean_PrimVar[1]
                     + tau_jacobian_j[2][2]*val_Mean_PrimVar[2]
                     + tau_jacobian_j[2][3]*val_Mean_PrimVar[3]
                     + heat_flux_jac_j[3];
    Jacobian_j[4][4] = heat_flux_jac_j[4];

    // CorrectJacobian(val_proj_vector, val_dS, val_Proj_Jac_Tensor_i, val_Proj_Jac_Tensor_j, Density, config);

    const su2double proj_flux_vel_i = val_Proj_Visc_Flux[1]*V_i[1] +
                                      val_Proj_Visc_Flux[2]*V_i[2] +
                                      val_Proj_Visc_Flux[3]*V_i[3];

    Jacobian_i[4][0] -= factor_i*proj_flux_vel_i;
    Jacobian_i[4][1] += factor_i*val_Proj_Visc_Flux[1];
    Jacobian_i[4][2] += factor_i*val_Proj_Visc_Flux[2];
    Jacobian_i[4][3] += factor_i*val_Proj_Visc_Flux[3];

    const su2double proj_flux_vel_j = val_Proj_Visc_Flux[1]*V_j[1] +
                                      val_Proj_Visc_Flux[2]*V_j[2] +
                                      val_Proj_Visc_Flux[3]*V_j[3];

    Jacobian_j[4][0] -= factor_j*proj_flux_vel_j;
    Jacobian_j[4][1] += factor_j*val_Proj_Visc_Flux[1];
    Jacobian_j[4][2] += factor_j*val_Proj_Visc_Flux[2];
    Jacobian_j[4][3] += factor_j*val_Proj_Visc_Flux[3];

  }
  
}

void CAvgGrad_Base::CorrectJacobian(const su2double val_proj_vector,
                                    const su2double val_dS,
                                    su2double **val_Proj_Jac_Tensor_i,
                                    su2double **val_Proj_Jac_Tensor_j,
                                    const su2double Density,
                                    const CConfig *config) {
  
  AD_BEGIN_PASSIVE

  /*--- Add contributions of GG gradients ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    const su2double weight_i = 0.5 / Volume_i;
    const su2double weight_j = 0.5 / Volume_j;
    
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        for (unsigned short jVar = 0; jVar < nVar; jVar++) {
          Jacobian_ic[iDim][iVar][jVar] -= weight_i*(Normal[iDim]/val_proj_vector
                                         - Edge_Vector[iDim]) * Jacobian_i[iVar][jVar];
          Jacobian_jc[iDim][iVar][jVar] += weight_j*(Normal[iDim]/val_proj_vector
                                         - Edge_Vector[iDim]) * Jacobian_j[iVar][jVar];
        }
      }
    }
  }
  
  /*--- TODO: add contributions of WLS gradients ---*/
  else if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {}
  
  if (!correct_gradient) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      for (unsigned short jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] = 0.;
        Jacobian_j[iVar][jVar] = 0.;
      }
    }
  }
  
  AD_END_PASSIVE

}

CAvgGrad_Flow::CAvgGrad_Flow(unsigned short val_nDim,
                             unsigned short val_nVar,
                             bool val_correct_grad,
                             const CConfig* config)
    : CAvgGrad_Base(val_nDim, val_nVar, val_nDim+3, val_correct_grad, config) { }

CNumerics::ResidualType<> CAvgGrad_Flow::ComputeResidual(const CConfig* config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+9);   AD::SetPreaccIn(V_j, nDim+9);
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);
  AD::SetPreaccIn(PrimVar_Grad_j, nDim+1, nDim);
  AD::SetPreaccIn(turb_ke_i); AD::SetPreaccIn(turb_ke_j);
  AD::SetPreaccIn(Normal, nDim);

  // if (TurbVar_Grad_i != NULL) {
  //   AD::SetPreaccIn(TurbVar_Grad_i, 1, nDim); AD::SetPreaccIn(TurbVar_Grad_j, 1, nDim);
  //   AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F1_j);
  // }
  
  unsigned short iVar, jVar, iDim;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Flux_Tensor[iVar] = 0.0;
    for (jVar = 0; jVar < nVar; jVar++) {
      Jacobian_i[iVar][jVar] = 0.0;
      Jacobian_j[iVar][jVar] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Jacobian_ic[iDim][iVar][jVar] = 0.0;
        Jacobian_jc[iDim][iVar][jVar] = 0.0;
      }
    }
  }

  /*--- Normalized normal vector ---*/

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  PrimVar_i = V_i;
  PrimVar_j = V_j;

  for (iVar = 0; iVar < nPrimVar; iVar++) {
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }


  /*--- Compute vector going from iPoint to jPoint ---*/

  su2double proj_vector_ij = 0.0;
  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }

  /*--- Laminar and Eddy viscosity ---*/

  Laminar_Viscosity_i = V_i[nDim+5]; Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6]; Eddy_Viscosity_j = V_j[nDim+6];

  /*--- Mean Viscosities and turbulent kinetic energy---*/

  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke = 0.5*(turb_ke_i + turb_ke_j);

  /*--- Mean gradient approximation ---*/

  for (iVar = 0; iVar < nDim+1; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
    }
  }

  // if (TurbVar_Grad_i != NULL) {
  //   for (iDim = 0; iDim < nDim; iDim++) {
  //     Mean_GradTurbVar[iDim] = 0.5*(TurbVar_Grad_i[0][iDim] + TurbVar_Grad_j[0][iDim]);
  //   }
  // }
  
  

  /*--- Projection of the mean gradient in the direction of the edge ---*/

  if (correct_gradient) {
    CorrectGradient(Mean_GradPrimVar, PrimVar_i, PrimVar_j, Edge_Vector,
                    dist_ij_2, nDim+1);
  }

  /*--- Wall shear stress values (wall functions) ---*/

  if (TauWall_i > 0.0 && TauWall_j > 0.0) Mean_TauWall = 0.5*(TauWall_i + TauWall_j);
  else if (TauWall_i > 0.0) Mean_TauWall = TauWall_i;
  else if (TauWall_j > 0.0) Mean_TauWall = TauWall_j;
  else Mean_TauWall = -1.0;

  /* --- If using UQ methodology, set Reynolds Stress tensor and perform perturbation--- */

  if (using_uq){
    SetReynoldsStressMatrix(Mean_turb_ke);
    SetPerturbedRSM(Mean_turb_ke, config);
  }

  /*--- Get projected flux tensor (viscous residual) ---*/

  SetStressTensor(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke,
         Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  if (config->GetQCR()) AddQCR(Mean_GradPrimVar);
  if (Mean_TauWall > 0) AddTauWall(Normal, Mean_TauWall);

  SetHeatFluxVector(Mean_GradPrimVar, Mean_Laminar_Viscosity,
                    Mean_Eddy_Viscosity);

  GetViscousProjFlux(Mean_PrimVar, Normal);

  /*--- Compute the implicit part ---*/

  if (implicit && dist_ij_2 != 0.0) {

    SetTauJacobian(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, proj_vector_ij, Area, UnitNormal);
    SetHeatFluxJacobian(Mean_PrimVar, Mean_Laminar_Viscosity,
                        Mean_Eddy_Viscosity, proj_vector_ij, Area, UnitNormal);
    GetViscousProjJacs(Mean_PrimVar, proj_vector_ij, Area, Proj_Flux_Tensor, Jacobian_i, Jacobian_j, config);
    SetLaminarViscosityJacobian(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, Normal, config);
    SetEddyViscosityJacobian(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, Normal, config);

  }

  AD::SetPreaccOut(Proj_Flux_Tensor, nVar);
  AD::EndPreacc();

  return ResidualType<>(Proj_Flux_Tensor, Jacobian_i, Jacobian_j, Jacobian_ic, Jacobian_jc);

}

void CAvgGrad_Flow::SetHeatFluxVector(const su2double* const *val_gradprimvar,
                                      const su2double val_laminar_viscosity,
                                      const su2double val_eddy_viscosity) {

  const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  const su2double heat_flux_factor = Cp * (val_laminar_viscosity/Prandtl_Lam + val_eddy_viscosity/Prandtl_Turb);

  /*--- Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure] ---*/

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    heat_flux_vector[iDim] = heat_flux_factor*val_gradprimvar[0][iDim];
  }
}

void CAvgGrad_Flow::SetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                                        const su2double val_laminar_viscosity,
                                        const su2double val_eddy_viscosity,
                                        const su2double val_proj_vector,
                                        const su2double val_area,
                                        const su2double *val_normal) {
  su2double sqvel, Density, Pressure;
  su2double phi, R_dTdu0, R_dTdu1, R_dTdu2;

  const su2double heat_flux_factor = val_laminar_viscosity/Prandtl_Lam + val_eddy_viscosity/Prandtl_Turb;
  const su2double cpoR = Gamma/Gamma_Minus_One; // cp over R
  const su2double conductivity_over_Rd = cpoR*heat_flux_factor*val_proj_vector/dist_ij_2;

  sqvel = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    sqvel += V_i[iDim+1]*V_i[iDim+1];
  }

  Density = V_i[nDim+2];
  Pressure = V_i[nDim+1];
  phi = Gamma_Minus_One/Density;

  /*--- R times partial derivatives of temp. ---*/

  R_dTdu0 = -Pressure/(Density*Density) + 0.5*sqvel*phi;
  R_dTdu1 = -phi*V_i[1];
  R_dTdu2 = -phi*V_i[2];

  heat_flux_jac_i[0] = -conductivity_over_Rd * R_dTdu0;
  heat_flux_jac_i[1] = -conductivity_over_Rd * R_dTdu1;
  heat_flux_jac_i[2] = -conductivity_over_Rd * R_dTdu2;

  if (nDim == 2) {

    const su2double R_dTdu3 = phi;
    heat_flux_jac_i[3] = -conductivity_over_Rd * R_dTdu3;

  } else {

    const su2double R_dTdu3 = -phi*V_i[3];
    const su2double R_dTdu4 = phi;
    heat_flux_jac_i[3] = -conductivity_over_Rd * R_dTdu3;
    heat_flux_jac_i[4] = -conductivity_over_Rd * R_dTdu4;

  }
  
  sqvel = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    sqvel += V_j[iDim+1]*V_j[iDim+1];
  }

  Density = V_j[nDim+2];
  Pressure = V_j[nDim+1];
  phi = Gamma_Minus_One/Density;

  /*--- R times partial derivatives of temp. ---*/

  R_dTdu0 = -Pressure/(Density*Density) + 0.5*sqvel*phi;
  R_dTdu1 = -phi*V_j[1];
  R_dTdu2 = -phi*V_j[2];

  heat_flux_jac_j[0] = conductivity_over_Rd * R_dTdu0;
  heat_flux_jac_j[1] = conductivity_over_Rd * R_dTdu1;
  heat_flux_jac_j[2] = conductivity_over_Rd * R_dTdu2;

  if (nDim == 2) {

    const su2double R_dTdu3 = phi;
    heat_flux_jac_j[3] = conductivity_over_Rd * R_dTdu3;

  } else {

    const su2double R_dTdu3 = -phi*V_j[3];
    const su2double R_dTdu4 = phi;
    heat_flux_jac_j[3] = conductivity_over_Rd * R_dTdu3;
    heat_flux_jac_j[4] = conductivity_over_Rd * R_dTdu4;

  }
}

void CAvgGrad_Flow::SetLaminarViscosityJacobian(const su2double *val_Mean_PrimVar,
                                                const su2double val_laminar_viscosity,
                                                const su2double val_eddy_viscosity,
                                                const su2double *val_normal,
                                                const CConfig   *config) {
  
  su2double stress_tensor[3][3];
  su2double heat_flux_vec[3];
  su2double div_vel = 0.0;
  
  const su2double WF_Factor = (Mean_TauWall > 0) ? Mean_TauWall/WallShearStress : su2double(1.0);
//  const su2double WF_Factor = 1.0;
  
  const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  const su2double Cv = Cp/Gamma;
  const su2double heat_flux_factor = Cp/Prandtl_Lam;
  const su2double muref = config->GetMu_RefND();
  const su2double Tref = config->GetMu_Temperature_RefND();
  const su2double Sref = config->GetMu_SND();
  
  for (unsigned short iDim = 0 ; iDim < nDim; iDim++)
    div_vel += Mean_GradPrimVar[iDim+1][iDim];
  
  for (unsigned short iDim = 0 ; iDim < nDim; iDim++) {
    heat_flux_vec[iDim] = heat_flux_factor*Mean_GradPrimVar[0][iDim];
    for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
      stress_tensor[iDim][jDim] = WF_Factor * ( Mean_GradPrimVar[jDim+1][iDim] + Mean_GradPrimVar[iDim+1][jDim]
                                            - TWO3*div_vel*delta[iDim][jDim] );
  }

  /*--- Jacobian wrt laminar viscosity ---*/
  
  su2double v2_i = 0., v2_j = 0.;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    v2_i += V_i[iDim+1]*V_i[iDim+1];
    v2_j += V_j[iDim+1]*V_j[iDim+1];
  }
  su2double T_i = V_i[0],
            T_j = V_j[0],
            r_i = V_i[nDim+2],
            r_j = V_j[nDim+2];
  
  su2double dmudT_i  = muref*(Tref+Sref)/pow(Tref,1.5) * (3.*Sref*sqrt(T_i) + pow(T_i,1.5))/(2.*pow((T_i+Sref),2.)),
            dmudT_j  = muref*(Tref+Sref)/pow(Tref,1.5) * (3.*Sref*sqrt(T_j) + pow(T_j,1.5))/(2.*pow((T_j+Sref),2.)),
            factor_i = dmudT_i/(r_i*Cv),
            factor_j = dmudT_j/(r_j*Cv);
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      for (unsigned short kDim = 0; kDim < nDim; kDim++) {
        Jacobian_i[iDim+1][kDim+1] += -0.5*factor_i*V_i[kDim+1]*stress_tensor[iDim][jDim]*Normal[jDim];
        Jacobian_i[nDim+1][kDim+1] += -0.5*factor_i*V_i[kDim+1]*stress_tensor[iDim][jDim]*val_Mean_PrimVar[iDim+1]*Normal[jDim];
        
        Jacobian_j[iDim+1][kDim+1] += -0.5*factor_j*V_j[kDim+1]*stress_tensor[iDim][jDim]*Normal[jDim];
        Jacobian_j[nDim+1][kDim+1] += -0.5*factor_j*V_j[kDim+1]*stress_tensor[iDim][jDim]*val_Mean_PrimVar[iDim+1]*Normal[jDim];
      }
      Jacobian_i[iDim+1][0]      += 0.5*factor_i*(v2_i/2. - Cv*T_i)*stress_tensor[iDim][jDim]*Normal[jDim];
      Jacobian_i[iDim+1][nDim+1] += 0.5*factor_i*stress_tensor[iDim][jDim]*Normal[jDim];
      Jacobian_i[nDim+1][nDim+1] += 0.5*factor_i*stress_tensor[iDim][jDim]*val_Mean_PrimVar[iDim+1]*Normal[jDim];
      Jacobian_i[nDim+1][iDim+1] += -0.5*factor_i*V_i[iDim+1]*heat_flux_vec[jDim]*Normal[jDim];
      
      Jacobian_j[iDim+1][0]      += 0.5*factor_j*(v2_j/2. - Cv*T_j)*stress_tensor[iDim][jDim]*Normal[jDim];
      Jacobian_j[iDim+1][nDim+1] += 0.5*factor_j*stress_tensor[iDim][jDim]*Normal[jDim];
      Jacobian_j[nDim+1][nDim+1] += 0.5*factor_j*stress_tensor[iDim][jDim]*val_Mean_PrimVar[iDim+1]*Normal[jDim];
      Jacobian_j[nDim+1][iDim+1] += -0.5*factor_j*V_j[iDim+1]*heat_flux_vec[jDim]*Normal[jDim];
    }
    Jacobian_i[nDim+1][0]      += 0.5*factor_i*(v2_i/2. - Cv*T_i)*heat_flux_vec[iDim]*Normal[iDim];
    Jacobian_i[nDim+1][nDim+1] += 0.5*factor_i*heat_flux_vec[iDim]*Normal[iDim];
    
    Jacobian_j[nDim+1][0]      += 0.5*factor_j*(v2_j/2. - Cv*T_j)*heat_flux_vec[iDim]*Normal[iDim];
    Jacobian_j[nDim+1][nDim+1] += 0.5*factor_j*heat_flux_vec[iDim]*Normal[iDim];
  }
}

void CAvgGrad_Flow::SetEddyViscosityJacobian(const su2double *val_Mean_PrimVar,
                                             const su2double val_laminar_viscosity,
                                             const su2double val_eddy_viscosity,
                                             const su2double *val_normal,
                                             const CConfig   *config) {
  
  if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)) {
    
    const su2double WF_Factor = (Mean_TauWall > 0) ? Mean_TauWall/WallShearStress : su2double(1.0);
//    const su2double WF_Factor = 1.0;
        
    const su2double a1 = 0.31;
    
    su2double stress_tensor[3][3];
    su2double heat_flux_vec[3];
    su2double div_vel = 0.0;
    const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
    const su2double heat_flux_factor = Cp/Prandtl_Turb;
    
    for (unsigned short iDim = 0 ; iDim < nDim; iDim++)
      div_vel += Mean_GradPrimVar[iDim+1][iDim];
    
    for (unsigned short iDim = 0 ; iDim < nDim; iDim++) {
      heat_flux_vec[iDim] = heat_flux_factor*Mean_GradPrimVar[0][iDim];
      for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
        stress_tensor[iDim][jDim] = WF_Factor * ( Mean_GradPrimVar[jDim+1][iDim] + Mean_GradPrimVar[iDim+1][jDim]
                                              - TWO3*div_vel*delta[iDim][jDim] );
    }
  
    /*--- Jacobian wrt eddy viscosity ---*/
    
    const su2double VorticityMag_i = sqrt(Vorticity_i[0]*Vorticity_i[0] +
                                          Vorticity_i[1]*Vorticity_i[1] +
                                          Vorticity_i[2]*Vorticity_i[2]);

    if (turb_omega_i > VorticityMag_i*F2_i/a1) {
      const su2double factor = turb_ke_i/turb_omega_i;
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        for (unsigned short jDim = 0; jDim < nDim; jDim++) {
          Jacobian_i[iDim+1][0] += 0.5*factor*stress_tensor[iDim][jDim]*Normal[jDim];
          Jacobian_i[nDim+1][0] += 0.5*factor*stress_tensor[iDim][jDim]*val_Mean_PrimVar[iDim+1]*Normal[jDim];
        }
        Jacobian_i[nDim+1][0] += 0.5*factor*heat_flux_vec[iDim]*Normal[iDim];
      }
    }

    const su2double VorticityMag_j = sqrt(Vorticity_j[0]*Vorticity_j[0] +
                                          Vorticity_j[1]*Vorticity_j[1] +
                                          Vorticity_j[2]*Vorticity_j[2]);

    if (turb_omega_j > VorticityMag_j*F2_j/a1) {
      const su2double factor = turb_ke_j/turb_omega_j;
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        for (unsigned short jDim = 0; jDim < nDim; jDim++) {
          Jacobian_j[iDim+1][0] += 0.5*factor*stress_tensor[iDim][jDim]*Normal[jDim];
          Jacobian_j[nDim+1][0] += 0.5*factor*stress_tensor[iDim][jDim]*val_Mean_PrimVar[iDim+1]*Normal[jDim];
        }
        Jacobian_j[nDim+1][0] += 0.5*factor*heat_flux_vec[iDim]*Normal[iDim];
      }
    }
  }
  
}

CAvgGradInc_Flow::CAvgGradInc_Flow(unsigned short val_nDim,
                                   unsigned short val_nVar,
                                   bool val_correct_grad, const CConfig* config)
    : CAvgGrad_Base(val_nDim, val_nVar, val_nDim+3, val_correct_grad, config) {

  energy   = config->GetEnergy_Equation();

}

CNumerics::ResidualType<> CAvgGradInc_Flow::ComputeResidual(const CConfig* config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+9);   AD::SetPreaccIn(V_j, nDim+9);
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(PrimVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(PrimVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(turb_ke_i); AD::SetPreaccIn(turb_ke_j);
  AD::SetPreaccIn(Normal, nDim);

  unsigned short iVar, jVar, iDim;

  /*--- Normalized normal vector ---*/

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  PrimVar_i = V_i;
  PrimVar_j = V_j;

  for (iVar = 0; iVar < nPrimVar; iVar++) {
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }

  /*--- Compute vector going from iPoint to jPoint ---*/

  su2double proj_vector_ij = 0.0;
  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }

  /*--- Density and transport properties ---*/

  Laminar_Viscosity_i    = V_i[nDim+4];  Laminar_Viscosity_j    = V_j[nDim+4];
  Eddy_Viscosity_i       = V_i[nDim+5];  Eddy_Viscosity_j       = V_j[nDim+5];
  Thermal_Conductivity_i = V_i[nDim+6];  Thermal_Conductivity_j = V_j[nDim+6];

  /*--- Mean transport properties ---*/

  Mean_Laminar_Viscosity    = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity       = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke              = 0.5*(turb_ke_i + turb_ke_j);
  Mean_Thermal_Conductivity = 0.5*(Thermal_Conductivity_i + Thermal_Conductivity_j);

  /*--- Mean gradient approximation ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    for (iDim = 0; iDim < nDim; iDim++)
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);

  /*--- Projection of the mean gradient in the direction of the edge ---*/

  if (correct_gradient && dist_ij_2 != 0.0) {
    CorrectGradient(Mean_GradPrimVar, PrimVar_i, PrimVar_j, Edge_Vector,
                    proj_vector_ij, nVar);
  }

  /*--- Get projected flux tensor (viscous residual) ---*/
  SetStressTensor(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);

  GetViscousIncProjFlux(Mean_GradPrimVar, Normal, Mean_Thermal_Conductivity);

  /*--- Implicit part ---*/

  if (implicit) {

    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    } else {

      const su2double dist_ij = sqrt(dist_ij_2);
      SetIncTauJacobian(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, dist_ij, UnitNormal);

      GetViscousIncProjJacs(Area, Jacobian_i, Jacobian_j);

      /*--- Include the temperature equation Jacobian. ---*/
      su2double proj_vector_ij = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        proj_vector_ij += (Coord_j[iDim]-Coord_i[iDim])*Normal[iDim];
      }
      proj_vector_ij = proj_vector_ij/dist_ij_2;
      Jacobian_i[nDim+1][nDim+1] = -Mean_Thermal_Conductivity*proj_vector_ij;
      Jacobian_j[nDim+1][nDim+1] =  Mean_Thermal_Conductivity*proj_vector_ij;
    }

  }

  if (!energy) {
    Proj_Flux_Tensor[nDim+1] = 0.0;
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[iVar][nDim+1] = 0.0;
        Jacobian_j[iVar][nDim+1] = 0.0;

        Jacobian_i[nDim+1][iVar] = 0.0;
        Jacobian_j[nDim+1][iVar] = 0.0;
      }
    }
  }

  AD::SetPreaccOut(Proj_Flux_Tensor, nVar);
  AD::EndPreacc();

  return ResidualType<>(Proj_Flux_Tensor, Jacobian_i, Jacobian_j);

}

void CAvgGradInc_Flow::GetViscousIncProjFlux(const su2double* const *val_gradprimvar,
                                             const su2double *val_normal,
                                             su2double val_thermal_conductivity) {

  /*--- Gradient of primitive variables -> [Pressure vel_x vel_y vel_z Temperature] ---*/

  if (nDim == 2) {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = val_thermal_conductivity*val_gradprimvar[nDim+1][0];

    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = val_thermal_conductivity*val_gradprimvar[nDim+1][1];

  } else {

    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][2];
    Flux_Tensor[4][0] = val_thermal_conductivity*val_gradprimvar[nDim+1][0];

    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][2];
    Flux_Tensor[4][1] = val_thermal_conductivity*val_gradprimvar[nDim+1][1];

    Flux_Tensor[0][2] = 0.0;
    Flux_Tensor[1][2] = tau[2][0];
    Flux_Tensor[2][2] = tau[2][1];
    Flux_Tensor[3][2] = tau[2][2];
    Flux_Tensor[4][2] = val_thermal_conductivity*val_gradprimvar[nDim+1][2];

  }

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Proj_Flux_Tensor[iVar] = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim] * val_normal[iDim];
  }

}

void CAvgGradInc_Flow::GetViscousIncProjJacs(su2double val_dS,
                                             su2double **val_Proj_Jac_Tensor_i,
                                             su2double **val_Proj_Jac_Tensor_j) {
  if (nDim == 2) {

    val_Proj_Jac_Tensor_i[0][0] = 0.0;
    val_Proj_Jac_Tensor_i[0][1] = 0.0;
    val_Proj_Jac_Tensor_i[0][2] = 0.0;
    val_Proj_Jac_Tensor_i[0][3] = 0.0;

    val_Proj_Jac_Tensor_i[1][0] = tau_jacobian_i[0][0];
    val_Proj_Jac_Tensor_i[1][1] = tau_jacobian_i[0][1];
    val_Proj_Jac_Tensor_i[1][2] = tau_jacobian_i[0][2];
    val_Proj_Jac_Tensor_i[1][3] = tau_jacobian_i[0][3];

    val_Proj_Jac_Tensor_i[2][0] = tau_jacobian_i[1][0];
    val_Proj_Jac_Tensor_i[2][1] = tau_jacobian_i[1][1];
    val_Proj_Jac_Tensor_i[2][2] = tau_jacobian_i[1][2];
    val_Proj_Jac_Tensor_i[2][3] = tau_jacobian_i[1][3];

    val_Proj_Jac_Tensor_i[3][0] = 0.0;
    val_Proj_Jac_Tensor_i[3][1] = 0.0;
    val_Proj_Jac_Tensor_i[3][2] = 0.0;
    val_Proj_Jac_Tensor_i[3][3] = 0.0;

  } else {

    val_Proj_Jac_Tensor_i[0][0] = 0.0;
    val_Proj_Jac_Tensor_i[0][1] = 0.0;
    val_Proj_Jac_Tensor_i[0][2] = 0.0;
    val_Proj_Jac_Tensor_i[0][3] = 0.0;
    val_Proj_Jac_Tensor_i[0][4] = 0.0;

    val_Proj_Jac_Tensor_i[1][0] = tau_jacobian_i[0][0];
    val_Proj_Jac_Tensor_i[1][1] = tau_jacobian_i[0][1];
    val_Proj_Jac_Tensor_i[1][2] = tau_jacobian_i[0][2];
    val_Proj_Jac_Tensor_i[1][3] = tau_jacobian_i[0][3];
    val_Proj_Jac_Tensor_i[1][4] = tau_jacobian_i[0][4];

    val_Proj_Jac_Tensor_i[2][0] = tau_jacobian_i[1][0];
    val_Proj_Jac_Tensor_i[2][1] = tau_jacobian_i[1][1];
    val_Proj_Jac_Tensor_i[2][2] = tau_jacobian_i[1][2];
    val_Proj_Jac_Tensor_i[2][3] = tau_jacobian_i[1][3];
    val_Proj_Jac_Tensor_i[2][4] = tau_jacobian_i[1][4];

    val_Proj_Jac_Tensor_i[3][0] = tau_jacobian_i[2][0];
    val_Proj_Jac_Tensor_i[3][1] = tau_jacobian_i[2][1];
    val_Proj_Jac_Tensor_i[3][2] = tau_jacobian_i[2][2];
    val_Proj_Jac_Tensor_i[3][3] = tau_jacobian_i[2][3];
    val_Proj_Jac_Tensor_i[3][4] = tau_jacobian_i[2][4];

    val_Proj_Jac_Tensor_i[4][0] = 0.0;
    val_Proj_Jac_Tensor_i[4][1] = 0.0;
    val_Proj_Jac_Tensor_i[4][2] = 0.0;
    val_Proj_Jac_Tensor_i[4][3] = 0.0;
    val_Proj_Jac_Tensor_i[4][4] = 0.0;

  }

  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    for (unsigned short jVar = 0; jVar < nVar; jVar++)
      val_Proj_Jac_Tensor_j[iVar][jVar] = -val_Proj_Jac_Tensor_i[iVar][jVar];

}

CGeneralAvgGrad_Flow::CGeneralAvgGrad_Flow(unsigned short val_nDim,
                                           unsigned short val_nVar,
                                           bool val_correct_grad,
                                           const CConfig* config)
    : CAvgGrad_Base(val_nDim, val_nVar, val_nDim+4, val_correct_grad, config) { }

void CGeneralAvgGrad_Flow::SetHeatFluxVector(const su2double* const *val_gradprimvar,
                                             const su2double val_laminar_viscosity,
                                             const su2double val_eddy_viscosity,
                                             const su2double val_thermal_conductivity,
                                             const su2double val_heat_capacity_cp) {

  const su2double heat_flux_factor = val_thermal_conductivity + val_heat_capacity_cp*val_eddy_viscosity/Prandtl_Turb;

  /*--- Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure] ---*/
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    heat_flux_vector[iDim] = heat_flux_factor*val_gradprimvar[0][iDim];
  }
}

void CGeneralAvgGrad_Flow::SetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                                               const su2double *val_Mean_SecVar,
                                               const su2double val_eddy_viscosity,
                                               const su2double val_thermal_conductivity,
                                               const su2double val_heat_capacity_cp,
                                               const su2double val_proj_vector,
                                               const su2double val_area) {
  /* Viscous flux Jacobians for arbitrary equations of state */

  //order of val_mean_primitives: T, vx, vy, vz, P, rho, ht
  //order of secondary:dTdrho_e, dTde_rho

  su2double sqvel = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
  }

  su2double rho = val_Mean_PrimVar[nDim+2];
  su2double P= val_Mean_PrimVar[nDim+1];
  su2double h= val_Mean_PrimVar[nDim+3];
  su2double dTdrho_e= val_Mean_SecVar[0];
  su2double dTde_rho= val_Mean_SecVar[1];

  su2double dTdu0= dTdrho_e + dTde_rho*(-(h-P/rho) + sqvel)*(1/rho);
  su2double dTdu1= dTde_rho*(-val_Mean_PrimVar[1])*(1/rho);
  su2double dTdu2= dTde_rho*(-val_Mean_PrimVar[2])*(1/rho);

  su2double total_conductivity = val_thermal_conductivity + val_heat_capacity_cp*val_eddy_viscosity/Prandtl_Turb;
  su2double factor2 = total_conductivity*val_proj_vector;

  heat_flux_jac_i[0] = factor2*dTdu0;
  heat_flux_jac_i[1] = factor2*dTdu1;
  heat_flux_jac_i[2] = factor2*dTdu2;

  if (nDim == 2) {

    su2double dTdu3= dTde_rho*(1/rho);
    heat_flux_jac_i[3] = factor2*dTdu3;

  } else {

    su2double dTdu3= dTde_rho*(-val_Mean_PrimVar[3])*(1/rho);
    su2double dTdu4= dTde_rho*(1/rho);
    heat_flux_jac_i[3] = factor2*dTdu3;
    heat_flux_jac_i[4] = factor2*dTdu4;

  }

}

CNumerics::ResidualType<> CGeneralAvgGrad_Flow::ComputeResidual(const CConfig* config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+9);   AD::SetPreaccIn(V_j, nDim+9);
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(S_i, 4); AD::SetPreaccIn(S_j, 4);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);
  AD::SetPreaccIn(PrimVar_Grad_j, nDim+1, nDim);
  AD::SetPreaccIn(turb_ke_i); AD::SetPreaccIn(turb_ke_j);
  AD::SetPreaccIn(Normal, nDim);

  unsigned short iVar, jVar, iDim;

  /*--- Normalized normal vector ---*/

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Mean primitive variables ---*/

  PrimVar_i = V_i;
  PrimVar_j = V_j;

  for (iVar = 0; iVar < nPrimVar; iVar++) {
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }

  /*--- Compute vector going from iPoint to jPoint ---*/

  su2double proj_vector_ij = 0.0;
  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }

  /*--- Laminar and Eddy viscosity ---*/

  Laminar_Viscosity_i = V_i[nDim+5];    Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6];       Eddy_Viscosity_j = V_j[nDim+6];
  Thermal_Conductivity_i = V_i[nDim+7]; Thermal_Conductivity_j = V_j[nDim+7];
  Cp_i = V_i[nDim+8]; Cp_j = V_j[nDim+8];

  /*--- Mean secondary variables ---*/

  for (iVar = 0; iVar < 2; iVar++) {
    Mean_SecVar[iVar] = 0.5*(S_i[iVar+2]+S_j[iVar+2]);
  }

  /*--- Mean Viscosities and turbulent kinetic energy---*/

  Mean_Laminar_Viscosity    = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity       = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke              = 0.5*(turb_ke_i + turb_ke_j);
  Mean_Thermal_Conductivity = 0.5*(Thermal_Conductivity_i + Thermal_Conductivity_j);
  Mean_Cp                   = 0.5*(Cp_i + Cp_j);

  /*--- Mean gradient approximation ---*/

  for (iVar = 0; iVar < nDim+1; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
    }
  }

  /*--- Projection of the mean gradient in the direction of the edge ---*/

  if (correct_gradient && dist_ij_2 != 0.0) {
    CorrectGradient(Mean_GradPrimVar, PrimVar_i, PrimVar_j, Edge_Vector,
                    proj_vector_ij, nDim+1);
  }

  /* --- If using UQ methodology, set Reynolds Stress tensor and perform perturbation--- */

  if (using_uq){
    SetReynoldsStressMatrix(Mean_turb_ke);
    SetPerturbedRSM(Mean_turb_ke, config);
  }

  /*--- Get projected flux tensor (viscous residual) ---*/

  SetStressTensor(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke,
                  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);

  SetHeatFluxVector(Mean_GradPrimVar, Mean_Laminar_Viscosity,
                    Mean_Eddy_Viscosity, Mean_Thermal_Conductivity, Mean_Cp);

  GetViscousProjFlux(Mean_PrimVar, Normal);

  /*--- Compute the implicit part ---*/

  if (implicit) {

    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    } else {
      // const su2double dist_ij = sqrt(dist_ij_2);

      SetTauJacobian(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, proj_vector_ij, Area, UnitNormal);

      SetHeatFluxJacobian(Mean_PrimVar, Mean_SecVar, Mean_Eddy_Viscosity,
                          Mean_Thermal_Conductivity, Mean_Cp, proj_vector_ij, Area);

      GetViscousProjJacs(Mean_PrimVar, proj_vector_ij, Area, Proj_Flux_Tensor, Jacobian_i, Jacobian_j, config);
    }

  }

  AD::SetPreaccOut(Proj_Flux_Tensor, nVar);
  AD::EndPreacc();

  return ResidualType<>(Proj_Flux_Tensor, Jacobian_i, Jacobian_j);

}
