/*!
 * \file flow_diffusion.cpp
 * \brief Implementation of numerics classes for discretization
 *        of viscous fluxes in fluid flow problems.
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

#include "../../../include/numerics/flow/flow_diffusion.hpp"
#include "../../../../Common/include/toolboxes/geometry_toolbox.hpp"

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

  TauWall_i = 0; TauWall_j = 0;

  Mean_PrimVar = new su2double [nPrimVar];

  Mean_GradPrimVar = new su2double* [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++)
    Mean_GradPrimVar[iVar] = new su2double [nDim];

  Proj_Mean_GradPrimVar_Edge = new su2double[val_nPrimVar];

  tau_jacobian_i = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    tau_jacobian_i[iDim] = new su2double [nVar];
  }

  heat_flux_jac_i = new su2double[nVar];

  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }

}

CAvgGrad_Base::~CAvgGrad_Base() {

  delete [] Mean_PrimVar;

  if (Mean_GradPrimVar != nullptr) {
    for (unsigned short iVar = 0; iVar < nPrimVar; iVar++)
      delete [] Mean_GradPrimVar[iVar];
    delete [] Mean_GradPrimVar;
  }

  delete [] Proj_Mean_GradPrimVar_Edge;

  if (tau_jacobian_i != nullptr) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      delete [] tau_jacobian_i[iDim];
    }
    delete [] tau_jacobian_i;
  }

  delete [] heat_flux_jac_i;

  if (Jacobian_i != nullptr) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete [] Jacobian_i[iVar];
      delete [] Jacobian_j[iVar];
    }
    delete [] Jacobian_i;
    delete [] Jacobian_j;
  }

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

  const su2double Density = val_primvar[nDim+2];

  /* --- If UQ methodology is used, use the perturbed Reynolds stress tensor
   * for the turbulent part of tau. Otherwise both the laminar and turbulent
   * parts of tau can be computed with the total viscosity. --- */

  if (sstParsedOptions.uq) {
    // laminar part
    ComputeStressTensor(nDim, tau, val_gradprimvar+1, val_laminar_viscosity);
    // add turbulent part which was perturbed
    for (unsigned short iDim = 0 ; iDim < nDim; iDim++)
      for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
        tau[iDim][jDim] += (-Density) * MeanPerturbedRSM[iDim][jDim];
  } else {
    const su2double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
    // turb_ke is not considered in the stress tensor, see #797
    ComputeStressTensor(nDim, tau, val_gradprimvar+1, total_viscosity, Density, su2double(0.0));
  }
}

void CAvgGrad_Base::AddTauWall(const su2double *UnitNormal,
                               const su2double TauWall) {

  /*--- Compute the wall shear stress as the magnitude of the
   tangential projection of the shear stress tensor. ---*/

  su2double TauTangent[MAXNDIM];
  GeometryToolbox::TangentProjection(nDim, tau, UnitNormal, TauTangent);

  su2double WallShearStress = GeometryToolbox::Norm(nDim, TauTangent);
  su2double Scale = TauWall / WallShearStress;

  /*--- Scale the stress tensor by the ratio of the wall shear stress
   (from wall functions) to the one computed above. ---*/

  for (auto iDim = 0u; iDim < nDim; iDim++)
    for (auto jDim = 0u; jDim < nDim; jDim++)
      tau[iDim][jDim] *= Scale;
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

  su2double Flux_Tensor[5][3];

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

CAvgGrad_Flow::CAvgGrad_Flow(unsigned short val_nDim,
                             unsigned short val_nVar,
                             bool val_correct_grad,
                             const CConfig* config)
    : CAvgGrad_Base(val_nDim, val_nVar, val_nDim+3, val_correct_grad, config) { }

CNumerics::ResidualType<> CAvgGrad_Flow::ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+9);   AD::SetPreaccIn(V_j, nDim+9);
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);
  AD::SetPreaccIn(PrimVar_Grad_j, nDim+1, nDim);
  AD::SetPreaccIn(turb_ke_i); AD::SetPreaccIn(turb_ke_j);
  AD::SetPreaccIn(TauWall_i); AD::SetPreaccIn(TauWall_j);
  AD::SetPreaccIn(Normal, nDim);

  unsigned short iVar, jVar, iDim;

  /*--- Normalized normal vector ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  PrimVar_i = V_i;
  PrimVar_j = V_j;

  for (iVar = 0; iVar < nPrimVar; iVar++) {
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
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

  /*--- Projection of the mean gradient in the direction of the edge ---*/

  if (correct_gradient && dist_ij_2 != 0.0) {
    CorrectGradient(Mean_GradPrimVar, PrimVar_i, PrimVar_j, Edge_Vector,
                    dist_ij_2, nDim+1);
  }

  /*--- Wall shear stress values (wall functions) only used if present for one but not both points (xor) ---*/

  const int scale = (TauWall_i > 0.0) ^ (TauWall_j > 0.0);
  Mean_TauWall = (max(TauWall_i,0.0) + max(TauWall_j,0.0)) * scale;

  /*--- If using UQ methodology, set Reynolds Stress tensor and perform perturbation ---*/

  if (sstParsedOptions.uq) {
    ComputePerturbedRSM(nDim, Eig_Val_Comp, uq_permute, uq_delta_b, uq_urlx,
                        Mean_GradPrimVar+1, Mean_PrimVar[nDim+2], Mean_Eddy_Viscosity,
                        Mean_turb_ke, MeanPerturbedRSM);
  }

  /*--- Get projected flux tensor (viscous residual) ---*/

  SetStressTensor(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke,
                  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  if (config->GetSAParsedOptions().qcr2000) AddQCR(nDim, &Mean_GradPrimVar[1], tau);
  if (Mean_TauWall > 0) AddTauWall(UnitNormal, Mean_TauWall);

  SetHeatFluxVector(Mean_GradPrimVar, Mean_Laminar_Viscosity,
                    Mean_Eddy_Viscosity);

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
      const su2double dist_ij = sqrt(dist_ij_2);
      SetTauJacobian(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, dist_ij, UnitNormal);

      SetHeatFluxJacobian(Mean_PrimVar, Mean_Laminar_Viscosity,
                          Mean_Eddy_Viscosity, dist_ij, UnitNormal);

      GetViscousProjJacs(Mean_PrimVar, Area, Proj_Flux_Tensor, Jacobian_i, Jacobian_j);
    }

  }

  AD::SetPreaccOut(Proj_Flux_Tensor, nVar);
  AD::EndPreacc();

  return ResidualType<>(Proj_Flux_Tensor, Jacobian_i, Jacobian_j);

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
                                        const su2double val_dist_ij,
                                        const su2double *val_normal) {
  su2double sqvel = 0.0;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
  }

  const su2double Density = val_Mean_PrimVar[nDim+2];
  const su2double Pressure = val_Mean_PrimVar[nDim+1];
  const su2double phi = Gamma_Minus_One/Density;

  /*--- R times partial derivatives of temp. ---*/

  const su2double R_dTdu0 = -Pressure/(Density*Density) + 0.5*sqvel*phi;
  const su2double R_dTdu1 = -phi*val_Mean_PrimVar[1];
  const su2double R_dTdu2 = -phi*val_Mean_PrimVar[2];

  const su2double heat_flux_factor = val_laminar_viscosity/Prandtl_Lam + val_eddy_viscosity/Prandtl_Turb;
  const su2double cpoR = Gamma/Gamma_Minus_One; // cp over R
  const su2double conductivity_over_Rd = cpoR*heat_flux_factor/val_dist_ij;

  heat_flux_jac_i[0] = conductivity_over_Rd * R_dTdu0;
  heat_flux_jac_i[1] = conductivity_over_Rd * R_dTdu1;
  heat_flux_jac_i[2] = conductivity_over_Rd * R_dTdu2;

  if (nDim == 2) {

    const su2double R_dTdu3 = phi;
    heat_flux_jac_i[3] = conductivity_over_Rd * R_dTdu3;

  } else {

    const su2double R_dTdu3 = -phi*val_Mean_PrimVar[3];
    const su2double R_dTdu4 = phi;
    heat_flux_jac_i[3] = conductivity_over_Rd * R_dTdu3;
    heat_flux_jac_i[4] = conductivity_over_Rd * R_dTdu4;

  }
}

CAvgGradInc_Flow::CAvgGradInc_Flow(unsigned short val_nDim,
                                   unsigned short val_nVar,
                                   bool val_correct_grad, const CConfig* config)
    : CAvgGrad_Base(val_nDim, val_nVar, val_nDim+3, val_correct_grad, config) {

  energy   = config->GetEnergy_Equation();

}

CNumerics::ResidualType<> CAvgGradInc_Flow::ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+9);   AD::SetPreaccIn(V_j, nDim+9);
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(PrimVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(PrimVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(turb_ke_i); AD::SetPreaccIn(turb_ke_j);
  AD::SetPreaccIn(TauWall_i); AD::SetPreaccIn(TauWall_j);
  AD::SetPreaccIn(Normal, nDim);

  unsigned short iVar, jVar, iDim;

  /*--- Normalized normal vector ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  PrimVar_i = V_i;
  PrimVar_j = V_j;

  for (iVar = 0; iVar < nPrimVar; iVar++) {
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
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
                    dist_ij_2, nVar);
  }

  /*--- Wall shear stress values (wall functions) only used if present for one but not both points (xor) ---*/

  const int scale = (TauWall_i > 0.0) ^ (TauWall_j > 0.0);
  Mean_TauWall = (max(TauWall_i,0.0) + max(TauWall_j,0.0)) * scale;

  /*--- If using UQ methodology, set Reynolds Stress tensor and perform perturbation ---*/

  if (sstParsedOptions.uq) {
    ComputePerturbedRSM(nDim, Eig_Val_Comp, uq_permute, uq_delta_b, uq_urlx,
                        Mean_GradPrimVar+1, Mean_PrimVar[nDim+2], Mean_Eddy_Viscosity,
                        Mean_turb_ke, MeanPerturbedRSM);
  }

  /*--- Get projected flux tensor (viscous residual) ---*/
  SetStressTensor(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke,
                  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  if (config->GetSAParsedOptions().qcr2000) AddQCR(nDim, &Mean_GradPrimVar[1], tau);
  if (Mean_TauWall > 0) AddTauWall(UnitNormal, Mean_TauWall);

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

  su2double Flux_Tensor[5][3];

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

    val_Proj_Jac_Tensor_i[1][0] = val_dS*tau_jacobian_i[0][0];
    val_Proj_Jac_Tensor_i[1][1] = val_dS*tau_jacobian_i[0][1];
    val_Proj_Jac_Tensor_i[1][2] = val_dS*tau_jacobian_i[0][2];
    val_Proj_Jac_Tensor_i[1][3] = val_dS*tau_jacobian_i[0][3];

    val_Proj_Jac_Tensor_i[2][0] = val_dS*tau_jacobian_i[1][0];
    val_Proj_Jac_Tensor_i[2][1] = val_dS*tau_jacobian_i[1][1];
    val_Proj_Jac_Tensor_i[2][2] = val_dS*tau_jacobian_i[1][2];
    val_Proj_Jac_Tensor_i[2][3] = val_dS*tau_jacobian_i[1][3];

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
                                               const su2double val_dist_ij) {
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
  su2double factor2 = total_conductivity/val_dist_ij;

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

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+9);   AD::SetPreaccIn(V_j, nDim+9);
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(S_i, 4); AD::SetPreaccIn(S_j, 4);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);
  AD::SetPreaccIn(PrimVar_Grad_j, nDim+1, nDim);
  AD::SetPreaccIn(turb_ke_i); AD::SetPreaccIn(turb_ke_j);
  AD::SetPreaccIn(TauWall_i); AD::SetPreaccIn(TauWall_j);
  AD::SetPreaccIn(Normal, nDim);

  unsigned short iVar, jVar, iDim;

  /*--- Normalized normal vector ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Mean primitive variables ---*/

  PrimVar_i = V_i;
  PrimVar_j = V_j;

  for (iVar = 0; iVar < nPrimVar; iVar++) {
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
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
                    dist_ij_2, nDim+1);
  }

  /*--- Wall shear stress values (wall functions) only used if present for one but not both points (xor) ---*/

  const int scale = (TauWall_i > 0.0) ^ (TauWall_j > 0.0);
  Mean_TauWall = (max(TauWall_i,0.0) + max(TauWall_j,0.0)) * scale;

  /*--- If using UQ methodology, set Reynolds Stress tensor and perform perturbation ---*/

  if (sstParsedOptions.uq) {
    ComputePerturbedRSM(nDim, Eig_Val_Comp, uq_permute, uq_delta_b, uq_urlx,
                        Mean_GradPrimVar+1, Mean_PrimVar[nDim+2], Mean_Eddy_Viscosity,
                        Mean_turb_ke, MeanPerturbedRSM);
  }

  /*--- Get projected flux tensor (viscous residual) ---*/

  SetStressTensor(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke,
                  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  if (config->GetSAParsedOptions().qcr2000) AddQCR(nDim, &Mean_GradPrimVar[1], tau);
  if (Mean_TauWall > 0) AddTauWall(UnitNormal, Mean_TauWall);

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
      const su2double dist_ij = sqrt(dist_ij_2);

      SetTauJacobian(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, dist_ij, UnitNormal);

      SetHeatFluxJacobian(Mean_PrimVar, Mean_SecVar, Mean_Eddy_Viscosity,
                          Mean_Thermal_Conductivity, Mean_Cp, dist_ij);

      GetViscousProjJacs(Mean_PrimVar, Area, Proj_Flux_Tensor, Jacobian_i, Jacobian_j);
    }

  }

  AD::SetPreaccOut(Proj_Flux_Tensor, nVar);
  AD::EndPreacc();

  return ResidualType<>(Proj_Flux_Tensor, Jacobian_i, Jacobian_j);

}
