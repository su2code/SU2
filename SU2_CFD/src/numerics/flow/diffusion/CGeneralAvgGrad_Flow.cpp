/*!
 * \file CGeneralAvgGrad_Flow.cpp
 * \brief Implementation of numerics class CGeneralAvgGrad_Flow.
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

#include "../../../../include/numerics/flow/diffusion/CGeneralAvgGrad_Flow.hpp"

CGeneralAvgGrad_Flow::CGeneralAvgGrad_Flow(unsigned short val_nDim,
                                           unsigned short val_nVar,
                                           bool val_correct_grad,
                                           CConfig *config)
    : CAvgGrad_Base(val_nDim, val_nVar, val_nDim+4, val_correct_grad, config) {

  Mean_SecVar = new su2double [2];

}

CGeneralAvgGrad_Flow::~CGeneralAvgGrad_Flow(void) {

  delete [] Mean_SecVar;

}

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

void CGeneralAvgGrad_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

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

  for (iVar = 0; iVar < nPrimVar; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
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

  /* --- If using UQ methodology, set Reynolds Stress tensor and perform perturbation--- */

  if (using_uq){
    SetReynoldsStressMatrix(Mean_turb_ke);
    SetPerturbedRSM(Mean_turb_ke, config);
  }

  /*--- Get projected flux tensor ---*/

  SetStressTensor(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke,
         Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);

  SetHeatFluxVector(Mean_GradPrimVar, Mean_Laminar_Viscosity,
                    Mean_Eddy_Viscosity, Mean_Thermal_Conductivity, Mean_Cp);

  GetViscousProjFlux(Mean_PrimVar, Normal);

  /*--- Update viscous residual ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];

  /*--- Compute the implicit part ---*/

  if (implicit) {

    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    } else {
      const su2double dist_ij = sqrt(dist_ij_2);
      SetTauJacobian(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
                     dist_ij, UnitNormal);
      SetHeatFluxJacobian(Mean_PrimVar, Mean_SecVar, Mean_Eddy_Viscosity,
                          Mean_Thermal_Conductivity, Mean_Cp, dist_ij);
      GetViscousProjJacs(Mean_PrimVar, Area, Proj_Flux_Tensor,
                         val_Jacobian_i, val_Jacobian_j);
    }

  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}
