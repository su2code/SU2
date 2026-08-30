/*!
 * \file pressure_based.cpp
 * \brief Implementations of fluxes for pressure-based solvers.
 * \author T. Aalbers
 * \version 8.5.0 "Harrier"
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

#include "../../../../include/numerics/flow/convection/pressure_based.hpp"

CPBConvection_Base::CPBConvection_Base(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  dynamic_grid = config->GetDynamic_Grid();
  energy = config->GetEnergy_Equation();
  variable_density = (config->GetVariable_Density_Model());
  
  AdvectedVelocity = new su2double [MAXNDIM];
  Flux = new su2double [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CPBConvection_Base::~CPBConvection_Base(void) {

  delete [] AdvectedVelocity;
  delete [] Flux;
    
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }

  delete [] Jacobian_i;
  delete [] Jacobian_j;
  
}

CNumerics::ResidualType<> CPBConvection_Base::ComputeResidual(const CConfig *config) {

  /*--- Primitive variables at point i and j ---*/

  Pressure_i    = V_i[0];             Pressure_j    = V_j[0];
  DensityInc_i  = V_i[nDim+2];        DensityInc_j  = V_j[nDim+2];
  Enthalpy_i    = V_i[nDim+3];        Enthalpy_j    = V_j[nDim+3];
  MeanPressure = 0.5 * (Pressure_i + Pressure_j);
  MeanDensity = 0.5 * (DensityInc_i + DensityInc_j);

  /*--- Find the velocity that is advected ---*/

  ComputeAdvectedQuantities();

  /*--- Set the flux vector. ---*/

  Flux[0] = MassFlux;
  for (iDim = 0; iDim < nDim; ++iDim) 
    Flux[1+iDim] = MassFlux * AdvectedVelocity[iDim];
  Flux[nDim+1] = MassFlux * AdvectedEnthalpy;

  /*--- Find Jacobian ---*/

  if (implicit) {

    for (jVar = 0; jVar < nVar; jVar++)
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[iVar][jVar] = 0.0;
        Jacobian_j[iVar][jVar] = 0.0;
      }

    /*--- We need the derivative of the equation of state to build the
    preconditioning matrix. For now, the only option is the ideal gas
    law, but in the future, dRhodT should be in the fluid model. ---*/

    dRhodh_i = 0.0; dRhodh_j = 0.0;
    if (variable_density) {
      Temperature_i = V_i[nDim+1];        Temperature_j = V_j[nDim+1];
      Cp_i          = V_i[nDim+8];        Cp_j          = V_j[nDim+8];

      dRhodh_i = -DensityInc_i / (Temperature_i * Cp_i);
      dRhodh_j = -DensityInc_j / (Temperature_j * Cp_j);
    }

    ComputeJacobianWeights();
    ComputeJacobian(DensityInc_i, &V_i[1], Enthalpy_i, dRhodh_i, weight_jacobian_i, Jacobian_i);
    ComputeJacobian(DensityInc_j, &V_j[1], Enthalpy_j, dRhodh_j, weight_jacobian_j, Jacobian_j);

  }

  /*--- Remove energy contributions if we aren't solving the energy equation. ---*/

  if (!energy) {
    Flux[nDim+1] = 0.0;
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[iVar][nDim+1] = 0.0;
        Jacobian_j[iVar][nDim+1] = 0.0;

        Jacobian_i[nDim+1][iVar] = 0.0;
        Jacobian_j[nDim+1][iVar] = 0.0;
      }
    }
  }

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}


void CPBConvection_Base::ComputeJacobian(su2double val_density, const su2double *val_velocity,
                                         su2double val_enthalpy, su2double val_dRhodh,
                                         su2double val_scale, su2double **val_Proj_Jac_Tensor) {

  su2double proj_vel = MassFlux / MeanDensity;

  /*--- Fill continuity parts ---*/

  val_Proj_Jac_Tensor[0][0] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    val_Proj_Jac_Tensor[0][iDim+1] = val_scale*(Normal[iDim] * (val_density));
  }

  /*--- Fill momentum parts ---*/

  for (jDim = 0; jDim < nDim; jDim++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      val_Proj_Jac_Tensor[iDim+1][jDim+1] = val_scale * val_density * (val_velocity[iDim] * Normal[jDim] + proj_vel * delta[iDim][jDim]);
    }
  }

  /*--- Fill enthalpy parts ---*/

  val_Proj_Jac_Tensor[nDim+1][0] = 0.0;
  val_Proj_Jac_Tensor[0][nDim+1] = val_scale * ((val_dRhodh) * proj_vel);
  val_Proj_Jac_Tensor[nDim+1][nDim+1] = val_scale*(((val_enthalpy)*(val_dRhodh) + (val_density))*proj_vel);
  for (iDim = 0; iDim < nDim; iDim++) {
    val_Proj_Jac_Tensor[nDim+1][iDim+1] = val_scale * ((val_enthalpy) * Normal[iDim] * (val_density));
    val_Proj_Jac_Tensor[iDim+1][nDim+1] = val_scale*((val_dRhodh) * val_velocity[iDim] * proj_vel);
  }

}

void CPBConvection_Central::ComputeAdvectedQuantities() {

  for (iDim = 0; iDim < nDim; iDim++) 
    AdvectedVelocity[iDim] = 0.5 * (V_i[iDim+1] + V_j[iDim+1]);

  AdvectedEnthalpy = 0.5 * (Enthalpy_i + Enthalpy_j);
  
}

void CPBConvection_Central::ComputeJacobianWeights() {

  weight_jacobian_i = weight_jacobian_j = 0.5;

}

void CPBConvection_Upwind::ComputeAdvectedQuantities() {

  bool Upw_i = (MassFlux>0);

  for (iDim = 0; iDim < nDim; iDim++)
    AdvectedVelocity[iDim] = Upw_i ? V_i[iDim+1] : V_j[iDim+1];

  AdvectedEnthalpy = (Upw_i) ? Enthalpy_i : Enthalpy_j;

}

void CPBConvection_Upwind::ComputeJacobianWeights() {

  bool Upw_i = (MassFlux>0);

  weight_jacobian_i = static_cast<su2double>(Upw_i);
  weight_jacobian_j = static_cast<su2double>(!Upw_i);

}
