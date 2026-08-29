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
  Pressure_i = V_i[0]; Pressure_j = V_j[0];
  DensityInc_i = V_i[nDim+2];  DensityInc_j = V_j[nDim+2];
  MeanPressure = 0.5 * (Pressure_i + Pressure_j);
  MeanDensity = 0.5 * (DensityInc_i + DensityInc_j);

  /*--- Find projected velocity (note that the massflux itself already should include grid movement) ---*/

  su2double ProjVelocity = MassFlux / MeanDensity;

  /*--- Find the velocity that is advected ---*/

  ComputeAdvectedVelocity();

  /*--- Set continuity (pressure) flux, only used by bounded scalar or coupled solver. ---*/

  Flux[0] = MassFlux;

  /*--- Set momentum flux. ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    Flux[iDim+1] = ProjVelocity * AdvectedVelocity[iDim];

  /*--- Find Jacobian ---*/

  if (implicit) {

    for (jVar = 0; jVar < nVar; jVar++)
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[iVar][jVar] = 0.0;
        Jacobian_j[iVar][jVar] = 0.0;
      }

    ComputeJacobian();

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

void CPBConvection_Central::ComputeAdvectedVelocity() {

  for (iDim = 0; iDim < nDim; iDim++) 
    AdvectedVelocity[iDim] = 0.5 * (V_i[iDim+1] + V_j[iDim+1]);
  
}

void CPBConvection_Central::ComputeJacobian() {
  
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {

      su2double dFdu = (AdvectedVelocity[iDim] * Normal[jDim] + 0.5 * MassFlux / MeanDensity * delta[iDim][jDim]);

      Jacobian_i[iDim+1][jDim+1] = 0.5 * dFdu;
      Jacobian_j[iDim+1][jDim+1] = 0.5 * dFdu;
    }
  }
}

void CPBConvection_Upwind::ComputeAdvectedVelocity() {

  bool Upw_i = (MassFlux>0);
  for (iDim = 0; iDim < nDim; iDim++) {
    if (Upw_i)
      AdvectedVelocity[iDim] = V_i[iDim+1];
    else
      AdvectedVelocity[iDim] = V_j[iDim+1];
  }

}

void CPBConvection_Upwind::ComputeJacobian() {

  bool Upw_i = (MassFlux>0);
  
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {

      su2double dFdu= (AdvectedVelocity[iDim] * Normal[jDim] + MassFlux / MeanDensity * delta[iDim][jDim]);

      Jacobian_i[iDim+1][jDim+1] = static_cast<su2double>(Upw_i) * dFdu;
      Jacobian_j[iDim+1][jDim+1] = static_cast<su2double>(!Upw_i) * dFdu;
    }
  }
}