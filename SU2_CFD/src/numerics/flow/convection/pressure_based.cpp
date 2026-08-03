/*!
 * \file centered.cpp
 * \brief Implementations of centered schemes.
 * \author F. Palacios, T. Economon
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
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"




CPBConvection_Base::CPBConvection_Base(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  dynamic_grid = config->GetDynamic_Grid();
  
  AdvectedVelocity = new su2double [nVar];
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

  DensityInc_i = V_i[nDim+2];  DensityInc_j = V_j[nDim+2];
  MeanDensity = 0.5*(DensityInc_i + DensityInc_j);

  /*--- Find mass flux (note that edgevelocity itself already included grid movement) ---*/

  MassFlux = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    MassFlux += MeanDensity * EdgeVelocity[iDim] * Normal[iDim];

  // if (dynamic_grid) {
  //   su2double ProjGridMassFlux = 0.0;
  //   for (iDim = 0; iDim < nDim; iDim++)
  //     ProjGridMassFlux += 0.5*(DensityInc_i*GridVel_i[iDim]+DensityInc_j*GridVel_j[iDim])*Normal[iDim];
  //   MassFlux -= ProjGridMassFlux;
  // }

  /*--- Find the velocity that is advected ---*/

  ComputeAdvectedVelocity();

  /*--- Find (momentum) flux. ---*/

  for (iVar = 0; iVar < nDim; iVar++) {
    Flux[iVar] = MassFlux * AdvectedVelocity[iVar];
  }

  /*--- Find Jacobian ---*/

  if (implicit)
    ComputeJacobian();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}

void CPBConvection_Central::ComputeAdvectedVelocity() {

  for (iVar = 0; iVar < nDim; iVar++) 
    AdvectedVelocity[iVar] = 0.5 * (V_i[iVar+1] + V_j[iVar+1]);
  
}

void CPBConvection_Central::ComputeJacobian() {
  
  for (iVar = 0; iVar < nDim; iVar++) {
    for (jVar = 0; jVar < nDim; jVar++) {

      su2double dFdMomentum= (AdvectedVelocity[iVar] * Normal[jVar] + 0.5 * MassFlux / MeanDensity * delta[iVar][jVar]);

      Jacobian_i[iVar][jVar] = 0.5*dFdMomentum;
      Jacobian_j[iVar][jVar] = 0.5*dFdMomentum;
    }
  }
}

void CPBConvection_Upwind::ComputeAdvectedVelocity() {

  bool Upw_i = (MassFlux>0);
  for (iVar = 0; iVar < nDim; iVar++) {
    if (Upw_i)
      AdvectedVelocity[iVar] = V_i[iVar+1];
    else
      AdvectedVelocity[iVar] = V_j[iVar+1];
  }
}

void CPBConvection_Upwind::ComputeJacobian() {

  bool Upw_i = (MassFlux>0);
  
  for (iVar = 0; iVar < nDim; iVar++) {
    for (jVar = 0; jVar < nDim; jVar++) {

      su2double dFdMomentum= (AdvectedVelocity[iVar] * Normal[jVar] + MassFlux / MeanDensity * delta[iVar][jVar]);

      Jacobian_i[iVar][jVar] = static_cast<su2double>(Upw_i) * dFdMomentum;
      Jacobian_j[iVar][jVar] = static_cast<su2double>(!Upw_i) * dFdMomentum;
    }
  }
}