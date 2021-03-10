/*!
 * \file scalar_convection.cpp
 * \brief This file contains the numerical methods for scalar transport eqns.
 * \author T. Economon, D. Mayer, N. Beishuizen
 * \version 7.1.0 "Blackbird"
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

#include "../../../include/numerics/flamelet/scalar_convection.hpp"

CUpwtransportedScalar::CUpwtransportedScalar(unsigned short val_nDim,
                       unsigned short val_nVar,
                       const CConfig* config) :
  CNumerics(val_nDim, val_nVar, config),
  // FIXME dan: this has to work for turbulence scalars and others. 
  //implicit(config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT)
  implicit(config->GetKind_TimeIntScheme_Scalar() == EULER_IMPLICIT),
  incompressible(config->GetKind_Regime() == INCOMPRESSIBLE),
  dynamic_grid(config->GetDynamic_Grid())
{
  Flux = new su2double [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CUpwtransportedScalar::~CUpwtransportedScalar(void) {

  delete [] Flux;
  if (Jacobian_i != nullptr) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete [] Jacobian_i[iVar];
      delete [] Jacobian_j[iVar];
    }
    delete [] Jacobian_i;
    delete [] Jacobian_j;
  }
}

CNumerics::ResidualType<> CUpwtransportedScalar::ComputeResidual(const CConfig* config) {

  unsigned short iDim;

  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);

  // FIXME daniel: this has to work for TurbVar and for flamelet scalars
  // AD::SetPreaccIn(TurbVar_i, nVar);  AD::SetPreaccIn(TurbVar_j, nVar);

  AD::SetPreaccIn(Scalar_i, nVar);  AD::SetPreaccIn(Scalar_j, nVar);

  if (dynamic_grid) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }

  ExtraADPreaccIn();

  Density_i = V_i[nDim+2];
  Density_j = V_j[nDim+2];

  q_ij = 0.0;
  if (dynamic_grid) {
    for (iDim = 0; iDim < nDim; iDim++) {
      su2double Velocity_i = V_i[iDim+1] - GridVel_i[iDim];
      su2double Velocity_j = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i+Velocity_j)*Normal[iDim];
    }
  }
  else {
    for (iDim = 0; iDim < nDim; iDim++) {
      su2double Velocity_i = V_i[iDim+1];
      su2double Velocity_j = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i+Velocity_j)*Normal[iDim];
    }
  }

  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));

  FinishResidualCalc(config);

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

CUpwtransportedScalar_General::CUpwtransportedScalar_General(unsigned short val_nDim,
                                       unsigned short val_nVar,
                                       CConfig *config)
: CUpwtransportedScalar(val_nDim, val_nVar, config) { }

CUpwtransportedScalar_General::~CUpwtransportedScalar_General(void) { }

void CUpwtransportedScalar_General::ExtraADPreaccIn() {
  AD::SetPreaccIn(Scalar_i, nVar);  AD::SetPreaccIn(Scalar_j, nVar);
  AD::SetPreaccIn(V_i, nDim+3); AD::SetPreaccIn(V_j, nDim+3);
  
}

void CUpwtransportedScalar_General::FinishResidualCalc(su2double *val_residual,
                                            su2double **val_Jacobian_i,
                                            su2double **val_Jacobian_j,
                                            CConfig *config) {
  
  unsigned short iVar, jVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = (a0*Density_i*Scalar_i[iVar] +
                          a1*Density_j*Scalar_j[iVar]);
    if (implicit) {
      for (jVar = 0; jVar < nVar; jVar++) {
        if (iVar == jVar) {
          val_Jacobian_i[iVar][jVar] = a0*Density_i;
          val_Jacobian_j[iVar][jVar] = a1*Density_j;
        } else {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
  }
}

void CUpwtransportedScalar_General::FinishResidualCalc(const CConfig* config) {

  unsigned short iVar, jVar;

  for (iVar = 0; iVar < nVar; iVar++) {

    Flux[iVar] = a0 * Density_i * Scalar_i[iVar]
               + a1 * Density_j * Scalar_j[iVar];

    if (implicit) {
      for (jVar = 0; jVar < nVar; jVar++) {
        if (iVar == jVar) {
          Jacobian_i[iVar][jVar] = a0*Density_i;
          Jacobian_j[iVar][jVar] = a1*Density_j;
        } else {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
  }
}