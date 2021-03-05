/*!
 * \file scalar_sources.cpp
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

#include "../../../include/numerics/flamelet/scalar_sources.hpp"

CSourcePieceWise_Scalar::CSourcePieceWise_Scalar(unsigned short val_nDim,
                                                 unsigned short val_nVar,
                                                 CConfig *config) :
CNumerics(val_nDim, val_nVar, config) {
  
  implicit       = (config->GetKind_TimeIntScheme_Scalar() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
}

CSourcePieceWise_Scalar::~CSourcePieceWise_Scalar(void) { }

void CSourcePieceWise_Scalar::ComputeResidual(su2double *val_residual,
                                              su2double **val_Jacobian_i,
                                              su2double **val_Jacobian_j,
                                              CConfig *config) {
  
  unsigned short iVar, jVar;
  
  Density_i = V_i[nDim+2];
  
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.0;
    if (implicit) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_i[iVar][jVar] = 0.0;
      }
    }
  }
  
}

CSourceAxisymmetric_Scalar::CSourceAxisymmetric_Scalar(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  energy   = config->GetEnergy_Equation();
  viscous  = config->GetViscous();
  
}

CSourceAxisymmetric_Scalar::~CSourceAxisymmetric_Scalar(void) { }

void CSourceAxisymmetric_Scalar::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config) {
  
  su2double yinv, Velocity_i[3];
  unsigned short iDim, iVar, jVar;
  
  if (Coord_i[1] > EPS) {
    yinv          = 1.0/Coord_i[1];
    Density_i     = V_i[nDim+2];
    
    /*--- Set primitive variables at points iPoint. ---*/
    
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity_i[iDim] = V_i[iDim+1];
    
    /*--- Inviscid component of the source term. ---*/
    
    for (iVar=0; iVar < nVar; iVar++)
      val_residual[iVar] = yinv*Volume*Density_i*scalar_i[iVar]*Velocity_i[1];
    
    if (implicit) {
      
      for (iVar=0; iVar < nVar; iVar++) {
        for (jVar=0; jVar < nVar; jVar++) {
          if (iVar == jVar) Jacobian_i[iVar][jVar] = Velocity_i[1];
          Jacobian_i[iVar][jVar] *= yinv*Volume*Density_i;
        }
      }
      
    }
    
    /*--- Add the viscous terms if necessary. ---*/
    
    if (viscous) {
      
      for (iVar=0; iVar < nVar; iVar++)
        val_residual[iVar] -= Volume*yinv*Diffusion_Coeff_i[iVar]*scalar_grad_i[iVar][1];
      
    }
    
  } else {
    
    for (iVar=0; iVar < nVar; iVar++)
      val_residual[iVar] = 0.0;
    
    if (implicit) {
      for (iVar=0; iVar < nVar; iVar++) {
        for (jVar=0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;
      }
    }
    
  }
  
}

