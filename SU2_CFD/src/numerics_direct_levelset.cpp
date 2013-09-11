/*!
 * \file numerics_direct_levelset.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.7
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/numerics_structure.hpp"
#include <limits>

CSourcePieceWise_FreeSurface::CSourcePieceWise_FreeSurface(unsigned short val_nDim, unsigned short val_nVar,
                                                           CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
}

CSourcePieceWise_FreeSurface::~CSourcePieceWise_FreeSurface(void) { }

void CSourcePieceWise_FreeSurface::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	unsigned short iVar;
  
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
	for (iVar = 0; iVar < nVar; iVar++)
		val_residual[iVar] = 0.0;
  
	if (nDim == 2) {
		val_residual[1] = -Volume * U_i[0] * AuxVar_Grad_i[0];
		val_residual[2] = -Volume * U_i[0] * AuxVar_Grad_i[1];
    
		if (implicit) {
			val_Jacobian_i[0][0] = 0.0;													val_Jacobian_i[0][1] = 0.0;		val_Jacobian_i[0][2] = 0.0;
			val_Jacobian_i[1][0] = -Volume*AuxVar_Grad_i[0];		val_Jacobian_i[1][1] = 0.0;		val_Jacobian_i[1][2] = 0.0;
			val_Jacobian_i[2][0] = -Volume*AuxVar_Grad_i[1];		val_Jacobian_i[2][1] = 0.0;		val_Jacobian_i[2][2] = 0.0;
		}
	}
  
  
	/*	if (nDim == 2) {
   val_residual[0] = 0.0; //(Volume/DensityInc_i) * ( ( U_i[1] ) * AuxVar_Grad_i[0] + ( U_i[2] ) * AuxVar_Grad_i[1] );
   val_residual[1] = (Volume/DensityInc_i) * ( ( U_i[1]*U_i[1] + U_i[0]/DensityInc_i ) * AuxVar_Grad_i[0] + (U_i[1]*U_i[2]) * AuxVar_Grad_i[1] );
   val_residual[2] = (Volume/DensityInc_i) * ( (U_i[1]*U_i[2]) * AuxVar_Grad_i[0] + ( U_i[2]*U_i[2] + U_i[0]/DensityInc_i ) * AuxVar_Grad_i[1] );
   
   if (implicit) {
   val_Jacobian_i[0][0] = 0.0;
   val_Jacobian_i[0][1] = 0.0; //(Volume/DensityInc_i) * AuxVar_Grad_i[0];
   val_Jacobian_i[0][2] = 0.0; //(Volume/DensityInc_i) * AuxVar_Grad_i[1];
   val_Jacobian_i[1][0] = (Volume/DensityInc_i) * (1.0/DensityInc_i) * AuxVar_Grad_i[0];
   val_Jacobian_i[1][1] = (Volume/DensityInc_i) * (2.0*U_i[1]*AuxVar_Grad_i[0]+U_i[2]*AuxVar_Grad_i[1]);
   val_Jacobian_i[1][2] = (Volume/DensityInc_i) * (U_i[1]*AuxVar_Grad_i[1]);
   val_Jacobian_i[2][0] = (Volume/DensityInc_i) * (1.0/DensityInc_i) * AuxVar_Grad_i[1];
   val_Jacobian_i[2][1] = (Volume/DensityInc_i) * (U_i[2]*AuxVar_Grad_i[0]);
   val_Jacobian_i[2][2] = (Volume/DensityInc_i) * (U_i[1]*AuxVar_Grad_i[0]+2.0*U_i[2]*AuxVar_Grad_i[1]);
   }
   
   } */
  
}

CSourcePieceWise_LevelSet::CSourcePieceWise_LevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
}

CSourcePieceWise_LevelSet::~CSourcePieceWise_LevelSet(void) { }

void CSourcePieceWise_LevelSet::ComputeResidual(double *val_residual, CConfig *config) {}

