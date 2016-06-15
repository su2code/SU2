/*!
 * \file numerics_adjoint_levelset.cpp
 * \brief This file contains all the convective term discretization.
 * \author F. Palacios
 * \version 4.2.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/numerics_structure.hpp"
#include <limits>

CUpwLin_AdjLevelSet::CUpwLin_AdjLevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_AdjLevelSet() == EULER_IMPLICIT);
	Velocity_i = new su2double [nDim];
	Velocity_j = new su2double [nDim];
  
}

CUpwLin_AdjLevelSet::~CUpwLin_AdjLevelSet(void) {
  
	delete [] Velocity_i;
	delete [] Velocity_j;
  
}

void CUpwLin_AdjLevelSet::ComputeResidual (su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii,
                                       su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) {
  
  unsigned short iDim;
	su2double proj_conv_flux_i = 0.0, proj_conv_flux_j = 0.0, proj_conv_flux_ij = 0.0;
  
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/DensityInc_i;
		Velocity_j[iDim] = U_j[iDim+1]/DensityInc_j;
		proj_conv_flux_i -= Velocity_i[iDim]*Normal[iDim]; // projection of convective flux at iPoint
		proj_conv_flux_j -= Velocity_j[iDim]*Normal[iDim]; // projection of convective flux at jPoint
	}
	proj_conv_flux_ij = 0.5*fabs(proj_conv_flux_i+proj_conv_flux_j); // projection of average convective flux
  
	val_residual_i[0] = 0.5*( proj_conv_flux_i*(LevelSetVar_i[0]+LevelSetVar_j[0])-proj_conv_flux_ij*(LevelSetVar_j[0]-LevelSetVar_i[0]));
	val_residual_j[0] = 0.5*(-proj_conv_flux_j*(LevelSetVar_j[0]+LevelSetVar_i[0])-proj_conv_flux_ij*(LevelSetVar_i[0]-LevelSetVar_j[0]));
  
	if (implicit) {
		val_Jacobian_ii[0][0] = 0.5*( proj_conv_flux_i+proj_conv_flux_ij);
		val_Jacobian_ij[0][0] = 0.5*( proj_conv_flux_i-proj_conv_flux_ij);
		val_Jacobian_ji[0][0] = 0.5*(-proj_conv_flux_j-proj_conv_flux_ij);
		val_Jacobian_jj[0][0] = 0.5*(-proj_conv_flux_j+proj_conv_flux_ij);
	}
  
}

CSourcePieceWise_AdjLevelSet::CSourcePieceWise_AdjLevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
}

CSourcePieceWise_AdjLevelSet::~CSourcePieceWise_AdjLevelSet(void) { }

void CSourcePieceWise_AdjLevelSet::ComputeResidual(su2double *val_residual, CConfig *config) {}
