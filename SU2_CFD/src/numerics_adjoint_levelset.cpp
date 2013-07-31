/*!
 * \file numerics_adjoint_levelset.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/numerics_structure.hpp"
#include <limits>

CUpwLin_AdjLevelSet::CUpwLin_AdjLevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_AdjLevelSet() == EULER_IMPLICIT);
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
  
}

CUpwLin_AdjLevelSet::~CUpwLin_AdjLevelSet(void) {
  
	delete [] Velocity_i;
	delete [] Velocity_j;
  
}

void CUpwLin_AdjLevelSet::ComputeResidual (double *val_residual_i, double *val_residual_j, double **val_Jacobian_ii,
                                       double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config)  {
  
  unsigned short iDim;
	double proj_conv_flux_i = 0.0, proj_conv_flux_j = 0.0, proj_conv_flux_ij = 0.0;
  
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

void CSourcePieceWise_AdjLevelSet::ComputeResidual(double *val_residual, CConfig *config) {}

