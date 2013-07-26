/*!
 * \file numerics_direct_levelset.cpp
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

CUpwLin_LevelSet::CUpwLin_LevelSet(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
  
}

CUpwLin_LevelSet::~CUpwLin_LevelSet(void) {
  
	delete [] Velocity_i;
	delete [] Velocity_j;
  
}

void CUpwLin_LevelSet::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                   double **val_JacobianMeanFlow_i, double **val_JacobianMeanFlow_j, CConfig *config) {
	unsigned short iDim; //, jDim;
	double a0, a1, q_ij, Velocity_i[3], Velocity_j[3]; //, dqij_dvi[3], dqij_dvj[3], dabsqij_dvi[3], dabsqij_dvj[3], da0_dvi[3], da0_dvj[3], da1_dvi[3], da1_dvj[3];
  
	q_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = V_i[iDim+1];
		Velocity_j[iDim] = V_i[iDim+1];
		q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
	}
  
	a0 = 0.5*(q_ij+fabs(q_ij)); a1 = 0.5*(q_ij-fabs(q_ij));
  
	val_residual[0] = a0*LevelSetVar_i[0]+a1*LevelSetVar_j[0];
  
	if (implicit) {
		val_Jacobian_i[0][0] = a0;
		val_Jacobian_j[0][0] = a1;
    
    //    for (iDim = 0; iDim < nDim; iDim++) {
    //      dqij_dvi[iDim] = 0.5 * Normal[iDim]/V_i[0];
    //      dqij_dvj[iDim] = 0.5 * Normal[iDim]/V_j[0];
    //      if ( q_ij >= 0.0 ) {
    //        dabsqij_dvi[iDim] = dqij_dvi[iDim];
    //        dabsqij_dvj[iDim] = dqij_dvj[iDim];
    //      }
    //      else {
    //        dabsqij_dvi[iDim] = -dqij_dvi[iDim];
    //        dabsqij_dvj[iDim] = -dqij_dvj[iDim];
    //      }
    //      da0_dvi[iDim] = 0.5 * (dqij_dvi[iDim] + dabsqij_dvi[iDim]);
    //      da1_dvi[iDim] = 0.5 * (dqij_dvi[iDim] - dabsqij_dvi[iDim]);
    //
    //      da0_dvj[iDim] = 0.5 * (dqij_dvj[iDim] + dabsqij_dvj[iDim]);
    //      da1_dvj[iDim] = 0.5 * (dqij_dvj[iDim] - dabsqij_dvj[iDim]);
    //    }
    //
    //    for (iDim = 0; iDim < nDim+1; iDim++) {
    //      for (jDim = 0; jDim < nDim+1; jDim++) {
    //        val_JacobianMeanFlow_i[iDim][jDim] = 0.0;
    //        val_JacobianMeanFlow_j[iDim][jDim] = 0.0;
    //      }
    //    }
    //
    //    val_JacobianMeanFlow_i[0][0] = 0.0; val_JacobianMeanFlow_j[0][0] = 0.0;
    //    for (iDim = 0; iDim < nDim; iDim++) {
    //      val_JacobianMeanFlow_i[0][iDim+1] = da0_dvi[iDim]*LevelSetVar_i[0]+da1_dvi[iDim]*LevelSetVar_j[0];
    //      val_JacobianMeanFlow_j[0][iDim+1] = da0_dvj[iDim]*LevelSetVar_i[0]+da1_dvj[iDim]*LevelSetVar_j[0];
    //    }
    
  }
}

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

