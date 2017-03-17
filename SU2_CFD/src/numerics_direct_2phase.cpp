/*!
 * \file numerics_direct_2phase.cpp
 * \brief This file contains all the convective term discretization.
 * \author F. Palacios, A. Bueno
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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



CUpw_2phaseHill::CUpw_2phaseHill(unsigned short val_nDim, unsigned short val_nVar,
                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_2phase() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  
}

CUpw_2phaseHill::~CUpw_2phaseHill(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  
}

void CUpw_2phaseHill::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(Two_phaseVar_i,2);
  AD::SetPreaccIn(Two_phaseVar_j,2);
  AD::SetPreaccIn(Normal, nDim);
  if (grid_movement) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+2);
    AD::SetPreaccIn(V_j, nDim+2);

    Density_i = V_i[nDim+1];
    Density_j = V_j[nDim+1];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+3);
    AD::SetPreaccIn(V_j, nDim+3);

    Density_i = V_i[nDim+2];
    Density_j = V_j[nDim+2];
  }
  
  q_ij = 0.0;
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  
  val_residual[0] = a0*Density_i*Two_phaseVar_i[0]+a1*Density_j*Two_phaseVar_j[0];
  val_residual[1] = a0*Density_i*Two_phaseVar_i[1]+a1*Density_j*Two_phaseVar_j[1];
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;    val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[1][0] = 0.0;    val_Jacobian_i[1][1] = a0;
    
    val_Jacobian_j[0][0] = a1;    val_Jacobian_j[0][1] = 0.0;
    val_Jacobian_j[1][0] = 0.0;    val_Jacobian_j[1][1] = a1;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}

void  CUpw_2phaseHill::ComputeResidual_2phase(su2double *Primitive, su2double *Residual, su2double **Jacobian_i) {

	unsigned short iVar;

	// for euler equations, mass transfer

	q_ij = 0;

	if (grid_movement) {
	    for (iDim = 0; iDim < nDim; iDim++) {
	      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
	      q_ij += Velocity_i[iDim]*Velocity_i[iDim];
	    }
	  }
	  else {
	    for (iDim = 0; iDim < nDim; iDim++) {
	      Velocity_i[iDim] = V_i[iDim+1];
	      q_ij += Velocity_i[iDim]*Velocity_i[iDim];
	    }
	  }

	Residual[0] = Primitive[sizeof Primitive -4] * Volume;

	for (iDim=0; iDim< nDim; iDim++) {
		Residual  [iDim+1] = Primitive[sizeof Primitive-4] * Primitive[iDim + 1] * Volume;
	}

	Residual[nDim] = Primitive[sizeof Primitive-4] * (Primitive[sizeof Primitive-3] + 0.5*q_ij) * Volume;

	for (iVar=0; iVar< nVar; iVar++) {
		for (iVar=0; iVar< nVar; iVar++) {
			Jacobian_i[iVar][iVar] = 0;
		}
	}

}

void CUpw_2phaseHill::ComputeResidual(su2double *Residual, su2double **Jacobian_i, CConfig *config) {
}





CUpw_2phaseQMOM::CUpw_2phaseQMOM(unsigned short val_nDim, unsigned short val_nVar,
                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit        = (config->GetKind_TimeIntScheme_2phase() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();

  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];

}

CUpw_2phaseQMOM::~CUpw_2phaseQMOM(void) {

  delete [] Velocity_i;
  delete [] Velocity_j;

}

void CUpw_2phaseQMOM::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
}

void  CUpw_2phaseQMOM::ComputeResidual_2phase(su2double *Primitive, su2double *Residual, su2double **Jacobian_i) {
}

void CUpw_2phaseQMOM::ComputeResidual(su2double *Residual, su2double **Jacobian_i, CConfig *config) {
}



