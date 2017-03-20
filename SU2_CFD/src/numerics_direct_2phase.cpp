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



CUpw_2phaseHill_Rus::CUpw_2phaseHill_Rus(unsigned short val_nDim, unsigned short val_nVar,
                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_2phase() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  
}

CUpw_2phaseHill_Rus::~CUpw_2phaseHill_Rus(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  
}

void CUpw_2phaseHill_Rus::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
	// compute the advective fluxes for the moments equations

  unsigned short iVar, jVar;

  AD::StartPreacc();
  AD::SetPreaccIn(Two_phaseVar_i,2);
  AD::SetPreaccIn(Two_phaseVar_j,2);
  AD::SetPreaccIn(Normal, nDim);

  if (grid_movement) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  
  AD::SetPreaccIn(V_i, nDim+3);
  AD::SetPreaccIn(V_j, nDim+3);

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
  
  a0 = 0.5*(1.0+fabs(q_ij));
  a1 = 0.5*(1.0-fabs(q_ij));
  
  val_residual[0] = a0*Two_phaseVar_i[0]+a1*Two_phaseVar_j[0];
  val_residual[1] = a0*Two_phaseVar_i[1]+a1*Two_phaseVar_j[1];
  val_residual[2] = a0*Two_phaseVar_i[2]+a1*Two_phaseVar_j[2];
  val_residual[3] = a0*Two_phaseVar_i[3]+a1*Two_phaseVar_j[3];
  
  if (implicit) {

	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			val_Jacobian_i[iVar][jVar] = 0.0;
		}
	}

    val_Jacobian_i[0][0] = a0;   val_Jacobian_j[0][0] = a1;
    val_Jacobian_i[1][1] = a0;   val_Jacobian_j[1][1] = a1;
    val_Jacobian_i[2][2] = a0;   val_Jacobian_j[2][2] = a1;
    val_Jacobian_i[3][3] = a0;   val_Jacobian_j[3][3] = a1;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}

void  CUpw_2phaseHill_Rus::ComputeResidual_HeatMassTransfer(su2double *Primitive, su2double *Residual, su2double **Jacobian_i) {

	unsigned short iVar;

	// compute the source terms for the governing equations of the continuum phase

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

void CUpw_2phaseHill_Rus::ComputeResidual(su2double *Residual, su2double **Jacobian_i, su2double *val_liquid_i, CConfig *config) {

	unsigned short iVar, jVar;

	su2double Density_mixture_i, Critical_radius, Nucleation_rate, Growth_rate;
	su2double P, T, rho, h, k, mu;

	// compute the source terms for the moments equations

	Density_mixture_i = val_liquid_i[8];
	Critical_radius = val_liquid_i[6];

	// retrieve the thermodynamic properties of the continuum phase primitive

	T   = V_i[0];
	P   = V_i[nDim+1];
	rho = V_i[nDim+2];
	h   = V_i[nDim+3];
	mu  = V_i[nDim+5];
	k   = V_i[nDim+7];

	// compute the nucleation rate the growth rate

	SetNucleationRate(P, T, rho, h, k, mu, val_liquid_i);
	SetGrowthRate(P, T, rho, h, k, mu, val_liquid_i);

	Nucleation_rate = GetNucleationRate();
	Growth_rate = GetGrowthRate();

	// compute the source terms

	Residual[0] = Density_mixture_i*Nucleation_rate;
	Residual[1] = Density_mixture_i*Nucleation_rate*Critical_radius + Density_mixture_i*Growth_rate*Two_phaseVar_i[0];
	Residual[2] = Density_mixture_i*Nucleation_rate*pow(Critical_radius,2) + 2.0*Density_mixture_i*Growth_rate*Two_phaseVar_i[1];
	Residual[3] = Density_mixture_i*Nucleation_rate*pow(Critical_radius,3) + 3.0*Density_mixture_i*Growth_rate*Two_phaseVar_i[2];

	for (iVar=0; iVar< nVar; iVar++) {
		Residual[iVar] = Residual[iVar]*Volume;
	}

	// compute the Jacobians of the source terms

	if (implicit) {

		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_i[iVar][jVar] = 0.0;
			}
		}

		Jacobian_i[1][0] = Density_mixture_i*Growth_rate*Volume;
		Jacobian_i[2][1] = 2.0*Density_mixture_i*Growth_rate*Volume;
		Jacobian_i[3][2] = 3.0*Density_mixture_i*Growth_rate*Volume;

	}

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



