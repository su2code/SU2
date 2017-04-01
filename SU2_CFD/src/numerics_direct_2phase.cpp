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

  unsigned short iVar, jVar, *Liquid_vec;
  su2double qi, qj;

  AD::StartPreacc();
  AD::SetPreaccIn(Two_phaseVar_i,nVar);
  AD::SetPreaccIn(Two_phaseVar_j,nVar);

  AD::SetPreaccIn(Primitive_Liquid,10);

  AD::SetPreaccIn(Normal, nDim);

  if (grid_movement) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  
  AD::SetPreaccIn(V_i, nDim+9);
  AD::SetPreaccIn(V_j, nDim+9);


//  if (Primitive_Liquid[4] > V_i[0]) {
	  qi = 0.0;
	  qj = 0.0;
	  if (grid_movement) {
		for (iDim = 0; iDim < nDim; iDim++) {
		  Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
		  Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
		  qi += Velocity_i[iDim]*Normal[iDim];
		  qj += Velocity_j[iDim]*Normal[iDim];
		}
	  }
	  else {
		for (iDim = 0; iDim < nDim; iDim++) {
		  Velocity_i[iDim] = V_i[iDim+1];
		  Velocity_j[iDim] = V_j[iDim+1];
		  qi += Velocity_i[iDim]*Normal[iDim];
		  qj += Velocity_j[iDim]*Normal[iDim];
		}
	  }


	  for (iVar=0; iVar<nVar; iVar++) {
		  val_residual[iVar] = 0.5* qi * Two_phaseVar_i[iVar] + 0.5 * qj* Two_phaseVar_j[iVar];
		  val_residual[iVar] = val_residual[iVar] - 0.5 * max(qi, qj) * (Two_phaseVar_j[iVar] - Two_phaseVar_i[iVar]);
	  }

//  } else {
//
//	  val_residual[0] = 0;
//	  val_residual[1] = 0;
//	  val_residual[2] = 0;
//	  val_residual[3] = 0;
//
//  }
  
  if (implicit) {

	for (iVar = 0; iVar < nVar; iVar++) {
		for (jVar = 0; jVar < nVar; jVar++) {
			val_Jacobian_i[iVar][jVar] = 0.0;
			val_Jacobian_j[iVar][jVar] = 0.0;
		}
	}

//	if (Primitive_Liquid[4] > V_i[0]) {
    val_Jacobian_i[0][0] = 0.5 * (qi + max(qi, qj));   val_Jacobian_j[0][0] = 0.5 * (qi - max(qi, qj));
    val_Jacobian_i[1][1] = 0.5 * (qi + max(qi, qj));   val_Jacobian_j[1][1] = 0.5 * (qi - max(qi, qj));
    val_Jacobian_i[2][2] = 0.5 * (qi + max(qi, qj));   val_Jacobian_j[2][2] = 0.5 * (qi - max(qi, qj));
    val_Jacobian_i[3][3] = 0.5 * (qi + max(qi, qj));   val_Jacobian_j[3][3] = 0.5 * (qi - max(qi, qj));
//	}
  }


  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}


CSourcePieceWise_Hill::CSourcePieceWise_Hill(unsigned short val_nDim, unsigned short val_nVar,
                                              CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);


}

CSourcePieceWise_Hill::~CSourcePieceWise_Hill(void) { }


void CSourcePieceWise_Hill::ComputeResidual(su2double *val_Residual, su2double **val_Jacobian_i, su2double *val_liquid_i, su2double *V, CConfig *config) {

	unsigned short iVar, jVar;

	su2double Density_mixture, Critical_radius, Nucleation_rate, Growth_rate;
	su2double P, T, rho, h, k, mu;

	// compute the source terms for the moments equations
	Critical_radius = val_liquid_i[6];
	Density_mixture = val_liquid_i[8];

	// retrieve the thermodynamic properties of the continuum phase primitive

	T   = V_i[0];
	P   = V_i[nDim+1];
	rho = V_i[nDim+2];
	h   = V_i[nDim+3];
	mu  = V_i[nDim+5];
	k   = V_i[nDim+7];

	// compute the nucleation rate and the growth rate

	if (val_liquid_i[4] > T) {

		Nucleation_rate = GetNucleation_Rate(P, T, rho, h, k, mu, val_liquid_i);
		Growth_rate = GetGrowth_Rate(P, T, rho, h, k, mu, val_liquid_i);

		// store G for source term euler
		val_liquid_i[9] = Growth_rate;

		// compute the source terms
		val_Residual[0] = Density_mixture * Nucleation_rate;
		val_Residual[1] = Density_mixture * Nucleation_rate*Critical_radius        +     Growth_rate*Two_phaseVar_i[0];
		val_Residual[2] = Density_mixture * Nucleation_rate*pow(Critical_radius,2) + 2.0*Growth_rate*Two_phaseVar_i[1];
		val_Residual[3] = Density_mixture * Nucleation_rate*pow(Critical_radius,3) + 3.0*Growth_rate*Two_phaseVar_i[2];


		for (iVar=0; iVar<nVar; iVar++) {
			val_Residual[iVar] = val_Residual[iVar] * Volume;
		}

	}	else 	{

		for (iVar=0; iVar<nVar; iVar++) {
			val_Residual[iVar] = 0;
		}
	}

	// compute the Jacobians of the source terms

	if (implicit) {

		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_i[iVar][jVar] = 0.0;
			}
		}

/*		if (val_liquid_i[4] > T) {
			val_Jacobian_i[1][0] =     Growth_rate* Volume;
			val_Jacobian_i[2][1] = 2.0*Growth_rate* Volume;
			val_Jacobian_i[3][2] = 3.0*Growth_rate* Volume;
		}
*/
	}
}


/*


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


*/


