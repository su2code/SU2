/*!
 * \file numerics_direct_turbulent.cpp
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

CUpwSca_TurbSA::CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                               CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement  = config->GetGrid_Movement();
	incompressible = config->GetIncompressible();
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
}

CUpwSca_TurbSA::~CUpwSca_TurbSA(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
}

void CUpwSca_TurbSA::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//
  
	//SU2_CPP2C START CUpwSca_TurbSA::ComputeResidual
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *U_i *U_j *TurbVar_i *TurbVar_j
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE *Normal
	//SU2_CPP2C CALL_LIST END
  
	//SU2_CPP2C DEFINE nDim
  
	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS INT SCALAR iDim
	//SU2_CPP2C VARS DOUBLE SCALAR Density_i Density_j q_ij a0 a1
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nDim Velocity_i Velocity_j
	//SU2_CPP2C DECL_LIST END
  
	//SU2_CPP2C COMMENT START
	if (incompressible) {
		Density_i = DensityInc_i;
		Density_j = DensityInc_j;
	}
	else {
    
		//SU2_CPP2C COMMENT END
		Density_i = U_i[0];
		Density_j = U_j[0];
		//SU2_CPP2C COMMENT START
	}
	//SU2_CPP2C COMMENT END
  
	q_ij = 0;
	//SU2_CPP2C COMMENT START
	if (rotating_frame) {
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i - RotVel_i[iDim];
			Velocity_j[iDim] = U_j[iDim+1]/Density_j - RotVel_j[iDim];
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
	} else if (grid_movement) {
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i - GridVel_i[iDim];
			Velocity_j[iDim] = U_j[iDim+1]/Density_j - GridVel_j[iDim];
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
	} else {
    //SU2_CPP2C COMMENT END
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i;
			Velocity_j[iDim] = U_j[iDim+1]/Density_j;
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
    //SU2_CPP2C COMMENT START
	}
	//SU2_CPP2C COMMENT END
  
	a0 = 0.5*(q_ij+fabs(q_ij));
	a1 = 0.5*(q_ij-fabs(q_ij));
	val_residual[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];
  
	//SU2_CPP2C COMMENT START
	if (implicit) {
		val_Jacobian_i[0][0] = a0;
		val_Jacobian_j[0][0] = a1;
	}
	//SU2_CPP2C COMMENT END
  
	//SU2_CPP2C END CUpwSca_TurbSA::ComputeResidual
  
}

CUpwSca_TurbSST::CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement  = config->GetGrid_Movement();
	incompressible = config->GetIncompressible();
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
}

CUpwSca_TurbSST::~CUpwSca_TurbSST(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
}

void CUpwSca_TurbSST::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	if (incompressible) {
		Density_i = DensityInc_i;
		Density_j = DensityInc_j;
	}
	else {
		Density_i = U_i[0];
		Density_j = U_j[0];
	}
  
	q_ij = 0;
	if (rotating_frame) {
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i - RotVel_i[iDim];
			Velocity_j[iDim] = U_j[iDim+1]/Density_j - RotVel_j[iDim];
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
	} else if (grid_movement) {
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i - GridVel_i[iDim];
			Velocity_j[iDim] = U_j[iDim+1]/Density_j - GridVel_j[iDim];
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
	} else {
		for (iDim = 0; iDim < nDim; iDim++) {
			Velocity_i[iDim] = U_i[iDim+1]/Density_i;
			Velocity_j[iDim] = U_j[iDim+1]/Density_j;
			q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
		}
	}
  
	a0 = 0.5*(q_ij+fabs(q_ij));
	a1 = 0.5*(q_ij-fabs(q_ij));
  
	val_residual[0] = a0*Density_i*TurbVar_i[0]+a1*Density_j*TurbVar_j[0];
	val_residual[1] = a0*Density_i*TurbVar_i[1]+a1*Density_j*TurbVar_j[1];
  
	if (implicit) {
		val_Jacobian_i[0][0] = a0;		val_Jacobian_i[0][1] = 0.0;
		val_Jacobian_i[1][0] = 0.0;		val_Jacobian_i[1][1] = a0;
    
		val_Jacobian_j[0][0] = a1;		val_Jacobian_j[0][1] = 0.0;
		val_Jacobian_j[1][0] = 0.0;		val_Jacobian_j[1][1] = a1;
	}
}

CAvgGrad_TurbSA::CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;
  
	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	incompressible = config->GetIncompressible();
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	sigma = 2./3.;
  
	Edge_Vector = new double [nDim];
	Proj_Mean_GradTurbVar_Kappa = new double [nVar];
	Proj_Mean_GradTurbVar_Edge = new double [nVar];
	Mean_GradTurbVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradTurbVar[iVar] = new double [nDim];
}

CAvgGrad_TurbSA::~CAvgGrad_TurbSA(void) {
	unsigned short iVar;
  
	delete [] Edge_Vector;
	delete [] Proj_Mean_GradTurbVar_Kappa;
	delete [] Proj_Mean_GradTurbVar_Edge;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradTurbVar[iVar];
	delete [] Mean_GradTurbVar;
}

void CAvgGrad_TurbSA::ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {
  
	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//
  
	//SU2_CPP2C START CAvgGrad_TurbSA::ComputeResidual
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *U_i *U_j Laminar_Viscosity_i Laminar_Viscosity_j *TurbVar_i *TurbVar_j **TurbVar_Grad_i **TurbVar_Grad_j
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE sigma *Coord_i *Coord_j *Normal
	//SU2_CPP2C CALL_LIST END
  
	//SU2_CPP2C DEFINE nDim nVar
  
	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS INT SCALAR iDim iVar
	//SU2_CPP2C VARS DOUBLE SCALAR Density_i Density_j nu_i nu_j nu_e
	//SU2_CPP2C VARS DOUBLE SCALAR dist_ij_2 proj_vector_ij
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nDim Edge_Vector
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nVar Proj_Mean_GradTurbVar_Kappa Proj_Mean_GradTurbVar_Edge
	//SU2_CPP2C VARS DOUBLE MATRIX SIZE=nVar SIZE=nDim Mean_GradTurbVar
	//SU2_CPP2C DECL_LIST END
  
	//SU2_CPP2C COMMENT START
	if (incompressible) {
		Density_i = DensityInc_i;
		Density_j = DensityInc_j;
	}
	else {
		//SU2_CPP2C COMMENT END
		Density_i = U_i[0];
		Density_j = U_j[0];
		//SU2_CPP2C COMMENT START
	}
	//SU2_CPP2C COMMENT END
  
	//SU2_CPP2C COMMENT START
	switch (config->GetKind_Turb_Model()) {
    case SA :
      //SU2_CPP2C COMMENT END
      
      /*--- Compute mean effective viscosity ---*/
      nu_i = Laminar_Viscosity_i/Density_i;
      nu_j = Laminar_Viscosity_j/Density_j;
      nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
      
      /*--- Compute vector going from iPoint to jPoint ---*/
      dist_ij_2 = 0; proj_vector_ij = 0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
        dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
        proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
      }
      proj_vector_ij = proj_vector_ij/dist_ij_2; // to normalize vectors
      
      /*--- Mean gradient approximation ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
        Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
          Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
        }
      }
      
      val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;
      //SU2_CPP2C COMMENT START
      /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
      if (implicit) {
        Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-nu_e*proj_vector_ij)/sigma;
        Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+nu_e*proj_vector_ij)/sigma;
      }
      //SU2_CPP2C COMMENT END
      
      //SU2_CPP2C COMMENT START
      break;
      
	}
	//SU2_CPP2C COMMENT END
  
	//SU2_CPP2C END CAvgGrad_TurbSA::ComputeResidual
  
}

CAvgGrad_TurbSST::CAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;
  
	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	incompressible = config->GetIncompressible();
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	sigma = 2./3.;
  
	Edge_Vector = new double [nDim];
	Proj_Mean_GradTurbVar_Kappa = new double [nVar];
	Proj_Mean_GradTurbVar_Edge = new double [nVar];
	Mean_GradTurbVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradTurbVar[iVar] = new double [nDim];
}

CAvgGrad_TurbSST::~CAvgGrad_TurbSST(void) {
	unsigned short iVar;
  
	delete [] Edge_Vector;
	delete [] Proj_Mean_GradTurbVar_Kappa;
	delete [] Proj_Mean_GradTurbVar_Edge;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradTurbVar[iVar];
	delete [] Mean_GradTurbVar;
}

void CAvgGrad_TurbSST::ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {
  
	/*--- Compute mean effective viscosity ---*/
	diff_i = Laminar_Viscosity_i + Eddy_Viscosity_i;
	diff_j = Laminar_Viscosity_j + Eddy_Viscosity_j;
	diff_e = 0.5*(diff_i + diff_j);
  
	/*--- Compute vector going from iPoint to jPoint ---*/
	dist_ij_2 = 0; proj_vector_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
		proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
	}
	proj_vector_ij = proj_vector_ij/dist_ij_2; // to normalize vectors
  
	/*--- Mean gradient approximation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
		Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
			Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
		}
	}
  
	val_residual[0] = diff_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;
  
	/*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
	if (implicit) {
		Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-diff_e*proj_vector_ij)/sigma;
		Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+diff_e*proj_vector_ij)/sigma;
	}
}

CAvgGradCorrected_TurbSA::CAvgGradCorrected_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;
  
	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	incompressible = config->GetIncompressible();
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	sigma = 2./3.;
  
	Edge_Vector = new double [nDim];
	Proj_Mean_GradTurbVar_Kappa = new double [nVar];
	Proj_Mean_GradTurbVar_Edge = new double [nVar];
	Proj_Mean_GradTurbVar_Corrected = new double [nVar];
	Mean_GradTurbVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradTurbVar[iVar] = new double [nDim];
}

CAvgGradCorrected_TurbSA::~CAvgGradCorrected_TurbSA(void) {
	unsigned short iVar;
  
	delete [] Edge_Vector;
	delete [] Proj_Mean_GradTurbVar_Kappa;
	delete [] Proj_Mean_GradTurbVar_Edge;
	delete [] Proj_Mean_GradTurbVar_Corrected;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradTurbVar[iVar];
	delete [] Mean_GradTurbVar;
}

void CAvgGradCorrected_TurbSA::ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {
  
	if (incompressible) {
		Density_i = DensityInc_i;
		Density_j = DensityInc_j;
	}
	else {
		Density_i = U_i[0];
		Density_j = U_j[0];
	}
  
	/*--- Compute mean effective viscosity ---*/
	nu_i = Laminar_Viscosity_i/Density_i;
	nu_j = Laminar_Viscosity_j/Density_j;
	nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
	/*--- Compute vector going from iPoint to jPoint ---*/
	dist_ij_2 = 0; proj_vector_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
		proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
	}
	proj_vector_ij = proj_vector_ij/dist_ij_2;
  
	/*--- Mean gradient approximation. Projection of the mean gradient
	 in the direction of the edge ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
		Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
			Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
			Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
		}
		Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
		Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
	}
  
	val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;
  
	/*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
	if (implicit) {
		Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
		Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
	}
  
}

CAvgGradCorrected_TurbSST::CAvgGradCorrected_TurbSST(unsigned short val_nDim, unsigned short val_nVar, double *constants,
                                                     CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;
  
	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	incompressible = config->GetIncompressible();
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	sigma_k1  = constants[0];
	sigma_om1 = constants[2];
	sigma_k2  = constants[1];
	sigma_om2 = constants[3];
  
	Edge_Vector = new double [nDim];
	Proj_Mean_GradTurbVar_Normal = new double [nVar];
	Proj_Mean_GradTurbVar_Edge = new double [nVar];
	Proj_Mean_GradTurbVar_Corrected = new double [nVar];
	Mean_GradTurbVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradTurbVar[iVar] = new double [nDim];
}

CAvgGradCorrected_TurbSST::~CAvgGradCorrected_TurbSST(void) {
	unsigned short iVar;
  
	delete [] Edge_Vector;
	delete [] Proj_Mean_GradTurbVar_Normal;
	delete [] Proj_Mean_GradTurbVar_Edge;
	delete [] Proj_Mean_GradTurbVar_Corrected;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradTurbVar[iVar];
	delete [] Mean_GradTurbVar;
}

void CAvgGradCorrected_TurbSST::ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {
  
	double sigma_kine_i, sigma_kine_j, sigma_omega_i, sigma_omega_j;
	double diff_i_kine, diff_i_omega, diff_j_kine, diff_j_omega;
  
	/*--- Compute the blended constant for the viscous terms ---*/
	sigma_kine_i  = F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2;
	sigma_kine_j  = F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2;
	sigma_omega_i = F1_i*sigma_om1 + (1.0 - F1_i)*sigma_om2;
	sigma_omega_j = F1_j*sigma_om1 + (1.0 - F1_j)*sigma_om2;
  
	/*--- Compute mean effective viscosity ---*/
	diff_i_kine  = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
	diff_j_kine  = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
	diff_i_omega = Laminar_Viscosity_i + sigma_omega_i*Eddy_Viscosity_i;
	diff_j_omega = Laminar_Viscosity_j + sigma_omega_j*Eddy_Viscosity_j;
  
	diff_kine  = 0.5*(diff_i_kine + diff_j_kine);    // Could instead use weighted average!
	diff_omega = 0.5*(diff_i_omega + diff_j_omega);
  
	/*--- Compute vector going from iPoint to jPoint ---*/
	dist_ij_2 = 0; proj_vector_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
		proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
	}
	proj_vector_ij = proj_vector_ij/dist_ij_2;
  
	/*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
		Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
			Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
			Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
		}
		Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
		Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
	}
  
	val_residual[0] = diff_kine*Proj_Mean_GradTurbVar_Corrected[0];
	val_residual[1] = diff_omega*Proj_Mean_GradTurbVar_Corrected[1];
  
	/*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
	if (implicit) {
		Jacobian_i[0][0] = -diff_kine*proj_vector_ij/U_i[0];		Jacobian_i[0][1] = 0.0;
		Jacobian_i[1][0] = 0.0;									    Jacobian_i[1][1] = -diff_omega*proj_vector_ij/U_i[0];
    
		Jacobian_j[0][0] = diff_kine*proj_vector_ij/U_j[0]; 		Jacobian_j[0][1] = 0.0;
		Jacobian_j[1][0] = 0.0;									    Jacobian_j[1][1] = diff_omega*proj_vector_ij/U_j[0];
	}
  
}

CSourcePieceWise_TurbSA::CSourcePieceWise_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	incompressible = config->GetIncompressible();
	//transition     = (config->GetKind_Trans_Model() == LM);
  transition = false; // Debugging, -AA
  rotating_frame = config->GetRotating_Frame();
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	/*--- Spalart-Allmaras closure constants ---*/
	cv1_3 = pow(7.1,3.0);
	k2 = pow(0.41,2.0);
	cb1 = 0.1355;
	cw2 = 0.3;
	cw3_6 = pow(2.0,6.0);
	sigma = 2./3.;
	cb2 = 0.622;
	cw1 = cb1/k2+(1+cb2)/sigma;
  
	/*--- LM transition model constants ---*/
	beta = 0.5;
	s1   = 2.0;
  
}

CSourcePieceWise_TurbSA::~CSourcePieceWise_TurbSA(void) { }

void CSourcePieceWise_TurbSA::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//
  
	//SU2_CPP2C START CSourcePieceWise_TurbSA::ComputeResidual
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *U_i **PrimVar_Grad_i Laminar_Viscosity_i *TurbVar_i **TurbVar_Grad_i
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE dist_i cv1_3 k2 cb1 cb2 cw1 cw2 cw3_6 Volume sigma TURB_EPS
	//SU2_CPP2C CALL_LIST END
  
	//SU2_CPP2C DEFINE nDim
  
	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS INT SCALAR iDim
	//SU2_CPP2C VARS DOUBLE SCALAR Density_i DivVelocity Vorticity dist_i_2
	//SU2_CPP2C VARS DOUBLE SCALAR nu Ji Ji_2 Ji_3 fv1 fv2 Omega Shat
	//SU2_CPP2C VARS DOUBLE SCALAR r g g_6 glim fw norm2_Grad
	//SU2_CPP2C DECL_LIST END
  
	//SU2_CPP2C COMMENT START
	if (incompressible) Density_i = DensityInc_i;
	else {
		//SU2_CPP2C COMMENT END
		Density_i = U_i[0];
		//SU2_CPP2C COMMENT START
	}
	//SU2_CPP2C COMMENT END
  
	val_residual[0] = 0.0;
  
	//SU2_CPP2C COMMENT START
  val_Jacobian_i[0][0] = 0.0;
	//SU2_CPP2C COMMENT END
  
	/*--- Computation of vorticity ---*/
	Vorticity = (PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1])*(PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1]);
	if (nDim == 3) Vorticity += ( (PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2])*(PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2]) +
                               (PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) );
	Omega = max(sqrt(Vorticity), 1.0e-10);
	dist_i = max(dist_i, 1.0e-10);
  
  /*--- Rotational correction term ---*/
  if (rotating_frame) {
    div = PrimVar_Grad_i[1][0] + PrimVar_Grad_i[2][1];
    if (nDim == 3) div += PrimVar_Grad_i[3][2];
    StrainMag = 0.0;
    // add diagonals
    StrainMag += pow(PrimVar_Grad_i[1][0] - 1.0/3.0*div,2.0);
    StrainMag += pow(PrimVar_Grad_i[2][1] - 1.0/3.0*div,2.0);
    if (nDim == 3) StrainMag += pow(PrimVar_Grad_i[3][2] - 1.0/3.0*div,2.0);
    // add off diagonals
    StrainMag += 2.0*pow(0.5*(PrimVar_Grad_i[1][1]+PrimVar_Grad_i[2][0]),2.0);
    if (nDim == 3) {
      StrainMag += 2.0*pow(0.5*(PrimVar_Grad_i[1][2]+PrimVar_Grad_i[3][0]),2.0);
      StrainMag += 2.0*pow(0.5*(PrimVar_Grad_i[2][2]+PrimVar_Grad_i[3][1]),2.0);
    }
    StrainMag = sqrt(2.0*StrainMag);
    Omega += 2.0*min(0.0,StrainMag-Omega);
  }
  
	if (dist_i > 0.0) {
    
		/*--- Production term ---*/
		dist_i_2 = dist_i*dist_i;
		nu = Laminar_Viscosity_i/Density_i;
		Ji = TurbVar_i[0]/nu;
		Ji_2 = Ji*Ji;
		Ji_3 = Ji_2*Ji;
		fv1 = Ji_3/(Ji_3+cv1_3);
		fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    S = Omega;
    inv_k2_d2 = 1.0/(k2*dist_i_2);
    
		Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
    inv_Shat = 1.0/max(Shat, 1.0e-10);
    
    /*--- Production term ---*/
		if (!transition) val_residual[0] += cb1*Shat*TurbVar_i[0]*Volume;
    else val_residual[0] += cb1*Shat*TurbVar_i[0]*Volume*intermittency;
    
		/*--- Destruction term ---*/
		r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
		g = r + cw2*(pow(r,6.0)-r);
		g_6 =	pow(g,6.0);
		glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
		fw = g*glim;
    
		if (!transition) val_residual[0] -= cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
		else val_residual[0] -= cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume*min(max(intermittency,0.1),1.0);
    
		/*--- Diffusion term ---*/
		norm2_Grad = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];
		val_residual[0] += cb2/sigma*norm2_Grad*Volume;
    
		//SU2_CPP2C COMMENT START
    
		/*--- Implicit part ---*/
    
    /*--- Production term ---*/
    dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
    dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
    if ( Shat <= 1.0e-10 ) dShat = 0.0;
    else dShat = (fv2+TurbVar_i[0]*dfv2)*inv_k2_d2;
    val_Jacobian_i[0][0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;
    
    /*--- Destruction term ---*/
    dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
    if (r == 10.0) dr = 0.0;
    dg = dr*(1.+cw2*(6.*pow(r,5.)-1.));
    dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
    val_Jacobian_i[0][0] -= cw1*(dfw*TurbVar_i[0] +	2.*fw)*TurbVar_i[0]/dist_i_2*Volume;
    
		//SU2_CPP2C COMMENT END
	}
	//SU2_CPP2C COMMENT START
	//SU2_CPP2C COMMENT END
  
	//SU2_CPP2C END CSourcePieceWise_TurbSA::ComputeResidual
  
}

CSourcePieceWise_TurbSST::CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar, double *constants,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	incompressible = config->GetIncompressible();
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	/*--- Closure constants ---*/
	beta_star     = constants[6];
	sigma_omega_1 = constants[2];
	sigma_omega_2 = constants[3];
	beta_1        = constants[4];
	beta_2        = constants[5];
	alfa_1        = constants[8];
	alfa_2        = constants[9];
	a1            = constants[7];
}

CSourcePieceWise_TurbSST::~CSourcePieceWise_TurbSST(void) { }

void CSourcePieceWise_TurbSST::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//
  
	//SU2_CPP2C START CSourcePieceWise_TurbSST::ComputeResidual
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *val_U_i
	//SU2_CPP2C OUTVARS *val_laminar_viscosity_i
	//SU2_CPP2C VARS DOUBLE Temperature_Ref Viscosity_Ref, Gamma_Minus_One
	//SU2_CPP2C CALL_LIST END
  
	//SU2_CPP2C DEFINE nDim NONE SA SST
  
	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C DECL_LIST END
  
	unsigned short iDim;
	double alfa_blended, beta_blended;
	double diverg, pk, pw, zeta;
  
	val_residual[0] = 0.0;
	val_residual[1] = 0.0;
	//SU2_CPP2C COMMENT START
  val_Jacobian_i[0][0] = 0.0;		val_Jacobian_i[0][1] = 0.0;
  val_Jacobian_i[1][0] = 0.0;		val_Jacobian_i[1][1] = 0.0;
	//SU2_CPP2C COMMENT END
  
	/*--- Computation of blended constants for the source terms---*/
	alfa_blended = F1_i*alfa_1 + (1.0 - F1_i)*alfa_2;
	beta_blended = F1_i*beta_1 + (1.0 - F1_i)*beta_2;
  
	if (dist_i > 0.0) {
		/*--- Production ---*/
		diverg = 0;
		for (iDim = 0; iDim < nDim; iDim++)
			diverg += PrimVar_Grad_i[iDim+1][iDim];
    
		pk = Eddy_Viscosity_i*StrainMag*StrainMag - 2.0/3.0*U_i[0]*TurbVar_i[0]*diverg;
		pk = min(pk,20.0*beta_star*U_i[0]*TurbVar_i[1]*TurbVar_i[0]);
		pk = max(pk,0.0);
    
		zeta = max(TurbVar_i[1],StrainMag*F2_i/a1);
		pw = StrainMag*StrainMag - 2.0/3.0*zeta*diverg;
		pw = max(pw,0.0);
    
		val_residual[0] += pk*Volume;
		val_residual[1] += alfa_blended*U_i[0]*pw*Volume;
    
		/*--- Dissipation ---*/
		val_residual[0] -= beta_star*U_i[0]*TurbVar_i[1]*TurbVar_i[0]*Volume;
		val_residual[1] -= beta_blended*U_i[0]*TurbVar_i[1]*TurbVar_i[1]*Volume;
    
		/*--- Cross diffusion ---*/
		val_residual[1] += (1.0 - F1_i)*CDkw*Volume;
    
		//SU2_CPP2C COMMENT START
		/*--- Implicit part ---*/
    val_Jacobian_i[0][0] = -beta_star*TurbVar_i[1]*Volume;		val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[1][0] = 0.0;                               val_Jacobian_i[1][1] = -2.0*beta_blended*TurbVar_i[1]*Volume;
		//SU2_CPP2C COMMENT END
	}
  
	//SU2_CPP2C END CSourcePieceWise_TurbSST::ComputeResidual
  
}
