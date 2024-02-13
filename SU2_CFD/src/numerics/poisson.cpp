/*!
 * \file numerics_direct_poisson.cpp
 * \brief This file contains all the convective term discretization.
 * \author F. Palacios
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../include/numerics/poisson.hpp"


CAvgGradCorrected_Poisson::CAvgGradCorrected_Poisson(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = false;
  implicit = (config->GetKind_TimeIntScheme_Poisson() == EULER_IMPLICIT);
  direct = false;

  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradPoissonVar_Edge = new su2double [nVar];
  Proj_Mean_GradPoissonVar_Kappa = new su2double [nVar];
  Proj_Mean_GradPoissonVar_Corrected = new su2double [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  Mean_GradPoissonVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nDim];
    Jacobian_j[iVar] = new su2double [nDim];
    Mean_GradPoissonVar[iVar] = new su2double [nDim];
  }

  Poisson_Coeff_i = 1.0;
  Poisson_Coeff_j = 1.0;
  
  Mom_Coeff_i = new su2double [nDim];
  Mom_Coeff_j = new su2double [nDim];

}


CAvgGradCorrected_Poisson::~CAvgGradCorrected_Poisson(void) {

	delete [] Edge_Vector;
	delete [] Proj_Mean_GradPoissonVar_Edge;
	delete [] Proj_Mean_GradPoissonVar_Kappa;
	delete [] Proj_Mean_GradPoissonVar_Corrected;
	delete [] Mom_Coeff_i;
	delete [] Mom_Coeff_j;
	
	for (iVar = 0; iVar < nVar; iVar++) {
       delete [] Mean_GradPoissonVar[iVar];
       delete [] Jacobian_i[iVar];
       delete [] Jacobian_j[iVar];
   }
    delete [] Mean_GradPoissonVar;
    delete [] Jacobian_i;
    delete [] Jacobian_j;
}

CNumerics::ResidualType<> CAvgGradCorrected_Poisson::ComputeResidual(const CConfig *config) {

  su2double Coeff_Mean;
  su2double *MomCoeffxNormal = new su2double[nDim];
  su2double *MomCoeffxNormalCorrected = new su2double[nDim];
  su2double  Mean_GradPoissonVar_Edge[MAXNDIM], GradPoisson[MAXNDIM];
  su2double  Mean_GradPoissonVar_Face[MAXNDIM][MAXNDIM], Num, Den;
   
  Poisson_Coeff_Mean = 1.0;//0.5*(Poisson_Coeff_i + Poisson_Coeff_j);

 /*--- Multiply the normal with the coefficients of the momentum equation
   *    and use it instead of the normal during the projection operation ---*/
   for (iDim = 0; iDim < nDim; iDim++) {
	   Coeff_Mean = 0.5*(Mom_Coeff_i[iDim] + Mom_Coeff_j[iDim]) ;
	   MomCoeffxNormal[iDim] = Coeff_Mean*Normal[iDim];
   }

  /*--- Compute vector going from iPoint to jPoint ---*/
  /*--- Minimum correction approach. ---*/
  /*dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*MomCoeffxNormal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;*/

  /*--- Over relaxed approach. ---*/
  dist_ij_2 = 0; proj_vector_ij = 0;
  Num = 0.0; Den = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    Num += MomCoeffxNormal[iDim]*MomCoeffxNormal[iDim];
    Den += Edge_Vector[iDim]*MomCoeffxNormal[iDim];
  }
  if (Den == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = Num/Den;
  
  /*--- Orthogonal correction approach. ---*/
  /*dist_ij_2 = 0; proj_vector_ij = 0;
  Num = 0.0; Den = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    Num += MomCoeffxNormal[iDim]*MomCoeffxNormal[iDim];
  }
  Den = dist_ij_2;
  if (Den == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = sqrt(Num/Den);*/

  /*--- Correct the face gradient for odd-even decoupling ---*/
  /*--- Steps are 
   * 1. Interpolate the gradient at the face -> du/ds|_in
   * 2. Find the projection of the interpolated gradient on the edge vector -> (du/ds|_in) . e
   * 3. Find the gradient as the difference between the neighboring nodes -> (u_j - u_i)/ds
   * 4. Correct the gradient at the face using du/ds = du/ds|_in + ( (u_j - u_i)/ds - (du/ds|_in . e) ) . e ---*/
  
  /*for (iVar = 0; iVar < nVar; iVar++) {
	  Mean_GradPoissonVar_Edge[iVar] = 0.0;	  
	  
	  GradPoisson[iVar] = (Poissonval_j - Poissonval_i)/sqrt(dist_ij_2);
	  
	  for (iDim = 0; iDim < nDim; iDim++) {
		  Mean_GradPoissonVar[iVar][iDim] = 0.5*(ConsVar_Grad_i[iVar][iDim] + ConsVar_Grad_j[iVar][iDim]);  
		  
		  Mean_GradPoissonVar_Edge[iVar] += Mean_GradPoissonVar[iVar][iDim]*Edge_Vector[iDim]/dist_ij_2;
       }
   }
   
   for (iVar = 0; iVar < nVar; iVar++) {
	   for (iDim = 0; iDim < nDim; iDim++) {
		   Mean_GradPoissonVar_Face[iVar][iDim] = Mean_GradPoissonVar[iVar][iDim] + 
                                        (GradPoisson[iVar] - Mean_GradPoissonVar_Edge[iVar])*Edge_Vector[iDim]/dist_ij_2;
       }
   }*/


  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradPoissonVar_Edge[iVar] = 0.0;
    Proj_Mean_GradPoissonVar_Kappa[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPoissonVar[iVar][iDim] = 0.5*(ConsVar_Grad_i[iVar][iDim] + ConsVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradPoissonVar_Kappa[iVar] += Mean_GradPoissonVar[iVar][iDim]*MomCoeffxNormal[iDim];
      Proj_Mean_GradPoissonVar_Edge[iVar] += Mean_GradPoissonVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradPoissonVar_Corrected[iVar] = Proj_Mean_GradPoissonVar_Kappa[iVar];
    Proj_Mean_GradPoissonVar_Corrected[iVar] -= Proj_Mean_GradPoissonVar_Edge[iVar]*proj_vector_ij -
    (Poissonval_j-Poissonval_i)*proj_vector_ij;
  }
  
  /*--- Jacobians for implicit scheme ---*/

  if (implicit) {
    Jacobian_i[0][0] = -Poisson_Coeff_Mean*proj_vector_ij;
    Jacobian_j[0][0] = Poisson_Coeff_Mean*proj_vector_ij;
  }

  /*AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();*/
  
  delete [] MomCoeffxNormal;
  delete [] MomCoeffxNormalCorrected;
  
  return ResidualType<>(Proj_Mean_GradPoissonVar_Corrected, Jacobian_i, Jacobian_j);
}

CAvgGrad_Poisson::CAvgGrad_Poisson(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = false;
  implicit = (config->GetKind_TimeIntScheme_Poisson() == EULER_IMPLICIT); // TODO: PBFlow: FIXED
  direct = false;

   
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradPoissonVar_Normal = new su2double [nVar];
  Proj_Mean_GradPoissonVar_Corrected = new su2double [nVar];
  Mom_Coeff_i = new su2double [nDim];
  Mom_Coeff_j = new su2double [nDim];
  
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  Mean_GradPoissonVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nDim];
    Jacobian_j[iVar] = new su2double [nDim];
    Mean_GradPoissonVar[iVar] = new su2double [nDim];
  }

  Poisson_Coeff_i = 1.0;
  Poisson_Coeff_j = 1.0;
}


CAvgGrad_Poisson::~CAvgGrad_Poisson(void) {
	
	unsigned int iDim;
	delete [] Edge_Vector;
	delete [] Proj_Mean_GradPoissonVar_Normal;
	delete [] Proj_Mean_GradPoissonVar_Corrected;
	delete [] Mom_Coeff_i;
	delete [] Mom_Coeff_j;

	for (iVar = 0; iVar < nVar; iVar++) {
       delete [] Mean_GradPoissonVar[iVar];
       delete [] Jacobian_i[iVar];
       delete [] Jacobian_j[iVar];
    }
    delete [] Mean_GradPoissonVar;
    delete [] Jacobian_i;
    delete [] Jacobian_j;
      
}

CNumerics::ResidualType<> CAvgGrad_Poisson::ComputeResidual(const CConfig *config) {


  su2double Coeff_Mean, Num, Den;
  su2double *MomCoeffxNormal = new su2double[nDim];
  su2double  Mean_GradPoissonVar_Edge[MAXNDIM], GradPoisson[MAXNDIM];
  su2double  Mean_GradPoissonVar_Face[MAXNDIM][MAXNDIM];
  
  
  
  Poisson_Coeff_Mean = 1.0;//0.5*(Poisson_Coeff_i + Poisson_Coeff_j);

  /*--- Multiply the normal with the coefficients of the momentum equation
   *    and use it instead of the normal during the projection operation ---*/
   for (iDim = 0; iDim < nDim; iDim++) {
	   Coeff_Mean = 0.5*(Mom_Coeff_i[iDim] + Mom_Coeff_j[iDim]) ;
	   MomCoeffxNormal[iDim] = Coeff_Mean*Normal[iDim];
   }

  /*--- Compute vector going from iPoint to jPoint ---*/
  /*--- Minimum correction approach. ---*/
  /*dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*MomCoeffxNormal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;*/

  /*--- Over relaxed approach. ---*/
  dist_ij_2 = 0; proj_vector_ij = 0;
  Num = 0.0; Den = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    Num += MomCoeffxNormal[iDim]*MomCoeffxNormal[iDim];
    Den += Edge_Vector[iDim]*MomCoeffxNormal[iDim];
  }
  if (Den == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = Num/Den;
  
  /*--- Orthogonal correction approach. ---*/
  /*dist_ij_2 = 0; proj_vector_ij = 0;
  Num = 0.0; Den = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    Num += MomCoeffxNormal[iDim]*MomCoeffxNormal[iDim];
  }
  Den = dist_ij_2;
  if (Den == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = sqrt(Num/Den);*/
  
  /*--- Correct the face gradient for odd-even decoupling ---*/
  /*--- Steps are 
   * 1. Interpolate the gradient at the face -> du/ds|_in
   * 2. Find the projection of the interpolated gradient on the edge vector -> (du/ds|_in) . e
   * 3. Find the gradient as the difference between the neighboring nodes -> (u_j - u_i)/ds
   * 4. Correct the gradient at the face using du/ds = du/ds|_in + ( (u_j - u_i)/ds - (du/ds|_in . e) ) . e ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
	  Mean_GradPoissonVar_Edge[iVar] = 0.0;	  
	  
	  GradPoisson[iVar] = (Poissonval_j - Poissonval_i)/sqrt(dist_ij_2);
	  
	  for (iDim = 0; iDim < nDim; iDim++) {
		  Mean_GradPoissonVar[iVar][iDim] = 0.5*(ConsVar_Grad_i[iVar][iDim] + ConsVar_Grad_j[iVar][iDim]);  
		  
		  Mean_GradPoissonVar_Edge[iVar] += Mean_GradPoissonVar[iVar][iDim]*Edge_Vector[iDim]/sqrt(dist_ij_2);
       }
   }
   
   for (iVar = 0; iVar < nVar; iVar++) {
	   for (iDim = 0; iDim < nDim; iDim++) {
		   Mean_GradPoissonVar_Face[iVar][iDim] = Mean_GradPoissonVar[iVar][iDim] + 
                                        (GradPoisson[iVar] - Mean_GradPoissonVar_Edge[iVar])*Edge_Vector[iDim]/sqrt(dist_ij_2);
       }
   }

  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradPoissonVar_Normal[iVar] = 0.0;
    Proj_Mean_GradPoissonVar_Corrected[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Proj_Mean_GradPoissonVar_Normal[iVar] += Mean_GradPoissonVar_Face[iVar][iDim]*MomCoeffxNormal[iDim];
    }
  }
  
  /*--- Jacobians for implicit scheme ---*/
  if (implicit) {
    Jacobian_i[0][0] = -Poisson_Coeff_Mean*proj_vector_ij;
    Jacobian_j[0][0] = Poisson_Coeff_Mean*proj_vector_ij;
  }
  
  delete [] MomCoeffxNormal;
  
  return ResidualType<>(Proj_Mean_GradPoissonVar_Normal, Jacobian_i, Jacobian_j);
}



CPressure_Poisson::CPressure_Poisson(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = true;
     
  Edge_Vector = new su2double [nDim];
  
  Poisson_Coeff_i = 1.0;
  Poisson_Coeff_j = 1.0;
  
  Mom_Coeff_i = new su2double [nDim];
  Mom_Coeff_j = new su2double [nDim];
  
}


CPressure_Poisson::~CPressure_Poisson(void) {
	
	unsigned int iDim;
	delete [] Edge_Vector;
	
	delete Mom_Coeff_i;
	delete Mom_Coeff_j;
}

void CPressure_Poisson::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {

  
  su2double Area;
  su2double UnitNormal[MAXNDIM];

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0; proj_vector_ij = 0;
  Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
    UnitNormal[iDim] = Normal[iDim]/Area;
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;

  
  Poisson_Coeff_Mean = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
     Poisson_Coeff_Mean += 0.5*(Mom_Coeff_i[iDim] + Mom_Coeff_j[iDim]);
   
   
  val_residual[0] = 0.5*Poisson_Coeff_Mean*(Poissonval_i + Poissonval_j);
    
  if (implicit){
	Jacobian_i[0][0] = Poisson_Coeff_Mean;
	Jacobian_j[0][0] = -Poisson_Coeff_Mean;
  }
}


CSource_PoissonFVM::CSource_PoissonFVM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSource_PoissonFVM::~CSource_PoissonFVM(void) { }

void CSource_PoissonFVM::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {
  
 if (config->GetKind_Incomp_System()==ENUM_INCOMP_SYSTEM::PRESSURE_BASED) 
    val_residual[0] = Source_Term;
 else 
   val_residual[0] = Source_Term*Volume; 
}