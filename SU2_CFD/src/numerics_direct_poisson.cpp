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

#include "../include/numerics_structure.hpp"
#include <limits>

CGalerkin_Flow::CGalerkin_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
}

CGalerkin_Flow::~CGalerkin_Flow(void) { }

void CGalerkin_Flow::ComputeResidual(su2double **val_stiffmatrix_elem, CConfig *config) {
  
  su2double a[4], b[4], c[4], d[4], Area, B_Matrix[4][4];
  unsigned short iVar, jVar;
  
  if (nDim == 2) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      a[iDim] = Coord_0[iDim]-Coord_2[iDim];
      b[iDim] = Coord_1[iDim]-Coord_2[iDim];
    }
    
    Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);  /* Norm of the normal component of area, area = 1/2*cross(a, b) */
    
    a[0] = 0.5 * (Coord_1[0]*Coord_2[1]-Coord_2[0]*Coord_1[1]) / Area;
    a[1] = 0.5 * (Coord_2[0]*Coord_0[1]-Coord_0[0]*Coord_2[1]) / Area;
    a[2] = 0.5 * (Coord_0[0]*Coord_1[1]-Coord_1[0]*Coord_0[1]) / Area;
    
    b[0] = 0.5 * (Coord_1[1]-Coord_2[1]) / Area;
    b[1] = 0.5 * (Coord_2[1]-Coord_0[1]) / Area;
    b[2] = 0.5 * (Coord_0[1]-Coord_1[1]) / Area;
    
    c[0] = 0.5 * (Coord_2[0]-Coord_1[0]) / Area;
    c[1] = 0.5 * (Coord_0[0]-Coord_2[0]) / Area;
    c[2] = 0.5 * (Coord_1[0]-Coord_0[0]) / Area;
    
    /* Compute the stiffness matrix, K & multiply it by the Area */
    
    val_stiffmatrix_elem[0][0] = Area * (b[0]*b[0]+c[0]*c[0]);
    val_stiffmatrix_elem[0][1] = Area * (b[0]*b[1]+c[0]*c[1]);
    val_stiffmatrix_elem[0][2] = Area * (b[0]*b[2]+c[0]*c[2]);
    val_stiffmatrix_elem[1][0] = Area * (b[0]*b[1]+c[0]*c[1]);
    val_stiffmatrix_elem[1][1] = Area * (b[1]*b[1]+c[1]*c[1]);
    val_stiffmatrix_elem[1][2] = Area * (b[1]*b[2]+c[1]*c[2]);
    val_stiffmatrix_elem[2][0] = Area * (b[0]*b[2]+c[0]*c[2]);
    val_stiffmatrix_elem[2][1] = Area * (b[1]*b[2]+c[1]*c[2]);
    val_stiffmatrix_elem[2][2] = Area * (b[2]*b[2]+c[2]*c[2]);
  }
  
  if (nDim == 3) {
    su2double Volume = 0.0;
    Volume -= Determinant_3x3(Coord_1[0], Coord_1[1], Coord_1[2], Coord_2[0], Coord_2[1], Coord_2[2], Coord_3[0], Coord_3[1], Coord_3[2]);
    Volume += Determinant_3x3(Coord_0[0], Coord_0[1], Coord_0[2], Coord_2[0], Coord_2[1], Coord_2[2], Coord_3[0], Coord_3[1], Coord_3[2]);
    Volume -= Determinant_3x3(Coord_0[0], Coord_0[1], Coord_0[2], Coord_1[0], Coord_1[1], Coord_1[2], Coord_3[0], Coord_3[1], Coord_3[2]);
    Volume += Determinant_3x3(Coord_0[0], Coord_0[1], Coord_0[2], Coord_1[0], Coord_1[1], Coord_1[2], Coord_2[0], Coord_2[1], Coord_2[2]);
    Volume = fabs(Volume / 6.0);
    
    a[0] = Determinant_3x3(Coord_1[0], Coord_1[1], Coord_1[2], Coord_2[0], Coord_2[1], Coord_2[2], Coord_3[0], Coord_3[1], Coord_3[2])/(6.0*Volume);
    b[0] = -Determinant_3x3(1.0, Coord_1[1], Coord_1[2],1.0, Coord_2[1], Coord_2[2],1.0, Coord_3[1], Coord_3[2])/(6.0*Volume);
    c[0] = -Determinant_3x3(Coord_1[0],1.0, Coord_1[2], Coord_2[0],1.0, Coord_2[2], Coord_3[0],1.0, Coord_3[2])/(6.0*Volume);
    d[0] = -Determinant_3x3(Coord_1[0], Coord_1[1],1.0, Coord_2[0], Coord_2[1],1.0, Coord_3[0], Coord_3[1],1.0)/(6.0*Volume);
    
    a[1] = -Determinant_3x3(Coord_2[0], Coord_2[1], Coord_2[2], Coord_3[0], Coord_3[1], Coord_3[2], Coord_0[0], Coord_0[1], Coord_0[2])/(6.0*Volume);
    b[1] = Determinant_3x3(1.0, Coord_2[1], Coord_2[2],1.0, Coord_3[1], Coord_3[2],1.0, Coord_0[1], Coord_0[2])/(6.0*Volume);
    c[1] = Determinant_3x3(Coord_2[0],1.0, Coord_2[2], Coord_3[0],1.0, Coord_3[2], Coord_0[0],1.0, Coord_0[2])/(6.0*Volume);
    d[1] = Determinant_3x3(Coord_2[0], Coord_2[1],1.0, Coord_3[0], Coord_3[1],1.0, Coord_0[0], Coord_0[1],1.0)/(6.0*Volume);
    
    a[2] = Determinant_3x3(Coord_3[0], Coord_3[1], Coord_3[2], Coord_0[0], Coord_0[1], Coord_0[2], Coord_1[0], Coord_1[1], Coord_1[2])/(6.0*Volume);
    b[2] = -Determinant_3x3(1.0, Coord_3[1], Coord_3[2],1.0, Coord_0[1], Coord_0[2],1.0, Coord_1[1], Coord_1[2])/(6.0*Volume);
    c[2] = -Determinant_3x3(Coord_3[0],1.0, Coord_3[2], Coord_0[0],1.0, Coord_0[2], Coord_1[0],1.0, Coord_1[2])/(6.0*Volume);
    d[2] = -Determinant_3x3(Coord_3[0], Coord_3[1],1.0, Coord_0[0], Coord_0[1],1.0, Coord_1[0], Coord_1[1],1.0)/(6.0*Volume);
    
    a[3] = -Determinant_3x3(Coord_0[0], Coord_0[1], Coord_0[2], Coord_1[0], Coord_1[1], Coord_1[2], Coord_2[0], Coord_2[1], Coord_2[2])/(6.0*Volume);
    b[3] = Determinant_3x3(1.0, Coord_0[1], Coord_0[2],1.0, Coord_1[1], Coord_1[2],1.0, Coord_2[1], Coord_2[2])/(6.0*Volume);
    c[3] = Determinant_3x3(Coord_0[0],1.0, Coord_0[2], Coord_1[0],1.0, Coord_1[2], Coord_2[0],1.0, Coord_2[2])/(6.0*Volume);
    d[3] = Determinant_3x3(Coord_0[0], Coord_0[1],1.0, Coord_1[0], Coord_1[1],1.0, Coord_2[0], Coord_2[1],1.0)/(6.0*Volume);
    
    /*--- Compute the B Matrix = grad N_j, dot grad N_i  ---*/
    B_Matrix[0][0] = b[0]*b[0] + c[0]*c[0] + d[0]*d[0];
    B_Matrix[0][1] = b[0]*b[1] + c[0]*c[1] + d[0]*d[1];
    B_Matrix[0][2] = b[0]*b[2] + c[0]*c[2] + d[0]*d[2];
    B_Matrix[0][3] = b[0]*b[3] + c[0]*c[3] + d[0]*d[3];
    
    B_Matrix[1][0] = b[1]*b[0] + c[1]*c[0] + d[1]*d[0];
    B_Matrix[1][1] = b[1]*b[1] + c[1]*c[1] + d[1]*d[1];
    B_Matrix[1][2] = b[1]*b[2] + c[1]*c[2] + d[1]*d[2];
    B_Matrix[1][3] = b[1]*b[3] + c[1]*c[3] + d[1]*d[3];
    
    B_Matrix[2][0] = b[2]*b[0] + c[2]*c[0] + d[2]*d[0];
    B_Matrix[2][1] = b[2]*b[1] + c[2]*c[1] + d[2]*d[1];
    B_Matrix[2][2] = b[2]*b[2] + c[2]*c[2] + d[2]*d[2];
    B_Matrix[2][3] = b[2]*b[3] + c[2]*c[3] + d[2]*d[3];
    
    B_Matrix[3][0] = b[3]*b[0] + c[3]*c[0] + d[3]*d[0];
    B_Matrix[3][1] = b[3]*b[1] + c[3]*c[1] + d[3]*d[1];
    B_Matrix[3][2] = b[3]*b[2] + c[3]*c[2] + d[3]*d[2];
    B_Matrix[3][3] = b[3]*b[3] + c[3]*c[3] + d[3]*d[3];
    
    
    /*--- Compute the BT.D.B Matrix (stiffness matrix) ---*/
    for (iVar = 0; iVar < 4; iVar++)
      for (jVar = 0; jVar < 4; jVar++)
        val_stiffmatrix_elem[iVar][jVar] = Volume * B_Matrix[iVar][jVar];
    
  }
}




CAvgGradCorrected_Poisson::CAvgGradCorrected_Poisson(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = false;
  implicit = (config->GetKind_TimeIntScheme_Poisson() == EULER_IMPLICIT);
  direct = false;

  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradPoissonVar_Edge = new su2double [nVar];
  Proj_Mean_GradPoissonVar_Kappa = new su2double [nVar];
  Proj_Mean_GradPoissonVar_Corrected = new su2double [nVar];
  Mean_GradPoissonVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradPoissonVar[iVar] = new su2double [nDim];

  ConsVar_Grad_i = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    ConsVar_Grad_i[iVar] = new su2double [nDim];
  
  ConsVar_Grad_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    ConsVar_Grad_j[iVar] = new su2double [nDim];

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
	
	for (iVar = 0; iVar < nVar; iVar++) 
       delete [] Mean_GradPoissonVar[iVar];
    delete [] Mean_GradPoissonVar;
    
    /*for (iVar = 0; iVar < nVar; iVar++) 
       if (ConsVar_Grad_i[iVar] != NULL) delete [] ConsVar_Grad_i[iVar];
    if (ConsVar_Grad_i != NULL) delete [] ConsVar_Grad_i;*/
    
    /*for (iVar = 0; iVar < nVar; iVar++)
       delete [] ConsVar_Grad_j[iVar];
    delete [] ConsVar_Grad_j;*/
    
    delete Mom_Coeff_i;
	delete Mom_Coeff_j;


}

void CAvgGradCorrected_Poisson::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {

  /*AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(Temp_i); AD::SetPreaccIn(Temp_j);
  AD::SetPreaccIn(ConsVar_Grad_i[0],nDim); AD::SetPreaccIn(ConsVar_Grad_j[0],nDim);
  AD::SetPreaccIn(Poisson_Coeff_i); AD::SetPreaccIn(Poisson_Coeff_j);*/

  Poisson_Coeff_Mean = 1.0;//0.5*(Poisson_Coeff_i + Poisson_Coeff_j);

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;

  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradPoissonVar_Edge[iVar] = 0.0;
    Proj_Mean_GradPoissonVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradPoissonVar_Corrected[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPoissonVar[iVar][iDim] = 0.5*(ConsVar_Grad_i[iVar][iDim] + ConsVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradPoissonVar_Kappa[iVar] += Mean_GradPoissonVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradPoissonVar_Edge[iVar] += Mean_GradPoissonVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradPoissonVar_Corrected[iVar] = (Poissonval_j-Poissonval_i)*proj_vector_ij;
    Proj_Mean_GradPoissonVar_Corrected[iVar] = Proj_Mean_GradPoissonVar_Corrected[iVar] + Proj_Mean_GradPoissonVar_Kappa[iVar] ;
    Proj_Mean_GradPoissonVar_Corrected[iVar] = Proj_Mean_GradPoissonVar_Corrected[iVar] - Proj_Mean_GradPoissonVar_Edge[iVar]*proj_vector_ij;
  }

  val_residual[0] = Poisson_Coeff_Mean*Proj_Mean_GradPoissonVar_Corrected[0];

  /*--- Jacobians for implicit scheme ---*/

  if (implicit) {
    Jacobian_i[0][0] = -Poisson_Coeff_Mean*proj_vector_ij;
    Jacobian_j[0][0] = Poisson_Coeff_Mean*proj_vector_ij;
  }

  /*AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();*/
}

CAvgGrad_Poisson::CAvgGrad_Poisson(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = false;
  implicit = (config->GetKind_TimeIntScheme_Poisson() == EULER_IMPLICIT);
  direct = false;

   
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradPoissonVar_Normal = new su2double [nVar];
  Proj_Mean_GradPoissonVar_Corrected = new su2double [nVar];
  Mom_Coeff_i = new su2double [nDim];
  Mom_Coeff_j = new su2double [nDim];
  
  Mean_GradPoissonVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
   Mean_GradPoissonVar[iVar] = new su2double [nDim];

  ConsVar_Grad_i = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    ConsVar_Grad_i[iVar] = new su2double [nDim];
  
  ConsVar_Grad_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    ConsVar_Grad_j[iVar] = new su2double [nDim];
  


  Poisson_Coeff_i = 1.0;
  Poisson_Coeff_j = 1.0;
}


CAvgGrad_Poisson::~CAvgGrad_Poisson(void) {
	
	unsigned int iDim;
	delete [] Edge_Vector;
	delete [] Proj_Mean_GradPoissonVar_Normal;
	delete [] Proj_Mean_GradPoissonVar_Corrected;
	/*delete [] Mom_Coeff_i;
	delete [] Mom_Coeff_j;
	
	for (iDim = 0; iVar < nVar; iVar++) 
      delete [] Mean_GradPoissonVar[iVar];
    
    delete [] Mean_GradPoissonVar;
    
    /*for (iDim = 0; iVar < nVar; iVar++) 
      delete [] ConsVar_Grad_i[iVar];
    
    delete [] ConsVar_Grad_i;
    
    for (iDim = 0; iVar < nVar; iVar++) 
      delete [] ConsVar_Grad_j[iVar];
   
    delete [] ConsVar_Grad_j;	*/
    
    

}

void CAvgGrad_Poisson::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {


  su2double Coeff_Mean;

  Poisson_Coeff_Mean = 1.0;//0.5*(Poisson_Coeff_i + Poisson_Coeff_j);

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;

  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradPoissonVar_Normal[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPoissonVar[iVar][iDim] = 0.5*(ConsVar_Grad_i[iVar][iDim] + ConsVar_Grad_j[iVar][iDim]);
      
      Coeff_Mean = 0.5*(Mom_Coeff_i[iDim] + Mom_Coeff_j[iDim]) ;
      
      Proj_Mean_GradPoissonVar_Normal[iVar] += Mean_GradPoissonVar[iVar][iDim]*Normal[iDim];
    }
  }

  val_residual[0] = Proj_Mean_GradPoissonVar_Normal[0];
  
  if (config->GetKind_Incomp_System() == PRESSURE_BASED) {
  Poisson_Coeff_Mean = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
     Poisson_Coeff_Mean += 0.5*Edge_Vector[iDim]*(Mom_Coeff_i[iDim] + Mom_Coeff_j[iDim])*Normal[iDim];
     
  Poisson_Coeff_Mean = Poisson_Coeff_Mean/dist_ij_2;
  if (dist_ij_2 == 0.0) cout<<"dist_ij is zero"<<endl;
  }
  else {
	  Poisson_Coeff_Mean = 1.0;
  }
  
  /*--- Jacobians for implicit scheme ---*/
  if (implicit) {
    Jacobian_i[0][0] = -Poisson_Coeff_Mean*proj_vector_ij;
    Jacobian_j[0][0] = Poisson_Coeff_Mean*proj_vector_ij;
  }

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
  su2double UnitNormal[3];

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
 
 su2double src_term;
 
 src_term = 0.0; //analytical solution, u = 1+2x^2+3y^2
 
 if (config->GetKind_Incomp_System()==PRESSURE_BASED) 
    val_residual[0] = Source_Term;
 else 
   val_residual[0] = Source_Term*Volume;
    
}


