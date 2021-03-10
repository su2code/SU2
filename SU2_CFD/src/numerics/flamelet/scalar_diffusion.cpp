/*!
 * \file scalar_diffusion.cpp
 * \brief This file contains the numerical methods for scalar transport eqns.
 * \author T. Economon, D. Mayer, N. Beishuizen
 * \version 7.1.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/flamelet/scalar_diffusion.hpp"

CAvgGradtransportedScalar::CAvgGradtransportedScalar(unsigned short val_nDim,
                               unsigned short val_nVar,
                               bool correct_grad,
                               CConfig *config)
: CNumerics(val_nDim, val_nVar, config), correct_gradient(correct_grad) {
  
  implicit       = (config->GetKind_TimeIntScheme_Scalar() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  Edge_Vector = new su2double[nDim];
  
  Proj_Mean_GradScalarVar_Normal = new su2double[nVar];
  Proj_Mean_GradScalarVar_Edge   = new su2double[nVar];
  Proj_Mean_GradScalarVar        = new su2double[nVar];
  
  Mean_GradScalarVar = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradScalarVar[iVar] = new su2double[nDim];
}

CAvgGradtransportedScalar::~CAvgGradtransportedScalar(void) {
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradScalarVar_Normal;
  delete [] Proj_Mean_GradScalarVar_Edge;
  delete [] Proj_Mean_GradScalarVar;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradScalarVar[iVar];
  delete [] Mean_GradScalarVar;
}

void CAvgGradtransportedScalar::ComputeResidual(su2double *val_residual,
                                     su2double **Jacobian_i,
                                     su2double **Jacobian_j,
                                     CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(scalar_grad_i, nVar, nDim);
  AD::SetPreaccIn(scalar_grad_j, nVar, nDim);
  AD::SetPreaccIn(Diffusion_Coeff_i, nVar);
  AD::SetPreaccIn(Diffusion_Coeff_j, nVar);
  if (correct_gradient) {
    AD::SetPreaccIn(scalar_i, nVar); AD::SetPreaccIn(scalar_j, nVar);
  }
  ExtraADPreaccIn();
  
  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+6); AD::SetPreaccIn(V_j, nDim+6);
    
    Density_i           = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+4];  Laminar_Viscosity_j = V_j[nDim+4];
    Eddy_Viscosity_i    = V_i[nDim+5];     Eddy_Viscosity_j = V_j[nDim+5];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);
    
    Density_i           = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i    = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradScalarVar_Normal[iVar] = 0.0;
    Proj_Mean_GradScalarVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradScalarVar[iVar][iDim] = 0.5*(scalar_grad_i[iVar][iDim] +
                                            scalar_grad_j[iVar][iDim]);
      Proj_Mean_GradScalarVar_Normal[iVar] += Mean_GradScalarVar[iVar][iDim] *
      Normal[iDim];
      if (correct_gradient)
        Proj_Mean_GradScalarVar_Edge[iVar] += Mean_GradScalarVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradScalarVar[iVar] = Proj_Mean_GradScalarVar_Normal[iVar];
    if (correct_gradient) {
      Proj_Mean_GradScalarVar[iVar] -= Proj_Mean_GradScalarVar_Edge[iVar]*proj_vector_ij -
      (scalar_j[iVar]-scalar_i[iVar])*proj_vector_ij;
    }
  }
  
  FinishResidualCalc(val_residual, Jacobian_i, Jacobian_j, config);
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}

CAvgGradtransportedScalar_General::CAvgGradtransportedScalar_General(unsigned short val_nDim,
                                               unsigned short val_nVar, bool correct_grad,
                                               CConfig *config)
: CAvgGradtransportedScalar(val_nDim, val_nVar, correct_grad, config) {
  
  //Mean_Diffusivity = new su2double[nVar];
}

CAvgGradtransportedScalar_General::~CAvgGradtransportedScalar_General(void) {
  
  //if (Mean_Diffusivity != NULL) delete [] Mean_Diffusivity;
  
}

void CAvgGradtransportedScalar_General::ExtraADPreaccIn() { }

void CAvgGradtransportedScalar_General::FinishResidualCalc(su2double *val_residual,
                                                su2double **Jacobian_i,
                                                su2double **Jacobian_j,
                                                CConfig *config) {
  
  unsigned short iVar, jVar;

  for (iVar = 0; iVar < nVar; iVar++) {
    
    /*--- Get the diffusion coefficient(s). ---*/
    
    //Mean_Diffusivity[iVar] = 0.5*(Diffusion_Coeff_i[iVar] + Diffusion_Coeff_j[iVar]);

    /*--- Compute the viscous residual. ---*/

    //val_residual[iVar] = Mean_Diffusivity[iVar]*Proj_Mean_GradScalarVar[iVar];
    val_residual[iVar] = 1.0e-4*Proj_Mean_GradScalarVar[iVar];
    
    /*--- Use TSL approx. to compute derivatives of the gradients. ---*/

    if (implicit) {
      for (jVar = 0; jVar < nVar; jVar++) {
        if (iVar == jVar) {
          //Jacobian_i[iVar][jVar] = -Mean_Diffusivity[iVar]*proj_vector_ij;
          //Jacobian_j[iVar][jVar] =  Mean_Diffusivity[iVar]*proj_vector_ij;
          Jacobian_i[iVar][jVar] = -1.0e-4*proj_vector_ij;
          Jacobian_j[iVar][jVar] =  1.0e-4*proj_vector_ij;
        } else {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
  }
}