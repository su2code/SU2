/*!
 * \file CAvgGrad_AdjTurb.cpp
 * \brief Implementation of numerics class CAvgGrad_AdjTurb.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/continuous_adjoint/CAvgGrad_AdjTurb.hpp"

CAvgGrad_AdjTurb::CAvgGrad_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  Edge_Vector = new su2double [nDim];
  Mean_GradTurbPsi = new su2double* [nVar];
  Proj_Mean_GradTurbPsi_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbPsi_Edge = new su2double [nVar];
  Proj_Mean_GradTurbPsi_Corrected = new su2double [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbPsi[iVar] = new su2double [nDim];
}

CAvgGrad_AdjTurb::~CAvgGrad_AdjTurb(void) {
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbPsi_Kappa;
  delete [] Proj_Mean_GradTurbPsi_Edge;
  delete [] Proj_Mean_GradTurbPsi_Corrected;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] Mean_GradTurbPsi[iVar];
  }
  delete [] Mean_GradTurbPsi;
}

void CAvgGrad_AdjTurb::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                                su2double **val_Jacobian_j, CConfig *config) {
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);
  
  su2double sigma = 2./3.;
  su2double nu_i, nu_j, nu_e;
  su2double dist_ij_2 = 0;
  su2double proj_vector_ij = 0;
  unsigned short iVar, iDim;
  
  /*--- Compute mean effective viscosity ---*/
  nu_i = Laminar_Viscosity_i/U_i[0];
  nu_j = Laminar_Viscosity_j/U_j[0];
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0])/sigma;
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbPsi_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbPsi[iVar][iDim] = 0.5*(TurbPsi_Grad_i[iVar][iDim] + TurbPsi_Grad_j[iVar][iDim]);
    }
    
    /*--- Projection of the corrected gradient ---*/
    Proj_Mean_GradTurbPsi_Corrected[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Proj_Mean_GradTurbPsi_Corrected[iVar] += Mean_GradTurbPsi[iVar][iDim]*Normal[iDim];
  }
  
  val_residual[0] = -nu_e*Proj_Mean_GradTurbPsi_Corrected[0];
  
  if (implicit) {
    val_Jacobian_i[0][0] =  nu_e*proj_vector_ij;
    val_Jacobian_j[0][0] = -nu_e*proj_vector_ij;
  }
  
}

void CAvgGrad_AdjTurb::ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                                                su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                                                su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) {
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);
  
  su2double sigma = 2./3.;
  su2double nu_i, nu_j, nu_e_i, nu_e_j;
  su2double dist_ij_2 = 0;
  su2double proj_vector_ij = 0;
  unsigned short iVar, iDim;
  
  /*--- Compute mean effective viscosity ---*/
  nu_i = Laminar_Viscosity_i/U_i[0];
  nu_j = Laminar_Viscosity_j/U_j[0];
  nu_e_i = (nu_i+TurbVar_i[0])/sigma;
  nu_e_j = (nu_j+TurbVar_j[0])/sigma;
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  proj_vector_ij = proj_vector_ij/dist_ij_2; // to normalize vectors
  
  /*--- Mean gradient approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbPsi_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbPsi[iVar][iDim] = 0.5*(TurbPsi_Grad_i[iVar][iDim] + TurbPsi_Grad_j[iVar][iDim]);
    }
    
    /*--- Projection of the corrected gradient ---*/
    Proj_Mean_GradTurbPsi_Corrected[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Proj_Mean_GradTurbPsi_Corrected[iVar] += Mean_GradTurbPsi[iVar][iDim]*Normal[iDim];
  }
  
  val_residual_i[0] = -nu_e_i*Proj_Mean_GradTurbPsi_Corrected[0];
  val_residual_j[0] =  nu_e_j*Proj_Mean_GradTurbPsi_Corrected[0];
  
  if (implicit) {
    val_Jacobian_ii[0][0] =  nu_e_i*proj_vector_ij;
    val_Jacobian_ij[0][0] = -nu_e_i*proj_vector_ij;
    val_Jacobian_ji[0][0] = -nu_e_j*proj_vector_ij;
    val_Jacobian_jj[0][0] =  nu_e_j*proj_vector_ij;
  }
  
}
