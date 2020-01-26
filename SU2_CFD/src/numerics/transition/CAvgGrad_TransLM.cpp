/*!
 * \file CAvgGrad_TransLM.cpp
 * \brief Implementation of numerics class CAvgGrad_TransLM.
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

#include "../../../include/numerics/transition/CAvgGrad_TransLM.hpp"

CAvgGrad_TransLM::CAvgGrad_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iVar;

  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  sigma = 2./3.;

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTransVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTransVar_Edge = new su2double [nVar];
  Mean_GradTransVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTransVar[iVar] = new su2double [nDim];
}

CAvgGrad_TransLM::~CAvgGrad_TransLM(void) {

  unsigned short iVar;

  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTransVar_Kappa;
  delete [] Proj_Mean_GradTransVar_Edge;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTransVar[iVar];
  delete [] Mean_GradTransVar;
}

void CAvgGrad_TransLM::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
 /*--- This section is commented out on 04/11/2016
       after review of the static scan ---*/
 // su2double *Density_Grad_i      = new su2double[nDim];
 // su2double *Density_Grad_j      = new su2double[nDim];
 // su2double *Conservative_Grad_i = new su2double[nDim];
 // su2double *Conservative_Grad_j = new su2double[nDim];
 // su2double *Primitive_Grad_i    = new su2double[nDim];
 // su2double *Primitive_Grad_j    = new su2double[nDim];
 //
 // /*--- Intermediate values for combining viscosities ---*/
 // su2double Inter_Viscosity_i, Inter_Viscosity_j, REth_Viscosity_i, REth_Viscosity_j, Inter_Viscosity_Mean, REth_Viscosity_Mean;
 //
 // /*--- Model constants---*/
 // su2double sigmaf       = 1.0;
 // su2double sigma_thetat = 2.0;
 //
 // /*--- Get density ---*/
 // Density_i = U_i[0];
 // Density_j = U_j[0];
 //
 // /*--- Construct combinations of viscosity ---*/
 // Inter_Viscosity_i    = (Laminar_Viscosity_i+Eddy_Viscosity_i/sigmaf);
 // Inter_Viscosity_j    = (Laminar_Viscosity_j+Eddy_Viscosity_j/sigmaf);
 // Inter_Viscosity_Mean = 0.5*(Inter_Viscosity_i+Inter_Viscosity_j);
 // REth_Viscosity_i     = sigma_thetat*(Laminar_Viscosity_i+Eddy_Viscosity_i);
 // REth_Viscosity_j     = sigma_thetat*(Laminar_Viscosity_j+Eddy_Viscosity_j);
 // REth_Viscosity_Mean  = 0.5*(REth_Viscosity_i+REth_Viscosity_j);
 //
  ///*--- Compute vector going from iPoint to jPoint ---*/
  //dist_ij_2 = 0; proj_vector_ij = 0;
  //for (iDim = 0; iDim < nDim; iDim++) {
  //  Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
  //  dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  //  proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  //}
  //proj_vector_ij = proj_vector_ij/dist_ij_2; // to normalize vectors
 //
  ///*--- Mean gradient approximation ---*/
  //for (iVar = 0; iVar < nVar; iVar++) {
  //  Proj_Mean_GradTransVar_Kappa[iVar] = 0.0;
  //  // Proj_Mean_GradTransVar_Edge[iVar] = 0.0;
  //  for (iDim = 0; iDim < nDim; iDim++) {
 //
 //     /* -- Compute primitive grad using chain rule -- */
 //     Density_Grad_i[iDim]      = ConsVar_Grad_i[0][iDim];
 //     Density_Grad_j[iDim]      = ConsVar_Grad_j[0][iDim];
 //     Conservative_Grad_i[iDim] = TransVar_Grad_i[iVar][iDim];
 //     Conservative_Grad_j[iDim] = TransVar_Grad_j[iVar][iDim];
 //     Primitive_Grad_i[iDim]    = 1./Density_i*(Conservative_Grad_i[iDim]-TransVar_i[iVar]*Density_Grad_i[iDim]);
 //     Primitive_Grad_j[iDim]    = 1./Density_j*(Conservative_Grad_j[iDim]-TransVar_j[iVar]*Density_Grad_j[iDim]);
 //
 //     /*--- Compute the average primitive gradient and project it in the normal direction ---*/
 //     Mean_GradTransVar[iVar][iDim] = 0.5*(Primitive_Grad_i[iDim] + Primitive_Grad_j[iDim]);
  //    Proj_Mean_GradTransVar_Kappa[iVar] += Mean_GradTransVar[iVar][iDim]*Normal[iDim];
  //  }
  //}
 //
  //val_residual[0] = Inter_Viscosity_Mean*Proj_Mean_GradTransVar_Kappa[0];
  //val_residual[1] = REth_Viscosity_Mean*Proj_Mean_GradTransVar_Kappa[1];
 //
  ///*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  //if (implicit) {
  //  Jacobian_i[0][0] = (0.5*Proj_Mean_GradTransVar_Kappa[0]-Inter_Viscosity_Mean*proj_vector_ij);
  //  Jacobian_j[0][0] = (0.5*Proj_Mean_GradTransVar_Kappa[0]+Inter_Viscosity_Mean*proj_vector_ij);
  //  Jacobian_i[1][1] = (0.5*Proj_Mean_GradTransVar_Kappa[1]-REth_Viscosity_Mean*proj_vector_ij);
  //  Jacobian_j[1][1] = (0.5*Proj_Mean_GradTransVar_Kappa[1]+REth_Viscosity_Mean*proj_vector_ij);
  //}
 //
 // /*--- Free locally allocated memory. For efficiency, these arrays
 //  should really be allocated/deallocated in the constructor/destructor. ---*/
 // delete [] Density_Grad_i;
 // delete [] Density_Grad_j;
 // delete [] Conservative_Grad_i;
 // delete [] Conservative_Grad_j;
 // delete [] Primitive_Grad_i;
 // delete [] Primitive_Grad_j;
 //
}
