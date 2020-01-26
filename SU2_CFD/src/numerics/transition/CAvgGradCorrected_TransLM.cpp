/*!
 * \file CAvgGradCorrected_TransLM.cpp
 * \brief Implementation of numerics class CAvgGradCorrected_TransLM.
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

#include "../../../include/numerics/transition/CAvgGradCorrected_TransLM.hpp"

CAvgGradCorrected_TransLM::CAvgGradCorrected_TransLM(unsigned short val_nDim, unsigned short val_nVar,
                                                     CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  unsigned short iVar;

  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  sigma = 2./3.;

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
}

CAvgGradCorrected_TransLM::~CAvgGradCorrected_TransLM(void) {

  unsigned short iVar;

  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
}

void CAvgGradCorrected_TransLM::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {

  //  switch (config->GetKind_Turb_Model()) {
  //  case SA :
  //    /*--- Compute mean effective viscosity ---*/
  //    nu_i = Laminar_Viscosity_i/U_i[0];
  //    nu_j = Laminar_Viscosity_j/U_j[0];
  //    nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  //
  //    /*--- Compute vector going from iPoint to jPoint ---*/
  //    dist_ij_2 = 0; proj_vector_ij = 0;
  //    for (iDim = 0; iDim < nDim; iDim++) {
  //      Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
  //      dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  //      proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  //    }
  //    proj_vector_ij = proj_vector_ij/dist_ij_2;
  //
  //    /*--- Mean gradient approximation. Projection of the mean gradient
  //       in the direction of the edge ---*/
  //    for (iVar = 0; iVar < nVar; iVar++) {
  //      Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
  //      Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
  //      for (iDim = 0; iDim < nDim; iDim++) {
  //        Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
  //        Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
  //        Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
  //      }
  //      Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
  //      Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
  //          (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  //    }
  //
  //    val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;
  //
  //    /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  //    if (implicit) {
  //      Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
  //      Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
  //    }
  //    break;
  //
  //  }
}
