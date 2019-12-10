/*!
 * \file numerics_direct_radiation.cpp
 * \brief Numerical methods for computing the viscous residual in the P1 equation.
 * \author Ruben Sanchez
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

#include "../../../include/numerics/CNumericsRadiation.hpp"
#include "../../../include/numerics/viscous/CAvgGradCorrected_P1.hpp"

CAvgGradCorrected_P1::CAvgGradCorrected_P1(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumericsRadiation(val_nDim, val_nVar, config) {

  // Initialization
  iVar = 0; iDim = 0;

  GammaP1 = 1.0 / (3.0*(Absorption_Coeff + Scattering_Coeff));

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradP1Var_Edge = new su2double [nVar];
  Proj_Mean_GradP1Var_Kappa = new su2double [nVar];
  Proj_Mean_GradP1Var_Corrected = new su2double [nVar];
  Mean_GradP1Var = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradP1Var[iVar] = new su2double [nDim];

}

CAvgGradCorrected_P1::~CAvgGradCorrected_P1(void) {

  delete [] Edge_Vector;
  delete [] Proj_Mean_GradP1Var_Edge;
  delete [] Proj_Mean_GradP1Var_Kappa;
  delete [] Proj_Mean_GradP1Var_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradP1Var[iVar];
  delete [] Mean_GradP1Var;

}

void CAvgGradCorrected_P1::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(RadVar_i,nVar); AD::SetPreaccIn(RadVar_j,nVar);
  AD::SetPreaccIn(RadVar_Grad_i,nVar,nDim); AD::SetPreaccIn(RadVar_Grad_j,nVar,nDim);

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij;

  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradP1Var_Edge[iVar] = 0.0;
    Proj_Mean_GradP1Var_Kappa[iVar] = 0.0;

    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradP1Var[iVar][iDim] = 0.5*(RadVar_Grad_i[iVar][iDim] + RadVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradP1Var_Kappa[iVar] += Mean_GradP1Var[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradP1Var_Edge[iVar] += Mean_GradP1Var[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradP1Var_Corrected[iVar] = Proj_Mean_GradP1Var_Kappa[iVar];
    Proj_Mean_GradP1Var_Corrected[iVar] -= Proj_Mean_GradP1Var_Edge[iVar]*proj_vector_ij -
                                           (RadVar_j[iVar]-RadVar_i[iVar])*proj_vector_ij;
  }

  val_residual[0] = GammaP1*Proj_Mean_GradP1Var_Corrected[0];

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  if (implicit) {
    Jacobian_i[0][0] = -GammaP1*proj_vector_ij;
    Jacobian_j[0][0] =  GammaP1*proj_vector_ij;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}
