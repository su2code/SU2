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
#include "../../../include/numerics/viscous/CAvgGrad_P1.hpp"

CAvgGrad_P1::CAvgGrad_P1(unsigned short val_nDim,
                         unsigned short val_nVar,
                         CConfig *config)
                         : CNumericsRadiation(val_nDim, val_nVar, config) {

  // Initialization
  iVar = 0; iDim = 0;

  GammaP1 = 1.0 / (3.0*(Absorption_Coeff + Scattering_Coeff));

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradP1Var = new su2double [nVar];
  Mean_GradP1Var = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradP1Var[iVar] = new su2double [nDim];

}

CAvgGrad_P1::~CAvgGrad_P1(void) {

  delete [] Edge_Vector;
  delete [] Proj_Mean_GradP1Var;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradP1Var[iVar];
  delete [] Mean_GradP1Var;

}

void CAvgGrad_P1::ComputeResidual(su2double *val_residual,
                                  su2double **Jacobian_i,
                                  su2double **Jacobian_j,
                                  CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(RadVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(RadVar_Grad_j, nVar, nDim);

  /*--- Mean gradient approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradP1Var[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      /*--- Average gradients at faces ---*/
      Mean_GradP1Var[iVar][iDim] = 0.5*(RadVar_Grad_i[iVar][iDim] +
                                        RadVar_Grad_j[iVar][iDim]);
      /*--- Project over edge (including area information) ---*/
      Proj_Mean_GradP1Var[iVar] += Mean_GradP1Var[iVar][iDim] *
                                          Normal[iDim];
    }
  }

  /*--- Compute mean effective viscosity ---*/

  val_residual[0] = GammaP1*Proj_Mean_GradP1Var[0];

  if (implicit) {

    /*--- Compute vector going from iPoint to jPoint ---*/
    dist_ij = 0.0; proj_vector_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      dist_ij        += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      proj_vector_ij += (Coord_j[iDim]-Coord_i[iDim])*Normal[iDim];
    }
    if (dist_ij == 0.0){
      Jacobian_i[0][0] = 0.0;
      Jacobian_j[0][0] = 0.0;
    }
    else{
      proj_vector_ij = proj_vector_ij/dist_ij;
      Jacobian_i[0][0] = -GammaP1*proj_vector_ij;
      Jacobian_j[0][0] =  GammaP1*proj_vector_ij;
    }

  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}
