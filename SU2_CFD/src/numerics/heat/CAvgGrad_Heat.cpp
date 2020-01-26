/*!
 * \file CAvgGrad_Heat.cpp
 * \brief Implementation of numerics class CAvgGrad_Heat.
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

#include "../../../include/numerics/heat/CAvgGrad_Heat.hpp"

CAvgGrad_Heat::CAvgGrad_Heat(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit        = (config->GetKind_TimeIntScheme_Heat() == EULER_IMPLICIT);

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradHeatVar_Normal = new su2double [nVar];
  Proj_Mean_GradHeatVar_Corrected = new su2double [nVar];
  Mean_GradHeatVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradHeatVar[iVar] = new su2double [nDim];

}

CAvgGrad_Heat::~CAvgGrad_Heat(void) {

  delete [] Edge_Vector;
  delete [] Proj_Mean_GradHeatVar_Normal;
  delete [] Proj_Mean_GradHeatVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradHeatVar[iVar];
  delete [] Mean_GradHeatVar;

}

void CAvgGrad_Heat::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(Temp_i); AD::SetPreaccIn(Temp_j);
  AD::SetPreaccIn(ConsVar_Grad_i[0],nDim); AD::SetPreaccIn(ConsVar_Grad_j[0],nDim);
  AD::SetPreaccIn(Thermal_Diffusivity_i); AD::SetPreaccIn(Thermal_Conductivity_j);

  Thermal_Diffusivity_Mean = 0.5*(Thermal_Diffusivity_i + Thermal_Diffusivity_j);

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
    Proj_Mean_GradHeatVar_Normal[iVar] = 0.0;
    Proj_Mean_GradHeatVar_Corrected[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradHeatVar[iVar][iDim] = 0.5*(ConsVar_Grad_i[iVar][iDim] + ConsVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradHeatVar_Normal[iVar] += Mean_GradHeatVar[iVar][iDim]*Normal[iDim];
    }
    Proj_Mean_GradHeatVar_Corrected[iVar] = Proj_Mean_GradHeatVar_Normal[iVar];
  }

  val_residual[0] = Thermal_Diffusivity_Mean*Proj_Mean_GradHeatVar_Corrected[0];

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {
    Jacobian_i[0][0] = -Thermal_Diffusivity_Mean*proj_vector_ij;
    Jacobian_j[0][0] = Thermal_Diffusivity_Mean*proj_vector_ij;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

