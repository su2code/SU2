/*!
 * \file radiation.cpp
 * \brief This file contains the implementation of the numerical
 *        methods for radiation.
 * \author Ruben Sanchez
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/numerics/radiation.hpp"

CNumericsRadiation::CNumericsRadiation(unsigned short val_nDim,
                         unsigned short val_nVar,
                         const CConfig *config)
                         : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Radiation() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);

  Absorption_Coeff = config->GetAbsorption_Coeff();
  Scattering_Coeff = config->GetScattering_Coeff();

  Absorption_Coeff = max(Absorption_Coeff,0.001);

  Temperature_Ref = config->GetTemperature_Ref();

}

CSourceP1::CSourceP1(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config) :
           CNumericsRadiation(val_nDim, val_nVar, config) { }

void CSourceP1::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {

  /*--- Retrieve the energy at the node i ---*/
  Energy_i = RadVar_i[0];

  /*--- Retrieve the temperature at the node i ---*/
  Temperature_i = V_i[nDim+1];

  /*--- Compute the blackbody intensity for gray media ---*/
  BlackBody_Intensity = 4.0*STEFAN_BOLTZMANN*pow(Temperature_i,4.0);

  /*--- Source term from black-body and energy contributions ---*/
  val_residual[0] = Absorption_Coeff * Volume * (BlackBody_Intensity - Energy_i);

  /*--- Contribution to the Jacobian ---*/
  if (implicit) {
    val_Jacobian_i[0][0] = - Absorption_Coeff * Volume;
  }

}

CAvgGradCorrected_P1::CAvgGradCorrected_P1(unsigned short val_nDim, unsigned short val_nVar,
                                           const CConfig *config) :
                      CNumericsRadiation(val_nDim, val_nVar, config) {

  GammaP1 = 1.0 / (3.0*(Absorption_Coeff + Scattering_Coeff));
}

void CAvgGradCorrected_P1::ComputeResidual(su2double *val_residual, su2double **Jacobian_i,
                                           su2double **Jacobian_j, CConfig *config) {
  constexpr int nVar = 1;

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(RadVar_i,nVar); AD::SetPreaccIn(RadVar_j,nVar);
  AD::SetPreaccIn(RadVar_Grad_i,nVar,nDim); AD::SetPreaccIn(RadVar_Grad_j,nVar,nDim);

  su2double NormalGrad[nVar], CorrectedGrad[nVar];

  auto proj_vector_ij = ComputeProjectedGradient(nDim, nVar, Normal, Coord_i, Coord_j, RadVar_Grad_i,
                                                 RadVar_Grad_j, true, RadVar_i, RadVar_j,
                                                 NormalGrad, CorrectedGrad);

  val_residual[0] = GammaP1*CorrectedGrad[0];

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  if (implicit) {
    Jacobian_i[0][0] = -GammaP1*proj_vector_ij;
    Jacobian_j[0][0] =  GammaP1*proj_vector_ij;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}
