/*!
 * \file turb_diffusion.cpp
 * \brief Implementation of numerics classes to compute viscous
 *        fluxes in scalar problems.
 * \author F. Palacios, T. Economon
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/scalar/scalar_diffusion.hpp"

CAvgGrad_Scalar::CAvgGrad_Scalar(unsigned short val_nDim, unsigned short val_nVar, bool correct_grad,
                                 const CConfig* config)
    : CNumerics(val_nDim, val_nVar, config),
      correct_gradient(correct_grad),
      implicit(config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT),
      incompressible(config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE) {
  Proj_Mean_GradScalarVar_Normal = new su2double[nVar]();
  Proj_Mean_GradScalarVar = new su2double[nVar]();

  Flux = new su2double[nVar]();
  Jacobian_i = new su2double*[nVar];
  Jacobian_j = new su2double*[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double[nVar]();
    Jacobian_j[iVar] = new su2double[nVar]();
  }
}

CAvgGrad_Scalar::~CAvgGrad_Scalar(void) {
  delete[] Proj_Mean_GradScalarVar_Normal;
  delete[] Proj_Mean_GradScalarVar;

  delete[] Flux;
  if (Jacobian_i != nullptr) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete[] Jacobian_i[iVar];
      delete[] Jacobian_j[iVar];
    }
    delete[] Jacobian_i;
    delete[] Jacobian_j;
  }
}

CNumerics::ResidualType<> CAvgGrad_Scalar::ComputeResidual(const CConfig* config) {
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim);
  AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(ScalarVar_Grad_j, nVar, nDim);
  if (correct_gradient) {
    AD::SetPreaccIn(ScalarVar_i, nVar);
    AD::SetPreaccIn(ScalarVar_j, nVar);
  }
  ExtraADPreaccIn();

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim + 6);
    AD::SetPreaccIn(V_j, nDim + 6);

    Density_i = V_i[nDim + 2];
    Density_j = V_j[nDim + 2];
    Laminar_Viscosity_i = V_i[nDim + 4];
    Laminar_Viscosity_j = V_j[nDim + 4];
    Eddy_Viscosity_i = V_i[nDim + 5];
    Eddy_Viscosity_j = V_j[nDim + 5];
  } else {
    AD::SetPreaccIn(V_i, nDim + 7);
    AD::SetPreaccIn(V_j, nDim + 7);

    Density_i = V_i[nDim + 2];
    Density_j = V_j[nDim + 2];
    Laminar_Viscosity_i = V_i[nDim + 5];
    Laminar_Viscosity_j = V_j[nDim + 5];
    Eddy_Viscosity_i = V_i[nDim + 6];
    Eddy_Viscosity_j = V_j[nDim + 6];
  }

  proj_vector_ij = ComputeProjectedGradient(nDim, nVar, Normal, Coord_i, Coord_j, ScalarVar_Grad_i, ScalarVar_Grad_j,
                                            correct_gradient, ScalarVar_i, ScalarVar_j, Proj_Mean_GradScalarVar_Normal,
                                            Proj_Mean_GradScalarVar);
  FinishResidualCalc(config);

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}
