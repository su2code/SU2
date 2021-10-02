/*!
 * \file turb_diffusion.cpp
 * \brief Implementation of numerics classes to compute viscous
 *        fluxes in turbulence problems.
 * \author T. Kattmann
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

#include "../../../include/numerics/species/species_diffusion.hpp"

CAvgGrad_Species::CAvgGrad_Species(unsigned short val_nDim, unsigned short val_nVar, bool correct_grad,
                                   const CConfig* config)
    : CAvgGrad_Scalar(val_nDim, val_nVar, correct_grad, config) {}

void CAvgGrad_Species::ExtraADPreaccIn() {}

void CAvgGrad_Species::FinishResidualCalc(const CConfig* config) {
  const su2double Sc_t = config->GetSchmidt_Number_Turbulent();
  const bool inc_rans = (config->GetKind_Solver() == INC_RANS) || (config->GetKind_Solver() == DISC_ADJ_INC_RANS);
  // const bool flame = (config->GetKind_Scalar_Model() == PROGRESS_VARIABLE);
  const bool flame = false;
  su2double Mass_Diffusivity_Lam;

  for (auto iVar = 0u; iVar < nVar; iVar++) {
    /*--- Compute the viscous residual. ---*/

    /* the general diffusion term for species transport is given by:
       (rho * D_{i,m} + mu_t/Sc_t ) * grad(Y_i))
       with D_{i,m} the mass diffusion coefficient of species i into the mixture m
     */

    if (flame)
      /* --- in case of combustion, Diffusion_Coeff from the lookup table is actually the complete diffusivity rho*D---
       */
      Mass_Diffusivity_Lam = 0.5 * (Diffusion_Coeff_i[iVar] + Diffusion_Coeff_j[iVar]);
    else
      /* --- in case of species transport, Diffusion_Coeff is the binary diffusion coefficient --- */
      Mass_Diffusivity_Lam = 0.5 * (Density_i * Diffusion_Coeff_i[iVar] + Density_j * Diffusion_Coeff_j[iVar]);

    su2double Mass_Diffusivity_Tur = 0.0;
    if (inc_rans) Mass_Diffusivity_Tur = 0.5 * (Eddy_Viscosity_i / Sc_t + Eddy_Viscosity_j / Sc_t);

    su2double Mass_Diffusivity = Mass_Diffusivity_Lam + Mass_Diffusivity_Tur;

    Flux[iVar] = Mass_Diffusivity * Proj_Mean_GradScalarVar[iVar];

    /*--- Use TSL approx. to compute derivatives of the gradients. ---*/

    if (implicit) {
      for (auto jVar = 0u; jVar < nVar; jVar++) {
        if (iVar == jVar) {
          su2double proj_on_rho = proj_vector_ij / Density_i;
          Jacobian_i[iVar][jVar] = -Mass_Diffusivity * proj_on_rho;
          proj_on_rho = proj_vector_ij / Density_j;
          Jacobian_j[iVar][jVar] = Mass_Diffusivity * proj_on_rho;
        } else {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
  }
}
