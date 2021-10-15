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

void CAvgGrad_Species::ExtraADPreaccIn() {
  AD::SetPreaccIn(Diffusion_Coeff_i, nVar);
  AD::SetPreaccIn(Diffusion_Coeff_j, nVar);
}

void CAvgGrad_Species::FinishResidualCalc(const CConfig* config) {
  for (auto iVar = 0u; iVar < nVar; iVar++) {
    const bool flame = false;  // const bool flame = (config->GetKind_Scalar_Model() == PROGRESS_VARIABLE);
    if (flame) {
      /*--- For combustion, Diffusion_Coeff from the lookup table is actually the complete diffusivity rho*D ---*/
      Density_i = 1.0;
      Density_j = 1.0;
    }
    /* --- in case of species transport, Diffusion_Coeff is the binary diffusion coefficient --- */
    const su2double Diffusivity_Lam = 0.5 * (Density_i * Diffusion_Coeff_i[iVar] + Density_j * Diffusion_Coeff_j[iVar]);

    const bool turbulence =
        (config->GetKind_Solver() == INC_RANS) ||
        (config->GetKind_Solver() == DISC_ADJ_INC_RANS);  // TODO TK:: this should be general for inc and comp
    const su2double Sc_t = config->GetSchmidt_Number_Turbulent();
    const su2double Diffusivity_Turb = turbulence ? 0.5 * (Eddy_Viscosity_i / Sc_t + Eddy_Viscosity_j / Sc_t) : 0.0;

    const su2double Diffusivity = Diffusivity_Lam + Diffusivity_Turb;

    Flux[iVar] = Diffusivity * Proj_Mean_GradScalarVar[iVar];

    /*--- Use TSL approx. to compute derivatives of the gradients. ---*/

    /*--- Off-diagonal entries are all zero. ---*/
    const su2double proj_on_rhoi = proj_vector_ij / Density_i;
    Jacobian_i[iVar][iVar] = -Diffusivity * proj_on_rhoi;

    const su2double proj_on_rhoj = proj_vector_ij / Density_j;
    Jacobian_j[iVar][iVar] = Diffusivity * proj_on_rhoj;

  }  // iVar
}
