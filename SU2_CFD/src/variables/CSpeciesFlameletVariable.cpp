/*!
 * \file CSpeciesFlameletVariable.cpp
 * \brief Definition of the variable fields for the flamelet class.
 * \author D. Mayer, T. Economon, N. Beishuizen
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/variables/CSpeciesFlameletVariable.hpp"

CSpeciesFlameletVariable::CSpeciesFlameletVariable(const su2double* species_inf, unsigned long npoint,
                                                   unsigned long ndim, unsigned long nvar, const CConfig* config)
    : CSpeciesVariable(species_inf, npoint, ndim, nvar, config) {
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution(iPoint, iVar) = species_inf[iVar];

  Solution_Old = Solution;

  /*--- Allocate and initialize solution for the dual time strategy ---*/
  bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));

  if (dual_time) {
    Solution_time_n = Solution;
    Solution_time_n1 = Solution;
  }

  /*--- Allocate residual structures ---*/

  Res_TruncError.resize(nPoint, nVar) = su2double(0.0);

  /* Allocate space for the source and scalars for visualization */
  const auto& flamelet_config_options = config->GetFlameletParsedOptions();
  source_scalar.resize(nPoint, flamelet_config_options.n_scalars) = su2double(0.0);
  lookup_scalar.resize(nPoint, flamelet_config_options.n_lookups) = su2double(0.0);
  source_pd.resize(nPoint, flamelet_config_options.n_control_vars) = su2double(0.0);
  table_misses.resize(nPoint) = 0;
  hull_miss_dcv1_.resize(nPoint) = su2double(0.0);
  hull_miss_dcv2_.resize(nPoint) = su2double(0.0);
  z_level_dist_.resize(nPoint)   = su2double(0.0);
  source_cons_jac.resize(nPoint, flamelet_config_options.n_user_scalars) = su2double(0.0);

  const bool beta_correction_active =
      flamelet_config_options.preferential_diffusion &&
      flamelet_config_options.pd_method != FLAMELET_PD_METHOD::SOURCE_TERM;
  const bool source_term_active =
      flamelet_config_options.preferential_diffusion &&
      flamelet_config_options.pd_method != FLAMELET_PD_METHOD::BETA_CORRECTION;

  /*--- Auxiliary variables: β scalars for the BETA_CORRECTION method, major species mass
   fractions (PREFERENTIAL_DIFFUSION_MAJOR_SPECIES) for the SOURCE_TERM method, whose gradients drive the Eq. (14)
   preferential diffusion fluxes (Schepers & van Oijen, C&F 280, 2025). Note that nAuxVar must
   be set: the SetAuxVar_Gradient_* routines iterate over GetnAuxVar() variables, so leaving it
   at its default of zero silently disables the aux-variable gradients. ---*/
  if (beta_correction_active) {
    nAuxVar = FLAMELET_PREF_DIFF_SCALARS::N_BETA_TERMS;
  } else if (source_term_active) {
    nAuxVar = flamelet_config_options.n_pd_major_species;
  }
  if (nAuxVar > 0) {
    AuxVar.resize(nPoint, nAuxVar) = su2double(0.0);
    Grad_AuxVar.resize(nPoint, nAuxVar, nDim, 0.0);
  }

  /*--- Eq. (14) flux coefficients for the SOURCE_TERM method: per control variable, the
   molecular coefficients D_{phi_k,i} for each major species plus one thermal (Soret)
   coefficient D^T_{phi_k}. ---*/
  if (source_term_active) {
    pd_terms_per_cv = FlameletPDTermsPerCV(flamelet_config_options.n_pd_major_species);
    pd_flux_coeff.resize(nPoint, flamelet_config_options.n_control_vars * pd_terms_per_cv) = su2double(0.0);
  }
}
