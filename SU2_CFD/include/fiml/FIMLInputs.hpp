/*!
 * \file FIMLInputs.hpp
 * \brief Catalogue of the inputs of the FIML neural network.
 * \version 8.6.0 "Harrier"
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

#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include "../../../Common/include/basic_types/datatype_structure.hpp"
#include "../../../Common/include/option_structure.hpp"

namespace FIML {

/*! \brief Number of inputs in the catalogue. */
constexpr unsigned short N_INPUTS = 9;

/*! \brief Scaling applied before the network: Box-Cox and z-score with fitted constants, or none. */
enum class TRANSFORM { BOX_COX, NONE };

/*!
 * \brief Local quantities of one point from which the inputs are computed.
 */
struct LocalState {
  su2double nu_tilde = 0.0;      /*!< \brief SA variable. */
  su2double density = 1.0;       /*!< \brief Density. */
  su2double mu_lam = 1.0;        /*!< \brief Laminar viscosity. */
  su2double mu_turb = 0.0;       /*!< \brief Eddy viscosity. */
  su2double strain = 0.0;        /*!< \brief Strain rate magnitude sqrt(2 S_ij S_ij). */
  su2double vorticity = 0.0;     /*!< \brief Vorticity magnitude. */
  su2double wall_dist = 0.0;     /*!< \brief Distance to the nearest wall. */
  passivedouble production = 0;  /*!< \brief SA production of the previous source evaluation. */
  passivedouble destruction = 0; /*!< \brief SA destruction of the previous source evaluation. */
  su2double tau_wall_ref = 1.0;  /*!< \brief Flat-plate wall shear stress of the free stream. */
};

/*!
 * \brief Name in the configuration and the network file, volume output field, and scaling of an input.
 */
struct InputDefinition {
  FIML_INPUT kind;
  const char* name;
  const char* field;
  const char* description;
  TRANSFORM transform;
};

/*!
 * \brief All inputs, in the order of FIML_INPUT.
 */
inline const std::array<InputDefinition, N_INPUTS>& Catalogue() {
  static const std::array<InputDefinition, N_INPUTS> catalogue = {{
      {FIML_INPUT::PROD_DEST, "PROD_DEST", "FIML_Prod_Dest", "SA production over destruction", TRANSFORM::BOX_COX},
      {FIML_INPUT::CHI, "CHI", "FIML_Chi", "Ratio of nu_tilde to the laminar viscosity", TRANSFORM::BOX_COX},
      {FIML_INPUT::DELTA, "DELTA", "FIML_Delta", "Eddy stress over the flat-plate wall stress", TRANSFORM::BOX_COX},
      {FIML_INPUT::STRAIN_VORT, "STRAIN_VORT", "FIML_Strain_Vort", "Strain rate over vorticity", TRANSFORM::BOX_COX},
      {FIML_INPUT::VISC_RATIO, "VISC_RATIO", "FIML_Visc_Ratio", "nu_t / (nu_t + 100 nu)", TRANSFORM::NONE},
      {FIML_INPUT::PROD_DEST_RATIO, "PROD_DEST_RATIO", "FIML_Prod_Dest_Ratio", "SA production / (production + destruction)",
       TRANSFORM::NONE},
      {FIML_INPUT::Q_CRIT, "Q_CRIT", "FIML_Q_Crit", "(Omega^2 - S^2) / (Omega^2 + S^2)", TRANSFORM::NONE},
      {FIML_INPUT::BL_RATIO, "BL_RATIO", "FIML_BL_Ratio", "1 / (1 + r_d)", TRANSFORM::NONE},
      {FIML_INPUT::FD, "FD", "FIML_Fd", "1 - tanh((8 r_d)^3)", TRANSFORM::NONE},
  }};
  return catalogue;
}

/*!
 * \brief Definition of one input.
 */
inline const InputDefinition& Definition(FIML_INPUT kind) { return Catalogue()[static_cast<unsigned short>(kind)]; }

/*!
 * \brief r_d of delayed DES, (nu_t + nu) / (kappa^2 d^2 |grad u|), limited to 10.
 */
inline su2double DistanceRatio(const LocalState& s) {
  const su2double grad_u = sqrt(0.5 * (s.strain * s.strain + s.vorticity * s.vorticity));
  const su2double denominator = 0.41 * 0.41 * s.wall_dist * s.wall_dist * grad_u;
  const su2double r_d = (s.mu_turb + s.mu_lam) / s.density / std::max(denominator, su2double(1e-30));
  return std::min(r_d, su2double(10.0));
}

/*!
 * \brief Shielding function 1 - tanh((C r_d)^3): 0 in attached boundary layers, 1 outside; C = 8 in delayed DES.
 */
inline su2double ShieldFunction(const LocalState& s, su2double coeff) {
  return 1.0 - tanh(pow(coeff * DistanceRatio(s), 3));
}

/*!
 * \brief Unscaled value of an input at a point.
 */
inline su2double Evaluate(FIML_INPUT kind, const LocalState& s) {
  switch (kind) {
    case FIML_INPUT::PROD_DEST:
      return s.production / (s.destruction + 1e-16);
    case FIML_INPUT::CHI:
      return s.nu_tilde * s.density / s.mu_lam;
    case FIML_INPUT::DELTA:
      return s.mu_turb * s.strain / (1.5 * s.tau_wall_ref);
    case FIML_INPUT::STRAIN_VORT:
      return s.strain / (s.vorticity + 1e-16);
    case FIML_INPUT::VISC_RATIO:
      return s.mu_turb / (s.mu_turb + 100.0 * s.mu_lam);
    case FIML_INPUT::PROD_DEST_RATIO:
      return s.production / (s.production + s.destruction + 1e-16);
    case FIML_INPUT::Q_CRIT: {
      const su2double omega2 = s.vorticity * s.vorticity, strain2 = s.strain * s.strain;
      return (omega2 - strain2) / (omega2 + strain2 + 1e-16);
    }
    case FIML_INPUT::BL_RATIO:
      return 1.0 / (1.0 + DistanceRatio(s));
    case FIML_INPUT::FD:
      return ShieldFunction(s, 8.0);
  }
  return 0.0;
}

}  // namespace FIML
