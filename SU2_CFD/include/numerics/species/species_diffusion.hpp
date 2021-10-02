/*!
 * \file species_diffusion.hpp
 * \brief Declarations of numerics classes for discretization of
 *        viscous fluxes in species problems.
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

#pragma once

#include "../scalar/scalar_diffusion.hpp"

/*!
 * \class CAvgGrad_Species
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 */
class CAvgGrad_Species final : public CAvgGrad_Scalar {
 private:
  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void) override;

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_grad - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Species(unsigned short val_nDim, unsigned short val_nVar, bool correct_grad, const CConfig* config);
};
