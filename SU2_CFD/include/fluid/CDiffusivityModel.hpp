/*!
 * \file CDiffusivityModel.hpp
 * \brief Interface class for defining mass diffusivity models.
 * \author T. Economon, C. Morales
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

#pragma once

#include "../../../Common/include/basic_types/datatype_structure.hpp"

/*!
 * \class CDiffusivityModel
 * \brief Defines a mass diffusivity model for species equations.
 * \author T. Economon
 */
class CDiffusivityModel {
 public:
  virtual ~CDiffusivityModel() = default;

  /*!
   * \brief Get mass diffusivity
   */
  su2double GetDiffusivity() const { return diff_; }

  /*!
   * \brief Set mass diffusivity
   */
  virtual void SetDiffusivity(su2double rho, su2double mu_lam, su2double cp, su2double kt) = 0;

 protected:
  su2double diff_{0.0}; /*!< \brief mass diffusivity. */
};
