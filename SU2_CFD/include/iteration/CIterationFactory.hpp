/*!
 * \file CIterationFactory.hpp
 * \brief Headers of the iteration classes used by SU2_CFD.
 *        Each CIteration class represents an available physics package.
 * \author F. Palacios, T. Economon
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

#include "CIteration.hpp"

/*!
 * \brief Creates a new iteration instance based on the current main solver.
 * \ingroup Drivers
 */
class CIterationFactory {
 public:
  /*!
   * \brief Deleted constructor to avoid creating instances of this class
   */
  CIterationFactory() = delete;

  /*!
   * \brief Create a new iteration instance based on the current main solver
   * \param[in] kindSolver - The kind of main solver we are running
   * \return               - Pointer to the allocated iteration instance
   */
  static CIteration* CreateIteration(MAIN_SOLVER kindSolver, const CConfig* config);
};
