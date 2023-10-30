/*!
 * \file CIntegrationFactory.hpp
 * \brief Headers of the CIntegrationFactory class
 * \author T. Albring
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

#include "../../include/solvers/CSolverFactory.hpp"

class CIntegration;
class CConfig;

/*!
 * \brief Factory for integration classes.
 * \ingroup Drivers
 */
class CIntegrationFactory{
 public:
  /*!
   * \brief Deleted constructor to avoid creating instances of this class
   */
  CIntegrationFactory() = delete;

  /*!
   * \brief Create the integration container based on the current main solver
   * \param[in] kindSolver       - The kind of main solver
   * \param[in] solver_container - The solver container
   * \return                  - Pointer to the allocated integration container
   */
  static CIntegration** CreateIntegrationContainer(MAIN_SOLVER kindSolver, const CSolver * const *solver_container);

  /*!
   * \brief Create a new integration instance based on the current sub solver
   * \param[in] integrationType  - The integration type
   * \return                     - Pointer to the allocated integration instance
   */
  static CIntegration* CreateIntegration(INTEGRATION_TYPE integrationType);

};
