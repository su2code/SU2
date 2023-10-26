/*!
 * \file COutputFactory.hpp
 * \brief Headers of the output class.
 * \author T.Albring
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

#include "../../../Common/include/option_structure.hpp"

class COutput;
class CConfig;

class COutputFactory{

public:

  /*!
   * \brief Deleted constructor to avoid creating instances of this class
   */
  COutputFactory() = delete;

  /*!
   * \brief Create the Output class based on the current main solver
   * \param[in] kindSolver       - The kind of main solver
   * \param[in] config           - Pointer to the config
   * \param[in] nDim                - Number of physical dimensions
   * \return                     - Pointer to the allocated Output
   */
  static COutput* CreateOutput(MAIN_SOLVER kindSolver, CConfig *config, int nDim);

  /*!
   * \brief Create a multizone output
   * \param driverConfig        - Pointer to the driver config
   * \param config_container    - Pointer to the config container
   * \param nDim                - Number of physical dimensions
   * \return                    - Pointer to the allocated multizone output
   */
  static COutput* CreateMultizoneOutput(CConfig *driverConfig, CConfig** config_container, int nDim);
};
