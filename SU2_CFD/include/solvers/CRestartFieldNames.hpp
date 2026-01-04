/*!
 * \file CRestartFieldNames.hpp
 * \brief Constants for restart file field names.
 * \note These names must match those used in output classes (CFlowIncOutput, CFlowCompOutput, etc.)
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include <string>

/*!
 * \namespace RestartFieldNames
 * \brief Namespace containing constants for restart file field names.
 * \note Field names in restart files match output field names (without quotes in output, with quotes in restart).
 */
namespace RestartFieldNames {
  /*--- Solution field names (must match output class field names) ---*/
  constexpr const char* PRESSURE = "Pressure";
  constexpr const char* TEMPERATURE = "Temperature";
  constexpr const char* DENSITY = "Density";
  constexpr const char* ENERGY = "Energy";
  
  /*--- Velocity component field names ---*/
  constexpr const char* VELOCITY_X = "Velocity_x";
  constexpr const char* VELOCITY_Y = "Velocity_y";
  constexpr const char* VELOCITY_Z = "Velocity_z";
  
  /*--- Momentum component field names ---*/
  constexpr const char* MOMENTUM_X = "Momentum_x";
  constexpr const char* MOMENTUM_Y = "Momentum_y";
  constexpr const char* MOMENTUM_Z = "Momentum_z";
  
  /*--- Grid velocity component field names ---*/
  constexpr const char* GRID_VELOCITY_X = "Grid_Velocity_x";
  constexpr const char* GRID_VELOCITY_Y = "Grid_Velocity_y";
  constexpr const char* GRID_VELOCITY_Z = "Grid_Velocity_z";
  
  /*!
   * \brief Get velocity field name for a given dimension index.
   * \param[in] iDim - Dimension index (0=x, 1=y, 2=z)
   * \return Field name string
   */
  inline std::string GetVelocityName(unsigned short iDim) {
    if (iDim == 0) return VELOCITY_X;
    if (iDim == 1) return VELOCITY_Y;
    return VELOCITY_Z;
  }
  
  /*!
   * \brief Get momentum field name for a given dimension index.
   * \param[in] iDim - Dimension index (0=x, 1=y, 2=z)
   * \return Field name string
   */
  inline std::string GetMomentumName(unsigned short iDim) {
    if (iDim == 0) return MOMENTUM_X;
    if (iDim == 1) return MOMENTUM_Y;
    return MOMENTUM_Z;
  }
  
  /*!
   * \brief Get grid velocity field name for a given dimension index.
   * \param[in] iDim - Dimension index (0=x, 1=y, 2=z)
   * \return Field name string
   */
  inline std::string GetGridVelocityName(unsigned short iDim) {
    if (iDim == 0) return GRID_VELOCITY_X;
    if (iDim == 1) return GRID_VELOCITY_Y;
    return GRID_VELOCITY_Z;
  }
}

