/*!
 * \file CFVMOutput.hpp
 * \brief  Headers of the Finite Volume Method output.
 * \author T. Kattmann
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

#include "COutput.hpp"

class CFVMOutput : public COutput{
 protected:
  /*!
   * \brief Constructor of the class
   */
  CFVMOutput(const CConfig *config, unsigned short nDim, bool femOutput);

  /*!
   * \brief Add Coordinates to output.
   */
  void AddCoordinates();

  /*!
   * \brief Load the coordinates.
   */
  template<class T>
  inline void LoadCoordinates(const T& Coord, const unsigned long iPoint) {
    SetVolumeOutputValue("COORD-X", iPoint, Coord[0]);
    SetVolumeOutputValue("COORD-Y", iPoint, Coord[1]);
    if (nDim == 3)
      SetVolumeOutputValue("COORD-Z", iPoint, Coord[2]);
  }

  /*!
   * \brief Add common FVM outputs.
   */
  void AddCommonFVMOutputs(const CConfig* config);

  /*!
   * \brief Load common FVM outputs.
   */
  void LoadCommonFVMOutputs(const CConfig* config, const CGeometry* geometry, unsigned long iPoint);
};
