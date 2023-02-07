/*!
 * \file CSU2ASCIIMeshReaderFVM.hpp
 * \brief Header file for the class CSU2ASCIIMeshReaderFVM.
 *        The implementations are in the <i>CSU2ASCIIMeshReaderFVM.cpp</i> file.
 * \author T. Economon
 * \version 7.5.1 "Blackbird"
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

#include <array>

#include "CSU2ASCIIMeshReaderBase.hpp"

/*!
 * \class CSU2ASCIIMeshReaderFVM
 * \brief Reads a native SU2 ASCII grid into linear partitions for the finite volume solver (FVM).
 * \author T. Economon
 */
class CSU2ASCIIMeshReaderFVM final: public CSU2ASCIIMeshReaderBase {
  
private:

  /*!
   * \brief Splits a single surface actuator disk boundary into two separate markers (repeated points).
   */
  void SplitActuatorDiskSurface();
  
public:

  /*!
   * \brief Constructor of the CSU2ASCIIMeshReaderFVM class.
   */
  CSU2ASCIIMeshReaderFVM(CConfig *val_config,
                         unsigned short val_iZone,
                         unsigned short val_nZone);

};
