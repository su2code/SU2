/*!
 * \file CSU2BinaryMeshReaderFVM.hpp
 * \brief Header file for the class CSU2BinaryMeshReaderFVM.
 *        The implementations are in the <i>CSU2BinaryMeshReaderFVM.cpp</i> file.
 * \author T. Econonmon, E. van der Weide
 * \version 8.2.0 "Harrier"
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

#include "CSU2BinaryMeshReaderBase.hpp"

/*!
 * \class CSU2BinaryMeshReaderFVM
 * \brief Reads a native SU2 binary grid into linear partitions for the finite volume solver (FVM).
 * \author T. Econonmon, E. van der Weide
 */
class CSU2BinaryMeshReaderFVM : public CSU2BinaryMeshReaderBase {
 private:
  /*!
   * \brief Splits a single surface actuator disk boundary into two separate markers (repeated points).
   */
  void SplitActuatorDiskSurface();

 public:
  /*!
   * \brief Constructor of the CSU2BinaryMeshReaderFVM class.
   */
  CSU2BinaryMeshReaderFVM(CConfig* val_config, unsigned short val_iZone, unsigned short val_nZone);

  /*!
   * \brief Destructor of the CSU2BinaryMeshReaderFVM class.
   */
  ~CSU2BinaryMeshReaderFVM(void) override;
};
