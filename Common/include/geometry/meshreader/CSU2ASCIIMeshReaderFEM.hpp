/*!
 * \file CSU2ASCIIMeshReaderFEM.hpp
 * \brief Header file for the class CSU2ASCIIMeshReaderFEM.
 *        The implementations are in the <i>CSU2ASCIIMeshReaderFEM.cpp</i> file.
 * \author T. Economon, E. van der Weide
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

#include <array>

#include "CSU2ASCIIMeshReaderBase.hpp"

/*!
 * \class CSU2ASCIIMeshReaderFEM
 * \brief Reads a native SU2 ASCII grid into linear partitions for the finite element solver (FEM).
 * \author T. Economon, E. van der Weide
 */
class CSU2ASCIIMeshReaderFEM : public CSU2ASCIIMeshReaderBase {
 private:
  /*!
   * \brief Reads the grid points from an SU2 zone into linear partitions across all ranks.
   */
  void ReadPointCoordinates();

  /*!
   * \brief Reads the interior volume elements from one section of an SU2 zone into linear partitions across all ranks.
   */
  void ReadVolumeElementConnectivity();

  /*!
   * \brief Reads the surface (boundary) elements from one section of an SU2 zone into linear partitions across all
   * ranks.
   */
  void ReadSurfaceElementConnectivity();

 public:
  /*!
   * \brief Constructor of the CSU2ASCIIMeshReaderFEM class.
   */
  CSU2ASCIIMeshReaderFEM(CConfig* val_config, unsigned short val_iZone, unsigned short val_nZone);

  /*!
   * \brief Destructor of the CSU2ASCIIMeshReaderFEM class.
   */
  ~CSU2ASCIIMeshReaderFEM(void) override;
};
