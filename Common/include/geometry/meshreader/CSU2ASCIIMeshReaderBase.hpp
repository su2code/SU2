/*!
 * \file CSU2ASCIIMeshReaderBase.hpp
 * \brief Header file for the class CSU2ASCIIMeshReaderBase.
 *        The implementations are in the <i>CSU2ASCIIMeshReaderBase.cpp</i> file.
 * \author T. Economon
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "CSU2MeshReaderBase.hpp"

/*!
 * \class CSU2ASCIIMeshReaderBase
 * \brief Base class for the reading of a native SU2 ASCII grid.
 * \author T. Economon
 */
class CSU2ASCIIMeshReaderBase : public CSU2MeshReaderBase {
 protected:
  enum class FileSection { POINTS, ELEMENTS, MARKERS }; /*!< \brief Different sections of the file. */
  std::array<FileSection, 3> SectionOrder{};            /*!< \brief Order of the sections in the file. */

  ifstream mesh_file; /*!< \brief File object for the SU2 ASCII mesh file. */

  /*!
   * \brief Reads all SU2 ASCII mesh metadata and checks for errors.
   * \param[in] single_pass - Try to read the contents together with the metadata if the order allows (points before
   * elements). \param[in,out] config - Problem configuration where some metadata is updated (e.g. AoA). \returns True
   * if single_pass was successful.
   */
  bool ReadMetadata(bool single_pass, CConfig* config);

  /*!
   * \brief Reads the grid points from an SU2 zone into linear partitions across all ranks.
   */
  virtual void ReadPointCoordinates(bool single_pass = false);

  /*!
   * \brief Reads the interior volume elements from one section of an SU2 zone into linear partitions across all ranks.
   */
  virtual void ReadVolumeElementConnectivity(bool single_pass = false);

  /*!
   * \brief Reads the surface (boundary) elements from the SU2 zone.
   */
  virtual void ReadSurfaceElementConnectivity(bool single_pass = false);

  /*!
   * \brief Helper function to find the current zone in an SU2 ASCII mesh object.
   */
  void FastForwardToMyZone();

 public:
  /*!
   * \brief Constructor of the CSU2ASCIIMeshReaderBase class.
   */
  CSU2ASCIIMeshReaderBase(CConfig* val_config, unsigned short val_iZone, unsigned short val_nZone);

  /*!
   * \brief Destructor of the CSU2ASCIIMeshReaderBase class.
   */
  virtual ~CSU2ASCIIMeshReaderBase(void) override;
};
