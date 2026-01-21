/*!
 * \file CSU2BinaryMeshReaderBase.hpp
 * \brief Header file for the class CSU2BinaryMeshReaderBase.
 *        The implementations are in the <i>CSU2BinaryMeshReaderBase.cpp</i> file.
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

#include "CSU2MeshReaderBase.hpp"
#include "../../../include/toolboxes/SwapBytes.hpp"

/*!
 * \class CSU2BinaryMeshReaderBase
 * \brief Base class for the reading of a native SU2 binary grid.
 * \author T. Economon, E. van der Weide
 */
class CSU2BinaryMeshReaderBase : public CSU2MeshReaderBase {
 protected:
  constexpr static int SU2_STRING_SIZE = 33; /*!< \brief Size of the strings in the SU2 binary mesh file. */

  FILE* mesh_file;    /*!< \brief File object for the SU2 binary mesh file. */
  bool swap_bytes;    /*!< \brief Whether or not byte swapping must be used. */
  int size_conn_type; /*!< \brief Size, in bytes of the connectivity type. */

  /*!
   * \brief Reads the connectivity type used in the binary file and check if
   *        byte swapping must be applied.
   */
  void ReadConnectivityType();

  /*!
   * \brief Reads all SU2 binary mesh metadata and checks for errors.
   * \param[in,out] config - Problem configuration where some metadata is updated (e.g. AoA).
   */
  void ReadMetadata(CConfig* config);

  /*!
   * \brief Reads the grid points from an SU2 zone into linear partitions across all ranks.
   */
  virtual void ReadPointCoordinates();

  /*!
   * \brief Reads the interior volume elements from one section of an SU2 zone into linear partitions across all ranks.
   */
  virtual void ReadVolumeElementConnectivity();

  /*!
   * \brief Reads the surface (boundary) elements from the SU2 zone.
   */
  virtual void ReadSurfaceElementConnectivity();

  /*!
   * \brief Helper function to find the current zone in an SU2 binary mesh object.
   */
  void FastForwardToMyZone();

  /*!
   * \brief Function to read one entity of the connectivity type from the binary file.
   * \return uint64_t version of the the data.
   */
  uint64_t ReadBinaryNEntities();

  /*!
   * \brief Template function to read data from the binary file.
   */
  template <typename T>
  void ReadBinaryData(T* data, const size_t nItems) {
    /*--- Read the actual data. ---*/
    auto ret = fread(data, sizeof(T), nItems, mesh_file);
    if (ret != nItems) SU2_MPI::Error(string("Error while reading the file ") + meshFilename, CURRENT_FUNCTION);

    /*--- Apply byte swapping, if needed. ---*/
    if (swap_bytes) SwapBytes((char*)data, sizeof(T), nItems);
  }

 private:
  /*!
   * \brief Read the meta data for a zone.
   */
  void ReadMetadataZone();

 public:
  /*!
   * \brief Constructor of the CSU2BinaryMeshReaderBase class.
   */
  CSU2BinaryMeshReaderBase(CConfig* val_config, unsigned short val_iZone, unsigned short val_nZone);

  /*!
   * \brief Destructor of the CSU2BinaryMeshReaderBase class.
   */
  ~CSU2BinaryMeshReaderBase(void) override;
};
