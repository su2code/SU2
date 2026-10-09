/*!
 * \file CParaviewBinaryFileWriter.hpp
 * \brief Headers fo paraview binary file writer class.
 * \author T. Albring
 * \version 8.5.0 "Harrier"
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

#include "CFileWriter.hpp"

class CParaviewBinaryFileWriter final: public CFileWriter{
  private:

  /*!
   * \brief Boolean storing whether we are on a big or little endian machine
   */
  bool bigEndian;

  static constexpr unsigned short NCOORDS = 3; /*!< \brief Points and vectors always have 3 components. */

  /*!
   * \brief Element types, in the order in which the cells are written.
   */
  static constexpr std::array<GEO_TYPE, 7> elemTypes = {LINE,       TRIANGLE, QUADRILATERAL, TETRAHEDRON,
                                                        HEXAHEDRON, PRISM,    PYRAMID};

  /*!
   * \brief Write the point coordinates.
   */
  void WritePoints();

  /*!
   * \brief Write the cells in the classic layout: the number of nodes followed by the node ids of each cell, Int32.
   * \param[in] GlobalCellStorage - Total size of that array.
   */
  void WriteCellsInt32(unsigned long GlobalCellStorage);

  /*!
   * \brief Write the cells in the layout of VTK >= 9.0: Int64 offsets and connectivity.
   */
  void WriteCellsInt64();

  /*!
   * \brief Write the type of each cell.
   */
  void WriteCellTypes();

  /*!
   * \brief Write the fields, as scalars or 3-component vectors.
   */
  void WritePointData();

  /*!
   * \brief Write 3 fields starting at firstVar as a vector (the third component is 0 in 2D).
   * \param[in] firstVar - Index of the first field in the data sorter.
   */
  void WriteVectorArray(unsigned short firstVar);

  /*!
   * \brief Write the point values of this rank (nComponents per point) at their place in the file.
   * \param[in,out] buffer - Values, byte-swapped in place to big endian.
   * \param[in] nComponents - Values per point.
   */
  void WritePointArray(vector<float>& buffer, unsigned short nComponents);

public:

  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] valDataSorter - The parallel sorted data to write
   */
  CParaviewBinaryFileWriter(CParallelDataSorter* valDataSorter);

  /*!
   * \brief Destructor
   */
  ~CParaviewBinaryFileWriter() override;

  /*!
   * \brief Write sorted data to file in paraview binary file format
   * \param[in] val_filename - The name of the file
   */
  void WriteData(string val_filename) override ;
};

