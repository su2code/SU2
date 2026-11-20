/*!
 * \file CParaviewPVTUFileWriter.hpp
 * \brief Headers of the Paraview parallel XML file writer (.pvtu and one .vtu piece per rank).
 * \author A. Rausa
 * \version 8.6.0 "Harrier"
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

#include <cstdint>

#include "CFileWriter.hpp"

/*!
 * \class CParaviewPVTUFileWriter
 * \brief Writes the data as a Paraview parallel unstructured grid: a .pvtu file that lists one .vtu
 *        piece per rank. Each rank writes its own piece, so no rank gathers data of other ranks and
 *        ParaView can read, filter and render the pieces in parallel (unlike a single .vtu file,
 *        which one process has to load as a whole).
 */
class CParaviewPVTUFileWriter final : public CFileWriter {
 private:
  /*!
   * \brief Boolean storing whether we are on a big or little endian machine.
   */
  bool bigEndian;

  /*!
   * \brief True to write the coordinates and fields in double precision instead of single.
   */
  bool doublePrecision;

  /*!
   * \brief Collect the points used by the elements of this rank and get their data from the ranks that own
   *        them, so that the piece of this rank is self-contained.
   * \param[out] usedPoints - Global indices of the points of the piece, sorted.
   * \param[out] pieceData - Fields of the points of the piece, field by field.
   */
  void CollectExternalPoints(vector<unsigned long>& usedPoints, vector<passivedouble>& pieceData) const;

  /*!
   * \brief Write the piece of this rank.
   * \param[in] val_filename - Name of the piece file.
   * \param[in] usedPoints - Global indices of the points of the piece, sorted.
   * \param[in] pieceData - Fields of the points of the piece, field by field.
   * \param[in] connInt64 - Write the connectivity and offsets as Int64 instead of Int32.
   * \returns True if the piece was written and closed without errors.
   */
  bool WritePiece(const string& val_filename, const vector<unsigned long>& usedPoints,
                  const vector<passivedouble>& pieceData, bool connInt64);

  /*!
   * \brief Write the .pvtu file that lists the pieces, only the master node writes it.
   * \param[in] val_filename - Name of the .pvtu file.
   * \param[in] pieceNames - Names of the pieces, relative to the .pvtu file.
   * \param[in] connInt64 - The connectivity and offsets of the pieces are Int64 instead of Int32.
   */
  bool WriteIndexFile(const string& val_filename, const vector<string>& pieceNames, bool connInt64);

 public:
  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] valDataSorter - The parallel sorted data to write.
   * \param[in] valDoublePrecision - Write the coordinates and fields in double precision.
   */
  CParaviewPVTUFileWriter(CParallelDataSorter* valDataSorter, bool valDoublePrecision = false);

  /*!
   * \brief Destructor
   */
  ~CParaviewPVTUFileWriter() override;

  /*!
   * \brief Write sorted data to the .pvtu file and its pieces.
   * \param[in] val_filename - The name of the file.
   */
  void WriteData(string val_filename) override;
};
