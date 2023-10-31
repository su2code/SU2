/*!
 * \file CFileReaderLUT.hpp
 * \brief reading lookup table for tabulated fluid properties
 * \author D. Mayer, T. Economon
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

#include <fstream>
#include <string>
#include <vector>

#include "../../Common/include/parallelization/mpi_structure.hpp"
#include "../../../Common/include/linear_algebra/blas_structure.hpp"
#include "../../../Common/include/toolboxes/CSquareMatrixCM.hpp"

/*!
 * \brief File reader for look up tables.
 * \ingroup LookUpInterp
 */
class CFileReaderLUT {
 protected:
  int rank;

  unsigned short table_dim = 2;
  std::string version_lut;
  std::string version_reader;
  unsigned long n_levels = 1;
  su2vector<unsigned long> n_points, n_triangles, n_hull_points;

  unsigned long n_variables;
  su2vector<su2double> table_levels;

  /*! \brief Holds the variable names stored in the table file.
   * Order is in sync with tableFlamelet.
   */
  su2vector<std::string> names_var;

  /*! \brief Holds all data stored in the table.
   * First index addresses the variable while second index addresses the point.
   */
  su2vector<su2activematrix> table_data;

  su2vector<su2matrix<unsigned long>> triangles;

  su2vector<su2vector<unsigned long>> hull;
  /*! \brief Searches for the position of flag in file_stream and
   *         sets the stream position of file_stream to that position.
   */
  void SkipToFlag(std::ifstream& file_stream, const std::string& current_line, const std::string& flag) const;

  /*! \brief Extracts the next non-empty characters from file_stream and stores them into line.
   */
  bool GetNextNonEmptyLine(std::ifstream& file_stream, std::string& line) const;

  /*! \brief Extracts characters from file_stream, removes trailing control characters,
   *         and stores them into line.
   */
  bool GetStrippedLine(std::ifstream& file_stream, std::string& line) const;

 public:
  /*! \brief Get table version as listed in input file.
   */
  inline const std::string& GetVersionLUT() const { return version_lut; }

  /*! \brief Get table reader version.
   */
  inline const std::string& GetVersionReader() const { return version_reader; }

  /*! \brief Get number of data points at specific table level.
   * \param[in] i_level - table level index.
   * \returns data point count at table level.
   */
  inline unsigned long GetNPoints(std::size_t i_level = 0) const { return n_points[i_level]; }

  /*! \brief Get number of triangles at specific table level.
   * \param[in] i_level - table level index.
   * \returns triangle count at table level.
   */
  inline unsigned long GetNTriangles(std::size_t i_level = 0) const { return n_triangles[i_level]; }

  /*! \brief Get number of hull points at specific table level.
   * \param[in] i_level - table level index.
   * \returns hull point count at table level.
   */
  inline unsigned long GetNHullPoints(std::size_t i_level = 0) const { return n_hull_points[i_level]; }

  /*! \brief Get number of variables for which data is stored in the table
   */
  inline unsigned long GetNVariables() const { return n_variables; }

  /*! \brief Get number of table levels.
   */
  inline unsigned long GetNLevels() const { return n_levels; }

  /*! \brief Get variable names for which data is stored in the table
   */
  inline const su2vector<std::string>& GetNamesVar() const { return names_var; }

  /*! \brief Get table data at a specific level.
   * \param[in] i_level - table level index.
   * \returns table data
   */
  inline const su2activematrix& GetTableData(std::size_t i_level = 0) const { return table_data[i_level]; }

  /*! \brief Get table connectivity at a specific level.
   * \param[in] i_level - table level index.
   * \returns data connectivity
   */
  inline const su2matrix<unsigned long>& GetTriangles(std::size_t i_level = 0) const { return triangles[i_level]; }

  /*! \brief Get hull node information at a specific table level.
   * \param[in] i_level - table level index.
   * \returns hull node indices.
   */
  inline const su2vector<unsigned long>& GetHull(std::size_t i_level = 0) const { return hull[i_level]; }

  /*! \brief Get table level value.
   * \param[in] i_level - table level index.
   * \returns value of the third controlling variable at table level.
   */
  inline su2double GetTableLevel(std::size_t i_level) const { return table_levels[i_level]; }

  /*! \brief Get table dimension
   */
  inline unsigned short GetTableDim() const { return table_dim; }

  /*! \brief Read LUT file and store information
   * \param[in] file_name - LUT input file name.
   */
  void ReadRawLUT(const std::string& file_name);
};
