/*!
 * \file CSTLFileWriter.hpp
 * \brief Headers fo the STL file writer class.
 * \author T. Kattmann, T. Albring
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "CFileWriter.hpp"
#include <set>


class CSTLFileWriter final : public CFileWriter{

  std::set<unsigned long> halo_nodes; /*!< \brief Solution of the problem. */
  vector<unsigned long> sorted_halo_nodes;

  vector<passivedouble> data_to_send;
  vector<passivedouble> halo_var_data;
  vector<int> num_values_to_send;
  vector<int> values_to_send_displacements;
  vector<int> values_to_receive_displacements;
  vector<unsigned long> nodes_to_send;
  vector<int> num_values_to_receive;
  vector<int> nodes_to_send_displacements;
  vector<int> nodes_to_receive_displacements;
  vector<int> num_nodes_to_send;
  size_t num_halo_nodes;
  vector<int> num_nodes_to_receive;

public:

  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names, file extension and dimension.
   * \param[in] fields - A list of field names
   * \param[in] nDim - Physical dimension
   * \param[in] fileName - The name of the file
   * \param[in] data_sorter - The parallel sorted data to write
   */
  CSTLFileWriter(vector<string> fields,
                 unsigned short nDim,
                 string fileName,
                 CParallelDataSorter* data_sorter);

  /*!
   * \brief Destructor
   */
  ~CSTLFileWriter() override;

  /*!
   * \brief Write sorted data to file in STL file format
   */
  void Write_Data() override;

  /*!
   * \brief Get the halo-node value of a global renumbered Point for a specific variable.
   * \param[in] global_node_number - global renumbered Point ID
   * \param[in] iVar - Variable number
   * \return Value of the halo-node variable
   */
  double GetHaloNodeValue(unsigned long global_node_number, unsigned short iVar);

};

