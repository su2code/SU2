/*!
 * \file CTecplotBinaryFileWriter.hpp
 * \brief Headers fo the tecplot binary writer class.
 * \author T. Albring
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

#include "CFileWriter.hpp"

#include <assert.h>

class CTecplotBinaryFileWriter final: public CFileWriter{

  unsigned long timeIter;  //!< Current value of the time iteration
  su2double timeStep;       //!< Current value of the time step

public:

  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] valDataSorter - The parallel sorted data to write
   * \param[in] valTimeIter - The current time iteration
   * \param[in] valTimeStep - The current physical time step value
   */
  CTecplotBinaryFileWriter(CParallelDataSorter* valDataSorter,
                           unsigned long valTimeIter, su2double valTimeStep);

  /*!
   * \brief Destructor
   */
  ~CTecplotBinaryFileWriter() override;

  /*!
   * \brief Write sorted data to file in tecplot binary file format
   * \param[in] val_filename - The name of the file
   */
  void WriteData(string val_filename) override ;

  /*!
   * \brief Calculate the partitioning of nodes to determine:
   * (a) For a given global node number, to which partition does it belong and what is its local node number; and
   * (b) How many nodes are held by each partition.
   */
  class NodePartitioner {
  public:
    /*!
     * \param[in] global_num_nodes - The total number of nodes being output
     * \param[in] num_ranks - The number of MPI ranks involved in the output
     */
    NodePartitioner(unsigned long global_num_nodes, int num_ranks)
      : m_num_ranks(num_ranks) {
      /* rank i has (1-based) global nodes m_node_range[i] + 1 through m_node_range[i + 1] */
      unsigned long nodes_per_rank = global_num_nodes / num_ranks;
      unsigned long num_extra_nodes = global_num_nodes - nodes_per_rank * num_ranks;
      m_node_range.resize(num_ranks + 1);
      m_node_range[0] = 0;
      for (int ii = 1; ii <= num_ranks; ii++) {
        m_node_range[ii] = m_node_range[ii - 1] + nodes_per_rank;
        if (num_extra_nodes > 0) {
          ++m_node_range[ii];
          --num_extra_nodes;
        }
      }
      assert(m_node_range[num_ranks] == global_num_nodes);
    }

    /*!
     * \brief Determine the MPI rank that owns a global node number and its corresponding local node number.
     * \param global_node_number[in] - The global node number; global node numbers are sequential across all MPI ranks.
     * \param owning_rank[out] - The MPI rank that owns (will output) the global node
     * \param node_number[out] - The rank-local node number for the given global node number
     */
    void GetOwningRankAndNodeNumber(unsigned long global_node_number, int &owning_rank, unsigned long &node_number)
    {
      owning_rank = static_cast<int>(global_node_number / m_node_range[1]);
      if (owning_rank >= m_num_ranks)
        owning_rank = m_num_ranks - 1;
      while(global_node_number > m_node_range[owning_rank + 1])
        ++owning_rank;
      while(global_node_number <= m_node_range[owning_rank])
        --owning_rank;
      node_number = global_node_number - m_node_range[owning_rank];
    }

    /*!
     * \brief Determine the number of nodes to be output by a particular rank
     * \param which_rank[in] - The MPI rank
     * \ret - The number of nodes that will be output by the give MPI rank.
     */
    int64_t GetRankNumNodes(int which_rank)
    {
      return static_cast<int64_t>(m_node_range[which_rank + 1] - m_node_range[which_rank]);
    }

  private:
    int m_num_ranks;
    vector<unsigned long> m_node_range;
  };

  int64_t GetHaloNodeNumber(unsigned long global_node_number, unsigned long last_local_node, vector<unsigned long> const &halo_node_list);

};

