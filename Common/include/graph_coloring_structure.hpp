/*!
 * \file graph_coloring_structure.hpp
 * \brief Include files and headers of the functions to carry out a
 *        coloring of a given graph. The functions are in the
 *        <i>graph_coloring_structure.cpp</i> file.
 * \author E. van der Weide
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

#include "./parallelization/mpi_structure.hpp"
#include "option_structure.hpp"
#include "CConfig.hpp"
#include <cstring>

using namespace std;

/*!
 * \class CGraphColoringStructure
 * \ingroup Graph
 * \brief Class, which provides distributed graph coloring algorithms.
 * \author: E. van der Weide
 * \version 8.0.0 "Harrier"
 */
class CGraphColoringStructure {
 public:
  /*!
   * \brief Function, which determines the colors for the vertices of the given graph.
   * \param[in]  config             - Definition of the particular problem.
   * \param[in]  nVerticesPerRank   - Number of vertices of the graph per MPI-rank,
                                      cumulative storage format.
   * \param[in]  entriesVertices    - The entries in the graph of the local vertices
                                      in global numbering.
   * \param[out] nGlobalColors      - Global number of colors in the graph.
   * \param[out] colorLocalVertices - The color of the local vertices of the graph.
   */
  void GraphVertexColoring(CConfig* config, const vector<unsigned long>& nVerticesPerRank,
                           const vector<vector<unsigned long> >& entriesVertices, int& nGlobalColors,
                           vector<int>& colorLocalVertices);
};
