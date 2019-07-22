/*!
 * \file graph_coloring_structure.hpp
 * \brief Include files and headers of the functions to carry out a
 *        coloring of a given graph. The functions are in the
 *        <i>graph_coloring_structure.cpp</i> file.
 * \author E. van der Weide
 * \version 6.2.0 "Falcon"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "./mpi_structure.hpp"
#include "option_structure.hpp"
#include "config_structure.hpp"
#include <cstring>

using namespace std;

/*!
 * \class CGraphColoringStructure
 * \brief Class, which provides graph coloring algorithms.
 * \author: E. van der Weide
 * \version 6.2.0 "Falcon"
 */
class CGraphColoringStructure {
public:
  /*!
   * \brief Constructor of the class. Nothing to be done.
   */
  CGraphColoringStructure(void);

  /*!
   * \brief Constructor of the class. Nothing to be done.
   */
  ~CGraphColoringStructure(void);

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
  void GraphVertexColoring(CConfig                              *config,
                           const vector<unsigned long>          &nVerticesPerRank,
                           const vector<vector<unsigned long> > &entriesVertices,
                           int                                  &nGlobalColors,
                           vector<int>                          &colorLocalVertices);
};
