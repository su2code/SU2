/*!
 * \file graph_vertex_coloring.cpp
 * \brief Function used to carry out the vertex coloring of a given graph.
 * \author E. van der Weide
 * \version 5.0.0 "Raven"
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

#include "../include/graph_vertex_coloring.hpp"

void GraphVertexColoring(const vector<unsigned long>          &nVerticesPerRank,
                         const vector<vector<unsigned long> > &neighborVertices,
                         int                                  &nGlobalColors,
                         vector<int>                          &colorLocalVertices) {

  /* Determine the number of ranks and the current rank. */
  int nRank  = SINGLE_NODE;
  int myRank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRank);
#endif

  /*--------------------------------------------------------------------------*/
  /* Step 1: Gather all the data of the graph on the master rank, because the */
  /*         greedy algorithm used here is a sequential algorithm. Moreover   */
  /*         the algorithm requires the neighbors of neighbors, which is      */
  /*         difficult to implement in an efficient way in parallel.          */
  /****************************************************************************/

  cout << "GraphVertexColoring not implemented yet" << endl;
  exit(1);
}
