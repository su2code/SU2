/*!
 * \file graph_coloring_structure.cpp
 * \brief Functions used to carry out the coloring of a given graph.
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

#include "../include/graph_coloring_structure.hpp"

/* Function, which determines the colors for the vertices of the given graph. */
void CGraphColoringStructure::GraphVertexColoring(CConfig* config, const vector<unsigned long>& nVerticesPerRank,
                                                  const vector<vector<unsigned long> >& entriesVertices,
                                                  int& nGlobalColors, vector<int>& colorLocalVertices) {
  /* Determine the number of ranks and the current rank. */
  int nRank = 1;
  int myRank = 0;

#ifdef HAVE_MPI
  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &myRank);
  SU2_MPI::Comm_size(SU2_MPI::GetComm(), &nRank);
#endif

  /*--- Determine the algorithm to use for the graph coloring. ---*/
  switch (config->GetKind_Matrix_Coloring()) {
    case GREEDY_COLORING: {
      /* Greedy algorithm, which is implemented sequentially.
         Make a distinction between the master rank and the other ranks. */
      if (myRank == 0) {
        /*--------------------------------------------------------------------*/
        /*             Master node, which does all the work.                  */
        /* Step 1: Create the global vector for the graph by gathering all the*/
        /*         data from the other ranks. This is done, because the greedy*/
        /*         algorithm used here is a sequential algorithm. Moreover the*/
        /*         algorithm requires the neighbors of neighbors, which would */
        /*         require a rather cumbersome communication step anyway.     */
        /**************************************************************************/

        /* Define the global vector and copy my data in it. */
        vector<vector<unsigned long> > entriesVert(nVerticesPerRank[nRank], vector<unsigned long>(0));

        for (unsigned long i = nVerticesPerRank[0]; i < nVerticesPerRank[1]; ++i) entriesVert[i] = entriesVertices[i];

          /* Receive the data from the other ranks. Only needed in parallel mode. */
#ifdef HAVE_MPI
        for (int rank = 1; rank < nRank; ++rank) {
          /* Determine the size of the message to be received. */
          SU2_MPI::Status status;
          SU2_MPI::Probe(rank, rank, SU2_MPI::GetComm(), &status);

          int sizeMess;
          SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

          /* Allocate the memory for the receive buffer and receive the message. */
          vector<unsigned long> recvBuf(sizeMess);
          SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG, rank, rank, SU2_MPI::GetComm(), &status);

          /* Store the data just received in the global vector for the graph. */
          unsigned long ii = 0;
          for (unsigned long i = nVerticesPerRank[rank]; i < nVerticesPerRank[rank + 1]; ++i) {
            const unsigned long nEntries = recvBuf[ii++];

            entriesVert[i].resize(nEntries);
            for (unsigned long j = 0; j < nEntries; ++j, ++ii) entriesVert[i][j] = recvBuf[ii];
          }
        }
#endif

        /**********************************************************************/
        /* Step 2: The greedy algorithm to determine the vertex colors.       */
        /**********************************************************************/

        /* Allocate the memory to store the color of the vertices. */
        vector<int> colorVertices(nVerticesPerRank[nRank]);

        /* Allocate and reserve the memory for vectors that are used
           in the greedy algorithm. */
        vector<unsigned long> flagColorStored(nVerticesPerRank[nRank], nVerticesPerRank[nRank]);
        vector<int> colorNeighbors;
        colorNeighbors.reserve(1000);

        /* Loop over the vertices of the graph. */
        for (unsigned long i = 0; i < nVerticesPerRank[nRank]; ++i) {
          /* Make sure that colorNeighbors is empty. */
          colorNeighbors.resize(0);

          /* Loop over the entries of this vertex. */
          for (unsigned long j = 0; j < entriesVert[i].size(); ++j) {
            const unsigned long jj = entriesVert[i][j];

            /* Add the color of jj if jj is less than i and if its color
               is not stored yet. */
            if (jj < i) {
              const int cJJ = colorVertices[jj];
              if (flagColorStored[cJJ] != i) {
                flagColorStored[cJJ] = i;
                colorNeighbors.push_back(cJJ);
              }
            }

            /* Loop over the entries of vertex jj. */
            for (unsigned long k = 0; k < entriesVert[jj].size(); ++k) {
              const unsigned long kk = entriesVert[jj][k];

              /* Add the color of kk if kk is less than i and if its color
                 is not stored yet. */
              if (kk < i) {
                const int cKK = colorVertices[kk];
                if (flagColorStored[cKK] != i) {
                  flagColorStored[cKK] = i;
                  colorNeighbors.push_back(cKK);
                }
              }
            }
          }

          /* Sort colorNeighbors in increasing order. */
          sort(colorNeighbors.begin(), colorNeighbors.end());

          /* Determine the color of this vertex. */
          int cVert;
          for (cVert = 0; cVert < (int)colorNeighbors.size(); ++cVert) {
            if (cVert != colorNeighbors[cVert]) break;
          }

          colorVertices[i] = cVert;
        }

        /* Check if the coloring is done correctly. */
        flagColorStored.assign(nVerticesPerRank[nRank], nVerticesPerRank[nRank]);
        for (unsigned long i = 0; i < entriesVert.size(); ++i) {
          for (unsigned long j = 0; j < entriesVert[i].size(); ++j) {
            const int cJJ = colorVertices[entriesVert[i][j]];
            if (flagColorStored[cJJ] == i) {
              cout << "In function GraphVertexColoring: Color " << cJJ
                   << " appears more than once in the stencil of vertex " << i << "." << endl;
              exit(1);
            }
            flagColorStored[cJJ] = i;
          }
        }

        /**********************************************************************/
        /* Step 3: Store the coloring information in the local vectors again. */
        /**********************************************************************/

        /* Store the data of the root rank. */
        colorLocalVertices.resize(nVerticesPerRank[1]);
        memcpy(colorLocalVertices.data(), colorVertices.data(), nVerticesPerRank[1] * sizeof(int));

        /* Send the color of the vertices to the other ranks. Only in parallel
           mode. Use blocking sends, because deadlock cannot occur. */
#ifdef HAVE_MPI
        for (int rank = 1; rank < nRank; ++rank) {
          int* sendBuf = colorVertices.data() + nVerticesPerRank[rank];
          unsigned long sizeMess = nVerticesPerRank[rank + 1] - nVerticesPerRank[rank];
          SU2_MPI::Send(sendBuf, sizeMess, MPI_INT, rank, rank + 1, SU2_MPI::GetComm());
        }
#endif
      }
#ifdef HAVE_MPI
      else {
        /* Not the master node. Communicate the data of my graph to the master
           node. First copy the data in a buffer.  */
        vector<unsigned long> sendBuf;
        for (unsigned long i = 0; i < entriesVertices.size(); ++i) {
          sendBuf.push_back(entriesVertices[i].size());
          sendBuf.insert(sendBuf.end(), entriesVertices[i].begin(), entriesVertices[i].end());
        }

        /* Send the data to the master node. A blocking send can be used,
           because there is no danger of deadlock here. */
        SU2_MPI::Send(sendBuf.data(), sendBuf.size(), MPI_UNSIGNED_LONG, 0, myRank, SU2_MPI::GetComm());

        /* Receive the data for the colors of my locally owned DOFs. */
        unsigned long nLocalVert = entriesVertices.size();
        colorLocalVertices.resize(nLocalVert);

        SU2_MPI::Status status;
        SU2_MPI::Recv(colorLocalVertices.data(), nLocalVert, MPI_INT, 0, myRank + 1, SU2_MPI::GetComm(), &status);
      }
#endif
      break;
    }

      /*-----------------------------------------------------------------------*/

    case NATURAL_COLORING: {
      /* Natural coloring, i.e. every DOF gets its own color. This is very
        inefficient and should only be used for debugging. */
      unsigned long nLocalVert = entriesVertices.size();
      colorLocalVertices.resize(nLocalVert);

      for (unsigned long i = 0; i < nLocalVert; ++i) colorLocalVertices[i] = (int)(i + nVerticesPerRank[myRank]);

      break;
    }
  }

  /* Determine the total (global) number of colors. */
  int nLocalColors = 0;
  for (unsigned long i = 0; i < colorLocalVertices.size(); ++i) nLocalColors = max(nLocalColors, colorLocalVertices[i]);
  ++nLocalColors;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalColors, &nGlobalColors, 1, MPI_INT, MPI_MAX, SU2_MPI::GetComm());
#else
  nGlobalColors = nLocalColors;
#endif
}
