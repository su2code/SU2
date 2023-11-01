/*!
 * \file CPastixWrapper.cpp
 * \brief An interface to the INRIA solver PaStiX
 *        (http://pastix.gforge.inria.fr/files/README-txt.html)
 * \author P. Gomes
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

#ifdef HAVE_PASTIX

#include "../../include/parallelization/mpi_structure.hpp"
#include "../../include/parallelization/omp_structure.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/linear_algebra/CPastixWrapper.hpp"

#include <numeric>

template <class ScalarType>
void CPastixWrapper<ScalarType>::Initialize(CGeometry* geometry, const CConfig* config) {
  using namespace PaStiX;

  if (isinitialized) return;  // only need to do this once

  unsigned long nVar = matrix.nVar, nPoint = matrix.nPoint, nPointDomain = matrix.nPointDomain;
  const unsigned long *row_ptr = matrix.rowptr, *col_ind = matrix.colidx;

  unsigned long iPoint, offset = 0, nNonZero = row_ptr[nPointDomain];

  /*--- Allocate ---*/

  nCols = pastix_int_t(nPointDomain);
  colptr.resize(nPointDomain + 1);
  rowidx.clear();
  rowidx.reserve(nNonZero);
  values.resize(nNonZero * nVar * nVar);
  loc2glb.resize(nPointDomain);
  perm.resize(nPointDomain);
  workvec.resize(nPointDomain * nVar);

  /*--- Set default parameter values ---*/

  pastix_int_t incomplete = iparm[IPARM_INCOMPLETE];

  iparm[IPARM_MODIFY_PARAMETER] = API_NO;
  Run();

  /*--- Customize important parameters ---*/

  switch (verb) {
    case 1:
      iparm[IPARM_VERBOSE] = API_VERBOSE_NO;
      break;
    case 2:
      iparm[IPARM_VERBOSE] = API_VERBOSE_YES;
      break;
    default:
      iparm[IPARM_VERBOSE] = API_VERBOSE_NOT;
      break;
  }
  iparm[IPARM_DOF_NBR] = pastix_int_t(nVar);
  iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
  iparm[IPARM_FREE_CSCPASTIX] = API_CSC_FREE;
  iparm[IPARM_CSCD_CORRECT] = API_NO;
  iparm[IPARM_RHSD_CHECK] = API_NO;
  iparm[IPARM_ORDERING] = API_ORDER_PTSCOTCH;
  iparm[IPARM_INCOMPLETE] = incomplete;
  iparm[IPARM_LEVEL_OF_FILL] = pastix_int_t(config->GetPastixFillLvl());
  iparm[IPARM_THREAD_NBR] = omp_get_max_threads();
#if defined(HAVE_MPI) && defined(HAVE_OMP)
  int comm_mode = MPI_THREAD_SINGLE;
  MPI_Query_thread(&comm_mode);
  if (comm_mode == MPI_THREAD_MULTIPLE)
    iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
  else
    iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_FUNNELED;
#endif

    /*--- Prepare sparsity structure ---*/

    /*--- We need it in global coordinates, i.e. shifted according to the position
     of the current rank in the linear partitioning space, and "unpacked" halo part.
     The latter forces us to re-sort the column indices of rows with halo points, which
     in turn requires blocks to be swapped accordingly. Moreover we need "pointer" and
     indices in Fortran-style numbering (start at 1), effectively the matrix is copied.
     Here we prepare the pointer and index part, and map the required swaps. ---*/

    /*--- 1 - Determine position in the linear partitioning ---*/

#ifdef HAVE_MPI
  vector<unsigned long> domain_sizes(mpi_size);
  MPI_Allgather(&nPointDomain, 1, MPI_UNSIGNED_LONG, domain_sizes.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
  for (int i = 0; i < mpi_rank; ++i) offset += domain_sizes[i];
#endif

  iota(loc2glb.begin(), loc2glb.end(), offset + 1);

  /*--- 2 - Communicate global indices of halo points to then renumber
   column indices from local to global when unpacking halos. ---*/

  vector<pastix_int_t> map(nPoint - nPointDomain, 0);

#ifdef HAVE_MPI
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) && (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      unsigned short MarkerS = iMarker, MarkerR = iMarker + 1;

      int sender = config->GetMarker_All_SendRecv(MarkerS) - 1;
      int recver = abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;

      unsigned long nVertexS = geometry->nVertex[MarkerS];
      unsigned long nVertexR = geometry->nVertex[MarkerR];

      /*--- Allocate Send/Receive buffers ---*/
      vector<unsigned long> Buffer_Recv(nVertexR), Buffer_Send(nVertexS);

      /*--- Prepare data to send ---*/
      for (unsigned long iVertex = 0; iVertex < nVertexS; iVertex++)
        Buffer_Send[iVertex] = geometry->vertex[MarkerS][iVertex]->GetNode() + offset;

      /*--- Send and Receive data ---*/
      MPI_Sendrecv(Buffer_Send.data(), nVertexS, MPI_UNSIGNED_LONG, sender, 0, Buffer_Recv.data(), nVertexR,
                   MPI_UNSIGNED_LONG, recver, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);

      /*--- Store received data---*/
      for (unsigned long iVertex = 0; iVertex < nVertexR; iVertex++)
        map[geometry->vertex[MarkerR][iVertex]->GetNode() - nPointDomain] = Buffer_Recv[iVertex];
    }
  }
#endif

  /*--- 3 - Copy, map the sparsity, and put it in Fortran numbering ---*/

  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    colptr[iPoint] = pastix_int_t(row_ptr[iPoint] + 1);

    unsigned long begin = row_ptr[iPoint], end = row_ptr[iPoint + 1], j;

    /*--- If last point of row is halo ---*/
    bool sort_required = (col_ind[end - 1] >= nPointDomain);

    if (sort_required) {
      unsigned long nnz_row = end - begin;

      sort_rows.push_back(iPoint);
      sort_order.push_back(vector<unsigned long>(nnz_row));

      /*--- Sort mapped indices ("first") and keep track of source ("second")
            for when we later need to swap blocks for these rows. ---*/

      vector<pair<pastix_int_t, unsigned long> > aux(nnz_row);

      for (j = begin; j < end; ++j) {
        if (col_ind[j] < nPointDomain)
          aux[j - begin].first = pastix_int_t(offset + col_ind[j] + 1);
        else
          aux[j - begin].first = pastix_int_t(map[col_ind[j] - nPointDomain] + 1);
        aux[j - begin].second = j;
      }
      sort(aux.begin(), aux.end());

      for (j = 0; j < nnz_row; ++j) {
        rowidx.push_back(aux[j].first);
        sort_order.back()[j] = aux[j].second;
      }
    } else {
      /*--- These are all internal, no need to go through map. ---*/
      for (j = begin; j < end; ++j) rowidx.push_back(pastix_int_t(offset + col_ind[j] + 1));
    }
  }
  colptr[nPointDomain] = pastix_int_t(nNonZero + 1);

  if (rowidx.size() != nNonZero) SU2_MPI::Error("Error during preparation of PaStiX data", CURRENT_FUNCTION);

  /*--- 4 - Perform ordering, symbolic factorization, and analysis steps ---*/

  if (mpi_rank == MASTER_NODE && verb > 0) cout << endl;

  iparm[IPARM_START_TASK] = API_TASK_ORDERING;
  iparm[IPARM_END_TASK] = API_TASK_ANALYSE;
  Run();

  if (mpi_rank == MASTER_NODE && verb > 0)
    cout << " +--------------------------------------------------------------------+" << endl;

  isinitialized = true;
}

template <class ScalarType>
void CPastixWrapper<ScalarType>::Factorize(CGeometry* geometry, const CConfig* config, unsigned short kind_fact) {
  using namespace PaStiX;

  /*--- Detect a possible change of settings between direct and adjoint that requires a reset ---*/
  if (isinitialized)
    if ((kind_fact == PASTIX_ILU) != (iparm[IPARM_INCOMPLETE] == API_YES)) {
      Clean();
      isinitialized = false;
      iter = 0;
    }

  verb = config->GetPastixVerbLvl();
  iparm[IPARM_INCOMPLETE] = (kind_fact == PASTIX_ILU);

  Initialize(geometry, config);

  /*--- Set some options that affect "compute" and could (one day) change during run ---*/

  switch (verb) {
    case 1:
      iparm[IPARM_VERBOSE] = API_VERBOSE_NO;
      break;
    case 2:
      iparm[IPARM_VERBOSE] = API_VERBOSE_YES;
      break;
    default:
      iparm[IPARM_VERBOSE] = API_VERBOSE_NOT;
      break;
  }

  if (kind_fact == PASTIX_LDLT || kind_fact == PASTIX_LDLT_P)
    iparm[IPARM_TRANSPOSE_SOLVE] = API_NO;  // symmetric so no need for slower transp. solve
  else
    iparm[IPARM_TRANSPOSE_SOLVE] = API_YES;  // negated due to CSR to CSC copy

  /*--- Is factorizing needed on this iteration? ---*/

  bool factorize = false;
  if (config->GetPastixFactFreq() != 0) factorize = (iter % config->GetPastixFactFreq() == 0);

  iter++;

  if (isfactorized && !factorize) return;  // No

  /*--- Yes ---*/

  if (mpi_rank == MASTER_NODE && verb > 0) {
    cout << endl;
    cout << " +--------------------------------------------------------------------+" << endl;
    cout << " +              PaStiX : Parallel Sparse matriX package               +" << endl;
    cout << " +--------------------------------------------------------------------+" << endl;
  }

  unsigned long i, j, k, iRow, begin, target, source, szBlk = matrix.nVar * matrix.nVar, nNonZero = values.size();

  /*--- Copy matrix values and swap blocks as required ---*/

  for (i = 0; i < nNonZero; ++i) values[i] = SU2_TYPE::GetValue(matrix.values[i]);

  for (i = 0; i < sort_rows.size(); ++i) {
    iRow = sort_rows[i];
    begin = matrix.rowptr[iRow];

    for (j = 0; j < sort_order[i].size(); ++j) {
      target = (begin + j) * szBlk;
      source = sort_order[i][j] * szBlk;

      for (k = 0; k < szBlk; ++k) values[target + k] = SU2_TYPE::GetValue(matrix.values[source + k]);
    }
  }

  /*--- Set factorization options ---*/

  switch (kind_fact) {
    case PASTIX_LDLT:
    case PASTIX_LDLT_P:
      iparm[IPARM_SYM] = API_SYM_YES;
      iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
      break;
    case PASTIX_LU:
    case PASTIX_LU_P:
    case PASTIX_ILU:
      iparm[IPARM_SYM] = API_SYM_NO;
      iparm[IPARM_FACTORIZATION] = API_FACT_LU;
      break;
    default:
      SU2_MPI::Error("Unknown type of PaStiX factorization.", CURRENT_FUNCTION);
      break;
  }

  /*--- Compute factorization ---*/

  iparm[IPARM_START_TASK] = API_TASK_NUMFACT;
  iparm[IPARM_END_TASK] = API_TASK_NUMFACT;
  Run();

  if (mpi_rank == MASTER_NODE && verb > 0)
    cout << " +--------------------------------------------------------------------+" << endl << endl;

  isfactorized = true;
}

#ifdef CODI_FORWARD_TYPE
template class CPastixWrapper<su2double>;
#else
template class CPastixWrapper<su2mixedfloat>;
#ifdef USE_MIXED_PRECISION
template class CPastixWrapper<passivedouble>;
#endif
#endif
#endif
