/*!
 * \file CPastixWrapper.cpp
 * \brief An interface to the INRIA solver PaStiX
 *        (http://pastix.gforge.inria.fr/files/README-txt.html)
 * \author P. Gomes
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

#ifdef HAVE_PASTIX

#include "../../include/parallelization/mpi_structure.hpp"
#include "../../include/parallelization/omp_structure.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/linear_algebra/CPastixWrapper.hpp"

#include <limits>
#include <numeric>

template <class ScalarType>
void CPastixWrapper<ScalarType>::Initialize(CGeometry* geometry, const CConfig* config) {
  if (isinitialized) return;  // only need to do this once

  const auto nVar = matrix.nVar, nPoint = matrix.nPoint, nPointDomain = matrix.nPointDomain;
  const auto *row_ptr = csr_row_ptr.data(), *col_ind = csr_col_ind.data();
  const auto nNonZero = row_ptr[nPointDomain];

  /*--- Allocate ---*/

  nCols = CheckedCast<pastix_int_t>(nPointDomain, "PaStiX local column count");
  colptr.resize(CheckedCast<size_t>(nPointDomain + 1, "PaStiX colptr size"));
  rowidx.clear();
  rowidx.reserve(CheckedCast<size_t>(nNonZero, "PaStiX rowidx size"));
  values.resize(CheckedCast<size_t>(CheckedMul(nNonZero, matrix.blkSz, "PaStiX value storage size"),
                                    "PaStiX value storage size"));
  loc2glb.resize(CheckedCast<size_t>(nPointDomain, "PaStiX loc2glb size"));
  perm.resize(CheckedCast<size_t>(nPointDomain, "PaStiX permutation size"));
  workvec.resize(
      CheckedCast<size_t>(CheckedMul(nPointDomain, nVar, "PaStiX work vector size"), "PaStiX work vector size"));

  /*--- Set default parameter values ---*/

  const auto incomplete = iparm[IPARM_INCOMPLETE];
  pastixInitParam(iparm, dparm);

  /*--- Customize important parameters ---*/

  switch (verb) {
    case 1:
      iparm[IPARM_VERBOSE] = PastixVerboseNo;
      break;
    case 2:
      iparm[IPARM_VERBOSE] = PastixVerboseYes;
      break;
    default:
      iparm[IPARM_VERBOSE] = PastixVerboseNot;
      break;
  }
  iparm[IPARM_ORDERING] = PastixOrderPtScotch;
  iparm[IPARM_INCOMPLETE] = incomplete;
  iparm[IPARM_LEVEL_OF_FILL] = static_cast<pastix_int_t>(config->GetPastixFillLvl());
  iparm[IPARM_THREAD_NBR] = omp_get_max_threads();

  pastixInit(&state, SU2_MPI::GetComm(), iparm, dparm);

  /*--- Prepare sparsity structure ---*/

  /*--- We need it in global coordinates, i.e. shifted according to the position
    of the current rank in the linear partitioning space, and "unpacked" halo part.
    The latter forces us to re-sort the column indices of rows with halo points, which
    in turn requires blocks to be swapped accordingly. Effectively the matrix is copied.
    Here we prepare the pointer and index part, and map the required swaps. ---*/

  /*--- 1 - Determine position in the linear partitioning ---*/

  su2_index_t offset = 0;
#ifdef HAVE_MPI
  vector<su2_index_t> domain_sizes(mpi_size);
  MPI_Allgather(&nPointDomain, 1, MPI_UINT64_T, domain_sizes.data(), 1, MPI_UINT64_T, SU2_MPI::GetComm());
  for (int i = 0; i < mpi_rank; ++i) offset = CheckedAdd(offset, domain_sizes[i], "PaStiX global offset");
#endif

  for (su2_index_t i = 0; i < nPointDomain; ++i) {
    loc2glb[static_cast<size_t>(i)] = CheckedCast<pastix_int_t>(CheckedAdd(offset, i + 1, "PaStiX global mapped index"),
                                                                "PaStiX global mapped index");
  }

  /*--- 2 - Communicate global indices of halo points to then renumber
   column indices from local to global when unpacking halos. ---*/

  vector<su2_index_t> map(CheckedCast<size_t>(nPoint - nPointDomain, "PaStiX halo map size"), 0);

#ifdef HAVE_MPI
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) && (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      unsigned short MarkerS = iMarker, MarkerR = iMarker + 1;

      int sender = config->GetMarker_All_SendRecv(MarkerS) - 1;
      int recver = abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;

      unsigned long nVertexS = geometry->nVertex[MarkerS];
      unsigned long nVertexR = geometry->nVertex[MarkerR];

      /*--- Allocate Send/Receive buffers ---*/
      vector<su2_index_t> Buffer_Recv(nVertexR), Buffer_Send(nVertexS);

      /*--- Prepare data to send ---*/
      for (unsigned long iVertex = 0; iVertex < nVertexS; iVertex++)
        Buffer_Send[iVertex] =
            CheckedAdd(geometry->vertex[MarkerS][iVertex]->GetNode(), offset, "PaStiX halo send index");

      /*--- Send and Receive data ---*/
      MPI_Sendrecv(Buffer_Send.data(), nVertexS, MPI_UINT64_T, sender, 0, Buffer_Recv.data(), nVertexR, MPI_UINT64_T,
                   recver, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);

      /*--- Store received data---*/
      for (unsigned long iVertex = 0; iVertex < nVertexR; iVertex++)
        map[geometry->vertex[MarkerR][iVertex]->GetNode() - nPointDomain] = Buffer_Recv[iVertex];
    }
  }
#endif

  /*--- 3 - Copy, map the sparsity, and put it in Fortran numbering ---*/

  for (su2_index_t iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    colptr[static_cast<size_t>(iPoint)] = CheckedCast<pastix_int_t>(
        CheckedAdd(row_ptr[iPoint], 1, "PaStiX row pointer entry"), "PaStiX row pointer entry");

    const su2_index_t begin = row_ptr[iPoint], end = row_ptr[iPoint + 1];

    /*--- If last point of row is halo ---*/
    const bool sort_required = (col_ind[end - 1] >= nPointDomain);

    if (sort_required) {
      const su2_index_t nnz_row = end - begin;

      sort_rows.push_back(iPoint);
      sort_order.emplace_back(CheckedCast<size_t>(nnz_row, "PaStiX row sort order size"));

      /*--- Sort mapped indices ("first") and keep track of source ("second")
            for when we later need to swap blocks for these rows. ---*/

      vector<pair<pastix_int_t, su2_index_t> > aux(CheckedCast<size_t>(nnz_row, "PaStiX row sort auxiliary size"));

      for (auto j = begin; j < end; ++j) {
        if (col_ind[j] < nPointDomain) {
          aux[j - begin].first = CheckedCast<pastix_int_t>(CheckedAdd(offset, col_ind[j] + 1, "PaStiX column index"),
                                                           "PaStiX column index");
        } else {
          aux[j - begin].first = CheckedCast<pastix_int_t>(
              CheckedAdd(map[col_ind[j] - nPointDomain], 1, "PaStiX halo column index"), "PaStiX halo column index");
        }
        aux[j - begin].second = j;
      }
      sort(aux.begin(), aux.end());

      for (su2_index_t j = 0; j < nnz_row; ++j) {
        rowidx.push_back(aux[j].first);
        sort_order.back()[j] = aux[j].second;
      }
    } else {
      /*--- These are all internal, no need to go through map. ---*/
      for (auto j = begin; j < end; ++j)
        rowidx.push_back(CheckedCast<pastix_int_t>(CheckedAdd(offset, col_ind[j] + 1, "PaStiX column index"),
                                                   "PaStiX column index"));
    }
  }
  colptr[static_cast<size_t>(nPointDomain)] =
      CheckedCast<pastix_int_t>(CheckedAdd(nNonZero, 1, "PaStiX final colptr entry"), "PaStiX final colptr entry");

  if (rowidx.size() != nNonZero) SU2_MPI::Error("Error during preparation of PaStiX data", CURRENT_FUNCTION);

  /*--- 4 - Perform ordering, symbolic factorization, and analysis steps ---*/

  spmInitDist(&spm, SU2_MPI::GetComm());
  spm.mtxtype = SpmGeneral;  // Despite being symmetric, we store the entire matrix.
  spm.flttype = std::is_same_v<su2mixedfloat, double> ? SpmDouble : SpmFloat;
  spm.fmttype = SpmCSC;
  spm.layout = SpmColMajor;
  spm.baseval = 1;

  spm.n = nCols;
  spm.nnz = CheckedCast<pastix_int_t>(nNonZero, "PaStiX nonzero count");
  spm.dof = CheckedCast<pastix_int_t>(nVar, "PaStiX block degree of freedom");

  spm.colptr = colptr.data();
  spm.rowptr = rowidx.data();
  spm.values = values.data();

  spm.replicated = static_cast<int>(mpi_size == 1);
  spm.loc2glob = mpi_size > 1 ? loc2glb.data() : nullptr;
  spmUpdateComputedFields(&spm);

  if (mpi_rank == MASTER_NODE && verb > 0) cout << endl;

  if (const auto rc = pastix_task_analyze(state, &spm); rc != PASTIX_SUCCESS) {
    SU2_MPI::Error("Error analyzing matrix: " + std::to_string(rc), CURRENT_FUNCTION);
  }

  if (mpi_rank == MASTER_NODE && verb > 0) cout << "+-------------------------------------------------+" << endl;

  isinitialized = true;
}

template <class ScalarType>
void CPastixWrapper<ScalarType>::AssembleValues() {
  const auto nDomain = matrix.nPointDomain;
  const auto blkSz = matrix.blkSz;
  const auto *d = matrix.d, *l = matrix.l, *u = matrix.u;
  for (su2_index_t iPoint = 0; iPoint < nDomain; ++iPoint) {
    auto* dst = values.data() + csr_row_ptr[iPoint] * blkSz;
    for (auto k = matrix.row_ptr_l[iPoint]; k < matrix.row_ptr_l[iPoint + 1]; ++k, dst += blkSz)
      for (auto b = 0ul; b < blkSz; ++b) dst[b] = SU2_TYPE::GetValue(l[k * blkSz + b]);
    for (auto b = 0ul; b < blkSz; ++b) dst[b] = SU2_TYPE::GetValue(d[iPoint * blkSz + b]);
    dst += blkSz;
    for (auto k = matrix.row_ptr_u[iPoint]; k < matrix.row_ptr_u[iPoint + 1]; ++k, dst += blkSz)
      for (auto b = 0ul; b < blkSz; ++b) dst[b] = SU2_TYPE::GetValue(u[k * blkSz + b]);
  }
}

template <class ScalarType>
void CPastixWrapper<ScalarType>::Factorize(CGeometry* geometry, const CConfig* config, unsigned short kind_fact) {
  /*--- Detect a possible change of settings between direct and adjoint that requires a reset ---*/
  if (isinitialized) {
    if ((kind_fact == PASTIX_ILU) != (iparm[IPARM_INCOMPLETE] == 1)) {
      Clean();
      iter = 0;
    }
  }
  verb = config->GetPastixVerbLvl();
  iparm[IPARM_INCOMPLETE] = (kind_fact == PASTIX_ILU);

  Initialize(geometry, config);

  /*--- Set some options that affect "compute" and could (one day) change during run ---*/

  switch (verb) {
    case 1:
      iparm[IPARM_VERBOSE] = PastixVerboseNo;
      break;
    case 2:
      iparm[IPARM_VERBOSE] = PastixVerboseYes;
      break;
    default:
      iparm[IPARM_VERBOSE] = PastixVerboseNot;
      break;
  }

  /*--- Is factorizing needed on this iteration? ---*/

  bool factorize = false;
  if (config->GetPastixFactFreq() != 0) factorize = (iter % config->GetPastixFactFreq() == 0);

  iter++;

  if (isfactorized && !factorize) return;  // No

  /*--- Yes: assemble LDU blocks into the flat CSR buffer ---*/
  AssembleValues();

  if (mpi_rank == MASTER_NODE && verb > 0) {
    cout << "\n+-------------------------------------------------+";
    cout << "\n+     PaStiX : Parallel Sparse matriX package     +" << endl;
  }

  const auto blkSz = matrix.blkSz;

  /*--- Permute blocks for rows with halo columns into global sorted order.
        AssembleValues wrote them in LDU order; copy to tmp then write back sorted. ---*/

  vector<su2mixedfloat> tmp;
  for (size_t i = 0; i < sort_rows.size(); ++i) {
    const auto iRow = sort_rows[i];
    /*--- colptr is 1-based Fortran numbering: row start = colptr[iRow] - 1. ---*/
    const auto begin = static_cast<su2_index_t>(colptr[iRow] - 1);
    const auto nnz_row = sort_order[i].size();

    tmp.assign(values.begin() + begin * blkSz, values.begin() + (begin + nnz_row) * blkSz);
    for (size_t j = 0; j < nnz_row; ++j) {
      const auto src_pos = sort_order[i][j] - begin;
      for (auto k = 0ul; k < blkSz; ++k) values[(begin + j) * blkSz + k] = tmp[src_pos * blkSz + k];
    }
  }

  /*--- Set factorization options ---*/

  switch (kind_fact) {
    case PASTIX_LDLT:
    case PASTIX_LDLT_P:
      iparm[IPARM_FACTORIZATION] = PastixFactLDLT;
      break;
    case PASTIX_LU:
    case PASTIX_LU_P:
    case PASTIX_ILU:
      iparm[IPARM_FACTORIZATION] = PastixFactLU;
      break;
    default:
      SU2_MPI::Error("Unknown type of PaStiX factorization.", CURRENT_FUNCTION);
      break;
  }

  /*--- Compute factorization ---*/

  if (const auto rc = pastix_task_numfact(state, &spm); rc != PASTIX_SUCCESS) {
    SU2_MPI::Error("Error factorizing matrix: " + std::to_string(rc), CURRENT_FUNCTION);
  }

  if (mpi_rank == MASTER_NODE && verb > 0) cout << "+-------------------------------------------------+\n" << endl;

  isfactorized = true;
}

#ifdef CODI_FORWARD_TYPE
template class CPastixWrapper<su2double>;
#else
template class CPastixWrapper<su2mixedfloat>;
#if defined(USE_MIXED_PRECISION) && !defined(USE_SINGLE_PRECISION)
template class CPastixWrapper<passivedouble>;
#endif
#endif
#endif
