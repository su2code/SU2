/*!
 * \file CPastixWrapper.hpp
 * \brief An interface to the INRIA solver PaStiX
 *        (http://pastix.gforge.inria.fr/files/README-txt.html)
 * \author P. Gomes
 * \version 8.4.0 "Harrier"
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

#ifdef HAVE_PASTIX

#ifdef CODI_FORWARD_TYPE
#error Cannot use PaStiX with forward mode AD
#endif

#include <pastix.h>
#include <spm.h>
#include <vector>

using namespace std;

class CConfig;
class CGeometry;

/*!
 * \class CPastixWrapper
 * \ingroup SpLinSys
 * \brief Wrapper class that converts between SU2 sparse format and PaStiX
 *        format and simplifies calls to the external solver.
 */
template <class ScalarType>
class CPastixWrapper {
 private:
  pastix_data_t* state{};        /*!< \brief Internal state of the solver. */
  spmatrix_t spm;                /*!< \brief Matrix format used by the solver. */
  pastix_int_t nCols;            /*!< \brief Local number of columns. */
  vector<pastix_int_t> colptr;   /*!< \brief Equiv. to our "row_ptr". */
  vector<pastix_int_t> rowidx;   /*!< \brief Equiv. to our "col_ind". */
  vector<su2mixedfloat> values;  /*!< \brief Equiv. to our "matrix". */
  vector<pastix_int_t> loc2glb;  /*!< \brief Global index of the columns held by this rank. */
  vector<pastix_int_t> perm;     /*!< \brief Ordering computed by PaStiX. */
  vector<su2mixedfloat> workvec; /*!< \brief RHS vector which then becomes the solution. */

  pastix_int_t iparm[IPARM_SIZE];  /*!< \brief Integer parameters for PaStiX. */
  passivedouble dparm[DPARM_SIZE]; /*!< \brief Floating point parameters for PaStiX. */

  struct {
    unsigned long nVar = 0;
    unsigned long nPoint = 0;
    unsigned long nPointDomain = 0;
    const unsigned long* rowptr = nullptr;
    const unsigned long* colidx = nullptr;
    const ScalarType* values = nullptr;

    unsigned long size_rhs() const { return nPointDomain * nVar; }
  } matrix; /*!< \brief Pointers and sizes of the input matrix. */

  bool issetup{};        /*!< \brief Signals that the matrix data has been provided. */
  bool isinitialized{};  /*!< \brief Signals that the sparsity pattern has been set. */
  bool isfactorized{};   /*!< \brief Signals that a factorization has been computed. */
  bool transpose{};      /*!< \brief Solve A^T x = b instead of A x = b. */
  unsigned long iter{};  /*!< \brief Number of times a factorization has been requested. */
  unsigned short verb{}; /*!< \brief Verbosity level. */
  const int mpi_size = SU2_MPI::GetSize();
  const int mpi_rank = SU2_MPI::GetRank();

  vector<unsigned long> sort_rows;          /*!< \brief List of rows with halo points. */
  vector<vector<unsigned long>> sort_order; /*!< \brief How each of those rows needs to be sorted. */

  /*!
   * \brief Run the "clean" task, releases all memory, leaves object in unusable state.
   */
  void Clean() {
    if (isinitialized) {
      iparm[IPARM_VERBOSE] = (verb > 0) ? PastixVerboseNo : PastixVerboseNot;
      pastixFinalize(&state);
      spm.colptr = nullptr;
      spm.rowptr = nullptr;
      spm.values = nullptr;
      if (mpi_size > 1) spm.loc2glob = nullptr;
      spmExit(&spm);
      isfactorized = false;
      isinitialized = false;
    }
  }

  /*!
   * \brief Initialize the matrix format that PaStiX requires.
   */
  void Initialize(CGeometry* geometry, const CConfig* config);

 public:
  CPastixWrapper() = default;

  /*--- Move or copy is not allowed. ---*/
  CPastixWrapper(CPastixWrapper&&) = delete;
  CPastixWrapper(const CPastixWrapper&) = delete;
  CPastixWrapper& operator=(CPastixWrapper&&) = delete;
  CPastixWrapper& operator=(const CPastixWrapper&) = delete;

  /*!
   * \brief Class destructor.
   */
  ~CPastixWrapper() { Clean(); }

  /*!
   * \brief Set matrix data, only once.
   * \param[in] nVar - DOF per point.
   * \param[in] nPoint - Total number of points including halos.
   * \param[in] nPointDomain - Number of internal points.
   * \param[in] rowptr - Array, where column index data starts for each matrix row.
   * \param[in] colidx - Non zeros column indices.
   * \param[in] values - Matrix coefficients.
   */
  void SetMatrix(unsigned long nVar, unsigned long nPoint, unsigned long nPointDomain, const unsigned long* rowptr,
                 const unsigned long* colidx, const ScalarType* values) {
    if (issetup) return;
    matrix.nVar = nVar;
    matrix.nPoint = nPoint;
    matrix.nPointDomain = nPointDomain;
    matrix.rowptr = rowptr;
    matrix.colidx = colidx;
    matrix.values = values;
    issetup = true;
  }

  /*!
   * \brief Factorize matrix.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] kind_fact - Type of factorization.
   */
  void Factorize(CGeometry* geometry, const CConfig* config, unsigned short kind_fact);

  /*!
   * \brief Request solves with the transposed matrix.
   * \param[in] transposed - Yes or no.
   */
  void SetTransposedSolve(bool transposed = true) { transpose = transposed; }

  /*!
   * \brief Runs the "solve" task for any rhs/sol with operator []
   * \param[in]  rhs - Right hand side of the linear system.
   * \param[out] sol - Solution of the system.
   */
  template <class T>
  void Solve(const T& rhs, T& sol) {
    if (!isfactorized) SU2_MPI::Error("The factorization has not been computed yet.", CURRENT_FUNCTION);

    /*--- Inverted logic due to CSR to CSC direct copy (PaStiX's A is our A^T). ---*/
    if (iparm[IPARM_FACTORIZATION] == PastixFactLDLT || transpose) {
      iparm[IPARM_TRANSPOSE_SOLVE] = PastixNoTrans;
    } else {
      iparm[IPARM_TRANSPOSE_SOLVE] = PastixTrans;
    }
    iparm[IPARM_VERBOSE] = PastixVerboseNot;

    for (auto i = 0ul; i < matrix.size_rhs(); ++i) workvec[i] = rhs[i];
    if (pastix_task_solve(state, matrix.size_rhs(), 1, workvec.data(), matrix.size_rhs()) != PASTIX_SUCCESS) {
      SU2_MPI::Error("Error solving linear system.", CURRENT_FUNCTION);
    }
    for (auto i = 0ul; i < matrix.size_rhs(); ++i) sol[i] = workvec[i];
  }
};
#endif
