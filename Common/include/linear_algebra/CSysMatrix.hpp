/*!
 * \file CSysMatrix.hpp
 * \brief Declaration of the block-sparse matrix class.
 *        The implemtation is in <i>CSysMatrix.cpp</i>.
 * \author F. Palacios, A. Bueno, T. Economon, P. Gomes
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

#pragma once

#include "../CConfig.hpp"
#include "CSysVector.hpp"
#include "CPastixWrapper.hpp"
#include "../toolboxes/graph_toolbox.hpp"

#include <cstdint>
#include <cstdlib>
#include <vector>
#include <cassert>

/*--- In forward mode the matrix is not of a built-in type. ---*/
#if defined(HAVE_MKL) && !defined(CODI_FORWARD_TYPE)
#include "mkl.h"
#ifndef __INTEL_MKL__
#error Could not determine the MKL version
#endif
/*--- JIT is only available since 2019. ---*/
#if __INTEL_MKL__ >= 2019
#define USE_MKL
/*---
 Lapack direct calls only seem to be created for Intel compilers, and it is not worthwhile
 making "getrf" and "getrs" compatible with AD since they are not used as often as "gemm".
---*/
#if defined(__INTEL_COMPILER) && defined(MKL_DIRECT_CALL_SEQ) && !defined(CODI_REVERSE_TYPE)
#define USE_MKL_LAPACK
#endif
template <class T>
struct mkl_jit_wrapper {
  using gemm_t = dgemm_jit_kernel_t;
  template <class... Ts>
  static void create_gemm(Ts&&... args) {
    mkl_jit_create_dgemm(args...);
  }
  static gemm_t get_gemm(void* jitter) { return mkl_jit_get_dgemm_ptr(jitter); }
};
template <>
struct mkl_jit_wrapper<float> {
  using gemm_t = sgemm_jit_kernel_t;
  template <class... Ts>
  static void create_gemm(Ts&&... args) {
    mkl_jit_create_sgemm(args...);
  }
  static gemm_t get_gemm(void* jitter) { return mkl_jit_get_sgemm_ptr(jitter); }
};
#else
#warning The current version of MKL does not support JIT gemm kernels
#endif
#endif

class CGeometry;

/*!
 * \brief Helper to communicate distributed vectors.
 * \ingroup SpLinSys
 */
struct CSysMatrixComms {
  /*!
   * \brief Routine to load a vector quantity into the data structures for MPI point-to-point
   *        communication and to launch non-blocking sends and recvs.
   * \param[in] x        - CSysVector holding the array of data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be communicated.
   */
  template <class T>
  static void Initiate(const CSysVector<T>& x, CGeometry* geometry, const CConfig* config,
                       MPI_QUANTITIES commType = MPI_QUANTITIES::SOLUTION_MATRIX);

  /*!
   * \brief Routine to complete the set of non-blocking communications launched by
   *        Initiate() and unpacking of the data in the vector.
   * \param[in] x        - CSysVector holding the array of data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] commType - Enumerated type for the quantity to be unpacked.
   */
  template <class T>
  static void Complete(CSysVector<T>& x, CGeometry* geometry, const CConfig* config,
                       MPI_QUANTITIES commType = MPI_QUANTITIES::SOLUTION_MATRIX);
};

/*!
 * \brief Reconstruct the float row-scale from a stored int8 binary exponent.
 *        The exponent \p e was packed as (e + 127) into the IEEE 754 biased-exponent field
 *        with a zero mantissa, giving an exact power of two: 2^e.
 *        This is the inverse of the encoding in EncodeQuantBlock.
 */
FORCEINLINE float DecodeQuantScale(int8_t e) noexcept {
  const uint32_t bits = static_cast<uint32_t>(std::max(0, static_cast<int>(e) + 127)) << 23;
  float scale;
  memcpy(&scale, &bits, sizeof(bits));
  return scale;
}

/*!
 * \brief Encode one nVar×nVar block into per-row int8 quantized storage.
 *        \p f(r,c) is called twice per entry (max-abs scan then encoding); it should be cheap.
 *        Stores a per-row scale exponent in \p qs and clamped int8 values in \p qv.
 */
template <class F>
FORCEINLINE void EncodeQuantBlock(const F& f, int8_t* __restrict qs, int8_t* __restrict qv,
                                  unsigned long nVar) noexcept {
  for (auto r = 0ul; r < nVar; ++r) {
    constexpr uint32_t eps_bits = 0x34000000u;
    uint32_t max_abs_bits = eps_bits;
    for (auto c = 0ul; c < nVar; ++c) {
      const float fv = SU2_TYPE::PassiveValue(f(r, c));
      uint32_t fb;
      memcpy(&fb, &fv, sizeof(fb));
      max_abs_bits = std::max(max_abs_bits, fb & 0x7FFFFFFFu);
    }
    const int e = std::min(127, std::max(-128, static_cast<int>(max_abs_bits >> 23) - 133));
    qs[r] = static_cast<int8_t>(e);
    const uint32_t inv_bits = static_cast<uint32_t>(127 - e) << 23;
    float inv_rscale;
    memcpy(&inv_rscale, &inv_bits, sizeof(inv_rscale));
    for (auto c = 0ul; c < nVar; ++c) {
      qv[r * nVar + c] =
          static_cast<int8_t>(std::max(-128.f, std::min(127.f, roundf(SU2_TYPE::PassiveValue(f(r, c)) * inv_rscale))));
    }
  }
}

/*!
 * \brief View of one matrix block, const-correct via the ScalarType template parameter.
 *        \c CBlockView<const ScalarType> is read-only; \c CBlockView<ScalarType> is mutable
 *        and exposes \c apply<Overwrite>(f) for writing with on-the-fly quantized encoding.
 *        Evaluates to \c false if the block is absent from the sparsity pattern.
 */
template <class ScalarType>
struct CBlockView {
  using QuantType = std::conditional_t<std::is_const_v<ScalarType>, const int8_t, int8_t>;

  ScalarType* ptr = nullptr;  ///< Full-precision block; non-null iff not quantized.
  QuantType* qs = nullptr;    ///< Per-row binary exponent; non-null iff quantized.
  QuantType* qv = nullptr;    ///< Quantized values (row-major); non-null iff quantized.
  unsigned long nVar = 0;

  /*! \brief False when the block is not present in the sparsity pattern. */
  explicit operator bool() const { return ptr || qs; }

  /*! \brief Return entry (row \p i, col \p j), decoding quantization if necessary. */
  std::remove_const_t<ScalarType> operator()(unsigned long i, unsigned long j) const {
    using T = std::remove_const_t<ScalarType>;
    if (ptr) return ptr[i * nVar + j];
    return static_cast<T>(qv[i * nVar + j] * DecodeQuantScale(qs[i]));
  }

  /*!
   * \brief Write the block from callable \p f(i,j).
   *        \p Overwrite=true overwrites (or quantizes for Q_LU_SGS off-diagonal blocks);
   *        \p Overwrite=false accumulates into non-quantized storage only — accumulating into
   *        quantized storage would require decode-accumulate-encode and is a silent no-op.
   *        Only enabled for mutable (non-const ScalarType) views.
   */
  template <bool Overwrite, class F, class S = ScalarType, su2enable_if<!std::is_const_v<S>> = 0>
  void apply(const F& f) const {
    if (ptr) {
      for (auto i = 0ul; i < nVar; ++i)
        for (auto j = 0ul; j < nVar; ++j) {
          if constexpr (Overwrite)
            ptr[i * nVar + j] = f(i, j);
          else
            ptr[i * nVar + j] += f(i, j);
        }
    } else if constexpr (Overwrite) {
      if (qs) EncodeQuantBlock(f, qs, qv, nVar);
    }
  }
};

/*!
 * \class CSysMatrix
 * \ingroup SpLinSys
 * \brief Main class for defining block-compressed-row-storage sparse matrices.
 */
template <class ScalarType>
class CSysMatrix {
 private:
  friend struct CSysMatrixComms;

  const int rank; /*!< \brief MPI Rank. */
  const int size; /*!< \brief MPI Size. */

  /*!< \brief Maximum number of variables the matrix can handle. The static
   * size is needed for fast, per-thread, static memory allocation. */
  enum : size_t { MAXNVAR = 20 };

  enum { OMP_MAX_SIZE_L = 8192 }; /*!< \brief Max. chunk size used in light parallel for loops. */
  enum { OMP_MAX_SIZE_H = 512 };  /*!< \brief Max. chunk size used in heavy parallel for loops. */
  enum { OMP_MIN_SIZE = 32 };     /*!< \brief Chunk size for finer grain operations. */
  unsigned long omp_light_size;   /*!< \brief Actual chunk size used in light loops (e.g. over non zeros). */
  unsigned long omp_heavy_size;   /*!< \brief Actual chunk size used in heavy loops (e.g. over rows). */
  unsigned long omp_num_parts;    /*!< \brief Number of threads used in thread-parallel LU_SGS and ILU. */
  unsigned long* omp_partitions;  /*!< \brief Point indexes of LU_SGS and ILU thread-parallel sub partitioning. */

  unsigned long nPoint;       /*!< \brief Number of points in the grid. */
  unsigned long nPointDomain; /*!< \brief Number of points in the grid (excluding halos). */
  unsigned long nVar;         /*!< \brief Number of variables (and rows of the blocks). */
  unsigned long nEqn;         /*!< \brief Number of equations (and columns of the blocks). */

  /*!
   * \brief Aggregates value arrays and sparse-structure pointers for an LDU-partitioned matrix.
   *        Each CSysMatrix holds three LDU<ScalarType> instances: the host matrix (mat), its
   *        device copy (gpu), and the ILU factorization (ilu). Ownership of the value arrays
   *        (d/l/u) and whether the pointers address host or device memory is managed by
   *        CSysMatrix. Also reused with T = QuantType to group the quantized scale/blocks
   *        storage (q_scale, q_blocks, d_q_scale, d_q_blocks) the same way; for those the pattern
   *        fields (row_ptr_l, col_ind_l, row_ptr_u, col_ind_u, nnz_l, nnz_u) are simply left
   *        unused, since the sparsity pattern is already available from mat/gpu.
   */
  template <class T>
  struct LDU {
    T* d = nullptr;                     /*!< \brief Diagonal block values. */
    T* l = nullptr;                     /*!< \brief Strictly-lower block values. */
    T* u = nullptr;                     /*!< \brief Strictly-upper block values. */
    const su2uint* row_ptr_l = nullptr; /*!< \brief Row pointers for L (geometry-owned or GPU copy). */
    const su2uint* col_ind_l = nullptr; /*!< \brief Column indices for L. */
    const su2uint* row_ptr_u = nullptr; /*!< \brief Row pointers for U. */
    const su2uint* col_ind_u = nullptr; /*!< \brief Column indices for U. */
    unsigned long nnz_l = 0;            /*!< \brief Number of L nonzeros. */
    unsigned long nnz_u = 0;            /*!< \brief Number of U nonzeros. */
  };

  LDU<ScalarType> mat;          /*!< \brief Host matrix (values owned via aligned_alloc; pattern from geometry). */
  LDU<ScalarType> gpu;          /*!< \brief Device matrix (all pointers to GPU memory). */
  LDU<ScalarType> ilu;          /*!< \brief ILU factorization, host (values owned; pattern from geometry). */
  LDU<ScalarType> gpu_ilu;      /*!< \brief ILU factorization, device (values and pattern in GPU memory). */
  ScalarType* d_invM = nullptr; /*!< \brief Device inverse diagonal blocks for the Jacobi preconditioner. */

  /*--- Quantized off-diagonal storage (used when quantized_mode == true). ---*/
  using QuantType = int8_t;

  /*! \brief Set by Initialize() when preconditioner == Q_LU_SGS, Q_JACOBI or Q_IDENTITY.
   *         mat.l and mat.u are NOT allocated; off-diagonal blocks live in q_scale/q_blocks
   *         below. Only the matrix-vector product (used by the Krylov solver and, for Q_LU_SGS,
   *         by the sweeps) reads the quantized blocks; the Jacobi preconditioner never touches
   *         them since it only applies the (full precision) inverse diagonal, and the identity
   *         preconditioner does not touch the matrix at all. */
#ifndef CODI_REVERSE_TYPE
  bool quantized_mode = false;
#else
  static constexpr bool quantized_mode = false;
#endif
  /*!< \brief Per-row exponents; .l/.u sized [nnz_l/u * nVar], .d [nPoint * nVar]
   *          (populated by QuantizeDiagonalBlocks(), always on the host, see below). */
  LDU<QuantType> q_scale;
  /*!< \brief Quantized block entries; .l/.u sized [nnz_l/u * nVar * nEqn], .d [nPoint * nVar * nEqn]. */
  LDU<QuantType> q_blocks;

  /*!< \brief Device mirrors of the quantized storage, only allocated when quantized_mode &&
   * useCuda (currently only reachable for Q_JACOBI/Q_IDENTITY, Q_LU_SGS stays host-only).
   * d_q_scale.l/.u and d_q_blocks.l/.u are plain device-side copies of q_scale.l/.u and
   * q_blocks.l/.u, uploaded *asynchronously* by HtDTransfer() so that transfer can overlap with
   * the host quantizing the diagonal; d_q_scale.d/d_q_blocks.d are that host result, uploaded (a
   * plain cudaMemcpy, not a kernel) once ready, at the end of QuantizeDiagonalBlocks(). */
  LDU<QuantType> d_q_scale;
  LDU<QuantType> d_q_blocks;

  bool useCuda = false; /*!< \brief Whether CUDA is enabled. */

  /*!< \brief Whether the inverse diagonal blocks are only needed on the device. False for the
   * Linelet preconditioner, which builds the Jacobi one but reads invM on the host. */
  bool jacobi_on_device = false;

  const su2uint* l_to_u_transp; /*!< \brief L-entry index -> U-entry index of its transpose. */
  const su2uint* u_to_l_transp; /*!< \brief U-entry index -> L-entry index of its transpose. */

  /*!
   * \brief Lookup table from edges to the L-index in the LDU split.
   * U-index == edge index by construction (edges are ordered 1:1 with the U pattern).
   * Therefore, edge_ptr_l == u_to_l_transp, but we keep a separate member for clarity.
   */
  const su2uint* edge_ptr_l;

  unsigned short ilu_fill_in; /*!< \brief Fill level for the ILU preconditioner. */

  /*!< \brief Level structure of the ILU dependency graph: rows within a level are independent,
   * rows in level k only depend on rows in levels < k. The same table drives the forward
   * (increasing level) and backward (decreasing level) substitution, because the U pattern is
   * the transpose of the L pattern. Used directly by the host/OMP substitution, and flattened
   * into ilu_level_ptr / d_ilu_level_idx below for the GPU triangular solves. */
  CCompressedSparsePatternUL levels_ilu;

  /*!< \brief Coloring of the (domain-only) ILU dependency graph, used only by the GPU iterative
   * factorization (see IluFactorColorKernel and ilu_color_ptr / d_ilu_color_idx below). The
   * host/OMP path and the GPU triangular solves use levels_ilu instead. */
  CCompressedSparsePatternUL color_ilu;

  /*!< \brief Number of colored Gauss-Seidel sweeps used to build the ILU factorization on the
   * device, see IluFactorColorKernel. Fixed (not adaptive) so the result is reproducible; set
   * from config in Initialize(). The triangular solves have no equivalent sweep count: they are
   * exact, one pass per level (see IluForwardKernel / IluBackwardKernel). */
  unsigned short ilu_gpu_sweeps = 1;

  vector<su2uint> ilu_color_ptr;      /*!< \brief Start of each color in d_ilu_color_idx, size nColors+1. */
  su2uint* d_ilu_color_idx = nullptr; /*!< \brief Row indices, grouped by color. */

  vector<su2uint> ilu_level_ptr;      /*!< \brief Start of each level in d_ilu_level_idx, size nLevels+1. */
  su2uint* d_ilu_level_idx = nullptr; /*!< \brief Row indices, grouped by level. */

  /*--- The per-color (factorization) and per-level (triangular solves) kernel launch sequences
   * are identical on every call: same grid/block sizes, same device pointers (all fixed members,
   * allocated once). Each is captured once into a CUDA graph and replayed to remove
   * host-side launch overhead without changing the parallelization. ---*/
  mutable struct CUgraphExec_st* ilu_build_graph_exec = nullptr;
  mutable struct CUgraphExec_st* ilu_apply_graph_exec = nullptr;
  mutable const ScalarType* ilu_apply_graph_vec = nullptr; /*!< \brief Pointers the apply graph
                                                            * was captured with, to detect when
                                                            * it must be recaptured. */
  mutable ScalarType* ilu_apply_graph_prod = nullptr;
  /*--- The legacy default stream cannot be captured into a graph. ---*/
  mutable struct CUstream_st* ilu_stream = nullptr;

  ScalarType* invM; /*!< \brief Inverse of (Jacobi) preconditioner. */

  /*--- Temporary (hence mutable) working memory used in the Linelet preconditioner, outer vector is for threads ---*/
  mutable vector<vector<const ScalarType*>>
      LineletUpper; /*!< \brief Pointers to the upper blocks of the tri-diag system (working memory). */
  mutable vector<vector<ScalarType>>
      LineletInvDiag; /*!< \brief Inverse of the diagonal blocks of the tri-diag system (working memory). */
  mutable vector<vector<ScalarType>>
      LineletVector; /*!< \brief Solution and RHS of the tri-diag system (working memory). */

#ifdef USE_MKL
  using gemm_t = typename mkl_jit_wrapper<ScalarType>::gemm_t;
  void* MatrixMatrixProductJitter;               /*!< \brief Jitter handle for MKL JIT based GEMM. */
  gemm_t MatrixMatrixProductKernel;              /*!< \brief MKL JIT based GEMM kernel. */
  void* MatrixVectorProductJitterBetaZero;       /*!< \brief Jitter handle for MKL JIT based GEMV. */
  gemm_t MatrixVectorProductKernelBetaZero;      /*!< \brief MKL JIT based GEMV kernel. */
  void* MatrixVectorProductJitterBetaOne;        /*!< \brief Jitter handle for MKL JIT based GEMV with BETA=1.0. */
  gemm_t MatrixVectorProductKernelBetaOne;       /*!< \brief MKL JIT based GEMV kernel with BETA=1.0. */
  void* MatrixVectorProductJitterAlphaMinusOne;  /*!< \brief Jitter handle for MKL JIT based GEMV with ALPHA=-1.0 and
                                                    BETA=1.0. */
  gemm_t MatrixVectorProductKernelAlphaMinusOne; /*!< \brief MKL JIT based GEMV kernel with ALPHA=-1.0 and BETA=1.0. */
#endif

#ifdef HAVE_PASTIX
  mutable CPastixWrapper<ScalarType> pastix_wrapper;
#endif

  /*!
   * \brief Handle type conversion for when we Set, Add, etc. blocks, preserving derivative information (if supported by
   * types).
   */
  template <class DstType, class SrcType, su2enable_if<std::is_arithmetic<DstType>::value> = 0>
  FORCEINLINE static DstType ActiveAssign(const SrcType& val) {
    return SU2_TYPE::GetValue(val);
  }

  template <class DstType, class SrcType, su2enable_if<!std::is_arithmetic<DstType>::value> = 0>
  FORCEINLINE static DstType ActiveAssign(const SrcType& val) {
    return val;
  }

  /*!
   * \brief Handle type conversion for when we Set, Add, etc. blocks, discarding derivative information.
   */
  template <class SrcType>
  FORCEINLINE static ScalarType PassiveAssign(const SrcType& val) {
    return SU2_TYPE::GetValue(val);
  }

  /*!
   * \brief Calculates the matrix-vector product: product = matrix*vector
   * \param[in] matrix
   * \param[in] vector
   * \param[out] product
   */
  void MatrixVectorProduct(const ScalarType* matrix, const ScalarType* vector, ScalarType* product) const;

  /*!
   * \brief Calculates the matrix-vector product: product += matrix*vector
   * \param[in] matrix
   * \param[in] vector
   * \param[in,out] product
   */
  void MatrixVectorProductAdd(const ScalarType* matrix, const ScalarType* vector, ScalarType* product) const;

  /*!
   * \brief Calculates the matrix-vector product: product -= matrix*vector
   * \param[in] matrix
   * \param[in] vector
   * \param[in,out] product
   */
  void MatrixVectorProductSub(const ScalarType* matrix, const ScalarType* vector, ScalarType* product) const;

  /*!
   * \brief Calculates the matrix-matrix product
   */
  void MatrixMatrixProduct(const ScalarType* matrix_a, const ScalarType* matrix_b, ScalarType* product) const;

  /*!
   * \brief Subtract b from a and store the result in c.
   */
  FORCEINLINE void VectorSubtraction(const ScalarType* a, const ScalarType* b, ScalarType* c) const {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) c[iVar] = a[iVar] - b[iVar];
  }

  /*!
   * \brief Subtract b from a and store the result in c.
   */
  FORCEINLINE void MatrixSubtraction(const ScalarType* a, const ScalarType* b, ScalarType* c) const {
    SU2_OMP_SIMD
    for (unsigned long iVar = 0; iVar < nVar * nEqn; iVar++) c[iVar] = a[iVar] - b[iVar];
  }

  /*!
   * \brief Copy matrix src into dst, transpose if required.
   */
  FORCEINLINE void MatrixCopy(const ScalarType* src, ScalarType* dst) const {
    SU2_OMP_SIMD
    for (auto iVar = 0ul; iVar < nVar * nEqn; ++iVar) dst[iVar] = src[iVar];
  }

  /*!
   * \brief Zero a matrix.
   */
  FORCEINLINE void ZeroMatrix(ScalarType* mat) const {
    SU2_OMP_SIMD
    for (auto iVar = 0ul; iVar < nVar * nEqn; ++iVar) mat[iVar] = 0;
  }

  /*!
   * \brief Solve a small (nVar x nVar) linear system using Gaussian elimination.
   * \param[in,out] matrix - On entry the system matrix, on exit the factorized matrix.
   * \param[in,out] vec - On entry the rhs, on exit the solution.
   */
  void GaussElimination(ScalarType* matrix, ScalarType* vec) const;

  /*!
   * \brief Invert a small dense matrix.
   * \param[in,out] matrix - On entry the system matrix, on exit the factorized matrix.
   * \param[out] inverse - the matrix inverse.
   */
  void MatrixInverse(ScalarType* matrix, ScalarType* inverse) const;

  /*!
   * \brief Performs the Gauss Elimination algorithm to solve the linear subsystem of the (i,i) subblock and rhs.
   * \param[in] block_i - Index of the (i,i) diagonal block.
   * \param[in] rhs - Right-hand-side of the linear system.
   * \return Solution of the linear system (overwritten on rhs).
   */
  inline void GaussElimination(unsigned long block_i, ScalarType* rhs) const;

  /*!
   * \brief Inverse diagonal block.
   * \param[in] block_i - Indexes of the block in the matrix-by-blocks structure.
   * \param[out] invBlock - Inverse block.
   */
  inline void InverseDiagonalBlock(unsigned long block_i, ScalarType* invBlock) const;

  /*!
   * \brief Invert diagonal block (Uii) of the ILU matrix in place.
   * \param[in] block_i - Index of the block to invert.
   * \return Inverted block.
   */
  inline const ScalarType* InvertDiagonalBlockILUMatrix(unsigned long block_i);

  /*!
   * \brief Returns the start of the ILU block or nullptr if (i,j) is not a nonzero.
   * \param[in] block_i/j - Indexes of the block in the matrix-by-blocks structure.
   */
  inline ScalarType* GetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j);

  /*!
   * \brief Performs the product of i-th row of the upper part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the upper part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \param[in] col_ub - Exclusive upper bound for column indices considered in multiplication.
   * \param[out] prod - Result of the product U(A)*vec.
   */
  inline void UpperProduct(const CSysVector<ScalarType>& vec, unsigned long row_i, unsigned long col_ub,
                           ScalarType* prod) const;

  /*!
   * \brief Performs the product of i-th row of the lower part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the lower part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \param[in] col_lb - Inclusive lower bound for column indices considered in multiplication.
   * \param[out] prod - Result of the product L(A)*vec.
   */
  inline void LowerProduct(const CSysVector<ScalarType>& vec, unsigned long row_i, unsigned long col_lb,
                           ScalarType* prod) const;

  /*!
   * \brief Performs the product of i-th row of the diagonal part of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the diagonal part of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return prod Result of the product D(A)*vec (stored at *prod_row_vector).
   */
  inline void DiagonalProduct(const CSysVector<ScalarType>& vec, unsigned long row_i, ScalarType* prod) const;

  /*!
   * \brief Performs the product of i-th row of a sparse matrix by a vector.
   * \param[in] vec - Vector to be multiplied by the row of the sparse matrix A.
   * \param[in] row_i - Row of the matrix to be multiplied by vector vec.
   * \return Result of the product (stored at *prod_row_vector).
   */
  void RowProduct(const CSysVector<ScalarType>& vec, unsigned long row_i, ScalarType* prod) const;

  /*!
   * \brief Computes product += A_k * vec using the quantized representation of block k.
   * \note Only valid after QuantizeDiagonalBlocks() has been called.
   * \param[in] k - Block index in the CSR flat storage.
   * \param[in] vec - Input vector (nEqn entries).
   * \param[in,out] prod - Accumulation output (nVar entries).
   */
  inline void QuantizedMatVecAdd(const QuantType* qs, const QuantType* qv, const ScalarType* vec,
                                 ScalarType* prod) const;

  /*! \brief Quantize one nVar×nVar block (row-major) into the int8 scale+value arrays.
   *         Called on the hot assembly path (SetBlocks/UpdateBlocks in Q_LU_SGS mode). */
  void QuantizeBlock(const ScalarType* blk, QuantType* qs, QuantType* qv) const;

  /*! \brief Full-row product using quantized L/D/U (Q_LU_SGS SpMV path). */
  inline void QuantizedRowProduct(const CSysVector<ScalarType>& vec, unsigned long row_i, ScalarType* prod) const;

  /*! \brief Upper-triangle product using quantized U (Q_LU_SGS backward sweep). */
  inline void QuantizedUpperProduct(const CSysVector<ScalarType>& vec, unsigned long row_i, unsigned long col_ub,
                                    ScalarType* prod) const;

  /*! \brief Lower-triangle product using quantized L (Q_LU_SGS forward sweep). */
  inline void QuantizedLowerProduct(const CSysVector<ScalarType>& vec, unsigned long row_i, unsigned long col_lb,
                                    ScalarType* prod) const;

  /*! \brief Diagonal product using quantized D (Q_LU_SGS backward sweep). */
  inline void QuantizedDiagonalProduct(const CSysVector<ScalarType>& vec, unsigned long row_i, ScalarType* prod) const;

  /*! \brief Gauss elimination on the quantized diagonal block: decodes q_blocks.d into a local
   *         ScalarType buffer and delegates to the scalar GaussElimination overload. */
  inline void QuantizedGaussElimination(unsigned long block_i, ScalarType* rhs) const;

  /*--- Hooks for GPU versions (implemented is in CSysMatrixGPU.cu). ---*/

  /*!
   * \brief Performs the product of a sparse matrix by a CSysVector on the device.
   */
  void MatrixVectorProductGPU(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod, CGeometry* geometry,
                              const CConfig* config) const;

  /*!
   * \brief Build the Jacobi preconditioner on the device, from the device copy of the matrix.
   * \note Requires the device matrix to be up to date, see HtDTransfer.
   */
  void BuildJacobiPreconditionerGPU();

  /*!
   * \brief Apply the Jacobi preconditioner on the GPU/device side.
   */
  void ComputeJacobiPreconditionerGPU(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                      CGeometry* geometry, const CConfig* config) const;

  /*!
   * \brief Build the ILU preconditioner on the device, from the device copy of the matrix.
   * \note Requires the device matrix to be up to date, see HtDTransfer.
   */
  void BuildILUPreconditionerGPU();

  /*!
   * \brief Apply the ILU preconditioner on the device.
   */
  void ComputeILUPreconditionerGPU(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod) const;

 public:
  /*!
   * \brief Constructor of the class.
   */
  CSysMatrix();

  /*!
   * \brief Destructor of the class.
   */
  ~CSysMatrix();

  /*!
   * \brief Initializes the sparse matrix.
   * \note The preconditioners require nVar == nEqn (square blocks).
   * \param[in] npoint - Number of points including halos.
   * \param[in] npointdomain - Number of points excluding halos.
   * \param[in] nvar - Number of variables (and rows of the blocks).
   * \param[in] neqn - Number of equations (and columns of the blocks).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] needTranspPtr - If the L/U transpose maps should be built, used for "SetDiagonalAsColumnSum".
   * \param[in] grad_mode - Gradient smoothing mode, only used to detect the right preconditioner type.
   * \param[in] allow_quant - Quantization is only possible with solvers that "set and forget" the off-diagonal
   *            blocks of the matrix. Solvers that perform multiple updates would lose too much information, so
   *            that pattern is not supported with quantization (the code will hit null pointers). It is up to
   *            the solver to declare whether it will "set and forget".
   */
  void Initialize(unsigned long npoint, unsigned long npointdomain, unsigned short nvar, unsigned short neqn,
                  bool EdgeConnect, CGeometry* geometry, const CConfig* config, bool needTranspPtr = false,
                  bool grad_mode = false, bool allow_quant = false);

  /*!
   * \brief Compresses off-diagonal blocks into quantized form for use with USE_QUANTIZATION.
   */
  void QuantizeDiagonalBlocks();

  /*!
   * \brief Sets to zero all the entries of the sparse matrix.
   */
  void SetValZero();

  /*!
   * \brief Sets to zero all the block diagonal entries of the sparse matrix.
   */
  void SetValDiagonalZero();

  /*!
   * \brief Performs the memory copy from host to device.
   * \param[in] trigger - boolean value that decides whether to conduct the transfer or not. True by default.
   */
  void HtDTransfer(bool trigger = true) const;

  /*!
   * \brief Get a pointer to the start of block "ij"
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \return Pointer to location in memory where the block starts.
   */
  FORCEINLINE const ScalarType* GetBlock(unsigned long block_i, unsigned long block_j) const {
    if (block_i == block_j) return &mat.d[block_i * nVar * nEqn];
    if (block_j < block_i) {
      for (auto index = mat.row_ptr_l[block_i]; index < mat.row_ptr_l[block_i + 1]; ++index)
        if (mat.col_ind_l[index] == block_j) return &mat.l[index * nVar * nEqn];
      return nullptr;
    }
    for (auto index = mat.row_ptr_u[block_i]; index < mat.row_ptr_u[block_i + 1]; ++index)
      if (mat.col_ind_u[index] == block_j) return &mat.u[index * nVar * nEqn];
    return nullptr;
  }

  /*!
   * \brief Get a pointer to the start of block "ij", non-const version.
   */
  FORCEINLINE ScalarType* GetBlock(unsigned long block_i, unsigned long block_j) {
    const CSysMatrix& const_this = *this;
    return const_cast<ScalarType*>(const_this.GetBlock(block_i, block_j));
  }

  /*!
   * \brief Read-only view of block (block_i, block_j). In Q_LU_SGS mode values are decoded
   *        on access inside CBlockView::operator()(i,j); no temporary copy is made.
   * \return A CBlockView<const ScalarType> that evaluates to false if the block is absent.
   */
  FORCEINLINE CBlockView<const ScalarType> GetBlockView(unsigned long block_i, unsigned long block_j) const {
#define GET_BLOCK_VIEW_IMPL                                                                                        \
  if (!quantized_mode || block_i == block_j) {                                                                     \
    return {GetBlock(block_i, block_j), nullptr, nullptr, nVar};                                                   \
  }                                                                                                                \
  if (block_j < block_i) {                                                                                         \
    for (auto k = mat.row_ptr_l[block_i]; k < mat.row_ptr_l[block_i + 1]; ++k)                                     \
      if (mat.col_ind_l[k] == block_j) return {nullptr, &q_scale.l[k * nVar], &q_blocks.l[k * nVar * nVar], nVar}; \
  } else {                                                                                                         \
    for (auto k = mat.row_ptr_u[block_i]; k < mat.row_ptr_u[block_i + 1]; ++k)                                     \
      if (mat.col_ind_u[k] == block_j) return {nullptr, &q_scale.u[k * nVar], &q_blocks.u[k * nVar * nVar], nVar}; \
  }                                                                                                                \
  return {}
    GET_BLOCK_VIEW_IMPL;
  }

  /*!
   * \overload Non const version of GetBlockView.
   */
  FORCEINLINE CBlockView<ScalarType> GetBlockView(unsigned long block_i, unsigned long block_j) {
    GET_BLOCK_VIEW_IMPL;
#undef GET_BLOCK_VIEW_IMPL
  }

  /*!
   * \brief Set the value of a scaled block in the sparse matrix.
   * \note This is an templated overload for C2Dcontainer specialization su2matrix.
   *       It assumes that MatrixType supports a member type Scalar and access operator(i, j).
   *       If the template param Overwrite is false we add to the block (bij += alpha*b).
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to set to A(i, j).
   * \param[in] alpha - Scale factor.
   */
  template <bool Overwrite = true, class MatrixType>
  inline void SetBlock(unsigned long block_i, unsigned long block_j, MatrixType& val_block,
                       std::decay_t<typename MatrixType::Scalar> alpha = 1.0) {
    auto view = GetBlockView(block_i, block_j);
    if (!view) return;
    view.template apply<Overwrite>(
        [&](unsigned long i, unsigned long j) { return PassiveAssign(alpha * val_block(i, j)); });
  }

  /*!
   * \overload val_block is a pointer instead of a matrix type.
   */
  template <bool Overwrite = true, class OtherType, su2enable_if<!is_pointer<OtherType>::value> = 0>
  inline void SetBlock(unsigned long block_i, unsigned long block_j, const OtherType* val_block,
                       std::decay_t<OtherType> alpha = 1.0) {
    auto view = GetBlockView(block_i, block_j);
    if (!view) return;
    view.template apply<Overwrite>(
        [&](unsigned long i, unsigned long j) { return PassiveAssign(alpha * val_block[i * nEqn + j]); });
  }

  /*!
   * \overload val_block is a double pointer instead of matrix type.
   */
  template <bool Overwrite = true, class OtherType>
  inline void SetBlock(unsigned long block_i, unsigned long block_j, const OtherType* const* val_block,
                       std::decay_t<OtherType> alpha = 1.0) {
    auto view = GetBlockView(block_i, block_j);
    if (!view) return;
    view.template apply<Overwrite>(
        [&](unsigned long i, unsigned long j) { return PassiveAssign(alpha * val_block[i][j]); });
  }

  /*!
   * \brief Add a scaled block (in flat format) to the sparse matrix (see SetBlock).
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to set to A(i, j).
   * \param[in] alpha - Scale factor.
   */
  template <class T, class OtherType = ScalarType>
  inline void AddBlock(unsigned long block_i, unsigned long block_j, const T& val_block, OtherType alpha = 1.0) {
    SetBlock<false>(block_i, block_j, val_block, alpha);
  }

  /*!
   * \brief Subtracts the specified block to the sparse matrix (see AddBlock).
   * \param[in] block_i - Row index.
   * \param[in] block_j - Column index.
   * \param[in] val_block - Block to subtract to A(i, j).
   */
  template <class T>
  inline void SubtractBlock(unsigned long block_i, unsigned long block_j, const T& val_block) {
    AddBlock(block_i, block_j, val_block, -1);
  }

  /*!
   * \brief Returns the 4 blocks ii, ij, ji, jj used by "UpdateBlocks".
   * \note This method assumes an FVM-type sparse pattern.
   * \param[in] edge - Index of edge that connects iPoint and jPoint.
   * \param[in] iPoint - Row to which we add the blocks.
   * \param[in] jPoint - Row from which we subtract the blocks.
   * \param[out] bii, bij, bji, bjj - Blocks of the matrix.
   */
  inline void GetBlocks(unsigned long iEdge, unsigned long iPoint, unsigned long jPoint, ScalarType*& bii,
                        ScalarType*& bij, ScalarType*& bji, ScalarType*& bjj) {
    const auto blkSz = nVar * nEqn;
    bii = &mat.d[iPoint * blkSz];
    bjj = &mat.d[jPoint * blkSz];
    bij = &mat.u[iEdge * blkSz];
    bji = &mat.l[edge_ptr_l[iEdge] * blkSz];
  }

  /*!
   * \brief Update 4 blocks ii, ij, ji, jj (add to i* sub from j*).
   * \note This method assumes an FVM-type sparse pattern.
   * \param[in] edge - Index of edge that connects iPoint and jPoint.
   * \param[in] iPoint - Row to which we add the blocks.
   * \param[in] jPoint - Row from which we subtract the blocks.
   * \param[in] block_i - Adds to ii, subs from ji.
   * \param[in] block_j - Adds to ij, subs from jj.
   * \param[in] scale - Scale blocks during update (axpy type op).
   */
  template <bool OverwriteOffDiag = false, class MatrixType, class OtherType = ScalarType>
  inline void UpdateBlocks(unsigned long iEdge, unsigned long iPoint, unsigned long jPoint, const MatrixType& block_i,
                           const MatrixType& block_j, OtherType scale = 1) {
    const auto blkSz = nVar * nEqn;
    auto* bii = &mat.d[iPoint * blkSz];
    auto* bjj = &mat.d[jPoint * blkSz];

    unsigned long iVar, jVar, offset = 0;

    if (quantized_mode) {
      assert(OverwriteOffDiag);
      /*--- Diagonal: full-precision accumulation. Off-diagonal: quantize on the fly. ---*/
      ScalarType bij_buf[MAXNVAR * MAXNVAR], bji_buf[MAXNVAR * MAXNVAR];
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++, ++offset) {
          bii[offset] += PassiveAssign(block_i[iVar][jVar] * scale);
          bjj[offset] -= PassiveAssign(block_j[iVar][jVar] * scale);
          bij_buf[offset] = PassiveAssign(block_j[iVar][jVar] * scale);
          bji_buf[offset] = -PassiveAssign(block_i[iVar][jVar] * scale);
        }
      QuantizeBlock(bij_buf, &q_scale.u[iEdge * nVar], &q_blocks.u[iEdge * blkSz]);
      const auto k_l = edge_ptr_l[iEdge];
      QuantizeBlock(bji_buf, &q_scale.l[k_l * nVar], &q_blocks.l[k_l * blkSz]);
      return;
    }

    auto* bij = &mat.u[iEdge * blkSz];
    auto* bji = &mat.l[edge_ptr_l[iEdge] * blkSz];
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nEqn; jVar++) {
        bii[offset] += PassiveAssign(block_i[iVar][jVar] * scale);
        bjj[offset] -= PassiveAssign(block_j[iVar][jVar] * scale);
        if constexpr (OverwriteOffDiag) {
          bij[offset] = PassiveAssign(block_j[iVar][jVar] * scale);
          bji[offset] = -PassiveAssign(block_i[iVar][jVar] * scale);
        } else {
          bij[offset] += PassiveAssign(block_j[iVar][jVar] * scale);
          bji[offset] -= PassiveAssign(block_i[iVar][jVar] * scale);
        }
        ++offset;
      }
    }
  }

  /*!
   * \brief Short-hand for the "subtractive" version (sub from i* add to j*) of UpdateBlocks.
   */
  template <class MatrixType>
  inline void UpdateBlocksSub(unsigned long iEdge, unsigned long iPoint, unsigned long jPoint,
                              const MatrixType& block_i, const MatrixType& block_j) {
    UpdateBlocks<false, MatrixType, ScalarType>(iEdge, iPoint, jPoint, block_i, block_j, -1);
  }

  /*!
   * \brief SIMD version, does the update for multiple edges and points.
   * \note Nothing is updated if the mask is 0.
   */
  template <class MatTypeSIMD, size_t N, class I, class F = ScalarType>
  FORCEINLINE void SetBlocks(simd::Array<I, N> iEdge, simd::Array<I, N> iPoint, simd::Array<I, N> jPoint,
                             const MatTypeSIMD& block_i, const MatTypeSIMD& block_j, simd::Array<F, N> mask = 1) {
    static_assert(MatTypeSIMD::StaticSize, "This method requires static size blocks.");
    static_assert(MatTypeSIMD::IsRowMajor, "Block storage is not compatible with matrix.");
    constexpr size_t blkSz = MatTypeSIMD::StaticSize;
    assert(blkSz == nVar * nEqn);

    /*--- "Transpose" the blocks, scale, and possibly convert types,
     * giving the compiler the chance to vectorize all of these. ---*/
    ScalarType blk_i[N][blkSz], blk_j[N][blkSz];

    for (size_t i = 0; i < blkSz; ++i) {
      SU2_OMP_SIMD_IF_NOT_AD
      for (size_t k = 0; k < N; ++k) {
        blk_i[k][i] = PassiveAssign(-mask[k] * block_i.data()[i][k]);
        blk_j[k][i] = PassiveAssign(mask[k] * block_j.data()[i][k]);
      }
    }

    /*--- Update one by one skipping if mask is 0. ---*/
    for (size_t k = 0; k < N; ++k) {
      if (mask[k] == 0) continue;

      auto bii = &mat.d[iPoint[k] * blkSz];
      auto bjj = &mat.d[jPoint[k] * blkSz];

      if (quantized_mode) {
        SU2_OMP_SIMD
        for (size_t i = 0; i < blkSz; ++i) {
          bii[i] -= blk_i[k][i];
          bjj[i] -= blk_j[k][i];
        }
        QuantizeBlock(blk_j[k], &q_scale.u[iEdge[k] * nVar], &q_blocks.u[iEdge[k] * blkSz]);
        const auto k_l = edge_ptr_l[iEdge[k]];
        QuantizeBlock(blk_i[k], &q_scale.l[k_l * nVar], &q_blocks.l[k_l * blkSz]);
      } else {
        auto bij = &mat.u[iEdge[k] * blkSz];
        auto bji = &mat.l[edge_ptr_l[iEdge[k]] * blkSz];
        /*--- Update, block i was negated during transpose in the
         * hope the assignments below become non-temporal stores. ---*/
        SU2_OMP_SIMD
        for (size_t i = 0; i < blkSz; ++i) {
          bii[i] -= blk_i[k][i];
          bjj[i] -= blk_j[k][i];
          bij[i] = blk_j[k][i];
          bji[i] = blk_i[k][i];
        }
      }
    }
  }

  /*!
   * \brief Sets 2 blocks ij and ji (add to i* sub from j*) associated with
   *        one edge of an FVM-type sparse pattern.
   * \note The parameter Overwrite allows completely writing over the
   *       current values held by the matrix (true), or updating them (false).
   * \param[in] edge - Index of edge that connects iPoint and jPoint.
   * \param[in] block_i - Subs from ji.
   * \param[in] block_j - Adds to ij.
   * \param[in] scale - Scale blocks during update (axpy type op).
   */
  template <class MatrixType, class OtherType = ScalarType, bool Overwrite = true>
  inline void SetBlocks(unsigned long iEdge, const MatrixType& block_i, const MatrixType& block_j,
                        OtherType scale = 1) {
    const auto blkSz = nVar * nEqn;
    unsigned long iVar, jVar, offset = 0;

    if (quantized_mode) {
      assert(Overwrite);
      ScalarType bij_buf[MAXNVAR * MAXNVAR], bji_buf[MAXNVAR * MAXNVAR];
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nEqn; jVar++, ++offset) {
          bij_buf[offset] = PassiveAssign(block_j[iVar][jVar] * scale);
          bji_buf[offset] = -PassiveAssign(block_i[iVar][jVar] * scale);
        }
      QuantizeBlock(bij_buf, &q_scale.u[iEdge * nVar], &q_blocks.u[iEdge * blkSz]);
      const auto k_l = edge_ptr_l[iEdge];
      QuantizeBlock(bji_buf, &q_scale.l[k_l * nVar], &q_blocks.l[k_l * blkSz]);
      return;
    }

    ScalarType* bij = &mat.u[iEdge * blkSz];
    ScalarType* bji = &mat.l[edge_ptr_l[iEdge] * blkSz];
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nEqn; jVar++) {
        bij[offset] = (Overwrite ? ScalarType(0) : bij[offset]) + PassiveAssign(block_j[iVar][jVar] * scale);
        bji[offset] = (Overwrite ? ScalarType(0) : bji[offset]) - PassiveAssign(block_i[iVar][jVar] * scale);
        ++offset;
      }
    }
  }

  /*!
   * \brief Short-hand for the "additive overwrite" version of SetBlocks.
   */
  template <class MatrixType, class OtherType = ScalarType>
  inline void UpdateBlocks(unsigned long iEdge, const MatrixType& block_i, const MatrixType& block_j,
                           OtherType scale = 1) {
    SetBlocks<MatrixType, OtherType, false>(iEdge, block_i, block_j, scale);
  }

  /*!
   * \brief Short-hand for the "subtractive" version (sub from i* add to j*) of SetBlocks.
   */
  template <class MatrixType>
  inline void UpdateBlocksSub(unsigned long iEdge, const MatrixType& block_i, const MatrixType& block_j) {
    SetBlocks<MatrixType, ScalarType, false>(iEdge, block_i, block_j, -1);
  }

  /*!
   * \brief SIMD version, does the update for multiple edges.
   * \note Nothing is updated if the mask is 0.
   */
  template <class MatTypeSIMD, size_t N, class I, class F = ScalarType>
  FORCEINLINE void SetBlocks(simd::Array<I, N> iEdge, const MatTypeSIMD& block_i, const MatTypeSIMD& block_j,
                             simd::Array<F, N> mask = 1) {
    static_assert(MatTypeSIMD::StaticSize, "This method requires static size blocks.");
    static_assert(MatTypeSIMD::IsRowMajor, "Block storage is not compatible with matrix.");
    constexpr size_t blkSz = MatTypeSIMD::StaticSize;
    assert(blkSz == nVar * nEqn);

    /*--- "Transpose" the blocks, scale, and possibly convert types,
     * giving the compiler the chance to vectorize all of these. ---*/
    ScalarType blk_i[N][blkSz], blk_j[N][blkSz];

    for (size_t i = 0; i < blkSz; ++i) {
      SU2_OMP_SIMD_IF_NOT_AD
      for (size_t k = 0; k < N; ++k) {
        blk_i[k][i] = PassiveAssign(-mask[k] * block_i.data()[i][k]);
        blk_j[k][i] = PassiveAssign(mask[k] * block_j.data()[i][k]);
      }
    }

    /*--- Update one by one skipping if mask is 0. ---*/
    for (size_t k = 0; k < N; ++k) {
      if (mask[k] == 0) continue;

      if (quantized_mode) {
        QuantizeBlock(blk_j[k], &q_scale.u[iEdge[k] * nVar], &q_blocks.u[iEdge[k] * blkSz]);
        const auto k_l = edge_ptr_l[iEdge[k]];
        QuantizeBlock(blk_i[k], &q_scale.l[k_l * nVar], &q_blocks.l[k_l * blkSz]);
      } else {
        ScalarType* bij = &mat.u[iEdge[k] * blkSz];
        ScalarType* bji = &mat.l[edge_ptr_l[iEdge[k]] * blkSz];
        /*--- Update, block i was negated during transpose in the
         * hope the assignments below become non-temporal stores. ---*/
        SU2_OMP_SIMD
        for (size_t i = 0; i < blkSz; ++i) {
          bij[i] = blk_j[k][i];
          bji[i] = blk_i[k][i];
        }
      }
    }
  }

  /*!
   * \brief Sets the specified block to the (i, i) subblock of the sparse matrix.
   *        Scales the input block by factor alpha. If the Overwrite parameter is
   *        false we update instead (bii += alpha*b).
   * \param[in] block_i - Diagonal index.
   * \param[in] val_block - Block to add to the diagonal of the matrix.
   * \param[in] alpha - Scale factor.
   */
  template <class OtherType, bool Overwrite = true, class T = ScalarType>
  inline void SetBlock2Diag(unsigned long block_i, const OtherType& val_block, T alpha = 1.0) {
    auto mat_ii = &mat.d[block_i * nVar * nEqn];

    for (auto iVar = 0ul; iVar < nVar; iVar++)
      for (auto jVar = 0ul; jVar < nEqn; jVar++) {
        *mat_ii = (Overwrite ? ScalarType(0) : *mat_ii) + PassiveAssign(alpha * val_block[iVar][jVar]);
        ++mat_ii;
      }
  }

  /*!
   * \brief Non overwrite version of SetBlock2Diag, also with scaling.
   */
  template <class OtherType, class T = ScalarType>
  inline void AddBlock2Diag(unsigned long block_i, const OtherType& val_block, T alpha = 1.0) {
    SetBlock2Diag<OtherType, false>(block_i, val_block, alpha);
  }

  /*!
   * \brief Short-hand to AddBlock2Diag with alpha = -1, i.e. subtracts from the current diagonal.
   */
  template <class OtherType>
  inline void SubtractBlock2Diag(unsigned long block_i, const OtherType& val_block) {
    AddBlock2Diag(block_i, val_block, -1.0);
  }

  /*!
   * \brief Adds the specified value to the diagonal of the (i, i) subblock
   *        of the matrix-by-blocks structure.
   * \param[in] block_i - Diagonal index.
   * \param[in] val_matrix - Value to add to the diagonal elements of A(i, i).
   */
  template <class OtherType>
  inline void AddVal2Diag(unsigned long block_i, OtherType val_matrix) {
    auto d = &mat.d[block_i * nVar * nVar];
    for (auto iVar = 0ul; iVar < nVar; iVar++) d[iVar * (nVar + 1)] += PassiveAssign(val_matrix);
  }

  /*!
   * \brief Adds the specified value to the diagonal of the (i, i) subblock
   *        of the matrix-by-blocks structure.
   * \param[in] block_i - Diagonal index.
   * \param[in] iVar - Variable index.
   * \param[in] val - Value to add to the diagonal elements of A(i, i).
   */
  template <class OtherType>
  inline void AddVal2Diag(unsigned long block_i, unsigned long iVar, OtherType val) {
    mat.d[block_i * nVar * nVar + iVar * (nVar + 1)] += PassiveAssign(val);
  }

  /*!
   * \brief Sets the specified value to the diagonal of the (i, i) subblock
   *        of the matrix-by-blocks structure.
   * \param[in] block_i - Diagonal index.
   * \param[in] val_matrix - Value to add to the diagonal elements of A(i, i).
   */
  template <class OtherType>
  inline void SetVal2Diag(unsigned long block_i, OtherType val_matrix) {
    /*--- Clear entire block before setting its diagonal. ---*/
    SU2_OMP_SIMD
    for (auto iVar = 0ul; iVar < nVar * nVar; iVar++) mat.d[block_i * nVar * nVar + iVar] = 0.0;

    AddVal2Diag(block_i, val_matrix);
  }

  /*!
   * \brief Deletes the values of a row of the sparse matrix.
   * \param[in] block_i - Index of the block.
   * \param[in] row - Row within the block.
   */
  void DeleteValsRowi(unsigned long block_i, unsigned long row);

  /*!
   * \brief Modifies this matrix (A) and a rhs vector (b) such that (A^-1 * b)_i = x_i.
   * \param[in] node_i - Index of the node for which to enforce the solution of all DOF's.
   * \param[in] x_i - Values to enforce (nVar sized).
   * \param[in,out] b - The rhs vector (b := b - A_{*,i} * x_i;  b_i = x_i).
   */
  template <class OtherType>
  void EnforceSolutionAtNode(unsigned long node_i, const OtherType* x_i, CSysVector<OtherType>& b);

  /*!
   * \brief Similar to EnforceSolutionAtNode, but for 0 projection in a given direction.
   */
  template <class OtherType>
  void EnforceZeroProjection(unsigned long node_i, const OtherType* n, CSysVector<OtherType>& b);

  /*!
   * \brief Sets the diagonal entries of the matrix as the sum of the blocks in the corresponding column.
   */
  void SetDiagonalAsColumnSum();

  /*!
   * \brief Transposes the matrix, any preconditioner that was computed may be invalid.
   */
  void TransposeInPlace();

  /*!
   * \brief Add a scaled sparse matrix to "this" (axpy-type operation, A = A+alpha*B).
   * \note Matrices must have the same sparse pattern.
   * \param[in] alpha - The scaling constant.
   * \param[in] B - Matrix being.
   */
  void MatrixMatrixAddition(ScalarType alpha, const CSysMatrix& B);

  /*!
   * \brief Performs the product of a sparse matrix by a CSysVector.
   * \param[in] vec - CSysVector to be multiplied by the sparse matrix A.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[out] prod - Result of the product.
   */
  void MatrixVectorProduct(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod, CGeometry* geometry,
                           const CConfig* config) const;

  /*!
   * \brief Build the Jacobi preconditioner.
   */
  void BuildJacobiPreconditioner();

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeJacobiPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod, CGeometry* geometry,
                                   const CConfig* config) const;

  /*!
   * \brief Build the ILU preconditioner.
   */
  void BuildILUPreconditioner();

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeILUPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod, CGeometry* geometry,
                                const CConfig* config) const;

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   */
  void ComputeLU_SGSPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod, CGeometry* geometry,
                                   const CConfig* config) const;

  /*!
   * \brief Build the Linelet preconditioner.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void BuildLineletPreconditioner(const CGeometry* geometry, const CConfig* config);

  /*!
   * \brief Multiply CSysVector by the preconditioner
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product A*vec.
   */
  void ComputeLineletPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod,
                                    CGeometry* geometry, const CConfig* config) const;

  /*!
   * \brief Compute the linear residual.
   * \param[in] sol - Solution (x).
   * \param[in] f - Right hand side (b).
   * \param[out] res - Residual (Ax-b).
   */
  void ComputeResidual(const CSysVector<ScalarType>& sol, const CSysVector<ScalarType>& f,
                       CSysVector<ScalarType>& res) const;

  /*!
   * \brief Factorize matrix using PaStiX.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] kind_fact - Type of factorization.
   */
  void BuildPastixPreconditioner(CGeometry* geometry, const CConfig* config, unsigned short kind_fact);

  /*!
   * \brief Apply the PaStiX factorization to CSysVec.
   * \param[in] vec - CSysVector to be multiplied by the preconditioner.
   * \param[out] prod - Result of the product M*vec.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputePastixPreconditioner(const CSysVector<ScalarType>& vec, CSysVector<ScalarType>& prod, CGeometry* geometry,
                                   const CConfig* config) const;
};
