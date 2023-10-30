/*!
 * \file CSysVector.hpp
 * \brief Declararion and inlines of the vector class used in the
 * solution of large, distributed, sparse linear systems.
 * \author P. Gomes, F. Palacios, J. Hicken, T. Economon
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

#include "../parallelization/mpi_structure.hpp"
#include "../parallelization/omp_structure.hpp"
#include "../parallelization/vectorization.hpp"
#include "vector_expressions.hpp"

/*!
 * \brief OpenMP worksharing construct used in CSysVector for loops.
 * \note The loop will only run in parallel if methods are called from a
 * parallel region (if not the results will still be correct).
 * Static schedule to reduce overhead, chunk size determined at initialization.
 * "nowait" clause is safe when calling CSysVector methods after each other
 * as the loop size is the same. Methods of other classes that operate on a
 * CSysVector and do not have the same work scheduling must use a
 * SU2_OMP_BARRIER before using the vector.
 */
#ifdef HAVE_OMP
#ifdef HAVE_OMP_SIMD
#define CSYSVEC_PARFOR SU2_OMP_FOR_(simd schedule(static, omp_chunk_size) SU2_NOWAIT)
#else
#define CSYSVEC_PARFOR SU2_OMP_FOR_(schedule(static, omp_chunk_size) SU2_NOWAIT)
#endif
#define END_CSYSVEC_PARFOR END_SU2_OMP_FOR
#else
#define CSYSVEC_PARFOR SU2_OMP_SIMD
#define END_CSYSVEC_PARFOR
#endif

/*!
 * \class CSysVector
 * \ingroup SpLinSys
 * \brief Class for holding and manipulating vectors needed by linear solvers.
 */
template <class ScalarType>
class CSysVector : public VecExpr::CVecExpr<CSysVector<ScalarType>, ScalarType> {
 private:
  enum { OMP_MAX_SIZE = 4096 }; /*!< \brief Maximum chunk size used in parallel for loops. */

  unsigned long omp_chunk_size = OMP_MAX_SIZE; /*!< \brief Static chunk size used in loops. */
  ScalarType* vec_val = nullptr;               /*!< \brief Storage, 64 byte aligned (do not use normal new/delete). */
  unsigned long nElm = 0;       /*!< \brief Total number of elements (or number elements on this processor). */
  unsigned long nElmDomain = 0; /*!< \brief Total number of elements without Ghost cells. */
  unsigned long nVar = 1;       /*!< \brief Number of elements in a block. */

  /*!
   * \brief Generic initialization from a scalar or array.
   * \note If val==nullptr vec_val is not initialized, only allocated.
   * \param[in] numBlk - Number of blocks locally.
   * \param[in] numBlkDomain - Number of blocks locally (without ghost cells).
   * \param[in] numVar - Number of variables in each block.
   * \param[in] val - Default value for elements.
   * \param[in] valIsArray - If true val is treated as array.
   * \param[in] errorIfParallel - Throw error if within parallel region (all ctors except the default one do this).
   */
  void Initialize(unsigned long numBlk, unsigned long numBlkDomain, unsigned long numVar, const ScalarType* val,
                  bool valIsArray, bool errorIfParallel = true);

  /*!
   * \brief Helper to unpack (transpose) a SIMD input block.
   */
  template <size_t N, size_t nVar, class VecTypeSIMD, class F>
  FORCEINLINE static void UnpackBlock(const VecTypeSIMD& in, simd::Array<F, N> mask, ScalarType out[][nVar]) {
    static_assert(VecTypeSIMD::StaticSize, "This method requires static size vectors.");
    for (size_t i = 0; i < nVar; ++i) {
      SU2_OMP_SIMD_IF_NOT_AD
      for (size_t k = 0; k < N; ++k) out[k][i] = mask[k] * in[i][k];
    }
  }

 public:
  static constexpr bool StoreAsRef = true; /*! \brief Required by CVecExpr. */

  /*!
   * \brief Default constructor of the class.
   */
  CSysVector() = default;

  /*!
   * \brief Destructor
   */
  ~CSysVector();

  /*!
   * \brief Construct from size and value.
   * \param[in] size - Number of elements locally.
   * \param[in] val - Default value for elements.
   */
  explicit CSysVector(unsigned long size, ScalarType val = 0.0) { Initialize(size, size, 1, &val, false); }

  /*!
   * \brief Construct from size and value (block version).
   * \param[in] numBlk - Number of blocks locally.
   * \param[in] numBlkDomain - Number of blocks locally (without ghost cells).
   * \param[in] numVar - Number of variables in each block.
   * \param[in] val - Default value for elements.
   */
  CSysVector(unsigned long numBlk, unsigned long numBlkDomain, unsigned long numVar, ScalarType val = 0.0) {
    Initialize(numBlk, numBlkDomain, numVar, &val, false);
  }

  /*!
   * \brief Construct from array.
   * \param[in] size - Number of elements locally.
   * \param[in] u_array - Vector stored as array being copied.
   */
  CSysVector(unsigned long size, const ScalarType* u_array) { Initialize(size, size, 1, u_array, true); }

  /*!
   * \brief Constructor from array (block version).
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain - number of blocks locally (without g cells)
   * \param[in] numVar - number of variables in each block
   * \param[in] u_array - vector stored as array being copied
   */
  CSysVector(unsigned long numBlk, unsigned long numBlkDomain, unsigned long numVar, const ScalarType* u_array) {
    Initialize(numBlk, numBlkDomain, numVar, u_array, true);
  }

  /*!
   * \brief Copy constructor of the class.
   * \note Not defined for expressions because we do not know their sizes.
   * \param[in] u - Vector being copied.
   */
  CSysVector(const CSysVector& u) { Initialize(u.GetNBlk(), u.GetNBlkDomain(), u.nVar, u.vec_val, true); }

  /*!
   * \brief Initialize the class with a scalar.
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain - number of blocks locally (without g cells)
   * \param[in] numVar - number of variables in each block
   * \param[in] val - default value for elements
   */
  void Initialize(unsigned long numBlk, unsigned long numBlkDomain, unsigned long numVar, ScalarType val = 0.0) {
    Initialize(numBlk, numBlkDomain, numVar, &val, false, false);
  }

  /*!
   * \brief Initialize the class with an array.
   * \note If ptr==nullptr no copy occurs.
   * \param[in] numBlk - number of blocks locally
   * \param[in] numBlkDomain - number of blocks locally (without g cells)
   * \param[in] numVar - number of variables in each block
   * \param[in] ptr - pointer to data with which to initialize the vector
   */
  void Initialize(unsigned long numBlk, unsigned long numBlkDomain, unsigned long numVar, const ScalarType* ptr) {
    Initialize(numBlk, numBlkDomain, numVar, ptr, true, false);
  }

  /*!
   * \brief Set our values (resizing if required) by copying from other, the derivative information is lost.
   * \param[in] other - source CSysVector
   */
  template <class T>
  void PassiveCopy(const CSysVector<T>& other) {
    /*--- This is a method and not the overload of an operator to make sure who
     * calls it knows the consequence to the derivative information (lost) ---*/

    /*--- check if self-assignment, otherwise perform deep copy ---*/
    if ((const void*)this == (const void*)&other) return;

    SU2_OMP_SAFE_GLOBAL_ACCESS(
        Initialize(other.GetNBlk(), other.GetNBlkDomain(), other.GetNVar(), nullptr, true, false);)

    CSYSVEC_PARFOR
    for (auto i = 0ul; i < nElm; i++) vec_val[i] = SU2_TYPE::GetValue(other[i]);
    END_CSYSVEC_PARFOR
  }

  /*!
   * \brief return the number of local elements in the CSysVector
   */
  inline unsigned long GetLocSize() const { return nElm; }

  /*!
   * \brief return the number of local elements in the CSysVector without ghost cells
   */
  inline unsigned long GetNElmDomain() const { return nElmDomain; }

  /*!
   * \brief return the number of variables at each block (typically number per node)
   */
  inline unsigned long GetNVar() const { return nVar; }

  /*!
   * \brief return the number of blocks (typically number of nodes locally)
   */
  inline unsigned long GetNBlk() const { return nElm / nVar; }

  /*!
   * \brief return the number of blocks (typically number of nodes locally)
   */
  inline unsigned long GetNBlkDomain() const { return nElmDomain / nVar; }

  /*!
   * \brief Access operator with assignment permitted.
   * \param[in] i - Local index to access.
   * \return Value at position i.
   */
  inline ScalarType& operator[](unsigned long i) { return vec_val[i]; }
  inline const ScalarType& operator[](unsigned long i) const { return vec_val[i]; }

  /*!
   * \brief Iterators for range for loops.
   */
  inline const ScalarType* begin() const { return vec_val; }
  inline const ScalarType* end() const { return vec_val + nElm; }

  /*!
   * \brief Access operator with assignment permitted block version.
   * \param[in] iPoint - Index of block.
   * \param[in] iVar - Index of variable.
   * \return Value at position (i,j).
   */
  inline ScalarType& operator()(unsigned long iPoint, unsigned long iVar) { return vec_val[iPoint * nVar + iVar]; }
  inline const ScalarType& operator()(unsigned long iPoint, unsigned long iVar) const {
    return vec_val[iPoint * nVar + iVar];
  }

  /*!
   * \brief Assignment operator from another vector.
   * \note Does not resize as it is meant for use in parallel.
   * \param[in] other - Another vector.
   */
  CSysVector& operator=(const CSysVector& other) {
    CSYSVEC_PARFOR
    for (auto i = 0ul; i < nElm; ++i) vec_val[i] = other.vec_val[i];
    END_CSYSVEC_PARFOR
    return *this;
  }

  /*!
   * \brief Compound assignement operations with scalars and expressions.
   * \param[in] val/expr - Scalar value or expression.
   */
#define MAKE_COMPOUND(OP)                                                 \
  CSysVector& operator OP(ScalarType val) {                               \
    CSYSVEC_PARFOR                                                        \
    for (auto i = 0ul; i < nElm; ++i) vec_val[i] OP val;                  \
    END_CSYSVEC_PARFOR                                                    \
    return *this;                                                         \
  }                                                                       \
  template <class T>                                                      \
  CSysVector& operator OP(const VecExpr::CVecExpr<T, ScalarType>& expr) { \
    CSYSVEC_PARFOR                                                        \
    for (auto i = 0ul; i < nElm; ++i) vec_val[i] OP expr.derived()[i];    \
    END_CSYSVEC_PARFOR                                                    \
    return *this;                                                         \
  }
  MAKE_COMPOUND(=)
  MAKE_COMPOUND(+=)
  MAKE_COMPOUND(-=)
  MAKE_COMPOUND(*=)
  MAKE_COMPOUND(/=)
#undef MAKE_COMPOUND

  /*!
   * \brief Sets to zero all the entries of the vector.
   */
  inline void SetValZero(void) { *this = ScalarType(0); }

  /*!
   * \brief Dot product between "this" and an expression.
   * \param[in] expr - Expression.
   * \return Result of dot product
   */
  template <class T>
  ScalarType dot(const VecExpr::CVecExpr<T, ScalarType>& expr) const {
    static ScalarType dotRes;
    /*--- All threads get the same "view" of the vectors and shared variable. ---*/
    SU2_OMP_SAFE_GLOBAL_ACCESS(dotRes = 0.0;)

    /*--- Local dot product for each thread. ---*/
    ScalarType sum = 0.0;

    CSYSVEC_PARFOR
    for (auto i = 0ul; i < nElmDomain; ++i) {
      sum += vec_val[i] * expr.derived()[i];
    }
    END_CSYSVEC_PARFOR

    /*--- Update shared variable with "our" partial sum. ---*/
    atomicAdd(sum, dotRes);

#ifdef HAVE_MPI
    /*--- Reduce across all mpi ranks, only master thread communicates. ---*/
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      sum = dotRes;
      const auto mpi_type = (sizeof(ScalarType) < sizeof(double)) ? MPI_FLOAT : MPI_DOUBLE;
      SelectMPIWrapper<ScalarType>::W::Allreduce(&sum, &dotRes, 1, mpi_type, MPI_SUM, SU2_MPI::GetComm());
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
#else
    /*--- Make view of result consistent across threads. ---*/
    SU2_OMP_BARRIER
#endif

    return dotRes;
  }

  /*!
   * \brief Squared L2 norm of the vector (via dot with self).
   * \return Squared L2 norm.
   */
  inline ScalarType squaredNorm() const { return dot(*this); }

  /*!
   * \brief L2 norm of the vector.
   * \return L2 norm.
   */
  inline ScalarType norm() const { return sqrt(squaredNorm()); }

  /*!
   * \brief Get pointer to a block.
   * \param[in] iPoint - Index of block.
   * \return Pointer to start of block.
   */
  inline ScalarType* GetBlock(unsigned long iPoint) { return &vec_val[iPoint * nVar]; }
  inline const ScalarType* GetBlock(unsigned long iPoint) const { return &vec_val[iPoint * nVar]; }

  /*!
   * \brief Set the values to zero for one block.
   * \param[in] iPoint - Index of the block being set to zero.
   */
  inline void SetBlock_Zero(unsigned long iPoint) {
    for (auto iVar = 0ul; iVar < nVar; iVar++) vec_val[iPoint * nVar + iVar] = 0.0;
  }

  /*!
   * \brief Set "block" to the vector.
   * \note Template param Overwrite can be set to false to update existing values.
   * \param[in] iPoint - index of the point where set the residual.
   * \param[in] block - Value to set to the residual.
   * \param[in] alpha - Scale factor (axpy-type operation).
   */
  template <class VectorType, bool Overwrite = true>
  FORCEINLINE void SetBlock(unsigned long iPoint, const VectorType& block, ScalarType alpha = 1) {
    if (Overwrite) {
      for (auto i = 0ul; i < nVar; ++i) vec_val[iPoint * nVar + i] = alpha * block[i];
    } else {
      for (auto i = 0ul; i < nVar; ++i) vec_val[iPoint * nVar + i] += alpha * block[i];
    }
  }

  /*!
   * \brief Add "block" to the vector, see SetBlock.
   */
  template <class VectorType>
  FORCEINLINE void AddBlock(unsigned long iPoint, const VectorType& block, ScalarType alpha = 1) {
    SetBlock<VectorType, false>(iPoint, block, alpha);
  }

  /*!
   * \brief Subtract "block" from the vector, see AddBlock.
   */
  template <class VectorType>
  FORCEINLINE void SubtractBlock(unsigned long iPoint, const VectorType& block) {
    AddBlock(iPoint, block, -1);
  }

  /*!
   * \brief Add to iPoint, subtract from jPoint.
   */
  template <class VectorType>
  FORCEINLINE void UpdateBlocks(unsigned long iPoint, unsigned long jPoint, const VectorType& block,
                                ScalarType alpha = 1) {
    AddBlock(iPoint, block, alpha);
    AddBlock(jPoint, block, -alpha);
  }

  /*!
   * \brief Vectorized version of SetBlock, sets multiple iPoint's.
   * \param[in] iPoint - SIMD integer, the positions to update.
   * \param[in] vector - Vector of SIMD scalars.
   * \param[in] mask - Optional scale factor (axpy type operation).
   * \note Nothing is updated if the mask is 0.
   */
  template <size_t N, class T, class VecTypeSIMD, class F = ScalarType>
  FORCEINLINE void SetBlock(simd::Array<T, N> iPoint, const VecTypeSIMD& vector, simd::Array<F, N> mask = 1) {
    /*--- "Transpose" and scale input vector. ---*/
    constexpr size_t nVar = VecTypeSIMD::StaticSize;
    assert(nVar == this->nVar);
    ScalarType vec[N][nVar];
    UnpackBlock(vector, mask, vec);

    /*--- Update one by one skipping if mask is 0. ---*/
    for (size_t k = 0; k < N; ++k) {
      if (mask[k] == 0) continue;
      SU2_OMP_SIMD
      for (size_t i = 0; i < nVar; ++i) vec_val[iPoint[k] * nVar + i] = vec[k][i];
    }
  }

  /*!
   * \brief Vectorized version of UpdateBlocks, updates multiple i/jPoint's.
   * \note See SIMD overload of SetBlock.
   */
  template <size_t N, class T, class VecTypeSIMD, class F = ScalarType>
  FORCEINLINE void UpdateBlocks(simd::Array<T, N> iPoint, simd::Array<T, N> jPoint, const VecTypeSIMD& vector,
                                simd::Array<F, N> mask = 1) {
    /*--- "Transpose" and scale input vector. ---*/
    constexpr size_t nVar = VecTypeSIMD::StaticSize;
    assert(nVar == this->nVar);
    ScalarType vec[N][nVar];
    UnpackBlock(vector, mask, vec);

    /*--- Update one by one skipping if mask is 0. ---*/
    for (size_t k = 0; k < N; ++k) {
      if (mask[k] == 0) continue;
      SU2_OMP_SIMD
      for (size_t i = 0; i < nVar; ++i) {
        vec_val[iPoint[k] * nVar + i] += vec[k][i];
        vec_val[jPoint[k] * nVar + i] -= vec[k][i];
      }
    }
  }
};

#undef CSYSVEC_PARFOR
#undef END_CSYSVEC_PARFOR
