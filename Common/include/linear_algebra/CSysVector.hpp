/*!
 * \file CSysVector.hpp
 * \brief Declararion and inlines of the vector class used in the
 * solution of large, distributed, sparse linear systems.
 * \author P. Gomes, F. Palacios, J. Hicken, T. Economon
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

#include <memory>
#include <vector>

#include "../code_config.hpp"
#include "../parallelization/mpi_structure.hpp"
#include "../parallelization/omp_structure.hpp"
#include "../parallelization/vectorization.hpp"
#include "vector_expressions.hpp"
#include "../../include/CConfig.hpp"

#ifdef __CUDACC__
#include "GPUComms.cuh"
#endif

template <class ScalarType>
class CSysVector;

/*!
 * \brief True for the plain floating-point scalar types the GPU vector kernels
 *        (CSysVectorGPU.cu) are instantiated for, in builds where those kernels
 *        exist at all. AD active types are never dispatched to the device: the
 *        tape/expression machinery those types pull in is not compatible with
 *        nvcc's device-code compilation, and device-resident autodiff is not
 *        supported by this GPU path.
 * \note In an AD build this is false even for su2mixedfloat, because the .cu
 *       translation units are not linked into the AD libraries (see
 *       SU2_ENABLE_CUDA_KERNELS in code_config.hpp).
 */
#ifdef SU2_ENABLE_CUDA_KERNELS
template <class ScalarType>
inline constexpr bool su2_gpu_capable_v = std::is_floating_point_v<ScalarType>;
#else
template <class ScalarType>
inline constexpr bool su2_gpu_capable_v = false;
#endif

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
 * \brief Brackets device work so that it is issued by a single thread with the whole team
 * synchronized before and after.
 * \note The GPU is one shared resource and the device path does not use OpenMP worksharing.
 * Issuing from one thread keeps kernel launches ordered on the default stream and, above
 * all, stops part of a team from entering a worksharing construct that the rest skipped.
 * Correctness relies on all threads reaching the same vector operations in the same order,
 * which is the assumption the "nowait" clause on CSYSVEC_PARFOR already makes. These
 * regions must not be nested; they are used by the operations of this class and by the
 * matrix-vector product and preconditioner wrappers, and by nothing above those.
 */
#define SU2_DEVICE_REGION(...) SU2_OMP_SAFE_GLOBAL_ACCESS(__VA_ARGS__)
#define BEGIN_SU2_DEVICE_REGION BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
#define END_SU2_DEVICE_REGION END_SU2_OMP_SAFE_GLOBAL_ACCESS

namespace VecExpr {

enum class DeviceAssignOp { Assign, Add, Subtract, Multiply, Divide };

/*!
 * \brief Whether vector expressions are currently evaluated on the device.
 * \note Defined in CSysVectorGPU.cu. This is a plain global, not per thread: every thread
 * of a team has to agree on it or they would split over the worksharing constructs below.
 * It is switched by CSysSolve at the same boundary that uploads and downloads the vectors
 * (HandleTemporariesIn/Out), and nowhere else, so it does not change while a solve runs.
 */
#ifdef SU2_ENABLE_CUDA_KERNELS
bool UseDeviceExpressions();
void SetUseDeviceExpressions(bool use);
#else
inline bool UseDeviceExpressions() { return false; }
inline void SetUseDeviceExpressions(bool) {}
#endif

/*!
 * \brief How a CSysVector is captured inside an expression: a bare pointer to whichever
 * storage the expression is going to be evaluated from.
 * \note Capturing by value (rather than a reference to the vector) is what makes an
 * arbitrary expression tree trivially copyable, and therefore passable by value to the
 * assignment kernel. The choice of storage is fixed when the expression is built, which is
 * sound because there is no fallback: while UseDeviceExpressions() holds, every expression
 * is evaluated by a kernel.
 */
template <class Scalar>
class CVectorView : public CVecExpr<CVectorView<Scalar>, Scalar> {
 private:
  const Scalar* data = nullptr;

 public:
  static constexpr bool StoreAsRef = false;

  CVectorView(const CSysVector<Scalar>& vector);

  SU2_CUDA_HOST_DEVICE FORCEINLINE const Scalar& operator[](size_t i) const { return data[i]; }
};

template <class Scalar>
struct store_type<CSysVector<Scalar>> {
  using type = CVectorView<Scalar>;
};

template <class Scalar>
struct store_type<const CSysVector<Scalar>> {
  using type = CVectorView<Scalar>;
};

template <DeviceAssignOp Op, class Scalar, class T>
void AssignDeviceExpression(Scalar* data, unsigned long size, const CVecExpr<T, Scalar>& expr);

#ifdef __CUDACC__
template <DeviceAssignOp Op, class Scalar, class T>
__global__ void DeviceAssignKernel(Scalar* data, unsigned long size, T expr) {
  const unsigned long i = static_cast<unsigned long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i >= size) return;

  if constexpr (Op == DeviceAssignOp::Assign) {
    data[i] = expr[i];
  } else if constexpr (Op == DeviceAssignOp::Add) {
    data[i] += expr[i];
  } else if constexpr (Op == DeviceAssignOp::Subtract) {
    data[i] -= expr[i];
  } else if constexpr (Op == DeviceAssignOp::Multiply) {
    data[i] *= expr[i];
  } else {
    data[i] /= expr[i];
  }
}

template <DeviceAssignOp Op, class Scalar, class T>
inline void AssignDeviceExpression(Scalar* data, unsigned long size, const CVecExpr<T, Scalar>& expr) {
  if (size == 0) return;
  constexpr unsigned block_size = 256;
  const auto grid_size = static_cast<unsigned>((size + block_size - 1) / block_size);
  DeviceAssignKernel<Op><<<grid_size, block_size>>>(data, size, expr.derived());
  gpuErrChk(cudaPeekAtLastError());
}
#endif

}  // namespace VecExpr

/*!
 * \class CSysVector
 * \ingroup SpLinSys
 * \brief Class for holding and manipulating vectors needed by linear solvers.
 */
template <class ScalarType>
class CSysVector : public VecExpr::CVecExpr<CSysVector<ScalarType>, ScalarType> {
 private:
  enum { OMP_MAX_SIZE = 4096 }; /*!< \brief Maximum chunk size used in parallel for loops. */

  /// NOTE: Update swap() if you add member variables.
  unsigned long omp_chunk_size = OMP_MAX_SIZE; /*!< \brief Static chunk size used in loops. */
  ScalarType* vec_val = nullptr;               /*!< \brief Storage, 64 byte aligned (do not use normal new/delete). */
  unsigned long nElm = 0;       /*!< \brief Total number of elements (or number elements on this processor). */
  unsigned long nElmDomain = 0; /*!< \brief Total number of elements without Ghost cells. */
  unsigned long nVar = 1;       /*!< \brief Number of elements in a block. */

  ScalarType* d_vec_val = nullptr; /*!< \brief Device Pointer to store the vector values on the GPU. */

#ifdef HAVE_OMP
  mutable std::unique_ptr<ScalarType[]>
      dot_scratch; /*!< \brief Stores partial sums for ordered reduction over OMP threads. */
#else
  mutable std::array<ScalarType, 1> dot_scratch;
#endif

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

  /*!
   * \brief Evaluates an expression into the device storage of this vector.
   * \note The kernel has to be instantiated for the expression type in CSysVectorGPU.cu,
   * a shape that is not in that list is an undefined symbol at link time.
   */
  template <VecExpr::DeviceAssignOp Op, class T>
  CSysVector& AssignDevice(const VecExpr::CVecExpr<T, ScalarType>& expr) {
#ifdef SU2_ENABLE_CUDA_KERNELS
    if constexpr (su2_gpu_capable_v<ScalarType>) {
      BEGIN_SU2_DEVICE_REGION {
        VecExpr::store_t<const T> stored_expr(expr.derived());
        VecExpr::AssignDeviceExpression<Op>(d_vec_val, nElm, stored_expr);
      }
      END_SU2_DEVICE_REGION
    }
#endif
    return *this;
  }

  template <VecExpr::DeviceAssignOp Op>
  CSysVector& AssignDevice(ScalarType val) {
#ifdef SU2_ENABLE_CUDA_KERNELS
    if constexpr (su2_gpu_capable_v<ScalarType>) {
      BEGIN_SU2_DEVICE_REGION
      VecExpr::AssignDeviceExpression<Op>(d_vec_val, nElm, VecExpr::Bcast<ScalarType>(val));
      END_SU2_DEVICE_REGION
    }
#endif
    return *this;
  }

  /*!
   * \brief GPU helper for `dot`.
   */
  ScalarType dotGPU(const CSysVector& other) const;

  /*!
   * \brief GPU helper for multiDot.
   */
  static su2matrix<ScalarType> multiDotGPU(const std::vector<CSysVector<ScalarType>>& V, size_t i0, size_t n,
                                           const std::vector<CSysVector<ScalarType>>& W, size_t m);

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
   * \brief Swap contents with another vector.
   */
  void swap(CSysVector& other) {
    std::swap(omp_chunk_size, other.omp_chunk_size);
    std::swap(vec_val, other.vec_val);
    std::swap(d_vec_val, other.d_vec_val);
    std::swap(nElm, other.nElm);
    std::swap(nElmDomain, other.nElmDomain);
    std::swap(nVar, other.nVar);
    std::swap(dot_scratch, other.dot_scratch);
  }

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
   * \brief Performs the memory copy from host to device.
   * \param[in] trigger - boolean value that decides whether to conduct the transfer or not. True by default.
   */
  void HtDTransfer(bool trigger = true) const;

  /*!
   * \brief Performs the memory copy from device to host.
   * \param[in] trigger - boolean value that decides whether to conduct the transfer or not. True by default.
   */
  void DtHTransfer(bool trigger = true) const;

  /*!
   * \brief return device pointer that points to the CSysVector values in GPU memory
   */
  inline ScalarType* GetDevicePointer() const { return d_vec_val; }

  /*!
   * \brief return host pointer that points to the CSysVector values, counterpart of
   * GetDevicePointer
   */
  inline const ScalarType* GetHostPointer() const { return vec_val; }

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
#ifdef SU2_ENABLE_CUDA_KERNELS
    if constexpr (su2_gpu_capable_v<ScalarType>) {
      if (VecExpr::UseDeviceExpressions()) return AssignDevice<VecExpr::DeviceAssignOp::Assign>(other);
    }
#endif
    CSYSVEC_PARFOR
    for (auto i = 0ul; i < nElm; ++i) vec_val[i] = other.vec_val[i];
    END_CSYSVEC_PARFOR
    return *this;
  }

  /*!
   * \brief Compound assignement operations with scalars and expressions.
   * \param[in] val/expr - Scalar value or expression.
   */
#define MAKE_COMPOUND(OP, ASSIGN_OP)                                             \
  CSysVector& operator OP(ScalarType val) {                                      \
    if constexpr (su2_gpu_capable_v<ScalarType>) {                               \
      if (VecExpr::UseDeviceExpressions()) return AssignDevice<ASSIGN_OP>(val);  \
    }                                                                            \
    CSYSVEC_PARFOR                                                               \
    for (auto i = 0ul; i < nElm; ++i) vec_val[i] OP val;                         \
    END_CSYSVEC_PARFOR                                                           \
    return *this;                                                                \
  }                                                                              \
  template <class T>                                                             \
  CSysVector& operator OP(const VecExpr::CVecExpr<T, ScalarType>& expr) {        \
    if constexpr (su2_gpu_capable_v<ScalarType>) {                               \
      if (VecExpr::UseDeviceExpressions()) return AssignDevice<ASSIGN_OP>(expr); \
    }                                                                            \
    CSYSVEC_PARFOR                                                               \
    for (auto i = 0ul; i < nElm; ++i) vec_val[i] OP expr.derived()[i];           \
    END_CSYSVEC_PARFOR                                                           \
    return *this;                                                                \
  }
  MAKE_COMPOUND(=, VecExpr::DeviceAssignOp::Assign)
  MAKE_COMPOUND(+=, VecExpr::DeviceAssignOp::Add)
  MAKE_COMPOUND(-=, VecExpr::DeviceAssignOp::Subtract)
  MAKE_COMPOUND(*=, VecExpr::DeviceAssignOp::Multiply)
  MAKE_COMPOUND(/=, VecExpr::DeviceAssignOp::Divide)
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
#ifdef SU2_ENABLE_CUDA_KERNELS
    if constexpr (su2_gpu_capable_v<ScalarType>) {
      using DeviceExpr = std::remove_cv_t<VecExpr::remove_reference_t<T>>;
      static_assert(std::is_same_v<DeviceExpr, CSysVector>,
                    "On the device the dot product needs a real device pointer (dotGPU takes a "
                    "materialized vector, not an expression template), so it only takes vectors. "
                    "Assign the expression to a vector first.");
      if (VecExpr::UseDeviceExpressions()) {
        /*--- dotGPU reduces over MPI, which has to happen once for the team, so the result
         * is published through the same scratch slot the host reduction below uses. ---*/
        SU2_DEVICE_REGION(dot_scratch[0] = dotGPU(expr.derived());)
        return dot_scratch[0];
      }
    }
#endif

    /*--- All threads get the same "view" of the vectors. ---*/
    SU2_OMP_BARRIER

    /*--- Local dot product for each thread. ---*/
    ScalarType sum = 0.0;

    CSYSVEC_PARFOR
    for (auto i = 0ul; i < nElmDomain; ++i) {
      sum += vec_val[i] * expr.derived()[i];
    }
    END_CSYSVEC_PARFOR

    dot_scratch[omp_get_thread_num()] = sum;

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      /*--- Reduce over all threads in an ordered way to ensure a deterministic result. ---*/
      for (int i = 1; i < omp_get_num_threads(); ++i) sum += dot_scratch[i];

      /*--- Reduce across all mpi ranks, only the master thread communicates. ---*/
      const auto mpi_type = (sizeof(ScalarType) < sizeof(double)) ? MPI_FLOAT : MPI_DOUBLE;
      SelectMPIWrapper<ScalarType>::W::Allreduce(&sum, &dot_scratch[0], 1, mpi_type, MPI_SUM, SU2_MPI::GetComm());
    }
    /*--- Make view of result consistent across threads. ---*/
    END_SU2_OMP_SAFE_GLOBAL_ACCESS

    return dot_scratch[0];
  }

  /*!
   * \brief Computes the product of V^T W efficiencly, where V and W are tall matrices stored as vectors of CSysVector.
   * \param[in] V - Tall matrix.
   * \param[in] i0 - First column of V to consider.
   * \param[in] n - Number of columns to consider from V starting at i0.
   * \param[in] W - Tall matrix.
   * \param[in] m - Number of columns to consider from W.
   * \return n by m matrix with the result of the product.
   */
  static const su2matrix<ScalarType>& multiDot(const std::vector<CSysVector>& V, size_t i0, size_t n,
                                               const std::vector<CSysVector>& W, size_t m);

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

namespace VecExpr {

template <class Scalar>
CVectorView<Scalar>::CVectorView(const CSysVector<Scalar>& vector)
    : data(UseDeviceExpressions() ? vector.GetDevicePointer() : vector.GetHostPointer()) {}

}  // namespace VecExpr

#undef CSYSVEC_PARFOR
#undef END_CSYSVEC_PARFOR
