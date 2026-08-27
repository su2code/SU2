/*!
 * \file util.hpp
 * \brief Generic auxiliary functions.
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

#pragma once

#include <type_traits>

#include "../../../Common/include/option_structure.hpp"
#include "../../../Common/include/parallelization/vectorization.hpp"
#include "../../../Common/include/containers/C2DContainer.hpp"
#include "../../../Common/include/linear_algebra/CSysVector.hpp"
#include "../../../Common/include/linear_algebra/CSysMatrix.hpp"

/*!
 * \enum UpdateType
 * \brief Ways to update vectors and system matrices.
 * COLORING is the typical i/j update, whereas for REDUCTION
 * the fluxes are stored and the matrix diagonal is not modified.
 */
enum class UpdateType { COLORING, REDUCTION };

#ifdef CODI_FORWARD_TYPE
using SparseMatrixType = CSysMatrix<su2double>;
#else
using SparseMatrixType = CSysMatrix<su2mixedfloat>;
#endif

/*!
 * \brief Alignment of the static containers backing a flux value type.
 * \note Yields the type's own alignment for a SIMD array, and the container default (0)
 *       for a plain scalar, which has no alignment of its own.
 */
template <class Type>
struct CAlignTraits {
  enum : size_t { Align = 0 };
};

template <class Scalar_t, size_t N>
struct CAlignTraits<simd::Array<Scalar_t, N>> {
  enum : size_t { Align = simd::Array<Scalar_t, N>::Align };
};

/*!
 * \brief Static vector and matrix types.
 * \note These should be used instead of C-style arrays.
 */
template <class Type, size_t Size>
using Vector = C2DContainer<unsigned long, Type, StorageType::ColumnMajor, CAlignTraits<Type>::Align, Size, 1>;

template <class Type, size_t Rows, size_t Cols>
using Matrix = C2DContainer<unsigned long, Type, StorageType::RowMajor, CAlignTraits<Type>::Align, Rows, Cols>;

/*!
 * \brief The flux value type and lane count that go with an index type.
 * \note There is exactly one floating type in play, su2double, active under AD; a plain
 *       integral index reads one of it, a lane-vector index reads a lane-vector of it. This
 *       is what lets every helper below deduce its value type from the index it is handed,
 *       instead of a caller naming it explicitly.
 */
template <class Int>
struct CValueTraits;

template <>
struct CValueTraits<unsigned long> {
  using Double = su2double;
  static constexpr size_t Size = 1;
};

template <size_t N>
struct CValueTraits<simd::Array<unsigned long, N>> {
  using Double = simd::Array<su2double, N>;
  static constexpr size_t Size = N;
};

/*!
 * \brief Index type and lane count that go with a flux value type, the converse of
 *        CValueTraits, used where the value type is already known (e.g. a class template
 *        parameter) and the index type is what needs deriving.
 */
template <class Double>
struct CLaneTraits;

template <>
struct CLaneTraits<su2double> {
  using Int = unsigned long;
  static constexpr size_t Size = 1;
};

template <class Scalar_t, size_t N>
struct CLaneTraits<simd::Array<Scalar_t, N>> {
  using Int = simd::Array<unsigned long, N>;
  static constexpr size_t Size = N;
};

/*!
 * \brief Constexpr version of max.
 */
inline constexpr size_t Max(size_t a, size_t b) { return a > b ? a : b; }

/*!
 * \brief Simple pair type for i/j variables.
 */
template <class T>
struct CPair {
  T i, j;
};

/*!
 * \brief Blocks a template parameter from participating in argument deduction.
 * \note Deduction never applies a user conversion, so a parameter typed plain Double would
 *       force a caller passing a bare su2double constant (kappa, a limiter ramp) to have
 *       already broadcast it. Wrapping the parameter type here defers Double entirely to
 *       the other, genuinely deduced arguments, and the broadcast then happens as an
 *       ordinary implicit conversion at the call.
 */
template <class T>
struct CIdentity {
  using type = T;
};
template <class T>
using CNonDeduced = typename CIdentity<T>::type;

/*!
 * \brief Equation count of a model whose value is only known at runtime.
 */
constexpr size_t Dynamic = size_t(-1);

/*!
 * \brief Backing size of the static arrays of a dynamic model.
 * \note The scalar numerics cap the equation count at this value and error above it, so a
 *       configuration that fits them fits these kernels.
 */
constexpr size_t MaxScalarVar = 8;

/*!
 * \brief Residual of one edge, accumulated by the convective and the diffusive terms.
 * \note flux_i and flux_j are the contributions to the rows of i and j. They are opposite
 *       for a conservative term and independent for a non-conservative one. The Jacobians
 *       map onto the ii, ij, ji and jj blocks of the edge. A dynamic model sizes the storage
 *       with the maximum and iterates to nVar, which the scheme sets from the solver.
 */
template <class Double, size_t nVar_>
struct EdgeResidual {
  /*!< \brief The Matrix<Double,1,1> a static nVar==1 model would otherwise need degenerates
   *   to vector-only indexing in C2DContainer (its RowMajor, one-row specialization), so the
   *   backing is never smaller than 2; the unused padding row/column is simply never visited,
   *   every loop here and in the model bounds itself to nVar, not Size. */
  static constexpr size_t Size = Max(2, (nVar_ == Dynamic) ? MaxScalarVar : nVar_);

  Vector<Double, Size> flux_i, flux_j;
  Matrix<Double, Size, Size> jac_ii, jac_ij, jac_ji, jac_jj;
  const size_t nVar;

  /*!
   * \brief Zero the terms of the equations in use, so that both terms can accumulate into them.
   * \note A static model zeroes its whole storage with constant trip counts; a dynamic one
   *       zeroes the leading nVar rows and columns and leaves the rest of the backing untouched.
   */
  FORCEINLINE explicit EdgeResidual(size_t nEqn = Size) : nVar(nEqn) {
    for (size_t iVar = 0; iVar < nVar; ++iVar) {
      flux_i(iVar) = 0.0;
      flux_j(iVar) = 0.0;
      for (size_t jVar = 0; jVar < nVar; ++jVar) {
        jac_ii(iVar, jVar) = 0.0;
        jac_ij(iVar, jVar) = 0.0;
        jac_ji(iVar, jVar) = 0.0;
        jac_jj(iVar, jVar) = 0.0;
      }
    }
  }
};

/*!
 * \brief Dot product.
 */
template <size_t nDim, class ForwardIterator, class T>
FORCEINLINE auto dot(ForwardIterator iterator, const T* ptr) -> typename std::decay<decltype(*ptr)>::type {
  typename std::decay<decltype(*ptr)>::type sum = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    sum += *(iterator++) * ptr[iDim];
  }
  return sum;
}

/*!
 * \overload Dot product.
 */
template <size_t nDim, class Double, class ForwardIterator>
FORCEINLINE Double dot(ForwardIterator iterator, const Vector<Double, nDim>& vector) {
  return dot<nDim>(iterator, vector.data());
}

/*!
 * \overload Dot product.
 */
template <size_t nDim, class Double>
FORCEINLINE Double dot(const Vector<Double, nDim>& a, const Vector<Double, nDim>& b) {
  return dot<nDim>(a.data(), b.data());
}

/*!
 * \brief Squared norm.
 */
template <size_t nDim, class ForwardIterator>
FORCEINLINE auto squaredNorm(ForwardIterator iterator) -> typename std::decay<decltype(*iterator)>::type {
  typename std::decay<decltype(*iterator)>::type sum = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    sum += pow(*(iterator++), 2);
  }
  return sum;
}

/*!
 * \overload Squared norm.
 */
template <size_t nDim, class Double>
FORCEINLINE Double squaredNorm(const Vector<Double, nDim>& vector) {
  return squaredNorm<nDim>(vector.data());
}

/*!
 * \brief Tangential projection.
 */
template <size_t nDim, class Double>
FORCEINLINE Vector<Double, nDim> tangentProjection(const Matrix<Double, nDim, nDim>& tensor,
                                                    const Vector<Double, nDim>& unitVector) {
  Vector<Double, nDim> proj;
  for (size_t iDim = 0; iDim < nDim; ++iDim) proj(iDim) = dot(tensor[iDim], unitVector);

  Double normalProj = dot(proj, unitVector);

  for (size_t iDim = 0; iDim < nDim; ++iDim) proj(iDim) -= normalProj * unitVector(iDim);

  return proj;
}

/*!
 * \brief Vector norm.
 */
template <size_t nDim, class Double>
FORCEINLINE Double norm(const Vector<Double, nDim>& vector) {
  return sqrt(squaredNorm(vector));
}

#ifndef CODI_REVERSE_TYPE
/*!
 * \brief Gather a single variable, from column iVar (0 by default) of row iPoint of a
 *        2D container, or from index iPoint of a 1D container.
 */
template <class Int, class Container, class Double = typename CValueTraits<Int>::Double>
FORCEINLINE Double gatherVariables(Int iPoint, const Container& vars, size_t iVar = 0) {
  return vars.template get<Vector<Double, 1>>(iPoint, iVar)(0);
}

/*!
 * \brief Gather nVar contiguous variables starting at column iVar (0 by default) of row
 *        iPoint of a 2D container.
 */
template <size_t nVar, class Int, class Container, class Double = typename CValueTraits<Int>::Double>
FORCEINLINE Vector<Double, nVar> gatherVariables(Int iPoint, const Container& vars, size_t iVar = 0) {
  return vars.template get<Vector<Double, nVar>>(iPoint, iVar);
}

/*!
 * \brief Gather an nRows x nCols block of a 3D container, from outer index iPoint and
 *        starting at middle index iRow.
 */
template <size_t nRows, size_t nCols, class Int, class Container, class Double = typename CValueTraits<Int>::Double>
FORCEINLINE Matrix<Double, nRows, nCols> gatherVariables(Int iPoint, const Container& vars, size_t iRow = 0) {
  return vars.template get<Matrix<Double, nRows, nCols>>(iPoint, iRow);
}
#else

namespace {
template <class Container, su2enable_if<Container::IsVector> = 0>
FORCEINLINE const su2double& get(const Container& vars, unsigned long iPoint) {
  return vars(iPoint);
}

/*--- When getting 1 variable from a matrix container, we assume it is the first. ---*/
template <class Container, su2enable_if<!Container::IsVector> = 0>
FORCEINLINE const su2double& get(const Container& vars, unsigned long iPoint, size_t iVar = 0) {
  return vars(iPoint, iVar);
}
}  // namespace

template <class Int, class Container, class Double = typename CValueTraits<Int>::Double>
FORCEINLINE Double gatherVariables(Int iPoint, const Container& vars, size_t iVar = 0) {
  Double x;
  for (size_t k = 0; k < CValueTraits<Int>::Size; ++k) {
    AD::SetPreaccIn(get(vars, iPoint[k], iVar));
    x[k] = get(vars, iPoint[k], iVar);
  }
  return x;
}

template <size_t nVar, class Int, class Container, class Double = typename CValueTraits<Int>::Double>
FORCEINLINE Vector<Double, nVar> gatherVariables(Int iPoint, const Container& vars, size_t iVar = 0) {
  Vector<Double, nVar> x;
  for (size_t i = 0; i < nVar; ++i) {
    for (size_t k = 0; k < CValueTraits<Int>::Size; ++k) {
      AD::SetPreaccIn(vars(iPoint[k], iVar + i));
      x[i][k] = vars(iPoint[k], iVar + i);
    }
  }
  return x;
}

template <size_t nRows, size_t nCols, class Int, class Container, class Double = typename CValueTraits<Int>::Double>
FORCEINLINE Matrix<Double, nRows, nCols> gatherVariables(Int iPoint, const Container& vars, size_t iRow = 0) {
  Matrix<Double, nRows, nCols> x;
  for (size_t i = 0; i < nRows; ++i) {
    for (size_t j = 0; j < nCols; ++j) {
      for (size_t k = 0; k < CValueTraits<Int>::Size; ++k) {
        AD::SetPreaccIn(vars(iPoint[k], iRow + i, j));
        x(i, j)[k] = vars(iPoint[k], iRow + i, j);
      }
    }
  }
  return x;
}
#endif

/*!
 * \brief Stop the AD preaccumulation.
 */
template <class Double, size_t nVar>
FORCEINLINE void stopPreacc(Vector<Double, nVar>& x) {
  AD::SetPreaccOut(x, nVar, CLaneTraits<Double>::Size);
  AD::EndPreacc();
}

/*!
 * \brief Distance vector, from point i to point j of one container.
 */
template <size_t nDim, class Int, class Container, class Double = typename CValueTraits<Int>::Double>
FORCEINLINE Vector<Double, nDim> distanceVector(Int iPoint, Int jPoint, const Container& coords) {
  return distanceVector<nDim>(iPoint, coords, jPoint, coords);
}

/*!
 * \brief Distance vector, from point i of one container to point j of another.
 * \note The two endpoints of a boundary flux read different containers, the solver's own
 *       and the marker's ghost one; the interior edge loop passes the same container twice.
 */
template <size_t nDim, class Int, class Container, class Double = typename CValueTraits<Int>::Double>
FORCEINLINE Vector<Double, nDim> distanceVector(Int iPoint, const Container& coords_i, Int jPoint,
                                                const Container& coords_j) {
  auto coord_i = gatherVariables<nDim>(iPoint, coords_i);
  auto coord_j = gatherVariables<nDim>(jPoint, coords_j);
  Vector<Double, nDim> vector_ij;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    vector_ij(iDim) = coord_j(iDim) - coord_i(iDim);
  }
  return vector_ij;
}

/*!
 * \brief Blended difference for U-MUSCL reconstruction.
 * \param[in] gradProj - Gradient projection at point i: dot(grad_i, vector_ij).
 * \param[in] delta - Centered difference: V_j - V_i.
 * \param[in] kappa - Blending parameter.
 * \return Blended difference for reconstruction from point i.
 */
template <class Double>
FORCEINLINE Double umusclProjection(const Double& gradProj, const Double& delta, const CNonDeduced<Double>& kappa) {
  /*-------------------------------------------------------------------*/
  /*--- The MUSCL kappa-scheme reconstruction is typically written: ---*/
  /*---     V_L = V_i + 0.25 * dV_ij^kap, where                     ---*/
  /*---     dV_ij^kap = (1-kappa) dV_ij^upw + (1+kappa) dV_ij^cen,  ---*/
  /*---     dV_ij^cen = V_j - V_i,                                  ---*/
  /*---     dV_ij^upw = 2 grad(Vi) dot vector_ij - dV_ij^cen.       ---*/
  /*--- To maintain proper scaling for edge limiters, the result of ---*/
  /*--- this function is 0.5 * dV_ij^kap.                           ---*/
  /*-------------------------------------------------------------------*/
  return (1.0 - kappa) * gradProj + kappa * delta;
}

/*!
 * \brief MUSCL reconstruction of the specified variable.
 * \note The result should be halved when added to i (or subtracted from j).
 * \note Reads its own row of the gradient container (rather than being handed an already
 *       gathered nVarGrad x nDim block, as it once was) so that a caller reconstructing a
 *       single variable, e.g. a scalar with nVar 1, never gathers a Matrix<Double,1,nDim>: that
 *       shape is the same RowMajor, one-row degeneracy that forces EdgeResidual's Size floor
 *       (see numerics/util.hpp), and here it would silently turn a row into a lone scalar
 *       instead of failing to compile, since Matrix<Double,1,nDim> still satisfies IsVector.
 */
template <size_t nDim, class Double, class Gradient_t, class Int = typename CLaneTraits<Double>::Int>
FORCEINLINE Double musclReconstruction(Int iPoint, const Gradient_t& gradient, size_t iRow,
                                       const Vector<Double, nDim>& vector_ij, const Double& delta,
                                       const CNonDeduced<Double>& kappa, const CNonDeduced<Double>& umusclRamp) {
  const auto grad = gatherVariables<nDim>(iPoint, gradient, iRow);
  const Double proj = dot(grad, vector_ij);
  return umusclRamp * umusclProjection(proj, delta, kappa);
}

/*!
 * \brief Unlimited reconstruction.
 * \param[in] iRow - Starting row of gradient to read, for reconstructing a slice of a
 *            larger set of gradients (e.g. only the velocity out of the primitives).
 */
template <size_t nVarGrad_ = 0, size_t nDim, class Double, class VarType, class Gradient_t>
FORCEINLINE void musclUnlimited(typename CLaneTraits<Double>::Int iPoint, typename CLaneTraits<Double>::Int jPoint,
                                const Vector<Double, nDim>& vector_ij, const Gradient_t& gradient, CPair<VarType>& V,
                                const CNonDeduced<Double>& kappa, const CNonDeduced<Double>& umusclRamp,
                                size_t iRow = 0) {
  constexpr auto nVarGrad = nVarGrad_ > 0 ? nVarGrad_ : VarType::nVar;

  for (size_t iVar = 0; iVar < nVarGrad; ++iVar) {
    /*--- Centered difference, needed for U-MUSCL projection ---*/
    const Double delta_ij = V.j.all(iVar) - V.i.all(iVar);

    /*--- U-MUSCL reconstructed variables ---*/
    const Double proj_i = musclReconstruction<nDim>(iPoint, gradient, iRow + iVar, vector_ij, delta_ij, kappa, umusclRamp);
    const Double proj_j = musclReconstruction<nDim>(jPoint, gradient, iRow + iVar, vector_ij, delta_ij, kappa, umusclRamp);

    /*--- Apply reconstruction: V_L = V_i + 0.5 * dV_ij^kap ---*/
    V.i.all(iVar) += 0.5 * proj_i;
    V.j.all(iVar) -= 0.5 * proj_j;
  }
}

/*!
 * \brief Limited reconstruction with point-based limiter.
 */
template <size_t nVarGrad_ = 0, size_t nDim, class Double, class VarType, class Limiter_t, class Gradient_t>
FORCEINLINE void musclPointLimited(typename CLaneTraits<Double>::Int iPoint,
                                   typename CLaneTraits<Double>::Int jPoint, const Vector<Double, nDim>& vector_ij,
                                   const Limiter_t& limiter, const Gradient_t& gradient, CPair<VarType>& V,
                                   const CNonDeduced<Double>& kappa, const CNonDeduced<Double>& umusclRamp,
                                   size_t iRow = 0) {
  constexpr auto nVarGrad = nVarGrad_ > 0 ? nVarGrad_ : VarType::nVar;

  auto lim_i = gatherVariables<nVarGrad>(iPoint, limiter, iRow);
  auto lim_j = gatherVariables<nVarGrad>(jPoint, limiter, iRow);

  for (size_t iVar = 0; iVar < nVarGrad; ++iVar) {
    /*--- Centered difference, needed for U-MUSCL projection ---*/
    const Double delta_ij = V.j.all(iVar) - V.i.all(iVar);

    /*--- U-MUSCL reconstructed variables ---*/
    const Double proj_i = musclReconstruction<nDim>(iPoint, gradient, iRow + iVar, vector_ij, delta_ij, kappa, umusclRamp);
    const Double proj_j = musclReconstruction<nDim>(jPoint, gradient, iRow + iVar, vector_ij, delta_ij, kappa, umusclRamp);

    /*--- Apply reconstruction: V_L = V_i + 0.5 * lim * dV_ij^kap ---*/
    V.i.all(iVar) += 0.5 * lim_i(iVar) * proj_i;
    V.j.all(iVar) -= 0.5 * lim_j(iVar) * proj_j;
  }
}

/*!
 * \brief Limited reconstruction with edge-based limiter.
 */
template <size_t nVarGrad_ = 0, size_t nDim, class Double, class VarType, class Gradient_t>
FORCEINLINE void musclEdgeLimited(typename CLaneTraits<Double>::Int iPoint,
                                  typename CLaneTraits<Double>::Int jPoint, const Vector<Double, nDim>& vector_ij,
                                  const Gradient_t& gradient, CPair<VarType>& V, const CNonDeduced<Double>& kappa,
                                  const CNonDeduced<Double>& umusclRamp, size_t iRow = 0) {
  constexpr auto nVarGrad = nVarGrad_ > 0 ? nVarGrad_ : VarType::nVar;

  for (size_t iVar = 0; iVar < nVarGrad; ++iVar) {
    /*--- Centered difference, needed for U-MUSCL projection and limiter ---*/
    const Double delta_ij = V.j.all(iVar) - V.i.all(iVar);
    const Double delta_ij_2 = pow(delta_ij, 2) + 1e-6;

    /*--- U-MUSCL reconstructed variables ---*/
    const Double proj_i = musclReconstruction<nDim>(iPoint, gradient, iRow + iVar, vector_ij, delta_ij, kappa, umusclRamp);
    const Double proj_j = musclReconstruction<nDim>(jPoint, gradient, iRow + iVar, vector_ij, delta_ij, kappa, umusclRamp);

    const Double lim_i = (delta_ij_2 + proj_i * delta_ij) / (pow(proj_i, 2) + delta_ij_2);
    const Double lim_j = (delta_ij_2 + proj_j * delta_ij) / (pow(proj_j, 2) + delta_ij_2);

    /*--- Apply reconstruction: V_L = V_i + 0.5 * lim * dV_ij^kap ---*/
    V.i.all(iVar) += 0.5 * lim_i * proj_i;
    V.j.all(iVar) -= 0.5 * lim_j * proj_j;
  }
}

/*!
 * \brief Reconstruct a slice of nVarGrad variables starting at column iRow, dispatching on the
 *        limiter type. This is the switch `reconstructPrimitives` used to perform inline; lifted
 *        here so both the flow and the scalar reconstructions call the same body.
 */
template <size_t nVarGrad_ = 0, size_t nDim, class Double, class VarType, class Limiter_t, class Gradient_t>
FORCEINLINE void reconstruct(typename CLaneTraits<Double>::Int iPoint, typename CLaneTraits<Double>::Int jPoint,
                             const Vector<Double, nDim>& vector_ij, const Gradient_t& gradient,
                             const Limiter_t& limiter, LIMITER limiterType, size_t iRow, CPair<VarType>& V,
                             const CNonDeduced<Double>& kappa, const CNonDeduced<Double>& umusclRamp) {
  switch (limiterType) {
    case LIMITER::NONE:
      musclUnlimited<nVarGrad_>(iPoint, jPoint, vector_ij, gradient, V, kappa, umusclRamp, iRow);
      break;
    case LIMITER::VAN_ALBADA_EDGE:
      musclEdgeLimited<nVarGrad_>(iPoint, jPoint, vector_ij, gradient, V, kappa, umusclRamp, iRow);
      break;
    default:
      musclPointLimited<nVarGrad_>(iPoint, jPoint, vector_ij, limiter, gradient, V, kappa, umusclRamp, iRow);
      break;
  }
}

/*!
 * \brief Update the matrix and right-hand-side of a linear system with one conservative flux.
 */
template <size_t nVar, class Double, class Int = typename CLaneTraits<Double>::Int>
FORCEINLINE void updateLinearSystem(Int iEdge, Int iPoint, Int jPoint, bool implicit, UpdateType updateType,
                                    Double updateMask, const Vector<Double, nVar>& flux,
                                    const Matrix<Double, nVar, nVar>& jac_i, const Matrix<Double, nVar, nVar>& jac_j,
                                    CSysVector<su2double>& vector, SparseMatrixType& matrix) {
  if (updateType == UpdateType::COLORING) {
    vector.UpdateBlocks(iPoint, jPoint, flux, updateMask);
    if (implicit) {
      auto wasActive = AD::BeginPassive();
      matrix.SetBlocks(iEdge, iPoint, jPoint, jac_i, jac_j, updateMask);
      AD::EndPassive(wasActive);
    }
  } else {
    vector.SetBlock(iEdge, flux, updateMask);
    if (implicit) {
      auto wasActive = AD::BeginPassive();
      matrix.SetBlocks(iEdge, jac_i, jac_j, updateMask);
      AD::EndPassive(wasActive);
    }
  }
}

/*!
 * \brief Update the matrix and right-hand-side of a linear system with two independent row
 *        contributions and four independent Jacobian blocks.
 * \note It carries a second CSysVector, the target of flux_j under UpdateType::REDUCTION and
 *       unused under COLORING, where both rows are written directly.
 */
template <size_t nVar, class Double, class Int = typename CLaneTraits<Double>::Int>
FORCEINLINE void updateLinearSystem(Int iEdge, Int iPoint, Int jPoint, bool implicit, UpdateType updateType,
                                    Double updateMask, const EdgeResidual<Double, nVar>& res,
                                    CSysVector<su2double>& vector, CSysVector<su2double>& vectorDiff,
                                    SparseMatrixType& matrix) {
  if (updateType == UpdateType::COLORING) {
    vector.AddBlock(iPoint, res.flux_i, updateMask);
    vector.AddBlock(jPoint, res.flux_j, updateMask);
    if (implicit) {
      auto wasActive = AD::BeginPassive();
      matrix.SetBlocks(iEdge, iPoint, jPoint, res.jac_ii, res.jac_ij, res.jac_ji, res.jac_jj, updateMask);
      AD::EndPassive(wasActive);
    }
  } else {
    vector.SetBlock(iEdge, res.flux_i, updateMask);
    vectorDiff.SetBlock(iEdge, res.flux_j, updateMask);
    if (implicit) {
      auto wasActive = AD::BeginPassive();
      matrix.SetOffDiagBlocks(iEdge, res.jac_ij, res.jac_ji, updateMask);
      AD::EndPassive(wasActive);
    }
  }
}

/*!
 * \brief Store the (scalar) mass flux of an edge, e.g. for "bounded scalar" transport equations.
 * \note No-op if "target" is null. As with CEdge's Nodes/Normal, edges within a SIMD group are
 * contiguous (coloring groups are multiples of the SIMD size), so this is a plain vectorized store
 * starting at iEdge[0], relying on "target" being padded to a multiple of the SIMD size.
 */
template <class Double, class Int = typename CLaneTraits<Double>::Int>
FORCEINLINE void updateEdgeMassFlux(Int iEdge, const Double& massFlux, su2activevector* target) {
  if (target) massFlux.store(&(*target)[iEdge[0]]);
}
