/*!
 * \file util.hpp
 * \brief Generic auxiliary functions.
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

#pragma once

#include "CNumericsSIMD.hpp"
#include "../../../Common/include/containers/C2DContainer.hpp"
#include "../../../Common/include/linear_algebra/CSysVector.hpp"
#include "../../../Common/include/linear_algebra/CSysMatrix.hpp"

/*!
 * \brief Static vector and matrix types.
 * \note These should be used instead of C-style arrays.
 */
template<class Type, size_t Size>
using Vector = C2DContainer<unsigned long, Type, StorageType::ColumnMajor, Type::Align, Size, 1>;

template<size_t Size> using VectorInt = Vector<Int, Size>;
template<size_t Size> using VectorDbl = Vector<Double, Size>;

template<class Type, size_t Rows, size_t Cols>
using Matrix = C2DContainer<unsigned long, Type, StorageType::RowMajor, Type::Align, Rows, Cols>;

template<size_t Rows, size_t Cols = Rows> using MatrixInt = Matrix<Int, Rows, Cols>;
template<size_t Rows, size_t Cols = Rows> using MatrixDbl = Matrix<Double, Rows, Cols>;

/*!
 * \brief Constexpr version of max.
 */
inline constexpr size_t Max(size_t a, size_t b) { return a>b? a : b; }

/*!
 * \brief Simple pair type for i/j variables.
 */
template<class T>
struct CPair {
  T i, j;
};

/*!
 * \brief Dot product.
 */
template<size_t nDim, class ForwardIterator, class T>
FORCEINLINE Double dot(ForwardIterator iterator, const T* ptr) {
  Double sum = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    sum += *(iterator++) * ptr[iDim];
  }
  return sum;
}

/*!
 * \overload Dot product.
 */
template<size_t nDim, class ForwardIterator>
FORCEINLINE Double dot(ForwardIterator iterator, const VectorDbl<nDim>& vector) {
  return dot<nDim>(iterator, vector.data());
}

/*!
 * \overload Dot product.
 */
template<size_t nDim>
FORCEINLINE Double dot(const VectorDbl<nDim>& a, const VectorDbl<nDim>& b) {
  return dot<nDim>(a.data(), b.data());
}

/*!
 * \brief Squared norm.
 */
template<size_t nDim, class ForwardIterator>
FORCEINLINE Double squaredNorm(ForwardIterator iterator) {
  Double sum = 0.0;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    sum += pow(*(iterator++),2);
  }
  return sum;
}

/*!
 * \overload Squared norm.
 */
template<size_t nDim>
FORCEINLINE Double squaredNorm(const VectorDbl<nDim>& vector) {
  return squaredNorm<nDim>(vector.data());
}

/*!
 * \brief Tangential projection.
 */
template<size_t nDim>
FORCEINLINE VectorDbl<nDim> tangentProjection(const MatrixDbl<nDim>& tensor,
                                              const VectorDbl<nDim>& unitVector) {
  VectorDbl<nDim> proj;
  for (size_t iDim = 0; iDim < nDim; ++iDim)
    proj(iDim) = dot(tensor[iDim], unitVector);

  Double normalProj = dot(proj, unitVector);

  for (size_t iDim = 0; iDim < nDim; ++iDim)
    proj(iDim) -= normalProj * unitVector(iDim);

  return proj;
}

/*!
 * \brief Vector norm.
 */
template<size_t nDim>
FORCEINLINE Double norm(const VectorDbl<nDim>& vector) { return sqrt(squaredNorm(vector)); }

#ifndef CODI_REVERSE_TYPE
/*!
 * \brief Gather a single variable from index iPoint of a 1D container.
 */
template<class Container>
FORCEINLINE Double gatherVariables(Int iPoint, const Container& vars) {
  return *vars.innerIter(iPoint);
}

/*!
 * \brief Gather a vector of variables (size nVar) from row iPoint of a 2D container.
 */
template<size_t nVar, class Container>
FORCEINLINE VectorDbl<nVar> gatherVariables(Int iPoint, const Container& vars) {
  return vars.template get<VectorDbl<nVar> >(iPoint);
}

/*!
 * \brief Gather a matrix of variables from outer index iPoint of a 3D container.
 */
template<size_t nRows, size_t nCols, class Container>
FORCEINLINE MatrixDbl<nRows,nCols> gatherVariables(Int iPoint, const Container& vars) {
  return vars.template get<MatrixDbl<nRows,nCols> >(iPoint);
}
#else

namespace {
  template<class Container, su2enable_if<Container::IsVector> = 0>
  FORCEINLINE const su2double& get(const Container& vars, unsigned long iPoint) { return vars(iPoint); }

  /*--- When getting 1 variable from a matrix container, we assume it is the first. ---*/
  template<class Container, su2enable_if<!Container::IsVector> = 0>
  FORCEINLINE const su2double& get(const Container& vars, unsigned long iPoint) { return vars(iPoint,0); }
}

template<class Container>
FORCEINLINE Double gatherVariables(Int iPoint, const Container& vars) {
  Double x;
  for (size_t k=0; k<Double::Size; ++k) {
    AD::SetPreaccIn(get(vars, iPoint[k]));
    x[k] = get(vars, iPoint[k]);
  }
  return x;
}

template<size_t nVar, class Container>
FORCEINLINE VectorDbl<nVar> gatherVariables(Int iPoint, const Container& vars) {
  VectorDbl<nVar> x;
  for (size_t i=0; i<nVar; ++i) {
    for (size_t k=0; k<Double::Size; ++k) {
      AD::SetPreaccIn(vars(iPoint[k],i));
      x[i][k] = vars(iPoint[k],i);
    }
  }
  return x;
}

template<size_t nRows, size_t nCols, class Container>
FORCEINLINE MatrixDbl<nRows,nCols> gatherVariables(Int iPoint, const Container& vars) {
  MatrixDbl<nRows,nCols> x;
  for (size_t i=0; i<nRows; ++i) {
    for (size_t j=0; j<nCols; ++j) {
      for (size_t k=0; k<Double::Size; ++k) {
        AD::SetPreaccIn(vars(iPoint[k],i,j));
        x(i,j)[k] = vars(iPoint[k],i,j);
      }
    }
  }
  return x;
}
#endif

/*!
 * \brief Stop the AD preaccumulation.
 */
template<size_t nVar>
FORCEINLINE void stopPreacc(VectorDbl<nVar>& x) {
  AD::SetPreaccOut(x, nVar, Double::Size);
  AD::EndPreacc();
}

/*!
 * \brief Distance vector, from point i to point j.
 */
template<size_t nDim, class Container>
FORCEINLINE VectorDbl<nDim> distanceVector(Int iPoint, Int jPoint,
                                           const Container& coords) {
  auto coord_i = gatherVariables<nDim>(iPoint, coords);
  auto coord_j = gatherVariables<nDim>(jPoint, coords);
  VectorDbl<nDim> vector_ij;
  for (size_t iDim = 0; iDim < nDim; ++iDim) {
    vector_ij(iDim) = coord_j(iDim) - coord_i(iDim);
  }
  return vector_ij;
}

/*!
 * \brief Update the matrix and right-hand-side of a linear system.
 */
template<size_t nVar>
FORCEINLINE void updateLinearSystem(Int iEdge,
                                    Int iPoint,
                                    Int jPoint,
                                    bool implicit,
                                    UpdateType updateType,
                                    Double updateMask,
                                    const VectorDbl<nVar>& flux,
                                    const MatrixDbl<nVar>& jac_i,
                                    const MatrixDbl<nVar>& jac_j,
                                    CSysVector<su2double>& vector,
                                    SparseMatrixType& matrix) {
  if (updateType == UpdateType::COLORING) {
    vector.UpdateBlocks(iPoint, jPoint, flux, updateMask);
    if(implicit) {
      auto wasActive = AD::BeginPassive();
      matrix.UpdateBlocks(iEdge, iPoint, jPoint, jac_i, jac_j, updateMask);
      AD::EndPassive(wasActive);
    }
  }
  else {
    vector.SetBlock(iEdge, flux, updateMask);
    if(implicit) {
      auto wasActive = AD::BeginPassive();
      matrix.SetBlocks(iEdge, jac_i, jac_j, updateMask);
      AD::EndPassive(wasActive);
    }
  }
}
