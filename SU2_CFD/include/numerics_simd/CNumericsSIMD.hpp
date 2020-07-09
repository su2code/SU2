/*!
 * \file CNumericsSIMD.hpp
 * \brief Vectorized (SIMD) numerics classes.
 * \author P. Gomes
 * \version 7.0.6 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/parallelization/vectorization.hpp"

/*!
 * \enum UpdateType
 * \brief Ways to update vectors and system matrices.
 * COLORING is the typical i/j update, whereas for REDUCTION
 * the fluxes are stored and the matrix diagonal is not modified.
 */
enum class UpdateType {COLORING, REDUCTION};

/*!
 * \brief Define Double and Int types, SIMD for normal mode, scalar for AD.
 */
#if !defined(CODI_REVERSE_TYPE) && !defined(CODI_FORWARD_TYPE)
/*--- Scalar types (SIMD). ---*/
using Double = simd::Array<su2double>;
using Int = simd::Array<unsigned long, Double::Size>;
constexpr size_t IntSize = Int::Size;
#else
/*--- No SIMD with AD. ---*/
using Double = su2double;
using Int = unsigned long;
constexpr size_t IntSize = 1;
#endif

template<class T> class CSysVector;
template<class T> class CSysMatrix;
class CConfig;
class CGeometry;
class CVariable;

#ifdef CODI_FORWARD_TYPE
using SparseMatrixType = CSysMatrix<su2double>;
#else
using SparseMatrixType = CSysMatrix<su2mixedfloat>;
#endif

/*!
 * \class CNumericsSIMD
 * \brief Base class to define the interface.
 * \note See CNumericsEmptyDecorator.
 */
class CNumericsSIMD {
public:
  /*!
   * \brief Interface for edge flux computation.
   * \param[in] iEdge - The edges for flux computation.
   * \param[in] config - Problem definitions.
   * \param[in] geometry - Problem geometry.
   * \param[in] solution - Solution variables.
   * \param[in] updateType - Type of update done on vector and matrix.
   * \param[in] updateMask - SIMD array of 1's and 0's, the latter prevent the update.
   * \param[in,out] vector - Target for the fluxes.
   * \param[in,out] matrix - Target for the flux Jacobians.
   * \note The update mask is used to handle "remainder" edges (nEdge mod simdSize).
   */
  virtual void ComputeFlux(Int iEdge,
                           const CConfig& config,
                           const CGeometry& geometry,
                           const CVariable& solution,
                           UpdateType updateType,
                           Double updateMask,
                           CSysVector<su2double>& vector,
                           SparseMatrixType& matrix) const = 0;

  /*! \brief Destructor of the class. */
  virtual ~CNumericsSIMD(void) = default;

  /*!
   * \brief Factory method.
   * \param[in] config - Problem definitions.
   * \param[in] nDim - 2D or 3D.
   * \param[in] iMesh - Grid index.
   */
  static CNumericsSIMD* CreateNumerics(const CConfig& config, int nDim, int iMesh);

};
