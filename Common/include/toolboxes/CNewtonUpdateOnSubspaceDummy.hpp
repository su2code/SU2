/*!
 * \file CNewtonUpdateOnSubspaceDummy.hpp
 * \brief
 * \note
 * \author O. Burghardt
 * \version 7.1.1 "Blackbird"
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

#include "CQuasiNewtonInvLeastSquares.hpp"

/*--- Both classes operate on a window of past corrected solutions (X) and an
 *    data structure (R) of similar size to construct a (quasi-) Newton scheme.
 *    Have to identify a proper common base class. ---*/
template<class Scalar_t>
class CNewtonUpdateOnSubspace : public CQuasiNewtonInvLeastSquares<Scalar_t> {

public:

  using Scalar = Scalar_t;
  using Index = typename su2matrix<Scalar>::Index;

  /*! \brief Default construction without allocation. */
  CNewtonUpdateOnSubspace() = default;

  /*! \brief Construction with allocation, see "resize". */
  CNewtonUpdateOnSubspace(Index nsample, Index npt, Index nvar, Index nptdomain = 0) {}

  /*!
   * \brief Resize the object.
   * \param[in] nsample - Number of samples used to build the FP history.
   * \param[in] nbasis - Dimension of basis the unstable space on which we apply the Newton update scheme.
   * \param[in] npt - Size of the solution including any halos.
   * \param[in] nvar - Number of solution variables.
   * \param[in] nptdomain - Local size (< npt), if 0 (default), MPI parallelization is skipped.
   */
  void resize(Index nsample, Index nbasis, Index npt, Index nvar, Index nptdomain = 0) {}

  /*! \brief Size of the object, the size of the subspace basis. */
  Index size() const { return 0; }

  /*! \brief Discard all history, keeping the current sample. */
  void reset() {}

  /*!
   * \brief Check for new basis vector and eventually append to basis.
   */
  bool checkBasis(su2double KrylovCriterionValue) { return false; }

  /*!
   * \brief Compute new projected subspace Jacobian and the inverse matrix for Newton steps.
   * \note To be used directly after basis dimension has been increased.
   */
  void computeProjectedJacobian(unsigned short iZone, su2matrix<int>& InputIndices, su2matrix<int>& OutputIndices) {}

  void computeProjections() {}

  void shift(bool outer) {}

  /*!
   * \brief Compute a new approximation.
   * \note To be used after storing the FP result.
   */
  void computeCorrection() {}
};
