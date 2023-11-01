/*!
 * \file CRadialBasisFunction.hpp
 * \brief Radial basis function interpolation.
 * \author Joel Ho, P. Gomes
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

#include "CInterpolator.hpp"
#include "../option_structure.hpp"
#include "../containers/C2DContainer.hpp"

/*!
 * \brief Radial basis function interpolation.
 * \ingroup Interfaces
 */
class CRadialBasisFunction final : public CInterpolator {
  static_assert(su2passivematrix::IsRowMajor, "This class relies on row major storage throughout.");

 private:
  unsigned long MinDonors = 0, AvgDonors = 0, MaxDonors = 0;
  passivedouble Density = 0.0, AvgCorrection = 0.0, MaxCorrection = 0.0;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - index of the donor zone.
   * \param[in] jZone - index of the target zone.
   */
  CRadialBasisFunction(CGeometry**** geometry_container, const CConfig* const* config, unsigned int iZone,
                       unsigned int jZone);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] config - Definition of the particular problem.
   */
  void SetTransferCoeff(const CConfig* const* config) override;

  /*!
   * \brief Print information about the interpolation.
   */
  void PrintStatistics(void) const override;

  /*!
   * \brief Compute the value of a radial basis function, this is static so it can be re-used.
   * \param[in] type - of radial basis function
   * \param[in] radius - the characteristic dimension
   * \param[in] dist - distance
   * \return value of the RBF.
   */
  static su2double Get_RadialBasisValue(RADIAL_BASIS type, const su2double radius, const su2double dist);

  /*!
   * \brief Compute the RBF "generator" matrix with or without polynomial terms.
   * \note Multiplying C_inv_trunc by a column vector gives specific coefficients for given "known values",
   * conversely, multiplying (on the left) by a row vector of polynomial and RBF values gives generic
   * interpolation coefficients for a given target evaluation point.
   * \param[in] type - Type of radial basis function.
   * \param[in] usePolynomial - Whether to use polynomial terms.
   * \param[in] radius - Normalizes point-to-point distance when computing RBF values.
   * \param[in] coords - Coordinates of the donor points.
   * \param[out] nPolynomial - Num of poly terms, -1 if !usePolynomial, nDim-1 if coords lie on plane, else nDim.
   * \param[out] keepPolynomialRow - Size nDim, signals which (if any) iDim was removed from polynomial term.
   * \param[out] C_inv_trunc - The generator matrix as described above.
   */
  static void ComputeGeneratorMatrix(RADIAL_BASIS type, bool usePolynomial, su2double radius,
                                     const su2activematrix& coords, int& nPolynomial, vector<int>& keepPolynomialRow,
                                     su2passivematrix& C_inv_trunc);

  /*!
   * \brief If the polynomial term is included in the interpolation, and the points lie on a plane, the matrix
   * becomes rank deficient and cannot be inverted. This method detects that condition and corrects it by
   * removing a row from P (the polynomial part of the interpolation matrix).
   * \param[in] max_diff_tol - Tolerance to detect whether points are on a plane.
   * \param[out] keep_row - Marks the dimensions of P kept.
   * \param[in,out] P - Polynomial part of the interpolation matrix, one row may be eliminated.
   * \return n_polynomial - Size of the polynomial part on exit (in practice nDim or nDim-1).
   */
  static int CheckPolynomialTerms(su2double max_diff_tol, vector<int>& keep_row, su2passivematrix& P);

 private:
  /*!
   * \brief Helper function, prunes (by setting to zero) small interpolation coefficients,
   * i.e. <= tolerance*max(abs(coeffs)). The vector is re-scaled such that sum(coeffs)==1.
   * \param[in] tolerance - Relative pruning tolerance.
   * \param[in] size - Size of the coefficient vector.
   * \param[in,out] coeffs - Iterator to start of vector of interpolation coefficients.
   * \return Number of non-zero coefficients after pruning and correction factor.
   */
  template <typename Float, typename Int, class ForwardIt>
  static pair<Int, Float> PruneSmallCoefficients(Float tolerance, Int size, ForwardIt coeffs) {
    /*--- Determine the pruning threshold. ---*/
    Float thresh = 0.0;
    auto end = coeffs;
    for (Int i = 0; i < size; ++i) thresh = max(thresh, fabs(*(end++)));
    thresh *= tolerance;

    /*--- Prune and count non-zeros. ---*/
    Int numNonZeros = 0;
    Float coeffSum = 0.0;
    for (auto it = coeffs; it != end; ++it) {
      if (fabs(*it) > thresh) {  // keep
        coeffSum += *it;
        ++numNonZeros;
      } else {
        *it = 0.0;
      }  // prune
    }

    /*--- Correct remaining coefficients, sum must be 1 for conservation. ---*/
    Float correction = 1.0 / coeffSum;
    while (coeffs != end) *(coeffs++) *= correction;

    return make_pair(numNonZeros, correction);
  }
};
