/*!
 * \file CRadialBasisFunction.hpp
 * \brief Radial basis function interpolation.
 * \author Joel Ho, P. Gomes
 * \version 7.0.2 "Blackbird"
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

#include "CInterpolator.hpp"
#include "../option_structure.hpp"
#include "../toolboxes/C2DContainer.hpp"

/*!
 * \brief Radial basis function interpolation.
 */
class CRadialBasisFunction final : public CInterpolator {
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - index of the donor zone
   * \param[in] jZone - index of the target zone
   */
  CRadialBasisFunction(CGeometry ****geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone);

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \param[in] config - Definition of the particular problem.
   */
  void Set_TransferCoeff(CConfig **config) override;

  /*!
   * \brief Compute the value of a radial basis function, this is static so it can be re-used.
   * \param[in] type - of radial basis function
   * \param[in] radius - the characteristic dimension
   * \param[in] dist - distance
   * \return value of the RBF.
   */
  static su2double Get_RadialBasisValue(ENUM_RADIALBASIS type, const su2double radius, const su2double dist);

private:
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
  void ComputeGeneratorMatrix(ENUM_RADIALBASIS type, bool usePolynomial, su2double radius,
                              const su2activematrix& coords, int& nPolynomial,
                              vector<int>& keepPolynomialRow, su2passivematrix& C_inv_trunc) const;

  /*!
   * \brief If the polynomial term is included in the interpolation, and the points lie on a plane, the matrix
   * becomes rank deficient and cannot be inverted. This method detects that condition and corrects it by
   * removing a row from P (the polynomial part of the interpolation matrix).
   * \param[in] max_diff_tol - Tolerance to detect whether points are on a plane.
   * \param[out] keep_row - Marks the dimensions of P kept.
   * \param[in,out] P - Polynomial part of the interpolation matrix, one row may be eliminated.
   * \return n_polynomial - Size of the polynomial part on exit (in practice nDim or nDim-1).
   */
  int CheckPolynomialTerms(su2double max_diff_tol, vector<int>& keep_row, su2passivematrix &P) const;

  /*!
   * \brief Prunes (by setting to zero) small interpolation coefficients, i.e.
   * <= tolerance*max(abs(coeffs)). The vector is re-scaled such that sum(coeffs)==1.
   * \param[in] tolerance - Relative pruning tolerance.
   * \param[in,out] coeffs - The vector of interpolation coefficients.
   * \return Number of non-zero coefficients after pruning.
   */
  int PruneSmallCoefficients(passivedouble tolerance, su2passivevector& coeffs) const;

};

/*!
 * \brief Helper class used by CRadialBasisFunction to compute the interpolation weights.
 * The matrix is symmetric but full storage is used as that gives much better performance
 * for some BLAS libraries (notably OpenBLAS). The code should be compiled with LAPACK
 * to use optimized matrix inversion and multiplication routines.
 */
class CSymmetricMatrix {
private:
  enum DecompositionType { NONE, CHOLESKY, LU };

  vector<passivedouble> val_vec, decomp_vec;
  vector<int> perm_vec;
  int sz = 0;
  bool initialized = false;
  DecompositionType decomposed = NONE;

  inline void CheckBounds(int i, int j) const {
    assert(initialized && "Matrix not initialized.");
    assert(i>=0 && i<sz && j>=0 && j<sz && "Index to access matrix out of bounds.");
  }

  inline int IdxFull(int i, int j) const {CheckBounds(i,j); return i*sz + j;}

  inline int IdxSym(int i, int j) const {return IdxFull(min(i,j), max(i,j));}

  inline passivedouble& decomp(int i, int j) { return decomp_vec[IdxFull(i,j)]; }

  // Not optimized dense matrix factorization and inversion for portability.
  void CholeskyDecompose();
  void LUDecompose();
  void CalcInv();
  // Matrix inversion using LAPACK routines (LDLT and LLT factorization).
  void CalcInv_sytri();
  void CalcInv_potri();

public:
  CSymmetricMatrix() = default;
  CSymmetricMatrix(int N) {Initialize(N);}

  void Initialize(int N);

  inline int GetSize() const { return sz; }

  inline passivedouble Get(int i, int j) const { return val_vec[IdxSym(i,j)]; }

  inline void Set(int i, int j, passivedouble val) { val_vec[IdxSym(i,j)] = val; }

  inline passivedouble& operator() (int i, int j) { return val_vec[IdxSym(i,j)]; }

  inline const passivedouble& operator() (int i, int j) const { return val_vec[IdxSym(i,j)]; }

  void MatVecMult(passivedouble *v) const;

  void MatMatMult(const char side, su2passivematrix& mat_in, su2passivematrix& mat_out);

  void Invert(const bool is_spd);

};
