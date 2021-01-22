/*!
 * \file CFEMStandardPyraVolumeSol.hpp
 * \brief Class for the FEM pyramid standard element for the volume solution
 *        The functions are in the <i>CFEMStandardPyraVolumeSol.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
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

#include "CFEMStandardPyraBase.hpp"

/*!
 * \class CFEMStandardPyraVolumeSol
 * \brief Class which defines the variables and methods for the
 *        pyramid standard element for the volume solution.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 */
class CFEMStandardPyraVolumeSol final: public CFEMStandardPyraBase {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardPyraVolumeSol() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly       - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact  - Polynomial order that must be integrated exactly
   *                              by the integration rule.
   * \param[in] val_locGridDOFs - Location of the grid DOFs (LGL or Equidistant).
   * \param[in] val_nVar        - Number of variables in the jitted gemm calls.
   */
  CFEMStandardPyraVolumeSol(const unsigned short val_nPoly,
                            const unsigned short val_orderExact,
                            const unsigned short val_locGridDOFs,
                            const unsigned short val_nVar);

  /*!
   * \brief Destructor.
   */
  ~CFEMStandardPyraVolumeSol();

  /*-----------------------------------------------------------------------------------*/
  /*---                     Public member functions.                                ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which determines the basis functions for the given parametric
   *        coordinates.
   * \param[in]  parCoor  - Double vector that contains the parametric coordinates
   *                        for which the basis functions must be determined.
   * \param[out] matBasis - Matrix that contains the values of the basis functions
   *                        in the given parametric coordinates.
   */
  void BasisFunctionsInPoints(const vector<vector<passivedouble> > &parCoor,
                              ColMajorMatrix<passivedouble>        &matBasis) override;

private:
  vector<passivedouble> rPyraSolDOFs; /*!< \brief Parametric r-coordinates of the pyramid solution DOFs. */
  vector<passivedouble> sPyraSolDOFs; /*!< \brief Parametric s-coordinates of the pyramid solution DOFs. */
  vector<passivedouble> tPyraSolDOFs; /*!< \brief Parametric t-coordinates of the pyramid solution DOFs. */

  ColMajorMatrix<passivedouble> legBasisInt;             /*!< \brief The values of the Legendre basis functions
                                                                     in the integration points. */
  vector<ColMajorMatrix<passivedouble> > derLegBasisInt; /*!< \brief The values of the derivatives of the Legendre
                                                                     basis functions in the integration points.
                                                                     It is a vector, because there are derivatives
                                                                     in three directions. */
  vector<ColMajorMatrix<passivedouble> > hesLegBasisInt; /*!< \brief The values of the 2nd derivatives of the Legendre
                                                                     basis functions in the integration points.
                                                                     It is a vector, because there are 6 2nd derivatives
                                                                     in three space dimensions. */

  ColMajorMatrix<passivedouble> legBasisSolDOFs;             /*!< \brief The values of the Legendre basis functions
                                                                         in the solution DOFs. */
  vector<ColMajorMatrix<passivedouble> > derLegBasisSolDOFs; /*!< \brief The values of the derivatives of the
                                                                         Legendre basis functions in the solution
                                                                         DOFs. It is a vector, because there
                                                                         are derivatives in three directions. */

  void *jitterDOFs2Int = nullptr;      /*!< \brief Pointer to the data for the jitted gemm function
                                                   to compute data in the integration points. */
  dgemm_jit_kernel_t gemmDOFs2Int;     /*!< \brief Pointer to the function to carry out jitterDOFs2Int. */

  void *jitterDOFs2SolDOFs = nullptr;  /*!< \brief Pointer to the data for the jitted gemm function
                                                   to compute data in the solution DOFs. */
  dgemm_jit_kernel_t gemmDOFs2SolDOFs; /*!< \brief Pointer to the function to carry out jitterDOFs2SolDOFs. */
};
