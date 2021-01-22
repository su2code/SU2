/*!
 * \file CFEMStandardQuadVolumeSol.hpp
 * \brief Class for the FEM quadrilateral standard element for the volume solution.
 *        The functions are in the <i>CFEMStandardQuadVolumeSol.cpp</i> file.
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

#include "CFEMStandardQuadBase.hpp"

/*!
 * \class CFEMStandardQuadVolumeSol
 * \brief Class which defines the variables and methods for the
 *        quadrilateral standard element for the volume solution.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 */
class CFEMStandardQuadVolumeSol final: public CFEMStandardQuadBase {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardQuadVolumeSol() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly       - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact  - Polynomial order that must be integrated exactly
   *                              by the integration rule.
   * \param[in] val_locGridDOFs - Location of the grid DOFs (LGL or Equidistant).
   * \param[in] val_nVar        - Number of variables in the jitted gemm calls (not used).
   */
  CFEMStandardQuadVolumeSol(const unsigned short val_nPoly,
                            const unsigned short val_orderExact,
                            const unsigned short val_locGridDOFs,
                            const unsigned short val_nVar);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFEMStandardQuadVolumeSol() = default;

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

  vector<passivedouble> rLineSolDOFs; /*!< \brief 1D parametric coordinates of the nodal solution DOFs. */

  ColMajorMatrix<passivedouble> legBasisLineInt;    /*!< \brief The values of the 1D Legendre basis functions
                                                                in the integration points. */
  ColMajorMatrix<passivedouble> derLegBasisLineInt; /*!< \brief The values of the derivatives of the 1D Legendre
                                                                basis functions in the integration points. */
  ColMajorMatrix<passivedouble> hesLegBasisLineInt; /*!< \brief The values of the 2nd derivatives of the 1D Legendre
                                                                basis functions in the integration points. */

  ColMajorMatrix<passivedouble> legBasisLineSolDOFs;    /*!< \brief The values of the 1D Legendre basis functions
                                                                    in the solution DOFs. */

  ColMajorMatrix<passivedouble> derLegBasisLineSolDOFs; /*!< \brief The values of the derivatives of the 1D Legendre
                                                                    basis functions in the solution DOFs. */

  TPI2D TensorProductDataVolIntPoints = nullptr; /*!< \brief Function pointer to carry out the tensor product
                                                             to compute the data in the volume integration points. */
  TPI2D TensorProductDataVolSolDOFs   = nullptr; /*!< \brief Function pointer to carry out the tensor product
                                                             to compute the data in the volume nodal solution DOFs. */
};
