/*!
 * \file CFEMStandardTriGrid.hpp
 * \brief Class for the FEM triangle standard element for the grid.
 *        The functions are in the <i>CFEMStandardTriGrid.cpp</i> file.
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

#include "CFEMStandardTriBase.hpp"

/*!
 * \class CFEMStandardTriGrid
 * \brief Class which defines the variables and methods for the
 *        triangle standard element for the grid.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 */
class CFEMStandardTriGrid final: public CFEMStandardTriBase {

public:
  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardTriGrid() = delete;

  /*!
   * \overload
   * \param[in] val_nPolyGrid   - Polynomial degree of the grid for this element.
   * \param[in] val_nPolyGrid   - Polynomial degree of the solution for this element.
   * \param[in] val_orderExact  - Polynomial order that must be integrated exactly
   *                              by the integration rule.
   * \param[in] val_locGridDOFs - Location of the grid DOFS, either LGL or equidistant.
   */
  CFEMStandardTriGrid(const unsigned short val_nPolyGrid,
                      const unsigned short val_nPolySol,
                      const unsigned short val_orderExact,
                      const unsigned short val_locGridDOFs);

  /*!
   * \brief Destructor.
   */
  ~CFEMStandardTriGrid();

  /*-----------------------------------------------------------------------------------*/
  /*---                  Inline public member functions.                            ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, that returns the number of solution DOFs.
   * \return The number of solution DOFs of the volume element.
   */
  inline unsigned short GetNSolDOFs(void) const override {return rTriangleSolDOFs.size();}

  /*!
   * \brief Function, that returns the padded number of solution DOFs.
   * \return The padded number of solution DOFs of the volume element.
   */
  inline unsigned short GetNSolDOFsPad(void) const override {
    return PaddedValue(rTriangleSolDOFs.size());
  }

  /*-----------------------------------------------------------------------------------*/
  /*---                     Public member functions.                                ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which computes the coordinates in the integration points.
   * \param[in]  notUsed    - Argument present to be consistent with the base class
   *                          function, which is overwritten.
   * \param[in]  matCoorDOF - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matCoorInt - Matrix that contains the coordinates of the integration
   *                          points.
   */
  void CoorIntPoints(const bool                notUsed,
                     ColMajorMatrix<su2double> &matCoorDOF,
                     ColMajorMatrix<su2double> &matCoorInt) override;

  /*!
   * \brief Function, which computes the coordinates of the solution DOFs. These are
   *        only different from the grid DOFs when different polynomials degrees
   *        are used for the grid and solution.
   * \param[in]  matCoorDOF    - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matCoorSolDOF - Matrix that contains the coordinates of the solution DOFs.
   */
  void CoorSolDOFs(ColMajorMatrix<su2double> &matCoorDOF,
                   ColMajorMatrix<su2double> &matCoorSolDOF) override;

  /*!
   * \brief Function, which computes the derivatives of the coordinates in the
   *        integration points.
   * \param[in]  notUsed    - Argument present to be consistent with the base class
   *                          function, which is overwritten.
   * \param[in]  matCoor    - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matDerCoor - Vector of matrices to store the derivatives of the coordinates.
   */
  void DerivativesCoorIntPoints(const bool                         notUsed,
                                ColMajorMatrix<su2double>          &matCoor,
                                vector<ColMajorMatrix<su2double> > &matDerCoor) override;

  /*!
   * \brief Function, which computes the 2nd derivatives of the coordinates in the
   *        integration points.
   * \param[in]  matCoor       - Matrix that contains the coordinates of the grid DOFs
   * \param[out] matDer2ndCoor - Vector of matrices to store the 2nd derivatives of the coordinates.
   */
  void Derivatives2ndCoorIntPoints(ColMajorMatrix<su2double>          &matCoor,
                                   vector<ColMajorMatrix<su2double> > &matDer2ndCoor) override;

  /*!
   * \brief Function, which computes the derivatives of the coordinates in the
   *        solution DOFs.
   * \param[in]  matCoor    - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matDerCoor - Vector of matrices to store the derivatives of the coordinates.
   */
  void DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                              vector<ColMajorMatrix<su2double> > &matDerCoor) override;

  /*!
   * \brief Function, which evaluates the coordinates and gradients for the given
   *        parametric coordinates.
   * \param[in]  matCoor - Matrix that contains the coordinates of the grid DOFs.
   * \param[in]  par     - Parametric coordinates for which the physical coordinates
   *                       and derivatives must be determined.
   * \param[out] x       - Physical coordinates to be determined.
   * \param[out] dxdpar  - Derivatives of the physical coordinates w.r.t. the parametric
   *                       coordinates.
   */
  void EvalCoorAndGradCoor(ColMajorMatrix<su2double> &matCoor,
                           const su2double           *par,
                           su2double                 x[3],
                           su2double                 dxdpar[3][3]) override;

  /*!
   * \brief Function, which interpolates the parametric coordinate within the
   *        given linear sub-element.
   * \param[in]  subElem - Sub-element in which the coordinate must be interpolated.
   * \param[in]  weights - Interpolation weights in the linear sub-element.
   * \param[out] parCoor - Parametric coordinates to be interpolated.
   */
  void InterpolCoorSubElem(const unsigned short subElem,
                           const su2double      *weights,
                           su2double            *parCoor) override;

private:

  vector<passivedouble> rTriangleDOFs; /*!< \brief Parametric r-coordinates of the triangle grid DOFs. */
  vector<passivedouble> sTriangleDOFs; /*!< \brief Parametric s-coordinates of the triangle grid DOFs. */

  vector<passivedouble> rTriangleSolDOFs; /*!< \brief Parametric r-coordinates of the triangle solution DOFs. */
  vector<passivedouble> sTriangleSolDOFs; /*!< \brief Parametric s-coordinates of the triangle solution DOFs. */

  ColMajorMatrix<passivedouble> lagBasisInt;             /*!< \brief The values of the Lagrangian basis functions
                                                                     in the integration points. */
  vector<ColMajorMatrix<passivedouble> > derLagBasisInt; /*!< \brief The values of the derivatives of the Lagrangian basis
                                                                     functions in the integration points. It is a vector,
                                                                     because there are derivatives in two directions. */
  vector<ColMajorMatrix<passivedouble> > hesLagBasisInt; /*!< \brief The values of the 2nd derivatives, Hessian, of the
                                                                     Lagrangian basis functions in the integration points
                                                                     for the point distribution in use. It is a vector,
                                                                     because there are 3 2nd derivatives in 2D. */

  ColMajorMatrix<passivedouble> lagBasisSolDOFs;             /*!< \brief The values of the Lagrangian basis functions
                                                                         in the nodal solution DOFs. */
  vector<ColMajorMatrix<passivedouble> > derLagBasisSolDOFs; /*!< \brief The values of the derivatives of the Lagrangian
                                                                         basis functions in the nodal solution DOFs.
                                                                         It is a vector, because there are derivatives
                                                                         in two directions. */

  void *jitterDOFs2Int = nullptr;      /*!< \brief Pointer to the data for the jitted gemm function
                                                   to compute data in the integration points. */
  dgemm_jit_kernel_t gemmDOFs2Int;     /*!< \brief Pointer to the function to carry out jitterDOFs2Int. */

  void *jitterDOFs2SolDOFs = nullptr;  /*!< \brief Pointer to the data for the jitted gemm function
                                                   to compute data in the solution DOFs. */
  dgemm_jit_kernel_t gemmDOFs2SolDOFs; /*!< \brief Pointer to the function to carry out jitterDOFs2SolDOFs. */
};
