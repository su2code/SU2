/*!
 * \file CFEMStandardQuadVolumeSol.hpp
 * \brief Class for the FEM quadrilateral standard element for the volume solution.
 *        The functions are in the <i>CFEMStandardQuadVolumeSol.cpp</i> file.
 * \author E. van der Weide
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

#include "CFEMStandardQuadBase.hpp"

/*!
 * \class CFEMStandardQuadVolumeSol
 * \brief Class which defines the variables and methods for the
 *        quadrilateral standard element for the volume solution.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
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
  /*---                  Inline public member functions.                            ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, that makes available the correction factor for the inviscid
   *        spectral radius.
   * \return The correction factor for the inviscid spectral radius.
   */
  inline passivedouble GetFactorInviscidSpectralRadius(void) const override {return factInviscidRad;}

  /*!
   * \brief Function, that makes available the correction factor for the viscous
   *        spectral radius.
   * \return The correction factor for the viscous spectral radius.
   */
  inline passivedouble GetFactorViscousSpectralRadius(void) const override {return factViscousRad;}

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

  /*!
   * \brief Function, that converts the modal form of the DOFs to the nodal form.
   * \param[in,out] solDOFs - On entry it contains the modal solution in the DOFs,
   *                          on exit it contains the nodal solution.
   */
  void ModalToNodal(ColMajorMatrix<su2double> &solDOFs) override;

  /*!
   * \brief Function, that converts the nodal form of the DOFs to the modal form.
   * \param[in,out] solDOFs - On entry it contains the nodal solution in the DOFs,
   *                          on exit it contains the modal solution.
   */
  void NodalToModal(ColMajorMatrix<su2double> &solDOFs) override;

  /*!
   * \brief Function that computes the gradients of the solution in integration points.
   * \param[in]  matSolDOF     - Matrix that contains the modal solution DOFs.
   * \param[out] matGradSolInt - Vector of matrices the contains the gradients of the
   *                             solution in the integration points.
   */
  void GradSolIntPoints(ColMajorMatrix<su2double>          &matSolDOF,
                        vector<ColMajorMatrix<su2double> > &matGradSolInt) override;

  /*!
   * \brief Function that computes the solution in integration points.
   * \param[in]  matSolDOF - Matrix that contains the modal solution DOFs.
   * \param[out] matSolInt - Matrix that contains the solution in the integration points.
   */
  void SolIntPoints(ColMajorMatrix<su2double> &matSolDOF,
                    ColMajorMatrix<su2double> &matSolInt) override;

  /*!
   * \brief Function that computes the solution in integration points
   *        from the padded modal solution.
   * \param[in]  matSolDOF - Matrix that contains the modal solution DOFs, the number
   *                         DOFs are padded.
   * \param[out] matSolInt - Matrix that contains the solution in the integration points.
   */
  void SolIntPointsDOFsPadded(ColMajorMatrix<su2double> &matSolDOF,
                              ColMajorMatrix<su2double> &matSolInt) override;

  /*!
   * \brief Function, that updates the residuals of the DOFs with the integral of the
   *        product of the given scalar data and the basis function. The integral is
   *        approximated by the weighted sum of the data in the integration points.
   * \param[in]     scalarDataInt - The scalar data in the integration points that must
   *                                be multiplied by the basis functions.
   * \param[in,out] resDOFs       - The residual of the DOFs that must be updated.
   */
  void ResidualBasisFunctions(ColMajorMatrix<su2double> &scalarDataInt,
                              ColMajorMatrix<su2double> &resDOFs) override;

  /*!
   * \brief Function, that updates the residuals of the DOFs with the integral of the
   *        dot product of the given vector data and the gradient of the basis function.
   *        The integral is approximated by the weighted sum of the data in the integration points.
   * \brief Virtual function, that, if used, must be overwritten by the derived class.
   * \param[in]     vectorDataInt - The vector data in the integration points that must
   *                                be multiplied by the gradient of the basis functions.
   * \param[in,out] resDOFs       - The residual of the DOFs that must be updated.
   */
  void ResidualGradientBasisFunctions(vector<ColMajorMatrix<su2double> > &vectorDataInt,
                                      ColMajorMatrix<su2double>          &resDOFs) override;

  /*!
   * \brief Function that makes available the value of the first (constant)
   *        basis function of this element.
   * \return - The value of the first (constant) basis function.
   */
  passivedouble ValBasis0(void) override;

private:
  passivedouble factInviscidRad = -1.0;  /*!< \brief Correction factor for the inviscid spectral radius
                                                     for the high order element. */
  passivedouble factViscousRad = -1.0;   /*!< \brief Correction factor for the viscous spectral radius
                                                     for the high order element. */

  vector<passivedouble> rLineSolDOFs; /*!< \brief 1D parametric coordinates of the nodal solution DOFs. */

  ColMajorMatrix<passivedouble> legBasisLineInt;    /*!< \brief The values of the 1D Legendre basis functions
                                                                in the integration points. */
  ColMajorMatrix<passivedouble> derLegBasisLineInt; /*!< \brief The values of the derivatives of the 1D Legendre
                                                                basis functions in the integration points. */
  ColMajorMatrix<passivedouble> hesLegBasisLineInt; /*!< \brief The values of the 2nd derivatives of the 1D Legendre
                                                                basis functions in the integration points. */

  ColMajorMatrix<passivedouble> legBasisLineIntTranspose;    /*!< \brief Transpose of legBasisLineInt. */
  ColMajorMatrix<passivedouble> derLegBasisLineIntTranspose; /*!< \brief Transpose of derLegBasisLineInt. */

  ColMajorMatrix<passivedouble> legBasisLineSolDOFs;    /*!< \brief The values of the 1D Legendre basis functions
                                                                    in the solution DOFs. */
  ColMajorMatrix<passivedouble> legBasisLineSolDOFsInv; /*!< \brief The inverse of legBasisLineSolDOFs. */

  ColMajorMatrix<passivedouble> derLegBasisLineSolDOFs; /*!< \brief The values of the derivatives of the 1D Legendre
                                                                    basis functions in the solution DOFs. */

  TPI2D TensorProductDataVolIntPoints = nullptr; /*!< \brief Function pointer to carry out the tensor product
                                                             to compute the data in the volume integration points. */
  TPI2D TensorProductDataVolSolDOFs   = nullptr; /*!< \brief Function pointer to carry out the tensor product
                                                             to compute the data in the volume nodal solution DOFs. */
  TPI2D TensorProductResVolDOFs       = nullptr; /*!< \brief Function pointer to carry out the tensor product
                                                             to update the residual of the DOFs for the volume term. */
};
