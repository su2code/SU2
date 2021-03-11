/*!
 * \file scalar_sources.hpp
 * \brief This file contains the numerical methods for scalar transport eqns.
 * \author T. Economon, D. Mayer, N. Beishuizen
 * \version 7.1.0 "Blackbird"
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

#include "../CNumerics.hpp"

/*!
 * \class CSourcePieceWise_Scalar
 * \brief Class for integrating the source terms of scalar transport equations.
 * \ingroup SourceDiscr
 * \author T. Economon
 */
class CSourcePieceWise_Scalar : public CNumerics {
private:
  bool incompressible;  /*!< \brief Flag defining compressibility. */
  bool implicit;        /*!< \brief Flag defining implicit scheme. */

public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_Scalar(unsigned short val_nDim,
                          unsigned short val_nVar,
                          CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_Scalar(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual,
                       su2double **val_Jacobian_i,
                       su2double **val_Jacobian_j,
                       CConfig *config);

};

/*!
 * \class CSourceAxisymmetric_Scalar
 * \brief Class for source term for solving scalar axisymmetric problems.
 * \ingroup SourceDiscr
 * \author T. Economon
 */
class CSourceAxisymmetric_Scalar : public CNumerics {
  bool implicit, /*!< \brief Implicit calculation. */
  viscous, /*!< \brief Viscous incompressible flows. */
  energy; /*!< \brief computation with the energy equation. */

  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceAxisymmetric_Scalar(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceAxisymmetric_Scalar(void);
  
  /*!
   * \brief Residual of the rotational frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config);
  
};