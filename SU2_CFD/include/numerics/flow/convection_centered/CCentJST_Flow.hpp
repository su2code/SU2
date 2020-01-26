/*!
 * \file CCentJST_Flow.hpp
 * \brief Delaration of numerics class CCentJST_Flow, the
 *        implementation is in the CCentJST_Flow.cpp file.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "CCentBase_Flow.hpp"

/*!
 * \class CCentJST_Flow
 * \brief Class for centered scheme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios
 */
class CCentJST_Flow : public CCentBase_Flow {
  
private:
  su2double Param_Kappa_2, Param_Kappa_4; /*!< \brief Artificial dissipation parameters. */
  su2double sc2, sc4;                     /*!< \brief Streching parameters. */
  su2double Epsilon_2, Epsilon_4;         /*!< \brief Artificial dissipation coefficients. */

  /*!
   * \brief JST second and forth order dissipation terms.
   * \param[in,out] val_residual - Pointer to the convective flux contribution to the residual.
   * \param[in,out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in,out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   */
  void DissipationTerm(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j);

  /*!
   * \brief Set input variables for AD preaccumulation.
   * \return true, as we will define inputs.
   */
  bool SetPreaccInVars(void);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CCentJST_Flow(void);

};
