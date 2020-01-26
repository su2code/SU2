/*!
 * \file CUpwRoe_Flow.hpp
 * \brief Delaration of numerics class CUpwRoe_Flow, the
 *        implementation is in the CUpwRoe_Flow.cpp file.
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

#include "CUpwRoeBase_Flow.hpp"

/*!
 * \class CUpwRoe_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the flow equations.
 * \ingroup ConvDiscr
 * \author A. Bueno, F. Palacios
 */
class CUpwRoe_Flow : public CUpwRoeBase_Flow {
private:
  /*!
   * \brief Add standard Roe dissipation to the flux.
   * \param[out] val_residual - Convective flux.
   * \param[out] val_Jacobian_i - Flux Jacobian wrt node i conservatives (implicit computation).
   * \param[out] val_Jacobian_j - Flux Jacobian wrt node j conservatives (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_low_dissipation - Use a low dissipation formulation.
   */
  CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config, bool val_low_dissipation);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwRoe_Flow(void);
  
};
