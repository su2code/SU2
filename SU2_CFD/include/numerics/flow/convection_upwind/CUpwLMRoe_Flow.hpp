/*!
 * \file CUpwLMRoe_Flow.hpp
 * \brief Delaration of numerics class CUpwLMRoe_Flow, the
 *        implementation is in the CUpwLMRoe_Flow.cpp file.
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
 * \class CUpwLMRoe_Flow
 * \brief Class for solving an approximate Riemann solver of LMRoe for the flow equations.
 * \ingroup ConvDiscr
 * \author E. Molina, A. Bueno, F. Palacios
 * \version 7.0.0 "Blackbird"
 */
class CUpwLMRoe_Flow : public CUpwRoeBase_Flow {
private:
  /*!
   * \brief Add LMRoe dissipation to the flux (low-Mach scheme).
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
   */
  CUpwLMRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwLMRoe_Flow(void);
};
