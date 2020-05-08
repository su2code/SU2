/*!
 * \file cusp.hpp
 * \brief Declaration of the CUSP numerics class.
 * \author F. Palacios, T. Economon
 * \version 7.0.4 "Blackbird"
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

#include "../../CNumerics.hpp"

/*!
 * \class CUpwCUSP_Flow
 * \brief Class for centered scheme - CUSP.
 * \ingroup ConvDiscr
 * \author F. Palacios
 */
class CUpwCUSP_Flow final : public CNumerics {

private:
  su2double Velocity_i[MAXNDIM], Velocity_j[MAXNDIM], *ProjFlux_i, *ProjFlux_j;
  bool implicit;

  su2double* Flux;        /*!< \brief The flux accross the face. */
  su2double** Jacobian_i; /*!< \brief The Jacobian w.r.t. point i after computation. */
  su2double** Jacobian_j; /*!< \brief The Jacobian w.r.t. point j after computation. */
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwCUSP_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwCUSP_Flow(void) override;

  /*!
   * \brief Compute the flow residual using a JST method.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};
