/*!
 * \file CSourceViscous_AdjFlow.hpp
 * \brief Delaration of numerics class CSourceViscous_AdjFlow, the
 *        implementation is in the CSourceViscous_AdjFlow.cpp file.
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

#include "../CNumerics.hpp"

/*!
 * \class CSourceViscous_AdjFlow
 * \brief Class for source term integration in adjoint problem.
 * \ingroup SourceDiscr
 * \author F. Palacios
 */
class CSourceViscous_AdjFlow : public CNumerics {
private:
  su2double *Velocity, *GradDensity, *GradInvDensity, *dPoDensity2, *alpha, *beta, *Sigma_5_vec;
  su2double **GradVel_o_Rho, **sigma, **Sigma_phi, **Sigma_5_Tensor, **Sigma;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceViscous_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourceViscous_AdjFlow(void);

  /*!
   * \brief Source term integration of the flow adjoint equation.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_residual, CConfig *config);

};
