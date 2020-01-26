/*!
 * \file CCentSca_Heat.hpp
 * \brief Delaration of numerics class CCentSca_Heat, the
 *        implementation is in the CCentSca_Heat.cpp file.
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
 * \class CCentSca_Heat
 * \brief Class for scalar centered scheme.
 * \ingroup ConvDiscr
 * \author O. Burghardt
 * \version 7.0.0 "Blackbird"
 */
class CCentSca_Heat : public CNumerics {

private:
  unsigned short iDim;             /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_Lapl,            /*!< \brief Diference of conservative variables and undivided laplacians. */
  *MeanVelocity, ProjVelocity,
  ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
  Param_Kappa_4,                   /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j,
  MeanLambda,                      /*!< \brief Local eingenvalues. */
  cte_0, cte_1;                    /*!< \brief Artificial dissipation values. */
  bool implicit,                   /*!< \brief Implicit calculation. */
  dynamic_grid;                    /*!< \brief Modification for grid movement. */


public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentSca_Heat(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CCentSca_Heat(void);

  /*!
   * \brief Compute the flow residual using a JST method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};
