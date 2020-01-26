/*!
 * \file CAvgGradInc_Flow.hpp
 * \brief Delaration of numerics class CAvgGradInc_Flow, the
 *        implementation is in the CAvgGradInc_Flow.cpp file.
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

#include "CAvgGrad_Base.hpp"

/*!
 * \class CAvgGradInc_Flow
 * \brief Class for computing viscous term using an average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, F. Palacios, T. Economon
 */
class CAvgGradInc_Flow : public CAvgGrad_Base {
private:
  su2double Mean_Thermal_Conductivity; /*!< \brief Mean value of the effective thermal conductivity. */
  bool energy;                         /*!< \brief computation with the energy equation. */

  /*
   * \brief Compute the projection of the viscous fluxes into a direction
   *
   * The viscous + turbulent stress tensor must be calculated before calling
   * this function.
   *
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] val_thermal_conductivity - Thermal conductivity.
   */
  void GetViscousIncProjFlux(const su2double* const *val_gradprimvar,
                             const su2double *val_normal,
                             su2double val_thermal_conductivity);

  /*!
   * \brief Compute the projection of the viscous Jacobian matrices.
   *
   * The Jacobian of the stress tensor must be calculated before calling
   * this function.
   *
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousIncProjJacs(su2double val_dS,
                             su2double **val_Proj_Jac_Tensor_i,
                             su2double **val_Proj_Jac_Tensor_j);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_correct_grad - Apply a correction to the gradient
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradInc_Flow(unsigned short val_nDim, unsigned short val_nVar,
                   bool val_correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradInc_Flow(void);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};
