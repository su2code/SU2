/*!
 * \file CGeneralAvgGrad_Flow.hpp
 * \brief Delaration of numerics class CGeneralAvgGrad_Flow, the
 *        implementation is in the CGeneralAvgGrad_Flow.cpp file.
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
 * \class CGeneralAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author M.Pini, S. Vitale
 */
class CGeneralAvgGrad_Flow : public CAvgGrad_Base {
private:
  su2double *Mean_SecVar,    /*!< \brief Mean secondary variables. */
  Mean_Thermal_Conductivity, /*!< \brief Mean value of the thermal conductivity. */
  Mean_Cp;                   /*!< \brief Mean value of the Cp. */

  /*!
   * \brief Compute the heat flux due to molecular and turbulent diffusivity
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   * \param[in] val_thermal_conductivity - Thermal Conductivity.
   * \param[in] val_heat_capacity_cp - Heat Capacity at constant pressure.
   */
  void SetHeatFluxVector(const su2double* const *val_gradprimvar,
                         su2double val_laminar_viscosity,
                         su2double val_eddy_viscosity,
                         su2double val_thermal_conductivity,
                         su2double val_heat_capacity_cp);

  /*!
   * \brief Compute the Jacobian of the heat flux vector
   *
   * This Jacobian is projected onto the normal vector, so it is of
   * dimension nVar.
   *
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_Mean_SecVar - Mean value of the secondary variables.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_thermal_conductivity - Value of the thermal conductivity.
   * \param[in] val_heat_capacity_cp - Value of the specific heat at constant pressure.
   * \param[in] val_dist_ij - Distance between the points.
   */
  void SetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                           const su2double *val_Mean_SecVar,
                           su2double val_eddy_viscosity,
                           su2double val_thermal_conductivity,
                           su2double val_heat_capacity_cp,
                           su2double val_dist_ij);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_correct_grad - Apply a correction to the gradient
   * \param[in] config - Definition of the particular problem.
   */
  CGeneralAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, bool val_correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CGeneralAvgGrad_Flow(void);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};
