/*!
 * \file radiation.hpp
 * \brief Declaration and inlines of the classes used to compute
 *        residual terms in radiation problems.
 * \author Ruben Sanchez
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "CNumerics.hpp"

class CNumericsRadiation : public CNumerics {
protected:
  bool implicit, incompressible;
  const su2double *RadVar_i;              /*!< \brief Vector of radiation variables at point i. */
  const su2double *RadVar_j;              /*!< \brief Vector of radiation variables at point j. */
  CMatrixView<const su2double> RadVar_Grad_i;  /*!< \brief Gradient of turbulent variables at point i. */
  CMatrixView<const su2double> RadVar_Grad_j;  /*!< \brief Gradient of turbulent variables at point j. */
  su2double Absorption_Coeff;             /*!< \brief Absorption coefficient. */
  su2double Scattering_Coeff;             /*!< \brief Scattering coefficient. */

  su2double Temperature_Ref;  /*!< \brief Reference temperature for redimensionalization of P1 solver. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CNumericsRadiation(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Set the value of the radiation variable.
   * \param[in] val_radvar_i - Value of the turbulent variable at point i.
   * \param[in] val_radvar_j - Value of the turbulent variable at point j.
   */
  inline void SetRadVar(const su2double *val_radvar_i, const su2double *val_radvar_j) final {
    RadVar_i = val_radvar_i;
    RadVar_j = val_radvar_j;
  }

  /*!
   * \brief Set the gradient of the radiation variables.
   * \param[in] val_radvar_grad_i - Gradient of the turbulent variable at point i.
   * \param[in] val_radvar_grad_j - Gradient of the turbulent variable at point j.
   */
  inline void SetRadVarGradient(CMatrixView<const su2double> val_radvar_grad_i,
                                CMatrixView<const su2double> val_radvar_grad_j) final {
    RadVar_Grad_i = val_radvar_grad_i;
    RadVar_Grad_j = val_radvar_grad_j;
  }

};


class CSourceP1 final : public CNumericsRadiation {
private:
  su2double BlackBody_Intensity;
  su2double Energy_i, Temperature_i;

public:
  /*!
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceP1(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Source term integration of the P1 model.
   * \param[out] val_residual - Pointer to the residual vector.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) override;

};

class CAvgGradCorrected_P1 final : public CNumericsRadiation {
 private:
  su2double GammaP1;  /*!< \brief P1 parameter */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_P1(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Compute the viscous residual of the P1 equation.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i,
                       su2double **Jacobian_j, CConfig *config) override;

};
