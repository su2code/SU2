/*!
 * \file turb_diffusion.hpp
 * \brief Declarations of numerics classes for discretization of
 *        viscous fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
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

#include "../scalar/scalar_diffusion.hpp"

/*!
 * \class CAvgGrad_TurbSA
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
template <class FlowIndices>
class CAvgGrad_TurbSA final : public CAvgGrad_Scalar<FlowIndices> {
private:
  using Base = CAvgGrad_Scalar<FlowIndices>;
  using Base::Laminar_Viscosity_i;
  using Base::Laminar_Viscosity_j;
  using Base::Density_i;
  using Base::Density_j;
  using Base::ScalarVar_i;
  using Base::ScalarVar_j;
  using Base::Proj_Mean_GradScalarVar;
  using Base::proj_vector_ij;
  using Base::Flux;
  using Base::Jacobian_i;
  using Base::Jacobian_j;

  const su2double sigma = 2.0/3.0;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override {}

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override {
    const bool implicit = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;

    /*--- Compute mean effective viscosity ---*/

    const su2double nu_i = Laminar_Viscosity_i/Density_i;
    const su2double nu_j = Laminar_Viscosity_j/Density_j;
    const su2double nu_e = 0.5*(nu_i+nu_j+ScalarVar_i[0]+ScalarVar_j[0]);

    Flux[0] = nu_e*Proj_Mean_GradScalarVar[0]/sigma;

    /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

    if (implicit) {
      Jacobian_i[0][0] = (0.5*Proj_Mean_GradScalarVar[0]-nu_e*proj_vector_ij)/sigma;
      Jacobian_j[0][0] = (0.5*Proj_Mean_GradScalarVar[0]+nu_e*proj_vector_ij)/sigma;
    }
  }

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_grad - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                  bool correct_grad, const CConfig* config)
    : CAvgGrad_Scalar<FlowIndices>(val_nDim, val_nVar, correct_grad, config) {}
};

/*!
 * \class CAvgGrad_TurbSA_Neg
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author F. Palacios
 */
template <class FlowIndices>
class CAvgGrad_TurbSA_Neg final : public CAvgGrad_Scalar<FlowIndices> {
private:
  using Base = CAvgGrad_Scalar<FlowIndices>;
  using Base::Laminar_Viscosity_i;
  using Base::Laminar_Viscosity_j;
  using Base::Density_i;
  using Base::Density_j;
  using Base::ScalarVar_i;
  using Base::ScalarVar_j;
  using Base::Proj_Mean_GradScalarVar;
  using Base::proj_vector_ij;
  using Base::Flux;
  using Base::Jacobian_i;
  using Base::Jacobian_j;

  const su2double sigma = 2.0/3.0;
  const su2double cn1 = 16.0;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override {}

  /*!
   * \brief SA-neg specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override {
    const bool implicit = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;

    /*--- Compute mean effective viscosity ---*/

    const su2double nu_i = Laminar_Viscosity_i/Density_i;
    const su2double nu_j = Laminar_Viscosity_j/Density_j;

    const su2double nu_ij = 0.5*(nu_i+nu_j);
    const su2double nu_tilde_ij = 0.5*(ScalarVar_i[0] + ScalarVar_j[0]);

    su2double nu_e;

    if (nu_tilde_ij > 0.0) {
      nu_e = nu_ij + nu_tilde_ij;
    }
    else {
      const su2double Xi = nu_tilde_ij/nu_ij;
      const su2double fn = (cn1 + Xi*Xi*Xi)/(cn1 - Xi*Xi*Xi);
      nu_e = nu_ij + fn*nu_tilde_ij;
    }

    Flux[0] = nu_e*Proj_Mean_GradScalarVar[0]/sigma;

    /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

    if (implicit) {
      Jacobian_i[0][0] = (0.5*Proj_Mean_GradScalarVar[0]-nu_e*proj_vector_ij)/sigma;
      Jacobian_j[0][0] = (0.5*Proj_Mean_GradScalarVar[0]+nu_e*proj_vector_ij)/sigma;
    }
  }

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_grad - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar,
                      bool correct_grad, const CConfig* config)
    : CAvgGrad_Scalar<FlowIndices>(val_nDim, val_nVar, correct_grad, config) {}
};

/*!
 * \class CAvgGrad_TurbSST
 * \brief Class for computing viscous term using average of gradient with correction (Menter SST turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
template <class FlowIndices>
class CAvgGrad_TurbSST final : public CAvgGrad_Scalar<FlowIndices> {
private:
  using Base = CAvgGrad_Scalar<FlowIndices>;
  using Base::Laminar_Viscosity_i;
  using Base::Laminar_Viscosity_j;
  using Base::Eddy_Viscosity_i;
  using Base::Eddy_Viscosity_j;
  using Base::Density_i;
  using Base::Density_j;
  using Base::ScalarVar_i;
  using Base::ScalarVar_j;
  using Base::Proj_Mean_GradScalarVar;
  using Base::proj_vector_ij;
  using Base::Flux;
  using Base::Jacobian_i;
  using Base::Jacobian_j;

  const su2double sigma_k1; /*!< \brief Constants for the viscous terms, k-w (1), k-eps (2)*/
  const su2double sigma_k2;
  const su2double sigma_om1;
  const su2double sigma_om2;

  su2double F1_i, F1_j; /*!< \brief Menter's first blending function */

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override {
    AD::SetPreaccIn(F1_i, F1_j);
  }

  /*!
   * \brief SST specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override {
    const bool implicit = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;

    /*--- Compute the blended constant for the viscous terms ---*/
    const su2double sigma_kine_i = F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2;
    const su2double sigma_kine_j = F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2;
    const su2double sigma_omega_i = F1_i*sigma_om1 + (1.0 - F1_i)*sigma_om2;
    const su2double sigma_omega_j = F1_j*sigma_om1 + (1.0 - F1_j)*sigma_om2;

    /*--- Compute mean effective dynamic viscosity ---*/
    const su2double diff_i_kine = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
    const su2double diff_j_kine = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
    const su2double diff_i_omega = Laminar_Viscosity_i + sigma_omega_i*Eddy_Viscosity_i;
    const su2double diff_j_omega = Laminar_Viscosity_j + sigma_omega_j*Eddy_Viscosity_j;

    const su2double diff_kine = 0.5*(diff_i_kine + diff_j_kine);
    const su2double diff_omega = 0.5*(diff_i_omega + diff_j_omega);

    Flux[0] = diff_kine*Proj_Mean_GradScalarVar[0];
    Flux[1] = diff_omega*Proj_Mean_GradScalarVar[1];

    /*--- For Jacobians -> Use of TSL (Thin Shear Layer) approx. to compute derivatives of the gradients ---*/
    if (implicit) {
      const su2double proj_on_rho_i = proj_vector_ij/Density_i;
      Jacobian_i[0][0] = -diff_kine*proj_on_rho_i;  Jacobian_i[0][1] = 0.0;
      Jacobian_i[1][0] = 0.0;                       Jacobian_i[1][1] = -diff_omega*proj_on_rho_i;

      const su2double proj_on_rho_j = proj_vector_ij/Density_j;
      Jacobian_j[0][0] = diff_kine*proj_on_rho_j;   Jacobian_j[0][1] = 0.0;
      Jacobian_j[1][0] = 0.0;                       Jacobian_j[1][1] = diff_omega*proj_on_rho_j;
    }
  }

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] constants - Constants of the model.
   * \param[in] correct_grad - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
                   const su2double* constants, bool correct_grad, const CConfig* config)
    : CAvgGrad_Scalar<FlowIndices>(val_nDim, val_nVar, correct_grad, config),
      sigma_k1(constants[0]),
      sigma_k2(constants[1]),
      sigma_om1(constants[2]),
      sigma_om2(constants[3]) {
  }

  /*!
   * \brief Sets value of first blending function.
   */
  void SetF1blending(su2double val_F1_i, su2double val_F1_j) override {
    F1_i = val_F1_i; F1_j = val_F1_j;
  }
};
