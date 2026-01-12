/*!
 * \file turb_diffusion.hpp
 * \brief Declarations of numerics classes for discretization of
 *        viscous fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
  const su2double cb2 = 0.622;

  const bool use_accurate_jacobians;

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

    /*--- First Term. Normal diffusion, and conservative part of the quadratic diffusion.
     * ||grad nu_t||^2 = div(nu_t grad nu_t) - nu_t div grad nu_t ---*/
    const su2double nu_i = Laminar_Viscosity_i/Density_i;
    const su2double nu_j = Laminar_Viscosity_j/Density_j;
    const su2double nu_e = 0.5 * (nu_i + nu_j + (1 + cb2) * (ScalarVar_i[0] + ScalarVar_j[0]));
    const su2double term_1 = nu_e;

    /* Second Term (quadratic diffusion, non conservative). */
    const su2double nu_tilde_i = ScalarVar_i[0];
    const su2double term_2 = cb2 * nu_tilde_i;

    const su2double diffusion_coefficient = term_1 - term_2;
    Flux[0] = diffusion_coefficient * Proj_Mean_GradScalarVar[0] / sigma;

    if (implicit) {
      /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
      Jacobian_i[0][0] = -diffusion_coefficient * proj_vector_ij / sigma;
      Jacobian_j[0][0] = diffusion_coefficient * proj_vector_ij / sigma;

      if (use_accurate_jacobians) {
        /*--- The diffusion coefficient is also a function of nu_t. ---*/
        const su2double dTerm1_dnut_i = (1 + cb2) * 0.5;
        const su2double dTerm1_dnut_j = (1 + cb2) * 0.5;

        const su2double dTerm2_dnut_i = cb2;
        const su2double dTerm2_dnut_j = 0.0;

        const su2double dDC_dnut_i = dTerm1_dnut_i - dTerm2_dnut_i;
        const su2double dDC_dnut_j = dTerm1_dnut_j - dTerm2_dnut_j;

        Jacobian_i[0][0] += dDC_dnut_i * Proj_Mean_GradScalarVar[0] / sigma;
        Jacobian_j[0][0] += dDC_dnut_j * Proj_Mean_GradScalarVar[0] / sigma;
      }
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
    : CAvgGrad_Scalar<FlowIndices>(val_nDim, val_nVar, correct_grad, config),
      use_accurate_jacobians(config->GetUse_Accurate_Turb_Jacobians()) {}
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
  const su2double cb2 = 0.622;

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

    const su2double nu_ij = 0.5 * (nu_i + nu_j);
    const su2double nu_tilde_i = ScalarVar_i[0];
    const su2double nu_tilde_j = ScalarVar_j[0];
    const su2double nu_tilde_ij = 0.5 * (nu_tilde_i + nu_tilde_j);

    /*--- Following Diskin's implementation from 10.2514/1.J064629, they propose a new fn function
     * to be evaluated at the cell to maintain positivity in the diffusion coefficient, which is
     * used in both terms. The new fn term averaged across the face reverts to the original fn
     * function. ---*/

    /*--- Second Term (LHS) ---*/
    const su2double zeta_i = ((1 + cb2) * nu_tilde_ij - cb2 * nu_tilde_i) / nu_ij;
    su2double fn_i = 1.0;
    if (zeta_i < 0.0) {
      fn_i = (cn1 + pow(zeta_i,3)) / (cn1 - pow(zeta_i,3));
    }

    const su2double term_1 = (nu_ij + (1 + cb2) * nu_tilde_ij * fn_i);
    const su2double term_2 = cb2 * nu_tilde_i * fn_i;
    Flux[0] = (term_1 - term_2) * Proj_Mean_GradScalarVar[0] / sigma;

    /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients
    * Exact Jacobians were tested on multiple cases but resulted in divergence of all
    * simulations, hence only frozen diffusion coefficient (approximate) Jacobians are used. ---*/

    if (implicit) {
      const su2double diffusion_coefficient = (term_1 - term_2);

      const su2double dGrad_dnut_i = -proj_vector_ij;
      const su2double dGrad_dnut_j = proj_vector_ij;

      Jacobian_i[0][0] = diffusion_coefficient * dGrad_dnut_i / sigma;
      Jacobian_j[0][0] = diffusion_coefficient * dGrad_dnut_j / sigma;
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
  const bool use_accurate_jacobians;

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
    const su2double diff_omega_T1 = 0.5*(diff_i_omega + diff_j_omega);

    /*--- We aim to treat the cross-diffusion as a diffusion term rather than a source term.
    * Re-writing the cross-diffusion contribution as λ/w ∇w ∇k, where λ = (2 (1- F1) ρ σ_ω2)
    * and expanding using the product rule for divergence theorem gives: ∇(w λ/w ∇k) - w ∇(λ/w ∇k). 
    * Discretising using FVM, gives: (λ)_ij ∇k - w_c (λ/w)_ij ∇k. where w_c is the cell centre value ---*/

    const su2double lambda_i = 2 * (1 - F1_i) * Density_i * sigma_omega_i;
    const su2double lambda_j = 2 * (1 - F1_j) * Density_j * sigma_omega_j;
    const su2double lambda_ij = 0.5 * (lambda_i + lambda_j);
    const su2double w_ij = 0.5 * (ScalarVar_i[1] + ScalarVar_j[1]);

    const su2double diff_omega_T2 = lambda_ij;
    
    const su2double diff_omega_T3 = -ScalarVar_i[1] * lambda_ij/w_ij;

    Flux[0] = diff_kine*Proj_Mean_GradScalarVar[0];
    Flux[1] = diff_omega_T1*Proj_Mean_GradScalarVar[1] + (diff_omega_T2 + diff_omega_T3)*Proj_Mean_GradScalarVar[0];

    /*--- For Jacobians -> Use of TSL (Thin Shear Layer) approx. to compute derivatives of the gradients ---*/
    if (implicit) {
      const su2double proj_on_rho_i = proj_vector_ij/Density_i;
      const su2double proj_on_rho_j = proj_vector_ij/Density_j;
      Jacobian_i[0][0] = -diff_kine*proj_on_rho_i;
      Jacobian_i[0][1] = 0.0;
      Jacobian_i[1][0] = (diff_omega_T2+diff_omega_T3)*-proj_on_rho_i;
      Jacobian_i[1][1] = -diff_omega_T1*proj_on_rho_i;

      Jacobian_j[0][0] = diff_kine*proj_on_rho_j;
      Jacobian_j[0][1] = 0.0;
      Jacobian_j[1][0] = (diff_omega_T2+diff_omega_T3)*proj_on_rho_j;
      Jacobian_j[1][1] = diff_omega_T1*proj_on_rho_j;

      if (use_accurate_jacobians) {
        Jacobian_i[0][0] = -diff_kine*proj_on_rho_i;
        Jacobian_i[0][1] = 0.0;
        Jacobian_i[1][0] = (diff_omega_T2 + diff_omega_T3)*-proj_on_rho_i;
        Jacobian_i[1][1] = -proj_on_rho_i * diff_omega_T1 - 2*lambda_ij*ScalarVar_j[1]/pow(ScalarVar_i[1]+ScalarVar_j[1],2) * Proj_Mean_GradScalarVar[0];

        Jacobian_j[0][0] = diff_kine*proj_on_rho_j;
        Jacobian_j[0][1] = 0.0;
        Jacobian_j[1][0] = (diff_omega_T2 + diff_omega_T3)*proj_on_rho_j;
        Jacobian_j[1][1] = proj_on_rho_j * diff_omega_T1 + 2*lambda_ij*ScalarVar_i[1]/pow(ScalarVar_i[1]+ScalarVar_j[1],2) * Proj_Mean_GradScalarVar[0];
      }      
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
      sigma_om2(constants[3]),
      use_accurate_jacobians(config->GetUse_Accurate_Turb_Jacobians()) {
  }

  /*!
   * \brief Sets value of first blending function.
   */
  void SetF1blending(su2double val_F1_i, su2double val_F1_j) override {
    F1_i = val_F1_i; F1_j = val_F1_j;
  }
};
