/*!
 * \file trans_diffusion.hpp
 * \brief Declarations of numerics classes for discretization of
 *        viscous fluxes in transition problems.
 * \author S. Kang
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


#include "../../scalar/scalar_diffusion.hpp"

/*!
 * \class CAvgGrad_TransLM
 * \brief Class for computing viscous term using average of gradient with correction (LM transition model).
 * \ingroup ViscDiscr
 * \author S. Kang.
 */
template <class FlowIndices>
class CAvgGrad_TransLM final : public CAvgGrad_Scalar<FlowIndices> {
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

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override {}

  /*!
   * \brief LM transition model specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override {
    const bool implicit = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;

    /*--- Compute mean effective dynamic viscosity ---*/
    const su2double diff_i_gamma = Laminar_Viscosity_i + Eddy_Viscosity_i;
    const su2double diff_j_gamma = Laminar_Viscosity_j + Eddy_Viscosity_j;
    const su2double diff_i_ReThetaT = 2.0*(Laminar_Viscosity_i + Eddy_Viscosity_i);
    const su2double diff_j_ReThetaT = 2.0*(Laminar_Viscosity_j + Eddy_Viscosity_j);

    const su2double diff_gamma = 0.5*(diff_i_gamma + diff_j_gamma);
    const su2double diff_ReThetaT = 0.5*(diff_i_ReThetaT + diff_j_ReThetaT);

    Flux[0] = diff_gamma*Proj_Mean_GradScalarVar[0];
    Flux[1] = diff_ReThetaT*Proj_Mean_GradScalarVar[1];

    /*--- For Jacobians -> Use of TSL (Thin Shear Layer) approx. to compute derivatives of the gradients ---*/
    if (implicit) {
      const su2double proj_on_rho_i = proj_vector_ij/Density_i;
      Jacobian_i[0][0] = -diff_gamma*proj_on_rho_i;  Jacobian_i[0][1] = 0.0;
      Jacobian_i[1][0] = 0.0;                        Jacobian_i[1][1] = -diff_ReThetaT*proj_on_rho_i;

      const su2double proj_on_rho_j = proj_vector_ij/Density_j;
      Jacobian_j[0][0] = diff_gamma*proj_on_rho_j;   Jacobian_j[0][1] = 0.0;
      Jacobian_j[1][0] = 0.0;                        Jacobian_j[1][1] = diff_ReThetaT*proj_on_rho_j;
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
  CAvgGrad_TransLM(unsigned short val_nDim, unsigned short val_nVar, bool correct_grad, const CConfig* config)
    : CAvgGrad_Scalar<FlowIndices>(val_nDim, val_nVar, correct_grad, config){
  }

};
