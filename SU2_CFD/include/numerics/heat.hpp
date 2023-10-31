/*!
 * \file heat.hpp
 * \brief Declarations of numerics classes for heat transfer problems.
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

#include "scalar/scalar_diffusion.hpp"
#include "scalar/scalar_convection.hpp"
#include "../variables/CIncEulerVariable.hpp"

/*!
 * \class CUpwSca_Heat
 * \brief Class for doing a scalar upwind solver for the heat convection equation.
 * \ingroup ConvDiscr
 * \author O. Burghardt.
 * \version 8.0.0 "Harrier"
 */
class CUpwSca_Heat final : public CUpwScalar<typename CIncEulerVariable::template CIndices<unsigned short>> {
 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_Heat(unsigned short val_nDim, const CConfig *config)
    : CUpwScalar<typename CIncEulerVariable::template CIndices<unsigned short>>(val_nDim, 1, config) {}

 private:
  /*!
   * \brief Adds extra variables to AD
   */
  void ExtraADPreaccIn(void) override {}

  /*!
   * \brief Heat-specific specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override {
    Flux[0] = a0 * ScalarVar_i[0] + a1 * ScalarVar_j[0];
    Jacobian_i[0][0] = a0;
    Jacobian_j[0][0] = a1;
  }
};

/*!
 * \class CAvgGrad_Heat
 * \brief Class for computing viscous term using average of gradients without correction (heat equation).
 * \ingroup ViscDiscr
 * \author O. Burghardt.
 * \version 8.0.0 "Harrier"
 */
class CAvgGrad_Heat final : public CAvgGrad_Scalar<CNoFlowIndices> {
 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] correct - Whether to correct the gradient.
   */
  CAvgGrad_Heat(unsigned short val_nDim, const CConfig *config, bool correct)
    : CAvgGrad_Scalar<CNoFlowIndices>(val_nDim, 1, correct, config) {}

 private:
  /*!
   * \brief Adds extra variables to AD
   */
  void ExtraADPreaccIn(void) override {
    AD::SetPreaccIn(*Diffusion_Coeff_i, *Diffusion_Coeff_j);
  }

  /*!
   * \brief Heat-specific specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override {
    const su2double Thermal_Diffusivity_Mean = 0.5 * (*Diffusion_Coeff_i + *Diffusion_Coeff_j);

    Flux[0] = Thermal_Diffusivity_Mean * Proj_Mean_GradScalarVar[0];

    /*--- Use TSL for Jacobians. ---*/
    Jacobian_i[0][0] = -Thermal_Diffusivity_Mean * proj_vector_ij;
    Jacobian_j[0][0] = Thermal_Diffusivity_Mean * proj_vector_ij;
  }
};
