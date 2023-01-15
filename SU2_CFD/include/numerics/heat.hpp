/*!
 * \file heat.hpp
 * \brief Declarations of numerics classes for heat transfer problems.
 * \author F. Palacios, T. Economon
 * \version 7.5.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

/*!
 * \class CUpwSca_Heat
 * \brief Class for doing a scalar upwind solver for the heat convection equation.
 * \ingroup ConvDiscr
 * \author O. Burghardt.
 * \version 7.5.0 "Blackbird"
 */
class CUpwSca_Heat : public CNumerics {
private:
  bool implicit, dynamic_grid;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_Heat(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config);

  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) override;
};

/*!
 * \class CAvgGrad_Heat
 * \brief Class for computing viscous term using average of gradients without correction (heat equation).
 * \ingroup ViscDiscr
 * \author O. Burghardt.
 * \version 7.5.0 "Blackbird"
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
