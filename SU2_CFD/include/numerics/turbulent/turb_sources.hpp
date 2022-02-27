/*!
 * \file turb_sources.hpp
 * \brief Delarations of numerics classes for integration of source
 *        terms in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.3.0 "Blackbird"
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

#include "../scalar/scalar_sources.hpp"
#include "turb_sources_new.hpp"

/*!
 * \class CSourcePieceWise_TurbSST
 * \brief Class for integrating the source terms of the Menter SST turbulence model equations.
 * \ingroup SourceDiscr
 * \author A. Campos.
 */
template <class FlowIndices>
class CSourcePieceWise_TurbSST final : public CNumerics {
private:
  const FlowIndices idx;  /*!< \brief Object to manage the access to the flow primitives. */

  su2double F1_i,
  F1_j,
  F2_i,
  F2_j;

  su2double alfa_1,
  alfa_2,
  beta_1,
  beta_2,
  sigma_k_1,
  sigma_k_2,
  sigma_w_1,
  sigma_w_2,
  beta_star,
  a1;

  su2double CDkw_i, CDkw_j;

  su2double kAmb, omegaAmb;

  su2double Residual[2];
  su2double* Jacobian_i[2];
  su2double Jacobian_Buffer[4]; /// Static storage for the Jacobian (which needs to be pointer for return type).

  bool sustaining_terms;
  bool axisymmetric;

  /*!
   * \brief A virtual member. Get strain magnitude based on perturbed reynolds stress matrix
   * \param[in] turb_ke: turbulent kinetic energy of the node
   */
  void SetPerturbedStrainMag(su2double turb_ke);

  /*!
   * \brief Add contribution due to axisymmetric formulation to 2D residual
   */
  inline void ResidualAxisymmetric(su2double alfa_blended, su2double zeta) {

    if (Coord_i[1] < EPS) return;

    su2double yinv, rhov, k, w;
    su2double sigma_k_i, sigma_w_i;
    su2double pk_axi, pw_axi, cdk_axi, cdw_axi;

    AD::SetPreaccIn(Coord_i[1]);

    yinv = 1.0/Coord_i[1];
    rhov = Density_i*V_i[2];
    k = ScalarVar_i[0];
    w = ScalarVar_i[1];

    /*--- Compute blended constants ---*/
    sigma_k_i = F1_i*sigma_k_1+(1.0-F1_i)*sigma_k_2;
    sigma_w_i = F1_i*sigma_w_1+(1.0-F1_i)*sigma_w_2;

    /*--- Production ---*/
    pk_axi = max(0.0,2.0/3.0*rhov*k*((2.0*yinv*V_i[2]-PrimVar_Grad_i[2][1]-PrimVar_Grad_i[1][0])/zeta-1.0));
    pw_axi = alfa_blended*zeta/k*pk_axi;

    /*--- Convection-Diffusion ---*/
    cdk_axi = rhov*k-(Laminar_Viscosity_i+sigma_k_i*Eddy_Viscosity_i)*ScalarVar_Grad_i[0][1];
    cdw_axi = rhov*w-(Laminar_Viscosity_i+sigma_w_i*Eddy_Viscosity_i)*ScalarVar_Grad_i[1][1];

    /*--- Add terms to the residuals ---*/
    Residual[0] += yinv*Volume*(pk_axi-cdk_axi);
    Residual[1] += yinv*Volume*(pw_axi-cdw_axi);

  }

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar, const su2double* constants,
                           su2double val_kine_Inf, su2double val_omega_Inf, const CConfig* config);

  /*!
   * \brief Set the value of the first blending function.
   * \param[in] val_F1_i - Value of the first blending function at point i.
   * \param[in] val_F1_j - Value of the first blending function at point j.
   */
  inline void SetF1blending(su2double val_F1_i, su2double val_F1_j) override {
    F1_i = val_F1_i;
    F1_j = val_F1_j;
  }

  /*!
   * \brief Set the value of the second blending function.
   * \param[in] val_F2_i - Value of the second blending function at point i.
   * \param[in] val_F2_j - Value of the second blending function at point j.
   */
  inline void SetF2blending(su2double val_F2_i, su2double val_F2_j) override {
    F2_i = val_F2_i;
    F2_j = val_F2_j;
  }

  /*!
   * \brief Set the value of the cross diffusion for the SST model.
   * \param[in] val_CDkw_i - Value of the cross diffusion at point i.
   * \param[in] val_CDkw_j - Value of the cross diffusion at point j.
   */
  inline void SetCrossDiff(su2double val_CDkw_i, su2double val_CDkw_j) override {
    CDkw_i = val_CDkw_i;
    CDkw_j = val_CDkw_j;
  }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

};
