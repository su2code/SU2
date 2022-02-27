/*!
 * \file turb_sources.hpp
 * \brief Numerics classes for integration of source terms in turbulence problems.
 * \version 7.3.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../scalar/scalar_sources.hpp"

/*!
 * \class CSAVariables
 * \brief Structure with SA common auxiliary functions and constants.
 */
struct CSAVariables {
  /*--- List of constants ---*/
  const su2double cv1_3 = pow(7.1, 3), k2 = pow(0.41, 2), cb1 = 0.1355, cw2 = 0.3, ct3 = 1.2, ct4 = 0.5,
                  cw3_6 = pow(2, 6), sigma = 2.0 / 3.0, cb2 = 0.622, cb2_sigma = cb2 / sigma,
                  cw1 = cb1 / k2 + (1 + cb2) / sigma, cr1 = 0.5;

  /*--- List of auxiliary functions ---*/
  su2double ft2, d_ft2, r, d_r, g, d_g, glim, fw, d_fw, Ji, d_Ji, S, Shat, d_Shat, fv1, d_fv1, fv2, d_fv2;

  /*--- List of helpers ---*/
  su2double Omega, dist_i_2, inv_k2_d2, inv_Shat, g_6, norm2_Grad, gamma_bc;
};

/*!
 * \class CSourceBase_TurbSA
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * The variables that are subject to change in each variation/correction have their own class.
 * \note Additional source terms (e.g. compressibility) are implemented as decorators.
 */
template <class FlowIndices, class Omega, class ft2, class ModVort, class rFunc, class SourceTerms>
class CSourceBase_TurbSA : public CNumerics {
 protected:
  su2double Gamma_BC = 0.0;

  /*--- Residual and Jacobian ---*/
  su2double Residual, *Jacobian_i;
  su2double Jacobian_Buffer; /*!< \brief Static storage for the Jacobian (which needs to be pointer for return type). */

  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */
  const bool rotating_frame = false;
  const bool transition = false;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceBase_TurbSA(unsigned short nDim, unsigned short, const CConfig* config)
      : CNumerics(nDim, 1, config),
        idx(nDim, config->GetnSpecies()),
        rotating_frame(config->GetRotating_Frame()),
        transition(config->GetKind_Trans_Model() == TURB_TRANS_MODEL::BC) {
    /*--- Setup the Jacobian pointer, we need to return su2double** but we know
     * the Jacobian is 1x1 so we use this trick to avoid heap allocation. ---*/
    Jacobian_i = &Jacobian_Buffer;
  }

  /*!
   * \brief Get the intermittency for the BC transition model.
   * \return Value of the intermittency.
   */
  su2double GetGammaBC() const final { return Gamma_BC; }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override {
    /// TODO: Fix AD preaccumulation.
    //    AD::StartPreacc();
    //    AD::SetPreaccIn(V_i, nDim+6);
    //    AD::SetPreaccIn(Vorticity_i, nDim);
    //    AD::SetPreaccIn(StrainMag_i);
    //    AD::SetPreaccIn(ScalarVar_i[0]);
    //    AD::SetPreaccIn(ScalarVar_Grad_i[0], nDim);
    //    AD::SetPreaccIn(Volume);
    //    AD::SetPreaccIn(dist_i);

    /*--- Common auxiliary variables and constants of the model. ---*/
    CSAVariables var;

    const auto& density = V_i[idx.Density()];
    const auto& laminar_viscosity = V_i[idx.LaminarViscosity()];

    Residual = 0.0;
    Jacobian_i[0] = 0.0;

    /*--- Evaluate Omega with a rotational correction term. ---*/

    Omega::get(Vorticity_i, nDim, PrimVar_Grad_i + idx.Velocity(), var);

    /// TODO: Make this one of the template parameters?
    if (rotating_frame) {
      var.Omega += 2.0 * min(0.0, StrainMag_i - var.Omega);
    }

    if (dist_i > 1e-10) {
      /*--- Vorticity ---*/
      var.S = var.Omega;

      var.dist_i_2 = pow(dist_i, 2);
      const su2double nu = laminar_viscosity / density;
      var.inv_k2_d2 = 1.0 / (var.k2 * var.dist_i_2);

      /*--- Modified values for roughness, roughness_i = 0 for smooth walls and Ji remains the same.
       * Ref: Aupoix, B. and Spalart, P. R., "Extensions of the Spalart-Allmaras Turbulence Model to Account for Wall
       * Roughness," International Journal of Heat and Fluid Flow, Vol. 24, 2003, pp. 454-462.
       * See https://turbmodels.larc.nasa.gov/spalart.html#sarough for detailed explanation. ---*/
      var.Ji = ScalarVar_i[0] / nu + var.cr1 * (roughness_i / (dist_i + EPS));
      var.d_Ji = 1.0 / nu;

      const su2double Ji_2 = pow(var.Ji, 2);
      const su2double Ji_3 = Ji_2 * var.Ji;

      var.fv1 = Ji_3 / (Ji_3 + var.cv1_3);
      var.d_fv1 = 3 * Ji_2 * var.cv1_3 / (nu * pow(Ji_3 + var.cv1_3, 2));

      /*--- Using a modified relation so as to not change the Shat that depends on fv2.
       * From NASA turb modeling resource and 2003 paper. ---*/
      var.fv2 = 1 - ScalarVar_i[0] / (nu + ScalarVar_i[0] * var.fv1);
      var.d_fv2 = -(1 / nu - Ji_2 * var.d_fv1) / pow(1 + var.Ji * var.fv1, 2);

      /*--- Compute ft2 term ---*/
      ft2::get(var);

      /*--- Compute modified vorticity ---*/
      ModVort::get(ScalarVar_i[0], nu, var);
      var.inv_Shat = 1.0 / var.Shat;

      /*--- Compute auxiliary function r ---*/
      rFunc::get(ScalarVar_i[0], var);

      var.g = var.r + var.cw2 * (pow(var.r, 6) - var.r);
      var.g_6 = pow(var.g, 6);
      var.glim = pow((1 + var.cw3_6) / (var.g_6 + var.cw3_6), 1.0 / 6.0);
      var.fw = var.g * var.glim;

      var.d_g = var.d_r * (1 + var.cw2 * (6 * pow(var.r, 5) - 1));
      var.d_fw = var.d_g * var.glim * (1 - var.g_6 / (var.g_6 + var.cw3_6));

      var.norm2_Grad = GeometryToolbox::SquaredNorm(nDim, ScalarVar_Grad_i[0]);

      /// TODO: Make this one of the template parameters?
      if (transition) {
        /*--- BC model constants (2020 revision). ---*/
        const su2double chi_1 = 0.002;
        const su2double chi_2 = 50.0;

        /*--- turbulence intensity is u'/U so we multiply by 100 to get percentage ---*/
        const su2double tu = 100.0 * config->GetTurbulenceIntensity_FreeStream();
        const su2double nu_t = ScalarVar_i[0] * var.fv1;

        const su2double re_v = density * var.dist_i_2 / laminar_viscosity * var.Omega;
        const su2double re_theta = re_v / 2.193;
        const su2double re_theta_t = 803.73 * pow(tu + 0.6067, -1.027);  // MENTER correlation
        // re_theta_t = 163.0 + exp(6.91-tu); //ABU-GHANNAM & SHAW correlation

        const su2double term1 = sqrt(max(re_theta - re_theta_t, 0.0) / (chi_1 * re_theta_t));
        const su2double term2 = sqrt(max((nu_t * chi_2) / nu, 0.0));
        const su2double term_exponential = (term1 + term2);

        Gamma_BC = 1.0 - exp(-term_exponential);
        var.gamma_bc = Gamma_BC;
      } else {
        /*--- Do not modify the production. ---*/
        var.gamma_bc = 1.0;
      }

      /*--- Compute production, destruction and cross production and jacobian ---*/
      su2double Production = 0.0, Destruction = 0.0, CrossProduction = 0.0;
      SourceTerms::get(ScalarVar_i[0], var, Production, Destruction, CrossProduction, Jacobian_i[0]);

      Residual = (Production - Destruction + CrossProduction) * Volume;
      Jacobian_i[0] *= Volume;
    }

    //    AD::SetPreaccOut(Residual);
    //    AD::EndPreacc();

    return ResidualType<>(&Residual, &Jacobian_i, nullptr);
  }
};

namespace detail {

/* =============================================================================
 * SPALART-ALLMARAS VARIATIONS AND CORRECTIONS
 * ============================================================================*/

/*!
 * \brief Strain rate classes.
 * \param[in] vorticity: Vorticity array.
 * \param[in] nDim: Problem dimension.
 * \param[in] velocity_grad: Velocity gradients.
 * \param[out] var: Common SA variables struct (to set Omega).
 */
namespace Omega {

/*! \brief Baseline. */
struct Bsl {
  template <class MatrixType>
  static void get(const su2double* vorticity, unsigned short, const MatrixType&, CSAVariables& var) {
    var.Omega = GeometryToolbox::Norm(3, vorticity);
  }
};

/*! \brief Edward. Here Omega is the Strain Rate. */
struct Edw {
  template <class MatrixType>
  static void get(const su2double*, unsigned short nDim, const MatrixType& velocity_grad, CSAVariables& var) {
    su2double Sbar = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      for (unsigned short jDim = 0; jDim < nDim; ++jDim) {
        Sbar += (velocity_grad[iDim][jDim] + velocity_grad[jDim][iDim]) * velocity_grad[iDim][jDim];
      }
    }
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      Sbar -= (2.0 / 3.0) * pow(velocity_grad[iDim][iDim], 2);
    }
    var.Omega = sqrt(max(Sbar, 0.0));
  }
};
}  // namespace Omega

/*!
 * \brief Classes to set the ft2 term and its derivative.
 * \param[in,out] var: Common SA variables struct.
 */
namespace ft2 {

/*! \brief No-ft2 term. */
struct Zero {
  static void get(CSAVariables& var) {
    var.ft2 = 0.0;
    var.d_ft2 = 0.0;
  }
};

/*! \brief Non-zero ft2 term according to the literature. */
struct Nonzero {
  static void get(CSAVariables& var) {
    const su2double xsi2 = pow(var.Ji, 2);
    var.ft2 = var.ct3 * exp(-var.ct4 * xsi2);
    var.d_ft2 = -2.0 * var.ct4 * var.Ji * var.ft2 * var.d_Ji;
  }
};
}  // namespace ft2

/*!
 * \brief Classes to compute the modified vorticity (\tilde{S}) and its derivative.
 * \param[in] nue: SA variable.
 * \param[in] nu: Laminar viscosity.
 * \param[in,out] var: Common SA variables struct.
 */
namespace ModVort {

/*! \brief Baseline. */
struct Bsl {
  static void get(const su2double& nue, const su2double& nu, CSAVariables& var) {
    const su2double Sbar = nue * var.fv2 * var.inv_k2_d2;
    var.Shat = var.S + Sbar;
    var.Shat = max(var.Shat, 1.0e-10);
    if (var.Shat <= 1.0e-10) {
      var.d_Shat = 0.0;
    } else {
      var.d_Shat = (var.fv2 + nue * var.d_fv2) * var.inv_k2_d2;
    }
  }
};

/*! \brief Edward. */
struct Edw {
  static void get(const su2double& nue, const su2double& nu, CSAVariables& var) {
    var.Shat = max(var.S * ((1.0 / max(var.Ji, 1.0e-16)) + var.fv1), 1.0e-16);
    var.Shat = max(var.Shat, 1.0e-10);
    if (var.Shat <= 1.0e-10) {
      var.d_Shat = 0.0;
    } else {
      var.d_Shat = -var.S * pow(var.Ji, -2) / nu + var.S * var.d_fv1;
    }
  }
};

/*! \brief Negative. */
struct Neg {
  static void get(const su2double& nue, const su2double& nu, CSAVariables& var) {
    if (nue > 0.0) {
      // Baseline solution
      Bsl::get(nue, nu, var);
    } else {
      var.Shat = 1.0e-10;
      var.d_Shat = 0.0;
    }
    /*--- Don't check whether Sbar <>= -cv2*S.
     * Steven R. Allmaras, Forrester T. Johnson and Philippe R. Spalart -
     * "Modifications and Clarifications for the Implementation of the Spalart-Allmaras Turbulence Model" eq. 12
     * No need for Sbar ---*/
  }
};
}  // namespace ModVort

/*!
 * \brief Auxiliary function r and its derivative.
 * \param[in] nue: SA variable.
 * \param[in,out] var: Common SA variables struct.
 */
namespace r {

/*! \brief Baseline. */
struct Bsl {
  static void get(const su2double& nue, CSAVariables& var) {
    var.r = min(nue * var.inv_Shat * var.inv_k2_d2, 10.0);
    var.d_r = (var.Shat - nue * var.d_Shat) * pow(var.inv_Shat, 2) * var.inv_k2_d2;
    if (var.r >= 10.0) var.d_r = 0.0;
  }
};

/*! \brief Edward. */
struct Edw {
  static void get(const su2double& nue, CSAVariables& var) {
    var.r = min(nue * var.inv_Shat * var.inv_k2_d2, 10.0);
    var.r = tanh(var.r) / tanh(1.0);

    var.d_r = (var.Shat - nue * var.d_Shat) * pow(var.inv_Shat, 2) * var.inv_k2_d2;
    var.d_r = (1 - pow(tanh(var.r), 2.0)) * (var.d_r) / tanh(1.0);
  }
};
}  // namespace r

/*!
 * \brief Source terms classes: production, destruction and cross-productions term and their derivative.
 * \param[in] nue: SA variable.
 * \param[in] var: Common SA variables struct.
 * \param[out] production: Production term.
 * \param[out] destruction: Destruction term.
 * \param[out] cross_production: CrossProduction term.
 * \param[out] jacobian: Derivative of the combined source term wrt nue.
 */
namespace SourceTerms {

/*! \brief Baseline (Original SA model). */
struct Bsl {
  static void get(const su2double& nue, const CSAVariables& var, su2double& production, su2double& destruction,
                  su2double& cross_production, su2double& jacobian) {
    ComputeProduction(nue, var, production, jacobian);
    ComputeDestruction(nue, var, destruction, jacobian);
    ComputeCrossProduction(nue, var, cross_production, jacobian);
  }

  static void ComputeProduction(const su2double& nue, const CSAVariables& var, su2double& production,
                                su2double& jacobian) {
    const su2double factor = var.gamma_bc * var.cb1;
    production = factor * (1.0 - var.ft2) * var.Shat * nue;
    jacobian += factor * (-var.Shat * nue * var.d_ft2 + (1.0 - var.ft2) * (nue * var.d_Shat + var.Shat));
  }

  static void ComputeDestruction(const su2double& nue, const CSAVariables& var, su2double& destruction,
                                 su2double& jacobian) {
    const su2double cb1_k2 = var.cb1 / var.k2;
    const su2double factor = var.cw1 * var.fw - cb1_k2 * var.ft2;
    destruction = factor * pow(nue, 2) / var.dist_i_2;
    jacobian -= ((var.cw1 * var.d_fw - cb1_k2 * var.d_ft2) * pow(nue, 2) + factor * 2 * nue) / var.dist_i_2;
  }

  static void ComputeCrossProduction(const su2double& nue, const CSAVariables& var, su2double& cross_production,
                                     su2double&) {
    cross_production = var.cb2_sigma * var.norm2_Grad;
    /*--- No contribution to the jacobian. ---*/
  }
};

/*! \brief Negative. */
struct Neg {
  static void get(const su2double& nue, const CSAVariables& var, su2double& production, su2double& destruction,
                  su2double& cross_production, su2double& jacobian) {
    if (nue > 0.0) {
      Bsl::get(nue, var, production, destruction, cross_production, jacobian);
    } else {
      ComputeProduction(nue, var, production, jacobian);
      ComputeDestruction(nue, var, destruction, jacobian);
      ComputeCrossProduction(nue, var, cross_production, jacobian);
    }
  }

  static void ComputeProduction(const su2double& nue, const CSAVariables& var, su2double& production,
                                su2double& jacobian) {
    const su2double dP_dnu = var.gamma_bc * var.cb1 * (1.0 - var.ct3) * var.S;
    production = dP_dnu * nue;
    jacobian += dP_dnu;
  }

  static void ComputeDestruction(const su2double& nue, const CSAVariables& var, su2double& destruction,
                                 su2double& jacobian) {
    const su2double dD_dnu = var.cw1 * nue / var.dist_i_2;
    destruction = dD_dnu * nue;
    jacobian -= 2 * dD_dnu;
  }

  static void ComputeCrossProduction(const su2double& nue, const CSAVariables& var, su2double& cross_production,
                                     su2double& jacobian) {
    Bsl::ComputeCrossProduction(nue, var, cross_production, jacobian);
  }
};
}  // namespace SourceTerms

/* =============================================================================
 * SPALART-ALLMARAS ADDITIONAL SOURCE TERMS DECORATORS
 * ============================================================================*/

/*!
 * \class CCompressibilityCorrection
 * \brief Mixing Layer Compressibility Correction (SA-comp).
 */
template <class ParentClass>
class CCompressibilityCorrection final : public ParentClass {
 private:
  using ParentClass::Gamma;
  using ParentClass::idx;
  using ParentClass::nDim;
  using ParentClass::PrimVar_Grad_i;
  using ParentClass::ScalarVar_i;
  using ParentClass::V_i;
  using ParentClass::Volume;

  using ResidualType = typename ParentClass::template ResidualType<>;

  const su2double c5 = 3.5;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCompressibilityCorrection(unsigned short nDim, unsigned short, const CConfig* config)
      : ParentClass(nDim, 0, config) {}

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType ComputeResidual(const CConfig* config) override {
    /*--- Residual from standard SA ---*/
    ParentClass::ComputeResidual(config);

    /*--- Compressibility Correction term ---*/
    const auto& pressure = V_i[idx.Pressure()];
    const auto& density = V_i[idx.Density()];
    const su2double sound_speed = sqrt(pressure * Gamma / density);
    su2double aux_cc = 0;
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      for (unsigned short jDim = 0; jDim < nDim; ++jDim) {
        aux_cc += pow(PrimVar_Grad_i[idx.Velocity() + iDim][jDim], 2);
      }
    }
    const su2double d_CompCorrection = 2.0 * c5 * ScalarVar_i[0] / pow(sound_speed, 2) * aux_cc * Volume;
    const su2double CompCorrection = 0.5 * ScalarVar_i[0] * d_CompCorrection;

    this->Residual -= CompCorrection;
    this->Jacobian_i[0] -= d_CompCorrection;

    return ResidualType(&this->Residual, &this->Jacobian_i, nullptr);
  }
};

}  // namespace detail

/* =============================================================================
 * SPALART-ALLMARAS CLASSES
 * ============================================================================*/

/// TODO: Factory method to create combinations of the different variations based on the config.

/*!
 * \class CSourcePieceWise_TurbSA
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 */
template <class FlowIndices>
using CSourcePieceWise_TurbSA = CSourceBase_TurbSA<FlowIndices, detail::Omega::Bsl, detail::ft2::Zero,
                                                   detail::ModVort::Bsl, detail::r::Bsl, detail::SourceTerms::Bsl>;

/*!
 * \class CSourcePieceWise_TurbSA_COMP
 * \brief Class for integrating the source terms of the Spalart-Allmaras model with compressibility correction.
 */
template <class FlowIndices>
using CSourcePieceWise_TurbSA_COMP = detail::CCompressibilityCorrection<CSourcePieceWise_TurbSA<FlowIndices>>;

/*!
 * \class CSourcePieceWise_TurbSA_E
 * \brief Class for integrating the source terms of the Spalart-Allmaras model with Edwards modification.
 * \note From NASA Turbulence model site. http://turbmodels.larc.nasa.gov/spalart.html
 * This form was developed primarily to improve the near-wall numerical behavior of the model (i.e. the goal was to
 * improve the convergence behavior). The reference is: Edwards, J. R. and Chandra, S. "Comparison of Eddy
 * Viscosity-Transport Turbulence Models for Three-Dimensional, Shock-Separated Flowfields," AIAA Journal, Vol. 34, No.
 * 4, 1996, pp. 756-763.
 */
template <class FlowIndices>
using CSourcePieceWise_TurbSA_E = CSourceBase_TurbSA<FlowIndices, detail::Omega::Edw, detail::ft2::Zero,
                                                     detail::ModVort::Edw, detail::r::Edw, detail::SourceTerms::Bsl>;

/*!
 * \class CSourcePieceWise_TurbSA_E_COMP
 * \brief Class for integrating the source terms of the Spalart-Allmaras model with Edwards modification
 *        and compressibility correction.
 */
template <class FlowIndices>
using CSourcePieceWise_TurbSA_E_COMP = detail::CCompressibilityCorrection<CSourcePieceWise_TurbSA_E<FlowIndices>>;

/*!
 * \class CSourcePieceWise_TurbSA_Neg
 * \brief Class for integrating the source terms of the negative Spalart-Allmaras model.
 */
template <class FlowIndices>
using CSourcePieceWise_TurbSA_Neg = CSourceBase_TurbSA<FlowIndices, detail::Omega::Bsl, detail::ft2::Zero,
                                                       detail::ModVort::Neg, detail::r::Bsl, detail::SourceTerms::Neg>;
