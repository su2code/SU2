/*!
 * \file turb_sources.hpp
 * \brief Numerics classes for integration of source terms in turbulence problems.
 * \version 7.3.1 "Blackbird"
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

#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../scalar/scalar_sources.hpp"

/*!
 * \class CSAVariables
 * \brief Structure with SA common auxiliary functions and constants.
 */
struct CSAVariables {
  /*--- List of constants ---*/
  const su2double cv1_3 = pow(7.1, 3);
  const su2double k2 = pow(0.41, 2);
  const su2double cb1 = 0.1355;
  const su2double cw2 = 0.3;
  const su2double ct3 = 1.2;
  const su2double ct4 = 0.5;
  const su2double cw3_6 = pow(2, 6);
  const su2double sigma = 2.0 / 3.0;
  const su2double cb2 = 0.622;
  const su2double cb2_sigma = cb2 / sigma;
  const su2double cw1 = cb1 / k2 + (1 + cb2) / sigma;
  const su2double cr1 = 0.5;

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
    const auto& density = V_i[idx.Density()];
    const auto& laminar_viscosity = V_i[idx.LaminarViscosity()];

    AD::StartPreacc();
    AD::SetPreaccIn(density, laminar_viscosity, StrainMag_i, ScalarVar_i[0], Volume, dist_i, roughness_i);
    AD::SetPreaccIn(Vorticity_i, 3);
    AD::SetPreaccIn(PrimVar_Grad_i + idx.Velocity(), nDim, nDim);
    AD::SetPreaccIn(ScalarVar_Grad_i[0], nDim);

    /*--- Common auxiliary variables and constants of the model. ---*/
    CSAVariables var;

    Residual = 0.0;
    Jacobian_i[0] = 0.0;

    /*--- Evaluate Omega with a rotational correction term. ---*/

    Omega::get(Vorticity_i, nDim, PrimVar_Grad_i + idx.Velocity(), var);

    /*--- Dacles-Mariani et. al. rotation correction ("-R"), this is applied by
     * default for rotating frame, but should be controled in the config. ---*/
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

      if (transition) {
        /*--- BC transition model (2020 revision). This should only be used with SA-noft2.
         * TODO: Consider making this part of the "SourceTerms" template. ---*/
        const su2double chi_1 = 0.002;
        const su2double chi_2 = 50.0;

        /*--- Turbulence intensity is u'/U so we multiply by 100 to get percentage. ---*/
        const su2double tu = 100.0 * config->GetTurbulenceIntensity_FreeStream();
        const su2double nu_t = ScalarVar_i[0] * var.fv1;

        const su2double re_v = density * var.dist_i_2 / laminar_viscosity * var.Omega;
        const su2double re_theta = re_v / 2.193;
        /*--- Menter correlation. ---*/
        const su2double re_theta_t = 803.73 * pow(tu + 0.6067, -1.027);

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

    AD::SetPreaccOut(Residual);
    AD::EndPreacc();

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
    /*--- The destruction when nue < 0 is added instead of the usual subtraction, hence the negative sign. ---*/
    const su2double dD_dnu = -var.cw1 * nue / var.dist_i_2;
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
/// See PR #1413, the combinations exposed should follow https://turbmodels.larc.nasa.gov/spalart.html

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

/*!
 * \class CSourcePieceWise_TurbSST
 * \brief Class for integrating the source terms of the Menter SST turbulence model equations.
 */
template <class FlowIndices>
class CSourcePieceWise_TurbSST final : public CNumerics {
 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */
  const bool sustaining_terms = false;
  const bool axisymmetric = false;

  /*--- Closure constants ---*/
  const su2double sigma_k_1, sigma_k_2, sigma_w_1, sigma_w_2, beta_1, beta_2, beta_star, a1, alfa_1, alfa_2;

  /*--- Ambient values for SST-SUST. ---*/
  const su2double kAmb, omegaAmb;

  su2double F1_i, F2_i, CDkw_i;
  su2double Residual[2];
  su2double* Jacobian_i[2];
  su2double Jacobian_Buffer[4];  /// Static storage for the Jacobian (which needs to be pointer for return type).

  /*!
   * \brief Get strain magnitude based on perturbed reynolds stress matrix.
   * \param[in] turb_ke: turbulent kinetic energy of the node.
   */
  inline su2double PerturbedStrainMag(su2double turb_ke) const {
    /*--- Compute norm of perturbed strain rate tensor. ---*/

    su2double perturbedStrainMag = 0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      for (unsigned short jDim = 0; jDim < nDim; jDim++) {
        su2double StrainRate_ij = MeanPerturbedRSM[iDim][jDim] - TWO3 * turb_ke * delta[iDim][jDim];
        StrainRate_ij = -StrainRate_ij * Density_i / (2 * Eddy_Viscosity_i);

        perturbedStrainMag += pow(StrainRate_ij, 2.0);
      }
    }
    return sqrt(2.0 * perturbedStrainMag);
  }

  /*!
   * \brief Add contribution due to axisymmetric formulation to 2D residual
   */
  inline void ResidualAxisymmetric(su2double alfa_blended, su2double zeta) {
    if (Coord_i[1] < EPS) return;

    AD::SetPreaccIn(Coord_i[1]);
    AD::SetPreaccIn(V_i[idx.Velocity() + 1]);

    const su2double yinv = 1.0 / Coord_i[1];
    const su2double rhov = Density_i * V_i[idx.Velocity() + 1];
    const su2double& k = ScalarVar_i[0];
    const su2double& w = ScalarVar_i[1];

    /*--- Compute blended constants ---*/
    const su2double sigma_k_i = F1_i * sigma_k_1 + (1.0 - F1_i) * sigma_k_2;
    const su2double sigma_w_i = F1_i * sigma_w_1 + (1.0 - F1_i) * sigma_w_2;

    /*--- Production ---*/
    const su2double pk_axi = max(
        0.0, 2.0 / 3.0 * rhov * k * ((2.0 * yinv * V_i[idx.Velocity() + 1] - PrimVar_Grad_i[idx.Velocity()+1][1] - PrimVar_Grad_i[idx.Velocity()][0]) / zeta - 1.0));
    const su2double pw_axi = alfa_blended * zeta / k * pk_axi;

    /*--- Convection-Diffusion ---*/
    const su2double cdk_axi = rhov * k - (Laminar_Viscosity_i + sigma_k_i * Eddy_Viscosity_i) * ScalarVar_Grad_i[0][1];
    const su2double cdw_axi = rhov * w - (Laminar_Viscosity_i + sigma_w_i * Eddy_Viscosity_i) * ScalarVar_Grad_i[1][1];

    /*--- Add terms to the residuals ---*/
    Residual[0] += yinv * Volume * (pk_axi - cdk_axi);
    Residual[1] += yinv * Volume * (pw_axi - cdw_axi);
  }

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] constants - SST model constants.
   * \param[in] val_kine_Inf - Freestream k, for SST with sustaining terms.
   * \param[in] val_omega_Inf - Freestream w, for SST with sustaining terms.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short, const su2double* constants, su2double val_kine_Inf,
                           su2double val_omega_Inf, const CConfig* config)
      : CNumerics(val_nDim, 2, config),
        idx(val_nDim, config->GetnSpecies()),
        sustaining_terms(config->GetKind_Turb_Model() == TURB_MODEL::SST_SUST),
        axisymmetric(config->GetAxisymmetric()),
        sigma_k_1(constants[0]),
        sigma_k_2(constants[1]),
        sigma_w_1(constants[2]),
        sigma_w_2(constants[3]),
        beta_1(constants[4]),
        beta_2(constants[5]),
        beta_star(constants[6]),
        a1(constants[7]),
        alfa_1(constants[8]),
        alfa_2(constants[9]),
        kAmb(val_kine_Inf),
        omegaAmb(val_omega_Inf) {
    /*--- "Allocate" the Jacobian using the static buffer. ---*/
    Jacobian_i[0] = Jacobian_Buffer;
    Jacobian_i[1] = Jacobian_Buffer + 2;
  }

  /*!
   * \brief Set the value of the first blending function.
   * \param[in] val_F1_i - Value of the first blending function at point i.
   * \param[in] Not used.
   */
  inline void SetF1blending(su2double val_F1_i, su2double) override {
    F1_i = val_F1_i;
  }

  /*!
   * \brief Set the value of the second blending function.
   * \param[in] val_F2_i - Value of the second blending function at point i.
   */
  inline void SetF2blending(su2double val_F2_i) override {
    F2_i = val_F2_i;
  }

  /*!
   * \brief Set the value of the cross diffusion for the SST model.
   * \param[in] val_CDkw_i - Value of the cross diffusion at point i.
   */
  inline void SetCrossDiff(su2double val_CDkw_i) override {
    CDkw_i = val_CDkw_i;
  }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override {
    AD::StartPreacc();
    AD::SetPreaccIn(StrainMag_i);
    AD::SetPreaccIn(ScalarVar_i, nVar);
    AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
    AD::SetPreaccIn(Volume);
    AD::SetPreaccIn(dist_i);
    AD::SetPreaccIn(F1_i);
    AD::SetPreaccIn(F2_i);
    AD::SetPreaccIn(CDkw_i);
    AD::SetPreaccIn(PrimVar_Grad_i, nDim + idx.Velocity(), nDim);
    AD::SetPreaccIn(Vorticity_i, 3);
    AD::SetPreaccIn(V_i[idx.Density()], V_i[idx.LaminarViscosity()], V_i[idx.EddyViscosity()]);

    Density_i = V_i[idx.Density()];
    Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];
    Eddy_Viscosity_i = V_i[idx.EddyViscosity()];

    Residual[0] = 0.0;
    Residual[1] = 0.0;
    Jacobian_i[0][0] = 0.0;
    Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;
    Jacobian_i[1][1] = 0.0;

    /*--- Computation of blended constants for the source terms ---*/

    const su2double alfa_blended = F1_i * alfa_1 + (1.0 - F1_i) * alfa_2;
    const su2double beta_blended = F1_i * beta_1 + (1.0 - F1_i) * beta_2;

    if (dist_i > 1e-10) {
      /*--- Production ---*/

      su2double diverg = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        diverg += PrimVar_Grad_i[iDim + idx.Velocity()][iDim];

      /*--- If using UQ methodolgy, calculate production using perturbed Reynolds stress matrix ---*/

      su2double StrainMag = StrainMag_i;

      if (using_uq) {
        ComputePerturbedRSM(nDim, Eig_Val_Comp, uq_permute, uq_delta_b, uq_urlx, PrimVar_Grad_i + idx.Velocity(),
                            Density_i, Eddy_Viscosity_i, ScalarVar_i[0], MeanPerturbedRSM);
        StrainMag = PerturbedStrainMag(ScalarVar_i[0]);
      }

      su2double pk = Eddy_Viscosity_i * pow(StrainMag, 2) - 2.0 / 3.0 * Density_i * ScalarVar_i[0] * diverg;
      pk = max(0.0, min(pk, 20.0 * beta_star * Density_i * ScalarVar_i[1] * ScalarVar_i[0]));

      const su2double VorticityMag = GeometryToolbox::Norm(3, Vorticity_i);
      const su2double zeta = max(ScalarVar_i[1], VorticityMag * F2_i / a1);
      su2double pw = alfa_blended * Density_i * max(pow(StrainMag, 2) - 2.0 / 3.0 * zeta * diverg, 0.0);

      /*--- Sustaining terms, if desired. Note that if the production terms are
            larger equal than the sustaining terms, the original formulation is
            obtained again. This is in contrast to the version in literature
            where the sustaining terms are simply added. This latter approach could
            lead to problems for very big values of the free-stream turbulence
            intensity. ---*/

      if (sustaining_terms) {
        const su2double sust_k = beta_star * Density_i * kAmb * omegaAmb;
        const su2double sust_w = beta_blended * Density_i * omegaAmb * omegaAmb;
        pk = max(pk, sust_k);
        pw = max(pw, sust_w);
      }

      /*--- Add the production terms to the residuals. ---*/

      Residual[0] += pk * Volume;
      Residual[1] += pw * Volume;

      /*--- Dissipation ---*/

      Residual[0] -= beta_star * Density_i * ScalarVar_i[1] * ScalarVar_i[0] * Volume;
      Residual[1] -= beta_blended * Density_i * ScalarVar_i[1] * ScalarVar_i[1] * Volume;

      /*--- Cross diffusion ---*/

      Residual[1] += (1.0 - F1_i) * CDkw_i * Volume;

      /*--- Contribution due to 2D axisymmetric formulation ---*/

      if (axisymmetric) ResidualAxisymmetric(alfa_blended, zeta);

      /*--- Implicit part ---*/

      Jacobian_i[0][0] = -beta_star * ScalarVar_i[1] * Volume;
      Jacobian_i[0][1] = -beta_star * ScalarVar_i[0] * Volume;
      Jacobian_i[1][0] = 0.0;
      Jacobian_i[1][1] = -2.0 * beta_blended * ScalarVar_i[1] * Volume;
    }

    AD::SetPreaccOut(Residual, nVar);
    AD::EndPreacc();

    return ResidualType<>(Residual, Jacobian_i, nullptr);
  }
};
