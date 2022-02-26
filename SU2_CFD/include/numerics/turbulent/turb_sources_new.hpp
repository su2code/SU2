/*!
 * \file turb_sources.hpp
 * \brief Numerics classes for integration of source terms in turbulence problems.
 * \author F. Palacios, T. Economon
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

#include "../scalar/scalar_sources.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

/*!
 * \class CommonSAVariables
 * \brief Structure with SA common auxiliary functions and constants.
 */
struct CommonSAVariables {
  /*--- List of constants ---*/
  const su2double cv1_3 = pow(7.1, 3.0), k2 = pow(0.41, 2.0), cb1 = 0.1355, cw2 = 0.3, ct3 = 1.2, ct4 = 0.5,
                  cw3_6 = pow(2.0, 6.0), sigma = 2. / 3., cb2 = 0.622, cb2_sigma = cb2 / sigma,
                  cw1 = cb1 / k2 + (1.0 + cb2) / sigma, cr1 = 0.5;

  /*--- List of auxiliary functions ---*/
  su2double ft2, d_ft2, r, d_r, g, d_g, glim, fw, d_fw, Ji, d_Ji, S, Shat, d_Shat, fv1, d_fv1, fv2, d_fv2;

  /*--- List of helpers ---*/
  su2double Omega, dist_i_2, inv_k2_d2, inv_Shat, g_6, norm2_Grad;
};

/*!
 * \class CSourceBase_TurbSA
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * The variables that are subject to change in each variation/correction have their own class.
 * \note Additional source terms (e.g. compressibility) are implemented as decorators.
 */
template <class FlowIndices, class OmegaType, class ft2Type, class ModVortType, class rType, class SourceTermsType>
class CSourceBase_TurbSA_ : public CNumerics {
 protected:
  su2double Gamma_BC = 0.0;
  su2double intermittency = 0.0;

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
  CSourceBase_TurbSA_(unsigned short nDim, unsigned short, const CConfig* config)
      : CNumerics(nDim, 1, config),
        idx(nDim, config->GetnSpecies()),
        rotating_frame(config->GetRotating_Frame()),
        transition(config->GetKind_Trans_Model() == TURB_TRANS_MODEL::BC) {}

  /*!
   * \brief Residual for source term integration.
   * \param[in] intermittency_in - Value of the intermittency.
   */
  void SetIntermittency(su2double intermittency_in) final { intermittency = intermittency_in; }

  /*!
   * \brief  Get the intermittency for the BC trans. model.
   * \return Value of the intermittency.
   */
  su2double GetGammaBC() const final { return Gamma_BC; }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override {
    /*--- Model common auxiliary and constant variables ---*/
    CommonSAVariables model_var;

    const auto& Density_i = V_i[idx.Density()];
    const auto& Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];

    Residual = 0.0;
    Jacobian_i[0] = 0.0;
    su2double Production = 0.0, Destruction = 0.0, CrossProduction = 0.0;

    /*--- Evaluate Omega ---*/
    if (nDim == 2) {
      OmegaType::template get<2>(Vorticity_i, PrimVar_Grad_i[idx.Velocity()], model_var);
    } else {
      OmegaType::template get<3>(Vorticity_i, PrimVar_Grad_i[idx.Velocity()], model_var);
    }

    /*--- Rotational correction term ---*/

    if (rotating_frame) {
      model_var.Omega += 2.0 * min(0.0, StrainMag_i - model_var.Omega);
    }

    if (dist_i > 1e-10) {
      /*--- Vorticity ---*/
      model_var.S = model_var.Omega;

      model_var.dist_i_2 = dist_i * dist_i;
      const su2double nu = Laminar_Viscosity_i / Density_i;
      model_var.inv_k2_d2 = 1.0 / (model_var.k2 * model_var.dist_i_2);

      /*--- Modified values for roughness ---*/
      /*--- Ref: Aupoix, B. and Spalart, P. R., "Extensions of the Spalart-Allmaras Turbulence Model to Account for Wall
       * Roughness," International Journal of Heat and Fluid Flow, Vol. 24, 2003, pp. 454-462. ---*/
      /* --- See https://turbmodels.larc.nasa.gov/spalart.html#sarough for detailed explanation. ---*/
      model_var.Ji =
          ScalarVar_i[0] / nu +
          model_var.cr1 * (roughness_i / (dist_i + EPS));  // roughness_i = 0 for smooth walls and Ji remains the
                                                           // same, changes only if roughness is specified.nue/nul
      model_var.d_Ji = 1.0 / nu;

      const su2double Ji_2 = pow(model_var.Ji, 2);
      const su2double Ji_3 = Ji_2 * model_var.Ji;

      model_var.fv1 = Ji_3 / (Ji_3 + model_var.cv1_3);
      model_var.d_fv1 = 3.0 * Ji_2 * model_var.cv1_3 / (nu * pow(Ji_3 + model_var.cv1_3, 2));

      /*--- Using a modified relation so as to not change the Shat that depends on fv2. ---*/
      model_var.fv2 =
          1.0 -
          ScalarVar_i[0] / (nu + ScalarVar_i[0] * model_var.fv1);  // From NASA turb modeling resource and 2003 paper
      model_var.d_fv2 = -(1.0 / nu - Ji_2 * model_var.d_fv1) / pow(1.0 + model_var.Ji * model_var.fv1, 2);

      /*--- Compute ft2 term ---*/
      ft2Type::get(model_var);

      /*--- Compute modified vorticity ---*/
      ModVortType::get(ScalarVar_i[0], nu, model_var);

      model_var.inv_Shat = 1.0 / model_var.Shat;

      /*--- Compute auxiliary function r ---*/
      rType::get(ScalarVar_i[0], model_var);

      model_var.g = model_var.r + model_var.cw2 * (pow(model_var.r, 6.0) - model_var.r);
      model_var.g_6 = pow(model_var.g, 6.0);
      model_var.glim = pow((1.0 + model_var.cw3_6) / (model_var.g_6 + model_var.cw3_6), 1.0 / 6.0);
      model_var.fw = model_var.g * model_var.glim;

      model_var.d_g = model_var.d_r * (1. + model_var.cw2 * (6.0 * pow(model_var.r, 5.0) - 1.0));
      model_var.d_fw = model_var.d_g * model_var.glim * (1. - model_var.g_6 / (model_var.g_6 + model_var.cw3_6));

      model_var.norm2_Grad = GeometryToolbox::SquaredNorm(nDim, ScalarVar_Grad_i[0]);

      /*--- Compute production, destruction and cross production and jacobian ---*/
      SourceTermsType::get(ScalarVar_i[0], model_var, Production, Destruction, CrossProduction, Jacobian_i[0]);

      /*--- Compute any necessary additional source term and jacobian contribution ---*/

      /*--- Residual ---*/
      Residual = Production - Destruction + CrossProduction;
      Residual *= Volume;

      /*--- Jacobian ---*/
      Jacobian_i[0] *= Volume;
    }

    //  AD::SetPreaccOut(Residual);
    //  AD::EndPreacc();

    return ResidualType<>(&Residual, &Jacobian_i, nullptr);
  }
};

/* =============================================================================
 * SPALART-ALLMARAS VARIANT CLASSES
 * ============================================================================*/

namespace detail {

/*!
 * \brief Strain rate classes.
 * \tparam nDim: Problem dimension.
 * \param[in] vorticity: Vorticity array.
 * \param[in] velocity_grad: Velocity gradients.
 * \param[out] model_var: Common SA variables struct (to set Omega).
 */
namespace Omega {

/*! \brief Baseline. */
struct Bsl {
  template <int nDim, class MatrixType>
  static void get(const su2double* vorticity, const MatrixType&, CommonSAVariables& model_var) {
    model_var.Omega = GeometryToolbox::Norm(3, vorticity);
  }
};

/*! \brief Edward. Here Omega is the Strain Rate. */
struct Edw {
  template <int nDim, class MatrixType>
  static void get(const su2double*, const MatrixType& velocity_grad, CommonSAVariables& model_var) {
    su2double Sbar = 0.0;
    for (int iDim = 0; iDim < nDim; ++iDim) {
      for (int jDim = 0; jDim < nDim; ++jDim) {
        Sbar += (velocity_grad[iDim][jDim] + velocity_grad[jDim][iDim]) * velocity_grad[iDim][jDim];
      }
    }
    for (int iDim = 0; iDim < nDim; ++iDim) {
      Sbar -= (2.0 / 3.0) * pow(velocity_grad[iDim][iDim], 2);
    }
    model_var.Omega = sqrt(max(Sbar, 0.0));
  }
};
}  // namespace Omega

/*!
 * \brief Classes to set the ft2 term and its derivative.
 * \param[in,out] model_var: Common SA variables struct.
 */
namespace ft2 {

/*! \brief SU2 baseline ft2 term value and its derivative (ft2=0.0). */
struct Bsl {
  static void get(CommonSAVariables& model_var) {
    model_var.ft2 = 0.0;
    model_var.d_ft2 = 0.0;
  }
};

/*! \brief Non-zero ft2 term according to the literature. */
struct Nonzero {
  static void get(CommonSAVariables& model_var) {
    const su2double xsi2 = pow(model_var.Ji, 2);
    model_var.ft2 = model_var.ct3 * exp(-model_var.ct4 * xsi2);
    model_var.d_ft2 = -2.0 * model_var.ct4 * model_var.Ji * model_var.ft2 * model_var.d_Ji;
  }
};
}  // namespace ft2

/*!
 * \brief Classes to compute the modified vorticity (\tilde{S}) and its derivative.
 * \param[in] nue: SA variable.
 * \param[in] nu: Laminar viscosity.
 * \param[in,out] model_var: Common SA variables struct.
 */
namespace ModVort {

/*! \brief Baseline. */
struct Bsl {
  static void get(const su2double& nue, const su2double& nu, CommonSAVariables& model_var) {
    const su2double Sbar = nue * model_var.fv2 * model_var.inv_k2_d2;
    model_var.Shat = model_var.S + Sbar;
    model_var.Shat = max(model_var.Shat, 1.0e-10);

    const su2double d_Sbar = (model_var.fv2 + nue * model_var.d_fv2) * model_var.inv_k2_d2;
    model_var.d_Shat = (model_var.Shat <= 1.0e-10) ? 0.0 : d_Sbar;
  }
};

/*! \brief Edward. */
struct Edw {
  static void get(const su2double& nue, const su2double& nu, CommonSAVariables& model_var) {
    model_var.Shat = max(model_var.S * ((1.0 / max(model_var.Ji, 1.0e-16)) + model_var.fv1), 1.0e-16);
    model_var.Shat = max(model_var.Shat, 1.0e-10);

    model_var.d_Shat =
        (model_var.Shat <= 1.0e-10) ? 0.0 : -model_var.S * pow(model_var.Ji, -2) / nu + model_var.S * model_var.d_fv1;
  }
};

/*! \brief Negative. */
struct Neg {
  static void get(const su2double& nue, const su2double& nu, CommonSAVariables& model_var) {
    if (nue > 0.0) {
      // Baseline solution
      Bsl::get(nue, nu, model_var);
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
 * \param[in,out] model_var: Common SA variables struct.
 */
namespace r {

/*! \brief Baseline. */
struct Bsl {
  static void get(const su2double& nue, CommonSAVariables& model_var) {
    model_var.r = min(nue * model_var.inv_Shat * model_var.inv_k2_d2, 10.0);
    model_var.d_r = (model_var.Shat - nue * model_var.d_Shat) * pow(model_var.inv_Shat, 2) * model_var.inv_k2_d2;
    if (model_var.r == 10.0) model_var.d_r = 0.0;
  }
};

/*! \brief Edward. */
struct Edw {
  static void get(const su2double& nue, CommonSAVariables& model_var) {
    model_var.r = min(nue * model_var.inv_Shat * model_var.inv_k2_d2, 10.0);
    model_var.r = tanh(model_var.r) / tanh(1.0);

    model_var.d_r = (model_var.Shat - nue * model_var.d_Shat) * pow(model_var.inv_Shat, 2) * model_var.inv_k2_d2;
    model_var.d_r = (1 - pow(tanh(model_var.r), 2.0)) * (model_var.d_r) / tanh(1.0);
  }
};
}  // namespace r

/*!
 * \brief Source terms classes: production, destruction and cross-productions term and their derivative.
 * \param[in] nue: SA variable.
 * \param[in] model_var: Common SA variables struct.
 * \param[out] production: Production term.
 * \param[out] destruction: Destruction term.
 * \param[out] cross_production: CrossProduction term.
 * \param[out] jacobian: Derivative of the combined source term wrt nue.
 */
namespace SourceTerms {

/*! \brief Baseline (Original SA model). */
struct Bsl {
  static void get(const su2double& nue, const CommonSAVariables& model_var, su2double& production,
                  su2double& destruction, su2double& cross_production, su2double& jacobian) {
    ComputeProduction(nue, model_var, production, jacobian);

    ComputeDestruction(nue, model_var, destruction, jacobian);

    ComputeCrossProduction(nue, model_var, cross_production, jacobian);
  }

  static void ComputeProduction(const su2double& nue, const CommonSAVariables& model_var, su2double& production,
                                su2double& jacobian) {
    production = model_var.cb1 * (1.0 - model_var.ft2) * model_var.Shat * nue;
    jacobian += model_var.cb1 * (-model_var.Shat * nue * model_var.d_ft2 +
                                 (1.0 - model_var.ft2) * (nue * model_var.d_Shat + model_var.Shat));
  }

  static void ComputeDestruction(const su2double& nue, const CommonSAVariables& model_var, su2double& destruction,
                                 su2double& jacobian) {
    destruction =
        (model_var.cw1 * model_var.fw - model_var.cb1 * model_var.ft2 / model_var.k2) * nue * nue / model_var.dist_i_2;
    jacobian -=
        (model_var.cw1 * model_var.d_fw - model_var.cb1 / model_var.k2 * model_var.d_ft2) * nue * nue /
            model_var.dist_i_2 +
        (model_var.cw1 * model_var.fw - model_var.cb1 * model_var.ft2 / model_var.k2) * 2.0 * nue / model_var.dist_i_2;
  }

  static void ComputeCrossProduction(const su2double& nue, const CommonSAVariables& model_var,
                                     su2double& cross_production, su2double&) {
    cross_production = model_var.cb2_sigma * model_var.norm2_Grad;
    /*--- No contribution to the jacobian. ---*/
  }
};

/*! \brief Negative. */
struct Neg {
  static void get(const su2double& nue, const CommonSAVariables& model_var, su2double& production,
                  su2double& destruction, su2double& cross_production, su2double& jacobian) {
    if (nue > 0.0) {
      Bsl::get(nue, model_var, production, destruction, cross_production, jacobian);
    } else {
      ComputeProduction(nue, model_var, production, jacobian);

      ComputeDestruction(nue, model_var, destruction, jacobian);

      ComputeCrossProduction(nue, model_var, cross_production, jacobian);
    }
  }

  static void ComputeProduction(const su2double& nue, const CommonSAVariables& model_var, su2double& production,
                                su2double& jacobian) {
    const su2double dP_dnu = model_var.cb1 * (1.0 - model_var.ct3) * model_var.S;
    production = dP_dnu * nue;
    jacobian += dP_dnu;
  }

  static void ComputeDestruction(const su2double& nue, const CommonSAVariables& model_var, su2double& destruction,
                                 su2double& jacobian) {
    const su2double dD_dnu = model_var.cw1 * nue / model_var.dist_i_2;
    destruction = dD_dnu * nue;
    jacobian -= 2 * dD_dnu;
  }

  static void ComputeCrossProduction(const su2double& nue, const CommonSAVariables& model_var,
                                     su2double& cross_production, su2double& jacobian) {
    Bsl::ComputeCrossProduction(nue, model_var, cross_production, jacobian);
  }
};
}  // namespace SourceTerms

}  // namespace detail

/* =============================================================================
 * SPALART-ALLMARAS ADDITIONAL SOURCE TERMS DECORATORS
 * ============================================================================*/

/*!
 * \class CompressiblityCorrection
 * \brief Mixing Layer Compressibility Correction (SA-comp).
 */
template <class ParentClass>
class CompressiblityCorrection final : public ParentClass {
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
  CompressiblityCorrection(unsigned short nDim, unsigned short, const CConfig* config) : ParentClass(nDim, 0, config) {}

  ResidualType ComputeResidual(const CConfig* config) override {
    /*--- Residual from standard SA ---*/
    ParentClass::ComputeResidual(config);

    /*--- Compressibility Correction term ---*/
    const auto& Pressure_i = V_i[idx.Pressure()];
    const auto& Density_i = V_i[idx.Density()];
    const su2double SoundSpeed_i = sqrt(Pressure_i * Gamma / Density_i);
    su2double aux_cc = 0;
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      for (unsigned short jDim = 0; jDim < nDim; ++jDim) {
        aux_cc += pow(PrimVar_Grad_i[idx.Velocity() + iDim][jDim], 2);
      }
    }
    const su2double d_CompCorrection = 2.0 * c5 * ScalarVar_i[0] / pow(SoundSpeed_i, 2) * aux_cc * Volume;
    const su2double CompCorrection = 0.5 * ScalarVar_i[0] * d_CompCorrection;

    this->Residual -= CompCorrection;
    this->Jacobian_i[0] -= d_CompCorrection;

    return ResidualType(&this->Residual, &this->Jacobian_i, nullptr);
  }
};
