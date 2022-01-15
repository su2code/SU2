/*!
 * \file turb_sources.hpp
 * \brief Delarations of numerics classes for integration of source
 *        terms in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.2.1 "Blackbird"
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

/*!
 * \class CSourcePieceWise_TurbSA
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * \brief The variables that are subject to change in each variation/correction have their own class. Additional source
 * terms are implemented as decorators. \ingroup SourceDiscr \author A. Bueno.
 */
template <class FlowIndices, class ft2_class, class ModVort_class, class r_class>
class CSourceBase_TurbSA : public CNumerics {
 protected:
  /*--- List of constants ---*/
  su2double cv1_3, k2, cb1, cw2, ct3, ct4, cw3_6, cb2_sigma, sigma, cb2, cw1, cr1;

  /*--- List of auxiliary functions ---*/
  su2double ft2, d_ft2, r, d_r, g, d_g, g_6, glim, fw, Ji, Ji_2, Ji_3, d_Ji, S, Omega, Shat, d_Shat, inv_Shat, fv1, fv2;

  su2double Gamma_BC = 0.0;
  su2double intermittency;

  /*--- Source term components ---*/
  su2double Production, Destruction, CrossProduction, AddSourceTerm;

  /*--- Residual and Jacobian ---*/
  su2double Residual, *Jacobian_i;

 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */

  su2double Jacobian_Buffer;  /// Static storage for the Jacobian (which needs to be pointer for return type).

 protected:
  const bool rotating_frame = false;
  bool roughwall = false;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceBase_TurbSA(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
      : CNumerics(val_nDim, val_nVar, config), idx(val_nDim, config->GetnSpecies()) {}

  /*!
   * \brief Residual for source term integration.
   * \param[in] intermittency_in - Value of the intermittency.
   */
  inline void SetIntermittency(su2double intermittency_in) final { intermittency = intermittency_in; }

  /*!
   * \brief Residual for source term integration.
   * \param[in] val_production - Value of the Production.
   */
  inline void SetProduction(su2double val_production) final { Production = val_production; }

  /*!
   * \brief Residual for source term integration.
   * \param[in] val_destruction - Value of the Destruction.
   */
  inline void SetDestruction(su2double val_destruction) final { Destruction = val_destruction; }

  /*!
   * \brief Residual for source term integration.
   * \param[in] val_crossproduction - Value of the CrossProduction.
   */
  inline void SetCrossProduction(su2double val_crossproduction) final { CrossProduction = val_crossproduction; }

  /*!
   * \brief ______________.
   */
  inline su2double GetProduction(void) const final { return Production; }

  /*!
   * \brief  Get the intermittency for the BC trans. model.
   * \return Value of the intermittency.
   */
  inline su2double GetGammaBC(void) const final { return Gamma_BC; }

  /*!
   * \brief  ______________.
   */
  inline su2double GetDestruction(void) const final { return Destruction; }

  /*!
   * \brief  ______________.
   */
  inline su2double GetCrossProduction(void) const final { return CrossProduction; }

  /*!
   * \brief  ______________.
   */
  inline su2double ComputeXsi(const su2double nue, const su2double nul, su2double xsi, su2double d_xsi) {
    xsi = nue / nul + cr1 * (roughness_i / (dist_i + EPS));  // roughness_i = 0 for smooth walls and Ji remains the
                                                             // same, changes only if roughness is specified.nue/nul
    d_xsi = 1.0 / nul;
  }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) final {
    // Set the boolean here depending on whether the point is closest to a rough wall or not.
    roughwall = (roughness_i > 0.0);

    Density_i = V_i[idx.Density()];
    Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];

    Residual = 0.0;
    Production = 0.0;
    Destruction = 0.0;
    CrossProduction = 0.0;
    Jacobian_i[0] = 0.0;

    /*--- Evaluate Omega ---*/

    Omega = sqrt(Vorticity_i[0] * Vorticity_i[0] + Vorticity_i[1] * Vorticity_i[1] + Vorticity_i[2] * Vorticity_i[2]);

    /*--- Rotational correction term ---*/

    if (rotating_frame) {
      Omega += 2.0 * min(0.0, StrainMag_i - Omega);
    }

    if (dist_i > 1e-10) {
      /*--- Production term ---*/

      dist_i_2 = dist_i * dist_i;
      nu = Laminar_Viscosity_i / Density_i;

      /*--- Modified values for roughness ---*/
      /*--- Ref: Aupoix, B. and Spalart, P. R., "Extensions of the Spalart-Allmaras Turbulence Model to Account for Wall
       * Roughness," International Journal of Heat and Fluid Flow, Vol. 24, 2003, pp. 454-462. ---*/
      /* --- See https://turbmodels.larc.nasa.gov/spalart.html#sarough for detailed explanation. ---*/

      ComputeXsi(ScalarVar_i[0], Laminar_Viscosity_i, &Ji, &d_Ji);
      Ji_2 = Ji * Ji;
      Ji_3 = Ji_2 * Ji;
      fv1 = Ji_3 / (Ji_3 + cv1_3);

      /*--- Using a modified relation so as to not change the Shat that depends on fv2. ---*/
      fv2 = 1.0 - ScalarVar_i[0] / (nu + ScalarVar_i[0] * fv1);  // From NASA turb modeling resource and 2003 paper

      ft2 = ct3 * exp(-ct4 * Ji_2);
      S = Omega;
      inv_k2_d2 = 1.0 / (k2 * dist_i_2);

      Shat = S + ScalarVar_i[0] * fv2 * inv_k2_d2;
      Shat = max(Shat, 1.0e-10);
      inv_Shat = 1.0 / Shat;

      //    Original SA model
      //    Production = cb1*(1.0-ft2)*Shat*ScalarVar_i[0]*Volume;

      if (transition) {
        /*--- BC model constants (2020 revision). ---*/
        const su2double chi_1 = 0.002;
        const su2double chi_2 = 50.0;

        /*--- turbulence intensity is u'/U so we multiply by 100 to get percentage ---*/
        su2double tu = 100.0 * config->GetTurbulenceIntensity_FreeStream();

        su2double nu_t = (ScalarVar_i[0] * fv1);  // S-A variable

        su2double re_v = ((Density_i * pow(dist_i, 2.)) / (Laminar_Viscosity_i)) * Omega;
        su2double re_theta = re_v / 2.193;
        su2double re_theta_t = (803.73 * pow((tu + 0.6067), -1.027));  // MENTER correlation
        // re_theta_t = 163.0 + exp(6.91-tu); //ABU-GHANNAM & SHAW correlation

        su2double term1 = sqrt(max(re_theta - re_theta_t, 0.) / (chi_1 * re_theta_t));
        su2double term2 = sqrt(max((nu_t * chi_2) / nu, 0.));
        su2double term_exponential = (term1 + term2);

        Gamma_BC = 1.0 - exp(-term_exponential);

        Production = Gamma_BC * cb1 * Shat * ScalarVar_i[0] * Volume;
      } else {
        Production = cb1 * Shat * ScalarVar_i[0] * Volume;
      }

      /*--- Destruction term ---*/

      r = min(ScalarVar_i[0] * inv_Shat * inv_k2_d2, 10.0);
      g = r + cw2 * (pow(r, 6.0) - r);
      g_6 = pow(g, 6.0);
      glim = pow((1.0 + cw3_6) / (g_6 + cw3_6), 1.0 / 6.0);
      fw = g * glim;

      //    Original SA model
      //    Destruction = (cw1*fw-cb1*ft2/k2)*ScalarVar_i[0]*ScalarVar_i[0]/dist_i_2*Volume;

      Destruction = cw1 * fw * ScalarVar_i[0] * ScalarVar_i[0] / dist_i_2 * Volume;

      /*--- Diffusion term ---*/

      norm2_Grad = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) norm2_Grad += ScalarVar_Grad_i[0][iDim] * ScalarVar_Grad_i[0][iDim];

      CrossProduction = cb2_sigma * norm2_Grad * Volume;

      Residual = Production - Destruction + CrossProduction;

      /*--- Implicit part, production term ---*/

      dfv1 = 3.0 * Ji_2 * cv1_3 / (nu * pow(Ji_3 + cv1_3, 2.));
      dfv2 = -(1 / nu - Ji_2 * dfv1) / pow(1. + Ji * fv1, 2.);
      if (Shat <= 1.0e-10)
        dShat = 0.0;
      else
        dShat = (fv2 + ScalarVar_i[0] * dfv2) * inv_k2_d2;

      if (transition) {
        Jacobian_i[0] += Gamma_BC * cb1 * (ScalarVar_i[0] * dShat + Shat) * Volume;
      } else {
        Jacobian_i[0] += cb1 * (ScalarVar_i[0] * dShat + Shat) * Volume;
      }

      /*--- Implicit part, destruction term ---*/

      dr = (Shat - ScalarVar_i[0] * dShat) * inv_Shat * inv_Shat * inv_k2_d2;
      if (r == 10.0) dr = 0.0;
      dg = dr * (1. + cw2 * (6.0 * pow(r, 5.0) - 1.0));
      dfw = dg * glim * (1. - g_6 / (g_6 + cw3_6));
      Jacobian_i[0] -= cw1 * (dfw * ScalarVar_i[0] + 2.0 * fw) * ScalarVar_i[0] / dist_i_2 * Volume;

      ft2_class::get(&ft2, &d_ft2);

      su2double pingu = ft2 * 3.14;

      ModVort_class::get();
    }
  }
};

/*------------------------------------------------------------------------------
| Compute the ft2-term and its derivative
| * \param[in] ct3 - constant ct3
| * \param[in] ct4 - constant ct4
| * \param[in] nue - nu_tilde
| * \param[in] nul - dynamic laminar viscosity
| * \param[out] ft2 - ft2 variable
| * \param[out] d_ft2 - derivative of ft2 w.r.t. nue variable
------------------------------------------------------------------------------*/

/*!
 * \brief SU2 baseline ft2 term value and its derivative. ft2=0.0
 */
template <class Base>
class ft2_Bsl : public Base {
  static double get(su2double* ft2, su2double* d_ft2) {
    *ft2 = 0.0;
    *d_ft2 = 0.0;
  }
};

/*!
 * \brief non-zero ft2 term according to the literature and its derivative.
 */
template <class Base>
class ft2_nonzero : public Base {
  using Base::ct3;
  using Base::ct4;
  using Base::d_xsi;
  using Base::xsi;

  static void get(su2double* ft2, su2double* d_ft2) {
    const su2double xsi2 = xsi * xsi;

    *ft2 = ct3 * exp(-ct4 * xsi2);
    *d_ft2 = -2.0 * ct4 * xsi * *ft2 * d_xsi;
  }
};

/* The above causes a cyclic dependency between the turbulence base and the ft2 and co. classes.
 * To break that you need to template the get method instead of the class, for example: */
struct ft2_nonzero {
  template <class Base>
  static void get(const Base& base, su2double& ft2, su2double& d_ft2) {
    const su2double xsi2 = pow(base.xsi, 2);
    ft2 = base.ct3 * exp(-base.ct4 * xsi2);
    d_ft2 = -2.0 * base.ct4 * base.xsi * ft2 * base.d_xsi;
  }
};
/* Now you pass the numerics class itself to "get" to access the member variables,
 * which you'll have to make public, for example in ComputeResidual:
 * ft2_class::get(*this, ft2, d_ft2);
 *
 * This is not a perfect solution because for performance we should cut down on "aux" class variables.
 * An alternative is to put such variables in a struct, which can then be a local variable in ComputeResidual. */
struct CommonVariables {
  su2double ft2, d_ft2, r, d_r, g, d_g, g_6, glim, fw, Ji, Ji_2, Ji_3, d_Ji, S, Omega, Shat, d_Shat, inv_Shat, fv1, fv2;
};
/*
 *  ResidualType<> ComputeResidual(const CConfig* config) final {
 *    CommonVariables data;
 *    ...
 *    data.Ji = ...
 *    ...
 *    ft2_class::get(data, ft2, d_ft2);
 */
/* This way you can also dispense with templating "get". */

/*------------------------------------------------------------------------------
| Compute the modified vorticity (\tilde{S}) and its derivative
| * \param[in] S - vorticity (Omega)
| * \param[in] PrimVar_Grad_i - Gradient of primitive variables at point i
| * \param[in] nul - dynamic laminar viscosity
| * \param[in] nue - nu_tilde
| * \param[in] fv1 - auxiliary function fv1
| * \param[in] d_fv1 - derivative of fv1 w.r.t. nue variable
| * \param[in] fv2 - auxiliary function fv2
| * \param[in] d_fv2 - derivative of fv2 w.r.t. nue variable
| * \param[in] inv_k2_d2 - 1/(k^2*d^2)
| * \param[out] Shat - modified vorticity (\tilde{S})
| * \param[out] d_Shat - derivative of Shat w.r.t. nue variable
------------------------------------------------------------------------------*/

/*!
 * \brief Baseline
 */
template <class Base>
class ModVort_Bsl : Base {
  using Base::d_fv2;
  using Base::fv2;
  using Base::inv_k2_d2;
  using Base::S;

  static constexpr su2double nue = Base::ScalarVar_i[0];

 public:
  static void get(su2double* Shat, su2double* d_Shat) {
    const su2double Sbar = nue * fv2 * inv_k2_d2;

    *Shat = S + Sbar;
    *Shat = max(*Shat, 1.0e-10);

    const su2double d_Sbar = (fv2 + nue * d_fv2) * inv_k2_d2;

    *d_Shat = (*Shat <= 1.0e-10) ? 0.0 : d_Sbar;
  }
};

/*!
 * \brief Edward
 */
template <class Base>
class ModVort_Edw {
  using Base::d_fv1;
  using Base::fv1;
  using Base::nDim;
  using Base::PrimVar_Grad_i;
  using Base::rotating_frame;
  using Base::StrainMag_i;

  static constexpr su2double xsi = Base::Ji;

  static constexpr su2double nul = Base::Laminar_Viscosity_i;
  static constexpr su2double nue = Base::ScalarVar_i[0];

 public:
  static void get(su2double* Shat, su2double* d_Shat) {
    unsigned short iDim, jDim;

    /*--- Evaluate Omega, here Omega is the Strain Rate ---*/

    su2double Sbar = 0.0;
    for (iDim = 0; iDim < nDim; ++iDim) {
      for (jDim = 0; jDim < nDim; ++jDim) {
        Sbar += (PrimVar_Grad_i[1 + iDim][jDim] + PrimVar_Grad_i[1 + jDim][iDim]) * (PrimVar_Grad_i[1 + iDim][jDim]);
      }
    }
    for (iDim = 0; iDim < nDim; ++iDim) {
      Sbar -= ((2.0 / 3.0) * pow(PrimVar_Grad_i[1 + iDim][iDim], 2.0));
    }

    const su2double Omega = sqrt(max(Sbar, 0.0));

    /*--- Rotational correction term ---*/
    if (rotating_frame) {
      Omega += 2.0 * min(0.0, StrainMag_i - Omega);
    }

    Base::S = Omega;

    *Shat = max(Omega * ((1.0 / max(xsi, 1.0e-16)) + fv1), 1.0e-16);
    *Shat = max(*Shat, 1.0e-10);

    *d_Shat = (*Shat <= 1.0e-10) ? 0.0 : -Omega * pow(xsi, -2.0) / nul + Omega * d_fv1;
  }
};

/*!
 * \brief Negative
 */
template <class Base>
class ModVort_Neg : Base {
  static constexpr su2double nue = Base::ScalarVar_i[0];

  static void get(su2double* Shat, su2double* d_Shat) {
    if (nue > 0.0) {
      // Don't check whether Sbar <>= -cv2*S.
      // Baseline solution
      ModVort_Bsl::get(Shat, d_Shat);
    }
    //    else {
    //      // No need for Sbar
    //    }
  }
};

/*------------------------------------------------------------------------------
| Compute the auxiliary function r and its derivative.
| * \param[in] Shat - modified vorticity (\tilde{S})
| * \param[in] d_Shat - derivative of Shat w.r.t. nue variable
| * \param[in] inv_Shat - inverse of modified vorticity (\tilde{S})
| * \param[in] nue - nu_tilde
| * \param[in] inv_k2_d2 - 1/(k^2*d^2)
| * \param[out] r - auxiliary function r
| * \param[out] dr - derivative of auxiliary function r
------------------------------------------------------------------------------*/

/*!
 * \brief Baseline
 */
template <class Base>
class r_Bsl {
  using Base::d_Shat;
  using Base::inv_k2_d2;
  using Base::inv_Shat;
  using Base::Shat;

  static constexpr su2double nue = Base::ScalarVar_i[0];

 public:
  static double get(su2double* r, su2double* d_r) {
    *r = min(nue * inv_Shat * inv_k2_d2, 10.0);

    *d_r = (Shat - nue * d_Shat) * inv_Shat * inv_Shat * inv_k2_d2;
  }
};

/*!
 * \brief Edward
 */
template <class Base>
class r_Edw {
  using Base::d_Shat;
  using Base::inv_k2_d2;
  using Base::inv_Shat;
  using Base::Shat;

  static constexpr su2double nue = Base::ScalarVar_i[0];

 public:
  static double get(su2double* r, su2double* d_r) {
    *r = min(nue * inv_Shat * inv_k2_d2, 10.0);
    *r = tanh(*r) / tanh(1.0);

    *d_r = (Shat - nue * d_Shat) * inv_Shat * inv_Shat * inv_k2_d2;
    *d_r = (1 - pow(tanh(*r), 2.0)) * (*d_r) / tanh(1.0);
  }
};

/*------------------------------------------------------------------------------
| Compute production term and its derivative.
| * \param[in] S - vorticity (Omega)
| * \param[in] PrimVar_Grad_i - Gradient of primitive variables at point i
| * \param[in] nul - dynamic laminar viscosity
| * \param[in] nue - nu_tilde
| * \param[in] fv1 - auxiliary function fv1
| * \param[in] d_fv1 - derivative of fv1 w.r.t. nue variable
| * \param[in] fv2 - auxiliary function fv2
| * \param[in] d_fv2 - derivative of fv2 w.r.t. nue variable
| * \param[in] inv_k2_d2 - 1/(k^2*d^2)
| * \param[out] Shat - modified vorticity (\tilde{S})
| * \param[out] d_Shat - derivative of Shat w.r.t. nue variable
------------------------------------------------------------------------------*/

/*!
 * \brief Baseline
 */
template <class Base>
class Production_Bsl {
  using Base::cb1;
  using Base::ct3;
  using Base::d_ft2;
  using Base::d_Shat;
  using Base::ft2;
  using Base::S;
  using Base::Shat;

  static constexpr su2double nue = Base::ScalarVar_i[0];

 public:
  static double get(su2double* Pr, su2double* d_Pr) {
    *Pr = cb1 * (1. - ft2) * Shat * nue;

    *d_Pr = cb1 * ((1. - ft2) * (Shat + nue * d_Shat) - d_ft2);
  }
};

/*!
 * \brief Negative
 */
template <class Base>
class Production_Neg {
  using Base::cb1;
  using Base::ct3;
  using Base::S;

  static constexpr su2double nue = Base::ScalarVar_i[0];

 public:
  static double get(su2double* Pr, su2double* d_Pr) {
    if (nue > 0.0) {
      // Baseline solution
      Production_Bsl::get(Pr, d_Pr);
    } else {
      *Pr = cb1 * (1.0 - ct3) * S * nue;

      *d_Pr = cb1 * (1.0 - ct3) * S;
    }
  }
};

/*------------------------------------------------------------------------------
| Compute destruction term and its derivative.
| * \param[in] S - vorticity (Omega)
| * \param[in] PrimVar_Grad_i - Gradient of primitive variables at point i
| * \param[in] nul - dynamic laminar viscosity
| * \param[in] nue - nu_tilde
| * \param[in] fv1 - auxiliary function fv1
| * \param[in] d_fv1 - derivative of fv1 w.r.t. nue variable
| * \param[in] fv2 - auxiliary function fv2
| * \param[in] d_fv2 - derivative of fv2 w.r.t. nue variable
| * \param[in] inv_k2_d2 - 1/(k^2*d^2)
| * \param[out] Shat - modified vorticity (\tilde{S})
| * \param[out] d_Shat - derivative of Shat w.r.t. nue variable
------------------------------------------------------------------------------*/

/*!
 * \brief Baseline
 */
template <class Base>
class Desctruction_Bsl {
  using Base::cb1;
  using Base::cw1;
  using Base::d_ft2;
  using Base::d_fw;
  using Base::d_Shat;
  using Base::dist_i_2;
  using Base::ft2;
  using Base::fw;
  using Base::k2;
  using Base::Shat;

  static constexpr su2double nue = Base::ScalarVar_i[0];

 public:
  static double get(su2double* De, su2double* d_De) {
    *De = (cw1 * fw - cb1 * ft2 / k2) * nue * nue / dist_i_2;

    *d_De = cw1 * (d_fw * nue + 2.0 * fw) * nue / dist_i_2;
  }
};

/*!
 * \brief Negative
 */
template <class Base>
class Destruction_Neg {
  using Base::cb1;
  using Base::cw1;
  using Base::d_ft2;
  using Base::d_fw;
  using Base::d_Shat;
  using Base::dist_i_2;
  using Base::ft2;
  using Base::fw;
  using Base::k2;
  using Base::Shat;

  static constexpr su2double nue = Base::ScalarVar_i[0];

 public:
  static double get(su2double* De, su2double* d_De) {
    if (nue > 0.0) {
      // Baseline solution
      Destruction_Bsl::get(De, d_De);
    } else {
      *De = -cw1 * nue * nue / dist_i_2;

      *d_De = -2.0 * cw1 * nue / dist_i_2;
    }
  }
};

/*------------------------------------------------------------------------------
| Compute cross production term and its derivative.
| * \param[in] S - vorticity (Omega)
| * \param[in] PrimVar_Grad_i - Gradient of primitive variables at point i
| * \param[in] nul - dynamic laminar viscosity
| * \param[in] nue - nu_tilde
| * \param[in] fv1 - auxiliary function fv1
| * \param[in] d_fv1 - derivative of fv1 w.r.t. nue variable
| * \param[in] fv2 - auxiliary function fv2
| * \param[in] d_fv2 - derivative of fv2 w.r.t. nue variable
| * \param[in] inv_k2_d2 - 1/(k^2*d^2)
| * \param[out] Shat - modified vorticity (\tilde{S})
| * \param[out] d_Shat - derivative of Shat w.r.t. nue variable
------------------------------------------------------------------------------*/

/*!
 * \brief Baseline
 */
template <class Base>
class CrossProduction_Bsl {
  using Base::cb2_sigma;
  using Base::norm2_Grad;

  static constexpr su2double nue = Base::ScalarVar_i[0];

  static double get(su2double* CrossProd, su2double* d_CrossProd) {
    *CrossProd = cb2_sigma * norm2_Grad;

    // No cross production influence in the Jacobian
    *d_CrossProd = 0.0;
  }
};
