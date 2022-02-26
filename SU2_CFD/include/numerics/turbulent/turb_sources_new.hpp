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
 * \class CommonSAVariables
 * \brief Structure with SA common auxiliary functions and constants.
 */
struct CommonSAVariables {

   /*--- List of constants ---*/
   const su2double
     cv1_3 = pow(7.1, 3.0),
     k2    = pow(0.41, 2.0),
     cb1   = 0.1355,
     cw2   = 0.3,
     ct3   = 1.2,
     ct4   = 0.5,
     cw3_6 = pow(2.0, 6.0),
     sigma = 2./3.,
     cb2   = 0.622,
     cb2_sigma = cb2/sigma,
     cw1 = cb1/k2+(1.0+cb2)/sigma,
     cr1 = 0.5;

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
template <class FlowIndices, class Omega_class, class ft2_class, class ModVort_class, class r_class, class SourceTerms_class>
class CSourceBase_TurbSA : public CNumerics {
 protected:

  su2double Gamma_BC = 0.0;
  su2double intermittency;

  /*--- Source term components ---*/
  su2double Production, Destruction, CrossProduction, AddSourceTerm;

  /*--- Residual and Jacobian ---*/
  su2double Residual, *Jacobian_i;

 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */
  su2double Jacobian_Buffer; /*!< \brief Static storage for the Jacobian (which needs to be pointer for return type). */

 protected:
  const bool rotating_frame = false;
  const bool transition = false;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceBase_TurbSA(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
      : CNumerics(val_nDim, val_nVar, config), idx(val_nDim, config->GetnSpecies()),
        rotating_frame(config->GetRotating_Frame()), transition(config->GetKind_Trans_Model() == TURB_TRANS_MODEL::BC) {}

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
    Production = 0.0;
    Destruction = 0.0;
    CrossProduction = 0.0;
    AddSourceTerm = 0.0;
    Jacobian_i[0] = 0.0;

    /*--- Evaluate Omega ---*/
    if (nDim == 2) {
      Omega_class::get<2>(Vorticity_i, PrimVar_Grad_i, idx.Velocity(), model_var);
    } else {
      Omega_class::get<3>(Vorticity_i, PrimVar_Grad_i, idx.Velocity(), model_var);
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
      model_var.Ji = ScalarVar_i[0] / nu + model_var.cr1 * (roughness_i / (dist_i + EPS));  // roughness_i = 0 for smooth walls and Ji remains the
                                                                                            // same, changes only if roughness is specified.nue/nul
      model_var.d_Ji = 1.0 / nu;

      const su2double Ji_2 = model_var.Ji * model_var.Ji;
      const su2double Ji_3 = Ji_2 * model_var.Ji;

      model_var.fv1 = Ji_3 / (Ji_3 + model_var.cv1_3);
      model_var.d_fv1 = 3.0 * Ji_2 * model_var.cv1_3 / (nu * pow(Ji_3 + model_var.cv1_3, 2.0));

      /*--- Using a modified relation so as to not change the Shat that depends on fv2. ---*/
      model_var.fv2 = 1.0 - ScalarVar_i[0] / (nu + ScalarVar_i[0] * model_var.fv1);  // From NASA turb modeling resource and 2003 paper
      model_var.d_fv2 = -(1.0 / nu - Ji_2 * model_var.d_fv1) / pow(1.0 + model_var.Ji * model_var.fv1, 2.0);

      /*--- Compute ft2 term ---*/
      ft2_class::get(model_var);

      /*--- Compute modified vorticity ---*/
      ModVort_class::get(ScalarVar_i[0], nu, model_var);

      model_var.inv_Shat = 1.0 / model_var.Shat;

      /*--- Compute auxiliary function r ---*/
      r_class::get(ScalarVar_i[0], model_var);

      model_var.g = model_var.r + model_var.cw2 * (pow(model_var.r, 6.0) - model_var.r);
      model_var.g_6 = pow(model_var.g, 6.0);
      model_var.glim = pow((1.0 + model_var.cw3_6) / (model_var.g_6 + model_var.cw3_6), 1.0 / 6.0);
      model_var.fw = model_var.g * model_var.glim;

      model_var.d_g  = model_var.d_r * (1. + model_var.cw2 * (6.0 * pow(model_var.r, 5.0) - 1.0));
      model_var.d_fw = model_var.d_g * model_var.glim * (1. - model_var.g_6 / (model_var.g_6 + model_var.cw3_6));

      model_var.norm2_Grad = GeometryToolbox::SquaredNorm(nDim, ScalarVar_Grad_i[0]);

      /*--- Compute production, destruction and cross production and jacobian ---*/
      SourceTerms_class::get(model_var, Production, Destruction, CrossProduction, Jacobian_i[0]);

      /*--- Compute any necessary additional source term and jacobian contribution ---*/

      /*--- Residual ---*/
      Residual = Production - Destruction + CrossProduction + AddSourceTerm;
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

/*------------------------------------------------------------------------------
| Strain rate
|
| \tparam nDim: problem dimension
| \param Vorticity_i: vorticity array
| \param PrimVar_Grad_i: primitive variables gradient at point i
| \param Velocity index
| \param model_var: Common SA variables struct
|
------------------------------------------------------------------------------*/

/*!
 * \brief Baseline
 */
struct Omega_Bsl {
   template<int nDim, class MatrixType>
   static void get(const su2double* Vorticity_i, const MatrixType& PrimVar_Grad_i, unsigned short vel_idx, CommonSAVariables &model_var) {
      model_var.Omega = GeometryToolbox::Norm(3, Vorticity_i);
   }
};

/*!
 * \brief Edward
 * Here Omega is the Strain Rate
 */
struct Omega_Edw {
  template<int nDim, class MatrixType>
  static void get(const su2double* Vorticity_i, const MatrixType& PrimVar_Grad_i, const unsigned short vel_idx,
                  CommonSAVariables &model_var) {
    su2double Sbar = 0.0;
    for(int iDim=0; iDim<nDim; ++iDim){
      for(int jDim=0; jDim<nDim; ++jDim){
        Sbar += (PrimVar_Grad_i[vel_idx+iDim][jDim]+
                 PrimVar_Grad_i[vel_idx+jDim][iDim]) * PrimVar_Grad_i[vel_idx+iDim][jDim];
      }
    }
    for(int iDim=0; iDim<nDim; ++iDim){
      Sbar -= (2.0/3.0) * pow(PrimVar_Grad_i[vel_idx+iDim][iDim], 2);
    }

    model_var.Omega = sqrt(max(Sbar, 0.0));
  }
};


/*------------------------------------------------------------------------------
| ft2-term and its derivative
|
| \param model_var: Common SA variables struct
|
------------------------------------------------------------------------------*/

/*!
 * \brief SU2 baseline ft2 term value and its derivative. ft2=0.0
 */
struct ft2_Bsl {
   static void get(CommonSAVariables &model_var) {
      model_var.ft2 = 0.0;
      model_var.d_ft2 = 0.0;
   }
};

/*!
 * \brief non-zero ft2 term according to the literature and its derivative.
 */
struct ft2_nonzero {
   static void get(CommonSAVariables &model_var) {
      const su2double xsi2 = pow(model_var.Ji, 2);
      model_var.ft2  = model_var.ct3 * exp(-model_var.ct4 * xsi2);
      model_var.d_ft2 = -2.0 * model_var.ct4 * model_var.Ji * model_var.ft2 * model_var.d_Ji;
   }
};


/*------------------------------------------------------------------------------
| Modified vorticity (\tilde{S}) and its derivative
|
| \param nue: SA variable
| \param nu: laminar viscosity
| \param model_var: Common SA varibales struct
|
------------------------------------------------------------------------------*/

/*!
 * \brief Baseline
 */
struct ModVort_Bsl {
   static void get(const su2double nue, const su2double nu, CommonSAVariables &model_var) {

      const su2double Sbar = nue * model_var.fv2 * model_var.inv_k2_d2;

      model_var.Shat = model_var.S + Sbar;
      model_var.Shat = max(model_var.Shat, 1.0e-10);

      const su2double d_Sbar = (model_var.fv2 + nue * model_var.d_fv2) * model_var.inv_k2_d2;

      model_var.d_Shat = (model_var.Shat <= 1.0e-10) ? 0.0 : d_Sbar;
   }
};

/*!
 * \brief Edward
 */
struct ModVort_Edw {
   static void get(const su2double nue, const su2double nu, CommonSAVariables &model_var) {

      model_var.Shat = max(model_var.S*((1.0/max(model_var.Ji,1.0e-16))+model_var.fv1),1.0e-16);
      model_var.Shat = max(model_var.Shat, 1.0e-10);

      model_var.d_Shat = (model_var.Shat <= 1.0e-10) ? 0.0 : -model_var.S*pow(model_var.Ji,-2.0)/nu + model_var.S*model_var.d_fv1;
   }
};

/*!
 * \brief Negative
 */
struct ModVort_Neg {
  static void get(const su2double nue, const su2double nu, CommonSAVariables &model_var) {
    if (nue > 0.0) {
       // Baseline solution
       ModVort_Bsl::get(nue, nu, model_var);
    }
    /*--- Don't check whether Sbar <>= -cv2*S.
     * Steven R. Allmaras, Forrester T. Johnson and Philippe R. Spalart -
     * "Modifications and Clarifications for the Implementation of the Spalart-Allmaras Turbulence Model" eq. 12
     * No need for Sbar ---*/
  }
};


/*------------------------------------------------------------------------------
| Auxiliary function r and its derivative.
|
| \param nue: SA variable
|
------------------------------------------------------------------------------*/

/*!
 * \brief Baseline
 */
struct r_Bsl {

   static void get(const su2double nue, CommonSAVariables &model_var) {

      model_var.r = min(nue * model_var.inv_Shat * model_var.inv_k2_d2, 10.0);
      model_var.d_r = (model_var.Shat - nue * model_var.d_Shat) * model_var.inv_Shat * model_var.inv_Shat * model_var.inv_k2_d2;
      if (model_var.r == 10.0) model_var.d_r = 0.0;
   }
};

/*!
 * \brief Edward
 */
struct r_Edw {

   static void get(const su2double nue, CommonSAVariables &model_var) {

      model_var.r = min(nue * model_var.inv_Shat * model_var.inv_k2_d2, 10.0);
      model_var.r = tanh(model_var.r) / tanh(1.0);

      model_var.d_r = (model_var.Shat - nue * model_var.d_Shat) * model_var.inv_Shat * model_var.inv_Shat * model_var.inv_k2_d2;
      model_var.d_r = (1 - pow(tanh(model_var.r), 2.0)) * (model_var.d_r) / tanh(1.0);
   }
};


/*------------------------------------------------------------------------------
| Compute source terms: production, destruction and cross-productions term and their derivatives.
|
| \brief get:
| \param nue: SA variable
|
------------------------------------------------------------------------------*/

/*!
 * \brief Baseline (Original SA model)
 */
struct SourceTerms_Bsl {

   static void get(const su2double nue, const CommonSAVariables &model_var, su2double& Production, su2double& Destruction,
                   su2double& CrossProduction, su2double& Jacobian) {

      ComputeProduction(nue, model_var, Production, Jacobian);

      ComputeDestruction(nue, model_var, Destruction, Jacobian);

      ComputeCrossProduction(nue, model_var, CrossProduction, Jacobian);
   }

   static void ComputeProduction(const su2double nue, const CommonSAVariables &model_var, su2double& Production, su2double& Jacobian) {
      Production = model_var.cb1 * (1.0 - model_var.ft2) * model_var.Shat * nue;
      Jacobian  += model_var.cb1 * (-model_var.Shat * nue * model_var.d_ft2 + (1.0 - model_var.ft2) * (nue * model_var.d_Shat + model_var.Shat));
   }

   static void ComputeDestruction(const su2double nue, const CommonSAVariables &model_var, su2double& Destruction, su2double& Jacobian) {
      Destruction = (model_var.cw1 * model_var.fw - model_var.cb1 * model_var.ft2 / model_var.k2) * nue * nue / model_var.dist_i_2;
      Jacobian   -= (model_var.cw1 * model_var.d_fw - model_var.cb1 / model_var.k2 * model_var.d_ft2) * nue * nue / model_var.dist_i_2 + (model_var.cw1 * model_var.fw - model_var.cb1 * model_var.ft2 / model_var.k2) * 2.0 * nue / model_var.dist_i_2;
   }

   static void ComputeCrossProduction(const su2double nue, const CommonSAVariables &model_var, su2double& CrossProduction, su2double& Jacobian) {
      CrossProduction = model_var.cb2_sigma * model_var.norm2_Grad;
      /*--- No contribution to the jacobian ---*/
   }
};

/*!
 * \brief Negative
 */
struct SourceTerms_Neg {

   static void get(const su2double nue, const CommonSAVariables &model_var, su2double& Production, su2double& Destruction,
                   su2double& CrossProduction, su2double& Jacobian) {
      if (nue > 0.0) {

         // Baseline solution
         SourceTerms_Bsl::get(nue, model_var, Production, Destruction, CrossProduction, Jacobian);

      } else {

         ComputeProduction(nue, model_var, Production, Jacobian);

         ComputeDestruction(nue, model_var, Destruction, Jacobian);

         ComputeCrossProduction(nue, model_var, CrossProduction, Jacobian);

      }
   }

   static void ComputeProduction(const su2double nue, const CommonSAVariables &model_var, su2double& Production, su2double& Jacobian) {
      Production = model_var.cb1 * (1.0 - model_var.ct3) * model_var.S * nue;
      Jacobian  += model_var.cb1 * (1.0 - model_var.ct3) * model_var.S;
   }

   static void ComputeDestruction(const su2double nue, const CommonSAVariables &model_var, su2double& Destruction, su2double& Jacobian) {
      Destruction = model_var.cw1 * nue * nue / model_var.dist_i_2;
      Jacobian   -= 2.0 * model_var.cw1 * nue / model_var.dist_i_2;
   }

   static void ComputeCrossProduction(const su2double nue, const CommonSAVariables &model_var, su2double& CrossProduction, su2double& Jacobian) {
      /*--- Same cross production as baseline. ---*/
      SourceTerms_Bsl::ComputeCrossProduction(nue, model_var, CrossProduction, Jacobian);
   }
};


/* =============================================================================
 * SPALART-ALLMARAS ADDITIONAL SOURCE TERMS DECORATORS
 * ============================================================================*/

/*!
 * \brief Mixing Layer Compressibility Correction (SA-comp) decorator
 */
template <class ParentClass>
class CompressiblityCorrection final : public ParentClass {
 private:
   const su2double c5 = 3.5;

 public:
     template <class... Ts>
     CompressiblityCorrection(const Ts&... args) : ParentClass(args...) {}

     ResidualType<> ComputeResidual(const CConfig* config) override {
        /*--- Residual from standard SA ---*/
        ParentClass::ComputeResidual(config);

        /*--- Compressibility Correction term ---*/
        const auto& Pressure_i = V_i[idx.Pressure()];
        const su2double SoundSpeed_i = sqrt(Pressure_i*Gamma/Density_i);
        su2double aux_cc=0;
        for(iDim=0;iDim<nDim;++iDim){
          for(jDim=0;jDim<nDim;++jDim){
            aux_cc += pow(PrimVar_Grad_i[idx.Velocity()+iDim][jDim], 2);
          }
        }

        const su2double d_CompCorrection = 2.0*c5*ScalarVar_i[0]/pow(SoundSpeed_i, 2)*aux_cc*Volume;
        const su2double CompCorrection = 0.5*ScalarVar_i[0]*d_CompCorrection;

        /*--- Decorator residual contribution ---*/
        Residual.residual[0] -= CompCorrection;
        Residual.jacobian_i[0][0] -= d_CompCorrection;


        return ResidualType<>(&Residual, &Jacobian_i, nullptr);
     }
};
