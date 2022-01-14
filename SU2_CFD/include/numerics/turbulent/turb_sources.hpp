/*!
 * \file turb_sources.hpp
 * \brief Delarations of numerics classes for integration of source
 *        terms in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.2.0 "Blackbird"
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
 * \brief The variables that are subject to change in each variation/correction have their own class. Additional source terms are implemented as decorators.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 */
template<class ft2_class, class ModVort_class, class r_class>
class CSourceBase_TurbSA : public CNumerics {
protected:
  /*--- List of constants ---*/
  su2double cv1_3,
            k2,
            cb1,
            cw2,
            ct3,
            ct4,
            cw3_6,
            cb2_sigma,
            sigma,
            cb2,
            cw1,
            cr1;

  /*--- List of auxiliary functions ---*/
  su2double ft2, d_ft2,
            r, d_r,
            g, d_g, g_6, glim,
            fw,
            Ji, Ji_2, Ji_3, d_Ji,
            S, Omega, Shat, d_Shat, inv_Shat,
            fv1, fv2;

  su2double Gamma_BC = 0.0;
  su2double intermittency;

  /*--- Source term components ---*/
  su2double Production, Destruction, CrossProduction, AddSourceTerm;

  /*--- Residual and Jacobian ---*/
  su2double Residual, *Jacobian_i;
  
private:
  su2double Jacobian_Buffer; /// Static storage for the Jacobian (which needs to be pointer for return type).

protected:
  const bool incompressible = false, rotating_frame = false;
  bool roughwall = false;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceBase_TurbSA(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

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

  inline su2double ComputeXsi(const su2double nue, const su2double nul, su2double xsi, su2double d_xsi) {
    xsi = nue/nul  + cr1*(roughness_i/(dist_i+EPS)); //roughness_i = 0 for smooth walls and Ji remains the same, changes only if roughness is specified.nue/nul
    d_xsi = 1.0/nul;
  }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) final {

     // Set the boolean here depending on whether the point is closest to a rough wall or not.
     roughwall = (roughness_i > 0.0);

     if (incompressible) {
       Density_i = V_i[nDim+2];
       Laminar_Viscosity_i = V_i[nDim+4];
     }
     else {
       Density_i = V_i[nDim+2];
       Laminar_Viscosity_i = V_i[nDim+5];
     }

     Residual        = 0.0;
     Production      = 0.0;
     Destruction     = 0.0;
     CrossProduction = 0.0;
     Jacobian_i[0]   = 0.0;

     /*--- Evaluate Omega ---*/

     Omega = sqrt(Vorticity_i[0]*Vorticity_i[0] + Vorticity_i[1]*Vorticity_i[1] + Vorticity_i[2]*Vorticity_i[2]);

     /*--- Rotational correction term ---*/

     if (rotating_frame) { Omega += 2.0*min(0.0, StrainMag_i-Omega); }

     if (dist_i > 1e-10) {

       /*--- Production term ---*/

       dist_i_2 = dist_i*dist_i;
       nu = Laminar_Viscosity_i/Density_i;

       /*--- Modified values for roughness ---*/
       /*--- Ref: Aupoix, B. and Spalart, P. R., "Extensions of the Spalart-Allmaras Turbulence Model to Account for Wall Roughness,"
        * International Journal of Heat and Fluid Flow, Vol. 24, 2003, pp. 454-462. ---*/
       /* --- See https://turbmodels.larc.nasa.gov/spalart.html#sarough for detailed explanation. ---*/

       ComputeXsi(ScalarVar_i[0], Laminar_Viscosity_i, &Ji, &d_Ji);
       Ji_2 = Ji*Ji;
       Ji_3 = Ji_2*Ji;
       fv1 = Ji_3/(Ji_3+cv1_3);

       /*--- Using a modified relation so as to not change the Shat that depends on fv2. ---*/
       fv2 = 1.0 - ScalarVar_i[0]/(nu+ScalarVar_i[0]*fv1);   // From NASA turb modeling resource and 2003 paper

       ft2 = ct3*exp(-ct4*Ji_2);
       S = Omega;
       inv_k2_d2 = 1.0/(k2*dist_i_2);

       Shat = S + ScalarVar_i[0]*fv2*inv_k2_d2;
       Shat = max(Shat, 1.0e-10);
       inv_Shat = 1.0/Shat;

   //    Original SA model
   //    Production = cb1*(1.0-ft2)*Shat*ScalarVar_i[0]*Volume;

       if (transition) {

         /*--- BC model constants (2020 revision). ---*/
         const su2double chi_1 = 0.002;
         const su2double chi_2 = 50.0;

         /*--- turbulence intensity is u'/U so we multiply by 100 to get percentage ---*/
         su2double tu = 100.0 * config->GetTurbulenceIntensity_FreeStream();

         su2double nu_t = (ScalarVar_i[0]*fv1); //S-A variable

         su2double re_v = ((Density_i*pow(dist_i,2.))/(Laminar_Viscosity_i))*Omega;
         su2double re_theta = re_v/2.193;
         su2double re_theta_t = (803.73 * pow((tu + 0.6067),-1.027)); //MENTER correlation
         //re_theta_t = 163.0 + exp(6.91-tu); //ABU-GHANNAM & SHAW correlation

         su2double term1 = sqrt(max(re_theta-re_theta_t,0.)/(chi_1*re_theta_t));
         su2double term2 = sqrt(max((nu_t*chi_2)/nu,0.));
         su2double term_exponential = (term1 + term2);

         Gamma_BC = 1.0 - exp(-term_exponential);

         Production = Gamma_BC*cb1*Shat*ScalarVar_i[0]*Volume;
       }
       else {
         Production = cb1*Shat*ScalarVar_i[0]*Volume;
       }

       /*--- Destruction term ---*/

       r = min(ScalarVar_i[0]*inv_Shat*inv_k2_d2,10.0);
       g = r + cw2*(pow(r,6.0)-r);
       g_6 =  pow(g,6.0);
       glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
       fw = g*glim;

   //    Original SA model
   //    Destruction = (cw1*fw-cb1*ft2/k2)*ScalarVar_i[0]*ScalarVar_i[0]/dist_i_2*Volume;

       Destruction = cw1*fw*ScalarVar_i[0]*ScalarVar_i[0]/dist_i_2*Volume;

       /*--- Diffusion term ---*/

       norm2_Grad = 0.0;
       for (iDim = 0; iDim < nDim; iDim++)
         norm2_Grad += ScalarVar_Grad_i[0][iDim]*ScalarVar_Grad_i[0][iDim];

       CrossProduction = cb2_sigma*norm2_Grad*Volume;

       Residual = Production - Destruction + CrossProduction;

       /*--- Implicit part, production term ---*/

       dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
       dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
       if ( Shat <= 1.0e-10 ) dShat = 0.0;
       else dShat = (fv2+ScalarVar_i[0]*dfv2)*inv_k2_d2;

       if (transition) {
         Jacobian_i[0] += Gamma_BC*cb1*(ScalarVar_i[0]*dShat+Shat)*Volume;
       }
       else {
         Jacobian_i[0] += cb1*(ScalarVar_i[0]*dShat+Shat)*Volume;
       }

       /*--- Implicit part, destruction term ---*/

       dr = (Shat-ScalarVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
       if (r == 10.0) dr = 0.0;
       dg = dr*(1.+cw2*(6.0*pow(r,5.0)-1.0));
       dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
       Jacobian_i[0] -= cw1*(dfw*ScalarVar_i[0] +  2.0*fw)*ScalarVar_i[0]/dist_i_2*Volume;

     ft2_class::get(&ft2, &d_ft2);

     su2double pingu = ft2*3.14;

     ModVort_class::get();
    
  };
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
template<class Base>
class ft2_Bsl : public Base  {
  static double get(su2double* ft2, su2double* d_ft2) {
    *ft2 = 0.0;
    *d_ft2 = 0.0;
  }
};

/*!
 * \brief non-zero ft2 term according to the literature and its derivative.
 */
template<class Base>
class ft2_nonzero : public Base {
    using Base::ct3;
    using Base::ct4;
    using Base::xsi;
    using Base::d_xsi;

    static void get(su2double* ft2, su2double* d_ft2) {

        const su2double xsi2 = xsi*xsi;

        *ft2 = ct3*exp(-ct4*xsi2);
        *d_ft2 = -2.0*ct4*xsi**ft2*d_xsi;
      }
};

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
template<class Base>
class ModVort_Bsl : Base {
    using Base::S;
    using Base::fv2;
    using Base::d_fv2;
    using Base::inv_k2_d2;

    static constexpr su2double nue = Base::ScalarVar_i[0];

public:
  static void get(su2double* Shat, su2double* d_Shat) {

    const su2double Sbar = nue*fv2*inv_k2_d2;

    *Shat = S + Sbar;
    *Shat = max(*Shat, 1.0e-10);

    const su2double d_Sbar = (fv2 + nue*d_fv2)*inv_k2_d2;
    
    *d_Shat = (*Shat <= 1.0e-10) ? 0.0 : d_Sbar;
  }
};

/*!
 * \brief Edward
 */
template<class Base>
class ModVort_Edw {
    using Base::nDim;
    using Base::PrimVar_Grad_i;
    using Base::fv1;
    using Base::d_fv1;
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
       for(iDim=0;iDim<nDim;++iDim){
         for(jDim=0;jDim<nDim;++jDim){
       Sbar+= (PrimVar_Grad_i[1+iDim][jDim]+PrimVar_Grad_i[1+jDim][iDim])*(PrimVar_Grad_i[1+iDim][jDim]);}}
       for(iDim=0;iDim<nDim;++iDim){
         Sbar-= ((2.0/3.0)*pow(PrimVar_Grad_i[1+iDim][iDim],2.0));}

       const su2double Omega = sqrt(max(Sbar,0.0));

       /*--- Rotational correction term ---*/
       if (rotating_frame) { Omega += 2.0*min(0.0, StrainMag_i-Omega); }

       Base::S = Omega;

       *Shat = max(Omega*((1.0/max(xsi,1.0e-16))+fv1),1.0e-16);
       *Shat = max(*Shat, 1.0e-10);

       *d_Shat = (*Shat <= 1.0e-10) ? 0.0 : -Omega*pow(xsi,-2.0)/nul + Omega*d_fv1;
    }
};

/*!
 * \brief Negative
 */
template<class Base>
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
template<class Base>
class r_Bsl {
   using Base::Shat;
   using Base::d_Shat;
   using Base::inv_Shat;
   using Base::inv_k2_d2;

   static constexpr su2double nue = Base::ScalarVar_i[0];

public:
   static double get(su2double* r, su2double* d_r){

       *r = min(nue*inv_Shat*inv_k2_d2,10.0);

       *d_r = (Shat-nue*d_Shat)*inv_Shat*inv_Shat*inv_k2_d2;
   }
};

/*!
 * \brief Edward
 */
template<class Base>
class r_Edw {
   using Base::Shat;
   using Base::d_Shat;
   using Base::inv_Shat;
   using Base::inv_k2_d2;

   static constexpr su2double nue = Base::ScalarVar_i[0];

public:
   static double get(su2double* r, su2double* d_r){

       *r = min(nue*inv_Shat*inv_k2_d2,10.0);
       *r = tanh(*r)/tanh(1.0);

       *d_r = (Shat-nue*d_Shat)*inv_Shat*inv_Shat*inv_k2_d2;
       *d_r = (1-pow(tanh(*r),2.0))*(*d_r)/tanh(1.0);
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
template<class Base>
class Production_Bsl{
   using Base::S;
   using Base::ft2;
   using Base::d_ft2;
   using Base::Shat;
   using Base::d_Shat;
   using Base::ct3;
   using Base::cb1;

   static constexpr su2double nue = Base::ScalarVar_i[0];

public:
  static double get(su2double* Pr, su2double* d_Pr){

    *Pr = cb1*(1. - ft2)*Shat*nue;

    *d_Pr = cb1*((1. - ft2)*(Shat + nue*d_Shat) - d_ft2);
  }
};

/*!
 * \brief Negative
 */
template<class Base>
class Production_Neg{
   using Base::S;
   using Base::ct3;
   using Base::cb1;

   static constexpr su2double nue = Base::ScalarVar_i[0];

public:
   static double get(su2double* Pr, su2double* d_Pr){

      if (nue > 0.0) {
         // Baseline solution
         Production_Bsl::get(Pr, d_Pr);
      }
      else {
         *Pr = cb1*(1.0-ct3)*S*nue;

         *d_Pr = cb1*(1.0-ct3)*S;
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
template<class Base>
class Desctruction_Bsl{
   using Base::Shat;
   using Base::d_Shat;
   using Base::fw;
   using Base::d_fw;
   using Base::ft2;
   using Base::d_ft2;
   using Base::k2;
   using Base::dist_i_2;
   using Base::cw1;
   using Base::cb1;

   static constexpr su2double nue = Base::ScalarVar_i[0];

public:
   static double get(su2double* De, su2double* d_De){

      *De = (cw1*fw-cb1*ft2/k2)*nue*nue/dist_i_2;

      *d_De = cw1*(d_fw*nue +  2.0*fw)*nue/dist_i_2;
   }
};

/*!
 * \brief Negative
 */
template<class Base>
class Destruction_Neg{
   using Base::Shat;
   using Base::d_Shat;
   using Base::fw;
   using Base::d_fw;
   using Base::ft2;
   using Base::d_ft2;
   using Base::k2;
   using Base::dist_i_2;
   using Base::cw1;
   using Base::cb1;

   static constexpr su2double nue = Base::ScalarVar_i[0];

public:
   static double get(su2double* De, su2double* d_De){

      if (nue > 0.0) {
         // Baseline solution
         Destruction_Bsl::get(De, d_De);
      }
      else {
         *De = -cw1*nue*nue/dist_i_2;

         *d_De = -2.0*cw1*nue/dist_i_2;
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
template<class Base>
class CrossProduction_Bsl{
   using Base::norm2_Grad;
   using Base::cb2_sigma;

   static constexpr su2double nue = Base::ScalarVar_i[0];

   static double get(su2double* CrossProd, su2double* d_CrossProd){

      *CrossProd = cb2_sigma*norm2_Grad;

      // No cross production influence in the Jacobian
      *d_CrossProd = 0.0;
   }
};





























  
/*!
 * \class CSourcePieceWise_TurbSST
 * \brief Class for integrating the source terms of the Menter SST turbulence model equations.
 * \ingroup SourceDiscr
 * \author A. Campos.
 */
class CSourcePieceWise_TurbSST final : public CNumerics {
private:
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

  su2double Residual[2],
  *Jacobian_i[2] = {nullptr},
  Jacobian_Buffer[4] = {0.0}; /// Static storage for the Jacobian (which needs to be pointer for return type).

  bool incompressible;
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
  inline void ResidualAxisymmetric(su2double alfa_blended, su2double zeta){

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
    pk_axi = max(0.0,2.0/3.0*rhov*k*(2.0/zeta*(yinv*V_i[2]-PrimVar_Grad_i[2][1]-PrimVar_Grad_i[1][0])-1.0));
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
