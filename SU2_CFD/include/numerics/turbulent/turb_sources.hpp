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
  /*--- List of constants and auxiliary functions ---*/
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

  su2double ComputeXsi(const su2double nue, const su2double nul) { return nue/nul  + cr1*(roughness_i/(dist_i+EPS)); //roughness_i = 0 for smooth walls and Ji remains the same, changes only if roughness is specified.nue/nul }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) final {
    
  };
};

/*------------------------------------------------------------------------------
| Define the Spalart-Allmaras variation/corrections.
| \author A. Bueno., E.Molina, F. Palacios
------------------------------------------------------------------------------*/
/*!
 * \brief SU2 baseline ft2 term value
 * \param[in] ft2 - ft2 variable
 * \param[in] d_ft2 - derivative of ft2 w.r.t. nue variable
 * \return Zero
 */
class ft2_bsl { static double get(su2double ft2, su2double d_ft2) { ft2 = 0.0; d_ft2 = 0.0; } };

/*!
 * \brief non-zero ft2 term according to the literature.
 * \param[in] ct3 - constant ct3
 * \param[in] ct4 - constant ct4
 * \param[in] nut - nu_tilde
 * \param[in] nul - dynamic laminar viscosity
 * \param[in] ft2 - ft2 variable
 * \param[in] d_ft2 - derivative of ft2 w.r.t. nue variable
 */
class ft2_nonzero {
  static void get(const su2double ct3, const su2double ct4, const su2double nue, const su2double nul, su2double ft2, su2double d_ft2) {
    const su2double xsi = ComputeXsi(nue, nul),
                    xsi2 = xsi*xsi;

    const su2double dxsi = 1.0/nul;
    
    ft2 = ct3*exp(-ct4*xsi2);
    d_ft2 = -2.0*ct4*xsi*ft2*dxsi;
  }
};

class ModVort_bsl {
  static double get(const su2double S, const su2double nue, const su2double fv2, const su2double inv_k2_d2, const su2double d_fv2, su2double Shat, su2double d_Shat){

    Shat = S + nue*fv2*inv_k2_d2;
    Shat = max(Shat, 1.0e-10);

    d_Shat = (Shat <= 1.0e-10) ? 0.0 : (fv2 + nue*d_fv2)*inv_k2_d2;
  }
};


class ModVort_Edw {
  static double get(const su2double nue, const su2double nul, const su2double fv1, const su2double inv_k2_d2, const su2double d_fv1, su2double Shat, su2double d_Shat){

    /*
      From NASA Turbulence model site. http://turbmodels.larc.nasa.gov/spalart.html
      This form was developed primarily to improve the near-wall numerical behavior of the model (i.e., the goal was to improve the convergence behavior). The reference is:
      Edwards, J. R. and Chandra, S. "Comparison of Eddy Viscosity-Transport Turbulence Models for Three-Dimensional, Shock-Separated Flowfields," AIAA Journal, Vol. 34, No. 4, 1996, pp. 756-763.
      In this modificaton Omega is replaced by Strain Rate
    */

    const su2double xsi = ComputeXsi(nue, nul);

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

    const su2double S = Omega;

    Shat = max(S*((1.0/max(xsi,1.0e-16))+fv1),1.0e-16);
    Shat = max(Shat, 1.0e-10);

    d_Shat = (Shat <= 1.0e-10) ? 0.0 : -S*pow(xsi,-2.0)/nul + S*d_fv1;
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
