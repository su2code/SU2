/*!
 * \file CTransLMVariable.hpp
 * \brief Declaration of the variables of the transition model.
 * \author F. Palacios, T. Economon, A. Rausa, M. Cerabona
 * \version 7.4.0 "Blackbird"
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

//
// Created by marcocera on 25/02/22.
//

#pragma once

#include "CTurbVariable.hpp"

class CTransLMVariable final : public CTurbVariable {

 protected:  
  VectorType F_length, /*!< \brief Transition length function from correlation. */
             F_onset,  /*!< \brief Transition onset function. */
             F_turb,   /*!< \brief Function for disabling destruction term of gamma. */
             F_length1,  /*!< \brief Preliminary F length function. */
             F_sublayer, /*!< \brief Sublayer detection function. */
             F_thetat,   /*!< \brief Blending function. */
             F_lambda;   /*!< \brief Empirical correlation function based on the Thwaites’ pressure gradient prameter. */

  VectorType reV,      /*!< \brief Vorticity Reynolds number. */
             rew,      /*!< \brief Dissipation Reynolds number. */
             rethetac, /*!< \brief Critical Reynolds number. */
             rethetat_eq;  /*!< \brief Predicted transition Reynolds number. */

  VectorType R_T; /*!< \brief Turbulent-to-laminar viscosity ratio. */

  VectorType T;   /*!< \brief Time scale. */

  VectorType delta_param,  /*!< \brief Boundary-layer thickness. */
             F_wake,       /*!< \brief Wake detection function. */
             lambda_theta, /*!< \brief Thwaites’ pressure gradient prameter. */
             Turb_Intens,  /*!< \brief Turbulence intensity. */
             dU_ds,        /*!< \brief Streamline derivative of the low velocity. */
             thetat;       /*!< \brief Transition momentum thickness. */
  
  VectorType gamma_sep,   /*!< \brief Separation induced intermittency. */
             gamma_eff,   /*!< \brief Effective intermittency. */
             FReattach;   /*!< \brief Flow reattachment function. */


  // Variables for LM2015 model
  VectorType F_thetat2,     /*!< \brief Blending function for LM2015 model. */
             ReThetat_SCF,  /*!< \brief Cross-flow induced transition reynolds number. */
             thetat_SCF;    /*!< \brief Cross-flow transition momentum thickness. */

  VectorType k, w;

  bool TurbKW, TurbSA;
  TURB_TRANS_MODEL TransModel;


  TURB_TRANS_CORRELATION Trans_Correlation;

 private:

  void ReThetaC_Correlations(unsigned long iPoint);
  void FLength_Correlations(unsigned long iPoint);

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] intermittency - Intermittency (gamma) (initialization value).
   * \param[in] Re_theta - theta Reynolds number (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] constants   // for now it is off
   * \param[in] config - Definition of the particular problem.
   */
  CTransLMVariable(su2double intermittency, su2double Re_theta, const su2double* constants, unsigned long npoint,
                   unsigned long ndim, unsigned long nvar, /*const su2double* constants*/ CConfig* config);

  /*!
     * \brief Destructor of the class.
   */
  ~CTransLMVariable() override = default;



  void Set_kAndw(unsigned long iPoint, su2double VorticityMag, su2double val_velocity, su2double val_eddy_viscosity, su2double StrainMag, su2double val_dist,
                 su2double val_viscosity, CConfig *config) override;

  void SetF_onset(unsigned long iPoint, su2double val_density, su2double StrainMag, su2double val_dist,
                  su2double val_viscosity, su2double *TurbVars, su2double val_eddy_viscosity, su2double VorticityMag, su2double val_velocity, CConfig *config) override;
  void SetF_length(su2double *Velocity, unsigned long iPoint, su2double val_density, su2double *TurbVars, su2double val_dist,
                   su2double val_viscosity, CConfig *config) override;
  void SetF_turb(unsigned long iPoint) override;
  void Setrethetat_eq(unsigned long iPoint, su2double *Velocity, su2double val_velocity, CMatrixView<const su2double> Velocity_Gradient,
                      su2double *TurbVars, su2double val_viscosity, su2double val_density, CConfig *config) override;
  void SetT(unsigned long iPoint, su2double val_velocity, su2double val_viscosity, su2double val_density, su2double val_localGridLength, su2double val_mu_T) override;
  void SetF_thetat(unsigned long iPoint, su2double val_density, su2double *TurbVars, su2double val_dist,
                   su2double val_viscosity, su2double val_vort, su2double val_velocity, su2double* constants, CConfig *config) override;
  void SetF_thetat_2(unsigned long iPoint, su2double val_dist) override;
  void SetReThetat_SCF(unsigned long iPoint, su2double val_dist, su2double val_density, su2double* val_velocity, su2double val_velocityMag, su2double hRoughness, su2double* val_vorticity, su2double val_viscosity, su2double val_eddy_viscosity) override;
  void Setgamma_sep(unsigned long iPoint, su2double val_density, su2double *TurbVars, su2double val_viscosity,
                    su2double val_dist, su2double StrainMag, su2double val_vort, su2double* constants) override;
  
  /*!
   * \brief Get the various quantities.
   */
  inline su2double Getrew(unsigned long iPoint) const override { return rew(iPoint); }
  inline su2double GetF_length1(unsigned long iPoint) const override { return F_length1(iPoint); }
  inline su2double GetF_sublayer(unsigned long iPoint) const override { return F_sublayer(iPoint); }
  inline su2double GetF_length(unsigned long iPoint) const override { return F_length(iPoint); }
  inline su2double GetreV(unsigned long iPoint) const override { return reV(iPoint); }
  inline su2double GetR_T(unsigned long iPoint) const override { return R_T(iPoint); }
  inline su2double GetF_onset(unsigned long iPoint) const override { return F_onset(iPoint); }
  inline su2double GetF_turb(unsigned long iPoint) const override { return F_turb(iPoint); }
  inline su2double GetT(unsigned long iPoint) const override { return T(iPoint); }
  inline su2double GetdU_ds(unsigned long iPoint) const override { return dU_ds(iPoint); }
  inline su2double GetTurb_Intens(unsigned long iPoint) const override { return Turb_Intens(iPoint); }
  inline su2double GetF_lambda(unsigned long iPoint) const override { return F_lambda(iPoint); }
  inline su2double Getrethetat_eq(unsigned long iPoint) const override { return rethetat_eq(iPoint); }
  inline su2double Getrethetac(unsigned long iPoint) const override { return rethetac(iPoint); }
  inline su2double Getthetat(unsigned long iPoint) const override { return thetat(iPoint); }
  inline su2double Getlambda_theta(unsigned long iPoint) const override { return lambda_theta(iPoint); }
  inline su2double Getdelta_param(unsigned long iPoint) const override { return delta_param(iPoint); }
  inline su2double GetF_wake(unsigned long iPoint) const override { return F_wake(iPoint); }
  inline su2double GetF_thetat(unsigned long iPoint) const override { return F_thetat(iPoint); }
  inline su2double GetF_Reattach(unsigned long iPoint) const override { return FReattach(iPoint); }
  inline su2double GetReThetat_SCF(unsigned long iPoint) const override { return ReThetat_SCF(iPoint); }
  inline su2double GetF_thetat_2(unsigned long iPoint) const override { return F_thetat2(iPoint); } 



  /*!
   * \brief Compute the correction for separation-induced transition.
   */
  inline void SetGammaSep(unsigned long iPoint, su2double val_gamma_sep) override {
    gamma_sep(iPoint) = val_gamma_sep;
  }

  /*!
   * \brief Correction for separation-induced transition.
   */
  inline void SetGammaEff(unsigned long iPoint) override {
    gamma_eff(iPoint) = max(Solution(iPoint, 0), gamma_sep(iPoint));
  }

  /*!
   * \brief Substitute gamma solution with effective gamma
   */
  inline void CorrectGamma(unsigned long iPoint) override {
    Solution(iPoint, 0) = gamma_eff(iPoint);
  }

  /*!
   * \brief Get intermittency value
   */
  inline su2double GetIntermittency(unsigned long iPoint) const override { return gamma_eff(iPoint); }

  /*!
   * \brief Get separation induced intermittency value
   */
  inline su2double GetGammaSep(unsigned long iPoint) const override { return gamma_sep(iPoint);}
  
  /*!
   * \brief Get effective intermittency value
   */
  inline su2double GetGammaEff(unsigned long iPoint) const override { return gamma_eff(iPoint);}
};