/*!
 * \file trans_sources.hpp
 * \brief Numerics classes for integration of source terms in transition problems.
 * \author A. Rausa, M. Cerabona
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
#pragma once

#include "../scalar/scalar_sources.hpp"

template <class FlowIndices>
class CSourcePiecewise_TransLM final : public CNumerics {
private:

    const FlowIndices idx;

 // costanti sono nel variable

    su2double gammaAmb,
              rethetaAmb;

    su2double Residual[2],
              *Jacobian_i[2] = {nullptr},
              Jacobian_Buffer[4] = {0.0};

    su2double F_length_i,
              F_onset_i,
              F_thetat_i,
              F_turb_i,
              rethetat_eq_i,
              T_param_i,
              F_wake_i,
              delta_param_i,
              ReThetat_SCF_i,
              F_thetat_2_i;

    su2double F_length_j,
              F_onset_j,
              F_thetat_j,
              F_turb_j,
              rethetat_eq_j,
              T_param_j,
              F_wake_j,
              delta_param_j,
              ReThetat_SCF_j,
              F_thetat_2_j;

    su2double c_a1,
              c_a2,
              c_e1,
              c_e2,
              c_thetat,
              c_CF;



    bool incompressible;
    TURB_TRANS_MODEL TransModel;

    // eventuali parti aggiuntive

public:
    /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   CSourcePiecewise_TransLM(unsigned short val_nDim, unsigned short val_nVar,
                               const su2double *constants, su2double val_gamma_Inf,
                               su2double val_rethetat_Inf, const CConfig* config):
                               CNumerics(val_nDim, val_nVar, config),
                               idx(val_nDim, config->GetnSpecies()),
                               incompressible(config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE),
                                                                               gammaAmb(val_gamma_Inf),
                                                                               rethetaAmb(val_rethetat_Inf),
                                                                               c_a1(constants[0]),
                                                                               c_a2(constants[1]),
                                                                               c_e1(constants[2]),
                                                                               c_e2(constants[3]),
                                                                               c_thetat(constants[4]),
                                                                               c_CF(constants[8]),
                                                                               TransModel(config->GetKind_Trans_Model()) {
     /*--- "Allocate" the Jacobian using the static buffer. ---*/
     Jacobian_i[0] = Jacobian_Buffer;
     Jacobian_i[1] = Jacobian_Buffer + 2;

   }


    /*!
    * \brief Residual for source term integration.
    * \param[in] config - Definition of the particular problem.
    * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
    */
    ResidualType<> ComputeResidual(const CConfig* config) override{


      unsigned short iDim;
      su2double P_gamma, P_rethetat, D_gamma, D_retheta_t;
      su2double VorticityMag = sqrt(Vorticity_i[0]*Vorticity_i[0] +
                                    Vorticity_i[1]*Vorticity_i[1] +
                                    Vorticity_i[2]*Vorticity_i[2]);


      Density_i = V_i[idx.Density()];
      Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];
      Eddy_Viscosity_i = V_i[idx.EddyViscosity()];



      Residual[0] = 0.0;            Residual[1] = 0.0;
      Jacobian_i[0][0] = 0.0;       Jacobian_i[0][1] = 0.0;
      Jacobian_i[1][0] = 0.0;       Jacobian_i[1][1] = 0.0;


    if (dist_i > 1e-10) {


      /*--- Production ---*/
      P_gamma = F_length_i * c_a1 * Density_i * StrainMag_i * pow(ScalarVar_i[0] * F_onset_i, 0.5) *
                (1 - c_e1 * ScalarVar_i[0]);
      P_rethetat = c_thetat * (Density_i / T_param_i) * (rethetat_eq_i - ScalarVar_i[1]) * (1.0 - F_thetat_i);
      D_gamma = c_a2 * Density_i * VorticityMag * ScalarVar_i[0] * F_turb_i * (c_e2 * ScalarVar_i[0] - 1);

      D_retheta_t = 0.0;
      if( TransModel == TURB_TRANS_MODEL::LM2015) {
        D_retheta_t =
            c_thetat * (Density_i / T_param_i) * c_CF * min(ReThetat_SCF_i - ScalarVar_i[1], 0.0) * F_thetat_2_i;
      }


      /*--- Add the production terms to the residuals ---*/
      Residual[0] += P_gamma * Volume;
      Residual[0] -= D_gamma * Volume;
      Residual[1] += P_rethetat * Volume;
      Residual[1] += D_retheta_t * Volume;



      /*--- Implicit part ---*/

      Jacobian_i[0][0] = (F_length_i*c_a1*StrainMag_i*sqrt(F_onset_i)*(0.5*pow(ScalarVar_i[0], -0.5) -1.5*c_e1*pow(ScalarVar_i[0], 0.5))
                          - c_a2 * VorticityMag*F_turb_i*(2.0*c_e2*ScalarVar_i[0]-1.0) )*Volume;
      Jacobian_i[0][1] = 0.0;
      Jacobian_i[1][0] = 0.0;
      Jacobian_i[1][1] = -(c_thetat/T_param_i)*(1.0-F_thetat_i)*Volume;

      if(TransModel == TURB_TRANS_MODEL::LM2015 && min(ReThetat_SCF_i-ScalarVar_i[1], 0.0) != 0.0)
        Jacobian_i[1][1] = -(c_thetat/T_param_i)*c_CF*F_thetat_2_i*Volume;

    }

      return ResidualType<>(Residual, Jacobian_i, nullptr);

   }

    inline void SetF_length(su2double val_F_length_i, su2double val_F_length_j) override {
      F_length_i = val_F_length_i;
      F_length_j = val_F_length_j;
    }

    inline void SetF_onset(su2double val_F_onset_i, su2double val_F_onset_j) override {
      F_onset_i = val_F_onset_i;
      F_onset_j = val_F_onset_j;
    }

    inline void SetF_thetat(su2double val_F_thetat_i, su2double val_F_thetat_j) override {
      F_thetat_i = val_F_thetat_i;
      F_thetat_j = val_F_thetat_j;
    }

    inline void SetF_turb(su2double val_F_turb_i, su2double val_F_turb_j) override {
      F_turb_i = val_F_turb_i;
      F_turb_j = val_F_turb_j;
    }

    inline void SetF_wake(su2double val_F_wake_i, su2double val_F_wake_j) override {
      F_wake_i = val_F_wake_i;
      F_wake_j = val_F_wake_j;
    }

    inline void Setrethetat_eq(su2double val_rethetat_eq_i, su2double val_rethetat_eq_j) override {
      rethetat_eq_i = val_rethetat_eq_i;
      rethetat_eq_j = val_rethetat_eq_j;
    }

    inline void SetT_param(su2double val_T_param_i, su2double val_T_param_j) override {
      T_param_i = val_T_param_i;
      T_param_j = val_T_param_j;
    }

    inline void Setdelta_param(su2double val_delta_param_i, su2double val_delta_param_j) override {
      delta_param_i = val_delta_param_i;
      delta_param_j = val_delta_param_j;
    }

    inline void SetReThetat_SCF(su2double val_ReThetat_SCF_i, su2double val_ReThetat_SCF_j) override {
      ReThetat_SCF_i = val_ReThetat_SCF_i;
      ReThetat_SCF_j = val_ReThetat_SCF_j;
    }

    inline void SetF_thetat_2(su2double val_F_thetat_2_i, su2double val_F_thetat_2_j) override {
      F_thetat_2_i = val_F_thetat_2_i;
      F_thetat_2_j = val_F_thetat_2_j;
    }


};
