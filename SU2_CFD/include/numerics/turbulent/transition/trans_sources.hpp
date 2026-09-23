/*!
 * \file trans_sources.hpp
 * \brief Numerics classes for integration of source terms in transition problems.
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../scalar/scalar_sources.hpp"
#include "./trans_correlations.hpp"

/*--- Lower limit of the RMS roughness height in the log(h / theta_t) of the Langtry et al. stationary cross-flow
 * correlation, which diverges for h = 0: 0.25 micrometers, the smallest roughness for which the correlation was
 * validated by Lee and Baeder (AIAA 2021-1532, Radeztsky et al. data) and the reference height h0 of Vallinayagam
 * Pillai and Lardeau (AIAA 2017-3159, Eq. 9). Smoother surfaces are treated as this one. HROUGHNESS is in meters. ---*/
constexpr passivedouble LANGTRY_CF_MIN_ROUGHNESS = 0.25e-6;

/*!
 * \class CSourcePieceWise_TranLM
 * \brief Class for integrating the source terms of the LM transition model equations.
 * \ingroup SourceDiscr
 * \author S. Kang.
 */
template <class FlowIndices>
class CSourcePieceWise_TransLM final : public CNumerics {
 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */

  const LM_ParsedOptions options;

  /*--- LM Closure constants ---*/
  const su2double c_e1 = 1.0;
  const su2double c_a1 = 2.0;
  const su2double c_e2 = 50.0;
  const su2double c_a2 = 0.06;
  const su2double sigmaf = 1.0;
  const su2double s1 = 2.0;
  const su2double c_theta = 0.03;
  const su2double c_CF = 0.6;
  const su2double sigmat = 2.0;

  TURB_FAMILY TurbFamily;
  su2double hRoughness;

  su2double Re_v_Here;
  su2double Corr_Rec_Here;
  su2double Prod_Here = 0.0;
  su2double Destr_Here = 0.0;
  su2double F_onset1_Here = 0.0;
  su2double F_onset2_Here = 0.0;
  su2double F_onset3_Here = 0.0;
  su2double F_onset_Here = 0.0;
  su2double lambda_theta_Here = 0.0;
  su2double duds_Here = 0.0;

  su2double IntermittencySep = 1.0;
  su2double IntermittencyEff = 1.0;

  su2double Residual[2];
  su2double* Jacobian_i[2];
  su2double Jacobian_Buffer[4];  // Static storage for the Jacobian (which needs to be pointer for return type).

  TransLMCorrelations TransCorrelations;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TransLM(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
      : CNumerics(val_nDim, 2, config), idx(val_nDim, config->GetnSpecies()), options(config->GetLMParsedOptions()){
    /*--- "Allocate" the Jacobian using the static buffer. ---*/
    Jacobian_i[0] = Jacobian_Buffer;
    Jacobian_i[1] = Jacobian_Buffer + 2;

    TurbFamily = TurbModelFamily(config->GetKind_Turb_Model());

    hRoughness = config->GethRoughness();

    TransCorrelations.SetOptions(options);

  }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override {
    /*--- ScalarVar[0] = k, ScalarVar[0] = w, TransVar[0] = gamma, and TransVar[1] = ReThetaT ---*/
    /*--- dU/dx = PrimVar_Grad[1][0] ---*/
    AD::StartPreacc();
    AD::SetPreaccIn(StrainMag_i);
    AD::SetPreaccIn(ScalarVar_i, nVar);
    AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
    AD::SetPreaccIn(TransVar_i, nVar);
    AD::SetPreaccIn(TransVar_Grad_i, nVar, nDim);
    AD::SetPreaccIn(Volume);
    AD::SetPreaccIn(dist_i);
    AD::SetPreaccIn(&V_i[idx.Velocity()], nDim);
    AD::SetPreaccIn(PrimVar_Grad_i, nDim + idx.Velocity(), nDim);
    AD::SetPreaccIn(Vorticity_i, 3);

    const su2double VorticityMag = max(GeometryToolbox::Norm(3, Vorticity_i), 1e-20);

    const su2double vel_u = V_i[idx.Velocity()];
    const su2double vel_v = V_i[1 + idx.Velocity()];
    const su2double vel_w = (nDim == 3) ? V_i[2 + idx.Velocity()] : 0.0;

    const su2double Velocity_Mag = max(sqrt(vel_u * vel_u + vel_v * vel_v + vel_w * vel_w), 1e-20);

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

    if (dist_i > 1e-10) {
      su2double Tu = 1.0;
      if (TurbFamily == TURB_FAMILY::KW) Tu = max(100.0 * sqrt(2.0 * ScalarVar_i[0] / 3.0) / Velocity_Mag, 0.027);
      if (TurbFamily == TURB_FAMILY::SA) Tu = config->GetTurbulenceIntensity_FreeStream() * 100;

      /*--- Corr_RetC correlation*/
      const su2double Corr_Rec = TransCorrelations.ReThetaC_Correlations(Tu, TransVar_i[1]);
       // AGGIUNTO PER DEBUG
      Corr_Rec_Here = Corr_Rec;

      /*--- F_length correlation*/
      const su2double Corr_F_length = TransCorrelations.FLength_Correlations(Tu, TransVar_i[1]);

      /*--- F_length ---*/
      su2double F_length = 0.0;
      if (TurbFamily == TURB_FAMILY::KW) {
        const su2double r_omega = Density_i * dist_i * dist_i * ScalarVar_i[1] / Laminar_Viscosity_i;
        const su2double f_sub = exp(-pow(r_omega / 200.0, 2));
        F_length = Corr_F_length * (1. - f_sub) + 40.0 * f_sub;
      }
      if (TurbFamily == TURB_FAMILY::SA) F_length = Corr_F_length;

      /*--- F_onset ---*/
      su2double R_t = 1.0;
      if (TurbFamily == TURB_FAMILY::KW) R_t = Density_i * ScalarVar_i[0] / Laminar_Viscosity_i / ScalarVar_i[1];
      if (TurbFamily == TURB_FAMILY::SA) R_t = Eddy_Viscosity_i / Laminar_Viscosity_i;

      const su2double Re_v = Density_i * dist_i * dist_i * StrainMag_i / Laminar_Viscosity_i;
      // AGGIUNTO PER DEBUG
      Re_v_Here = Re_v;

      const su2double F_onset1 = Re_v / (2.193 * Corr_Rec);
      su2double F_onset2 = 1.0;
      su2double F_onset3 = 1.0;
      if (TurbFamily == TURB_FAMILY::KW) {
        F_onset2 = min(max(F_onset1, pow(F_onset1, 4.0)), 2.0);
        F_onset3 = max(1.0 - pow(R_t / 2.5, 3.0), 0.0);
      }
      if (TurbFamily == TURB_FAMILY::SA) {
        F_onset2 = min(max(F_onset1, pow(F_onset1, 4.0)), 4.0);
        F_onset3 = max(2.0 - pow(R_t / 2.5, 3.0), 0.0);
      }
      const su2double F_onset = max(F_onset2 - F_onset3, 0.0);
      // AGGIUNTO PER DEBUG
      F_onset1_Here = F_onset1;
      F_onset2_Here = F_onset2;
      F_onset3_Here = F_onset3;
      F_onset_Here = F_onset;

      /*-- Gradient of velocity magnitude ---*/

      su2double dU_dx = 0.5 / Velocity_Mag * (2. * vel_u * PrimVar_Grad_i[1][0] + 2. * vel_v * PrimVar_Grad_i[2][0]);
      if (nDim == 3) dU_dx += 0.5 / Velocity_Mag * (2. * vel_w * PrimVar_Grad_i[3][0]);

      su2double dU_dy = 0.5 / Velocity_Mag * (2. * vel_u * PrimVar_Grad_i[1][1] + 2. * vel_v * PrimVar_Grad_i[2][1]);
      if (nDim == 3) dU_dy += 0.5 / Velocity_Mag * (2. * vel_w * PrimVar_Grad_i[3][1]);

      su2double dU_dz = 0.0;
      if (nDim == 3)
        dU_dz =
            0.5 / Velocity_Mag *
            (2. * vel_u * PrimVar_Grad_i[1][2] + 2. * vel_v * PrimVar_Grad_i[2][2] + 2. * vel_w * PrimVar_Grad_i[3][2]);

      su2double du_ds = vel_u / Velocity_Mag * dU_dx + vel_v / Velocity_Mag * dU_dy;
      if (nDim == 3) du_ds += vel_w / Velocity_Mag * dU_dz;

      /*-- Calculate blending function f_theta --*/
      su2double time_scale = 500.0 * Laminar_Viscosity_i / Density_i / Velocity_Mag / Velocity_Mag;
      if (options.CrossFlow)
        time_scale = min(time_scale,
                         Density_i * LocalGridLength_i * LocalGridLength_i / (Laminar_Viscosity_i + Eddy_Viscosity_i));
      const su2double theta_bl = TransVar_i[1] * Laminar_Viscosity_i / Density_i / Velocity_Mag;
      const su2double delta_bl = 7.5 * theta_bl;
      const su2double delta = 50.0 * VorticityMag * dist_i / Velocity_Mag * delta_bl + 1e-20;

      su2double f_wake = 0.0;
      if (TurbFamily == TURB_FAMILY::KW) {
        const su2double re_omega = Density_i * ScalarVar_i[1] * dist_i * dist_i / Laminar_Viscosity_i;
        f_wake = exp(-pow(re_omega / (1.0e+05), 2));
      }
      if (TurbFamily == TURB_FAMILY::SA) f_wake = 1.0;

      const su2double var1 = (TransVar_i[0] - 1.0 / c_e2) / (1.0 - 1.0 / c_e2);
      const su2double var2 = 1.0 - pow(var1, 2.0);
      const su2double f_theta = min(max(f_wake * exp(-pow(dist_i / delta, 4)), var2), 1.0);
      const su2double f_turb = exp(-pow(R_t / 4, 4));

      su2double f_theta_2 = 0.0;
      if (options.CrossFlow)
        f_theta_2 = min(f_wake * exp(-pow(dist_i / delta, 4.0)), 1.0);

      /*--- Corr_Ret correlation*/
      const su2double Corr_Ret_lim = 20.0;
      su2double f_lambda = 1.0;

      su2double Retheta_Error = 200.0, Retheta_old = 1.0;
      su2double lambda = 0.0;
      su2double Corr_Ret = 20.0;

      for (int iter = 0; iter < 100; iter++) {
        su2double theta = Corr_Ret * Laminar_Viscosity_i / Density_i / Velocity_Mag;
        lambda = Density_i * theta * theta / Laminar_Viscosity_i * du_ds;
        lambda = min(max(-0.1, lambda), 0.1);

        if (lambda <= 0.0) {
          f_lambda = 1. - (-12.986 * lambda - 123.66 * lambda * lambda - 405.689 * lambda * lambda * lambda) *
                              exp(-pow(Tu / 1.5, 1.5));
        } else {
          f_lambda = 1. + 0.275 * (1. - exp(-35. * lambda)) * exp(-Tu / 0.5);
        }

        if (Tu <= 1.3) {
          Corr_Ret = f_lambda * (1173.51 - 589.428 * Tu + 0.2196 / Tu / Tu);
        } else {
          Corr_Ret = 331.5 * f_lambda * pow(Tu - 0.5658, -0.671);
        }
        Corr_Ret = max(Corr_Ret, Corr_Ret_lim);

        Retheta_Error = fabs(Retheta_old - Corr_Ret) / Retheta_old;

        if (Retheta_Error < 0.0000001) {
          break;
        }

        Retheta_old = Corr_Ret;
      }

      // DEBUG
      lambda_theta_Here = lambda;
      duds_Here = du_ds;

      /*-- Corr_RetT_SCF Correlations--*/
      su2double ReThetat_SCF = 0.0;
      if (options.CrossFlow) {
        su2double VelocityNormalized[3];
        VelocityNormalized[0] = vel_u / Velocity_Mag;
        VelocityNormalized[1] = vel_v / Velocity_Mag;
        if (nDim == 3) VelocityNormalized[2] = vel_w / Velocity_Mag;

        su2double StreamwiseVort = 0.0;
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          StreamwiseVort += VelocityNormalized[iDim] * Vorticity_i[iDim];
        }
        StreamwiseVort = abs(StreamwiseVort);

        const su2double H_CF = StreamwiseVort * dist_i / Velocity_Mag;
        const su2double DeltaH_CF = H_CF * (1.0 + min(Eddy_Viscosity_i / Laminar_Viscosity_i, 0.4));
        const su2double DeltaH_CF_Minus = max(-1.0 * (0.1066 - DeltaH_CF), 0.0);
        const su2double DeltaH_CF_Plus = max(0.1066 - DeltaH_CF, 0.0);
        const su2double fDeltaH_CF_Minus = 75.0 * tanh(DeltaH_CF_Minus / 0.0125);
        const su2double fDeltaH_CF_Plus = 6200 * DeltaH_CF_Plus + 50000 * DeltaH_CF_Plus * DeltaH_CF_Plus;

        const su2double toll = 1e-5;
        su2double error = toll + 1.0;
        su2double thetat_SCF = 0.0;
        su2double rethetat_SCF_old = 20.0;
        const int nMax = 100;

        int iter;
        for (iter = 0; iter < nMax && error > toll; iter++) {
          thetat_SCF = rethetat_SCF_old * Laminar_Viscosity_i / (Density_i * (Velocity_Mag / 0.82));
          thetat_SCF = max(1e-20, thetat_SCF);

          ReThetat_SCF = -35.088 * log(max(hRoughness, su2double(LANGTRY_CF_MIN_ROUGHNESS)) / thetat_SCF) + 319.51 +
                         fDeltaH_CF_Plus - fDeltaH_CF_Minus;

          error = abs(ReThetat_SCF - rethetat_SCF_old) / rethetat_SCF_old;

          rethetat_SCF_old = ReThetat_SCF;
        }
      }

      /*-- production term of Intermeittency(Gamma) --*/
      const su2double Pg =
          F_length * c_a1 * Density_i * StrainMag_i * sqrt(F_onset * TransVar_i[0]) * (1.0 - c_e1 * TransVar_i[0]);

      /*-- destruction term of Intermeittency(Gamma) --*/
      const su2double Dg = c_a2 * Density_i * VorticityMag * TransVar_i[0] * f_turb * (c_e2 * TransVar_i[0] - 1.0);

      // DEBUG
      Prod_Here = Pg;
      Destr_Here = Dg;

      /*-- production term of ReThetaT --*/
      const su2double PRethetat = c_theta * Density_i / time_scale * (Corr_Ret - TransVar_i[1]) * (1.0 - f_theta);

      /*-- destruction term of ReThetaT --*/
      // It should not be with the minus sign but I put for consistency
      su2double DRethetat = 0.0;
      if (options.CrossFlow)
        DRethetat = -c_theta * (Density_i / time_scale) * c_CF * min(ReThetat_SCF - TransVar_i[1], 0.0) * f_theta_2;

      /*--- Source ---*/
      Residual[0] += (Pg - Dg) * Volume;
      Residual[1] += (PRethetat - DRethetat) * Volume;

      /*--- Implicit part ---*/
      Jacobian_i[0][0] = (F_length * c_a1 * StrainMag_i * sqrt(F_onset) *
                              (0.5 * pow(TransVar_i[0], -0.5) - 1.5 * c_e1 * pow(TransVar_i[0], 0.5)) -
                          c_a2 * VorticityMag * f_turb * (2.0 * c_e2 * TransVar_i[0] - 1.0)) *
                         Volume;
      Jacobian_i[0][1] = 0.0;
      Jacobian_i[1][0] = 0.0;
      Jacobian_i[1][1] = -c_theta / time_scale * (1.0 - f_theta) * Volume;
      if (options.CrossFlow && ReThetat_SCF - TransVar_i[1] < 0)
        Jacobian_i[1][1] += (c_theta / time_scale) * c_CF * f_theta_2 * Volume;
    }

    AD::SetPreaccOut(Residual, nVar);
    AD::EndPreacc();

    return ResidualType<>(Residual, Jacobian_i, nullptr);
  }

  inline su2double GetRe_v() override {return Re_v_Here;} 
  inline su2double GetCorr_Rec() override {return Corr_Rec_Here;} 
  inline su2double GetProd() override {return Prod_Here;} 
  inline su2double GetDestr() override {return Destr_Here;} 
  inline su2double GetF_onset1() override {return F_onset1_Here;} 
  inline su2double GetF_onset2() override {return F_onset2_Here;} 
  inline su2double GetF_onset3() override {return F_onset3_Here;} 
  inline su2double GetF_onset() override {return F_onset_Here;} 
  inline su2double GetLambda_theta() override {return lambda_theta_Here;} 
  inline su2double Getduds() override {return duds_Here;} 

};

/*!
 * \class CSourcePieceWise_TranSLM
 * \brief Class for integrating the source terms of the Simplified LM transition model equations.
 * \ingroup SourceDiscr
 * \author S. Kang.
 */
template <class FlowIndices>
class CSourcePieceWise_TransSLM final : public CNumerics {
 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */

  const LM_ParsedOptions options;

  /*--- LM Closure constants ---*/
  const su2double c_e1 = 1.0;
  const su2double c_a1 = 2.0;
  const su2double c_e2 = 50.0;
  const su2double c_a2 = 0.06;
  const su2double sigmaf = 1.0;
  const su2double s1 = 2.0;
  const su2double c_theta = 0.03;
  const su2double c_CF = 0.6;
  const su2double sigmat = 2.0;

  TURB_FAMILY TurbFamily;
  su2double hRoughness;

  su2double IntermittencySep = 1.0;
  su2double IntermittencyEff = 1.0;

  su2double Re_t;
  su2double Corr_Rec = 1.0;
  su2double AuxVar = 0.0;  /*!< \brief Wall-normal derivative of the wall-normal velocity (Menter correlation). */
  su2double CrossFlowPsi = 0.0;  /*!< \brief Wall-normal change of the vorticity direction times the wall distance. */
  su2double F2;
  su2double Tu_Here = 0.0;
  su2double duds_Here = 0.0;
  su2double lambda_theta_Here = 0.0;
  su2double Re_v_Here = 0.0;
  su2double Prod_Here = 0.0;
  su2double Destr_Here = 0.0;
  su2double F_onset1_Here = 0.0;
  su2double F_onset2_Here = 0.0;
  su2double F_onset3_Here = 0.0;
  su2double F_onset_Here = 0.0;

  su2double Residual;
  su2double* Jacobian_i;
  su2double Jacobian_Buffer;  // Static storage for the Jacobian (which needs to be pointer for return type).

  TransLMCorrelations TransCorrelations;

  /*!
   * \brief Magnitude of the streamwise vorticity, |U/|U| . Omega|.
   */
  su2double StreamwiseVorticity(su2double vel_u, su2double vel_v, su2double vel_w, su2double Velocity_Mag) const {
    const su2double velocity[3] = {vel_u, vel_v, vel_w};
    su2double streamwiseVort = 0.0;
    for (auto iDim = 0u; iDim < nDim; iDim++) streamwiseVort += velocity[iDim] / Velocity_Mag * Vorticity_i[iDim];
    return fabs(streamwiseVort);
  }

  /*!
   * \brief Stationary cross-flow Reynolds number of Langtry et al., Lee and Baeder (AIAA 2021-1532), Eqs. 36-42.
   *        Eq. 36 is implicit in theta_t and is solved by fixed-point iteration, as in the LM2015 option of
   *        CSourcePieceWise_TransLM.
   */
  su2double StationaryCrossFlowReynolds(su2double vel_u, su2double vel_v, su2double vel_w, su2double Velocity_Mag,
                                        su2double R_t, const CConfig* config) const {
    const su2double H_CF = StreamwiseVorticity(vel_u, vel_v, vel_w, Velocity_Mag) * dist_i / Velocity_Mag;  // Eq. 37
    const su2double DeltaH_CF = H_CF * (1.0 + min(R_t, 0.4));                                            // Eq. 38
    const su2double DeltaH_CF_Plus = max(0.1066 - DeltaH_CF, 0.0);                                       // Eq. 39
    const su2double fDeltaH_CF_Plus = 6200 * DeltaH_CF_Plus + 50000 * DeltaH_CF_Plus * DeltaH_CF_Plus;   // Eq. 40
    const su2double DeltaH_CF_Minus = max(-(0.1066 - DeltaH_CF), 0.0);                                   // Eq. 41
    const su2double fDeltaH_CF_Minus = 75.0 * tanh(DeltaH_CF_Minus / 0.0125);                            // Eq. 42

    const su2double h = max(config->GethRoughness(), su2double(LANGTRY_CF_MIN_ROUGHNESS));
    su2double reScf = 20.0, error = 1.0;
    for (int iter = 0; iter < 100 && error > 1e-5; iter++) {
      const su2double theta_t = max(1e-20, reScf * Laminar_Viscosity_i / (Density_i * Velocity_Mag / 0.82));
      const su2double reScfNew = -35.088 * log(h / theta_t) + 319.51 + fDeltaH_CF_Plus - fDeltaH_CF_Minus;
      error = fabs(reScfNew - reScf) / fabs(reScf);
      reScf = reScfNew;
    }
    return max(reScf, 1e-20);
  }

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TransSLM(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
      : CNumerics(val_nDim, 1, config), idx(val_nDim, config->GetnSpecies()), options(config->GetLMParsedOptions()){
    /*--- "Allocate" the Jacobian using the static buffer. ---*/
    Jacobian_i = &Jacobian_Buffer;

    TurbFamily = TurbModelFamily(config->GetKind_Turb_Model());

    hRoughness = config->GethRoughness();

    TransCorrelations.SetOptions(options);
    
  }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override {
    /*--- ScalarVar[0] = k, ScalarVar[0] = w, TransVar[0] = gamma ---*/
    /*--- dU/dx = PrimVar_Grad[1][0] ---*/
    AD::StartPreacc();
    AD::SetPreaccIn(StrainMag_i);
    /*--- ScalarVar_i holds the turbulence variables (k and omega for SST), nVar is the single transition variable. ---*/
    const unsigned short nVarTurb = (TurbFamily == TURB_FAMILY::KW) ? 2 : 1;
    AD::SetPreaccIn(ScalarVar_i, nVarTurb);
    AD::SetPreaccIn(ScalarVar_Grad_i, nVarTurb, nDim);
    AD::SetPreaccIn(AuxVar);
    AD::SetPreaccIn(CrossFlowPsi);
    AD::SetPreaccIn(TransVar_i, nVar);
    AD::SetPreaccIn(TransVar_Grad_i, nVar, nDim);
    AD::SetPreaccIn(Volume);
    AD::SetPreaccIn(dist_i);
    AD::SetPreaccIn(&V_i[idx.Velocity()], nDim);
    AD::SetPreaccIn(PrimVar_Grad_i, nDim + idx.Velocity(), nDim);
    AD::SetPreaccIn(Vorticity_i, 3);

    const su2double VorticityMag = GeometryToolbox::Norm(3, Vorticity_i);

    const su2double vel_u = V_i[idx.Velocity()];
    const su2double vel_v = V_i[1 + idx.Velocity()];
    const su2double vel_w = (nDim == 3) ? V_i[2 + idx.Velocity()] : 0.0;

    const su2double Velocity_Mag = max(sqrt(vel_u * vel_u + vel_v * vel_v + vel_w * vel_w), 1e-20);

    AD::SetPreaccIn(V_i[idx.Density()], V_i[idx.LaminarViscosity()], V_i[idx.EddyViscosity()]);

    Density_i = V_i[idx.Density()];
    Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];
    Eddy_Viscosity_i = V_i[idx.EddyViscosity()];

    Residual = 0.0;
    Jacobian_i[0] = 0.0;

    if (dist_i > 1e-10) {
      su2double Tu_L = 1.0;
      // Local value of the Turbulence intensity that makes it galileian invariant. Look at Eq. 7 in https://doi.org/10.1007/s10494-015-9622-4
      if (TurbFamily == TURB_FAMILY::KW) Tu_L = min(100.0 * sqrt(2.0 * ScalarVar_i[0] / 3.0) / (ScalarVar_i[1]*dist_i), 100.0);
      if (TurbFamily == TURB_FAMILY::SA) Tu_L = config->GetTurbulenceIntensity_FreeStream() * 100;

      Tu_Here = Tu_L;

      /*--- F_length ---*/
      const su2double F_length = 100.0;  // Menter et al. (2015), Eq. 6, and Lee and Baeder (2021), Eq. 7.

      /*--- F_onset ---*/
      su2double R_t = 1.0;
      if (TurbFamily == TURB_FAMILY::KW) R_t = Density_i * ScalarVar_i[0] / (Laminar_Viscosity_i * ScalarVar_i[1]);
      /*--- Lee and Baeder (AIAA 2021-1532), Eq. 6: k and omega are replaced by the viscosity ratio for SA. ---*/
      if (TurbFamily == TURB_FAMILY::SA) R_t = Eddy_Viscosity_i / Laminar_Viscosity_i;

      /*--- Menter et al. (2015), Eqs. 12-13, and Lee and Baeder (2021), Eq. 9. ---*/
      const su2double lambda_theta = max(min(-7.57e-3 * AuxVar * dist_i * dist_i * Density_i / Laminar_Viscosity_i + 0.0128, 1.0), -1.0);

      duds_Here = AuxVar;
      lambda_theta_Here = lambda_theta;

      /*--- Critical Reynolds number, Menter et al. (2015), Eq. 14, and Lee and Baeder (2021), Eqs. 10-13. ---*/
      Re_t = TransCorrelations.ReThetaC_Correlations_SLM(Tu_L, lambda_theta, dist_i, VorticityMag, Velocity_Mag,
                                                          TurbFamily == TURB_FAMILY::SA);
      Corr_Rec = Re_t;

      const su2double Re_v = Density_i * dist_i * dist_i * StrainMag_i / Laminar_Viscosity_i;
      Re_v_Here = Re_v;

      /*--- Menter et al. (2015), Eqs. 4-5, and Lee and Baeder (2021), Eqs. 4-5 (the same for SST and SA). ---*/
      su2double F_onset1 = Re_v / (2.2 * Corr_Rec);

      if (options.CrossFlow && TurbFamily == TURB_FAMILY::SA) {
        /*--- Langtry et al. stationary cross-flow criterion, Lee and Baeder (2021), Eqs. 36-43. ---*/
        const su2double Re_scf = StationaryCrossFlowReynolds(vel_u, vel_v, vel_w, Velocity_Mag, R_t, config);
        const su2double c_scf = 0.94;
        F_onset1 = max(F_onset1, c_scf * Re_v / (2.2 * Re_scf));
      }

      const su2double F_onset2 = min(F_onset1, 2.0);
      const su2double F_onset3 = max(1.0 - pow(R_t / 3.5, 3.0), 0.0);
      su2double F_onset = max(F_onset2 - F_onset3, 0.0);

      F_onset1_Here = F_onset1;
      F_onset2_Here = F_onset2;
      F_onset3_Here = F_onset3;
      F_onset_Here = F_onset;

      if (options.CrossFlow && TurbFamily == TURB_FAMILY::SA) {
        /*--- Menter and Smirnov C1-based cross-flow criterion, Lee and Baeder (2021), Eqs. 25-35. ---*/
        const su2double lambda_CF = min(max(-7.57e-3 * AuxVar * dist_i * dist_i * Density_i / Laminar_Viscosity_i + 0.0174, 0.0), 0.0477);
        const su2double g_CF = min(max(27864.0 * pow(lambda_CF, 3) - 1962.0 * pow(lambda_CF, 2) + 54.3 * lambda_CF + 1.0, 1.0), 2.3);
        const su2double G_CF = 0.684 / g_CF;
        const su2double C_RSF = 1.35;  // Calibrated for SA by Lee and Baeder (1.0 in the original model).
        const su2double T_C1 = C_RSF / 150.0 * G_CF * CrossFlowPsi * Re_v;
        const su2double F_onset_CF = min(max(100.0 * (T_C1 - 1.0), 0.0), 1.0);
        F_onset = max(F_onset, F_onset_CF);
      }

      if (options.CrossFlow && TurbFamily == TURB_FAMILY::KW) {

        /*--- Vallinayagam Pillai and Lardeau, "Accounting crossflow effects in one-equation local correlation-based
         * transition model", AIAA 2017-3159 (equation numbers below). ---*/

        // Shape factor, Eq. 3 with k = 0.25 - lambda (lambda_theta_L of the gamma model), limited to 2.7 above which
        // the crossflow criterion does not apply (text after Eq. 3)
        const su2double k = 0.25 - lambda_theta;
        const su2double FirstTerm = 4.14 * k;
        const su2double SecondTerm = 83.5 * pow(k, 2.0);
        const su2double ThirdTerm = 854.0 * pow(k, 3.0);
        const su2double ForthTerm = 3337.0 * pow(k, 4.0);
        const su2double FifthTerm = 4576.0 * pow(k, 5.0);
        const su2double H = min(2.0 + FirstTerm - SecondTerm + ThirdTerm - ForthTerm + FifthTerm, 2.7);

        // Critical crossflow Reynolds number, Eq. 2
        su2double Re_Crit_CF = 0.0;
        if(H < 2.3) {
          Re_Crit_CF = 150.0;
        } else {
          // Eq. 2 is printed with a minus sign, which gives -150 at H = 2.3 and negative critical Reynolds numbers;
          // the positive sign makes it continuous with the value 150 below H = 2.3
          Re_Crit_CF = (300.0/PI_NUMBER) * atan(0.106/(pow(H-2.3, 2.05)));
        }

        const su2double H_CF = StreamwiseVorticity(vel_u, vel_v, vel_w, Velocity_Mag) * dist_i / Velocity_Mag;

        // Crossflow strength H_cf, Eqs. 4-6, and Delta_H_cf, Eq. 8
        const su2double Delta_H_CF = H_CF * (1.0 + min(Eddy_Viscosity_i / Laminar_Viscosity_i, 0.4));

        // Roughness, Eq. 9, h0 = 0.25 micrometers (HROUGHNESS in meters)
        const su2double h_0 = 0.25e-6;
        const su2double C_r = 2.0 - pow(0.5, config->GethRoughness()/h_0);

        // f_cf, Eqs. 7 and 10, with C_cf = 1
        const su2double C_CF = 1.0;
        const su2double f_CF = (C_CF * C_r * Delta_H_CF * Corr_Rec) / Re_Crit_CF;
        const su2double F_onset_CF = min(max(0.0, f_CF - 1.0), 1.0);

        // Eqs. 17-18
        F_onset = max(F_onset, F_onset_CF);

      }

      /*--- Menter et al. (2015), Eq. 5, and Lee and Baeder (2021), Eq. 6. ---*/
      const su2double f_turb = exp(-pow(R_t / 2, 4));

      /*--- Production and destruction of the intermittency, Menter et al. (2015), Eqs. 2-3, and Lee and Baeder
       * (2021), Eqs. 2-3, written as gamma*(P - D) with (Eqs. 22-23) ---*/
      const su2double P = F_length * StrainMag_i * (1.0 - TransVar_i[0]) * F_onset;
      const su2double D = c_a2 * VorticityMag * f_turb * (c_e2 * TransVar_i[0] - 1.0);
      const su2double Pg = Density_i * TransVar_i[0] * P;
      const su2double Dg = Density_i * TransVar_i[0] * D;

      Prod_Here = Pg;
      Destr_Here = Dg;

      /*--- Source ---*/
      Residual += (Pg - Dg) * Volume;

      /*--- Implicit part, per unit rho*gamma. ---*/
      if (TurbFamily == TURB_FAMILY::SA) {
        /*--- Positivity of the implicit operator, Lee and Baeder (2021), Eq. 24. ---*/
        const su2double dP = -F_length * StrainMag_i * F_onset;
        const su2double dD = c_a2 * VorticityMag * f_turb * c_e2;
        Jacobian_i[0] = -(max(D - P, 0.0) + max(dD - dP, 0.0) * TransVar_i[0]) * Volume;
      } else {
        Jacobian_i[0] = (F_length * StrainMag_i * F_onset * (1 - 2*TransVar_i[0]) -
                         c_a2 * VorticityMag * f_turb * (2.0 * c_e2 * TransVar_i[0] - 1.0)) * Volume;
      }

    }

    AD::SetPreaccOut(Residual);
    AD::EndPreacc();

    return ResidualType<>(&Residual, &Jacobian_i, nullptr);
    
  }

  inline su2double GetRe_t() override {return Re_t;}
  inline su2double GetCorr_Rec() override {return Corr_Rec;} 
  inline su2double GetTu() override {return Tu_Here;}
  inline su2double GetLambda_theta() override {return lambda_theta_Here;} 
  inline su2double Getduds() override {return duds_Here;} 
  inline su2double GetRe_v() override {return Re_v_Here;} 
  inline su2double GetProd() override {return Prod_Here;} 
  inline su2double GetDestr() override {return Destr_Here;} 
  inline su2double GetF_onset1() override {return F_onset1_Here;} 
  inline su2double GetF_onset2() override {return F_onset2_Here;} 
  inline su2double GetF_onset3() override {return F_onset3_Here;} 
  inline su2double GetF_onset() override {return F_onset_Here;} 
  inline void SetAuxVar(su2double val_AuxVar) override { AuxVar = val_AuxVar;}
  inline void SetCrossFlowStrength(su2double val_Psi) override { CrossFlowPsi = val_Psi; }
  // non serve più
  inline void SetF2(su2double val_F2) override { F2 = val_F2;}

};