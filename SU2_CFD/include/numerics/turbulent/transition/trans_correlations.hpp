/*!
 * \file trans_correlations.hpp
 * \brief Numerics class for the LM model's correlation functions.
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

/*!
 * \class TransLMCorrelations
 * \brief Class for LM model's correlation functions.
 * \ingroup SourceDiscr
 * \author A. Rausa.
 */
class TransLMCorrelations {
 private:

  LM_ParsedOptions options;

 public:

  /*!
   * \brief Set LM options.
   * \param[in] val_options - LM options structure.
   */
  void SetOptions(const LM_ParsedOptions val_options){
    options = val_options;
  }

  /*!
   * \brief Compute Re_theta_c from correlations.
   * \param[in] Tu - Turbulence intensity.
   * \param[in] Re_theta_t - Re_theta_t (TransVar[1]).
   * \param[out] rethetac - Corrected value for Re_theta.
   */
  su2double ReThetaC_Correlations(const su2double Tu, const su2double Re_theta_t) const {

    su2double rethetac = 0.0;

    switch (options.Correlation) {
      case TURB_TRANS_CORRELATION::MALAN: {
        rethetac = min(0.615 * Re_theta_t + 61.5, Re_theta_t);
        break;
      }

      case TURB_TRANS_CORRELATION::SULUKSNA: {
        rethetac = min(0.1 * exp(-0.0022 * Re_theta_t + 12), 300.0);
        break;
      }

      case TURB_TRANS_CORRELATION::KRAUSE: {
        rethetac = 0.91 * Re_theta_t + 5.32;
        break;
      }

      case TURB_TRANS_CORRELATION::KRAUSE_HYPER: {
        const su2double FirstTerm = -0.042 * pow(Tu, 3);
        const su2double SecondTerm = 0.4233 * pow(Tu, 2);
        const su2double ThirdTerm = 0.0118 * pow(Tu, 1);
        rethetac = Re_theta_t / (FirstTerm + SecondTerm + ThirdTerm + 1.0744);
        break;
      }

      case TURB_TRANS_CORRELATION::MEDIDA_BAEDER: {
        const su2double FirstTerm = 4.45 * pow(Tu, 3);
        const su2double SecondTerm = 5.7 * pow(Tu, 2);
        const su2double ThirdTerm = 1.37 * pow(Tu, 1);
        rethetac = (FirstTerm - SecondTerm + ThirdTerm + 0.585) * Re_theta_t;
        break;
      }

      case TURB_TRANS_CORRELATION::MEDIDA: {
        rethetac = 0.62 * Re_theta_t;
        break;
      }

      case TURB_TRANS_CORRELATION::MENTER_LANGTRY: {
        if (Re_theta_t <= 1870) {
          const su2double FirstTerm = (-396.035 * pow(10, -2));
          const su2double SecondTerm = (10120.656 * pow(10, -4)) * Re_theta_t;
          const su2double ThirdTerm = (-868.230 * pow(10, -6)) * pow(Re_theta_t, 2);
          const su2double ForthTerm = (696.506 * pow(10, -9)) * pow(Re_theta_t, 3);
          const su2double FifthTerm = (-174.105 * pow(10, -12)) * pow(Re_theta_t, 4);
          rethetac = FirstTerm + SecondTerm + ThirdTerm + ForthTerm + FifthTerm;
        } else {
          rethetac = Re_theta_t - (593.11 + 0.482 * (Re_theta_t - 1870.0));
        }

        break;
      }
      case TURB_TRANS_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }

    return rethetac;
  }

  /*!
   * \brief Compute FLength from correlations.
   * \param[in] Tu - Turbulence intensity.
   * \param[in] Re_theta_t - Re_theta_t (TransVar[1]).
   * \param[out] F_length1 - Value for the F_length1 variable.
   */
  su2double FLength_Correlations(const su2double Tu, const su2double Re_theta_t) const {
    su2double F_length1 = 0.0;

    switch (options.Correlation) {
      case TURB_TRANS_CORRELATION::MALAN: {
        F_length1 = min(exp(7.168 - 0.01173 * Re_theta_t) + 0.5, 300.0);
        break;
      }

      case TURB_TRANS_CORRELATION::SULUKSNA: {
        const su2double FirstTerm = -pow(0.025 * Re_theta_t, 2) + 1.47 * Re_theta_t - 120.0;
        F_length1 = min(max(FirstTerm, 125.0), Re_theta_t);
        break;
      }

      case TURB_TRANS_CORRELATION::KRAUSE: {
        F_length1 = 3.39 * Re_theta_t + 55.03;
        break;
      }

      case TURB_TRANS_CORRELATION::KRAUSE_HYPER: {
        if (Tu <= 1.) {
          F_length1 = log(Re_theta_t + 1) / Tu;
        } else {
          const su2double FirstTerm = 0.2337 * pow(Tu, 2);
          const su2double SecondTerm = -1.3493 * pow(Tu, 1);
          F_length1 = log(Re_theta_t + 1) * (FirstTerm + SecondTerm + 2.1449);
        }
        break;
      }

      case TURB_TRANS_CORRELATION::MEDIDA_BAEDER: {
        const su2double FirstTerm = 0.171 * pow(Tu, 2);
        const su2double SecondTerm = 0.0083 * pow(Tu, 1);
        F_length1 = (FirstTerm - SecondTerm + 0.0306);
        break;
      }

      case TURB_TRANS_CORRELATION::MEDIDA: {
        F_length1 = 40;
        break;
      }

      case TURB_TRANS_CORRELATION::MENTER_LANGTRY: {
        if (Re_theta_t < 400) {
          F_length1 = 39.8189 + (-119.270 * pow(10, -4)) * Re_theta_t +
                      (-132.567 * pow(10, -6)) * Re_theta_t * Re_theta_t;
        } else if (Re_theta_t < 596) {
          F_length1 = 263.404 + (-123.939 * pow(10, -2)) * Re_theta_t +
                      (194.548 * pow(10, -5)) * pow(Re_theta_t, 2) +
                      (-101.695 * pow(10, -8)) * pow(Re_theta_t, 3);
        } else if (Re_theta_t < 1200) {
          F_length1 = 0.5 - (3.0 * pow(10, -4)) * (Re_theta_t - 596.0);
        } else {
          F_length1 = 0.3188;
        }
        break;
      }
      case TURB_TRANS_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }

    return F_length1;
  }
};
