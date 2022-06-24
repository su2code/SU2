/*!
 * \file CNonPolytropicPengRobinson.hpp
 * \brief Defines the non-polytropic Peng-Robinson model.
 * \author B. Fuentes Monjas
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

#include "CPengRobinson.hpp"

/*!
 * \class CNonPolytropicPengRobinson
 * \brief Child class for defining the non-polytropic Peng-Robinson model.
 * \author: B. Fuentes Monjas
 */
class CNonPolytropicPengRobinson : public CPengRobinson {

 public:
  /*!
   * \brief Constructor of the class.
   */
  CNonPolytropicPengRobinson(su2double gamma, su2double R, su2double Pstar, su2double Tstar, su2double w, const CConfig *config);

  /*!
   * \brief Set the temperature polynomial coefficients for variable Cp.
   * \param[in] config - configuration container for the problem.
   */
  virtual void SetCpModel(const CConfig* config) override;

  /*!
   * \brief Set the Dimensionless State using Density and Internal Energy
   * \param[in] rho - first thermodynamic variable.
   * \param[in] e - second thermodynamic variable.
   */
  virtual void SetTDState_rhoe(su2double rho, su2double e) override;

  /*!
   * \brief Set the Dimensionless State using Density and Temperature
   * \param[in] rho - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  virtual void SetTDState_rhoT(su2double rho, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Temperature
   * \param[in] P - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  virtual void SetTDState_PT(su2double P, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Enthalpy and Entropy
   * \param[in] h - first thermodynamic variable.
   * \param[in] s - second thermodynamic variable.
   */
  virtual void SetTDState_hs(su2double h, su2double s) override;

  /*!
   * \brief Set the Dimensionless Energy using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  virtual void SetEnergy_Prho(su2double P, su2double rho) override;

   /*!
   * \brief Set the Dimensionless State using Pressure and Entropy
   * \param[in] P - first thermodynamic variable.
   * \param[in] s - second thermodynamic variable.
   */
  virtual void SetTDState_Ps(su2double P, su2double s) override;

 private:

  /*!
   * \brief Internal function for the calculation of the integral of CP(IG)/T over T
   */
  su2double ComputeIntegralCp0_T_DT(su2double T);

  /*!
   * \brief Internal function for the calculation of the integral of CP(IG) over T
   */
  su2double ComputeIntegralCp0_DT(su2double T);

  /*!
   * \brief Internal function for the calculation of CV(IG) at the TD state Temperature
   */
  su2double GetCv0();

  /*!
   * \brief Internal function for the implicit call hs.
   */
  virtual su2double T_v_h(su2double v, su2double h) override;

  array<su2double, N_POLY_COEFFS> coeffs_; /*!< \brief Polynomial coefficients for heat capacity as a function of temperature. */
};