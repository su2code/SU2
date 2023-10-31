/*!
 * \file CPengRobinson.hpp
 * \brief Defines the Peng-Robinson model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
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

#include "CIdealGas.hpp"

/*!
 * \class CPengRobinson
 * \brief Child class for defining the Peng-Robinson model.
 * \author: S.Vitale, G. Gori
 */
class CPengRobinson final : public CIdealGas {
 protected:
  su2double a{0.0};         /*!< \brief model parameter. */
  su2double b{0.0};         /*!< \brief model parameter. */
  su2double k{0.0};         /*!< \brief model parameter (computed with acentric factor). */
  su2double Zed{0.0};       /*!< \brief compressibility factor. */
  su2double TstarCrit{0.0}; /*!< \brief Critical temperature. */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CPengRobinson(su2double gamma, su2double R, su2double Pstar, su2double Tstar, su2double w);

  /*!
   * \brief Set the Dimensionless State using Density and Internal Energy
   * \param[in] rho - first thermodynamic variable.
   * \param[in] e - second thermodynamic variable.
   */
  void SetTDState_rhoe(su2double rho, su2double e) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Temperature
   * \param[in] P - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  void SetTDState_PT(su2double P, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetTDState_Prho(su2double P, su2double rho) override;

  /*!
   * \brief Set the Dimensionless Energy using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetEnergy_Prho(su2double P, su2double rho) override;

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("hs").
   * \param[in] th1 - first thermodynamic variable (h).
   * \param[in] th2 - second thermodynamic variable (s).
   *
   */
  void SetTDState_hs(su2double h, su2double s) override;

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
   * \param[in] th1 - first thermodynamic variable (rho).
   * \param[in] th2 - second thermodynamic variable (T).
   *
   */
  void SetTDState_rhoT(su2double rho, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Entropy
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (s).
   */
  void SetTDState_Ps(su2double P, su2double s) override;

  /*!
   * \brief compute some derivatives of enthalpy and entropy needed for subsonic inflow BC
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (v).
   *
   */
  void ComputeDerivativeNRBC_Prho(su2double P, su2double rho) override;

 private:
  /*!
   * \brief Internal model parameter.
   */
  su2double alpha2(su2double T) const;

  /*!
   * \brief Internal function for the implicit call hs.
   */
  su2double T_v_h(su2double v, su2double h);

  /*!
   * \brief Internal function for the implicit call Ps.
   */
  su2double T_P_rho(su2double P, su2double rho);
};
