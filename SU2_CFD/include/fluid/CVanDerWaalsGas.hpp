/*!
 * \file CVanDerWaalsGas.hpp
 * \brief Declaration of the Polytropic Van der Waals model.
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
 * \class CVanDerWaalsGas
 * \brief Child class for defining the Van der Waals model.
 * \author: S.Vitale, M.Pini
 */
class CVanDerWaalsGas final : public CIdealGas {
 protected:
  su2double a{0.0};   /*!< \brief model parameter. */
  su2double b{0.0};   /*!< \brief model parameter. */
  su2double Zed{0.0}; /*!< \brief compressibility factor. */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CVanDerWaalsGas(su2double gamma, su2double R, su2double Pstar, su2double Tstar);

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
   * \brief Set the Dimensionless Internal Energy using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetEnergy_Prho(su2double P, su2double rho) override;

  /*!
   * \brief Set the Dimensionless state using Enthalpy and Entropy
   * \param[in] h - first thermodynamic variable (h).
   * \param[in] s - second thermodynamic variable (s).
   *
   */
  void SetTDState_hs(su2double h, su2double s) override;

  /*!
   * \brief Set the Dimensionless state using Density and Temperature
   * \param[in] rho - first thermodynamic variable (rho).
   * \param[in] T - second thermodynamic variable (T).
   *
   */
  void SetTDState_rhoT(su2double rho, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Entropy
   * \param[in] P - first thermodynamic variable (P).
   * \param[in] s - second thermodynamic variable (s).
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
};
