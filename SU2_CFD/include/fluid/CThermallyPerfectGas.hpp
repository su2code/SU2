/*!
 * \file CThermallyPerfectGas.hpp
 * \brief Defines the ideal gas model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
 * \version 7.0.8 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "CFluidModel.hpp"

#if defined(HAVE_MPP) && !defined(CODI_REVERSE_TYPE) && !defined(CODI_FORWARD_TYPE)
#include "mutation++.h"

/*!
 * \class CThermallyPerfectGas
 * \brief Child class for defining the ideal gas model.
 * \author: S.Vitale, M.Pini.
 */
class CThermallyPerfectGas : public CFluidModel {
 protected:
  std::unique_ptr<Mutation::Mixture> mix; /*!< \brief Pointer to object Mixture from Mutation++ library. */
  su2double Gas_Constant{0.0};            /*!< \brief Gas Constant. */
  bool ComputeEntropy{true};              /*!< \brief Whether or not to compute entropy. */
  int ns, i;
  vector<su2double> rhos, energy, PT, temp, MolarMass, Cv_ks, Cvtrs, Cvves, rhoenergy;
  const su2double *cs;
  su2double RuSI{UNIVERSAL_GAS_CONSTANT};/*!< \brief Universal gas constant [J/(mol*K)] */
  su2double Ru{1000.0*RuSI};             /*!< \brief Universal gas constant [J/(kmol*K)] */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CThermallyPerfectGas(const CConfig* config, su2double R, bool CompEntropy = true);

  /*!
   * \brief Set the Dimensionless State using Density and Internal Energy
   * \param[in] rho - first thermodynamic variable.
   * \param[in] e - second thermodynamic variable.
   */
  void SetTDState_rhoe(su2double rho, su2double e) override;

  /*!
   * \brief Set the Dimensionless State using Pressure  and Temperature
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
  void SetEnergy_Prho(su2double P, su2double rho, su2double T) override;
  void SetEnergy_Prho(su2double P, su2double rho) override;

  /*!
   * \brief Set the Dimensionless State using Enthalpy and Entropy
   * \param[in] th1 - first thermodynamic variable (h).
   * \param[in] th2 - second thermodynamic variable (s).
   *
   */
  void SetTDState_hs(su2double h, su2double s) override;

  /*!
   * \brief Set the Dimensionless State using Density and Temperature
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

  //su2double ComputeGamma(vector<su2double> rhos);

};
#endif
