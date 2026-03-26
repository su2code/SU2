/*!
 * \file CNasaGas.hpp
 * \brief Defines the NASA polynomial gas model.
 * \author Om Giri
 * \version 8.4.0 "Harrier"
 */

#pragma once

#include "CIdealGas.hpp"
#include <array>

/*!
 * \class CNasaGas
 * \brief Child class for defining the NASA polynomial gas model.
 */
class CNasaGas : public CIdealGas {
 protected:
  array<su2double, 7> CpLowPolyCoefficientsND{{0.0}};  /*!< \brief NASA low temp polynomial coefficients. */
  array<su2double, 7> CpHighPolyCoefficientsND{{0.0}}; /*!< \brief NASA high temp polynomial coefficients. */
  su2double TransitionTemperatureND{0.0};              /*!< \brief NASA transition temperature. */

 public:
  /*!
   * \brief Constructor of the class.
   */
  CNasaGas(su2double gamma, su2double R, bool CompEntropy = true);

  /*!
   * \brief Set specific heat Cp model.
   */
  void SetCpModel(const CConfig* config, su2double val_Temperature_Ref) override;

  /*!
   * \brief Set the Dimensionless State using Density and Internal Energy
   */
  void SetTDState_rhoe(su2double rho, su2double e) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Temperature
   */
  void SetTDState_PT(su2double P, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Density
   */
  void SetTDState_Prho(su2double P, su2double rho) override;

  /*!
   * \brief Set the Dimensionless State using Density and Temperature
   */
  void SetTDState_rhoT(su2double rho, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Entropy
   */
  void SetTDState_Ps(su2double P, su2double s) override;

 protected:
  /*!
   * \brief Calculate Cp, Enthalpy and Entropy from Temperature.
   */
  void ComputeProperties_T(su2double T);
};
