#pragma once

#include <vector>
#include <memory>

#include "CFluidModel.hpp"

class CFluidScalar final : public CFluidModel {
private:
  unsigned short n_scalars = 0;                       /*!< \brief Number of transported scalars. */
  unsigned short n_species_mixture = 0;               /*!< \brief Number of species in mixture. */
  su2double Gas_Constant = 0.0;                       /*!< \brief Specific gas constant. */
  su2double Gamma = 0.0;                              /*!< \brief Ratio of specific heats of the gas. */
  su2double Pressure_Thermodynamic = 0.0;             /*!< \brief Constant pressure thermodynamic. */

  bool wilke;
  bool davidson;

  std::vector<su2double> massFractions;               /*!< \brief Mass fractions of all species. */
  std::vector<su2double> moleFractions;               /*!< \brief Mole fractions of all species. */
  std::vector<su2double> molarMasses;                 /*!< \brief Molar masses of all species. */
  std::vector<su2double> laminarViscosity;            /*!< \brief Laminar viscosity of all species. */
  std::vector<su2double> specificHeat;                /*!< \brief Specific heat of all species. */
  std::vector<su2double> laminarthermalConductivity;  /*!< \brief Laminar thermal conductivity of all species. */

  static const int ARRAYSIZE = 100;
  std::unique_ptr<CViscosityModel> LaminarViscosityPointers[ARRAYSIZE];
  std::unique_ptr<CConductivityModel> ThermalConductivityPointers[ARRAYSIZE];

  /*!
   * \brief Convert mass fractions to mole fractions.
   * \param[in] val_scalars - Scalar mass fraction.
   */
  std::vector<su2double>& massToMoleFractions(const su2double * const val_scalars);

  /*!
   * \brief Wilke mixing law for mixture viscosity.
   * \param[in] val_scalars - Scalar mass fraction.
   */
  su2double wilkeViscosity(const su2double * const val_scalars);

  /*!
   * \brief Davidson mixing law for mixture viscosity.
   * \param[in] val_scalars - Scalar mass fraction.
   */
  su2double davidsonViscosity(const su2double * const val_scalars);

  /*!
   * \brief Wilke mixing law for mixture thermal conductivity.
   * \param[in] val_scalars - Scalar mass fraction.
   */
  su2double wilkeConductivity(const su2double * const val_scalars);

public:
  /*!
   * \brief Constructor of the class.
   */
  CFluidScalar(CConfig *config, const su2double value_pressure_operating);

  /*!
   * \brief Set viscosity model.
   */
  void SetLaminarViscosityModel(const CConfig* config) override;

  /*!
   * \brief Set thermal conductivity model.
   */
  void SetThermalConductivityModel(const CConfig* config) override;

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] val_temperature - Temperature value at the point.
   * \param[in] val_scalars - Scalar mass fraction.
   */
  unsigned long SetTDState_T(const su2double val_temperature, su2double * const val_scalars) override;

  /*!
   * \brief Get fluid dynamic viscosity.
   */
  inline su2double GetLaminarViscosity() override {return Mu; }

  /*!
   * \brief Get fluid thermal conductivity.
   */
  inline su2double GetThermalConductivity() override { return Kt; }
};
