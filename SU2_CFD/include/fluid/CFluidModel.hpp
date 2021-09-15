/*!
 * \file CFluidModel.hpp
 * \brief Defines the main fluid model class for thermophysical properties.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna, T. Economon
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include <cmath>
#include <memory>

#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/basic_types/datatype_structure.hpp"
#include "CConductivityModel.hpp"
#include "CViscosityModel.hpp"
#include "CDiffusivityModel.hpp"
#include "../numerics/CLookUpTable.hpp"

using namespace std;

/*!
 * \class CFluidModel
 * \brief Main class for defining the Thermo-Physical Model
 * \author: S.Vitale, G.Gori, M.Pini, T. Economon
 */
class CFluidModel {
 protected:
  su2double StaticEnergy{0.0}; /*!< \brief Internal Energy. */
  su2double Entropy{0.0};      /*!< \brief Entropy. */
  su2double Density{0.0};      /*!< \brief Density. */
  su2double Pressure{0.0};     /*!< \brief Pressure. */
  su2double SoundSpeed2{0.0};  /*!< \brief SpeedSound. */
  su2double Temperature{0.0};  /*!< \brief Temperature. */
  su2double dPdrho_e{0.0};     /*!< \brief DpDd_e. */
  su2double dPde_rho{0.0};     /*!< \brief DpDe_d. */
  su2double dTdrho_e{0.0};     /*!< \brief DTDd_e. */
  su2double dTde_rho{0.0};     /*!< \brief DTDe_d. */
  su2double dhdrho_P{0.0};     /*!< \brief DhDrho_p. */
  su2double dhdP_rho{0.0};     /*!< \brief DhDp_rho. */
  su2double dsdrho_P{0.0};     /*!< \brief DsDrho_p. */
  su2double dsdP_rho{0.0};     /*!< \brief DsDp_rho. */
  su2double Cp{0.0};           /*!< \brief Specific Heat Capacity at constant pressure. */
  su2double Cv{0.0};           /*!< \brief Specific Heat Capacity at constant volume. */
  su2double Mu{0.0};           /*!< \brief Laminar viscosity. */
  su2double Mu_Turb{0.0};      /*!< \brief Eddy viscosity provided by a turbulence model. */
  su2double dmudrho_T{0.0};    /*!< \brief Partial derivative of viscosity w.r.t. density. */
  su2double dmudT_rho{0.0};    /*!< \brief Partial derivative of viscosity w.r.t. temperature. */
  su2double Kt{0.0};           /*!< \brief Thermal conductivity. */
  su2double dktdrho_T{0.0};    /*!< \brief Partial derivative of conductivity w.r.t. density. */
  su2double dktdT_rho{0.0};    /*!< \brief Partial derivative of conductivity w.r.t. temperature. */
  su2double mass_diffusivity{0.0};

  unique_ptr<CViscosityModel> LaminarViscosity;        /*!< \brief Laminar Viscosity Model */
  unique_ptr<CConductivityModel> ThermalConductivity;  /*!< \brief Thermal Conductivity Model */
  unique_ptr<CDiffusivityModel> MassDiffusivity;       /*!< \brief Mass Diffusivity Model */

  CLookUpTable* look_up_table;

 public:
  virtual ~CFluidModel() {}

  /*!
   * \brief Get fluid pressure.
   */
  su2double GetPressure() const { return Pressure; }

  /*!
   * \brief Get fluid temperature.
   */
  su2double GetTemperature() const { return Temperature; }

  /*!
   * \brief Get fluid entropy.
   */
  su2double GetEntropy() const { return Entropy; }

  /*!
   * \brief Get fluid internal energy.
   */
  su2double GetStaticEnergy() const { return StaticEnergy; }

  /*!
   * \brief Get fluid density.
   */
  su2double GetDensity() const { return Density; }

  /*!
   * \brief Get fluid speed of sound.
   */
  su2double GetSoundSpeed() const { return sqrt(SoundSpeed2); }

  /*!
   * \brief Get fluid speed of sound squared.
   */
  su2double GetSoundSpeed2() const { return SoundSpeed2; }

  /*!
   * \brief Get fluid specific heat at constant pressure.
   */
  su2double GetCp() const { return Cp; }

  /*!
   * \brief Get fluid specific heat at constant volume.
   */
  su2double GetCv() const { return Cv; }

  /*!
   * \brief Get fluid mean specific heat capacity at constant pressure.
   */
  su2double ComputeMeanSpecificHeatCp(const std::vector<su2double> specific_heat_cp, const su2double * const val_scalars) {
    su2double val_scalars_sum = 0.0;
    unsigned short n_scalars = specific_heat_cp.size() - 1;

    for (int i_scalar = 0; i_scalar < n_scalars; i_scalar++){
      Cp += specific_heat_cp[i_scalar] * val_scalars[i_scalar];
      val_scalars_sum += val_scalars[i_scalar];
    }
    return Cp += specific_heat_cp[n_scalars]*(1 - val_scalars_sum);
  }

  /*!
   * \brief Get fluid mean molecular weight.
   */
  static su2double ComputeMeanMolecularWeight(const std::vector<su2double> molar_masses, const su2double * const val_scalars) {
    su2double OneOverMeanMolecularWeight = 0.0;
    su2double val_scalars_sum = 0.0;
    unsigned short n_scalars = molar_masses.size() - 1;

    for (int i_scalar = 0; i_scalar < n_scalars; i_scalar++){
      OneOverMeanMolecularWeight += val_scalars[i_scalar]/(molar_masses[i_scalar]/1000);
      val_scalars_sum += val_scalars[i_scalar];
    }
    OneOverMeanMolecularWeight += (1 - val_scalars_sum)/(molar_masses[n_scalars]/1000);
    return 1/OneOverMeanMolecularWeight;
  }

  /*!
   * \brief Get the source term of the transported scalar
   */
  virtual inline su2double* GetScalarSources(){ return nullptr; }

    /*!
   * \brief Get the source term of the transported scalar
   * \param[in] val_ix - Index of the scalar.
   */
  virtual inline su2double GetScalarSources(int val_ix){ return 37; }
  
  /*!
  * \brief Get the number of transported scalars
  */
  virtual inline unsigned short GetNScalars() {return 73; }

  /*!
   * \brief Get the looked up scalar field for combustion
   */
  virtual inline su2double GetScalarLookups(int){ return 37; }


  virtual CLookUpTable* GetLookUpTable() {return look_up_table; }

  virtual inline unsigned long GetEnthFromTemp(su2double *enthalpy, 
                                               su2double  val_prog, 
                                               su2double  val_temp) { return 73; }

  virtual inline pair<su2double, su2double> GetTableLimitsEnth() { return make_pair(73,37); }
  virtual inline pair<su2double, su2double> GetTableLimitsProg() { return make_pair(73,37); }

  virtual inline su2double GetdDensitydPV() { return 0.0; }
  virtual inline su2double GetdSourcePVdPV() { return 0.0; }
  virtual inline su2double GetdDensitydEnth() { return 0.0; }

  /*!
   * \brief Get fluid dynamic viscosity.
   */
  virtual inline su2double GetLaminarViscosity() {
      LaminarViscosity->SetViscosity(Temperature, Density);
      Mu = LaminarViscosity->GetViscosity();
      LaminarViscosity->SetDerViscosity(Temperature, Density);
      dmudrho_T = LaminarViscosity->Getdmudrho_T();
      dmudT_rho = LaminarViscosity->GetdmudT_rho();
    return Mu;
  }

  /*!
   * \brief Get fluid thermal conductivity.
   */
  virtual inline su2double GetThermalConductivity() {
    ThermalConductivity->SetConductivity(Temperature, Density, Mu, Mu_Turb, Cp);
    Kt = ThermalConductivity->GetConductivity();
    ThermalConductivity->SetDerConductivity(Temperature, Density, dmudrho_T, dmudT_rho, Cp);
    dktdrho_T = ThermalConductivity->Getdktdrho_T();
    dktdT_rho = ThermalConductivity->GetdktdT_rho();
    return Kt;
  }

/*!
   * \brief Get fluid mass diffusivity.
   */
  virtual inline su2double GetMassDiffusivity() {
    MassDiffusivity->SetDiffusivity(Temperature, Density, Mu, Mu_Turb, Cp, Kt);
    mass_diffusivity = MassDiffusivity->GetDiffusivity();
    //ThermalConductivity->SetDerConductivity(Temperature, Density, dmudrho_T, dmudT_rho, Cp);
    //dktdrho_T = ThermalConductivity->Getdktdrho_T();
    //dktdT_rho = ThermalConductivity->GetdktdT_rho();
    return mass_diffusivity;
  }
  
  /*!
   * \brief Get fluid pressure partial derivative.
   */
  su2double GetdPdrho_e() const { return dPdrho_e; }

  /*!
   * \brief Get fluid pressure partial derivative.
   */
  su2double GetdPde_rho() const { return dPde_rho; }

  /*!
   * \brief Get fluid temperature partial derivative.
   */
  su2double GetdTdrho_e() const { return dTdrho_e; }

  /*!
   * \brief Get fluid temperature partial derivative.
   */
  su2double GetdTde_rho() const { return dTde_rho; }

  /*!
   * \brief Get fluid enthalpy partial derivative.
   */
  su2double Getdhdrho_P() const { return dhdrho_P; }

  /*!
   * \brief Get fluid enthalpy partial derivative.
   */
  su2double GetdhdP_rho() const { return dhdP_rho; }

  /*!
   * \brief Get fluid entropy partial derivative.
   */
  su2double Getdsdrho_P() const { return dsdrho_P; }

  /*!
   * \brief Get fluid entropy partial derivative.
   */
  su2double GetdsdP_rho() const { return dsdP_rho; }

  /*!
   * \brief Get fluid dynamic viscosity partial derivative.
   */
  su2double Getdmudrho_T() { return LaminarViscosity->Getdmudrho_T(); }

  /*!
   * \brief Get fluid dynamic viscosity partial derivative.
   */
  su2double GetdmudT_rho() { return LaminarViscosity->GetdmudT_rho(); }

  /*!
   * \brief Get fluid thermal conductivity partial derivative.
   */
  su2double Getdktdrho_T() const { return dktdrho_T; }

  /*!
   * \brief Get fluid thermal conductivity partial derivative.
   */
  su2double GetdktdT_rho() const { return dktdT_rho; }

  /*!
   * \brief Set specific heat Cp model.
   */
  virtual void SetCpModel(const CConfig* config) {}

  /*!
   * \brief Set viscosity model.
   */
  virtual void SetLaminarViscosityModel(const CConfig* config);

  /*!
   * \brief Set thermal conductivity model.
   */
  virtual void SetThermalConductivityModel(const CConfig* config);

  /*!
   * \brief Set mass diffusivity model.
   */
  void SetMassDiffusivityModel(const CConfig* config);

  /*!
   * \brief Set the state needed for the mass diffusivity model.
   */
  /* nijso: why do we need this? put in constructor or delete */
  void SetDiffusivityState(su2double val_temperature, su2double val_Density,
                           su2double val_Mu, su2double val_Mu_Turb, su2double val_Cp) {
  Temperature = val_temperature;
  Density = val_Density;
  Mu = val_Mu;
  Mu_Turb = val_Mu_Turb;
  Cp = val_Cp;
}

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("e, rho").
   * \param[in] rho - first thermodynamic variable.
   * \param[in] e - second thermodynamic variable.
   */
  virtual void SetTDState_rhoe(su2double rho, su2double e) {}

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("PT").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (T).
   */
  virtual void SetTDState_PT(su2double P, su2double T) {}

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (v).
   */
  virtual void SetTDState_Prho(su2double P, su2double rho) {}

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (v).
   *
   */
  virtual void SetEnergy_Prho(su2double P, su2double rho) {}

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("hs").
   * \param[in] th1 - first thermodynamic variable (h).
   * \param[in] th2 - second thermodynamic variable (s).
   *
   */
  virtual void SetTDState_hs(su2double h, su2double s) {}

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
   * \param[in] th1 - first thermodynamic variable (rho).
   * \param[in] th2 - second thermodynamic variable (T).
   *
   */
  virtual void SetTDState_rhoT(su2double rho, su2double T) {}

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (s).
   */
  virtual void SetTDState_Ps(su2double P, su2double s) {}

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (v).
   *
   */
  virtual void ComputeDerivativeNRBC_Prho(su2double P, su2double rho) {}

  /*!
   * \brief Virtual member.
   * \param[in] T - Temperature value at the point.
   */
  virtual unsigned long SetTDState_T(su2double val_Temperature, su2double *val_scalars = nullptr) {return 0;}

  virtual unsigned long SetScalarSources(su2double *val_scalars) {return 0;}

  virtual unsigned long SetScalarLookups(su2double *val_scalars) {return 0;}

  /*!
   * \brief Set fluid eddy viscosity provided by a turbulence model needed for computing effective thermal conductivity.
   */
  void SetEddyViscosity(su2double val_Mu_Turb) { Mu_Turb = val_Mu_Turb; }

};
