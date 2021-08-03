#pragma once

#include <vector>
#include <memory>

#include "CFluidModel.hpp"

class CFluidScalar : public CFluidModel {
// class CFluidScalar final : public CFluidModel {

protected:
  unsigned short n_species_mixture = 0.0;
  su2double Gas_Constant = 0.0; 
  su2double Gamma = 0.0; 
  
  bool wilke;
  bool davidson; 

  std::vector<su2double> massFractions; 
  std::vector<su2double> moleFractions;  
  std::vector<su2double> molarMasses;
  std::vector<su2double> laminarViscosity;
  std::vector<su2double> specificHeat; 
  std::vector<su2double> laminarthermalConductivity; 
  
  unique_ptr<CViscosityModel> LaminarViscosityPointers[100];  // How to fix this?
  unique_ptr<CConductivityModel> ThermalConductivityPointers[100];  

 public:

  CFluidScalar(CConfig *config, su2double value_pressure_operating);

  ~CFluidScalar() {}; // Remove. 

  void SetLaminarViscosityModel(const CConfig* config);
  // void SetLaminarViscosityModel(const CConfig* config) override;

  void SetThermalConductivityModel(const CConfig* config);
  // void SetThermalConductivityModel(const CConfig* config) override;

  unsigned long SetTDState_T(su2double val_temperature, su2double *val_scalars);
  // unsigned long SetTDState_T(su2double val_temperature, su2double *val_scalars) override;

  std::vector<su2double> massToMoleFractions(su2double* val_scalars); // Make private. 

  su2double wilkeViscosity(su2double* val_scalars); // Make private. 

  su2double davidsonViscosity(su2double* val_scalars); // Make private. 

  su2double wilkeConductivity(su2double *val_scalars); // Make private. 

  inline su2double GetLaminarViscosity() {return Mu; }
  // inline su2double GetLaminarViscosity() const {return Mu; }

  inline su2double GetThermalConductivity() { return Kt; }
  // inline su2double GetThermalConductivity() const { return Kt; }
};
