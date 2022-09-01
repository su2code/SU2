#pragma once

#include "../../Common/include/containers/CLookUpTable.hpp"
#include "CFluidModel.hpp"


class CFluidFlamelet : public CFluidModel {

protected:

  int rank;

  unsigned short n_scalars;
  unsigned short n_lookups;
  unsigned short n_table_sources; 

  vector<string> table_scalar_names;    /*!< \brief vector to store names of scalar variables.   */
  vector<string> table_source_names;    /*!< \brief vector to store names of scalar source variables.   */
  vector<string> table_lookup_names;    /*!< \brief vector to store names of look up variables.   */

  su2double mass_diffusivity;
  //su2double source_energy;
  su2double dDensitydPV;
  su2double dSourcePVdPV;
  su2double dDensitydEnth;

  vector<su2double> source_scalar;
  vector<su2double> lookup_scalar;

  CLookUpTable *look_up_table;

 public:
  CFluidFlamelet(CConfig *config, su2double value_pressure_operating);

  ~CFluidFlamelet();

  unsigned long SetTDState_T(su2double val_temperature, su2double *val_scalars);

  unsigned long SetScalarSources(su2double *val_scalars);

  unsigned long SetScalarLookups(su2double *val_scalars);

  void SetTDState_prog_enth(su2double val_prog, su2double val_enth);

  unsigned long GetEnthFromTemp(su2double *enthalpy, su2double val_prog, su2double val_temp);

  inline CLookUpTable* GetLookUpTable() {return look_up_table; }

  //inline su2double GetSourceEnergy() { return source_energy; }

  inline su2double GetMassDiffusivity() { return mass_diffusivity; }

  inline su2double GetThermalConductivity() { return Kt; }

  inline su2double GetLaminarViscosity() { return Mu; }

  inline pair<su2double, su2double> GetTableLimitsEnth() { return look_up_table->GetTableLimitsEnth(); }

  inline pair<su2double, su2double> GetTableLimitsProg() { return look_up_table->GetTableLimitsProg(); }

  inline su2double GetdDensitydPV() { return dDensitydPV; }

  inline su2double GetdSourcePVdPV() { return dSourcePVdPV; }

  inline su2double GetdDensitydEnth() { return dDensitydEnth; }

  inline su2double GetScalarSources(int i_scalar) { return source_scalar.at(i_scalar); }

  inline su2double *GetScalarSources(){ return &source_scalar[0]; }

  inline unsigned short GetNScalars() { return n_scalars; }

  inline su2double GetScalarLookups(int i_scalar) { return lookup_scalar.at(i_scalar); }

};
