/*!
 * \file CCoolProp.cpp
 * \brief Source of the fluid model from CoolProp.
 * \author P. Yan, G. Gori, A. Guardone
 * \version 7.4.0 "Blackbird"
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

#include "../../include/fluid/CCoolProp.hpp"

#ifdef USE_COOLPROP
#include "CoolProp.h"
#include "AbstractState.h"

CCoolProp::CCoolProp(string fluidname) : CFluidModel() {
  fluid_entity = std::unique_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS",fluidname));
  Gas_Constant = fluid_entity->gas_constant()/fluid_entity->molar_mass();
  Pressure_Critical = fluid_entity->p_critical();
  Temperature_Critical = fluid_entity->T_critical();
  acentric_factor = fluid_entity->acentric_factor();

  su2double p = 5e4;
  su2double T = 250;
  su2double p_max = 2e6;
  su2double T_max = 600;
  su2double delta_p = (p_max - p)/500;
  su2double delta_T = (T_max - T)/500;
  cout << "rho, e, s, dsde_rho, dsdrho_e, d2sde2, d2sdedrho, d2sdrho2, pressure, temperature, soundspeed2" << endl;
  for(unsigned long i=0;i<500;i++){
    p = 5e4;
    for(unsigned long j=0; j<500;j++){
      fluid_entity->update(CoolProp::PT_INPUTS, p, T);
      cout << fluid_entity->rhomass() << ", "
           << fluid_entity->umass() << ", "
           << fluid_entity->smass() << ", "
           << fluid_entity->first_partial_deriv(CoolProp::iSmass, CoolProp::iUmass, CoolProp::iDmass) << ", "
           << fluid_entity->first_partial_deriv(CoolProp::iSmass, CoolProp::iDmass, CoolProp::iUmass) << ", "
           << fluid_entity->second_partial_deriv(CoolProp::iSmass, CoolProp::iUmass, CoolProp::iDmass, CoolProp::iUmass, CoolProp::iDmass) << ", "
           << fluid_entity->second_partial_deriv(CoolProp::iSmass, CoolProp::iUmass, CoolProp::iDmass, CoolProp::iDmass, CoolProp::iUmass) << ", "
           << fluid_entity->second_partial_deriv(CoolProp::iSmass, CoolProp::iDmass, CoolProp::iUmass, CoolProp::iDmass, CoolProp::iUmass) << ", "
           << p << ", "
           << T << ", "
           << pow(fluid_entity->speed_sound(), 2) << endl;
      p += delta_p;
    }
    T += delta_T;

  }
}

CCoolProp::~CCoolProp() {}

void CCoolProp::SetTDState_rhoe(su2double rho, su2double e) {
  //cout<<"p "<<Pressure<<"Pc "<<Pressure_Critical<<"T "<<Temperature<<"Tc"<<Temperature_Critical<<endl;
  Density = rho;
  StaticEnergy = e;
  fluid_entity->update(CoolProp::DmassUmass_INPUTS, Density, StaticEnergy);
  Cp = fluid_entity->cpmass();
  Cv = fluid_entity->cvmass();
  Gamma = Cp / Cv;
  Pressure = fluid_entity->p();
  Temperature = fluid_entity->T();
  Entropy = fluid_entity->smass();
  dPdrho_e = fluid_entity->first_partial_deriv(CoolProp::iP, CoolProp::iDmass, CoolProp::iUmass);
  dPde_rho = fluid_entity->first_partial_deriv(CoolProp::iP, CoolProp::iUmass, CoolProp::iDmass);
  dTdrho_e = fluid_entity->first_partial_deriv(CoolProp::iT, CoolProp::iDmass, CoolProp::iUmass);
  dTde_rho = fluid_entity->first_partial_deriv(CoolProp::iT, CoolProp::iUmass, CoolProp::iDmass);
  if (fluid_entity->phase() == 6) {
      fluid_entity->specify_phase(CoolProp::iphase_gas);
      SetTDState_PT(Pressure,Temperature);
  }
  else{
      SoundSpeed2 = pow(fluid_entity->speed_sound(), 2);
  }
}

void CCoolProp::SetTDState_PT(su2double P, su2double T) {
  fluid_entity->update(CoolProp::PT_INPUTS, P, T);
  su2double rho = fluid_entity->rhomass();
  su2double e = fluid_entity->umass();
  SetTDState_rhoe(rho, e);
}

void CCoolProp::SetTDState_Prho(su2double P, su2double rho) {
  fluid_entity->update(CoolProp::DmassP_INPUTS, rho, P);
  su2double e = fluid_entity->umass();
  SetTDState_rhoe(rho, e);
}

void CCoolProp::SetEnergy_Prho(su2double P, su2double rho) {
  fluid_entity->update(CoolProp::DmassP_INPUTS, rho, P);
  StaticEnergy = fluid_entity->umass();
}

void CCoolProp::SetTDState_hs(su2double h, su2double s) {
  fluid_entity->update(CoolProp::HmassSmass_INPUTS, h, s);
  su2double rho = fluid_entity->rhomass();
  su2double e = fluid_entity->umass();
  SetTDState_rhoe(rho, e);
}

void  CCoolProp::SetTDState_Ps(su2double P, su2double s) {
  fluid_entity->update(CoolProp::PSmass_INPUTS, P, s);
  su2double Rho = fluid_entity->rhomass();
  su2double e = fluid_entity->umass();
  SetTDState_rhoe(Rho, e);
}

void CCoolProp::SetTDState_rhoT(su2double rho, su2double T) {
  fluid_entity->update(CoolProp::DmassT_INPUTS, rho, T);
  su2double Rho = fluid_entity->rhomass();
  su2double e = fluid_entity->umass();
  SetTDState_rhoe(Rho, e);
}

void  CCoolProp::ComputeDerivativeNRBC_Prho(su2double P, su2double rho) {
  SetTDState_Prho(P, rho);
  dhdrho_P = fluid_entity->first_partial_deriv(CoolProp::iHmass,CoolProp::iDmass,CoolProp::iP);
  dhdP_rho = fluid_entity->first_partial_deriv(CoolProp::iHmass,CoolProp::iP,CoolProp::iDmass);
  dsdP_rho = fluid_entity->first_partial_deriv(CoolProp::iSmass,CoolProp::iP,CoolProp::iDmass);
  dsdrho_P = fluid_entity->first_partial_deriv(CoolProp::iSmass,CoolProp::iDmass,CoolProp::iP);
}

#else
CCoolProp::CCoolProp(string fluidname) {
  SU2_MPI::Error("SU2 was not compiled with CoolProp (-Denable-coolprop=true). Note that CoolProp cannot be used with directdiff or autodiff", CURRENT_FUNCTION);
}
#endif