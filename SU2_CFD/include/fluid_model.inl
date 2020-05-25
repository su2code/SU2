/*!
 * \file fluid_model.inl
 * \brief In-Line subroutines of the <i>fluid_model.hpp</i> file.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 7.0.4 "Blackbird"
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

inline void CFluidModel::SetCpModel(CConfig *config) { }
inline su2double CFluidModel::GetPressure () const { return Pressure; }
inline su2double CFluidModel::GetSoundSpeed () const { return sqrt(SoundSpeed2); }
inline su2double CFluidModel::GetSoundSpeed2 () const { return SoundSpeed2; }
inline su2double CFluidModel::GetDensity () const { return Density; }
inline su2double CFluidModel::GetEntropy () const { return Entropy; }
inline su2double CFluidModel::GetStaticEnergy () const { return StaticEnergy; }
inline su2double CFluidModel::GetTemperature () const { return Temperature; }
inline su2double CFluidModel::GetCp () const { return Cp; }
inline su2double CFluidModel::GetCv () const { return Cv; }
inline su2double CFluidModel::GetdPdrho_e () const { return dPdrho_e; }
inline su2double CFluidModel::GetdPde_rho () const { return dPde_rho; }
inline su2double CFluidModel::GetdTdrho_e () const { return dTdrho_e; }
inline su2double CFluidModel::GetdTde_rho () const { return dTde_rho; }
inline su2double CFluidModel::Getdhdrho_P () const {return dhdrho_P;}
inline su2double CFluidModel::GetdhdP_rho () const {return dhdP_rho;}
inline su2double CFluidModel::Getdsdrho_P () const {return dsdrho_P;}
inline su2double CFluidModel::GetdsdP_rho () const {return dsdP_rho;}

inline su2double CFluidModel::GetLaminarViscosity () {
        LaminarViscosity->SetViscosity(Temperature, Density);
        Mu = LaminarViscosity->GetViscosity();
        LaminarViscosity->SetDerViscosity(Temperature, Density);
        dmudrho_T= LaminarViscosity->Getdmudrho_T(); 
        dmudT_rho= LaminarViscosity->GetdmudT_rho();  
        return Mu;
}

inline su2double CFluidModel::Getdmudrho_T () {
        return LaminarViscosity->Getdmudrho_T();
}

inline su2double CFluidModel::GetdmudT_rho () {
        return LaminarViscosity->GetdmudT_rho();
}

inline su2double CFluidModel::GetThermalConductivity () {
        ThermalConductivity->SetConductivity(Temperature, Density, Mu, Mu_Turb, Cp);
        Kt = ThermalConductivity->GetConductivity();
        ThermalConductivity->SetDerConductivity(Temperature, Density, dmudrho_T, dmudT_rho, Cp);
        dktdrho_T= ThermalConductivity->Getdktdrho_T(); 
        dktdT_rho= ThermalConductivity->GetdktdT_rho();
        return Kt;
}


inline su2double CFluidModel::Getdktdrho_T () const {
        return dktdrho_T;
}

inline su2double CFluidModel::GetdktdT_rho () const {
        return dktdT_rho;
}

inline void CFluidModel::SetTDState_rhoe (su2double rho, su2double e ) { }
inline void CFluidModel::SetTDState_PT (su2double P, su2double T ) { }
inline void CFluidModel::SetTDState_Prho (su2double P, su2double rho ) { }
inline void CFluidModel::SetTDState_hs (su2double h, su2double s ) { }
inline void CFluidModel::SetTDState_rhoT (su2double rho, su2double T ) { }
inline void CFluidModel::SetEnergy_Prho (su2double P, su2double rho ) { }
inline void CFluidModel::SetTDState_Ps (su2double P, su2double s ) { }
inline void CFluidModel::ComputeDerivativeNRBC_Prho (su2double P, su2double rho ){ }
inline void CFluidModel::SetTDState_T (su2double val_Temperature) { }
inline void CFluidModel::SetEddyViscosity (su2double val_Mu_Turb) { Mu_Turb = val_Mu_Turb; }
