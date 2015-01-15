/*!
 * \file fluid_model.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 3.2.7.1 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (fpalacios@stanford.edu).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
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

inline double CFluidModel::GetPressure () { return Pressure; }
inline double CFluidModel::GetSoundSpeed () { return sqrt(SoundSpeed2); }
inline double CFluidModel::GetSoundSpeed2 () { return SoundSpeed2; }
inline double CFluidModel::GetDensity () { return Density; }
inline double CFluidModel::GetEntropy () { return Entropy; }
inline double CFluidModel::GetStaticEnergy () { return StaticEnergy; }
inline double CFluidModel::GetTemperature () { return Temperature; }
inline double CFluidModel::GetCp () { return Cp; }
inline double CFluidModel::GetdPdrho_e () { return dPdrho_e; }
inline double CFluidModel::GetdPde_rho () { return dPde_rho; }
inline double CFluidModel::GetdTdrho_e () { return dTdrho_e; }
inline double CFluidModel::GetdTde_rho () { return dTde_rho; }

inline double CFluidModel::GetLaminarViscosity () {
        LaminarViscosity->SetViscosity(Temperature, Density);
        Mu = LaminarViscosity->GetViscosity();
        LaminarViscosity->SetDerViscosity(Temperature, Density);
        dmudrho_T= LaminarViscosity->Getdmudrho_T(); 
        dmudT_rho= LaminarViscosity->GetdmudT_rho();  
        return Mu;
}

inline double CFluidModel::Getdmudrho_T () {
        return LaminarViscosity->Getdmudrho_T();
}

inline double CFluidModel::GetdmudT_rho () {
        return LaminarViscosity->GetdmudT_rho();
}

inline double CFluidModel::GetThermalConductivity () {
        ThermalConductivity->SetConductivity(Temperature, Density, Mu, Cp);
        Kt = ThermalConductivity->GetConductivity();
        ThermalConductivity->SetDerConductivity(Temperature, Density, dmudrho_T, dmudT_rho, Cp);
        dktdrho_T= ThermalConductivity->Getdktdrho_T(); 
        dktdT_rho= ThermalConductivity->GetdktdT_rho();
        return Kt;
}


inline double CFluidModel::Getdktdrho_T () {
        return dktdrho_T;
}

inline double CFluidModel::GetdktdT_rho () {
        return dktdT_rho;
}

inline void CFluidModel::SetTDState_rhoe (double rho, double e ) { }
inline void CFluidModel::SetTDState_PT (double P, double T ) { }
inline void CFluidModel::SetTDState_Prho (double P, double rho ) { }
inline void CFluidModel::SetTDState_hs (double h, double s ) { }
inline void CFluidModel::SetTDState_rhoT (double rho, double T ) { }
inline void CFluidModel::SetEnergy_Prho (double P, double rho ) { }
