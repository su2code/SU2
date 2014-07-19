/*!
 * \file gas_model.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
 * \author S.Vitale, M.Pini, G.Gori, A.Guardone, P.Colonna
 * \version 3.2.0 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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
inline double CFluidModel::GetDpDd_e () { return DpDd_e; }
inline double CFluidModel::GetDpDe_d () { return DpDe_d; }
inline double CFluidModel::GetLaminarViscosity (double T, double rho) {
        DynamicViscosity->SetViscosity(T, rho);
        return DynamicViscosity->GetViscosity();
}
inline void CFluidModel::SetTDState_rhoe (double rho, double e ) { }
inline void CFluidModel::SetTDState_PT (double P, double T ) { }
inline void CFluidModel::SetTDState_Prho (double P, double rho ) { }
inline void CFluidModel::SetTDState_hs (double h, double s ) { }
//inline void CFluidModel::SetTDState_Ps (double P, double s ) { }
inline void CFluidModel::SetEnergy_Prho (double P, double rho ) { }
