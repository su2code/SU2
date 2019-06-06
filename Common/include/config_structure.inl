/*!
 * \file config_structure.inl
 * \brief In-Line subroutines of the <i>config_structure.hpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

inline su2double CConfig::GetHTP_Axis(unsigned short val_index) { return HTP_Axis[val_index]; }

inline su2double CConfig::GetCFL_AdaptParam(unsigned short val_index) { return CFL_AdaptParam[val_index]; }

inline bool CConfig::GetCFL_Adapt(void) { return CFL_Adapt; }

inline bool CConfig::GetHB_Precondition(void) { return HB_Precondition; }

inline void CConfig::SetInflow_Mach(unsigned short val_imarker, su2double val_fanface_mach) { Inflow_Mach[val_imarker] = val_fanface_mach; }

inline void CConfig::SetInflow_Pressure(unsigned short val_imarker, su2double val_fanface_pressure) { Inflow_Pressure[val_imarker] = val_fanface_pressure; }

inline void CConfig::SetInflow_MassFlow(unsigned short val_imarker, su2double val_fanface_massflow) { Inflow_MassFlow[val_imarker] = val_fanface_massflow; }

inline void CConfig::SetInflow_ReverseMassFlow(unsigned short val_imarker, su2double val_fanface_reversemassflow) { Inflow_ReverseMassFlow[val_imarker] = val_fanface_reversemassflow; }

inline void CConfig::SetInflow_TotalPressure(unsigned short val_imarker, su2double val_fanface_totalpressure) { Inflow_TotalPressure[val_imarker] = val_fanface_totalpressure; }

inline void CConfig::SetInflow_Temperature(unsigned short val_imarker, su2double val_fanface_temperature) { Inflow_Temperature[val_imarker] = val_fanface_temperature; }

inline void CConfig::SetInflow_TotalTemperature(unsigned short val_imarker, su2double val_fanface_totaltemperature) { Inflow_TotalTemperature[val_imarker] = val_fanface_totaltemperature; }

inline void CConfig::SetInflow_RamDrag(unsigned short val_imarker, su2double val_fanface_ramdrag) { Inflow_RamDrag[val_imarker] = val_fanface_ramdrag; }

inline void CConfig::SetInflow_Force(unsigned short val_imarker, su2double val_fanface_force) { Inflow_Force[val_imarker] = val_fanface_force; }

inline void CConfig::SetInflow_Power(unsigned short val_imarker, su2double val_fanface_power) { Inflow_Power[val_imarker] = val_fanface_power; }

inline su2double CConfig::GetInflow_ReverseMassFlow(unsigned short val_imarker) { return Inflow_ReverseMassFlow[val_imarker]; }

inline void CConfig::SetExhaust_Pressure(unsigned short val_imarker, su2double val_exhaust_pressure) { Exhaust_Pressure[val_imarker] = val_exhaust_pressure; }

inline void CConfig::SetExhaust_Temperature(unsigned short val_imarker, su2double val_exhaust_temp) { Exhaust_Temperature[val_imarker] = val_exhaust_temp; }

inline void CConfig::SetExhaust_MassFlow(unsigned short val_imarker, su2double val_exhaust_massflow) { Exhaust_MassFlow[val_imarker] = val_exhaust_massflow; }

inline void CConfig::SetExhaust_TotalPressure(unsigned short val_imarker, su2double val_exhaust_totalpressure) { Exhaust_TotalPressure[val_imarker] = val_exhaust_totalpressure; }

inline void CConfig::SetExhaust_TotalTemperature(unsigned short val_imarker, su2double val_exhaust_totaltemp) { Exhaust_TotalTemperature[val_imarker] = val_exhaust_totaltemp; }

inline void CConfig::SetExhaust_GrossThrust(unsigned short val_imarker, su2double val_exhaust_grossthrust) { Exhaust_GrossThrust[val_imarker] = val_exhaust_grossthrust; }

inline void CConfig::SetExhaust_Force(unsigned short val_imarker, su2double val_exhaust_force) { Exhaust_Force[val_imarker] = val_exhaust_force; }

inline void CConfig::SetExhaust_Power(unsigned short val_imarker, su2double val_exhaust_power) { Exhaust_Power[val_imarker] = val_exhaust_power; }

inline void CConfig::SetEngine_Mach(unsigned short val_imarker, su2double val_engine_mach) { Engine_Mach[val_imarker] = val_engine_mach; }

inline void CConfig::SetEngine_Force(unsigned short val_imarker, su2double val_engine_force) { Engine_Force[val_imarker] = val_engine_force; }

inline void CConfig::SetEngine_Power(unsigned short val_imarker, su2double val_engine_power) { Engine_Power[val_imarker] = val_engine_power; }

inline void CConfig::SetEngine_NetThrust(unsigned short val_imarker, su2double val_engine_netthrust) { Engine_NetThrust[val_imarker] = val_engine_netthrust; }

inline void CConfig::SetEngine_GrossThrust(unsigned short val_imarker, su2double val_engine_grossthrust) { Engine_GrossThrust[val_imarker] = val_engine_grossthrust; }

inline void CConfig::SetEngine_Area(unsigned short val_imarker, su2double val_engine_area) { Engine_Area[val_imarker] = val_engine_area; }

inline void CConfig::SetActDisk_DeltaPress(unsigned short val_imarker, su2double val_actdisk_deltapress) { ActDisk_DeltaPress[val_imarker] = val_actdisk_deltapress; }

inline void CConfig::SetActDisk_Power(unsigned short val_imarker, su2double val_actdisk_power) { ActDisk_Power[val_imarker] = val_actdisk_power; }

inline void CConfig::SetActDisk_MassFlow(unsigned short val_imarker, su2double val_actdisk_massflow) { ActDisk_MassFlow[val_imarker] = val_actdisk_massflow; }

inline void CConfig::SetActDisk_Mach(unsigned short val_imarker, su2double val_actdisk_mach) { ActDisk_Mach[val_imarker] = val_actdisk_mach; }

inline void CConfig::SetActDisk_Force(unsigned short val_imarker, su2double val_actdisk_force) { ActDisk_Force[val_imarker] = val_actdisk_force; }

inline void CConfig::SetOutlet_MassFlow(unsigned short val_imarker, su2double val_massflow) { Outlet_MassFlow[val_imarker] = val_massflow; }

inline void CConfig::SetOutlet_Density(unsigned short val_imarker, su2double val_density) { Outlet_Density[val_imarker] = val_density; }

inline void CConfig::SetOutlet_Area(unsigned short val_imarker, su2double val_area) { Outlet_Area[val_imarker] = val_area; }

inline void CConfig::SetSurface_DC60(unsigned short val_imarker, su2double val_surface_distortion) { Surface_DC60[val_imarker] = val_surface_distortion; }

inline void CConfig::SetSurface_MassFlow(unsigned short val_imarker, su2double val_surface_massflow) { Surface_MassFlow[val_imarker] = val_surface_massflow; }

inline void CConfig::SetSurface_Mach(unsigned short val_imarker, su2double val_surface_mach) { Surface_Mach[val_imarker] = val_surface_mach; }

inline void CConfig::SetSurface_Temperature(unsigned short val_imarker, su2double val_surface_temperature) { Surface_Temperature[val_imarker] = val_surface_temperature; }

inline void CConfig::SetSurface_Pressure(unsigned short val_imarker, su2double val_surface_pressure) { Surface_Pressure[val_imarker] = val_surface_pressure; }

inline void CConfig::SetSurface_Density(unsigned short val_imarker, su2double val_surface_density) { Surface_Density[val_imarker] = val_surface_density; }

inline void CConfig::SetSurface_Enthalpy(unsigned short val_imarker, su2double val_surface_enthalpy) { Surface_Enthalpy[val_imarker] = val_surface_enthalpy; }

inline void CConfig::SetSurface_NormalVelocity(unsigned short val_imarker, su2double val_surface_normalvelocity) { Surface_NormalVelocity[val_imarker] = val_surface_normalvelocity; }

inline void CConfig::SetSurface_Uniformity(unsigned short val_imarker, su2double val_surface_streamwiseuniformity) { Surface_Uniformity[val_imarker] = val_surface_streamwiseuniformity; }

inline void CConfig::SetSurface_SecondaryStrength(unsigned short val_imarker, su2double val_surface_secondarystrength) { Surface_SecondaryStrength[val_imarker] = val_surface_secondarystrength; }

inline void CConfig::SetSurface_SecondOverUniform(unsigned short val_imarker, su2double val_surface_secondaryoverstream) { Surface_SecondOverUniform[val_imarker] = val_surface_secondaryoverstream; }

inline void CConfig::SetSurface_MomentumDistortion(unsigned short val_imarker, su2double val_surface_momentumdistortion) { Surface_MomentumDistortion[val_imarker] = val_surface_momentumdistortion; }

inline void CConfig::SetSurface_TotalTemperature(unsigned short val_imarker, su2double val_surface_totaltemperature) { Surface_TotalTemperature[val_imarker] = val_surface_totaltemperature; }

inline void CConfig::SetSurface_TotalPressure(unsigned short val_imarker, su2double val_surface_totalpressure) { Surface_TotalPressure[val_imarker] = val_surface_totalpressure; }

inline void CConfig::SetSurface_PressureDrop(unsigned short val_imarker, su2double val_surface_pressuredrop) { Surface_PressureDrop[val_imarker] = val_surface_pressuredrop; }

inline void CConfig::SetSurface_IDC(unsigned short val_imarker, su2double val_surface_distortion) { Surface_IDC[val_imarker] = val_surface_distortion; }

inline void CConfig::SetSurface_IDC_Mach(unsigned short val_imarker, su2double val_surface_distortion) { Surface_IDC_Mach[val_imarker] = val_surface_distortion; }

inline void CConfig::SetSurface_IDR(unsigned short val_imarker, su2double val_surface_distortion) { Surface_IDR[val_imarker] = val_surface_distortion; }

inline void CConfig::SetActDisk_DeltaTemp(unsigned short val_imarker, su2double val_actdisk_deltatemp) { ActDisk_DeltaTemp[val_imarker] = val_actdisk_deltatemp; }

inline void CConfig::SetActDisk_TotalPressRatio(unsigned short val_imarker, su2double val_actdisk_pressratio) { ActDisk_TotalPressRatio[val_imarker] = val_actdisk_pressratio; }

inline void CConfig::SetActDisk_TotalTempRatio(unsigned short val_imarker, su2double val_actdisk_tempratio) { ActDisk_TotalTempRatio[val_imarker] = val_actdisk_tempratio; }

inline void CConfig::SetActDisk_StaticPressRatio(unsigned short val_imarker, su2double val_actdisk_pressratio) { ActDisk_StaticPressRatio[val_imarker] = val_actdisk_pressratio; }

inline void CConfig::SetActDisk_StaticTempRatio(unsigned short val_imarker, su2double val_actdisk_tempratio) { ActDisk_StaticTempRatio[val_imarker] = val_actdisk_tempratio; }

inline void CConfig::SetActDisk_NetThrust(unsigned short val_imarker, su2double val_actdisk_netthrust) { ActDisk_NetThrust[val_imarker] = val_actdisk_netthrust; }

inline void CConfig::SetActDisk_BCThrust(unsigned short val_imarker, su2double val_actdisk_bcthrust) { ActDisk_BCThrust[val_imarker] = val_actdisk_bcthrust; }

inline void CConfig::SetActDisk_BCThrust_Old(unsigned short val_imarker, su2double val_actdisk_bcthrust_old) { ActDisk_BCThrust_Old[val_imarker] = val_actdisk_bcthrust_old; }

inline void CConfig::SetActDisk_GrossThrust(unsigned short val_imarker, su2double val_actdisk_grossthrust) { ActDisk_GrossThrust[val_imarker] = val_actdisk_grossthrust; }

inline void CConfig::SetActDisk_Area(unsigned short val_imarker, su2double val_actdisk_area) { ActDisk_Area[val_imarker] = val_actdisk_area; }

inline void CConfig::SetActDiskInlet_ReverseMassFlow(unsigned short val_imarker, su2double val_actdisk_area) { ActDisk_ReverseMassFlow[val_imarker] = val_actdisk_area; }

inline su2double CConfig::GetActDisk_DeltaPress(unsigned short val_imarker) { return ActDisk_DeltaPress[val_imarker]; }

inline su2double CConfig::GetActDisk_DeltaTemp(unsigned short val_imarker) { return ActDisk_DeltaTemp[val_imarker]; }

inline su2double CConfig::GetActDisk_TotalPressRatio(unsigned short val_imarker) { return ActDisk_TotalPressRatio[val_imarker]; }

inline su2double CConfig::GetActDisk_TotalTempRatio(unsigned short val_imarker) { return ActDisk_TotalTempRatio[val_imarker]; }

inline su2double CConfig::GetActDisk_StaticPressRatio(unsigned short val_imarker) { return ActDisk_StaticPressRatio[val_imarker]; }

inline su2double CConfig::GetActDisk_StaticTempRatio(unsigned short val_imarker) { return ActDisk_StaticTempRatio[val_imarker]; }

inline su2double CConfig::GetActDisk_Power(unsigned short val_imarker) { return ActDisk_Power[val_imarker]; }

inline su2double CConfig::GetActDisk_MassFlow(unsigned short val_imarker) { return ActDisk_MassFlow[val_imarker]; }

inline su2double CConfig::GetActDisk_Mach(unsigned short val_imarker) { return ActDisk_Mach[val_imarker]; }

inline su2double CConfig::GetActDisk_Force(unsigned short val_imarker) { return ActDisk_Force[val_imarker]; }

inline su2double CConfig::GetSurface_MassFlow(unsigned short val_imarker) { return Surface_MassFlow[val_imarker]; }

inline su2double CConfig::GetSurface_Mach(unsigned short val_imarker) { return Surface_Mach[val_imarker]; }

inline su2double CConfig::GetSurface_Temperature(unsigned short val_imarker) { return Surface_Temperature[val_imarker]; }

inline su2double CConfig::GetSurface_Pressure(unsigned short val_imarker) { return Surface_Pressure[val_imarker]; }

inline su2double CConfig::GetSurface_Density(unsigned short val_imarker) { return Surface_Density[val_imarker]; }

inline su2double CConfig::GetSurface_Enthalpy(unsigned short val_imarker) { return Surface_Enthalpy[val_imarker]; }

inline su2double CConfig::GetSurface_NormalVelocity(unsigned short val_imarker) { return Surface_NormalVelocity[val_imarker]; }

inline su2double CConfig::GetSurface_Uniformity(unsigned short val_imarker) { return Surface_Uniformity[val_imarker]; }

inline su2double CConfig::GetSurface_SecondaryStrength(unsigned short val_imarker) { return Surface_SecondaryStrength[val_imarker]; }

inline su2double CConfig::GetSurface_SecondOverUniform(unsigned short val_imarker) { return Surface_SecondOverUniform[val_imarker]; }

inline su2double CConfig::GetSurface_MomentumDistortion(unsigned short val_imarker) { return Surface_MomentumDistortion[val_imarker]; }

inline su2double CConfig::GetSurface_TotalTemperature(unsigned short val_imarker) { return Surface_TotalTemperature[val_imarker]; }

inline su2double CConfig::GetSurface_TotalPressure(unsigned short val_imarker) { return Surface_TotalPressure[val_imarker]; }

inline su2double CConfig::GetSurface_PressureDrop(unsigned short val_imarker) { return Surface_PressureDrop[val_imarker]; }

inline su2double CConfig::GetSurface_DC60(unsigned short val_imarker) { return Surface_DC60[val_imarker]; }

inline su2double CConfig::GetSurface_IDC(unsigned short val_imarker) { return Surface_IDC[val_imarker]; }

inline su2double CConfig::GetSurface_IDC_Mach(unsigned short val_imarker) { return Surface_IDC_Mach[val_imarker]; }

inline su2double CConfig::GetSurface_IDR(unsigned short val_imarker) { return Surface_IDR[val_imarker]; }

inline su2double CConfig::GetActDisk_NetThrust(unsigned short val_imarker) { return ActDisk_NetThrust[val_imarker]; }

inline su2double CConfig::GetActDisk_BCThrust(unsigned short val_imarker) { return ActDisk_BCThrust[val_imarker]; }

inline su2double CConfig::GetActDisk_BCThrust_Old(unsigned short val_imarker) { return ActDisk_BCThrust_Old[val_imarker]; }

inline su2double CConfig::GetActDisk_GrossThrust(unsigned short val_imarker) { return ActDisk_GrossThrust[val_imarker]; }

inline su2double CConfig::GetActDisk_Area(unsigned short val_imarker) { return ActDisk_Area[val_imarker]; }

inline su2double CConfig::GetActDisk_ReverseMassFlow(unsigned short val_imarker) { return ActDisk_ReverseMassFlow[val_imarker]; }

inline void CConfig::SetActDiskInlet_Pressure(unsigned short val_imarker, su2double val_actdisk_press) { ActDiskInlet_Pressure[val_imarker] = val_actdisk_press; }

inline void CConfig::SetActDiskInlet_TotalPressure(unsigned short val_imarker, su2double val_actdisk_totalpress) { ActDiskInlet_TotalPressure[val_imarker] = val_actdisk_totalpress; }

inline void CConfig::SetActDiskInlet_RamDrag(unsigned short val_imarker, su2double val_actdisk_ramdrag) { ActDiskInlet_RamDrag[val_imarker] = val_actdisk_ramdrag; }

inline void CConfig::SetActDiskInlet_Force(unsigned short val_imarker, su2double val_actdisk_force) { ActDiskInlet_Force[val_imarker] = val_actdisk_force; }

inline void CConfig::SetActDiskInlet_Power(unsigned short val_imarker, su2double val_actdisk_power) { ActDiskInlet_Power[val_imarker] = val_actdisk_power; }

inline void CConfig::SetActDiskInlet_Temperature(unsigned short val_imarker, su2double val_actdisk_temp) { ActDiskInlet_Temperature[val_imarker] = val_actdisk_temp; }

inline void CConfig::SetActDiskInlet_TotalTemperature(unsigned short val_imarker, su2double val_actdisk_totaltemp) { ActDiskInlet_TotalTemperature[val_imarker] = val_actdisk_totaltemp; }

inline void CConfig::SetActDiskInlet_MassFlow(unsigned short val_imarker, su2double val_actdisk_massflow) { ActDiskInlet_MassFlow[val_imarker] = val_actdisk_massflow; }

inline void CConfig::SetActDiskOutlet_Pressure(unsigned short val_imarker, su2double val_actdisk_press) { ActDiskOutlet_Pressure[val_imarker] = val_actdisk_press; }

inline void CConfig::SetActDiskOutlet_TotalPressure(unsigned short val_imarker, su2double val_actdisk_totalpress) { ActDiskOutlet_TotalPressure[val_imarker] = val_actdisk_totalpress; }

inline void CConfig::SetActDiskOutlet_GrossThrust(unsigned short val_imarker, su2double val_actdisk_grossthrust) { ActDiskOutlet_GrossThrust[val_imarker] = val_actdisk_grossthrust; }

inline void CConfig::SetActDiskOutlet_Force(unsigned short val_imarker, su2double val_actdisk_force) { ActDiskOutlet_Force[val_imarker] = val_actdisk_force; }

inline void CConfig::SetActDiskOutlet_Power(unsigned short val_imarker, su2double val_actdisk_power) { ActDiskOutlet_Power[val_imarker] = val_actdisk_power; }

inline void CConfig::SetActDiskOutlet_Temperature(unsigned short val_imarker, su2double val_actdisk_temp) { ActDiskOutlet_Temperature[val_imarker] = val_actdisk_temp; }

inline void CConfig::SetActDiskOutlet_TotalTemperature(unsigned short val_imarker, su2double val_actdisk_totaltemp) { ActDiskOutlet_TotalTemperature[val_imarker] = val_actdisk_totaltemp; }

inline void CConfig::SetActDiskOutlet_MassFlow(unsigned short val_imarker, su2double val_actdisk_massflow) { ActDiskOutlet_MassFlow[val_imarker] = val_actdisk_massflow; }

inline su2double CConfig::GetEngine_Mach(unsigned short val_imarker) { return Engine_Mach[val_imarker]; }

inline su2double CConfig::GetEngine_Force(unsigned short val_imarker) { return Engine_Force[val_imarker]; }

inline su2double CConfig::GetEngine_Power(unsigned short val_imarker) { return Engine_Power[val_imarker]; }

inline su2double CConfig::GetEngine_NetThrust(unsigned short val_imarker) { return Engine_NetThrust[val_imarker]; }

inline su2double CConfig::GetEngine_GrossThrust(unsigned short val_imarker) { return Engine_GrossThrust[val_imarker]; }

inline su2double CConfig::GetEngine_Area(unsigned short val_imarker) { return Engine_Area[val_imarker]; }

inline unsigned short CConfig::GetnZone(void) { return nZone; }

inline unsigned short CConfig::GetiZone(void) { return iZone; }

inline unsigned short CConfig::GetKind_SU2(void) { return Kind_SU2; }

inline unsigned short CConfig::GetRef_NonDim(void) { return Ref_NonDim; }

inline unsigned short CConfig::GetRef_Inc_NonDim(void) { return Ref_Inc_NonDim; }

inline void CConfig::SetKind_SU2(unsigned short val_kind_su2) { Kind_SU2 = val_kind_su2 ; }

inline bool CConfig::GetContinuous_Adjoint(void) { return ContinuousAdjoint; }

inline bool CConfig::GetViscous(void) { return Viscous; }

inline unsigned long CConfig::GetnExtIter(void) { return nExtIter; }

inline unsigned short CConfig::GetnTimeInstances(void) { return nTimeInstances; }

inline su2double CConfig::GetHarmonicBalance_Period(void) { return HarmonicBalance_Period; }

inline void CConfig::SetExtIter(unsigned long val_iter) { ExtIter = val_iter; }

inline void CConfig::SetExtIter_OffSet(unsigned long val_iter) { ExtIter_OffSet = val_iter; }

inline void CConfig::SetOuterIter(unsigned long val_iter) { OuterIter = val_iter; }

inline void CConfig::SetIntIter(unsigned long val_iter) { IntIter = val_iter; }

inline unsigned long CConfig::GetExtIter(void) { return ExtIter; }

inline unsigned long CConfig::GetExtIter_OffSet(void) { return ExtIter_OffSet; }

inline unsigned long CConfig::GetOuterIter(void) { return OuterIter; }

inline unsigned long CConfig::GetIntIter(void) { return IntIter; }

inline unsigned long CConfig::GetUnst_nIntIter(void) { return Unst_nIntIter; }

inline unsigned long CConfig::GetDyn_nIntIter(void) { return Dyn_nIntIter; }

inline long CConfig::GetUnst_RestartIter(void) { return Unst_RestartIter; }

inline long CConfig::GetUnst_AdjointIter(void) { return Unst_AdjointIter; }

inline void CConfig::SetPhysicalTime(su2double val_t) { PhysicalTime = val_t; }

inline su2double CConfig::GetPhysicalTime(void) { return PhysicalTime; }

inline bool CConfig::GetReorientElements(void) { return ReorientElements; }

inline unsigned long CConfig::GetIter_Avg_Objective(void) { return Iter_Avg_Objective ; }

inline long CConfig::GetDyn_RestartIter(void) { return Dyn_RestartIter; }

inline string CConfig::GetPlaneTag(unsigned short index) { return PlaneTag[index]; }

inline su2double CConfig::GetEA_IntLimit(unsigned short index) { return EA_IntLimit[index]; }

inline su2double CConfig::GetEA_ScaleFactor(void) { return EA_ScaleFactor; }

inline su2double CConfig::GetAdjointLimit(void) { return AdjointLimit; }

inline su2double *CConfig::GetHold_GridFixed_Coord(void) { return Hold_GridFixed_Coord; }

inline su2double *CConfig::GetSubsonicEngine_Cyl(void) { return SubsonicEngine_Cyl; }

inline su2double *CConfig::GetSubsonicEngine_Values(void) { return SubsonicEngine_Values; }

inline su2double *CConfig::GetDistortionRack(void) { return DistortionRack; }

inline unsigned short CConfig::GetAnalytical_Surface(void) { return Analytical_Surface; }

inline unsigned short CConfig::GetGeo_Description(void) { return Geo_Description; }

inline su2double CConfig::GetDualVol_Power(void) { return DualVol_Power; }

inline bool CConfig::GetExtraOutput(void) { return ExtraOutput; }

inline long CConfig::GetExtraHeatOutputZone(void) { return ExtraHeatOutputZone; }

inline su2double CConfig::GetRefArea(void) { return RefArea; }

inline su2double CConfig::GetThermalDiffusivity(void) { return Thermal_Diffusivity; }

inline su2double CConfig::GetThermalDiffusivity_Solid(void) { return Thermal_Diffusivity_Solid; }

inline su2double CConfig::GetTemperature_Freestream_Solid(void) { return Temperature_Freestream_Solid;  }

inline su2double CConfig::GetElasticyMod(unsigned short id_val) { return ElasticityMod[id_val]; }

inline bool CConfig::GetDE_Effects(void) { return DE_Effects; }

inline unsigned short CConfig::GetnElectric_Constant(void) { return nElectric_Constant; }

inline su2double CConfig::GetElectric_Constant(unsigned short iVar) { return Electric_Constant[iVar]; }

inline su2double CConfig::GetKnowles_B(void) { return Knowles_B; }

inline su2double CConfig::GetKnowles_N(void) { return Knowles_N; }

inline unsigned short CConfig::GetDV_FEA(void) { return Kind_DV_FEA; }

inline unsigned long CConfig::GetRefNode_ID(void) { return refNodeID; }

inline su2double CConfig::GetRefNode_Displacement(unsigned short val_coeff) { return RefNode_Displacement[val_coeff]; }

inline su2double CConfig::GetRefNode_Penalty(void) { return RefNode_Penalty; }

inline bool CConfig::GetRefGeom(void) { return RefGeom; }

inline string CConfig::GetRefGeom_FEMFileName(void) { return RefGeom_FEMFileName; }

inline unsigned short CConfig::GetRefGeom_FileFormat(void) { return RefGeom_FileFormat; }

inline unsigned short CConfig::GetElas2D_Formulation(void) { return Kind_2DElasForm; }

inline su2double CConfig::GetPoissonRatio(unsigned short id_val) { return PoissonRatio[id_val]; }

inline su2double CConfig::GetMaterialDensity(unsigned short id_val) { return MaterialDensity[id_val]; }

inline unsigned short CConfig::GetnElasticityMod(void) { return nElasticityMod; }

inline unsigned short CConfig::GetnPoissonRatio(void) { return nPoissonRatio; }

inline unsigned short CConfig::GetnMaterialDensity(void) { return nMaterialDensity; }

inline unsigned short CConfig::GetMaterialCompressibility(void) { return Kind_Material_Compress; }

inline unsigned short CConfig::GetMaterialModel(void) { return Kind_Material; }

inline unsigned short CConfig::GetGeometricConditions(void) { return Kind_Struct_Solver; }

inline bool CConfig::GetPrestretch(void) { return Prestretch; }

inline bool CConfig::Add_CrossTerm(void) { return addCrossTerm; }

inline void CConfig::Set_CrossTerm(bool needCrossTerm) { addCrossTerm = needCrossTerm; }

inline string CConfig::GetFEA_FileName(void) { return FEA_FileName; }

inline string CConfig::GetPrestretch_FEMFileName(void) { return Prestretch_FEMFileName; }

inline su2double CConfig::GetRefLength(void) { return RefLength; }

inline su2double CConfig::GetRefElemLength(void) { return RefElemLength; }

inline su2double CConfig::GetRefSharpEdges(void) { return RefSharpEdges; }

inline su2double CConfig::GetDomainVolume(void) { return DomainVolume; }

inline void CConfig::SetRefArea(su2double val_area) { RefArea = val_area; }

inline void CConfig::SetSemiSpan(su2double val_semispan) { SemiSpan = val_semispan; }

inline void CConfig::SetDomainVolume(su2double val_volume) { DomainVolume = val_volume; }

inline void CConfig::SetnExtIter(unsigned long val_niter) { nExtIter = val_niter; }

inline su2double CConfig::GetMach(void) { return Mach; }

inline su2double CConfig::GetGamma(void) { return Gamma; }

inline su2double CConfig::GetStations_Bounds(unsigned short val_var) { return Stations_Bounds[val_var]; }

inline su2double CConfig::GetFFD_Axis(unsigned short val_var) { return FFD_Axis[val_var]; }

inline su2double CConfig::GetBulk_Modulus(void) { return Bulk_Modulus; }

inline su2double CConfig::GetBeta_Factor(void) { return Beta_Factor; }

inline su2double CConfig::GetGas_Constant(void) { return Gas_Constant; }

inline su2double CConfig::GetGas_ConstantND(void) { return Gas_ConstantND; }

inline su2double CConfig::GetMolecular_Weight(void) { return Molecular_Weight; }

inline su2double CConfig::GetSpecific_Heat_Cp(void) { return Specific_Heat_Cp; }

inline su2double CConfig::GetSpecific_Heat_Cp_Solid(void) { return Specific_Heat_Cp_Solid; }

inline su2double CConfig::GetSpecific_Heat_CpND(void) { return Specific_Heat_CpND; }

inline su2double CConfig::GetSpecific_Heat_Cv(void) { return Specific_Heat_Cv; }

inline su2double CConfig::GetSpecific_Heat_CvND(void) { return Specific_Heat_CvND; }

inline su2double CConfig::GetThermal_Expansion_Coeff(void) { return Thermal_Expansion_Coeff; }

inline su2double CConfig::GetThermal_Expansion_CoeffND(void) { return Thermal_Expansion_CoeffND; }

inline su2double CConfig::GetInc_Density_Ref(void) { return Inc_Density_Ref; }

inline su2double CConfig::GetInc_Velocity_Ref(void) { return Inc_Velocity_Ref; }

inline su2double CConfig::GetInc_Temperature_Ref(void) { return Inc_Temperature_Ref; }

inline su2double CConfig::GetInc_Density_Init(void) { return Inc_Density_Init; }

inline su2double* CConfig::GetInc_Velocity_Init(void) { return Inc_Velocity_Init; }

inline su2double CConfig::GetInc_Temperature_Init(void) { return Inc_Temperature_Init; }

inline su2double CConfig::GetHeat_Flux_Ref(void) { return Heat_Flux_Ref; }

inline su2double CConfig::GetWallTemperature(void) { return Wall_Temperature; }

inline su2double CConfig::GetGas_Constant_Ref(void) { return Gas_Constant_Ref; }

inline su2double CConfig::GetTemperature_FreeStream(void) { return Temperature_FreeStream; }

inline su2double CConfig::GetEnergy_FreeStream(void) { return Energy_FreeStream; }

inline su2double CConfig::GetViscosity_FreeStream(void) { return Viscosity_FreeStream; }

inline su2double CConfig::GetDensity_FreeStream(void) { return Density_FreeStream; }

inline su2double CConfig::GetDensity_Solid(void) { return Density_Solid; }

inline su2double CConfig::GetModVel_FreeStream(void) { return ModVel_FreeStream; }

inline su2double CConfig::GetModVel_FreeStreamND(void) { return ModVel_FreeStreamND; }

inline su2double CConfig::GetPressure_FreeStream(void) { return Pressure_FreeStream; }

inline su2double CConfig::GetPressure_Thermodynamic(void) { return Pressure_Thermodynamic; }

inline su2double CConfig::GetTemperature_ve_FreeStream(void) { return Temperature_ve_FreeStream; }

inline su2double CConfig::GetPrandtl_Lam(void) { return Prandtl_Lam; }

inline su2double CConfig::GetPrandtl_Turb(void) { return Prandtl_Turb; }

inline su2double CConfig::GetThermalConductivity_Solid(void) { return Thermal_Conductivity_Solid; }

inline su2double CConfig::GetLength_Ref(void) { return Length_Ref; }

inline su2double CConfig::GetPressure_Ref(void) { return Pressure_Ref; }

inline su2double CConfig::GetTemperature_Ref(void) { return Temperature_Ref; }

inline su2double CConfig::GetDensity_Ref(void) { return Density_Ref; }

inline su2double CConfig::GetVelocity_Ref(void) { return Velocity_Ref; }

inline su2double CConfig::GetEnergy_Ref(void) { return Energy_Ref; }

inline su2double CConfig::GetTime_Ref(void) { return Time_Ref; }

inline su2double CConfig::GetViscosity_Ref(void) { return Viscosity_Ref; }

inline su2double CConfig::GetConductivity_Ref(void) { return Conductivity_Ref; }

inline su2double CConfig::GetHighlite_Area(void) { return Highlite_Area; }

inline su2double CConfig::GetFan_Poly_Eff(void) { return Fan_Poly_Eff; }

inline su2double CConfig::GetOmega_Ref(void) { return Omega_Ref; }

inline su2double CConfig::GetForce_Ref(void) { return Force_Ref; }

inline su2double CConfig::GetPressure_FreeStreamND(void) { return Pressure_FreeStreamND; }

inline su2double CConfig::GetPressure_ThermodynamicND(void) { return Pressure_ThermodynamicND; }

inline su2double CConfig::GetTemperature_FreeStreamND(void) { return Temperature_FreeStreamND; }

inline su2double CConfig::GetDensity_FreeStreamND(void) { return Density_FreeStreamND; }

inline su2double* CConfig::GetVelocity_FreeStreamND(void) { return Velocity_FreeStreamND; }

inline su2double* CConfig::GetVelocity_FreeStream(void) { return Velocity_FreeStream; }

inline su2double CConfig::GetEnergy_FreeStreamND(void) { return Energy_FreeStreamND; }

inline su2double CConfig::GetViscosity_FreeStreamND(void) { return Viscosity_FreeStreamND; }

inline su2double CConfig::GetTke_FreeStreamND(void) { return Tke_FreeStreamND; }

inline su2double CConfig::GetOmega_FreeStreamND(void) { return Omega_FreeStreamND; }

inline su2double CConfig::GetTke_FreeStream(void) { return Tke_FreeStream; }

inline su2double CConfig::GetOmega_FreeStream(void) { return Omega_FreeStream; }

inline su2double CConfig::GetNuFactor_FreeStream(void) { return NuFactor_FreeStream; }

inline su2double CConfig::GetNuFactor_Engine(void) { return NuFactor_Engine; }

inline su2double CConfig::GetSecondaryFlow_ActDisk(void) { return SecondaryFlow_ActDisk; }

inline su2double CConfig::GetInitial_BCThrust(void) { return Initial_BCThrust; }

inline void CConfig::SetInitial_BCThrust(su2double val_bcthrust) { Initial_BCThrust = val_bcthrust; }

inline su2double CConfig::GetIntermittency_FreeStream(void) { return Intermittency_FreeStream; }

inline su2double CConfig::GetTurbulenceIntensity_FreeStream(void) { return TurbulenceIntensity_FreeStream; }

inline su2double CConfig::GetTurb2LamViscRatio_FreeStream(void) { return Turb2LamViscRatio_FreeStream;}

inline su2double* CConfig::GetMassFrac_FreeStream(void) { return MassFrac_FreeStream; }

inline su2double CConfig::GetLength_Reynolds(void) { return Length_Reynolds; }

inline unsigned short CConfig::GetnStartUpIter(void) { return nStartUpIter; }

inline su2double *CConfig::GetRefOriginMoment(unsigned short val_marker) {
    if(val_marker < nMarker_Monitoring) {
      RefOriginMoment[0] = RefOriginMoment_X[val_marker];
      RefOriginMoment[1] = RefOriginMoment_Y[val_marker];
      RefOriginMoment[2] = RefOriginMoment_Z[val_marker];
    }
    return RefOriginMoment;
}

inline su2double CConfig::GetRefOriginMoment_X(unsigned short val_marker) { return RefOriginMoment_X[val_marker]; }

inline su2double CConfig::GetRefOriginMoment_Y(unsigned short val_marker) { return RefOriginMoment_Y[val_marker]; }

inline su2double CConfig::GetRefOriginMoment_Z(unsigned short val_marker) { return RefOriginMoment_Z[val_marker]; }

inline void CConfig::SetRefOriginMoment_X(unsigned short val_marker, su2double val_origin) { RefOriginMoment_X[val_marker] = val_origin; }

inline void CConfig::SetRefOriginMoment_Y(unsigned short val_marker, su2double val_origin) { RefOriginMoment_Y[val_marker] = val_origin; }

inline void CConfig::SetRefOriginMoment_Z(unsigned short val_marker, su2double val_origin) { RefOriginMoment_Z[val_marker] = val_origin; }

inline su2double CConfig::GetChargeCoeff(void) { return ChargeCoeff; }

inline su2double CConfig::GetVenkat_LimiterCoeff(void) { return Venkat_LimiterCoeff; }

inline unsigned long CConfig::GetLimiterIter(void) { return LimiterIter; }

inline su2double CConfig::GetAdjSharp_LimiterCoeff(void) { return AdjSharp_LimiterCoeff; }

inline su2double CConfig::GetReynolds(void) { return Reynolds; }

inline su2double CConfig::GetFroude(void) { return Froude; }

inline void CConfig::SetPressure_FreeStreamND(su2double val_pressure_freestreamnd) { Pressure_FreeStreamND = val_pressure_freestreamnd; }

inline void CConfig::SetPressure_FreeStream(su2double val_pressure_freestream) { Pressure_FreeStream = val_pressure_freestream; }

inline void CConfig::SetPressure_ThermodynamicND(su2double val_pressure_thermodynamicnd) { Pressure_ThermodynamicND = val_pressure_thermodynamicnd; }

inline void CConfig::SetPressure_Thermodynamic(su2double val_pressure_thermodynamic) { Pressure_Thermodynamic = val_pressure_thermodynamic; }

inline void CConfig::SetDensity_FreeStreamND(su2double val_density_freestreamnd) { Density_FreeStreamND = val_density_freestreamnd; }

inline void CConfig::SetDensity_FreeStream(su2double val_density_freestream) { Density_FreeStream = val_density_freestream; }

inline void CConfig::SetViscosity_FreeStream(su2double val_viscosity_freestream) { Viscosity_FreeStream = val_viscosity_freestream; }

inline void CConfig::SetModVel_FreeStream(su2double val_modvel_freestream) { ModVel_FreeStream = val_modvel_freestream; }

inline void CConfig::SetModVel_FreeStreamND(su2double val_modvel_freestreamnd) { ModVel_FreeStreamND = val_modvel_freestreamnd; }

inline void CConfig::SetTemperature_FreeStream(su2double val_temperature_freestream) { Temperature_FreeStream = val_temperature_freestream; }

inline void CConfig::SetTemperature_FreeStreamND(su2double val_temperature_freestreamnd) { Temperature_FreeStreamND = val_temperature_freestreamnd; }

inline void CConfig::SetGas_ConstantND(su2double val_gas_constantnd) { Gas_ConstantND = val_gas_constantnd; }

inline void CConfig::SetSpecific_Heat_CpND(su2double val_specific_heat_cpnd) { Specific_Heat_CpND = val_specific_heat_cpnd; }

inline void CConfig::SetSpecific_Heat_CvND(su2double val_specific_heat_cvnd) { Specific_Heat_CvND = val_specific_heat_cvnd; }

inline void CConfig::SetThermal_Expansion_CoeffND(su2double val_thermal_expansion_coeffnd) { Thermal_Expansion_CoeffND = val_thermal_expansion_coeffnd; }

inline void CConfig::SetVelocity_FreeStream(su2double val_velocity_freestream, unsigned short val_dim) { Velocity_FreeStream[val_dim] = val_velocity_freestream; }

inline void CConfig::SetVelocity_FreeStreamND(su2double val_velocity_freestreamnd, unsigned short val_dim) { Velocity_FreeStreamND[val_dim] = val_velocity_freestreamnd; }

inline void CConfig::SetViscosity_FreeStreamND(su2double val_viscosity_freestreamnd) { Viscosity_FreeStreamND = val_viscosity_freestreamnd; }

inline void CConfig::SetTke_FreeStreamND(su2double val_tke_freestreamnd) { Tke_FreeStreamND = val_tke_freestreamnd; }

inline void CConfig::SetOmega_FreeStreamND(su2double val_omega_freestreamnd) { Omega_FreeStreamND = val_omega_freestreamnd; }

inline void CConfig::SetTke_FreeStream(su2double val_tke_freestream) { Tke_FreeStream = val_tke_freestream; }

inline void CConfig::SetOmega_FreeStream(su2double val_omega_freestream) { Omega_FreeStream = val_omega_freestream; }

inline void CConfig::SetEnergy_FreeStreamND(su2double val_energy_freestreamnd) { Energy_FreeStreamND = val_energy_freestreamnd; }

inline void CConfig::SetEnergy_FreeStream(su2double val_energy_freestream) { Energy_FreeStream = val_energy_freestream; }

inline void CConfig::SetTotal_UnstTimeND(su2double val_total_unsttimend) { Total_UnstTimeND = val_total_unsttimend; }

inline void CConfig::SetFroude(su2double val_froude) { Froude = val_froude; }

inline void CConfig::SetReynolds(su2double val_reynolds) { Reynolds = val_reynolds; }

inline void CConfig::SetMach(su2double val_mach) { Mach = val_mach; }

inline void CConfig::SetLength_Ref(su2double val_length_ref) { Length_Ref = val_length_ref; }

inline void CConfig::SetVelocity_Ref(su2double val_velocity_ref) { Velocity_Ref = val_velocity_ref; }

inline void CConfig::SetPressure_Ref(su2double val_pressure_ref) { Pressure_Ref = val_pressure_ref; }

inline void CConfig::SetDensity_Ref(su2double val_density_ref) { Density_Ref = val_density_ref; }

inline void CConfig::SetTemperature_Ref(su2double val_temperature_ref) { Temperature_Ref = val_temperature_ref; }

inline void CConfig::SetTime_Ref(su2double val_time_ref) { Time_Ref = val_time_ref; }

inline void CConfig::SetOmega_Ref(su2double val_omega_ref) { Omega_Ref = val_omega_ref; }

inline void CConfig::SetForce_Ref(su2double val_force_ref) { Force_Ref = val_force_ref; }

inline void CConfig::SetGas_Constant_Ref(su2double val_gas_constant_ref) { Gas_Constant_Ref = val_gas_constant_ref; }

inline void CConfig::SetGas_Constant(su2double val_gas_constant) { Gas_Constant = val_gas_constant; }

inline void CConfig::SetSpecific_Heat_Cp(su2double val_specific_heat_cp) { Specific_Heat_Cp = val_specific_heat_cp; }

inline void CConfig::SetSpecific_Heat_Cv(su2double val_specific_heat_cv) { Specific_Heat_Cv = val_specific_heat_cv; }

inline void CConfig::SetThermal_Expansion_Coeff(su2double val_thermal_expansion_coeff) { Thermal_Expansion_Coeff = val_thermal_expansion_coeff; }

inline void CConfig::SetHeat_Flux_Ref(su2double val_heat_flux_ref) { Heat_Flux_Ref = val_heat_flux_ref; }

inline void CConfig::SetViscosity_Ref(su2double val_viscosity_ref) { Viscosity_Ref = val_viscosity_ref; }

inline void CConfig::SetConductivity_Ref(su2double val_conductivity_ref) { Conductivity_Ref = val_conductivity_ref; }

inline void CConfig::SetEnergy_Ref(su2double val_energy_ref) { Energy_Ref = val_energy_ref; }

inline void CConfig::SetThermalDiffusivity_Solid(su2double val_thermal_diffusivity) { Thermal_Diffusivity_Solid = val_thermal_diffusivity; }

inline su2double CConfig::GetAoA(void) { return AoA; }

inline void CConfig::SetAoA(su2double val_AoA) { AoA = val_AoA; }

inline void CConfig::SetAoA_Offset(su2double val_AoA_offset) { AoA_Offset = val_AoA_offset; }

inline void CConfig::SetAoS_Offset(su2double val_AoS_offset) { AoS_Offset = val_AoS_offset; }

inline void CConfig::SetAoA_Sens(su2double val_AoA_sens) { AoA_Sens = val_AoA_sens; }

inline void CConfig::SetAoS(su2double val_AoS) { AoS = val_AoS; }

inline su2double CConfig::GetAoS(void) { return AoS; }

inline su2double CConfig::GetAoA_Offset(void) { return AoA_Offset; }

inline su2double CConfig::GetAoS_Offset(void) { return AoS_Offset; }

inline su2double CConfig::GetAoA_Sens(void) { return AoA_Sens; }

inline unsigned short CConfig::GetnMGLevels(void) { return nMGLevels; }

inline void CConfig::SetMGLevels(unsigned short val_nMGLevels) { nMGLevels = val_nMGLevels; }

inline unsigned short CConfig::GetFinestMesh(void) { return FinestMesh; }

inline void CConfig::SetFinestMesh(unsigned short val_finestmesh) { FinestMesh = val_finestmesh; }

inline void CConfig::SubtractFinestMesh(void) { FinestMesh = FinestMesh-1; }

inline unsigned short CConfig::GetDesign_Variable(unsigned short val_dv) { return Design_Variable[val_dv]; }

inline bool CConfig::GetBuffet_Monitoring(void) { return Buffet_Monitoring; }

inline su2double CConfig::GetBuffet_k(void) { return Buffet_k; }

inline su2double CConfig::GetBuffet_lambda(void) { return Buffet_lambda; }

inline unsigned short CConfig::GetConvCriteria(void) { return ConvCriteria; }

inline unsigned short CConfig::GetMGCycle(void) { return MGCycle; }

inline unsigned short CConfig::GetGeometryMode(void) { return GeometryMode; }

inline su2double CConfig::GetCFL(unsigned short val_mesh) {	return CFL[val_mesh]; }

inline su2double CConfig::GetCFL_Solid(void) {	return CFLSolid; }

inline void CConfig::SetCFL(unsigned short val_mesh, su2double val_cfl) { CFL[val_mesh] = val_cfl; }

inline su2double CConfig::GetUnst_CFL(void) {	return Unst_CFL; }

inline su2double CConfig::GetMax_DeltaTime(void) {	return Max_DeltaTime; }

inline su2double CConfig::GetParamDV(unsigned short val_dv, unsigned short val_param) {	return ParamDV[val_dv][val_param]; }

inline su2double CConfig::GetCoordFFDBox(unsigned short val_ffd, unsigned short val_index) {	return CoordFFDBox[val_ffd][val_index]; }

inline unsigned short CConfig::GetDegreeFFDBox(unsigned short val_ffd, unsigned short val_index) {	return DegreeFFDBox[val_ffd][val_index]; }

inline string CConfig::GetFFDTag(unsigned short val_dv) {	return FFDTag[val_dv]; }

inline string CConfig::GetTagFFDBox(unsigned short val_ffd) {	return TagFFDBox[val_ffd]; }

inline unsigned short CConfig::GetnDV(void) {	return nDV; }

inline unsigned short CConfig::GetnDV_Value(unsigned short iDV) {	return nDV_Value[iDV]; }

inline unsigned short CConfig::GetnFFDBox(void) {	return nFFDBox; }

inline unsigned short CConfig::GetFFD_Continuity(void) { return FFD_Continuity; }

inline unsigned short CConfig::GetFFD_CoordSystem(void) { return FFD_CoordSystem; }

inline unsigned short CConfig::GetnRKStep(void) { return nRKStep; }

inline unsigned short CConfig::GetnLevels_TimeAccurateLTS(void) { return nLevels_TimeAccurateLTS; }

inline void CConfig::SetnLevels_TimeAccurateLTS(unsigned short val_nLevels) {nLevels_TimeAccurateLTS = val_nLevels;}

inline unsigned short CConfig::GetnTimeDOFsADER_DG(void) { return nTimeDOFsADER_DG; }

inline su2double *CConfig::GetTimeDOFsADER_DG(void) { return TimeDOFsADER_DG; }

inline unsigned short CConfig::GetnTimeIntegrationADER_DG(void) { return nTimeIntegrationADER_DG; }

inline su2double *CConfig::GetTimeIntegrationADER_DG(void) { return TimeIntegrationADER_DG; }

inline su2double *CConfig::GetWeightsIntegrationADER_DG(void) { return WeightsIntegrationADER_DG; }

inline su2double CConfig::Get_Alpha_RKStep(unsigned short val_step) { return RK_Alpha_Step[val_step]; }

inline su2double CConfig::GetLocationStations(unsigned short val_section) { return LocationStations[val_section]; }

inline su2double CConfig::GetNacelleLocation(unsigned short val_index) { return NacelleLocation[val_index]; }

inline unsigned short CConfig::GetnFFD_Fix_IDir(void) { return nFFD_Fix_IDir; }

inline unsigned short CConfig::GetnFFD_Fix_JDir(void) { return nFFD_Fix_JDir; }

inline unsigned short CConfig::GetnFFD_Fix_KDir(void) { return nFFD_Fix_KDir; }

inline short CConfig::GetFFD_Fix_IDir(unsigned short val_index) { return FFD_Fix_IDir[val_index]; }

inline short CConfig::GetFFD_Fix_JDir(unsigned short val_index) { return FFD_Fix_JDir[val_index]; }

inline short CConfig::GetFFD_Fix_KDir(unsigned short val_index) { return FFD_Fix_KDir[val_index]; }

inline unsigned short CConfig::GetMG_PreSmooth(unsigned short val_mesh) {	
	if (nMG_PreSmooth == 0) return 1;
	else return MG_PreSmooth[val_mesh]; 
}

inline unsigned short CConfig::GetMG_PostSmooth(unsigned short val_mesh) { 
	if (nMG_PostSmooth == 0) return 0;
	else return MG_PostSmooth[val_mesh];
}

inline unsigned short CConfig::GetMG_CorrecSmooth(unsigned short val_mesh) { 
	if (nMG_CorrecSmooth == 0) return 0;
	else return MG_CorrecSmooth[val_mesh]; 
}

inline unsigned long CConfig::GetWrt_Sol_Freq(void) { return Wrt_Sol_Freq; }

inline unsigned long CConfig::GetWrt_Sol_Freq_DualTime(void) { return Wrt_Sol_Freq_DualTime; }

inline unsigned long CConfig::GetWrt_Con_Freq(void) { return Wrt_Con_Freq; }

inline void CConfig::SetWrt_Con_Freq(unsigned long val_freq) { Wrt_Con_Freq = val_freq; }

inline unsigned long CConfig::GetWrt_Con_Freq_DualTime(void) { return Wrt_Con_Freq_DualTime; }

inline bool CConfig::GetWrt_Unsteady(void) { return Wrt_Unsteady; }

inline unsigned short CConfig::GetKind_Solver(void) { return Kind_Solver; }

inline void CConfig::SetKind_Solver(unsigned short val_solver) { Kind_Solver = val_solver; }

inline unsigned short CConfig::GetKind_MZSolver(void) { return Kind_MZSolver; }

inline unsigned short CConfig::GetKind_Regime(void) { return Kind_Regime; }

inline unsigned short CConfig::GetSystemMeasurements(void) { return SystemMeasurements; }

inline unsigned short CConfig::GetKind_GasModel(void) { return Kind_GasModel; }

inline unsigned short CConfig::GetKind_FluidModel(void) { return Kind_FluidModel; }

inline unsigned short CConfig::GetKind_FreeStreamOption(void) { return Kind_FreeStreamOption; } 

inline unsigned short CConfig::GetKind_DensityModel(void) { return Kind_DensityModel; } 

inline bool CConfig::GetEnergy_Equation(void) { return Energy_Equation; } 

inline unsigned short CConfig::GetKind_InitOption(void) { return Kind_InitOption; }

inline su2double CConfig::GetPressure_Critical(void) { return Pressure_Critical; }

inline su2double CConfig::GetTemperature_Critical(void) { return Temperature_Critical; }

inline su2double CConfig::GetAcentric_Factor(void) { return Acentric_Factor; }

inline unsigned short CConfig::GetKind_ViscosityModel(void) { return Kind_ViscosityModel; }

inline unsigned short CConfig::GetKind_ConductivityModel(void) { return Kind_ConductivityModel; }

inline unsigned short CConfig::GetKind_ConductivityModel_Turb(void) { return Kind_ConductivityModel_Turb; }

inline su2double CConfig::GetMu_Constant(void) { return Mu_Constant; }

inline su2double CConfig::GetMu_ConstantND(void) { return Mu_ConstantND; }

inline su2double CConfig::GetKt_Constant(void) { return Kt_Constant; }

inline su2double CConfig::GetKt_ConstantND(void) { return Kt_ConstantND; }

inline su2double CConfig::GetMu_Ref(void) { return Mu_Ref; }

inline su2double CConfig::GetMu_RefND(void) { return Mu_RefND; }

inline su2double CConfig::GetMu_Temperature_Ref(void) { return Mu_Temperature_Ref; }

inline su2double CConfig::GetMu_Temperature_RefND(void) { return Mu_Temperature_RefND; }

inline su2double CConfig::GetMu_S(void) { return Mu_S; }

inline su2double CConfig::GetMu_SND(void) { return Mu_SND; }

inline unsigned short CConfig::GetnPolyCoeffs(void) { return nPolyCoeffs; }

inline su2double CConfig::GetCp_PolyCoeff(unsigned short val_index) { return CpPolyCoefficients[val_index]; }

inline su2double CConfig::GetCp_PolyCoeffND(unsigned short val_index) { return CpPolyCoefficientsND[val_index]; }

inline su2double CConfig::GetMu_PolyCoeff(unsigned short val_index) { return MuPolyCoefficients[val_index]; }

inline su2double CConfig::GetMu_PolyCoeffND(unsigned short val_index) { return MuPolyCoefficientsND[val_index]; }

inline su2double* CConfig::GetMu_PolyCoeffND(void) { return MuPolyCoefficientsND; }

inline su2double CConfig::GetKt_PolyCoeff(unsigned short val_index) { return KtPolyCoefficients[val_index]; }

inline su2double CConfig::GetKt_PolyCoeffND(unsigned short val_index) { return KtPolyCoefficientsND[val_index]; }

inline su2double* CConfig::GetKt_PolyCoeffND(void) { return KtPolyCoefficientsND; }

inline void CConfig::SetMu_ConstantND(su2double mu_const) { Mu_ConstantND = mu_const; }

inline void CConfig::SetMu_RefND(su2double mu_ref) { Mu_RefND = mu_ref; }

inline void CConfig::SetMu_Temperature_RefND(su2double mu_Tref) {Mu_Temperature_RefND = mu_Tref; }

inline void CConfig::SetMu_SND(su2double mu_s) {Mu_SND = mu_s; }

inline void CConfig::SetKt_ConstantND(su2double kt_const) { Kt_ConstantND = kt_const; }

inline void CConfig::SetCp_PolyCoeffND(su2double val_coeff, unsigned short val_index) { CpPolyCoefficientsND[val_index] = val_coeff; }

inline void CConfig::SetMu_PolyCoeffND(su2double val_coeff, unsigned short val_index) { MuPolyCoefficientsND[val_index] = val_coeff; }

inline void CConfig::SetKt_PolyCoeffND(su2double val_coeff, unsigned short val_index) { KtPolyCoefficientsND[val_index] = val_coeff; }

inline unsigned short CConfig::GetKind_GridMovement(unsigned short val_iZone) { return Kind_GridMovement[val_iZone]; }

inline unsigned short CConfig::GetKind_GridMovement(void) { return Kind_GridMovement[0]; }

inline void CConfig::SetKind_GridMovement(unsigned short val_iZone, unsigned short motion_Type) { Kind_GridMovement[val_iZone] = motion_Type; }

inline su2double CConfig::GetMach_Motion(void) { return Mach_Motion; }

inline su2double CConfig::GetMotion_Origin_X(unsigned short val_iZone) { return Motion_Origin_X[val_iZone]; }

inline su2double CConfig::GetMotion_Origin_Y(unsigned short val_iZone) { return Motion_Origin_Y[val_iZone]; }

inline su2double CConfig::GetMotion_Origin_Z(unsigned short val_iZone) { return Motion_Origin_Z[val_iZone]; }

inline void CConfig::SetMotion_Origin_X(unsigned short val_iZone, su2double val_origin) { Motion_Origin_X[val_iZone] = val_origin; }

inline void CConfig::SetMotion_Origin_Y(unsigned short val_iZone, su2double val_origin) { Motion_Origin_Y[val_iZone] = val_origin; }

inline void CConfig::SetMotion_Origin_Z(unsigned short val_iZone, su2double val_origin) { Motion_Origin_Z[val_iZone] = val_origin; }

inline su2double CConfig::GetTranslation_Rate_X(unsigned short val_iZone) { return  Translation_Rate_X[val_iZone]; }

inline su2double CConfig::GetTranslation_Rate_Y(unsigned short val_iZone) { return  Translation_Rate_Y[val_iZone]; }

inline su2double CConfig::GetTranslation_Rate_Z(unsigned short val_iZone) { return  Translation_Rate_Z[val_iZone]; }

inline su2double CConfig::GetRotation_Rate_X(unsigned short val_iZone) { return  Rotation_Rate_X[val_iZone]; }

inline su2double CConfig::GetRotation_Rate_Y(unsigned short val_iZone) { return  Rotation_Rate_Y[val_iZone]; }

inline su2double CConfig::GetRotation_Rate_Z(unsigned short val_iZone) { return  Rotation_Rate_Z[val_iZone]; }

inline su2double CConfig::GetFinalRotation_Rate_Z(unsigned short val_iZone) { return  FinalRotation_Rate_Z[val_iZone]; }

inline void CConfig::SetRotation_Rate_Z(su2double newRotation_Rate_Z, unsigned short val_iZone) { Rotation_Rate_Z[val_iZone] = newRotation_Rate_Z; }

inline su2double CConfig::GetPitching_Omega_X(unsigned short val_iZone) { return  Pitching_Omega_X[val_iZone]; }

inline su2double CConfig::GetPitching_Omega_Y(unsigned short val_iZone) { return  Pitching_Omega_Y[val_iZone]; }

inline su2double CConfig::GetPitching_Omega_Z(unsigned short val_iZone) { return  Pitching_Omega_Z[val_iZone]; }

inline su2double CConfig::GetPitching_Ampl_X(unsigned short val_iZone) { return  Pitching_Ampl_X[val_iZone]; }

inline su2double CConfig::GetPitching_Ampl_Y(unsigned short val_iZone) { return  Pitching_Ampl_Y[val_iZone]; }

inline su2double CConfig::GetPitching_Ampl_Z(unsigned short val_iZone) { return  Pitching_Ampl_Z[val_iZone]; }

inline su2double CConfig::GetPitching_Phase_X(unsigned short val_iZone) { return  Pitching_Phase_X[val_iZone]; }

inline su2double CConfig::GetPitching_Phase_Y(unsigned short val_iZone) { return  Pitching_Phase_Y[val_iZone]; }

inline su2double CConfig::GetPitching_Phase_Z(unsigned short val_iZone) { return  Pitching_Phase_Z[val_iZone]; }

inline su2double CConfig::GetPlunging_Omega_X(unsigned short val_iZone) { return  Plunging_Omega_X[val_iZone]; }

inline su2double CConfig::GetPlunging_Omega_Y(unsigned short val_iZone) { return  Plunging_Omega_Y[val_iZone]; }

inline su2double CConfig::GetPlunging_Omega_Z(unsigned short val_iZone) { return  Plunging_Omega_Z[val_iZone]; }

inline su2double CConfig::GetPlunging_Ampl_X(unsigned short val_iZone) { return  Plunging_Ampl_X[val_iZone]; }

inline su2double CConfig::GetPlunging_Ampl_Y(unsigned short val_iZone) { return  Plunging_Ampl_Y[val_iZone]; }

inline su2double CConfig::GetPlunging_Ampl_Z(unsigned short val_iZone) { return  Plunging_Ampl_Z[val_iZone]; }

inline su2double* CConfig::GetOmega_HB(void) { return  Omega_HB; }

inline unsigned short CConfig::GetMoveMotion_Origin(unsigned short val_marker) { return MoveMotion_Origin[val_marker]; }

inline su2double CConfig::GetminTurkelBeta() { return  Min_Beta_RoeTurkel; }

inline su2double CConfig::GetmaxTurkelBeta() { return  Max_Beta_RoeTurkel; }

inline unsigned short CConfig::GetKind_Gradient_Method(void) { return Kind_Gradient_Method; }

inline unsigned short CConfig::GetKind_Linear_Solver(void) { return Kind_Linear_Solver; }

inline unsigned short CConfig::GetKind_Deform_Linear_Solver(void) { return Kind_Deform_Linear_Solver; }

inline void CConfig::SetKind_Deform_Linear_Solver_Prec(unsigned short val_kind_prec) { Kind_Deform_Linear_Solver_Prec = val_kind_prec; }

inline unsigned short CConfig::GetKind_Linear_Solver_Prec(void) { return Kind_Linear_Solver_Prec; }

inline void CConfig::SetKind_Linear_Solver_Prec(unsigned short val_kind_prec) { Kind_Linear_Solver_Prec = val_kind_prec; }

inline su2double CConfig::GetLinear_Solver_Error(void) { return Linear_Solver_Error; }

inline su2double CConfig::GetDeform_Linear_Solver_Error(void) { return Deform_Linear_Solver_Error; }

inline unsigned long CConfig::GetLinear_Solver_Iter(void) { return Linear_Solver_Iter; }

inline unsigned long CConfig::GetDeform_Linear_Solver_Iter(void) { return Deform_Linear_Solver_Iter; }

inline unsigned short CConfig::GetLinear_Solver_ILU_n(void) { return Linear_Solver_ILU_n; }

inline unsigned long CConfig::GetLinear_Solver_Restart_Frequency(void) { return Linear_Solver_Restart_Frequency; }

inline su2double CConfig::GetLinear_Solver_Smoother_Relaxation(void) const { return Linear_Solver_Smoother_Relaxation; }

inline su2double CConfig::GetRelaxation_Factor_Flow(void) { return Relaxation_Factor_Flow; }

inline su2double CConfig::GetRelaxation_Factor_AdjFlow(void) { return Relaxation_Factor_AdjFlow; }

inline su2double CConfig::GetRelaxation_Factor_Turb(void) { return Relaxation_Factor_Turb; }

inline su2double CConfig::GetRelaxation_Factor_CHT(void) { return Relaxation_Factor_CHT; }

inline su2double CConfig::GetRoe_Kappa(void) { return Roe_Kappa; }

inline su2double CConfig::GetSemiSpan(void) { return SemiSpan; }

inline unsigned short CConfig::GetKind_AdjTurb_Linear_Solver(void) { return Kind_AdjTurb_Linear_Solver; }

inline unsigned short CConfig::GetKind_AdjTurb_Linear_Prec(void) { return Kind_AdjTurb_Linear_Prec; }

inline unsigned short CConfig::GetKind_DiscAdj_Linear_Solver(void) { return Kind_DiscAdj_Linear_Solver; }

inline unsigned short CConfig::GetKind_DiscAdj_Linear_Prec(void) { return Kind_DiscAdj_Linear_Prec; }

inline unsigned short CConfig::GetKind_Deform_Linear_Solver_Prec(void) { return Kind_Deform_Linear_Solver_Prec; }

inline void CConfig::SetKind_AdjTurb_Linear_Prec(unsigned short val_kind_prec) { Kind_AdjTurb_Linear_Prec = val_kind_prec; }

inline su2double CConfig::GetAdjTurb_Linear_Error(void) { return AdjTurb_Linear_Error; }

inline su2double CConfig::GetEntropyFix_Coeff(void) { return EntropyFix_Coeff; }

inline unsigned short CConfig::GetAdjTurb_Linear_Iter(void) { return AdjTurb_Linear_Iter; }

inline su2double CConfig::GetCFLRedCoeff_AdjTurb(void) { return CFLRedCoeff_AdjTurb; }

inline unsigned long CConfig::GetGridDef_Linear_Iter(void) { return GridDef_Linear_Iter; }

inline unsigned long CConfig::GetGridDef_Nonlinear_Iter(void) { return GridDef_Nonlinear_Iter; }

inline bool CConfig::GetDeform_Output(void) { return Deform_Output; }

inline su2double CConfig::GetDeform_Coeff(void) { return Deform_Coeff; }

inline su2double CConfig::GetDeform_Limit(void) { return Deform_Limit; }

inline su2double CConfig::GetDeform_ElasticityMod(void) { return Deform_ElasticityMod; }

inline su2double CConfig::GetDeform_PoissonRatio(void) { return Deform_PoissonRatio; }

inline unsigned short CConfig::GetDeform_Stiffness_Type(void) { return Deform_Stiffness_Type; }

inline bool CConfig::GetVisualize_Volume_Def(void) { return Visualize_Volume_Def; }

inline bool CConfig::GetVisualize_Surface_Def(void) { return Visualize_Surface_Def; }

inline bool CConfig::GetFFD_Symmetry_Plane(void) { return FFD_Symmetry_Plane; }

inline unsigned short CConfig::GetKind_Adaptation(void) { return Kind_Adaptation; }

inline su2double CConfig::GetNew_Elem_Adapt(void) { return New_Elem_Adapt; }

inline unsigned short CConfig::GetKind_TimeIntScheme(void) { return Kind_TimeNumScheme; }

inline unsigned short CConfig::GetKind_ConvNumScheme(void) { return Kind_ConvNumScheme; }

inline unsigned short CConfig::GetKind_Centered(void) { return Kind_Centered; }

inline unsigned short CConfig::GetKind_Upwind(void) { return Kind_Upwind; }

inline bool CConfig::GetMUSCL(void) { return MUSCL; }

inline bool CConfig::GetMUSCL_Flow(void) { return MUSCL_Flow; }

inline bool CConfig::GetMUSCL_Turb(void) { return MUSCL_Turb; }

inline bool CConfig::GetMUSCL_Heat(void) { return MUSCL_Heat; }

inline bool CConfig::GetMUSCL_AdjFlow(void) { return MUSCL_AdjFlow; }

inline bool CConfig::GetMUSCL_AdjTurb(void) { return MUSCL_AdjTurb; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Flow(void) { return Kind_TimeIntScheme_Flow; }

inline unsigned short CConfig::GetKind_ADER_Predictor(void) { return Kind_ADER_Predictor; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Heat(void) { return Kind_TimeIntScheme_Heat; }

inline unsigned short CConfig::GetKind_TimeStep_Heat(void) { return Kind_TimeStep_Heat; }

inline unsigned short CConfig::GetKind_TimeIntScheme_FEA(void) { return Kind_TimeIntScheme_FEA; }

inline unsigned short CConfig::GetKind_SpaceIteScheme_FEA(void) { return Kind_SpaceIteScheme_FEA; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Flow(void) { return Kind_ConvNumScheme_Flow; }

inline unsigned short CConfig::GetKind_ConvNumScheme_FEM_Flow(void) { return Kind_ConvNumScheme_FEM_Flow; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Template(void) { return Kind_ConvNumScheme_Template; }

inline unsigned short CConfig::GetKind_Centered_Flow(void) { return Kind_Centered_Flow; }

inline unsigned short CConfig::GetKind_SlopeLimit(void) { return Kind_SlopeLimit; }

inline unsigned short CConfig::GetKind_SlopeLimit_Flow(void) { return Kind_SlopeLimit_Flow; }

inline unsigned short CConfig::GetKind_SlopeLimit_Turb(void) { return Kind_SlopeLimit_Turb; }

inline unsigned short CConfig::GetKind_SlopeLimit_AdjTurb(void) { return Kind_SlopeLimit_AdjTurb; }

inline unsigned short CConfig::GetKind_SlopeLimit_AdjFlow(void) { return Kind_SlopeLimit_AdjFlow; }

inline unsigned short CConfig::GetKind_Upwind_Flow(void) { return Kind_Upwind_Flow; }

inline unsigned short CConfig::GetKind_FEM_Flow(void) { return Kind_FEM_Flow; }

inline unsigned short CConfig::GetKind_FEM_DG_Shock(void) { return Kind_FEM_DG_Shock; }

inline unsigned short CConfig::GetKind_Matrix_Coloring(void) { return Kind_Matrix_Coloring; }

inline su2double CConfig::GetKappa_1st_Flow(void) { return Kappa_1st_Flow; }

inline su2double CConfig::GetKappa_2nd_Flow(void) { return Kappa_2nd_Flow; }

inline su2double CConfig::GetKappa_4th_Flow(void) { return Kappa_4th_Flow; }

inline su2double CConfig::GetKappa_2nd_Heat(void) { return Kappa_2nd_Heat; }

inline su2double CConfig::GetKappa_4th_Heat(void) { return Kappa_4th_Heat; }

inline unsigned short CConfig::GetKind_TimeIntScheme_AdjFlow(void) { return Kind_TimeIntScheme_AdjFlow; }

inline unsigned short CConfig::GetKind_ConvNumScheme_AdjFlow(void) { return Kind_ConvNumScheme_AdjFlow; }

inline unsigned short CConfig::GetKind_Centered_AdjFlow(void) { return Kind_Centered_AdjFlow; }

inline unsigned short CConfig::GetKind_Upwind_AdjFlow(void) { return Kind_Upwind_AdjFlow; }

inline su2double CConfig::GetKappa_1st_AdjFlow(void) { return Kappa_1st_AdjFlow; }

inline su2double CConfig::GetKappa_2nd_AdjFlow(void) { return Kappa_2nd_AdjFlow; }

inline su2double CConfig::GetKappa_4th_AdjFlow(void) { return Kappa_4th_AdjFlow; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Turb(void) { return Kind_TimeIntScheme_Turb; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Turb(void) { return Kind_ConvNumScheme_Turb; }

inline unsigned short CConfig::GetKind_Centered_Turb(void) { return Kind_Centered_Turb; }

inline unsigned short CConfig::GetKind_Upwind_Turb(void) {	return Kind_Upwind_Turb; }

inline unsigned short CConfig::GetKind_TimeIntScheme_AdjTurb(void) { return Kind_TimeIntScheme_AdjTurb; }

inline unsigned short CConfig::GetKind_ConvNumScheme_AdjTurb(void) { return Kind_ConvNumScheme_AdjTurb; }

inline unsigned short CConfig::GetKind_Centered_AdjTurb(void) { return Kind_Centered_AdjTurb; }

inline unsigned short CConfig::GetKind_Upwind_AdjTurb(void) { return Kind_Upwind_AdjTurb; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Heat(void) {	return Kind_ConvNumScheme_Heat; }

inline unsigned short CConfig::GetKind_Inlet(void) { return Kind_Inlet; }

inline bool CConfig::GetInlet_Profile_From_File(void) { return Inlet_From_File; }

inline string CConfig::GetInlet_FileName(void) { return Inlet_Filename; }

inline su2double CConfig::GetInlet_Profile_Matching_Tolerance(void) { return Inlet_Matching_Tol; }

inline unsigned short CConfig::GetnInc_Inlet(void) { return nInc_Inlet;}

inline bool CConfig::GetInc_Inlet_UseNormal(void) { return Inc_Inlet_UseNormal;}

inline su2double CConfig::GetInc_Inlet_Damping(void) { return Inc_Inlet_Damping; }

inline su2double CConfig::GetInc_Outlet_Damping(void) { return Inc_Outlet_Damping; }

inline unsigned short CConfig::GetKind_Engine_Inflow(void) { return Kind_Engine_Inflow; }

inline unsigned short CConfig::GetKind_ActDisk(void) { return Kind_ActDisk; }

inline su2double* CConfig::GetFreeStreamTurboNormal(void){return FreeStreamTurboNormal;}

inline unsigned short CConfig::GetKind_AverageProcess(void) { return Kind_AverageProcess; }

inline unsigned short CConfig::GetKind_PerformanceAverageProcess(void) { return Kind_PerformanceAverageProcess; }

inline void CConfig::SetKind_AverageProcess(unsigned short new_AverageProcess) {Kind_AverageProcess = new_AverageProcess; }

inline void CConfig::SetKind_PerformanceAverageProcess(unsigned short new_AverageProcess) {Kind_PerformanceAverageProcess = new_AverageProcess; }

inline su2double CConfig::GetRampRotatingFrame_Coeff(unsigned short iCoeff) { return RampRotatingFrame_Coeff[iCoeff];}

inline bool CConfig::GetRampRotatingFrame(void) { return RampRotatingFrame;}

inline su2double CConfig::GetRampOutletPressure_Coeff(unsigned short iCoeff) { return RampOutletPressure_Coeff[iCoeff];}

inline su2double CConfig::GetFinalOutletPressure(void) { return  FinalOutletPressure; }

inline su2double CConfig::GetMonitorOutletPressure(void) { return MonitorOutletPressure; }

inline void CConfig::SetMonitotOutletPressure(su2double newMonPres) {MonitorOutletPressure = newMonPres;}

inline bool CConfig::GetRampOutletPressure(void) { return RampOutletPressure;}

inline su2double CConfig::GetMixedout_Coeff(unsigned short iCoeff) { return Mixedout_Coeff[iCoeff];}

inline su2double CConfig::GetExtraRelFacGiles(unsigned short iCoeff) { return ExtraRelFacGiles[iCoeff];}

inline su2double CConfig::GetAverageMachLimit(void) { return AverageMachLimit;}

inline unsigned short CConfig::GetKind_MixingPlaneInterface(void) { return Kind_MixingPlaneInterface;}

inline unsigned short CConfig::GetKind_TurboMachinery(unsigned short val_iZone) { return Kind_TurboMachinery[val_iZone]; }

inline unsigned short CConfig::GetKind_SpanWise(void) { return Kind_SpanWise; }

inline bool CConfig::GetBoolMixingPlaneInterface(void) { return (nMarker_MixingPlaneInterface !=0);}

inline bool CConfig::GetBoolTurbomachinery(void) { return (nMarker_Turbomachinery !=0);}

inline bool CConfig::GetBoolZoneSpecific(void) { return ZoneSpecific_Problem;}

inline bool CConfig::GetBoolTurbMixingPlane(void) { return turbMixingPlane;}

inline bool CConfig::GetSpatialFourier(void){return SpatialFourier;}

inline su2double CConfig::GetnBlades(unsigned short val_iZone) { return nBlades[val_iZone];}

inline void CConfig::SetnBlades(unsigned short val_iZone, su2double nblades) { nBlades[val_iZone] = nblades;}

inline bool CConfig::GetBoolGiles(void) { return (nMarker_Giles!=0);}

inline bool CConfig::GetBoolRiemann(void) { return (nMarker_Riemann!=0);}

inline unsigned short CConfig::GetnMarker_MixingPlaneInterface(void) { return nMarker_MixingPlaneInterface;}

inline unsigned short CConfig::GetnMarker_Turbomachinery(void) { return nMarker_Turbomachinery;}

inline unsigned short CConfig::GetnMarker_Shroud(void) { return nMarker_Shroud;}

inline string CConfig::GetMarker_Shroud(unsigned short val_marker){return Marker_Shroud[val_marker];}

inline unsigned short CConfig::GetnMarker_TurboPerformance(void) { return nMarker_TurboPerformance;}

inline unsigned short CConfig::Get_nSpanWiseSections_User(void) { return nSpanWiseSections_User;}

inline unsigned short CConfig::GetnSpanWiseSections(void) { return nSpanWiseSections;}

inline void CConfig::SetnSpanWiseSections(unsigned short nSpan) {nSpanWiseSections = nSpan;}

inline void CConfig::SetnSpanMaxAllZones(unsigned short val_nSpna_max) { nSpanMaxAllZones = val_nSpna_max;}

inline unsigned short CConfig::GetnSpanMaxAllZones(void) { return nSpanMaxAllZones;}

inline void CConfig::SetnSpan_iZones(unsigned short nSpan, unsigned short iZone) {nSpan_iZones[iZone] = nSpan;}

inline unsigned short CConfig::GetnSpan_iZones(unsigned short iZone) { return nSpan_iZones[iZone];}

inline string CConfig::GetMarker_TurboPerf_BoundIn(unsigned short index) { return Marker_TurboBoundIn[index];}

inline string CConfig::GetMarker_TurboPerf_BoundOut(unsigned short index) { return Marker_TurboBoundOut[index];}

inline string CConfig::GetMarker_PerBound(unsigned short val_marker){return Marker_PerBound[val_marker];}

inline unsigned short CConfig::GetnLocationStations(void) { return nLocationStations; }

inline unsigned short CConfig::GetnWingStations(void) { return nWingStations; }

inline su2double CConfig::GetGeo_Waterline_Location(void) { return Geo_Waterline_Location; }

inline void CConfig::SetKind_TimeIntScheme(unsigned short val_kind_timeintscheme) { Kind_TimeNumScheme = val_kind_timeintscheme; }

inline unsigned short CConfig::GetKind_ObjFunc(void) { return Kind_ObjFunc[0]; }

inline unsigned short CConfig::GetKind_ObjFunc(unsigned short val_obj) { return Kind_ObjFunc[val_obj]; }

inline su2double CConfig::GetWeight_ObjFunc(unsigned short val_obj) { return Weight_ObjFunc[val_obj]; }

inline void CConfig::SetWeight_ObjFunc(unsigned short val_obj, su2double val) {Weight_ObjFunc[val_obj] = val; }

inline su2double CConfig::GetCoeff_ObjChainRule(unsigned short iVar) { return Obj_ChainRuleCoeff[iVar]; }

inline unsigned short CConfig::GetKind_SensSmooth(void) { return Kind_SensSmooth; }

inline unsigned short CConfig::GetUnsteady_Simulation(void) { return Unsteady_Simulation; }

inline bool CConfig::GetRestart(void) {	return Restart; }

inline bool CConfig::GetWrt_Binary_Restart(void) {	return Wrt_Binary_Restart; }

inline bool CConfig::GetRead_Binary_Restart(void) {	return Read_Binary_Restart; }

inline bool CConfig::GetRestart_Flow(void) { return Restart_Flow; }

inline bool CConfig::GetEquivArea(void) { return EquivArea; }

inline bool CConfig::GetInvDesign_Cp(void) { return InvDesign_Cp; }

inline bool CConfig::GetInvDesign_HeatFlux(void) { return InvDesign_HeatFlux; }

inline void CConfig::SetnMarker_All(unsigned short val_nmarker) { nMarker_All = val_nmarker; }

inline string CConfig::GetMarker_All_TagBound(unsigned short val_marker) { return Marker_All_TagBound[val_marker]; }

inline string CConfig::GetMarker_ActDiskInlet_TagBound(unsigned short val_marker) { return Marker_ActDiskInlet[val_marker]; }

inline string CConfig::GetMarker_ActDiskOutlet_TagBound(unsigned short val_marker) { return Marker_ActDiskOutlet[val_marker]; }

inline string CConfig::GetMarker_Outlet_TagBound(unsigned short val_marker) { return Marker_Outlet[val_marker]; }

inline string CConfig::GetMarker_EngineInflow_TagBound(unsigned short val_marker) { return Marker_EngineInflow[val_marker]; }

inline string CConfig::GetMarker_EngineExhaust_TagBound(unsigned short val_marker) { return Marker_EngineExhaust[val_marker]; }

inline string CConfig::GetMarker_Monitoring_TagBound(unsigned short val_marker) { return Marker_Monitoring[val_marker]; }

inline string CConfig::GetMarker_HeatFlux_TagBound(unsigned short val_marker) { return Marker_HeatFlux[val_marker]; }

inline string CConfig::GetMarker_Moving_TagBound(unsigned short val_marker) { return Marker_Moving[val_marker]; }

inline string CConfig::GetMarker_PyCustom_TagBound(unsigned short val_marker){ return Marker_PyCustom[val_marker]; }

inline string CConfig::GetMarker_Analyze_TagBound(unsigned short val_marker) { return Marker_Analyze[val_marker]; }

inline short CConfig::GetMarker_All_TagBound(string val_tag) {
	for (unsigned short iMarker = 0; iMarker < nMarker_All; iMarker++) {
		if (val_tag == Marker_All_TagBound[iMarker]) return iMarker;
	}
	return -1;
}

inline unsigned short CConfig::GetMarker_All_KindBC(unsigned short val_marker) { return Marker_All_KindBC[val_marker]; }

inline void CConfig::SetMarker_All_KindBC(unsigned short val_marker, unsigned short val_boundary) { Marker_All_KindBC[val_marker] = val_boundary; }

inline void CConfig::SetMarker_All_TagBound(unsigned short val_marker, string val_index) { Marker_All_TagBound[val_marker] = val_index; }

inline void CConfig::SetMarker_All_Monitoring(unsigned short val_marker, unsigned short val_monitoring) { Marker_All_Monitoring[val_marker] = val_monitoring; }

inline void CConfig::SetMarker_All_GeoEval(unsigned short val_marker, unsigned short val_geoeval) { Marker_All_GeoEval[val_marker] = val_geoeval; }

inline void CConfig::SetMarker_All_Designing(unsigned short val_marker, unsigned short val_designing) { Marker_All_Designing[val_marker] = val_designing; }

inline void CConfig::SetMarker_All_Plotting(unsigned short val_marker, unsigned short val_plotting) { Marker_All_Plotting[val_marker] = val_plotting; }

inline void CConfig::SetMarker_All_Analyze(unsigned short val_marker, unsigned short val_analyze) { Marker_All_Analyze[val_marker] = val_analyze; }

inline void CConfig::SetMarker_All_ZoneInterface(unsigned short val_marker, unsigned short val_fsiinterface) { Marker_All_ZoneInterface[val_marker] = val_fsiinterface; }

inline void CConfig::SetMarker_All_Turbomachinery(unsigned short val_marker, unsigned short val_turbo) { Marker_All_Turbomachinery[val_marker] = val_turbo; }

inline void CConfig::SetMarker_All_TurbomachineryFlag(unsigned short val_marker, unsigned short val_turboflag) { Marker_All_TurbomachineryFlag[val_marker] = val_turboflag; }

inline void CConfig::SetMarker_All_MixingPlaneInterface(unsigned short val_marker, unsigned short val_mixpla_interface) { Marker_All_MixingPlaneInterface[val_marker] = val_mixpla_interface; }

inline void CConfig::SetMarker_All_DV(unsigned short val_marker, unsigned short val_DV) { Marker_All_DV[val_marker] = val_DV; }

inline void CConfig::SetMarker_All_Moving(unsigned short val_marker, unsigned short val_moving) { Marker_All_Moving[val_marker] = val_moving; }

inline void CConfig::SetMarker_All_PyCustom(unsigned short val_marker, unsigned short val_PyCustom) { Marker_All_PyCustom[val_marker] = val_PyCustom; }

inline void CConfig::SetMarker_All_PerBound(unsigned short val_marker, short val_perbound) { Marker_All_PerBound[val_marker] = val_perbound; }

inline short CConfig::GetMarker_All_PerBound(unsigned short val_marker) { return Marker_All_PerBound[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Monitoring(unsigned short val_marker) { return Marker_All_Monitoring[val_marker]; }

inline unsigned short CConfig::GetMarker_All_GeoEval(unsigned short val_marker) { return Marker_All_GeoEval[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Designing(unsigned short val_marker) { return Marker_All_Designing[val_marker]; }

inline short CConfig::GetMarker_All_SendRecv(unsigned short val_marker) { return Marker_All_SendRecv[val_marker]; }

inline void CConfig::SetMarker_All_SendRecv(unsigned short val_marker, short val_index) { Marker_All_SendRecv[val_marker] = val_index; }

inline unsigned short CConfig::GetMarker_All_Plotting(unsigned short val_marker) { return Marker_All_Plotting[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Analyze(unsigned short val_marker) { return Marker_All_Analyze[val_marker]; }

inline unsigned short CConfig::GetMarker_All_ZoneInterface(unsigned short val_marker) { return Marker_All_ZoneInterface[val_marker]; }

inline unsigned short CConfig::GetMarker_n_ZoneInterface(void) { return nMarker_ZoneInterface; }

inline unsigned short CConfig::GetMarker_All_Turbomachinery(unsigned short val_marker) { return Marker_All_Turbomachinery[val_marker]; }

inline unsigned short CConfig::GetMarker_All_TurbomachineryFlag(unsigned short val_marker) { return Marker_All_TurbomachineryFlag[val_marker]; }

inline unsigned short CConfig::GetMarker_All_MixingPlaneInterface(unsigned short val_marker) { return Marker_All_MixingPlaneInterface[val_marker]; }

inline unsigned short CConfig::GetMarker_All_DV(unsigned short val_marker) { return Marker_All_DV[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Moving(unsigned short val_marker) { return Marker_All_Moving[val_marker]; }

inline unsigned short CConfig::GetMarker_All_PyCustom(unsigned short val_marker) { return Marker_All_PyCustom[val_marker];}

inline unsigned short CConfig::GetnMarker_All(void) { return nMarker_All; }

inline unsigned short CConfig::GetnMarker_Max(void) { return nMarker_Max; }

inline unsigned short CConfig::GetnMarker_EngineInflow(void) {	return nMarker_EngineInflow; }

inline unsigned short CConfig::GetnMarker_EngineExhaust(void) { return nMarker_EngineExhaust; }

inline unsigned short CConfig::GetnMarker_InterfaceBound(void) { return nMarker_InterfaceBound; }

inline unsigned short CConfig::GetnMarker_Fluid_InterfaceBound(void) { return nMarker_Fluid_InterfaceBound; }

inline unsigned short CConfig::GetnMarker_Monitoring(void) { return nMarker_Monitoring; }

inline unsigned short CConfig::GetnMarker_DV(void) { return nMarker_DV; }

inline unsigned short CConfig::GetnMarker_Moving(void) { return nMarker_Moving; }

inline unsigned short CConfig::GetnMarker_PyCustom(void) { return nMarker_PyCustom; }

inline unsigned short CConfig::GetnMarker_Analyze(void) { return nMarker_Analyze; }

inline unsigned short CConfig::GetnMarker_NearFieldBound(void) { return nMarker_NearFieldBound; }

inline unsigned short CConfig::GetnMarker_ActDiskInlet(void) { return nMarker_ActDiskInlet; }

inline unsigned short CConfig::GetnMarker_ActDiskOutlet(void) { return nMarker_ActDiskOutlet; }

inline unsigned short CConfig::GetnMarker_Outlet(void) { return nMarker_Outlet; }

inline unsigned short CConfig::GetnMarker_Periodic(void) { return nMarker_PerBound; }

inline unsigned short CConfig::GetnMarker_HeatFlux(void) { return nMarker_HeatFlux; }

inline unsigned short CConfig::GetnObj(void) { return nObj;}

inline string CConfig::GetMesh_FileName(void) { return Mesh_FileName; }

inline string CConfig::GetMesh_Out_FileName(void) { return Mesh_Out_FileName; }

inline unsigned short CConfig::GetMesh_FileFormat(void) { return Mesh_FileFormat; }

inline unsigned short CConfig::GetOutput_FileFormat(void) { return Output_FileFormat; }

inline unsigned short CConfig::GetActDisk_Jump(void) { return ActDisk_Jump; }

inline string CConfig::GetConv_FileName(void) { return Conv_FileName; }

inline string CConfig::GetConv_FileName_FSI(void) { return Conv_FileName_FSI; }

inline string CConfig::GetBreakdown_FileName(void) { return Breakdown_FileName; }

inline string CConfig::GetSolution_FlowFileName(void) { return Solution_FlowFileName; }

inline string CConfig::GetSolution_AdjFileName(void) { return Solution_AdjFileName; }

inline string CConfig::GetSolution_FEMFileName(void) { return Solution_FEMFileName; }

inline string CConfig::GetSolution_AdjFEMFileName(void) { return Solution_AdjFEMFileName; }

inline string CConfig::GetFlow_FileName(void) { return Flow_FileName; }

inline string CConfig::GetStructure_FileName(void) { return Structure_FileName; }

inline string CConfig::GetSurfStructure_FileName(void) { return SurfStructure_FileName; }

inline string CConfig::GetAdjStructure_FileName(void) { return Structure_FileName; }

inline string CConfig::GetAdjSurfStructure_FileName(void) { return SurfStructure_FileName; }

inline string CConfig::GetSurfHeat_FileName(void) { return SurfHeat_FileName; }

inline string CConfig::GetHeat_FileName(void) { return Heat_FileName; }

inline string CConfig::GetRestart_FlowFileName(void) { return Restart_FlowFileName; }

inline string CConfig::GetRestart_HeatFileName(void) { return Restart_HeatFileName; }

inline string CConfig::GetRestart_AdjFileName(void) { return Restart_AdjFileName; }

inline string CConfig::GetRestart_FEMFileName(void) { return Restart_FEMFileName; }

inline string CConfig::GetRestart_AdjFEMFileName(void) { return Restart_AdjFEMFileName; }

inline string CConfig::GetAdj_FileName(void) { return Adj_FileName; }

inline string CConfig::GetObjFunc_Grad_FileName(void) { return ObjFunc_Grad_FileName; }

inline string CConfig::GetObjFunc_Value_FileName(void) { return ObjFunc_Value_FileName; }

inline string CConfig::GetSurfFlowCoeff_FileName(void) { return SurfFlowCoeff_FileName; }

inline string CConfig::GetSurfAdjCoeff_FileName(void) { return SurfAdjCoeff_FileName; }

inline string CConfig::GetSurfSens_FileName(void) { return SurfSens_FileName; }

inline string CConfig::GetVolSens_FileName(void) { return VolSens_FileName; }

inline unsigned short CConfig::GetResidual_Criteria_FEM(void) { return Res_FEM_CRIT; }

inline unsigned short CConfig::GetResidual_Func_Flow(void) { return Residual_Func_Flow; }

inline unsigned short CConfig::GetCauchy_Func_Flow(void) { return Cauchy_Func_Flow; }

inline unsigned short CConfig::GetCauchy_Func_AdjFlow(void) { return Cauchy_Func_AdjFlow; }

inline unsigned short CConfig::GetCauchy_Elems(void) { return Cauchy_Elems; }

inline unsigned long CConfig::GetStartConv_Iter(void) { return StartConv_Iter; }

inline su2double CConfig::GetCauchy_Eps(void) { return Cauchy_Eps; }

inline su2double CConfig::GetDelta_UnstTimeND(void) { return Delta_UnstTimeND; }

inline su2double CConfig::GetTotal_UnstTimeND(void) { return Total_UnstTimeND; }

inline su2double CConfig::GetDelta_UnstTime(void) { return Delta_UnstTime; }

inline su2double CConfig::GetCurrent_UnstTime(void) { return Current_UnstTime; }

inline void CConfig::SetDelta_UnstTimeND(su2double val_delta_unsttimend) { Delta_UnstTimeND = val_delta_unsttimend; }

inline su2double CConfig::GetTotal_UnstTime(void) { return Total_UnstTime; }

inline bool CConfig::GetSubsonicEngine(void) { return SubsonicEngine; }

inline bool CConfig::GetActDisk_DoubleSurface(void) { return ActDisk_DoubleSurface; }

inline bool CConfig::GetEngine_HalfModel(void) { return Engine_HalfModel; }

inline bool CConfig::GetActDisk_SU2_DEF(void) { return ActDisk_SU2_DEF; }

inline su2double CConfig::GetDV_Value(unsigned short val_dv, unsigned short val_value) { return DV_Value[val_dv][val_value]; }

inline void CConfig::SetDV_Value(unsigned short val_dv, unsigned short val_ind, su2double val) { DV_Value[val_dv][val_ind] = val; }

inline su2double CConfig::GetOrderMagResidual(void) { return OrderMagResidual; }

inline su2double CConfig::GetMinLogResidual(void) { return MinLogResidual; }

inline su2double CConfig::GetDamp_Engine_Inflow(void) { return Damp_Engine_Inflow; }

inline su2double CConfig::GetDamp_Engine_Exhaust(void) { return Damp_Engine_Exhaust; }

inline su2double CConfig::GetDamp_Res_Restric(void) { return Damp_Res_Restric; }

inline su2double CConfig::GetDamp_Correc_Prolong(void) { return Damp_Correc_Prolong; }

inline su2double CConfig::GetPosition_Plane(void) { return Position_Plane; }

inline su2double CConfig::GetWeightCd(void) { return WeightCd; }

inline su2double CConfig::GetdCD_dCL(void) { return dCD_dCL; }

inline su2double CConfig::GetdCMx_dCL(void) { return dCMx_dCL; }

inline su2double CConfig::GetdCMy_dCL(void) { return dCMy_dCL; }

inline su2double CConfig::GetdCMz_dCL(void) { return dCMz_dCL; }

inline void CConfig::SetdCD_dCL(su2double val_dcd_dcl) { dCD_dCL = val_dcd_dcl; }

inline void CConfig::SetdCMx_dCL(su2double val_dcmx_dcl) { dCMx_dCL = val_dcmx_dcl; }

inline void CConfig::SetdCMy_dCL(su2double val_dcmy_dcl) { dCMy_dCL = val_dcmy_dcl; }

inline void CConfig::SetdCMz_dCL(su2double val_dcmz_dcl) { dCMz_dCL = val_dcmz_dcl; }

inline void CConfig::SetdCL_dAlpha(su2double val_dcl_dalpha) { dCL_dAlpha = val_dcl_dalpha; }

inline void CConfig::SetdCM_diH(su2double val_dcm_dhi) { dCM_diH = val_dcm_dhi; }

inline su2double CConfig::GetdCD_dCMy(void) { return dCD_dCMy; }

inline void CConfig::SetdCD_dCMy(su2double val_dcd_dcmy) { dCD_dCMy = val_dcd_dcmy; }

inline su2double CConfig::GetCL_Target(void) { return CL_Target; }

inline su2double CConfig::GetCM_Target(void) { return CM_Target; }

inline su2double CConfig::GetFixAzimuthalLine(void) { return FixAzimuthalLine; }

inline su2double CConfig::GetCFLRedCoeff_Turb(void) { return CFLRedCoeff_Turb; }

inline bool CConfig::GetGrid_Movement(void) { return Grid_Movement; }

inline bool CConfig::GetRotating_Frame(void) { return Rotating_Frame; }

inline bool CConfig::GetAxisymmetric(void) { return Axisymmetric; }

inline bool CConfig::GetAdaptBoundary(void) { return AdaptBoundary; }

inline bool CConfig::GetPoissonSolver(void) { return PoissonSolver; }

inline bool CConfig::Low_Mach_Preconditioning(void) { return Low_Mach_Precon; }

inline bool CConfig::Low_Mach_Correction(void) { return Low_Mach_Corr; }

inline bool CConfig::GetGravityForce(void) { return GravityForce; }

inline bool CConfig::GetBody_Force(void) { return Body_Force; }

inline su2double* CConfig::GetBody_Force_Vector(void) { return Body_Force_Vector; }

inline bool CConfig::GetSmoothNumGrid(void) { return SmoothNumGrid; }

inline void CConfig::SetSmoothNumGrid(bool val_smoothnumgrid) { SmoothNumGrid = val_smoothnumgrid; }

inline unsigned short CConfig::GetKind_Turb_Model(void) { return Kind_Turb_Model; }

inline unsigned short CConfig::GetKind_Trans_Model(void) { return Kind_Trans_Model; }

inline unsigned short CConfig::GetKind_SGS_Model(void) { return Kind_SGS_Model; }

inline bool CConfig::GetFrozen_Visc_Cont(void) { return Frozen_Visc_Cont; }

inline bool CConfig::GetFrozen_Visc_Disc(void) { return Frozen_Visc_Disc; }

inline bool CConfig::GetFrozen_Limiter_Disc(void){ return Frozen_Limiter_Disc; }

inline bool CConfig::GetInconsistent_Disc(void){ return Inconsistent_Disc; }

inline bool CConfig::GetSens_Remove_Sharp(void) { return Sens_Remove_Sharp; }

inline bool CConfig::GetWrite_Conv_FSI(void) { return Write_Conv_FSI; }

inline bool CConfig::GetHold_GridFixed(void) { return Hold_GridFixed; }

inline su2double CConfig::GetCyclic_Pitch(void) { return Cyclic_Pitch; }

inline su2double CConfig::GetCollective_Pitch(void) { return Collective_Pitch; }

inline string CConfig::GetDV_Filename(void) { return DV_Filename; }

inline string CConfig::GetDV_Sens_Filename(void) { return DV_Sens_Filename; }

inline string CConfig::GetDV_Unordered_Sens_Filename(void) { return DV_Unordered_Sens_Filename; }

inline bool CConfig::GetLow_MemoryOutput(void) { return Low_MemoryOutput; }

inline bool CConfig::GetWrt_Output(void) { return Wrt_Output; }

inline bool CConfig::GetWrt_Vol_Sol(void) { return Wrt_Vol_Sol; }

inline bool CConfig::GetWrt_Srf_Sol(void) { return Wrt_Srf_Sol; }

inline bool CConfig::GetWrt_Csv_Sol(void) { return Wrt_Csv_Sol; }

inline bool CConfig::GetWrt_Crd_Sol(void) { return Wrt_Crd_Sol; }

inline bool CConfig::GetWrt_Residuals(void) { return Wrt_Residuals; }

inline bool CConfig::GetWrt_Limiters(void) { return Wrt_Limiters; }

inline bool CConfig::GetWrt_Surface(void) { return Wrt_Surface; }

inline bool CConfig::GetWrt_SharpEdges(void) { return Wrt_SharpEdges; }

inline bool CConfig::GetWrt_Halo(void) { return Wrt_Halo; }

inline bool CConfig::GetWrt_Performance(void) { return Wrt_Performance; }

inline bool CConfig::GetWrt_InletFile(void) { return Wrt_InletFile; }

inline void CConfig::SetWrt_InletFile(bool val_wrt_inletfile) { Wrt_InletFile = val_wrt_inletfile; }

inline bool CConfig::GetWrt_Slice(void) { return Wrt_Slice; }

inline bool CConfig::GetWrt_Projected_Sensitivity(void) { return Wrt_Projected_Sensitivity; }

inline unsigned short CConfig::GetSensitivity_Format(void) { return Sensitivity_FileFormat; }

inline bool CConfig::GetPlot_Section_Forces(void) { return Plot_Section_Forces; }

inline vector<vector<su2double> > CConfig::GetAeroelastic_np1(unsigned short iMarker) { return Aeroelastic_np1[iMarker]; }

inline vector<vector<su2double> > CConfig::GetAeroelastic_n(unsigned short iMarker) { return Aeroelastic_n[iMarker]; }

inline vector<vector<su2double> > CConfig::GetAeroelastic_n1(unsigned short iMarker) { return Aeroelastic_n1[iMarker]; }

inline void CConfig::SetAeroelastic_np1(unsigned short iMarker, vector<vector<su2double> > solution) {Aeroelastic_np1[iMarker] = solution;}

inline su2double CConfig::GetAeroelastic_plunge(unsigned short val_marker) { return Aeroelastic_plunge[val_marker]; }

inline su2double CConfig::GetAeroelastic_pitch(unsigned short val_marker) { return Aeroelastic_pitch[val_marker]; }

inline void CConfig::SetAeroelastic_plunge(unsigned short val_marker, su2double val) {Aeroelastic_plunge[val_marker] = val; }

inline void CConfig::SetAeroelastic_pitch(unsigned short val_marker, su2double val) {Aeroelastic_pitch[val_marker] = val; }

inline void CConfig::SetAeroelastic_n1(void) {
        Aeroelastic_n1 = Aeroelastic_n;
}

inline void CConfig::SetAeroelastic_n(void) {
        Aeroelastic_n = Aeroelastic_np1;
}

inline su2double CConfig::GetAeroelastic_Flutter_Speed_Index(void) { return FlutterSpeedIndex; }

inline su2double CConfig::GetAeroelastic_Frequency_Plunge(void) { return PlungeNaturalFrequency; }

inline su2double CConfig::GetAeroelastic_Frequency_Pitch(void) { return PitchNaturalFrequency; }

inline su2double CConfig::GetAeroelastic_Airfoil_Mass_Ratio(void) { return AirfoilMassRatio; }

inline su2double CConfig::GetAeroelastic_CG_Location(void) { return CG_Location; }

inline su2double CConfig::GetAeroelastic_Radius_Gyration_Squared(void) { return RadiusGyrationSquared; }

inline unsigned short CConfig::GetAeroelasticIter(void) { return AeroelasticIter; }

inline bool CConfig::GetWind_Gust(void) { return Wind_Gust; }

inline bool CConfig::GetAeroelastic_Simulation(void) { return Aeroelastic_Simulation; }

inline unsigned short CConfig::GetGust_Type(void) { return Gust_Type; }

inline unsigned short CConfig::GetGust_Dir(void) { return Gust_Dir; }

inline su2double CConfig::GetGust_WaveLength(void) { return Gust_WaveLength; }

inline su2double CConfig::GetGust_Periods(void) { return Gust_Periods; }

inline su2double CConfig::GetGust_Ampl(void) { return Gust_Ampl; }

inline su2double CConfig::GetGust_Begin_Time(void) { return Gust_Begin_Time; }

inline su2double CConfig::GetGust_Begin_Loc(void) { return Gust_Begin_Loc; }

inline unsigned short CConfig::GetnFFD_Iter(void) { return nFFD_Iter; }

inline su2double CConfig::GetFFD_Tol(void) { return FFD_Tol; }

inline su2double CConfig::GetOpt_LineSearch_Bound(void) {return Opt_LineSearch_Bound; }

inline su2double CConfig::GetOpt_RelaxFactor(void) {return Opt_RelaxFactor; }

inline void CConfig::SetOpt_RelaxFactor(su2double val_scale) {Opt_RelaxFactor = val_scale; }

inline long CConfig::GetVisualize_CV(void) { return Visualize_CV; }

inline bool CConfig::GetWall_Functions(void) { return Wall_Functions; }

inline bool CConfig::GetFixed_CL_Mode(void) { return Fixed_CL_Mode; }

inline bool CConfig::GetFixed_CM_Mode(void) { return Fixed_CM_Mode; }

inline bool CConfig::GetEval_dOF_dCX(void) { return Eval_dOF_dCX; }

inline bool CConfig::GetDiscard_InFiles(void) { return Discard_InFiles; }

inline su2double CConfig::GetTarget_CL(void) { return Target_CL; }

inline su2double CConfig::GetdCL_dAlpha(void) { return dCL_dAlpha; }

inline su2double CConfig::GetdCM_diH(void) {return dCM_diH; }

inline unsigned long CConfig::GetIter_Fixed_NetThrust(void) {return Iter_Fixed_NetThrust; }

inline unsigned long CConfig::GetIter_Fixed_CL(void) { return Iter_Fixed_CL; }

inline unsigned long CConfig::GetUpdate_Alpha(void) {return Update_Alpha; }

inline unsigned long CConfig::GetIter_dCL_dAlpha(void) {return Iter_dCL_dAlpha; }

inline bool CConfig::GetUpdate_AoA(void) { return Update_AoA; }

inline bool CConfig::GetUpdate_BCThrust_Bool(void) { return Update_BCThrust_Bool; }

inline void CConfig::SetUpdate_AoA(bool val_update) { Update_AoA = val_update; }

inline unsigned long CConfig::GetUpdate_BCThrust(void) {return Update_BCThrust; }

inline void CConfig::SetUpdate_BCThrust_Bool(bool val_update) { Update_BCThrust_Bool = val_update; }

inline su2double CConfig::GetdNetThrust_dBCThrust(void) {return dNetThrust_dBCThrust; }

inline void CConfig::SetNonphysical_Points(unsigned long val_nonphys_points) { Nonphys_Points = val_nonphys_points; }

inline unsigned long CConfig::GetNonphysical_Points(void) { return Nonphys_Points; }

inline void CConfig::SetNonphysical_Reconstr(unsigned long val_nonphys_reconstr) { Nonphys_Reconstr = val_nonphys_reconstr; }

inline unsigned long CConfig::GetNonphysical_Reconstr(void) { return Nonphys_Reconstr; }

inline unsigned short CConfig::GetConsole_Output_Verb(void) { return Console_Output_Verb; }

inline unsigned short CConfig::GetKind_Average(void) { return Kind_Average; }

inline unsigned short CConfig::GetnIterFSI(void) { return nIterFSI; }

inline unsigned short CConfig::GetnIterFSI_Ramp(void) { return nIterFSI_Ramp; }

inline su2double CConfig::GetAitkenStatRelax(void) { return AitkenStatRelax; }

inline su2double CConfig::GetAitkenDynMaxInit(void) { return AitkenDynMaxInit; }

inline su2double CConfig::GetAitkenDynMinInit(void) { return AitkenDynMinInit; }

inline bool CConfig::GetDeadLoad(void) { return DeadLoad; }

inline bool CConfig::GetPseudoStatic(void) { return PseudoStatic; }

inline bool CConfig::GetSteadyRestart(void) { return SteadyRestart; }

inline unsigned short CConfig::GetDynamic_Analysis(void) { return Dynamic_Analysis; }

inline su2double CConfig::GetDelta_DynTime(void) { return Delta_DynTime; }

inline su2double CConfig::GetTotal_DynTime(void) { return Total_DynTime; }

inline su2double CConfig::GetCurrent_DynTime(void) { return Current_DynTime; }

inline unsigned short CConfig::GetiInst(void) { return iInst; }

inline void CConfig::SetiInst(unsigned short val_iInst) { iInst = val_iInst; }

inline bool CConfig::GetWrt_Dynamic(void) { return Wrt_Dynamic; }

inline su2double CConfig::GetNewmark_beta(void) { return Newmark_beta; }

inline su2double CConfig::GetNewmark_gamma(void) { return Newmark_gamma; }

inline unsigned short CConfig::GetnIntCoeffs(void) { return nIntCoeffs; }

inline su2double CConfig::Get_Int_Coeffs(unsigned short val_coeff) { return Int_Coeffs[val_coeff]; }

inline unsigned short CConfig::GetnElectric_Field(void) { return nElectric_Field; }

inline unsigned short CConfig::GetnDim_Electric_Field(void) { return nDim_Electric_Field; }

inline su2double CConfig::Get_Electric_Field_Mod(unsigned short val_coeff) { return Electric_Field_Mod[val_coeff]; }

inline void CConfig::Set_Electric_Field_Mod(unsigned short val_coeff, su2double val_el_field) { Electric_Field_Mod[val_coeff] = val_el_field; }

inline su2double* CConfig::Get_Electric_Field_Dir(void) { return Electric_Field_Dir; }

inline bool CConfig::GetRamp_Load(void) { return Ramp_Load; }

inline su2double CConfig::GetRamp_Time(void) { return Ramp_Time; }

inline bool CConfig::GetRampAndRelease_Load(void) { return RampAndRelease; }

inline bool CConfig::GetSine_Load(void) { return Sine_Load; }

inline su2double* CConfig::GetLoad_Sine(void) { return SineLoad_Coeff; }

inline su2double CConfig::GetRefGeom_Penalty(void) { return RefGeom_Penalty; }

inline su2double CConfig::GetTotalDV_Penalty(void) { return DV_Penalty; }

inline bool CConfig::GetPredictor(void) { return Predictor; }

inline unsigned short CConfig::GetPredictorOrder(void) { return Pred_Order; }

inline bool CConfig::GetRelaxation(void) { return Relaxation; }

inline bool CConfig::GetIncrementalLoad(void) { return IncrementalLoad; }

inline unsigned long CConfig::GetNumberIncrements(void) { return IncLoad_Nincrements; }

inline su2double CConfig::GetIncLoad_Criteria(unsigned short val_var) { return IncLoad_Criteria[val_var]; }

inline bool CConfig::GetEulerPersson(void) { return EulerPersson; }

inline void CConfig::SetEulerPersson(bool val_EulerPersson) { EulerPersson = val_EulerPersson; }

inline bool CConfig::GetFSI_Simulation(void) { return FSI_Problem; }

inline void CConfig::SetFSI_Simulation(bool FSI_sim) { FSI_Problem = FSI_sim; }

inline void CConfig::SetMultizone_Problem(bool MZ_problem) { Multizone_Problem = MZ_problem; }

inline bool CConfig::GetMultizone_Problem(void) { return Multizone_Problem; }

inline unsigned short CConfig::GetnID_DV(void) { return nID_DV; }

inline unsigned short CConfig::GetKindInterpolation(void) { return Kind_Interpolation; }

inline unsigned short CConfig::GetKindRadialBasisFunction(void) { return Kind_RadialBasisFunction; }

inline bool CConfig::GetRadialBasisFunctionPolynomialOption(void) {return RadialBasisFunction_PolynomialOption; }

inline su2double CConfig::GetRadialBasisFunctionParameter(void) {return RadialBasisFunction_Parameter; }

inline bool CConfig::GetConservativeInterpolation(void) { return ConservativeInterpolation; }

inline unsigned short CConfig::GetRelaxation_Method_FSI(void) { return Kind_BGS_RelaxMethod; }

inline su2double CConfig::GetOrderMagResidualFSI(void) { return OrderMagResidualFSI; }

inline su2double CConfig::GetMinLogResidualFSI(void) { return MinLogResidualFSI; }

inline su2double CConfig::GetOrderMagResidual_BGS_F(void) { return OrderMagResidual_BGS_F; }

inline su2double CConfig::GetMinLogResidual_BGS_F(void) { return MinLogResidual_BGS_F; }

inline su2double CConfig::GetOrderMagResidual_BGS_S(void) { return OrderMagResidual_BGS_S; }

inline su2double CConfig::GetMinLogResidual_BGS_S(void) { return MinLogResidual_BGS_S; }

inline su2double CConfig::GetResidual_FEM_UTOL(void) { return Res_FEM_UTOL; }

inline su2double CConfig::GetResidual_FEM_RTOL(void) { return Res_FEM_RTOL; }

inline su2double CConfig::GetResidual_FEM_ETOL(void) { return Res_FEM_ETOL; }

inline su2double CConfig::GetCriteria_FEM_ADJ(void) { return Res_FEM_ADJ; }

inline unsigned short CConfig::GetDynamic_LoadTransfer(void) { return Dynamic_LoadTransfer; }

inline unsigned short CConfig::GetDirectDiff() { return DirectDiff;}

inline bool CConfig::GetDiscrete_Adjoint() { return DiscreteAdjoint;}

inline unsigned short CConfig::GetRiemann_Solver_FEM(void) {return Riemann_Solver_FEM;}

inline su2double CConfig::GetQuadrature_Factor_Straight(void) {return Quadrature_Factor_Straight;}

inline su2double CConfig::GetQuadrature_Factor_Curved(void) {return Quadrature_Factor_Curved;}

inline su2double CConfig::GetQuadrature_Factor_Time_ADER_DG(void) {return Quadrature_Factor_Time_ADER_DG;}

inline su2double CConfig::GetTheta_Interior_Penalty_DGFEM(void) {return Theta_Interior_Penalty_DGFEM;}

inline unsigned short CConfig::GetSizeMatMulPadding(void) {return sizeMatMulPadding;}

inline bool CConfig::GetCompute_Entropy(void) {return Compute_Entropy;}

inline bool CConfig::GetUse_Lumped_MassMatrix_DGFEM(void) {return Use_Lumped_MassMatrix_DGFEM;}

inline bool CConfig::GetJacobian_Spatial_Discretization_Only(void) {return Jacobian_Spatial_Discretization_Only;}

inline bool CConfig::GetWeakly_Coupled_Heat(void) { return Weakly_Coupled_Heat; }

inline bool CConfig::GetIntegrated_HeatFlux(void) { return Integrated_HeatFlux; }

inline bool CConfig::GetAD_Mode(void) { return AD_Mode;}

inline bool CConfig::GetAD_Preaccumulation(void) {return AD_Preaccumulation;}

inline unsigned short CConfig::GetFFD_Blending(void){return FFD_Blending;}

inline su2double* CConfig::GetFFD_BSplineOrder(){return FFD_BSpline_Order;}

inline void CConfig::SetMax_Vel2(su2double val_max_vel2) { Max_Vel2 = val_max_vel2; }

inline su2double CConfig::GetMax_Vel2(void) { return Max_Vel2; }

inline void CConfig::SetRestart_Bandwidth_Agg(su2double val_restart_bandwidth_sum) { Restart_Bandwidth_Agg = val_restart_bandwidth_sum; }

inline su2double CConfig::GetRestart_Bandwidth_Agg(void) { return Restart_Bandwidth_Agg; }

inline unsigned long CConfig::GetWrt_Surf_Freq_DualTime(void) { return Wrt_Surf_Freq_DualTime; }

inline unsigned short CConfig::GetKind_HybridRANSLES(void) {return Kind_HybridRANSLES; }

inline unsigned short CConfig::GetKind_RoeLowDiss(void) {return Kind_RoeLowDiss; }

inline su2double CConfig::GetConst_DES(void) {return Const_DES; }

inline bool CConfig::GetQCR(void) {return QCR;}

inline bool CConfig::GetCompute_Average(void) {return Compute_Average;}

inline unsigned short CConfig::GetVerification_Solution(void) {return Kind_Verification_Solution;}

inline ofstream* CConfig::GetHistFile(void) { return ConvHistFile; }

inline void CConfig::SetHistFile(ofstream *HistFile) { ConvHistFile = HistFile; }

inline unsigned short CConfig::GetComm_Level(void) { return Comm_Level; }

inline bool CConfig::GetTopology_Optimization(void) const { return topology_optimization; }

inline string CConfig::GetTopology_Optim_FileName(void) const { return top_optim_output_file; }

inline su2double CConfig::GetSIMP_Exponent(void) const { return simp_exponent; }

inline su2double CConfig::GetSIMP_MinStiffness(void) const { return simp_minimum_stiffness; }
  
inline unsigned short CConfig::GetTopology_Optim_Num_Kernels(void) const { return top_optim_nKernel; }
  
inline void CConfig::GetTopology_Optim_Kernel(const unsigned short iKernel, unsigned short &type,
                                              su2double &param, su2double &radius) const {
  type = top_optim_kernels[iKernel];
  param = top_optim_kernel_params[iKernel];
  radius = top_optim_filter_radius[iKernel];
}

inline void CConfig::GetTopology_Optim_Projection(unsigned short &type, su2double &param) const {
  type = top_optim_proj_type;  param = top_optim_proj_param;
}

inline string CConfig::GetConfigFilename(unsigned short index) { return Config_Filenames[index]; }

inline unsigned short CConfig::GetnConfigFiles(void) { return nConfig_Files; }

inline unsigned short CConfig::GetnMarker_ZoneInterface(void) { return nMarker_ZoneInterface; }

inline string CConfig::GetMarkerTag_ZoneInterface(unsigned short val_iMarker) { return Marker_ZoneInterface[val_iMarker]; }

inline bool CConfig::GetTime_Domain(void) { return Time_Domain; }

inline unsigned long CConfig::GetnInner_Iter(void) { return Inner_Iter; }

inline unsigned long CConfig::GetnOuter_Iter(void) { return Outer_Iter; }

inline unsigned long CConfig::GetnTime_Iter(void) { return Time_Iter; }

inline unsigned long CConfig::GetnIter(void) { return Iter; }

inline unsigned long CConfig::GetRestart_Iter(void) { return Restart_Iter; }

inline su2double CConfig::GetTime_Step(void) { return Time_Step; }

inline su2double CConfig::GetMax_Time(void) { return Max_Time; }

inline bool CConfig::GetMultizone_Mesh(void) { return Multizone_Mesh; }

inline bool CConfig::GetMultizone_Residual(void) { return Multizone_Residual; }

inline bool CConfig::GetSinglezone_Driver(void) { return SinglezoneDriver; }

inline bool CConfig::GetSpecial_Output(void) { return SpecialOutput; }

inline bool CConfig::GetWrt_ForcesBreakdown(void) { return Wrt_ForcesBreakdown; }

inline bool CConfig::GetUsing_UQ(void) { return using_uq; }

inline su2double CConfig::GetUQ_Delta_B(void) { return uq_delta_b; }

inline unsigned short CConfig::GetEig_Val_Comp(void) {return eig_val_comp; }

inline su2double CConfig::GetUQ_URLX(void) {return uq_urlx; }

inline bool CConfig::GetUQ_Permute(void) { return uq_permute; }
