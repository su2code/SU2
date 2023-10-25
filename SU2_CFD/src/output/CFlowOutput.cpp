/*!
 * \file CFlowOutput.cpp
 * \brief Common functions for flow output.
 * \author R. Sanchez
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include <sstream>
#include <string>
#include <sstream>
#include <iomanip>

#include "../../include/output/CFlowOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/solvers/CSolver.hpp"
#include "../../include/variables/CPrimitiveIndices.hpp"
#include "../../include/fluid/CCoolProp.hpp"

CFlowOutput::CFlowOutput(const CConfig *config, unsigned short nDim, bool fem_output) :
  CFVMOutput(config, nDim, fem_output),
  lastInnerIter(curInnerIter) {
}

// The "AddHistoryOutput(" must not be split over multiple lines to ensure proper python parsing
// clang-format off
void CFlowOutput::AddAnalyzeSurfaceOutput(const CConfig *config){

  /// DESCRIPTION: Average mass flow
  AddHistoryOutput("SURFACE_MASSFLOW",         "Avg_Massflow",              ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average mass flow on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Mach number
  AddHistoryOutput("SURFACE_MACH",             "Avg_Mach",                  ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average mach number on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Temperature
  AddHistoryOutput("SURFACE_STATIC_TEMPERATURE","Avg_Temp",                 ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average temperature on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Pressure
  AddHistoryOutput("SURFACE_STATIC_PRESSURE",  "Avg_Press",                 ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average pressure on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Density
  AddHistoryOutput("AVG_DENSITY",              "Avg_Density",               ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average density on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Enthalpy
  AddHistoryOutput("AVG_ENTHALPY",             "Avg_Enthalpy",              ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average enthalpy on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average velocity in normal direction of the surface
  AddHistoryOutput("AVG_NORMALVEL",            "Avg_NormalVel",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average normal velocity on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Flow uniformity
  AddHistoryOutput("SURFACE_UNIFORMITY",       "Uniformity",                ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total flow uniformity on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Secondary strength
  AddHistoryOutput("SURFACE_SECONDARY",        "Secondary_Strength",        ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total secondary strength on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Momentum distortion
  AddHistoryOutput("SURFACE_MOM_DISTORTION",   "Momentum_Distortion",       ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total momentum distortion on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Secondary over uniformity
  AddHistoryOutput("SURFACE_SECOND_OVER_UNIFORM","Secondary_Over_Uniformity",ScreenOutputFormat::SCIENTIFIC,"FLOW_COEFF", "Total secondary over uniformity on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average total temperature
  AddHistoryOutput("SURFACE_TOTAL_TEMPERATURE","Avg_TotalTemp",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average total temperature all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average total pressure
  AddHistoryOutput("SURFACE_TOTAL_PRESSURE",   "Avg_TotalPress",            ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average total pressure on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Pressure drop
  if (config->GetnMarker_Analyze() >= 2) {
    AddHistoryOutput("SURFACE_PRESSURE_DROP",    "Pressure_Drop",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total pressure drop on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  } else if (rank == MASTER_NODE) {
    cout << "\nWARNING: SURFACE_PRESSURE_DROP can only be computed for at least 2 surfaces (outlet, inlet, ...)\n" << endl;
  }
  if (config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT) {
    /// DESCRIPTION: Average Species
    for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
      AddHistoryOutput("SURFACE_SPECIES_" + std::to_string(iVar), "Avg_Species_" + std::to_string(iVar), ScreenOutputFormat::FIXED, "SPECIES_COEFF", "Total average species " + std::to_string(iVar) + " on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
    }
    /// DESCRIPTION: Species Variance
    AddHistoryOutput("SURFACE_SPECIES_VARIANCE", "Species_Variance", ScreenOutputFormat::SCIENTIFIC, "SPECIES_COEFF", "Total species variance, measure for mixing quality. On all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  }
  /// END_GROUP

  /// BEGIN_GROUP: AERO_COEFF_SURF, DESCRIPTION: Surface values on non-solid markers.
  vector<string> Marker_Analyze;
  for (unsigned short iMarker_Analyze = 0; iMarker_Analyze < config->GetnMarker_Analyze(); iMarker_Analyze++){
    Marker_Analyze.push_back(config->GetMarker_Analyze_TagBound(iMarker_Analyze));
  }

  /// DESCRIPTION: Average mass flow
  AddHistoryOutputPerSurface("SURFACE_MASSFLOW",         "Avg_Massflow",              ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Mach number
  AddHistoryOutputPerSurface("SURFACE_MACH",             "Avg_Mach",                  ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Temperature
  AddHistoryOutputPerSurface("SURFACE_STATIC_TEMPERATURE","Avg_Temp",                 ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Pressure
  AddHistoryOutputPerSurface("SURFACE_STATIC_PRESSURE",  "Avg_Press",                 ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Density
  AddHistoryOutputPerSurface("AVG_DENSITY",              "Avg_Density",               ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Enthalpy
  AddHistoryOutputPerSurface("AVG_ENTHALPY",             "Avg_Enthalpy",              ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average velocity in normal direction of the surface
  AddHistoryOutputPerSurface("AVG_NORMALVEL",            "Avg_NormalVel",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Flow uniformity
  AddHistoryOutputPerSurface("SURFACE_UNIFORMITY",       "Uniformity",                ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Secondary strength
  AddHistoryOutputPerSurface("SURFACE_SECONDARY",        "Secondary_Strength",        ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Momentum distortion
  AddHistoryOutputPerSurface("SURFACE_MOM_DISTORTION",   "Momentum_Distortion",       ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Secondary over uniformity
  AddHistoryOutputPerSurface("SURFACE_SECOND_OVER_UNIFORM","Secondary_Over_Uniformity",ScreenOutputFormat::SCIENTIFIC,"FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average total temperature
  AddHistoryOutputPerSurface("SURFACE_TOTAL_TEMPERATURE","Avg_TotalTemp",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average total pressure
  AddHistoryOutputPerSurface("SURFACE_TOTAL_PRESSURE",   "Avg_TotalPress",            ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  if (config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT) {
    /// DESCRIPTION: Average Species
    for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
      AddHistoryOutputPerSurface("SURFACE_SPECIES_" + std::to_string(iVar), "Avg_Species_" + std::to_string(iVar), ScreenOutputFormat::FIXED, "SPECIES_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
    }
    /// DESCRIPTION: Species Variance
    AddHistoryOutputPerSurface("SURFACE_SPECIES_VARIANCE", "Species_Variance", ScreenOutputFormat::SCIENTIFIC, "SPECIES_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  }
}
// clang-format on

void CFlowOutput::SetAnalyzeSurface(const CSolver* const*solver, const CGeometry *geometry, CConfig *config, bool output){

  unsigned short iDim, iMarker, iMarker_Analyze;
  unsigned long iVertex, iPoint;
  su2double Mach = 0.0, Pressure, Temperature = 0.0, TotalPressure = 0.0, TotalTemperature = 0.0,
  Enthalpy, Velocity[3] = {0.0}, TangVel[3], Vector[3], Velocity2, MassFlow, Density, Area,
  SoundSpeed, Vn, Vn2, Vtang2, Weight = 1.0;

  const su2double Gas_Constant      = config->GetGas_ConstantND();
  const su2double Gamma             = config->GetGamma();
  const unsigned short nMarker      = config->GetnMarker_All();
  const unsigned short nDim         = geometry->GetnDim();
  const unsigned short Kind_Average = config->GetKind_Average();

  const bool compressible   = config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE;
  const bool incompressible = config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE;
  const bool energy         = config->GetEnergy_Equation();
  const bool streamwisePeriodic = (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE);
  const bool species        = config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT;
  const auto nSpecies       = config->GetnSpecies();

  const bool axisymmetric               = config->GetAxisymmetric();
  const unsigned short nMarker_Analyze  = config->GetnMarker_Analyze();

  const auto flow_nodes = solver[FLOW_SOL]->GetNodes();
  const CVariable* species_nodes = species ? solver[SPECIES_SOL]->GetNodes() : nullptr;

  vector<su2double> Surface_MassFlow          (nMarker,0.0);
  vector<su2double> Surface_Mach              (nMarker,0.0);
  vector<su2double> Surface_Temperature       (nMarker,0.0);
  vector<su2double> Surface_Density           (nMarker,0.0);
  vector<su2double> Surface_Enthalpy          (nMarker,0.0);
  vector<su2double> Surface_NormalVelocity    (nMarker,0.0);
  vector<su2double> Surface_StreamVelocity2   (nMarker,0.0);
  vector<su2double> Surface_TransvVelocity2   (nMarker,0.0);
  vector<su2double> Surface_Pressure          (nMarker,0.0);
  vector<su2double> Surface_TotalTemperature  (nMarker,0.0);
  vector<su2double> Surface_TotalPressure     (nMarker,0.0);
  vector<su2double> Surface_VelocityIdeal     (nMarker,0.0);
  vector<su2double> Surface_Area              (nMarker,0.0);
  vector<su2double> Surface_MassFlow_Abs      (nMarker,0.0);
  su2activematrix Surface_Species(nMarker, nSpecies);
  Surface_Species = su2double(0.0);

  su2double  Tot_Surface_MassFlow          = 0.0;
  su2double  Tot_Surface_Mach              = 0.0;
  su2double  Tot_Surface_Temperature       = 0.0;
  su2double  Tot_Surface_Density           = 0.0;
  su2double  Tot_Surface_Enthalpy          = 0.0;
  su2double  Tot_Surface_NormalVelocity    = 0.0;
  su2double  Tot_Surface_StreamVelocity2   = 0.0;
  su2double  Tot_Surface_TransvVelocity2   = 0.0;
  su2double  Tot_Surface_Pressure          = 0.0;
  su2double  Tot_Surface_TotalTemperature  = 0.0;
  su2double  Tot_Surface_TotalPressure     = 0.0;
  su2double  Tot_Momentum_Distortion       = 0.0;
  su2double  Tot_SecondOverUniformity      = 0.0;
  vector<su2double> Tot_Surface_Species(nSpecies,0.0);

  /*--- Compute the numerical fan face Mach number, and the total area of the inflow ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if (config->GetMarker_All_Analyze(iMarker) == YES) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->nodes->GetDomain(iPoint)) {

          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);

          const su2double AxiFactor = GetAxiFactor(axisymmetric, *geometry->nodes, iPoint, iMarker);

          Density = flow_nodes->GetDensity(iPoint);
          Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0; Vn = 0.0; Vtang2 = 0.0;

          for (iDim = 0; iDim < nDim; iDim++) {
            Area += (Vector[iDim] * AxiFactor) * (Vector[iDim] * AxiFactor);
            Velocity[iDim] = flow_nodes->GetVelocity(iPoint,iDim);
            Velocity2 += Velocity[iDim] * Velocity[iDim];
            Vn += Velocity[iDim] * Vector[iDim] * AxiFactor;
            MassFlow += Vector[iDim] * AxiFactor * Density * Velocity[iDim];
          }

          Area       = sqrt (Area);
          if (AxiFactor == 0.0) Vn = 0.0; else Vn /= Area;
          Vn2        = Vn * Vn;
          Pressure   = flow_nodes->GetPressure(iPoint);
          /*--- Use recovered pressure here as pressure difference between in and outlet is zero otherwise  ---*/
          if(streamwisePeriodic) Pressure = flow_nodes->GetStreamwise_Periodic_RecoveredPressure(iPoint);
          SoundSpeed = flow_nodes->GetSoundSpeed(iPoint);

          for (iDim = 0; iDim < nDim; iDim++) {
            TangVel[iDim] = Velocity[iDim] - Vn*Vector[iDim]*AxiFactor/Area;
            Vtang2       += TangVel[iDim]*TangVel[iDim];
          }

          if (incompressible){
            if (config->GetVariable_Density_Model()) {
              Mach = sqrt(flow_nodes->GetVelocity2(iPoint))/
              sqrt(flow_nodes->GetSpecificHeatCp(iPoint)*config->GetPressure_ThermodynamicND()/(flow_nodes->GetSpecificHeatCv(iPoint)*flow_nodes->GetDensity(iPoint)));
            } else {
              Mach = sqrt(flow_nodes->GetVelocity2(iPoint))/
              sqrt(config->GetBulk_Modulus()/(flow_nodes->GetDensity(iPoint)));
            }
            Temperature       = flow_nodes->GetTemperature(iPoint);
            Enthalpy          = flow_nodes->GetSpecificHeatCp(iPoint)*Temperature;
            TotalTemperature  = Temperature + 0.5*Velocity2/flow_nodes->GetSpecificHeatCp(iPoint);
            TotalPressure     = Pressure + 0.5*Density*Velocity2;
          }
          else{
            Mach              = sqrt(Velocity2)/SoundSpeed;
            Temperature       = Pressure / (Gas_Constant * Density);
            Enthalpy          = flow_nodes->GetEnthalpy(iPoint);
            TotalTemperature  = Temperature * (1.0 + Mach * Mach * 0.5 * (Gamma - 1.0));
            TotalPressure     = Pressure * pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), Gamma / (Gamma - 1.0));
          }

          /*--- Compute the mass Surface_MassFlow ---*/

          Surface_Area[iMarker]             += Area;
          Surface_MassFlow[iMarker]         += MassFlow;
          Surface_MassFlow_Abs[iMarker]     += abs(MassFlow);

          if (Kind_Average == AVERAGE_MASSFLUX) Weight = abs(MassFlow);
          else if (Kind_Average == AVERAGE_AREA) Weight = abs(Area);
          else Weight = 1.0;

          Surface_Mach[iMarker]             += Mach*Weight;
          Surface_Temperature[iMarker]      += Temperature*Weight;
          Surface_Density[iMarker]          += Density*Weight;
          Surface_Enthalpy[iMarker]         += Enthalpy*Weight;
          Surface_NormalVelocity[iMarker]   += Vn*Weight;
          Surface_Pressure[iMarker]         += Pressure*Weight;
          Surface_TotalTemperature[iMarker] += TotalTemperature*Weight;
          Surface_TotalPressure[iMarker]    += TotalPressure*Weight;
          if (species)
            for (unsigned short iVar = 0; iVar < nSpecies; iVar++)
              Surface_Species(iMarker, iVar) += species_nodes->GetSolution(iPoint, iVar)*Weight;

          /*--- For now, always used the area to weight the uniformities. ---*/

          Weight = abs(Area);

          Surface_StreamVelocity2[iMarker]   += Vn2*Weight;
          Surface_TransvVelocity2[iMarker]   += Vtang2*Weight;

        }
      }
    }
  }

  /*--- Copy to the appropriate structure ---*/

  vector<su2double> Surface_MassFlow_Local          (nMarker_Analyze,0.0);
  vector<su2double> Surface_Mach_Local              (nMarker_Analyze,0.0);
  vector<su2double> Surface_Temperature_Local       (nMarker_Analyze,0.0);
  vector<su2double> Surface_Density_Local           (nMarker_Analyze,0.0);
  vector<su2double> Surface_Enthalpy_Local          (nMarker_Analyze,0.0);
  vector<su2double> Surface_NormalVelocity_Local    (nMarker_Analyze,0.0);
  vector<su2double> Surface_StreamVelocity2_Local   (nMarker_Analyze,0.0);
  vector<su2double> Surface_TransvVelocity2_Local   (nMarker_Analyze,0.0);
  vector<su2double> Surface_Pressure_Local          (nMarker_Analyze,0.0);
  vector<su2double> Surface_TotalTemperature_Local  (nMarker_Analyze,0.0);
  vector<su2double> Surface_TotalPressure_Local     (nMarker_Analyze,0.0);
  vector<su2double> Surface_Area_Local              (nMarker_Analyze,0.0);
  vector<su2double> Surface_MassFlow_Abs_Local      (nMarker_Analyze,0.0);
  su2activematrix Surface_Species_Local(nMarker_Analyze,nSpecies);
  Surface_Species_Local = su2double(0.0);

  vector<su2double> Surface_MassFlow_Total          (nMarker_Analyze,0.0);
  vector<su2double> Surface_Mach_Total              (nMarker_Analyze,0.0);
  vector<su2double> Surface_Temperature_Total       (nMarker_Analyze,0.0);
  vector<su2double> Surface_Density_Total           (nMarker_Analyze,0.0);
  vector<su2double> Surface_Enthalpy_Total          (nMarker_Analyze,0.0);
  vector<su2double> Surface_NormalVelocity_Total    (nMarker_Analyze,0.0);
  vector<su2double> Surface_StreamVelocity2_Total   (nMarker_Analyze,0.0);
  vector<su2double> Surface_TransvVelocity2_Total   (nMarker_Analyze,0.0);
  vector<su2double> Surface_Pressure_Total          (nMarker_Analyze,0.0);
  vector<su2double> Surface_TotalTemperature_Total  (nMarker_Analyze,0.0);
  vector<su2double> Surface_TotalPressure_Total     (nMarker_Analyze,0.0);
  vector<su2double> Surface_Area_Total              (nMarker_Analyze,0.0);
  vector<su2double> Surface_MassFlow_Abs_Total      (nMarker_Analyze,0.0);
  su2activematrix Surface_Species_Total(nMarker_Analyze,nSpecies);
  Surface_Species_Total = su2double(0.0);

  vector<su2double> Surface_MomentumDistortion_Total (nMarker_Analyze,0.0);

  /*--- Compute the numerical fan face Mach number, mach number, temperature and the total area ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if (config->GetMarker_All_Analyze(iMarker) == YES)  {

      for (iMarker_Analyze= 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {

        /*--- Add the Surface_MassFlow, and Surface_Area to the particular boundary ---*/

        if (config->GetMarker_All_TagBound(iMarker) == config->GetMarker_Analyze_TagBound(iMarker_Analyze)) {
          Surface_MassFlow_Local[iMarker_Analyze]          += Surface_MassFlow[iMarker];
          Surface_Mach_Local[iMarker_Analyze]              += Surface_Mach[iMarker];
          Surface_Temperature_Local[iMarker_Analyze]       += Surface_Temperature[iMarker];
          Surface_Density_Local[iMarker_Analyze]           += Surface_Density[iMarker];
          Surface_Enthalpy_Local[iMarker_Analyze]          += Surface_Enthalpy[iMarker];
          Surface_NormalVelocity_Local[iMarker_Analyze]    += Surface_NormalVelocity[iMarker];
          Surface_StreamVelocity2_Local[iMarker_Analyze]   += Surface_StreamVelocity2[iMarker];
          Surface_TransvVelocity2_Local[iMarker_Analyze]   += Surface_TransvVelocity2[iMarker];
          Surface_Pressure_Local[iMarker_Analyze]          += Surface_Pressure[iMarker];
          Surface_TotalTemperature_Local[iMarker_Analyze]  += Surface_TotalTemperature[iMarker];
          Surface_TotalPressure_Local[iMarker_Analyze]     += Surface_TotalPressure[iMarker];
          Surface_Area_Local[iMarker_Analyze]              += Surface_Area[iMarker];
          Surface_MassFlow_Abs_Local[iMarker_Analyze]      += Surface_MassFlow_Abs[iMarker];
          for (unsigned short iVar = 0; iVar < nSpecies; iVar++)
            Surface_Species_Local(iMarker_Analyze, iVar) += Surface_Species(iMarker, iVar);
        }

      }

    }

  }

  auto Allreduce = [](const vector<su2double>& src, vector<su2double>& dst) {
    SU2_MPI::Allreduce(src.data(), dst.data(), src.size(), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  };

  auto Allreduce_su2activematrix = [](const su2activematrix& src, su2activematrix& dst) {
    SU2_MPI::Allreduce(src.data(), dst.data(), src.size(), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  };

  Allreduce(Surface_MassFlow_Local, Surface_MassFlow_Total);
  Allreduce(Surface_Mach_Local, Surface_Mach_Total);
  Allreduce(Surface_Temperature_Local, Surface_Temperature_Total);
  Allreduce(Surface_Density_Local, Surface_Density_Total);
  Allreduce(Surface_Enthalpy_Local, Surface_Enthalpy_Total);
  Allreduce(Surface_NormalVelocity_Local, Surface_NormalVelocity_Total);
  Allreduce(Surface_StreamVelocity2_Local, Surface_StreamVelocity2_Total);
  Allreduce(Surface_TransvVelocity2_Local, Surface_TransvVelocity2_Total);
  Allreduce(Surface_Pressure_Local, Surface_Pressure_Total);
  Allreduce(Surface_TotalTemperature_Local, Surface_TotalTemperature_Total);
  Allreduce(Surface_TotalPressure_Local, Surface_TotalPressure_Total);
  Allreduce(Surface_Area_Local, Surface_Area_Total);
  Allreduce(Surface_MassFlow_Abs_Local, Surface_MassFlow_Abs_Total);
  Allreduce_su2activematrix(Surface_Species_Local, Surface_Species_Total);

  /*--- Compute the value of Surface_Area_Total, and Surface_Pressure_Total, and
   set the value in the config structure for future use ---*/

  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {

    if (Kind_Average == AVERAGE_MASSFLUX) Weight = Surface_MassFlow_Abs_Total[iMarker_Analyze];
    else if (Kind_Average == AVERAGE_AREA) Weight = abs(Surface_Area_Total[iMarker_Analyze]);
    else Weight = 1.0;

    if (Weight != 0.0) {
      Surface_Mach_Total[iMarker_Analyze]             /= Weight;
      Surface_Temperature_Total[iMarker_Analyze]      /= Weight;
      Surface_Density_Total[iMarker_Analyze]          /= Weight;
      Surface_Enthalpy_Total[iMarker_Analyze]         /= Weight;
      Surface_NormalVelocity_Total[iMarker_Analyze]   /= Weight;
      Surface_Pressure_Total[iMarker_Analyze]         /= Weight;
      Surface_TotalTemperature_Total[iMarker_Analyze] /= Weight;
      Surface_TotalPressure_Total[iMarker_Analyze]    /= Weight;
      for (unsigned short iVar = 0; iVar < nSpecies; iVar++)
        Surface_Species_Total(iMarker_Analyze, iVar) /= Weight;
    }
    else {
      Surface_Mach_Total[iMarker_Analyze]             = 0.0;
      Surface_Temperature_Total[iMarker_Analyze]      = 0.0;
      Surface_Density_Total[iMarker_Analyze]          = 0.0;
      Surface_Enthalpy_Total[iMarker_Analyze]         = 0.0;
      Surface_NormalVelocity_Total[iMarker_Analyze]   = 0.0;
      Surface_Pressure_Total[iMarker_Analyze]         = 0.0;
      Surface_TotalTemperature_Total[iMarker_Analyze] = 0.0;
      Surface_TotalPressure_Total[iMarker_Analyze]    = 0.0;
      for (unsigned short iVar = 0; iVar < nSpecies; iVar++)
        Surface_Species_Total(iMarker_Analyze, iVar) = 0.0;
    }

    /*--- Compute flow uniformity parameters separately (always area for now). ---*/

    Area = fabs(Surface_Area_Total[iMarker_Analyze]);

    /*--- The definitions for Distortion and Uniformity Parameters are taken as defined by Banko, Andrew J., et al. in section 3.2 of
    https://www.sciencedirect.com/science/article/pii/S0142727X16301412 ------*/

    if (Area != 0.0) {
      Surface_MomentumDistortion_Total[iMarker_Analyze] = Surface_StreamVelocity2_Total[iMarker_Analyze]/(Surface_NormalVelocity_Total[iMarker_Analyze]*Surface_NormalVelocity_Total[iMarker_Analyze]*Area) - 1.0;
      Surface_StreamVelocity2_Total[iMarker_Analyze] /= Area;
      Surface_TransvVelocity2_Total[iMarker_Analyze] /= Area;
    }
    else {
      Surface_MomentumDistortion_Total[iMarker_Analyze] = 0.0;
      Surface_StreamVelocity2_Total[iMarker_Analyze]    = 0.0;
      Surface_TransvVelocity2_Total[iMarker_Analyze]    = 0.0;
    }

  }

  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {

    su2double MassFlow = Surface_MassFlow_Total[iMarker_Analyze] * config->GetDensity_Ref() * config->GetVelocity_Ref();
    if (us_units) MassFlow *= 32.174;
    SetHistoryOutputPerSurfaceValue("SURFACE_MASSFLOW", MassFlow, iMarker_Analyze);
    Tot_Surface_MassFlow += MassFlow;
    config->SetSurface_MassFlow(iMarker_Analyze, MassFlow);

    su2double Mach = Surface_Mach_Total[iMarker_Analyze];
    SetHistoryOutputPerSurfaceValue("SURFACE_MACH", Mach, iMarker_Analyze);
    Tot_Surface_Mach += Mach;
    config->SetSurface_Mach(iMarker_Analyze, Mach);

    su2double Temperature = Surface_Temperature_Total[iMarker_Analyze] * config->GetTemperature_Ref();
    SetHistoryOutputPerSurfaceValue("SURFACE_STATIC_TEMPERATURE", Temperature, iMarker_Analyze);
    Tot_Surface_Temperature += Temperature;
    config->SetSurface_Temperature(iMarker_Analyze, Temperature);

    su2double Pressure = Surface_Pressure_Total[iMarker_Analyze] * config->GetPressure_Ref();
    SetHistoryOutputPerSurfaceValue("SURFACE_STATIC_PRESSURE", Pressure, iMarker_Analyze);
    Tot_Surface_Pressure += Pressure;
    config->SetSurface_Pressure(iMarker_Analyze, Pressure);

    su2double Density = Surface_Density_Total[iMarker_Analyze] * config->GetDensity_Ref();
    SetHistoryOutputPerSurfaceValue("AVG_DENSITY", Density, iMarker_Analyze);
    Tot_Surface_Density += Density;
    config->SetSurface_Density(iMarker_Analyze, Density);

    su2double Enthalpy = Surface_Enthalpy_Total[iMarker_Analyze];
    SetHistoryOutputPerSurfaceValue("AVG_ENTHALPY", Enthalpy, iMarker_Analyze);
    Tot_Surface_Enthalpy += Enthalpy;
    config->SetSurface_Enthalpy(iMarker_Analyze, Enthalpy);

    su2double NormalVelocity = Surface_NormalVelocity_Total[iMarker_Analyze] * config->GetVelocity_Ref();
    SetHistoryOutputPerSurfaceValue("AVG_NORMALVEL", NormalVelocity, iMarker_Analyze);
    Tot_Surface_NormalVelocity += NormalVelocity;
    config->SetSurface_NormalVelocity(iMarker_Analyze, NormalVelocity);

    su2double Uniformity = sqrt(Surface_StreamVelocity2_Total[iMarker_Analyze]) * config->GetVelocity_Ref();
    SetHistoryOutputPerSurfaceValue("SURFACE_UNIFORMITY", Uniformity, iMarker_Analyze);
    Tot_Surface_StreamVelocity2 += Uniformity;
    config->SetSurface_Uniformity(iMarker_Analyze, Uniformity);

    su2double SecondaryStrength = sqrt(Surface_TransvVelocity2_Total[iMarker_Analyze]) * config->GetVelocity_Ref();
    SetHistoryOutputPerSurfaceValue("SURFACE_SECONDARY", SecondaryStrength, iMarker_Analyze);
    Tot_Surface_TransvVelocity2 += SecondaryStrength;
    config->SetSurface_SecondaryStrength(iMarker_Analyze, SecondaryStrength);

    su2double MomentumDistortion = Surface_MomentumDistortion_Total[iMarker_Analyze];
    SetHistoryOutputPerSurfaceValue("SURFACE_MOM_DISTORTION", MomentumDistortion, iMarker_Analyze);
    Tot_Momentum_Distortion += MomentumDistortion;
    config->SetSurface_MomentumDistortion(iMarker_Analyze, MomentumDistortion);

    su2double SecondOverUniform = SecondaryStrength/Uniformity;
    SetHistoryOutputPerSurfaceValue("SURFACE_SECOND_OVER_UNIFORM", SecondOverUniform, iMarker_Analyze);
    Tot_SecondOverUniformity += SecondOverUniform;
    config->SetSurface_SecondOverUniform(iMarker_Analyze, SecondOverUniform);

    su2double TotalTemperature = Surface_TotalTemperature_Total[iMarker_Analyze] * config->GetTemperature_Ref();
    SetHistoryOutputPerSurfaceValue("SURFACE_TOTAL_TEMPERATURE", TotalTemperature, iMarker_Analyze);
    Tot_Surface_TotalTemperature += TotalTemperature;
    config->SetSurface_TotalTemperature(iMarker_Analyze, TotalTemperature);

    su2double TotalPressure = Surface_TotalPressure_Total[iMarker_Analyze] * config->GetPressure_Ref();
    SetHistoryOutputPerSurfaceValue("SURFACE_TOTAL_PRESSURE", TotalPressure, iMarker_Analyze);
    Tot_Surface_TotalPressure += TotalPressure;
    config->SetSurface_TotalPressure(iMarker_Analyze, TotalPressure);

    if (species) {
      for (unsigned short iVar = 0; iVar < nSpecies; iVar++) {
        su2double Species = Surface_Species_Total(iMarker_Analyze, iVar);
        SetHistoryOutputPerSurfaceValue("SURFACE_SPECIES_" + std::to_string(iVar), Species, iMarker_Analyze);
        Tot_Surface_Species[iVar] += Species;
        if (iVar == 0)
          config->SetSurface_Species_0(iMarker_Analyze, Species);
      }
    }
  }

  /*--- Compute the average static pressure drop between two surfaces. Note
   that this assumes we have two surfaces being analyzed and that the outlet
   is first followed by the inlet. This is because we may also want to choose
   outlet values (temperature, uniformity, etc.) for our design problems,
   which require the outlet to be listed first. This is a simple first version
   that could be generalized to a different orders/lists/etc. ---*/

  if (nMarker_Analyze >= 2) {
    su2double PressureDrop = (Surface_Pressure_Total[1] - Surface_Pressure_Total[0]) * config->GetPressure_Ref();
    for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
      config->SetSurface_PressureDrop(iMarker_Analyze, PressureDrop);
    }
    SetHistoryOutputValue("SURFACE_PRESSURE_DROP", PressureDrop);
  }
  SetHistoryOutputValue("SURFACE_MASSFLOW", Tot_Surface_MassFlow);
  SetHistoryOutputValue("SURFACE_MACH", Tot_Surface_Mach);
  SetHistoryOutputValue("SURFACE_STATIC_TEMPERATURE", Tot_Surface_Temperature);
  SetHistoryOutputValue("SURFACE_STATIC_PRESSURE", Tot_Surface_Pressure);
  SetHistoryOutputValue("AVG_DENSITY", Tot_Surface_Density);
  SetHistoryOutputValue("AVG_ENTHALPY", Tot_Surface_Enthalpy);
  SetHistoryOutputValue("AVG_NORMALVEL", Tot_Surface_NormalVelocity);
  SetHistoryOutputValue("SURFACE_UNIFORMITY", Tot_Surface_StreamVelocity2);
  SetHistoryOutputValue("SURFACE_SECONDARY", Tot_Surface_TransvVelocity2);
  SetHistoryOutputValue("SURFACE_MOM_DISTORTION", Tot_Momentum_Distortion);
  SetHistoryOutputValue("SURFACE_SECOND_OVER_UNIFORM", Tot_SecondOverUniformity);
  SetHistoryOutputValue("SURFACE_TOTAL_TEMPERATURE", Tot_Surface_TotalTemperature);
  SetHistoryOutputValue("SURFACE_TOTAL_PRESSURE", Tot_Surface_TotalPressure);
  if (species) {
    for (unsigned short iVar = 0; iVar < nSpecies; iVar++)
      SetHistoryOutputValue("SURFACE_SPECIES_" + std::to_string(iVar), Tot_Surface_Species[iVar]);

    SetAnalyzeSurfaceSpeciesVariance(solver, geometry, config, Surface_Species_Total, Surface_MassFlow_Abs_Total,
                                      Surface_Area_Total);
  }

  if ((rank == MASTER_NODE) && !config->GetDiscrete_Adjoint() && output) {

    cout.precision(6);
    cout.setf(ios::scientific, ios::floatfield);
    cout << endl << "Computing surface mean values." << endl << endl;

    for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
      cout << "Surface "<< config->GetMarker_Analyze_TagBound(iMarker_Analyze) << ":" << endl;

      if (nDim == 3) { if (si_units) cout << setw(20) << "Area (m^2): "; else cout << setw(20) << "Area (ft^2): "; }
      else { if (si_units) cout << setw(20) << "Area (m): "; else cout << setw(20) << "Area (ft): "; }

      if (si_units)      cout << setw(15) << fabs(Surface_Area_Total[iMarker_Analyze]);
      else if (us_units) cout << setw(15) << fabs(Surface_Area_Total[iMarker_Analyze])*12.0*12.0;

      cout << endl;

      su2double MassFlow = config->GetSurface_MassFlow(iMarker_Analyze);
      if (si_units)      cout << setw(20) << "Mf (kg/s): " << setw(15) << MassFlow;
      else if (us_units) cout << setw(20) << "Mf (lbs/s): " << setw(15) << MassFlow;

      su2double NormalVelocity = config->GetSurface_NormalVelocity(iMarker_Analyze);
      if (si_units)      cout << setw(20) << "Vn (m/s): " << setw(15) << NormalVelocity;
      else if (us_units) cout << setw(20) << "Vn (ft/s): " << setw(15) << NormalVelocity;

      cout << endl;

      su2double Uniformity = config->GetSurface_Uniformity(iMarker_Analyze);
      if (si_units)      cout << setw(20) << "Uniformity (m/s): " << setw(15) << Uniformity;
      else if (us_units) cout << setw(20) << "Uniformity (ft/s): " << setw(15) << Uniformity;

      su2double SecondaryStrength = config->GetSurface_SecondaryStrength(iMarker_Analyze);
      if (si_units)      cout << setw(20) << "Secondary (m/s): " << setw(15) << SecondaryStrength;
      else if (us_units) cout << setw(20) << "Secondary (ft/s): " << setw(15) << SecondaryStrength;

      cout << endl;

      su2double MomentumDistortion = config->GetSurface_MomentumDistortion(iMarker_Analyze);
      cout << setw(20) << "Mom. Distortion: " << setw(15) << MomentumDistortion;

      su2double SecondOverUniform = config->GetSurface_SecondOverUniform(iMarker_Analyze);
      cout << setw(20) << "Second/Uniform: " << setw(15) << SecondOverUniform;

      cout << endl;

      su2double Pressure = config->GetSurface_Pressure(iMarker_Analyze);
      if (si_units)      cout << setw(20) << "P (Pa): " << setw(15) << Pressure;
      else if (us_units) cout << setw(20) << "P (psf): " << setw(15) << Pressure;

      su2double TotalPressure = config->GetSurface_TotalPressure(iMarker_Analyze);
      if (si_units)      cout << setw(20) << "PT (Pa): " << setw(15) <<TotalPressure;
      else if (us_units) cout << setw(20) << "PT (psf): " << setw(15) <<TotalPressure;

      cout << endl;

      su2double Mach = config->GetSurface_Mach(iMarker_Analyze);
      cout << setw(20) << "Mach: " << setw(15) << Mach;

      su2double Density = config->GetSurface_Density(iMarker_Analyze);
      if (si_units)      cout << setw(20) << "Rho (kg/m^3): " << setw(15) << Density;
      else if (us_units) cout << setw(20) << "Rho (lb/ft^3): " << setw(15) << Density*32.174;

      cout << endl;

      if (compressible || energy) {
        su2double Temperature = config->GetSurface_Temperature(iMarker_Analyze);
        if (si_units)      cout << setw(20) << "T (K): " << setw(15) << Temperature;
        else if (us_units) cout << setw(20) << "T (R): " << setw(15) << Temperature;

        su2double TotalTemperature = config->GetSurface_TotalTemperature(iMarker_Analyze);
        if (si_units)      cout << setw(20) << "TT (K): " << setw(15) << TotalTemperature;
        else if (us_units) cout << setw(20) << "TT (R): " << setw(15) << TotalTemperature;

        cout << endl;
      }

    }
    cout.unsetf(ios_base::floatfield);

  }

  std::cout << std::resetiosflags(std::cout.flags());
}

void CFlowOutput::SetAnalyzeSurfaceSpeciesVariance(const CSolver* const*solver, const CGeometry *geometry,
                                                    CConfig *config, const su2activematrix& Surface_Species_Total,
                                                    const vector<su2double>& Surface_MassFlow_Abs_Total,
                                                    const vector<su2double>& Surface_Area_Total) {

  const unsigned short nMarker      = config->GetnMarker_All();
  const unsigned short Kind_Average = config->GetKind_Average();

  const bool species        = config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT;
  const auto nSpecies       = config->GetnSpecies();

  const bool axisymmetric               = config->GetAxisymmetric();
  const unsigned short nMarker_Analyze  = config->GetnMarker_Analyze();

  const auto flow_nodes = solver[FLOW_SOL]->GetNodes();
  const CVariable* species_nodes = species ? solver[SPECIES_SOL]->GetNodes() : nullptr;

  /*--- Compute Variance of species on the analyze markers. This is done after the rest as the average species value is
   * necessary. The variance is computed for all species together and not for each species alone. ---*/
  vector<su2double> Surface_SpeciesVariance(nMarker,0.0);
  su2double Tot_Surface_SpeciesVariance = 0.0;

  /*--- sum += (Yj_i - mu_Yj)^2 * weight_i with i representing the node and j the species. ---*/
  for (unsigned short iMarker = 0; iMarker < nMarker; iMarker++) {

    if (config->GetMarker_All_Analyze(iMarker) == YES) {

      /*--- Find iMarkerAnalyze to iMarker. As SpeciesAvg is accessed via iMarkerAnalyze. ---*/
      unsigned short iMarker_Analyze_Stored = std::numeric_limits<unsigned short>::max();
      for (unsigned short iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++)
        if (config->GetMarker_All_TagBound(iMarker) == config->GetMarker_Analyze_TagBound(iMarker_Analyze))
          iMarker_Analyze_Stored = iMarker_Analyze;

      for (unsigned long iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->nodes->GetDomain(iPoint)) {

          const su2double AxiFactor = GetAxiFactor(axisymmetric, *geometry->nodes, iPoint, iMarker);

          su2double Vector[3];
          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
          const su2double Density = flow_nodes->GetDensity(iPoint);
          su2double Area = 0.0;
          su2double MassFlow = 0.0;

          for (unsigned short iDim = 0; iDim < nDim; iDim++) {
            Area += (Vector[iDim] * AxiFactor) * (Vector[iDim] * AxiFactor);
            MassFlow += Vector[iDim] * AxiFactor * Density * flow_nodes->GetVelocity(iPoint,iDim);
          }
          Area= sqrt(Area);

          su2double Weight;
          if (Kind_Average == AVERAGE_MASSFLUX) Weight = abs(MassFlow);
          else if (Kind_Average == AVERAGE_AREA) Weight = abs(Area);
          else Weight = 1.0;

          for (unsigned short iVar = 0; iVar < nSpecies; iVar++)
            Surface_SpeciesVariance[iMarker] += pow(species_nodes->GetSolution(iPoint, iVar) - Surface_Species_Total(iMarker_Analyze_Stored, iVar), 2) * Weight;
        }
      }
    }
  }

  /*--- MPI Communication ---*/
  vector<su2double> Surface_SpeciesVariance_Local(nMarker_Analyze,0.0);
  vector<su2double> Surface_SpeciesVariance_Total(nMarker_Analyze,0.0);

  for (unsigned short iMarker = 0; iMarker < nMarker; iMarker++) {

    if (config->GetMarker_All_Analyze(iMarker) == YES)  {

      for (unsigned short iMarker_Analyze= 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {

        /*--- Add the Surface_MassFlow, and Surface_Area to the particular boundary ---*/

        if (config->GetMarker_All_TagBound(iMarker) == config->GetMarker_Analyze_TagBound(iMarker_Analyze)) {
          Surface_SpeciesVariance_Local[iMarker_Analyze] += Surface_SpeciesVariance[iMarker];
        }
      }
    }
  }

  auto Allreduce = [](const vector<su2double>& src, vector<su2double>& dst) {
    SU2_MPI::Allreduce(src.data(), dst.data(), src.size(), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  };
  Allreduce(Surface_SpeciesVariance_Local, Surface_SpeciesVariance_Total);

  /*--- Divide quantity by weight. ---*/
  for (unsigned short iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {

    su2double Weight;
    if (Kind_Average == AVERAGE_MASSFLUX) Weight = Surface_MassFlow_Abs_Total[iMarker_Analyze];
    else if (Kind_Average == AVERAGE_AREA) Weight = abs(Surface_Area_Total[iMarker_Analyze]);
    else Weight = 1.0;

    if (Weight != 0.0) {
      Surface_SpeciesVariance_Total[iMarker_Analyze] /= Weight;
    }
    else {
      Surface_SpeciesVariance[iMarker_Analyze] = 0.0;
    }
  }

  /*--- Set values on markers ---*/
  for (unsigned short iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
    su2double SpeciesVariance = Surface_SpeciesVariance_Total[iMarker_Analyze];
    SetHistoryOutputPerSurfaceValue("SURFACE_SPECIES_VARIANCE", SpeciesVariance, iMarker_Analyze);
    config->SetSurface_Species_Variance(iMarker_Analyze, SpeciesVariance);
    Tot_Surface_SpeciesVariance += SpeciesVariance;
  }
  SetHistoryOutputValue("SURFACE_SPECIES_VARIANCE", Tot_Surface_SpeciesVariance);
}

void CFlowOutput::ConvertVariableSymbolsToIndices(const CPrimitiveIndices<unsigned long>& idx, const bool allowSkip,
                                                  CustomOutput& output) const {
  const auto nameToIndex = PrimitiveNameToIndexMap(idx);

  std::stringstream knownVariables;
  for (const auto& items : nameToIndex) {
    knownVariables << items.first + '\n';
  }
  knownVariables << "TURB[0,1,...]\nRAD[0,1,...]\nSPECIES[0,1,...]\nSCALAR[0,1,...]\n";

  auto IndexOfVariable = [](const map<std::string, unsigned long>& nameToIndex, const std::string& var) {
    /*--- Primitives of the flow solver. ---*/
    const auto flowOffset = FLOW_SOL * CustomOutput::MAX_VARS_PER_SOLVER;
    const auto it = nameToIndex.find(var);
    if (it != nameToIndex.end()) return flowOffset + it->second;

    /*--- Index-based (no name) access to variables of other solvers. ---*/
    auto GetIndex = [](const std::string& s, int nameLen) {
      /*--- Extract an int from "name[int]", nameLen is the length of "name". ---*/
      return std::stoi(std::string(s.begin() + nameLen + 1, s.end() - 1));
    };

    if (var.rfind("SPECIES", 0) == 0) return SPECIES_SOL * CustomOutput::MAX_VARS_PER_SOLVER + GetIndex(var, 7);
    if (var.rfind("SCALAR", 0) == 0) return SPECIES_SOL * CustomOutput::MAX_VARS_PER_SOLVER + GetIndex(var, 6);
    if (var.rfind("TURB", 0) == 0) return TURB_SOL * CustomOutput::MAX_VARS_PER_SOLVER + GetIndex(var, 4);
    if (var.rfind("RAD", 0) == 0) return RAD_SOL * CustomOutput::MAX_VARS_PER_SOLVER + GetIndex(var, 3);
    return CustomOutput::NOT_A_VARIABLE;
  };

  output.otherOutputs.clear();
  output.varIndices.clear();
  output.varIndices.reserve(output.varSymbols.size());

  for (const auto& var : output.varSymbols) {
    output.varIndices.push_back(IndexOfVariable(nameToIndex, var));

    if (output.type == OperationType::FUNCTION && output.varIndices.back() != CustomOutput::NOT_A_VARIABLE) {
      SU2_MPI::Error("Custom outputs of type 'Function' cannot reference solver variables.", CURRENT_FUNCTION);
    }
    /*--- Symbol is a valid solver variable. ---*/
    if (output.varIndices.back() < CustomOutput::NOT_A_VARIABLE) continue;

    /*--- An index above NOT_A_VARIABLE is not valid with current solver settings. ---*/
    if (output.varIndices.back() > CustomOutput::NOT_A_VARIABLE) {
      SU2_MPI::Error("Inactive solver variable (" + var + ") used in function " + output.name + "\n"
                      "E.g. this may only be a variable of the compressible solver.", CURRENT_FUNCTION);
    }

    /*--- An index equal to NOT_A_VARIABLE may refer to a history output. ---*/
    output.varIndices.back() += output.otherOutputs.size();
    output.otherOutputs.push_back(GetPtrToHistoryOutput(var));
    if (output.otherOutputs.back() == nullptr) {
      if (!allowSkip) {
        SU2_MPI::Error("Invalid history output or solver variable (" + var + ") used in function " + output.name +
                       "\nValid solvers variables:\n" + knownVariables.str(), CURRENT_FUNCTION);
      } else {
        if (rank == MASTER_NODE) {
          std::cout << "Info: Ignoring function " + output.name + " because it may be used by the primal/adjoint "
                       "solver.\n      If the function is ignored twice it is invalid." << std::endl;
        }
        output.skip = true;
        break;
      }
    }
  }
}

void CFlowOutput::SetCustomOutputs(const CSolver* const* solver, const CGeometry *geometry, const CConfig *config) {

  const bool adjoint = config->GetDiscrete_Adjoint();
  const bool axisymmetric = config->GetAxisymmetric();
  const auto* flowNodes = su2staticcast_p<const CFlowVariable*>(solver[FLOW_SOL]->GetNodes());

  for (auto& output : customOutputs) {
    if (output.skip) continue;

    if (output.varIndices.empty()) {
      const bool allowSkip = adjoint && (output.type == OperationType::FUNCTION);

      /*--- Setup indices for the symbols in the expression. ---*/
      const auto primIdx = CPrimitiveIndices<unsigned long>(config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE,
          config->GetNEMOProblem(), nDim, config->GetnSpecies());
      ConvertVariableSymbolsToIndices(primIdx, allowSkip, output);
      if (output.skip) continue;

      /*--- Convert marker names to their index (if any) in this rank. Or probe locations to nearest points. ---*/

      if (output.type != OperationType::PROBE) {
        output.markerIndices.clear();
        for (const auto& marker : output.markers) {
          for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); ++iMarker) {
            if (config->GetMarker_All_TagBound(iMarker) == marker) {
              output.markerIndices.push_back(iMarker);
              continue;
            }
          }
        }
      } else {
        if (output.markers.size() != nDim) {
          SU2_MPI::Error("Wrong number of coordinates to specify probe " + output.name, CURRENT_FUNCTION);
        }
        su2double coord[3] = {};
        for (auto iDim = 0u; iDim < nDim; ++iDim) coord[iDim] = std::stod(output.markers[iDim]);
        su2double minDist = std::numeric_limits<su2double>::max();
        unsigned long minPoint = 0;
        for (auto iPoint = 0ul; iPoint < geometry->GetnPointDomain(); ++iPoint) {
          const su2double dist = GeometryToolbox::SquaredDistance(nDim, coord, geometry->nodes->GetCoord(iPoint));
          if (dist < minDist) {
            minDist = dist;
            minPoint = iPoint;
          }
        }
        /*--- Decide which rank owns the probe. ---*/
        su2double globMinDist;
        SU2_MPI::Allreduce(&minDist, &globMinDist, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
        output.iPoint = fabs(minDist - globMinDist) < EPS ? minPoint : CustomOutput::PROBE_NOT_OWNED;
        if (output.iPoint != CustomOutput::PROBE_NOT_OWNED) {
          std::cout << "Probe " << output.name << " is using global point "
                    << geometry->nodes->GetGlobalIndex(output.iPoint)
                    << ", distance from target location is " << sqrt(minDist) << std::endl;
        }
      }
    }

    if (output.type == OperationType::FUNCTION) {
      auto Functor = [&](unsigned long i) {
        /*--- Functions only reference other history outputs. ---*/
        return *output.otherOutputs[i - CustomOutput::NOT_A_VARIABLE];
      };
      SetHistoryOutputValue(output.name, output.Eval(Functor));
      continue;
    }

    /*--- Prepares the functor that maps symbol indices to values at a given point
     * (see ConvertVariableSymbolsToIndices). ---*/

    auto MakeFunctor = [&](unsigned long iPoint) {
      /*--- This returns another lambda that captures iPoint by value. ---*/
      return [&, iPoint](unsigned long i) {
        if (i < CustomOutput::NOT_A_VARIABLE) {
          const auto solIdx = i / CustomOutput::MAX_VARS_PER_SOLVER;
          const auto varIdx = i % CustomOutput::MAX_VARS_PER_SOLVER;
          if (solIdx == FLOW_SOL) {
            return flowNodes->GetPrimitive(iPoint, varIdx);
          }
          return solver[solIdx]->GetNodes()->GetSolution(iPoint, varIdx);
        } else {
          return *output.otherOutputs[i - CustomOutput::NOT_A_VARIABLE];
        }
      };
    };

    if (output.type == OperationType::PROBE) {
      su2double value = std::numeric_limits<su2double>::max();
      if (output.iPoint != CustomOutput::PROBE_NOT_OWNED) {
        value = output.Eval(MakeFunctor(output.iPoint));
      }
      su2double tmp = value;
      SU2_MPI::Allreduce(&tmp, &value, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      SetHistoryOutputValue(output.name, value);
      continue;
    }

    /*--- Surface integral of the expression. ---*/

    std::array<su2double, 2> integral = {0.0, 0.0};

    SU2_OMP_PARALLEL {
      std::array<su2double, 2> local_integral = {0.0, 0.0};

      for (const auto iMarker : output.markerIndices) {

        SU2_OMP_FOR_(schedule(static) SU2_NOWAIT)
        for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; ++iVertex) {
          const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (!geometry->nodes->GetDomain(iPoint)) continue;

          const auto* normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          su2double weight = 1.0;
          if (output.type == OperationType::MASSFLOW_AVG || output.type == OperationType::MASSFLOW_INT) {
            weight = flowNodes->GetDensity(iPoint) * flowNodes->GetProjVel(iPoint, normal);
          } else {
            weight = GeometryToolbox::Norm(nDim, normal);
          }
          weight *= GetAxiFactor(axisymmetric, *geometry->nodes, iPoint, iMarker);
          local_integral[1] += weight;
          local_integral[0] += weight * output.Eval(MakeFunctor(iPoint));
        }
        END_SU2_OMP_FOR
      }

      SU2_OMP_CRITICAL {
        integral[0] += local_integral[0];
        integral[1] += local_integral[1];
      }
      END_SU2_OMP_CRITICAL
    }
    END_SU2_OMP_PARALLEL

    const auto local = integral;
    SU2_MPI::Allreduce(local.data(), integral.data(), 2, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    if (output.type == OperationType::AREA_AVG || output.type == OperationType::MASSFLOW_AVG) {
      integral[0] /= integral[1];
    }
    SetHistoryOutputValue(output.name, integral[0]);
  }
}

// The "AddHistoryOutput(" must not be split over multiple lines to ensure proper python parsing
// clang-format off
void CFlowOutput::AddHistoryOutputFields_ScalarRMS_RES(const CConfig* config) {

  switch (TurbModelFamily(config->GetKind_Turb_Model())) {
    case TURB_FAMILY::SA:
      /// DESCRIPTION: Root-mean square residual of nu tilde (SA model).
      AddHistoryOutput("RMS_NU_TILDE", "rms[nu]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of nu tilde (SA model).", HistoryFieldType::RESIDUAL);
      break;

    case TURB_FAMILY::KW:
      /// DESCRIPTION: Root-mean square residual of kinetic energy (SST model).
      AddHistoryOutput("RMS_TKE", "rms[k]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of kinetic energy (SST model).", HistoryFieldType::RESIDUAL);
      /// DESCRIPTION: Root-mean square residual of the dissipation (SST model).
      AddHistoryOutput("RMS_DISSIPATION", "rms[w]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of dissipation (SST model).", HistoryFieldType::RESIDUAL);
      break;

    case TURB_FAMILY::NONE: break;
  }
  switch (config->GetKind_Trans_Model()) {

    case TURB_TRANS_MODEL::LM:
      /// DESCRIPTION: Root-mean square residual of the intermittency (LM model).
      AddHistoryOutput("RMS_INTERMITTENCY", "rms[LM_1]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of intermittency (LM model).", HistoryFieldType::RESIDUAL);
      /// DESCRIPTION: Root-mean square residual of the momentum thickness Reynolds number (LM model).
      AddHistoryOutput("RMS_RE_THETA_T", "rms[LM_2]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of momentum thickness Reynolds number (LM model).", HistoryFieldType::RESIDUAL);
      break;

    case TURB_TRANS_MODEL::NONE: break;
  }

  switch (config->GetKind_Species_Model()) {
    case SPECIES_MODEL::SPECIES_TRANSPORT: {
      for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
        AddHistoryOutput("RMS_SPECIES_" + std::to_string(iVar), "rms[rho*Y_" + std::to_string(iVar)+"]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of transported species.", HistoryFieldType::RESIDUAL);
      }
      break;
    }
    case SPECIES_MODEL::FLAMELET: {
      /*--- Controlling variable transport. ---*/
      for (auto iCV = 0u; iCV < config->GetNControlVars(); iCV++){
        const auto& CV_name = config->GetControllingVariableName(iCV);
        AddHistoryOutput("RMS_"+CV_name, "rms["+CV_name+"]",ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean squared residual of " + CV_name + " controlling variable equation.", HistoryFieldType::RESIDUAL);
      }

      /*--- auxiliary species transport ---*/
      for (auto i_scalar = 0u; i_scalar < config->GetNUserScalars(); i_scalar++){
        const auto& scalar_name = config->GetUserScalarName(i_scalar);
        AddHistoryOutput("RMS_"+scalar_name, "rms["+scalar_name+"]", ScreenOutputFormat::FIXED  , "RMS_RES", "Root-mean squared residual of the "+scalar_name+" mass fraction equation." , HistoryFieldType::RESIDUAL);
      }
      break;
    }
    case SPECIES_MODEL::NONE: break;
  }
}

void CFlowOutput::AddHistoryOutputFields_ScalarMAX_RES(const CConfig* config) {

  switch (TurbModelFamily(config->GetKind_Turb_Model())) {
    case TURB_FAMILY::SA:
      /// DESCRIPTION: Maximum residual of nu tilde (SA model).
      AddHistoryOutput("MAX_NU_TILDE", "max[nu]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of nu tilde (SA model).", HistoryFieldType::RESIDUAL);
      break;

    case TURB_FAMILY::KW:
      /// DESCRIPTION: Maximum residual of kinetic energy (SST model).
      AddHistoryOutput("MAX_TKE", "max[k]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of kinetic energy (SST model).", HistoryFieldType::RESIDUAL);
      /// DESCRIPTION: Maximum residual of the dissipation (SST model).
      AddHistoryOutput("MAX_DISSIPATION", "max[w]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of dissipation (SST model).", HistoryFieldType::RESIDUAL);
      break;

    case TURB_FAMILY::NONE:
      break;
  }

  switch (config->GetKind_Trans_Model()) {

    case TURB_TRANS_MODEL::LM:
      /// DESCRIPTION: Maximum residual of the intermittency (LM model).
      AddHistoryOutput("MAX_INTERMITTENCY", "max[LM_1]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the intermittency (LM model).", HistoryFieldType::RESIDUAL);
      /// DESCRIPTION: Maximum residual of the momentum thickness Reynolds number (LM model).
      AddHistoryOutput("MAX_RE_THETA_T", "max[LM_2]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the momentum thickness Reynolds number (LM model).", HistoryFieldType::RESIDUAL);
      break;

    case TURB_TRANS_MODEL::NONE:
      break;
  }

  switch (config->GetKind_Species_Model()) {
    case SPECIES_MODEL::SPECIES_TRANSPORT: {
      for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
        AddHistoryOutput("MAX_SPECIES_" + std::to_string(iVar), "max[rho*Y_" + std::to_string(iVar)+"]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of transported species.", HistoryFieldType::RESIDUAL);
      }
      break;
    }
    case SPECIES_MODEL::FLAMELET: {
      /*--- Controlling variable transport. ---*/
      for (auto iCV=0u; iCV < config->GetNControlVars(); iCV++){
        const auto& cv_name = config->GetControllingVariableName(iCV);
        AddHistoryOutput("MAX_" + cv_name, "max[" + cv_name + "]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the " + cv_name + " equation.", HistoryFieldType::RESIDUAL);
      }

      /*--- auxiliary species transport ---*/
      for (auto i_scalar = 0u; i_scalar < config->GetNUserScalars(); i_scalar++){
        const auto& scalar_name = config->GetUserScalarName(i_scalar);
        AddHistoryOutput("MAX_" + scalar_name, "max[" + scalar_name + "]", ScreenOutputFormat::FIXED  , "MAX_RES", "Maximum residual of the " + scalar_name + " mass fraction equation." , HistoryFieldType::RESIDUAL);
      }
      break;
    }
    case SPECIES_MODEL::NONE: break;
  }
}

void CFlowOutput::AddHistoryOutputFields_ScalarBGS_RES(const CConfig* config) {
  if (!multiZone) return;

  switch (TurbModelFamily(config->GetKind_Turb_Model())) {
    case TURB_FAMILY::SA:
      /// DESCRIPTION: Maximum residual of nu tilde (SA model).
      AddHistoryOutput("BGS_NU_TILDE", "bgs[nu]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of nu tilde (SA model).", HistoryFieldType::RESIDUAL);
      break;

    case TURB_FAMILY::KW:
      /// DESCRIPTION: Maximum residual of kinetic energy (SST model).
      AddHistoryOutput("BGS_TKE", "bgs[k]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of kinetic energy (SST model).", HistoryFieldType::RESIDUAL);
      /// DESCRIPTION: Maximum residual of the dissipation (SST model).
      AddHistoryOutput("BGS_DISSIPATION", "bgs[w]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of dissipation (SST model).", HistoryFieldType::RESIDUAL);
      break;

    case TURB_FAMILY::NONE: break;
  }

  switch (config->GetKind_Trans_Model()) {
    case TURB_TRANS_MODEL::LM:
      /// DESCRIPTION: Maximum residual of the intermittency (LM model).
      AddHistoryOutput("BGS_INTERMITTENCY", "bgs[LM_1]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the intermittency (LM model).", HistoryFieldType::RESIDUAL);
      /// DESCRIPTION: Maximum residual of the momentum thickness Reynolds number (LM model).
      AddHistoryOutput("BGS_RE_THETA_T", "bgs[LM_2]",  ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the momentum thickness Reynolds number (LM model).", HistoryFieldType::RESIDUAL);
      break;

    case TURB_TRANS_MODEL::NONE: break;
  }

  switch (config->GetKind_Species_Model()) {
    case SPECIES_MODEL::SPECIES_TRANSPORT: {
      for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
        AddHistoryOutput("BGS_SPECIES_" + std::to_string(iVar), "bgs[rho*Y_" + std::to_string(iVar)+"]", ScreenOutputFormat::FIXED, "BGS_RES", "Maximum residual of transported species.", HistoryFieldType::RESIDUAL);
      }
      break;
    }
    case SPECIES_MODEL::FLAMELET: {
      /*--- Controlling variable transport. ---*/
      for (auto iCV=0u; iCV < config->GetNControlVars(); iCV++){
        const auto& cv_name = config->GetControllingVariableName(iCV);
        AddHistoryOutput("BGS_" + cv_name, "bgs[" + cv_name + "]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the " + cv_name + " controlling variable equation.", HistoryFieldType::RESIDUAL);
      }

      /*--- auxiliary species transport ---*/
      for (auto i_scalar = 0u; i_scalar < config->GetNUserScalars(); i_scalar++){
        const auto& scalar_name = config->GetUserScalarName(i_scalar);
        AddHistoryOutput("BGS_"+scalar_name, "bgs["+scalar_name+"]", ScreenOutputFormat::FIXED  , "BGS_RES", "BGS residual of the "+scalar_name+" mass fraction equation." , HistoryFieldType::RESIDUAL);
      }
      break;
    }
    case SPECIES_MODEL::NONE: break;
  }
}

void CFlowOutput::AddHistoryOutputFieldsScalarLinsol(const CConfig* config) {
  if (config->GetKind_Turb_Model() != TURB_MODEL::NONE) {
    AddHistoryOutput("LINSOL_ITER_TURB", "LinSolIterTurb", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver for turbulence solver.");
    AddHistoryOutput("LINSOL_RESIDUAL_TURB", "LinSolResTurb", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver for turbulence solver.");
  }

  if (config->GetKind_Trans_Model() != TURB_TRANS_MODEL::NONE) {
    AddHistoryOutput("LINSOL_ITER_TRANS", "LinSolIterTrans", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver for transition solver.");
    AddHistoryOutput("LINSOL_RESIDUAL_TRANS", "LinSolResTrans", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver for transition solver.");
  }

  switch (config->GetKind_Species_Model()) {
    case SPECIES_MODEL::SPECIES_TRANSPORT: {
      AddHistoryOutput("LINSOL_ITER_SPECIES", "LinSolIterSpecies", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver for species solver.");
      AddHistoryOutput("LINSOL_RESIDUAL_SPECIES", "LinSolResSpecies", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver for species solver.");
      break;
    }
    case SPECIES_MODEL::FLAMELET: {
      AddHistoryOutput("LINSOL_ITER_FLAMELET", "LinSolIterSpecies", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver for flamelet solver.");
      AddHistoryOutput("LINSOL_RESIDUAL_FLAMELET", "LinSolResSpecies", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver for flamelet solver.");
      break;
    }
    case SPECIES_MODEL::NONE: break;
  }
}
// clang-format on

void CFlowOutput::LoadHistoryDataScalar(const CConfig* config, const CSolver* const* solver) {

  switch (TurbModelFamily(config->GetKind_Turb_Model())) {
    case TURB_FAMILY::SA:
      SetHistoryOutputValue("RMS_NU_TILDE", log10(solver[TURB_SOL]->GetRes_RMS(0)));
      SetHistoryOutputValue("MAX_NU_TILDE", log10(solver[TURB_SOL]->GetRes_Max(0)));
      if (multiZone) {
        SetHistoryOutputValue("BGS_NU_TILDE", log10(solver[TURB_SOL]->GetRes_BGS(0)));
      }
      break;

    case TURB_FAMILY::KW:
      SetHistoryOutputValue("RMS_TKE", log10(solver[TURB_SOL]->GetRes_RMS(0)));
      SetHistoryOutputValue("RMS_DISSIPATION",log10(solver[TURB_SOL]->GetRes_RMS(1)));
      SetHistoryOutputValue("MAX_TKE", log10(solver[TURB_SOL]->GetRes_Max(0)));
      SetHistoryOutputValue("MAX_DISSIPATION", log10(solver[TURB_SOL]->GetRes_Max(1)));
      if (multiZone) {
        SetHistoryOutputValue("BGS_TKE", log10(solver[TURB_SOL]->GetRes_BGS(0)));
        SetHistoryOutputValue("BGS_DISSIPATION", log10(solver[TURB_SOL]->GetRes_BGS(1)));
      }
      break;

    case TURB_FAMILY::NONE: break;
  }

  if (config->GetKind_Turb_Model() != TURB_MODEL::NONE) {
    SetHistoryOutputValue("LINSOL_ITER_TURB", solver[TURB_SOL]->GetIterLinSolver());
    SetHistoryOutputValue("LINSOL_RESIDUAL_TURB", log10(solver[TURB_SOL]->GetResLinSolver()));
  }

  switch (config->GetKind_Trans_Model()) {
    case TURB_TRANS_MODEL::LM:
      SetHistoryOutputValue("RMS_INTERMITTENCY", log10(solver[TRANS_SOL]->GetRes_RMS(0)));
      SetHistoryOutputValue("RMS_RE_THETA_T",log10(solver[TRANS_SOL]->GetRes_RMS(1)));
      SetHistoryOutputValue("MAX_INTERMITTENCY", log10(solver[TRANS_SOL]->GetRes_Max(0)));
      SetHistoryOutputValue("MAX_RE_THETA_T", log10(solver[TRANS_SOL]->GetRes_Max(1)));
      if (multiZone) {
        SetHistoryOutputValue("BGS_INTERMITTENCY", log10(solver[TRANS_SOL]->GetRes_BGS(0)));
        SetHistoryOutputValue("BGS_RE_THETA_T", log10(solver[TRANS_SOL]->GetRes_BGS(1)));
      }
      SetHistoryOutputValue("LINSOL_ITER_TRANS", solver[TRANS_SOL]->GetIterLinSolver());
      SetHistoryOutputValue("LINSOL_RESIDUAL_TRANS", log10(solver[TRANS_SOL]->GetResLinSolver()));
      break;

    case TURB_TRANS_MODEL::NONE: break;
  }

  switch(config->GetKind_Species_Model()) {
    case SPECIES_MODEL::SPECIES_TRANSPORT: {
      for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
        SetHistoryOutputValue("RMS_SPECIES_" + std::to_string(iVar), log10(solver[SPECIES_SOL]->GetRes_RMS(iVar)));
        SetHistoryOutputValue("MAX_SPECIES_" + std::to_string(iVar), log10(solver[SPECIES_SOL]->GetRes_Max(iVar)));
        if (multiZone) {
          SetHistoryOutputValue("BGS_SPECIES_" + std::to_string(iVar), log10(solver[SPECIES_SOL]->GetRes_BGS(iVar)));
        }
      }
      SetHistoryOutputValue("LINSOL_ITER_SPECIES", solver[SPECIES_SOL]->GetIterLinSolver());
      SetHistoryOutputValue("LINSOL_RESIDUAL_SPECIES", log10(solver[SPECIES_SOL]->GetResLinSolver()));
      break;
    }

    case SPECIES_MODEL::FLAMELET: {
      /*--- Controlling variable transport. ---*/
      for (auto iCV=0u; iCV < config->GetNControlVars(); iCV++){
        const auto& cv_name = config->GetControllingVariableName(iCV);
        SetHistoryOutputValue("RMS_" + cv_name, log10(solver[SPECIES_SOL]->GetRes_RMS(iCV)));
        SetHistoryOutputValue("MAX_" + cv_name, log10(solver[SPECIES_SOL]->GetRes_Max(iCV)));
      }
      /*--- auxiliary species transport ---*/
      for (unsigned short iReactant=0; iReactant<config->GetNUserScalars(); iReactant++){
        const auto& species_name = config->GetUserScalarName(iReactant);
        SetHistoryOutputValue("RMS_" + species_name, log10(solver[SPECIES_SOL]->GetRes_RMS(config->GetNControlVars() + iReactant)));
        SetHistoryOutputValue("MAX_" + species_name, log10(solver[SPECIES_SOL]->GetRes_Max(config->GetNControlVars() + iReactant)));
        if (multiZone) {
          SetHistoryOutputValue("BGS_" + species_name, log10(solver[SPECIES_SOL]->GetRes_BGS(config->GetNControlVars() + iReactant)));
        }
      }

      SetHistoryOutputValue("LINSOL_ITER_FLAMELET", solver[SPECIES_SOL]->GetIterLinSolver());
      SetHistoryOutputValue("LINSOL_RESIDUAL_FLAMELET", log10(solver[SPECIES_SOL]->GetResLinSolver()));
      break;
    }

    case SPECIES_MODEL::NONE: break;
  }
}

void CFlowOutput::SetVolumeOutputFieldsScalarSolution(const CConfig* config){
  /*--- Only place outputs of the "SOLUTION" group here. ---*/

  switch (TurbModelFamily(config->GetKind_Turb_Model())) {
    case TURB_FAMILY::SA:
      AddVolumeOutput("NU_TILDE", "Nu_Tilde", "SOLUTION", "Spalart-Allmaras variable");
      break;

    case TURB_FAMILY::KW:
      AddVolumeOutput("TKE", "Turb_Kin_Energy", "SOLUTION", "Turbulent kinetic energy");
      AddVolumeOutput("DISSIPATION", "Omega", "SOLUTION", "Rate of dissipation");
      break;

    case TURB_FAMILY::NONE:
      break;
  }

  switch (config->GetKind_Trans_Model()) {
    case TURB_TRANS_MODEL::LM:
      AddVolumeOutput("INTERMITTENCY", "LM_gamma", "SOLUTION", "LM intermittency");
      AddVolumeOutput("RE_THETA_T", "LM_Re_t", "SOLUTION", "LM RE_THETA_T");
      break;

    case TURB_TRANS_MODEL::NONE:
      break;
  }

  switch (config->GetKind_Species_Model()) {
    case SPECIES_MODEL::SPECIES_TRANSPORT:
      for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++){
        AddVolumeOutput("SPECIES_" + std::to_string(iVar), "Species_" + std::to_string(iVar), "SOLUTION", "Species_" + std::to_string(iVar) + " mass fraction");
      }
      break;
    case SPECIES_MODEL::FLAMELET:
      /*--- Controlling variables. ---*/
      for (auto iCV=0u; iCV<config->GetNControlVars(); iCV++) {
        const auto& cv_name = config->GetControllingVariableName(iCV);
        AddVolumeOutput(cv_name, cv_name, "SOLUTION", cv_name + " solution.");
      }
      /*--- auxiliary species ---*/
      for (auto iReactant=0u; iReactant<config->GetNUserScalars(); iReactant++) {
        const auto& species_name = config->GetUserScalarName(iReactant);
        AddVolumeOutput(species_name, species_name, "SOLUTION", species_name + "Mass fraction solution");
      }

      break;
    case SPECIES_MODEL::NONE:
      break;
  }
}

void CFlowOutput::SetVolumeOutputFieldsScalarResidual(const CConfig* config) {
  /*--- Only place outputs of the "RESIDUAL" group here. ---*/

  switch (TurbModelFamily(config->GetKind_Turb_Model())){
    case TURB_FAMILY::SA:
      AddVolumeOutput("RES_NU_TILDE", "Residual_Nu_Tilde", "RESIDUAL", "Residual of the Spalart-Allmaras variable");
      break;

    case TURB_FAMILY::KW:
      AddVolumeOutput("RES_TKE", "Residual_TKE", "RESIDUAL", "Residual of turbulent kinetic energy");
      AddVolumeOutput("RES_DISSIPATION", "Residual_Omega", "RESIDUAL", "Residual of the rate of dissipation");
      break;

    case TURB_FAMILY::NONE:
      break;
  }

  switch (config->GetKind_Species_Model()) {
    case SPECIES_MODEL::SPECIES_TRANSPORT:
      for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++){
        AddVolumeOutput("RES_SPECIES_" + std::to_string(iVar), "Residual_Species_" + std::to_string(iVar), "RESIDUAL", "Residual of the transported species " + std::to_string(iVar));
      }
      break;
    case SPECIES_MODEL::FLAMELET:
      /*--- Residuals for controlling variable transport equations. ---*/
      for (auto iCV=0u; iCV<config->GetNControlVars(); iCV++) {
        const auto& cv_name = config->GetControllingVariableName(iCV);
        AddVolumeOutput("RES_"+cv_name, "Residual_"+cv_name, "RESIDUAL", "Residual of " + cv_name + " controlling variable.");
      }
      /*--- residuals for auxiliary species transport equations ---*/
      for (unsigned short iReactant=0; iReactant<config->GetNUserScalars(); iReactant++){
        const auto& species_name = config->GetUserScalarName(iReactant);
        AddVolumeOutput("RES_" + species_name, "Residual_" + species_name, "RESIDUAL", "Residual of the " + species_name + " equation");
      }
      break;
    case SPECIES_MODEL::NONE:
      break;
  }

  switch (config->GetKind_Trans_Model()) {
    case TURB_TRANS_MODEL::LM:
      AddVolumeOutput("RES_INTERMITTENCY", "Residual_LM_intermittency", "RESIDUAL", "Residual of LM intermittency");
      AddVolumeOutput("RES_RE_THETA_T", "Residual_LM_RE_THETA_T", "RESIDUAL", "Residual of LM RE_THETA_T");
      break;

    case TURB_TRANS_MODEL::NONE:
      break;
  }
}


void CFlowOutput::SetVolumeOutputFieldsScalarLimiter(const CConfig* config) {
  /*--- Only place outputs of the "SOLUTION" group for species transport here. ---*/


  if (config->GetKind_SlopeLimit_Turb() != LIMITER::NONE) {
    switch (TurbModelFamily(config->GetKind_Turb_Model())) {
      case TURB_FAMILY::SA:
        AddVolumeOutput("LIMITER_NU_TILDE", "Limiter_Nu_Tilde", "LIMITER", "Limiter value of the Spalart-Allmaras variable");
        break;

      case TURB_FAMILY::KW:
        AddVolumeOutput("LIMITER_TKE", "Limiter_TKE", "LIMITER", "Limiter value of turb. kinetic energy");
        AddVolumeOutput("LIMITER_DISSIPATION", "Limiter_Omega", "LIMITER", "Limiter value of dissipation rate");
        break;

      case TURB_FAMILY::NONE:
        break;
    }
  }

  if (config->GetKind_SlopeLimit_Species() != LIMITER::NONE) {
    switch (config->GetKind_Species_Model()) {
      case SPECIES_MODEL::SPECIES_TRANSPORT:
        for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++)
          AddVolumeOutput("LIMITER_SPECIES_" + std::to_string(iVar), "Limiter_Species_" + std::to_string(iVar), "LIMITER", "Limiter value of the transported species " + std::to_string(iVar));
      break;
      case SPECIES_MODEL::FLAMELET:
        /*--- Limiter for controlling variables transport. ---*/
        for (auto iCV=0u; iCV < config->GetNControlVars(); iCV++) {
          const auto& cv_name = config->GetControllingVariableName(iCV);
          AddVolumeOutput("LIMITER_" + cv_name, "Limiter_" + cv_name, "LIMITER", "Limiter of " + cv_name + " controlling variable.");
        }
        /*--- limiter for auxiliary species transport ---*/
        for (unsigned short iReactant=0; iReactant < config->GetNUserScalars(); iReactant++) {
          const auto& species_name = config->GetUserScalarName(iReactant);
          AddVolumeOutput("LIMITER_" + species_name, "LIMITER_" + species_name, "LIMITER", "Limiter value for the " + species_name + " equation");
        }
      break;
      default:
        break;
    }
  }
}

void CFlowOutput::SetVolumeOutputFieldsScalarPrimitive(const CConfig* config) {
  /*--- Only place outputs of the "PRIMITIVE" group for scalar transport here. ---*/

  switch (config->GetKind_Species_Model()) {
    case SPECIES_MODEL::SPECIES_TRANSPORT:
      for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++){
        AddVolumeOutput("DIFFUSIVITY_" + std::to_string(iVar), "Diffusivity_" + std::to_string(iVar), "PRIMITIVE", "Diffusivity of the transported species " + std::to_string(iVar));
      }
      break;
    default:
      break;
  }

  switch (config->GetKind_Trans_Model()) {
    case TURB_TRANS_MODEL::LM:
      AddVolumeOutput("INTERMITTENCY_SEP", "LM_gamma_sep", "PRIMITIVE", "LM intermittency");
      AddVolumeOutput("INTERMITTENCY_EFF", "LM_gamma_eff", "PRIMITIVE", "LM RE_THETA_T");
      AddVolumeOutput("TURB_INDEX", "Turb_index", "PRIMITIVE", "Turbulence index");
      break;

    case TURB_TRANS_MODEL::NONE:
      break;
  }

  if (config->GetKind_Turb_Model() != TURB_MODEL::NONE) {
    AddVolumeOutput("EDDY_VISCOSITY", "Eddy_Viscosity", "PRIMITIVE", "Turbulent eddy viscosity");
  }

}

void CFlowOutput::SetVolumeOutputFieldsScalarSource(const CConfig* config) {
  /*--- Only place outputs of the "SOURCE" group for scalar transport here. ---*/

  switch (config->GetKind_Species_Model()) {
    case SPECIES_MODEL::FLAMELET:
      for (auto iCV=0u; iCV < config->GetNControlVars(); iCV++) {
        const auto& cv_source_name = config->GetControllingVariableSourceName(iCV);
        const auto& cv_name = config->GetControllingVariableName(iCV);
        if (cv_source_name.compare("NULL") != 0)
          AddVolumeOutput("SOURCE_"+cv_name, "Source_" + cv_name, "SOURCE", "Source " + cv_name);
      }
      /*--- no source term for enthalpy ---*/
      /*--- auxiliary species source terms ---*/
      for (auto iReactant=0u; iReactant<config->GetNUserScalars(); iReactant++) {
        const auto& species_name = config->GetUserScalarName(iReactant);
        AddVolumeOutput("SOURCE_" + species_name, "Source_" + species_name, "SOURCE", "Source " + species_name);
      }
      break;
    default:
      break;
  }
}


void CFlowOutput::SetVolumeOutputFieldsScalarLookup(const CConfig* config) {
  /*--- Only place outputs of the "LOOKUP" group for scalar transport here. ---*/

  switch (config->GetKind_Species_Model()) {
    case SPECIES_MODEL::FLAMELET:
      for (auto i_lookup = 0u; i_lookup < config->GetNLookups(); ++i_lookup) {
        string strname1 = "lookup_" + config->GetLookupName(i_lookup);
        AddVolumeOutput(config->GetLookupName(i_lookup), strname1,"LOOKUP", config->GetLookupName(i_lookup));
      }
      AddVolumeOutput("TABLE_MISSES"       , "Table_misses"       , "LOOKUP", "Lookup table misses");
      break;
    default:
      break;
  }
}

void CFlowOutput::SetVolumeOutputFieldsScalarMisc(const CConfig* config) {
  /*--- Only place outputs of the group for scalar transport here that do not fit in other categories. ---*/

  if (config->GetSAParsedOptions().bc) {
    AddVolumeOutput("INTERMITTENCY", "gamma_BC", "INTERMITTENCY", "Intermittency");
  }

  // Hybrid RANS-LES
  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES) {
    AddVolumeOutput("DES_LENGTHSCALE", "DES_LengthScale", "DDES", "DES length scale value");
    AddVolumeOutput("WALL_DISTANCE", "Wall_Distance", "DDES", "Wall distance value");
  }

  if (config->GetViscous()) {
    if (nDim == 3) {
      AddVolumeOutput("VORTICITY_X", "Vorticity_x", "VORTEX_IDENTIFICATION", "x-component of the vorticity vector");
      AddVolumeOutput("VORTICITY_Y", "Vorticity_y", "VORTEX_IDENTIFICATION", "y-component of the vorticity vector");
      AddVolumeOutput("VORTICITY_Z", "Vorticity_z", "VORTEX_IDENTIFICATION", "z-component of the vorticity vector");
    } else {
      AddVolumeOutput("VORTICITY", "Vorticity", "VORTEX_IDENTIFICATION", "Value of the vorticity");
    }
    AddVolumeOutput("Q_CRITERION", "Q_Criterion", "VORTEX_IDENTIFICATION", "Value of the Q-Criterion");
  }

  // Timestep info
  AddVolumeOutput("DELTA_TIME", "Delta_Time", "TIMESTEP", "Value of the local timestep for the flow variables");
  AddVolumeOutput("CFL", "CFL", "TIMESTEP", "Value of the local CFL for the flow variables");
  if (config->GetKind_Turb_Model() != TURB_MODEL::NONE)
  {
    AddVolumeOutput("TURB_DELTA_TIME", "Turb_Delta_Time", "TIMESTEP", "Value of the local timestep for the turbulence variables");
    AddVolumeOutput("TURB_CFL", "Turb_CFL", "TIMESTEP", "Value of the local CFL for the turbulence variables");
  }
}

void CFlowOutput::LoadVolumeDataScalar(const CConfig* config, const CSolver* const* solver, const CGeometry* geometry,
                                        const unsigned long iPoint) {
  const auto* turb_solver = solver[TURB_SOL];
  const auto* trans_solver = solver[TRANS_SOL];
  const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();
  const auto* Node_Turb = (config->GetKind_Turb_Model() != TURB_MODEL::NONE) ? turb_solver->GetNodes() : nullptr;
  const auto* Node_Trans = (config->GetKind_Trans_Model() != TURB_TRANS_MODEL::NONE) ? trans_solver->GetNodes() : nullptr;
  const auto* Node_Geo = geometry->nodes;

  SetVolumeOutputValue("DELTA_TIME", iPoint, Node_Flow->GetDelta_Time(iPoint));
  SetVolumeOutputValue("CFL", iPoint, Node_Flow->GetLocalCFL(iPoint));

  if (config->GetViscous()) {
    if (nDim == 3){
      SetVolumeOutputValue("VORTICITY_X", iPoint, Node_Flow->GetVorticity(iPoint)[0]);
      SetVolumeOutputValue("VORTICITY_Y", iPoint, Node_Flow->GetVorticity(iPoint)[1]);
      SetVolumeOutputValue("VORTICITY_Z", iPoint, Node_Flow->GetVorticity(iPoint)[2]);
    } else {
      SetVolumeOutputValue("VORTICITY", iPoint, Node_Flow->GetVorticity(iPoint)[2]);
    }
    SetVolumeOutputValue("Q_CRITERION", iPoint, GetQCriterion(Node_Flow->GetVelocityGradient(iPoint)));
  }

  const bool limiter = (config->GetKind_SlopeLimit_Turb() != LIMITER::NONE);

  switch (TurbModelFamily(config->GetKind_Turb_Model())) {
    case TURB_FAMILY::SA:
      SetVolumeOutputValue("NU_TILDE", iPoint, Node_Turb->GetSolution(iPoint, 0));
      SetVolumeOutputValue("RES_NU_TILDE", iPoint, turb_solver->LinSysRes(iPoint, 0));
      if (limiter) {
        SetVolumeOutputValue("LIMITER_NU_TILDE", iPoint, Node_Turb->GetLimiter(iPoint, 0));
      }
      break;

    case TURB_FAMILY::KW:
      SetVolumeOutputValue("TKE", iPoint, Node_Turb->GetSolution(iPoint, 0));
      SetVolumeOutputValue("DISSIPATION", iPoint, Node_Turb->GetSolution(iPoint, 1));
      SetVolumeOutputValue("RES_TKE", iPoint, turb_solver->LinSysRes(iPoint, 0));
      SetVolumeOutputValue("RES_DISSIPATION", iPoint, turb_solver->LinSysRes(iPoint, 1));
      if (limiter) {
        SetVolumeOutputValue("LIMITER_TKE", iPoint, Node_Turb->GetLimiter(iPoint, 0));
        SetVolumeOutputValue("LIMITER_DISSIPATION", iPoint, Node_Turb->GetLimiter(iPoint, 1));
      }
      break;

    case TURB_FAMILY::NONE: break;
  }

  /*--- If we got here a turbulence model is being used, therefore there is eddy viscosity. ---*/
  if (config->GetKind_Turb_Model() != TURB_MODEL::NONE) {
    SetVolumeOutputValue("EDDY_VISCOSITY", iPoint, Node_Flow->GetEddyViscosity(iPoint));
    SetVolumeOutputValue("TURB_DELTA_TIME", iPoint, Node_Turb->GetDelta_Time(iPoint));
    SetVolumeOutputValue("TURB_CFL", iPoint, Node_Turb->GetLocalCFL(iPoint));
  }

  if (config->GetSAParsedOptions().bc) {
    SetVolumeOutputValue("INTERMITTENCY", iPoint, Node_Turb->GetIntermittencyEff(iPoint));
  }

  switch (config->GetKind_Trans_Model()) {
    case TURB_TRANS_MODEL::LM:
      SetVolumeOutputValue("INTERMITTENCY", iPoint, Node_Trans->GetSolution(iPoint, 0));
      SetVolumeOutputValue("RE_THETA_T", iPoint, Node_Trans->GetSolution(iPoint, 1));
      SetVolumeOutputValue("INTERMITTENCY_SEP", iPoint, Node_Trans->GetIntermittencySep(iPoint));
      SetVolumeOutputValue("INTERMITTENCY_EFF", iPoint, Node_Trans->GetIntermittencyEff(iPoint));
      SetVolumeOutputValue("TURB_INDEX", iPoint, Node_Turb->GetTurbIndex(iPoint));
      SetVolumeOutputValue("RES_INTERMITTENCY", iPoint, trans_solver->LinSysRes(iPoint, 0));
      SetVolumeOutputValue("RES_RE_THETA_T", iPoint, trans_solver->LinSysRes(iPoint, 1));
      break;

    case TURB_TRANS_MODEL::NONE: break;
  }

  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES) {
    SetVolumeOutputValue("DES_LENGTHSCALE", iPoint, Node_Flow->GetDES_LengthScale(iPoint));
    SetVolumeOutputValue("WALL_DISTANCE", iPoint, Node_Geo->GetWall_Distance(iPoint));
  }

  switch (config->GetKind_Species_Model()) {

    case SPECIES_MODEL::SPECIES_TRANSPORT: {
      const auto Node_Species = solver[SPECIES_SOL]->GetNodes();
      for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
        SetVolumeOutputValue("SPECIES_" + std::to_string(iVar), iPoint, Node_Species->GetSolution(iPoint, iVar));
        SetVolumeOutputValue("RES_SPECIES_" + std::to_string(iVar), iPoint, solver[SPECIES_SOL]->LinSysRes(iPoint, iVar));
        SetVolumeOutputValue("DIFFUSIVITY_"+ std::to_string(iVar), iPoint, Node_Species->GetDiffusivity(iPoint,iVar));
        if (config->GetKind_SlopeLimit_Species() != LIMITER::NONE)
          SetVolumeOutputValue("LIMITER_SPECIES_" + std::to_string(iVar), iPoint, Node_Species->GetLimiter(iPoint, iVar));
      }
      break;
    }

    case SPECIES_MODEL::FLAMELET: {
      const auto Node_Species = solver[SPECIES_SOL]->GetNodes();

      /*--- Controlling variables transport equations. ---*/
      for (auto iCV=0u; iCV < config->GetNControlVars(); iCV++) {
        const auto& cv_name = config->GetControllingVariableName(iCV);
        SetVolumeOutputValue(cv_name, iPoint, Node_Species->GetSolution(iPoint, iCV));
        SetVolumeOutputValue("RES_" + cv_name, iPoint, solver[SPECIES_SOL]->LinSysRes(iPoint, iCV));
        const auto& source_name = config->GetControllingVariableSourceName(iCV);
        if (source_name.compare("NULL") != 0)
          SetVolumeOutputValue("SOURCE_" + cv_name, iPoint, Node_Species->GetScalarSources(iPoint)[iCV]);
      }
      /*--- auxiliary species transport equations ---*/
      for (unsigned short i_scalar=0; i_scalar<config->GetNUserScalars(); i_scalar++) {
        const auto& scalar_name = config->GetUserScalarName(i_scalar);
        SetVolumeOutputValue(scalar_name, iPoint, Node_Species->GetSolution(iPoint, config->GetNControlVars() + i_scalar));
        SetVolumeOutputValue("SOURCE_" + scalar_name, iPoint, Node_Species->GetScalarSources(iPoint)[config->GetNControlVars() + i_scalar]);
        SetVolumeOutputValue("RES_" + scalar_name, iPoint, solver[SPECIES_SOL]->LinSysRes(iPoint, config->GetNControlVars() + i_scalar));
      }

      if (config->GetKind_SlopeLimit_Species() != LIMITER::NONE) {
        /*--- Limiter for controlling variable transport equations. ---*/
        for (auto iCV=0u; iCV<config->GetNControlVars(); iCV++) {
          const auto& cv_name = config->GetControllingVariableName(iCV);
          SetVolumeOutputValue("LIMITER_" + cv_name, iPoint, Node_Species->GetLimiter(iPoint, iCV));
        }
        /*--- limiter for auxiliary species transport equations ---*/
        for (unsigned short i_scalar=0; i_scalar<config->GetNUserScalars(); i_scalar++) {
          const auto& scalar_name = config->GetUserScalarName(i_scalar);
          SetVolumeOutputValue("LIMITER_" + scalar_name, iPoint, Node_Species->GetLimiter(iPoint, config->GetNControlVars() + i_scalar));
        }
      }

      /*--- variables that we look up from the LUT ---*/
      for (int i_lookup = 0; i_lookup < config->GetNLookups(); ++i_lookup) {
        if (config->GetLookupName(i_lookup)!="NULL")
          SetVolumeOutputValue(config->GetLookupName(i_lookup), iPoint, Node_Species->GetScalarLookups(iPoint)[i_lookup]);
      }

      SetVolumeOutputValue("TABLE_MISSES", iPoint, Node_Species->GetTableMisses(iPoint));

      break;
    }
    case SPECIES_MODEL::NONE: break;
  }
}

void CFlowOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){

  if (!config->GetViscous_Wall(iMarker)) return;

  const auto heat_sol = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE) &&
                         config->GetWeakly_Coupled_Heat() ? HEAT_SOL : FLOW_SOL;

  SetVolumeOutputValue("SKIN_FRICTION-X", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 0));
  SetVolumeOutputValue("SKIN_FRICTION-Y", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 1));
  if (nDim == 3)
    SetVolumeOutputValue("SKIN_FRICTION-Z", iPoint, solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 2));
  SetVolumeOutputValue("HEAT_FLUX", iPoint, solver[heat_sol]->GetHeatFlux(iMarker, iVertex));
  SetVolumeOutputValue("Y_PLUS", iPoint, solver[FLOW_SOL]->GetYPlus(iMarker, iVertex));
}

void CFlowOutput::AddAerodynamicCoefficients(const CConfig* config) {

  /// BEGIN_GROUP: AERO_COEFF, DESCRIPTION: Sum of the aerodynamic coefficients and forces on all surfaces (markers) set with MARKER_MONITORING.
  /// DESCRIPTION: Reference force for aerodynamic coefficients
  AddHistoryOutput("REFERENCE_FORCE", "RefForce", ScreenOutputFormat::FIXED, "AERO_COEFF", "Reference force used to compute aerodynamic coefficients", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Drag coefficient
  AddHistoryOutput("DRAG",       "CD",   ScreenOutputFormat::FIXED, "AERO_COEFF", "Total drag coefficient on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Lift coefficient
  AddHistoryOutput("LIFT",       "CL",   ScreenOutputFormat::FIXED, "AERO_COEFF", "Total lift coefficient on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Sideforce coefficient
  AddHistoryOutput("SIDEFORCE",  "CSF",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total sideforce coefficient on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Moment around the x-axis
  AddHistoryOutput("MOMENT_X",   "CMx",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total momentum x-component on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Moment around the y-axis
  AddHistoryOutput("MOMENT_Y",   "CMy",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total momentum y-component on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Moment around the z-axis
  AddHistoryOutput("MOMENT_Z",   "CMz",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total momentum z-component on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Force in x direction
  AddHistoryOutput("FORCE_X",    "CFx",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total force x-component on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Force in y direction
  AddHistoryOutput("FORCE_Y",    "CFy",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total force y-component on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Force in z direction
  AddHistoryOutput("FORCE_Z",    "CFz",  ScreenOutputFormat::FIXED, "AERO_COEFF", "Total force z-component on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Lift-to-drag ratio
  AddHistoryOutput("EFFICIENCY", "CEff", ScreenOutputFormat::FIXED, "AERO_COEFF", "Total lift-to-drag ratio on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  /// END_GROUP

  /// BEGIN_GROUP: AERO_COEFF_SURF, DESCRIPTION: Aerodynamic coefficients and forces per surface.
  vector<string> Marker_Monitoring;
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++){
    Marker_Monitoring.push_back(config->GetMarker_Monitoring_TagBound(iMarker_Monitoring));
  }
  /// DESCRIPTION: Drag coefficient
  AddHistoryOutputPerSurface("DRAG_ON_SURFACE",       "CD",   ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Lift coefficient
  AddHistoryOutputPerSurface("LIFT_ON_SURFACE",       "CL",   ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Sideforce coefficient
  AddHistoryOutputPerSurface("SIDEFORCE_ON_SURFACE",  "CSF",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Moment around the x-axis
  AddHistoryOutputPerSurface("MOMENT-X_ON_SURFACE",   "CMx",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Moment around the y-axis
  AddHistoryOutputPerSurface("MOMENT-Y_ON_SURFACE",   "CMy",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Moment around the z-axis
  AddHistoryOutputPerSurface("MOMENT-Z_ON_SURFACE",   "CMz",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Force in x direction
  AddHistoryOutputPerSurface("FORCE-X_ON_SURFACE",    "CFx",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Force in y direction
  AddHistoryOutputPerSurface("FORCE-Y_ON_SURFACE",    "CFy",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Force in z direction
  AddHistoryOutputPerSurface("FORCE-Z_ON_SURFACE",    "CFz",  ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Lift-to-drag ratio
  AddHistoryOutputPerSurface("EFFICIENCY_ON_SURFACE", "CEff", ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
  /// END_GROUP

  /// DESCRIPTION: Angle of attack
  AddHistoryOutput("AOA", "AoA", ScreenOutputFormat::FIXED, "AOA", "Angle of attack");

  AddHistoryOutput("COMBO", "ComboObj", ScreenOutputFormat::SCIENTIFIC, "COMBO", "Combined obj. function value.", HistoryFieldType::COEFFICIENT);
}

void CFlowOutput::SetAerodynamicCoefficients(const CConfig* config, const CSolver* flow_solver){

  SetHistoryOutputValue("REFERENCE_FORCE", flow_solver->GetAeroCoeffsReferenceForce());
  SetHistoryOutputValue("DRAG", flow_solver->GetTotal_CD());
  SetHistoryOutputValue("LIFT", flow_solver->GetTotal_CL());
  if (nDim == 3)
    SetHistoryOutputValue("SIDEFORCE", flow_solver->GetTotal_CSF());
  if (nDim == 3){
    SetHistoryOutputValue("MOMENT_X", flow_solver->GetTotal_CMx());
    SetHistoryOutputValue("MOMENT_Y", flow_solver->GetTotal_CMy());
  }
  SetHistoryOutputValue("MOMENT_Z", flow_solver->GetTotal_CMz());
  SetHistoryOutputValue("FORCE_X", flow_solver->GetTotal_CFx());
  SetHistoryOutputValue("FORCE_Y", flow_solver->GetTotal_CFy());
  if (nDim == 3)
    SetHistoryOutputValue("FORCE_Z", flow_solver->GetTotal_CFz());
  SetHistoryOutputValue("EFFICIENCY", flow_solver->GetTotal_CEff());

  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    SetHistoryOutputPerSurfaceValue("DRAG_ON_SURFACE", flow_solver->GetSurface_CD(iMarker_Monitoring), iMarker_Monitoring);
    SetHistoryOutputPerSurfaceValue("LIFT_ON_SURFACE", flow_solver->GetSurface_CL(iMarker_Monitoring), iMarker_Monitoring);
    if (nDim == 3)
      SetHistoryOutputPerSurfaceValue("SIDEFORCE_ON_SURFACE", flow_solver->GetSurface_CSF(iMarker_Monitoring), iMarker_Monitoring);
    if (nDim == 3){
      SetHistoryOutputPerSurfaceValue("MOMENT-X_ON_SURFACE", flow_solver->GetSurface_CMx(iMarker_Monitoring), iMarker_Monitoring);
      SetHistoryOutputPerSurfaceValue("MOMENT-Y_ON_SURFACE", flow_solver->GetSurface_CMy(iMarker_Monitoring), iMarker_Monitoring);
    }
    SetHistoryOutputPerSurfaceValue("MOMENT-Z_ON_SURFACE", flow_solver->GetSurface_CMz(iMarker_Monitoring), iMarker_Monitoring);
    SetHistoryOutputPerSurfaceValue("FORCE-X_ON_SURFACE", flow_solver->GetSurface_CFx(iMarker_Monitoring), iMarker_Monitoring);
    SetHistoryOutputPerSurfaceValue("FORCE-Y_ON_SURFACE", flow_solver->GetSurface_CFy(iMarker_Monitoring), iMarker_Monitoring);
    if (nDim == 3)
      SetHistoryOutputPerSurfaceValue("FORCE-Z_ON_SURFACE", flow_solver->GetSurface_CFz(iMarker_Monitoring), iMarker_Monitoring);

    SetHistoryOutputPerSurfaceValue("EFFICIENCY_ON_SURFACE", flow_solver->GetSurface_CEff(iMarker_Monitoring), iMarker_Monitoring);
    if (config->GetAeroelastic_Simulation()){
      SetHistoryOutputPerSurfaceValue("PITCH", config->GetAeroelastic_pitch(iMarker_Monitoring), iMarker_Monitoring);
      SetHistoryOutputPerSurfaceValue("PLUNGE", config->GetAeroelastic_plunge(iMarker_Monitoring), iMarker_Monitoring);
    }
  }

  SetHistoryOutputValue("AOA", config->GetAoA());
}

void CFlowOutput::AddHeatCoefficients(const CConfig* config) {

  if (!config->GetViscous()) return;

  /// BEGIN_GROUP: HEAT, DESCRIPTION: Heat coefficients on all surfaces set with MARKER_MONITORING.
  /// DESCRIPTION: Total heatflux
  AddHistoryOutput("TOTAL_HEATFLUX", "HF", ScreenOutputFormat::SCIENTIFIC, "HEAT", "Total heatflux on all surfaces set with MARKER_MONITORING.", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Maximal heatflux
  AddHistoryOutput("MAXIMUM_HEATFLUX", "maxHF", ScreenOutputFormat::SCIENTIFIC, "HEAT", "Maximum heatflux across all surfaces set with MARKER_MONITORING.", HistoryFieldType::COEFFICIENT);

  vector<string> Marker_Monitoring;
  Marker_Monitoring.reserve(config->GetnMarker_Monitoring());
for (auto iMarker = 0u; iMarker < config->GetnMarker_Monitoring(); iMarker++) {
    Marker_Monitoring.push_back(config->GetMarker_Monitoring_TagBound(iMarker));
  }
  /// DESCRIPTION:  Total heatflux
  AddHistoryOutputPerSurface("TOTAL_HEATFLUX_ON_SURFACE", "HF", ScreenOutputFormat::SCIENTIFIC, "HEAT_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION:  Total heatflux
  AddHistoryOutputPerSurface("MAXIMUM_HEATFLUX_ON_SURFACE", "maxHF", ScreenOutputFormat::SCIENTIFIC, "HEAT_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
  /// END_GROUP
}

void CFlowOutput::SetHeatCoefficients(const CConfig* config, const CSolver* flow_solver) {

  if (!config->GetViscous()) return;

  SetHistoryOutputValue("TOTAL_HEATFLUX", flow_solver->GetTotal_HeatFlux());
  SetHistoryOutputValue("MAXIMUM_HEATFLUX", flow_solver->GetTotal_MaxHeatFlux());

  for (auto iMarker = 0u; iMarker < config->GetnMarker_Monitoring(); iMarker++) {
    SetHistoryOutputPerSurfaceValue("TOTAL_HEATFLUX_ON_SURFACE", flow_solver->GetSurface_HF_Visc(iMarker), iMarker);
    SetHistoryOutputPerSurfaceValue("MAXIMUM_HEATFLUX_ON_SURFACE", flow_solver->GetSurface_MaxHF_Visc(iMarker), iMarker);
  }
}

void CFlowOutput::AddRotatingFrameCoefficients() {
  /// BEGIN_GROUP: ROTATING_FRAME, DESCRIPTION: Coefficients related to a rotating frame of reference.
  /// DESCRIPTION: Merit
  AddHistoryOutput("FIGURE_OF_MERIT", "CMerit", ScreenOutputFormat::SCIENTIFIC, "ROTATING_FRAME", "Thrust over torque", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: CT
  AddHistoryOutput("THRUST", "CT", ScreenOutputFormat::SCIENTIFIC, "ROTATING_FRAME", "Thrust coefficient", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: CQ
  AddHistoryOutput("TORQUE", "CQ", ScreenOutputFormat::SCIENTIFIC, "ROTATING_FRAME", "Torque coefficient", HistoryFieldType::COEFFICIENT);
  /// END_GROUP
}

void CFlowOutput::SetRotatingFrameCoefficients(const CSolver* flow_solver) {

  SetHistoryOutputValue("THRUST", flow_solver->GetTotal_CT());
  SetHistoryOutputValue("TORQUE", flow_solver->GetTotal_CQ());
  SetHistoryOutputValue("FIGURE_OF_MERIT", flow_solver->GetTotal_CMerit());
}

void CFlowOutput::AddCpInverseDesignOutput(){

  AddHistoryOutput("INVERSE_DESIGN_PRESSURE", "Cp_Diff", ScreenOutputFormat::FIXED, "CP_DIFF", "Cp difference for inverse design", HistoryFieldType::COEFFICIENT);
}

void CFlowOutput::SetCpInverseDesign(CSolver *solver, const CGeometry *geometry, const CConfig *config){

  /*--- Prepare to read the surface pressure files (CSV) ---*/

  const auto surfCp_filename = config->GetUnsteady_FileName("TargetCp", curTimeIter, ".dat");

  /*--- Read the surface pressure file, on the first inner iteration. ---*/

  ifstream Surface_file;
  Surface_file.open(surfCp_filename);

  if (!Surface_file.good()) {
    solver->SetTotal_CpDiff(0.0);
    SetHistoryOutputValue("INVERSE_DESIGN_PRESSURE", 0.0);
    return;
  }

  if ((config->GetInnerIter() == 0) || config->GetDiscrete_Adjoint()) {
    string text_line;
    getline(Surface_file, text_line);

    while (getline(Surface_file, text_line)) {
      /*--- remove commas ---*/
      for (auto& c : text_line) if (c == ',') c = ' ';
      stringstream point_line(text_line);

      /*--- parse line ---*/
      unsigned long iPointGlobal;
      su2double XCoord, YCoord, ZCoord=0, Pressure, PressureCoeff;

      point_line >> iPointGlobal >> XCoord >> YCoord;
      if (nDim == 3) point_line >> ZCoord;
      point_line >> Pressure >> PressureCoeff;

      const auto iPoint = geometry->GetGlobal_to_Local_Point(iPointGlobal);

      /*--- If the point is on this rank set the Cp to associated vertices
       *    (one point may be shared by multiple vertices). ---*/
      if (iPoint >= 0) {
        bool set = false;
        for (auto iMarker = 0u; iMarker < geometry->GetnMarker(); ++iMarker) {
          const auto iVertex = geometry->nodes->GetVertex(iPoint, iMarker);

          if (iVertex >= 0) {
            solver->SetCPressureTarget(iMarker, iVertex, PressureCoeff);
            set = true;
          }
        }
        if (!set)
          cout << "WARNING: In file " << surfCp_filename << ", point " << iPointGlobal << " is not a vertex." << endl;
      }
    }
  }

  /*--- Compute the pressure difference. ---*/

  su2double PressDiff = 0.0;

  for (auto iMarker = 0u; iMarker < geometry->GetnMarker(); ++iMarker) {

    const auto Boundary = config->GetMarker_All_KindBC(iMarker);

    if (config->GetSolid_Wall(iMarker) || (Boundary == NEARFIELD_BOUNDARY)) {
      for (auto iVertex = 0ul; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (!geometry->nodes->GetDomain(iPoint)) continue;

        const auto Cp = solver->GetCPressure(iMarker, iVertex);
        const auto CpTarget = solver->GetCPressureTarget(iMarker, iVertex);

        const auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        const auto Area = GeometryToolbox::Norm(nDim, Normal);

        PressDiff += Area * pow(CpTarget-Cp, 2);
      }
    }
  }
  su2double tmp = PressDiff;
  SU2_MPI::Allreduce(&tmp, &PressDiff, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

  /*--- Update the total Cp difference coeffient. ---*/

  solver->SetTotal_CpDiff(PressDiff);
  SetHistoryOutputValue("INVERSE_DESIGN_PRESSURE", PressDiff);

}

void CFlowOutput::AddNearfieldInverseDesignOutput(){

  AddHistoryOutput("EQUIVALENT_AREA", "CEquiv_Area", ScreenOutputFormat::SCIENTIFIC, "EQUIVALENT_AREA", "Equivalent area", HistoryFieldType::COEFFICIENT);
}

void CFlowOutput::SetNearfieldInverseDesign(CSolver *solver, const CGeometry *geometry, const CConfig *config){

  ofstream EquivArea_file;
  su2double auxXCoord, auxYCoord, auxZCoord, InverseDesign = 0.0, DeltaX,
    Coord_i, Coord_j, jp1Coord, *Coord = nullptr, MeanFunction,
    *Face_Normal = nullptr, auxArea, auxPress, jFunction, jp1Function;
  unsigned long iPoint, auxPoint, auxDomain;
  ofstream NearFieldEA_file; ifstream TargetEA_file;

  const su2double XCoordBegin_OF = config->GetEA_IntLimit(0);
  const su2double XCoordEnd_OF = config->GetEA_IntLimit(1);

  const su2double AoA = -(config->GetAoA()*PI_NUMBER/180.0);
  const su2double EAScaleFactor = config->GetEA_ScaleFactor(); // The EA Obj. Func. should be ~ force based Obj. Func.

  const su2double Mach  = config->GetMach();
  const su2double Gamma = config->GetGamma();
  const su2double Beta = sqrt(Mach*Mach-1.0);
  const su2double R_Plane = fabs(config->GetEA_IntLimit(2));
  const su2double Pressure_Inf = config->GetPressure_FreeStreamND();

  const su2double factor = 4.0*sqrt(2.0*Beta*R_Plane) / (Gamma*Pressure_Inf*Mach*Mach);

  if (rank == MASTER_NODE) cout << "Writing Equivalent Area files." << endl ;

  vector<unsigned long> Buffer_Receive_nVertex;
  if (rank == MASTER_NODE) {
    Buffer_Receive_nVertex.resize(size);
  }

  /*--- Compute the total number of points of the near-field ghost nodes ---*/

  unsigned long nLocalVertex_NearField = 0;
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (unsigned long iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->nodes->GetCoord(iPoint);

        if (geometry->nodes->GetDomain(iPoint))
          if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0))
            nLocalVertex_NearField ++;
      }

  /*--- Send Near-Field vertex information --*/
  unsigned long MaxLocalVertex_NearField, nVertex_NearField;

  SU2_MPI::Allreduce(&nLocalVertex_NearField, &nVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());
  SU2_MPI::Gather(&nLocalVertex_NearField, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex.data(), 1, MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());

  vector<su2double> Buffer_Send_Xcoord          (MaxLocalVertex_NearField, 0.0);
  vector<su2double> Buffer_Send_Ycoord          (MaxLocalVertex_NearField, 0.0);
  vector<su2double> Buffer_Send_Zcoord          (MaxLocalVertex_NearField, 0.0);
  vector<unsigned long> Buffer_Send_IdPoint     (MaxLocalVertex_NearField, 0);
  vector<su2double> Buffer_Send_Pressure        (MaxLocalVertex_NearField, 0.0);
  vector<su2double> Buffer_Send_FaceArea        (MaxLocalVertex_NearField, 0.0);

  vector<su2double> Buffer_Receive_Xcoord;
  vector<su2double> Buffer_Receive_Ycoord;
  vector<su2double> Buffer_Receive_Zcoord;
  vector<unsigned long> Buffer_Receive_IdPoint;
  vector<su2double> Buffer_Receive_Pressure;
  vector<su2double> Buffer_Receive_FaceArea;

  if (rank == MASTER_NODE) {
    Buffer_Receive_Xcoord.resize(size*MaxLocalVertex_NearField);
    Buffer_Receive_Ycoord.resize(size*MaxLocalVertex_NearField);
    Buffer_Receive_Zcoord.resize(size*MaxLocalVertex_NearField);
    Buffer_Receive_IdPoint.resize(size*MaxLocalVertex_NearField);
    Buffer_Receive_Pressure.resize(size*MaxLocalVertex_NearField);
    Buffer_Receive_FaceArea.resize(size*MaxLocalVertex_NearField);
  }

  const auto nBuffer_Xcoord = MaxLocalVertex_NearField;
  const auto nBuffer_Ycoord = MaxLocalVertex_NearField;
  const auto nBuffer_Zcoord = MaxLocalVertex_NearField;
  const auto nBuffer_IdPoint = MaxLocalVertex_NearField;
  const auto nBuffer_Pressure = MaxLocalVertex_NearField;
  const auto nBuffer_FaceArea = MaxLocalVertex_NearField;


  /*--- Copy coordinates, index points, and pressures to the auxiliar vector --*/

  nLocalVertex_NearField = 0;
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (unsigned long iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->nodes->GetCoord(iPoint);

        if (geometry->nodes->GetDomain(iPoint))
          if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) {
            Buffer_Send_IdPoint[nLocalVertex_NearField] = iPoint;
            Buffer_Send_Xcoord[nLocalVertex_NearField] = geometry->nodes->GetCoord(iPoint, 0);
            Buffer_Send_Ycoord[nLocalVertex_NearField] = geometry->nodes->GetCoord(iPoint, 1);
            if (nDim == 3) {
              Buffer_Send_Zcoord[nLocalVertex_NearField] = geometry->nodes->GetCoord(iPoint, 2);
            }
            Buffer_Send_Pressure[nLocalVertex_NearField] = solver->GetNodes()->GetPressure(iPoint);
            Buffer_Send_FaceArea[nLocalVertex_NearField] = fabs(Face_Normal[nDim-1]);
            nLocalVertex_NearField++;
          }
      }

  /*--- Send all the information --*/

  SU2_MPI::Gather(Buffer_Send_Xcoord.data(), nBuffer_Xcoord, MPI_DOUBLE, Buffer_Receive_Xcoord.data(), nBuffer_Xcoord, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
  SU2_MPI::Gather(Buffer_Send_Ycoord.data(), nBuffer_Ycoord, MPI_DOUBLE, Buffer_Receive_Ycoord.data(), nBuffer_Ycoord, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
  SU2_MPI::Gather(Buffer_Send_Zcoord.data(), nBuffer_Zcoord, MPI_DOUBLE, Buffer_Receive_Zcoord.data(), nBuffer_Zcoord, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
  SU2_MPI::Gather(Buffer_Send_IdPoint.data(), nBuffer_IdPoint, MPI_UNSIGNED_LONG, Buffer_Receive_IdPoint.data(), nBuffer_IdPoint, MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());
  SU2_MPI::Gather(Buffer_Send_Pressure.data(), nBuffer_Pressure, MPI_DOUBLE, Buffer_Receive_Pressure.data(), nBuffer_Pressure, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
  SU2_MPI::Gather(Buffer_Send_FaceArea.data(), nBuffer_FaceArea, MPI_DOUBLE, Buffer_Receive_FaceArea.data(), nBuffer_FaceArea, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());

  if (rank == MASTER_NODE) {

    vector<su2double> Xcoord(nVertex_NearField);
    vector<su2double> Ycoord(nVertex_NearField);
    vector<su2double> Zcoord(nVertex_NearField);
    vector<short> AzimuthalAngle(nVertex_NearField);
    vector<unsigned long> IdPoint(nVertex_NearField);
    vector<unsigned long> IdDomain(nVertex_NearField);
    vector<su2double> Pressure(nVertex_NearField);
    vector<su2double> FaceArea(nVertex_NearField);
    vector<su2double> EquivArea(nVertex_NearField);
    vector<su2double> TargetArea(nVertex_NearField);
    vector<su2double> NearFieldWeight(nVertex_NearField);
    vector<su2double> Weight(nVertex_NearField);

    nVertex_NearField = 0;
    for (int iProcessor = 0; iProcessor < size; iProcessor++) {
      for (unsigned long iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
        Xcoord[nVertex_NearField] = Buffer_Receive_Xcoord[iProcessor*MaxLocalVertex_NearField+iVertex];
        Ycoord[nVertex_NearField] = Buffer_Receive_Ycoord[iProcessor*MaxLocalVertex_NearField+iVertex];

        if (nDim == 2) {
          AzimuthalAngle[nVertex_NearField] = 0;
        }

        if (nDim == 3) {
          Zcoord[nVertex_NearField] = Buffer_Receive_Zcoord[iProcessor*MaxLocalVertex_NearField+iVertex];

          /*--- Rotate the nearfield cylinder  ---*/

          su2double YcoordRot = Ycoord[nVertex_NearField];
          su2double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA);

          /*--- Compute the Azimuthal angle ---*/

          su2double AngleDouble = fabs(atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER);

          /*--- Fix an azimuthal line due to misalignments of the near-field ---*/

          su2double FixAzimuthalLine = config->GetFixAzimuthalLine();

          if ((AngleDouble >= FixAzimuthalLine - 0.1) && (AngleDouble <= FixAzimuthalLine + 0.1))
            AngleDouble = FixAzimuthalLine - 0.1;

          const auto AngleInt = SU2_TYPE::Short(floor(AngleDouble + 0.5));

          if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
          else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
        }

        if (AzimuthalAngle[nVertex_NearField] <= 60) {
          IdPoint[nVertex_NearField] = Buffer_Receive_IdPoint[iProcessor*MaxLocalVertex_NearField+iVertex];
          Pressure[nVertex_NearField] = Buffer_Receive_Pressure[iProcessor*MaxLocalVertex_NearField+iVertex];
          FaceArea[nVertex_NearField] = Buffer_Receive_FaceArea[iProcessor*MaxLocalVertex_NearField+iVertex];
          IdDomain[nVertex_NearField] = iProcessor;
          nVertex_NearField++;
        }
      }
    }


    vector<short> PhiAngleList;
    vector<short>::iterator IterPhiAngleList;

    for (unsigned long iVertex = 0; iVertex < nVertex_NearField; iVertex++)
      PhiAngleList.push_back(AzimuthalAngle[iVertex]);

    sort( PhiAngleList.begin(), PhiAngleList.end());
    IterPhiAngleList = unique( PhiAngleList.begin(), PhiAngleList.end());
    PhiAngleList.resize( IterPhiAngleList - PhiAngleList.begin() );

    /*--- Create vectors and distribute the values among the different PhiAngle queues ---*/

    vector<vector<su2double> > Xcoord_PhiAngle(PhiAngleList.size());
    vector<vector<su2double> > Ycoord_PhiAngle(PhiAngleList.size());
    vector<vector<su2double> > Zcoord_PhiAngle(PhiAngleList.size());
    vector<vector<unsigned long> > IdPoint_PhiAngle(PhiAngleList.size());
    vector<vector<unsigned long> > IdDomain_PhiAngle(PhiAngleList.size());
    vector<vector<su2double> > Pressure_PhiAngle(PhiAngleList.size());
    vector<vector<su2double> > FaceArea_PhiAngle(PhiAngleList.size());
    vector<vector<su2double> > EquivArea_PhiAngle(PhiAngleList.size());
    vector<vector<su2double> > TargetArea_PhiAngle(PhiAngleList.size());
    vector<vector<su2double> > NearFieldWeight_PhiAngle(PhiAngleList.size());
    vector<vector<su2double> > Weight_PhiAngle(PhiAngleList.size());

    /*--- Distribute the values among the different PhiAngles ---*/

    for (unsigned long iVertex = 0; iVertex < nVertex_NearField; iVertex++)
      for (unsigned short iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        if (AzimuthalAngle[iVertex] == PhiAngleList[iPhiAngle]) {
          Xcoord_PhiAngle[iPhiAngle].push_back(Xcoord[iVertex]);
          Ycoord_PhiAngle[iPhiAngle].push_back(Ycoord[iVertex]);
          Zcoord_PhiAngle[iPhiAngle].push_back(Zcoord[iVertex]);
          IdPoint_PhiAngle[iPhiAngle].push_back(IdPoint[iVertex]);
          IdDomain_PhiAngle[iPhiAngle].push_back(IdDomain[iVertex]);
          Pressure_PhiAngle[iPhiAngle].push_back(Pressure[iVertex]);
          FaceArea_PhiAngle[iPhiAngle].push_back(FaceArea[iVertex]);
          EquivArea_PhiAngle[iPhiAngle].push_back(EquivArea[iVertex]);
          TargetArea_PhiAngle[iPhiAngle].push_back(TargetArea[iVertex]);
          NearFieldWeight_PhiAngle[iPhiAngle].push_back(NearFieldWeight[iVertex]);
          Weight_PhiAngle[iPhiAngle].push_back(Weight[iVertex]);
        }

    /*--- Order the arrays (x Coordinate, Pressure, Point, and Domain) ---*/

    for (unsigned long iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (unsigned long iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++)
        for (unsigned long jVertex = 0; jVertex < Xcoord_PhiAngle[iPhiAngle].size() - 1 - iVertex; jVertex++)
          if (Xcoord_PhiAngle[iPhiAngle][jVertex] > Xcoord_PhiAngle[iPhiAngle][jVertex+1]) {
            auxXCoord = Xcoord_PhiAngle[iPhiAngle][jVertex]; Xcoord_PhiAngle[iPhiAngle][jVertex] = Xcoord_PhiAngle[iPhiAngle][jVertex+1]; Xcoord_PhiAngle[iPhiAngle][jVertex+1] = auxXCoord;
            auxYCoord = Ycoord_PhiAngle[iPhiAngle][jVertex]; Ycoord_PhiAngle[iPhiAngle][jVertex] = Ycoord_PhiAngle[iPhiAngle][jVertex+1]; Ycoord_PhiAngle[iPhiAngle][jVertex+1] = auxYCoord;
            auxZCoord = Zcoord_PhiAngle[iPhiAngle][jVertex]; Zcoord_PhiAngle[iPhiAngle][jVertex] = Zcoord_PhiAngle[iPhiAngle][jVertex+1]; Zcoord_PhiAngle[iPhiAngle][jVertex+1] = auxZCoord;
            auxPress = Pressure_PhiAngle[iPhiAngle][jVertex]; Pressure_PhiAngle[iPhiAngle][jVertex] = Pressure_PhiAngle[iPhiAngle][jVertex+1]; Pressure_PhiAngle[iPhiAngle][jVertex+1] = auxPress;
            auxArea = FaceArea_PhiAngle[iPhiAngle][jVertex]; FaceArea_PhiAngle[iPhiAngle][jVertex] = FaceArea_PhiAngle[iPhiAngle][jVertex+1]; FaceArea_PhiAngle[iPhiAngle][jVertex+1] = auxArea;
            auxPoint = IdPoint_PhiAngle[iPhiAngle][jVertex]; IdPoint_PhiAngle[iPhiAngle][jVertex] = IdPoint_PhiAngle[iPhiAngle][jVertex+1]; IdPoint_PhiAngle[iPhiAngle][jVertex+1] = auxPoint;
            auxDomain = IdDomain_PhiAngle[iPhiAngle][jVertex]; IdDomain_PhiAngle[iPhiAngle][jVertex] = IdDomain_PhiAngle[iPhiAngle][jVertex+1]; IdDomain_PhiAngle[iPhiAngle][jVertex+1] = auxDomain;
          }


    /*--- Check that all the azimuth lists have the same size ---*/

    auto nVertex = Xcoord_PhiAngle[0].size();
    for (unsigned long iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      auto nVertex_aux = Xcoord_PhiAngle[iPhiAngle].size();
      if (nVertex_aux != nVertex) cout <<"Be careful! One azimuth list is shorter than the other.\n";
      nVertex = min(nVertex, nVertex_aux);
    }

    /*--- Compute equivalent area distribution at each azimuth angle ---*/

    for (unsigned long iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      EquivArea_PhiAngle[iPhiAngle][0] = 0.0;
      for (unsigned long iVertex = 1; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        EquivArea_PhiAngle[iPhiAngle][iVertex] = 0.0;

        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][iVertex]*sin(AoA);

        for (unsigned long jVertex = 0; jVertex < iVertex-1; jVertex++) {

          Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex]*sin(AoA);
          jp1Coord = Xcoord_PhiAngle[iPhiAngle][jVertex+1]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex+1]*sin(AoA);

          jFunction = factor*(Pressure_PhiAngle[iPhiAngle][jVertex] - Pressure_Inf)*sqrt(Coord_i-Coord_j);
          jp1Function = factor*(Pressure_PhiAngle[iPhiAngle][jVertex+1] - Pressure_Inf)*sqrt(Coord_i-jp1Coord);

          DeltaX = (jp1Coord-Coord_j);
          MeanFunction = 0.5*(jp1Function + jFunction);
          EquivArea_PhiAngle[iPhiAngle][iVertex] += DeltaX * MeanFunction;
        }
      }
    }

    /*--- Create a file with the equivalent area distribution at each azimuthal angle ---*/

    NearFieldEA_file.precision(config->GetOutput_Precision());

    NearFieldEA_file.open("Equivalent_Area.dat", ios::out);
    NearFieldEA_file << "TITLE = \"Equivalent Area evaluation at each azimuthal angle\"" << "\n";

    if (config->GetSystemMeasurements() == US)
      NearFieldEA_file << "VARIABLES = \"Height (in) at r="<< R_Plane*12.0 << " in. (cyl. coord. system)\"";
    else
      NearFieldEA_file << "VARIABLES = \"Height (m) at r="<< R_Plane << " m. (cylindrical coordinate system)\"";

    for (unsigned long iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      if (config->GetSystemMeasurements() == US)
        NearFieldEA_file << ", \"Equivalent Area (ft<sup>2</sup>), <greek>F</greek>= " << PhiAngleList[iPhiAngle] << " deg.\"";
      else
        NearFieldEA_file << ", \"Equivalent Area (m<sup>2</sup>), <greek>F</greek>= " << PhiAngleList[iPhiAngle] << " deg.\"";
    }

    NearFieldEA_file << "\n";
    for (unsigned long iVertex = 0; iVertex < EquivArea_PhiAngle[0].size(); iVertex++) {

      su2double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
      su2double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);

      if (config->GetSystemMeasurements() == US)
        NearFieldEA_file << scientific << (XcoordRot - XcoordRot_init) * 12.0;
      else
        NearFieldEA_file << scientific << (XcoordRot - XcoordRot_init);

      for (unsigned long iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
        NearFieldEA_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex];
      }

      NearFieldEA_file << "\n";

    }
    NearFieldEA_file.close();


    /*--- Read target equivalent area from the configuration file,
     this first implementation requires a complete table (same as the original
     EA table). so... no interpolation. ---*/

    vector<vector<su2double> > TargetArea_PhiAngle_Trans;
    TargetEA_file.open("TargetEA.dat", ios::in);

    if (TargetEA_file.fail()) {
      /*--- Set the table to 0 ---*/
      for (unsigned long iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        for (unsigned long iVertex = 0; iVertex < TargetArea_PhiAngle[iPhiAngle].size(); iVertex++)
          TargetArea_PhiAngle[iPhiAngle][iVertex] = 0.0;
    }
    else {

      /*--- skip header lines ---*/

      string line;
      getline(TargetEA_file, line);
      getline(TargetEA_file, line);

      while (TargetEA_file) {

        string line;
        getline(TargetEA_file, line);
        istringstream is(line);
        vector<su2double> row;
        unsigned short iter = 0;

        while (is.good()) {
          string token;
          getline(is, token,',');

          istringstream js(token);

          su2double data;
          js >> data;

          /*--- The first element in the table is the coordinate (in or m)---*/

          if (iter != 0) row.push_back(data);
          iter++;

        }
        TargetArea_PhiAngle_Trans.push_back(row);
      }

      for (unsigned long iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        for (unsigned long iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++)
          TargetArea_PhiAngle[iPhiAngle][iVertex] = TargetArea_PhiAngle_Trans[iVertex][iPhiAngle];

    }

    /*--- Divide by the number of Phi angles in the nearfield ---*/

    su2double PhiFactor = 1.0/su2double(PhiAngleList.size());

    /*--- Evaluate the objective function ---*/

    InverseDesign = 0;
    for (unsigned long iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (unsigned long iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;
        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];

        su2double Difference = EquivArea_PhiAngle[iPhiAngle][iVertex]-TargetArea_PhiAngle[iPhiAngle][iVertex];
        su2double percentage = fabs(Difference)*100/fabs(TargetArea_PhiAngle[iPhiAngle][iVertex]);

        if ((percentage < 0.1) || (Coord_i < XCoordBegin_OF) || (Coord_i > XCoordEnd_OF)) Difference = 0.0;

        InverseDesign += EAScaleFactor*PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*Difference*Difference;
      }

    /*--- Evaluate the weight of the nearfield pressure (adjoint input) ---*/

    for (unsigned long iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      for (unsigned long iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];
        NearFieldWeight_PhiAngle[iPhiAngle][iVertex] = 0.0;
        for (unsigned long jVertex = iVertex; jVertex < EquivArea_PhiAngle[iPhiAngle].size(); jVertex++) {
          Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex];
          Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;

          su2double Difference = EquivArea_PhiAngle[iPhiAngle][jVertex]-TargetArea_PhiAngle[iPhiAngle][jVertex];
          su2double percentage = fabs(Difference)*100/fabs(TargetArea_PhiAngle[iPhiAngle][jVertex]);

          if ((percentage < 0.1) || (Coord_j < XCoordBegin_OF) || (Coord_j > XCoordEnd_OF)) Difference = 0.0;

          NearFieldWeight_PhiAngle[iPhiAngle][iVertex] += EAScaleFactor*PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*2.0*Difference*factor*sqrt(Coord_j-Coord_i);
        }
      }
    }

    /*--- Write the Nearfield pressure at each Azimuthal PhiAngle ---*/

    EquivArea_file.precision(config->GetOutput_Precision());

    EquivArea_file.open("nearfield_flow.dat", ios::out);
    EquivArea_file << "TITLE = \"Equivalent Area evaluation at each azimuthal angle\"" << "\n";

    if (config->GetSystemMeasurements() == US)
      EquivArea_file << "VARIABLES = \"Height (in) at r="<< R_Plane*12.0 << " in. (cyl. coord. system)\",\"Equivalent Area (ft<sup>2</sup>)\",\"Target Equivalent Area (ft<sup>2</sup>)\",\"Cp\"" << "\n";
    else
      EquivArea_file << "VARIABLES = \"Height (m) at r="<< R_Plane << " m. (cylindrical coordinate system)\",\"Equivalent Area (m<sup>2</sup>)\",\"Target Equivalent Area (m<sup>2</sup>)\",\"Cp\"" << "\n";

    for (unsigned long iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      EquivArea_file << fixed << "ZONE T= \"<greek>F</greek>=" << PhiAngleList[iPhiAngle] << " deg.\"" << "\n";
      for (unsigned long iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++) {

        su2double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
        su2double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);

        if (config->GetSystemMeasurements() == US)
          EquivArea_file << scientific << (XcoordRot - XcoordRot_init) * 12.0;
        else
          EquivArea_file << scientific << (XcoordRot - XcoordRot_init);

        EquivArea_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex]
        << ", " << TargetArea_PhiAngle[iPhiAngle][iVertex] << ", " << (Pressure_PhiAngle[iPhiAngle][iVertex]-Pressure_Inf)/Pressure_Inf << "\n";
      }
    }

    EquivArea_file.close();

  }

  /*--- Send the value of the NearField coefficient to all the processors ---*/

  SU2_MPI::Bcast(&InverseDesign, 1, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());

  /*--- Store the value of the NearField coefficient ---*/

  solver->SetTotal_CEquivArea(InverseDesign);
  SetHistoryOutputValue("EQUIVALENT_AREA", InverseDesign);

}

void CFlowOutput::WriteAdditionalFiles(CConfig *config, CGeometry *geometry, CSolver **solver_container){

  if (config->GetFixed_CL_Mode() ||
      (config->GetKind_Streamwise_Periodic() == ENUM_STREAMWISE_PERIODIC::MASSFLOW)){
    WriteMetaData(config);
  }

  if (config->GetWrt_ForcesBreakdown()){
    WriteForcesBreakdown(config, solver_container[FLOW_SOL]);
  }

}

void CFlowOutput::WriteMetaData(const CConfig *config){

  ofstream meta_file;

  string filename = "flow";

  filename = config->GetFilename(filename, ".meta", curTimeIter);

  /*--- All processors open the file. ---*/

  if (rank == MASTER_NODE) {
    cout << "Writing Flow Meta-Data file: " << filename << endl;

    meta_file.open(filename.c_str(), ios::out);
    meta_file.precision(15);

    if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST || config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)
      meta_file <<"ITER= " << curTimeIter + 1 << endl;
    else
      meta_file <<"ITER= " << curInnerIter + config->GetExtIter_OffSet() + 1 << endl;

    if (config->GetFixed_CL_Mode()){
      meta_file <<"AOA= " << config->GetAoA() - config->GetAoA_Offset() << endl;
      meta_file <<"SIDESLIP_ANGLE= " << config->GetAoS() - config->GetAoS_Offset() << endl;
      meta_file <<"DCD_DCL_VALUE= " << config->GetdCD_dCL() << endl;
      if (nDim==3){
        meta_file <<"DCMX_DCL_VALUE= " << config->GetdCMx_dCL() << endl;
        meta_file <<"DCMY_DCL_VALUE= " << config->GetdCMy_dCL() << endl;
      }
      meta_file <<"DCMZ_DCL_VALUE= " << config->GetdCMz_dCL() << endl;
    }
    meta_file <<"INITIAL_BCTHRUST= " << config->GetInitial_BCThrust() << endl;


    if (( config->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_EULER ||
          config->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES ||
          config->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_RANS )) {
      meta_file << "SENS_AOA=" << GetHistoryFieldValue("SENS_AOA") * PI_NUMBER / 180.0 << endl;
    }

    if(config->GetKind_Streamwise_Periodic() == ENUM_STREAMWISE_PERIODIC::MASSFLOW) {
      meta_file << "STREAMWISE_PERIODIC_PRESSURE_DROP=" << GetHistoryFieldValue("STREAMWISE_DP") << endl;
    }
  }

  meta_file.close();
}

void CFlowOutput::WriteForcesBreakdown(const CConfig* config, const CSolver* flow_solver) const {
  // clang-format off
  if (rank != MASTER_NODE) return;

  const bool compressible = (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);
  const bool incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
  const bool unsteady = config->GetTime_Domain();
  const bool viscous = config->GetViscous();
  const bool dynamic_grid = config->GetDynamic_Grid();
  const bool gravity = config->GetGravityForce();
  const TURB_MODEL Kind_Turb_Model = config->GetKind_Turb_Model();
  const bool turbulent = Kind_Turb_Model != TURB_MODEL::NONE;
  const TURB_TRANS_MODEL Kind_Trans_Model = config->GetKind_Trans_Model();
  const bool transition = Kind_Trans_Model != TURB_TRANS_MODEL::NONE;
  const bool fixed_cl = config->GetFixed_CL_Mode();
  const auto Kind_Solver = config->GetKind_Solver();
  const auto Ref_NonDim = config->GetRef_NonDim();
  const auto nMonitoring = config->GetnMarker_Monitoring();

  auto fileName = config->GetBreakdown_FileName();
  if (unsteady) {
    const auto lastindex = fileName.find_last_of('.');
    const auto ext = fileName.substr(lastindex, fileName.size());
    fileName = fileName.substr(0, lastindex);
    fileName = config->GetFilename(fileName, ext, curTimeIter);
  }

  /*--- Output the mean flow solution using only the master node ---*/

  cout << "\nWriting the forces breakdown file (" << fileName << ")." << endl;

  vector<su2double> Surface_CL(nMonitoring);
  vector<su2double> Surface_CD(nMonitoring);
  vector<su2double> Surface_CSF(nMonitoring);
  vector<su2double> Surface_CEff(nMonitoring);
  vector<su2double> Surface_CFx(nMonitoring);
  vector<su2double> Surface_CFy(nMonitoring);
  vector<su2double> Surface_CFz(nMonitoring);
  vector<su2double> Surface_CMx(nMonitoring);
  vector<su2double> Surface_CMy(nMonitoring);
  vector<su2double> Surface_CMz(nMonitoring);

  vector<su2double> Surface_CL_Inv(nMonitoring);
  vector<su2double> Surface_CD_Inv(nMonitoring);
  vector<su2double> Surface_CSF_Inv(nMonitoring);
  vector<su2double> Surface_CEff_Inv(nMonitoring);
  vector<su2double> Surface_CFx_Inv(nMonitoring);
  vector<su2double> Surface_CFy_Inv(nMonitoring);
  vector<su2double> Surface_CFz_Inv(nMonitoring);
  vector<su2double> Surface_CMx_Inv(nMonitoring);
  vector<su2double> Surface_CMy_Inv(nMonitoring);
  vector<su2double> Surface_CMz_Inv(nMonitoring);

  vector<su2double> Surface_CL_Visc(nMonitoring);
  vector<su2double> Surface_CD_Visc(nMonitoring);
  vector<su2double> Surface_CSF_Visc(nMonitoring);
  vector<su2double> Surface_CEff_Visc(nMonitoring);
  vector<su2double> Surface_CFx_Visc(nMonitoring);
  vector<su2double> Surface_CFy_Visc(nMonitoring);
  vector<su2double> Surface_CFz_Visc(nMonitoring);
  vector<su2double> Surface_CMx_Visc(nMonitoring);
  vector<su2double> Surface_CMy_Visc(nMonitoring);
  vector<su2double> Surface_CMz_Visc(nMonitoring);

  vector<su2double> Surface_CL_Mnt(nMonitoring);
  vector<su2double> Surface_CD_Mnt(nMonitoring);
  vector<su2double> Surface_CSF_Mnt(nMonitoring);
  vector<su2double> Surface_CEff_Mnt(nMonitoring);
  vector<su2double> Surface_CFx_Mnt(nMonitoring);
  vector<su2double> Surface_CFy_Mnt(nMonitoring);
  vector<su2double> Surface_CFz_Mnt(nMonitoring);
  vector<su2double> Surface_CMx_Mnt(nMonitoring);
  vector<su2double> Surface_CMy_Mnt(nMonitoring);
  vector<su2double> Surface_CMz_Mnt(nMonitoring);

  /*--- Flow solution coefficients ---*/

  const auto Total_CL = flow_solver->GetTotal_CL();
  const auto Total_CD = flow_solver->GetTotal_CD();
  const auto Total_CSF = flow_solver->GetTotal_CSF();
  const auto Total_CEff = flow_solver->GetTotal_CEff();
  const auto Total_CMx = flow_solver->GetTotal_CMx();
  const auto Total_CMy = flow_solver->GetTotal_CMy();
  const auto Total_CMz = flow_solver->GetTotal_CMz();
  const auto Total_CFx = flow_solver->GetTotal_CFx();
  const auto Total_CFy = flow_solver->GetTotal_CFy();
  const auto Total_CFz = flow_solver->GetTotal_CFz();

  su2double Total_CoPx = 0.0, Total_CoPy = 0.0, Total_CoPz = 0.0;
  if (nDim == 2) {
    Total_CoPx = flow_solver->GetTotal_CoPx() / flow_solver->GetTotal_CFy();
    Total_CoPy = flow_solver->GetTotal_CoPy() / flow_solver->GetTotal_CFx();
  } else {
    Total_CoPx = flow_solver->GetTotal_CoPx() / flow_solver->GetTotal_CFz();
    Total_CoPz = flow_solver->GetTotal_CoPz() / flow_solver->GetTotal_CFx();
  }
  if (us_units) {
    Total_CoPx *= 12.0;
    Total_CoPy *= 12.0;
    Total_CoPz *= 12.0;
  }

  /*--- Flow inviscid solution coefficients ---*/

  const auto Inv_CL = flow_solver->GetAllBound_CL_Inv();
  const auto Inv_CD = flow_solver->GetAllBound_CD_Inv();
  const auto Inv_CSF = flow_solver->GetAllBound_CSF_Inv();
  const auto Inv_CEff = flow_solver->GetAllBound_CEff_Inv();
  const auto Inv_CMx = flow_solver->GetAllBound_CMx_Inv();
  const auto Inv_CMy = flow_solver->GetAllBound_CMy_Inv();
  const auto Inv_CMz = flow_solver->GetAllBound_CMz_Inv();
  const auto Inv_CFx = flow_solver->GetAllBound_CFx_Inv();
  const auto Inv_CFy = flow_solver->GetAllBound_CFy_Inv();
  const auto Inv_CFz = flow_solver->GetAllBound_CFz_Inv();

  /*--- Flow viscous solution coefficients ---*/

  const auto Visc_CL = flow_solver->GetAllBound_CL_Visc();
  const auto Visc_CD = flow_solver->GetAllBound_CD_Visc();
  const auto Visc_CSF = flow_solver->GetAllBound_CSF_Visc();
  const auto Visc_CEff = flow_solver->GetAllBound_CEff_Visc();
  const auto Visc_CMx = flow_solver->GetAllBound_CMx_Visc();
  const auto Visc_CMy = flow_solver->GetAllBound_CMy_Visc();
  const auto Visc_CMz = flow_solver->GetAllBound_CMz_Visc();
  const auto Visc_CFx = flow_solver->GetAllBound_CFx_Visc();
  const auto Visc_CFy = flow_solver->GetAllBound_CFy_Visc();
  const auto Visc_CFz = flow_solver->GetAllBound_CFz_Visc();

  /*--- Flow momentum solution coefficients ---*/

  const auto Mnt_CL = flow_solver->GetAllBound_CL_Mnt();
  const auto Mnt_CD = flow_solver->GetAllBound_CD_Mnt();
  const auto Mnt_CSF = flow_solver->GetAllBound_CSF_Mnt();
  const auto Mnt_CEff = flow_solver->GetAllBound_CEff_Mnt();
  const auto Mnt_CMx = flow_solver->GetAllBound_CMx_Mnt();
  const auto Mnt_CMy = flow_solver->GetAllBound_CMy_Mnt();
  const auto Mnt_CMz = flow_solver->GetAllBound_CMz_Mnt();
  const auto Mnt_CFx = flow_solver->GetAllBound_CFx_Mnt();
  const auto Mnt_CFy = flow_solver->GetAllBound_CFy_Mnt();
  const auto Mnt_CFz = flow_solver->GetAllBound_CFz_Mnt();

  /*--- Look over the markers being monitored and get the desired values ---*/

  for (auto iMarker = 0u; iMarker < nMonitoring; iMarker++) {
    Surface_CL[iMarker] = flow_solver->GetSurface_CL(iMarker);
    Surface_CD[iMarker] = flow_solver->GetSurface_CD(iMarker);
    Surface_CSF[iMarker] = flow_solver->GetSurface_CSF(iMarker);
    Surface_CEff[iMarker] = flow_solver->GetSurface_CEff(iMarker);
    Surface_CMx[iMarker] = flow_solver->GetSurface_CMx(iMarker);
    Surface_CMy[iMarker] = flow_solver->GetSurface_CMy(iMarker);
    Surface_CMz[iMarker] = flow_solver->GetSurface_CMz(iMarker);
    Surface_CFx[iMarker] = flow_solver->GetSurface_CFx(iMarker);
    Surface_CFy[iMarker] = flow_solver->GetSurface_CFy(iMarker);
    Surface_CFz[iMarker] = flow_solver->GetSurface_CFz(iMarker);

    Surface_CL_Inv[iMarker] = flow_solver->GetSurface_CL_Inv(iMarker);
    Surface_CD_Inv[iMarker] = flow_solver->GetSurface_CD_Inv(iMarker);
    Surface_CSF_Inv[iMarker] = flow_solver->GetSurface_CSF_Inv(iMarker);
    Surface_CEff_Inv[iMarker] = flow_solver->GetSurface_CEff_Inv(iMarker);
    Surface_CMx_Inv[iMarker] = flow_solver->GetSurface_CMx_Inv(iMarker);
    Surface_CMy_Inv[iMarker] = flow_solver->GetSurface_CMy_Inv(iMarker);
    Surface_CMz_Inv[iMarker] = flow_solver->GetSurface_CMz_Inv(iMarker);
    Surface_CFx_Inv[iMarker] = flow_solver->GetSurface_CFx_Inv(iMarker);
    Surface_CFy_Inv[iMarker] = flow_solver->GetSurface_CFy_Inv(iMarker);
    Surface_CFz_Inv[iMarker] = flow_solver->GetSurface_CFz_Inv(iMarker);
    Surface_CL_Visc[iMarker] = flow_solver->GetSurface_CL_Visc(iMarker);
    Surface_CD_Visc[iMarker] = flow_solver->GetSurface_CD_Visc(iMarker);
    Surface_CSF_Visc[iMarker] = flow_solver->GetSurface_CSF_Visc(iMarker);
    Surface_CEff_Visc[iMarker] = flow_solver->GetSurface_CEff_Visc(iMarker);
    Surface_CMx_Visc[iMarker] = flow_solver->GetSurface_CMx_Visc(iMarker);
    Surface_CMy_Visc[iMarker] = flow_solver->GetSurface_CMy_Visc(iMarker);
    Surface_CMz_Visc[iMarker] = flow_solver->GetSurface_CMz_Visc(iMarker);
    Surface_CFx_Visc[iMarker] = flow_solver->GetSurface_CFx_Visc(iMarker);
    Surface_CFy_Visc[iMarker] = flow_solver->GetSurface_CFy_Visc(iMarker);
    Surface_CFz_Visc[iMarker] = flow_solver->GetSurface_CFz_Visc(iMarker);

    Surface_CL_Mnt[iMarker] = flow_solver->GetSurface_CL_Mnt(iMarker);
    Surface_CD_Mnt[iMarker] = flow_solver->GetSurface_CD_Mnt(iMarker);
    Surface_CSF_Mnt[iMarker] = flow_solver->GetSurface_CSF_Mnt(iMarker);
    Surface_CEff_Mnt[iMarker] = flow_solver->GetSurface_CEff_Mnt(iMarker);
    Surface_CMx_Mnt[iMarker] = flow_solver->GetSurface_CMx_Mnt(iMarker);
    Surface_CMy_Mnt[iMarker] = flow_solver->GetSurface_CMy_Mnt(iMarker);
    Surface_CMz_Mnt[iMarker] = flow_solver->GetSurface_CMz_Mnt(iMarker);
    Surface_CFx_Mnt[iMarker] = flow_solver->GetSurface_CFx_Mnt(iMarker);
    Surface_CFy_Mnt[iMarker] = flow_solver->GetSurface_CFy_Mnt(iMarker);
    Surface_CFz_Mnt[iMarker] = flow_solver->GetSurface_CFz_Mnt(iMarker);
  }

  /*--- Write file name with extension ---*/

  ofstream file;
  file.open(fileName);

  file << "\n";
  file << "-------------------------------------------------------------------------\n";
  file << "|    ___ _   _ ___                                                      |\n";
  file << "|   / __| | | |_  )   Release 8.0.0 \"Harrier\"                           |\n";
  file << "|   \\__ \\ |_| |/ /                                                      |\n";
  file << "|   |___/\\___//___|   Suite (Computational Fluid Dynamics Code)         |\n";
  file << "|                                                                       |\n";
  file << "-------------------------------------------------------------------------\n";
  file << "| SU2 Project Website: https://su2code.github.io                        |\n";
  file << "|                                                                       |\n";
  file << "| The SU2 Project is maintained by the SU2 Foundation                   |\n";
  file << "| (http://su2foundation.org)                                            |\n";
  file << "-------------------------------------------------------------------------\n";
  file << "| Copyright 2012-2023, SU2 Contributors                                 |\n";
  file << "|                                                                       |\n";
  file << "| SU2 is free software; you can redistribute it and/or                  |\n";
  file << "| modify it under the terms of the GNU Lesser General Public            |\n";
  file << "| License as published by the Free Software Foundation; either          |\n";
  file << "| version 2.1 of the License, or (at your option) any later version.    |\n";
  file << "|                                                                       |\n";
  file << "| SU2 is distributed in the hope that it will be useful,                |\n";
  file << "| but WITHOUT ANY WARRANTY; without even the implied warranty of        |\n";
  file << "| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |\n";
  file << "| Lesser General Public License for more details.                       |\n";
  file << "|                                                                       |\n";
  file << "| You should have received a copy of the GNU Lesser General Public      |\n";
  file << "| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |\n";
  file << "-------------------------------------------------------------------------\n";

  file.precision(6);

  file << "\n\nProblem definition:\n\n";

  switch (Kind_Solver) {
    case MAIN_SOLVER::EULER:
    case MAIN_SOLVER::INC_EULER:
      if (compressible) file << "Compressible Euler equations.\n";
      if (incompressible) file << "Incompressible Euler equations.\n";
      break;
    case MAIN_SOLVER::NAVIER_STOKES:
    case MAIN_SOLVER::INC_NAVIER_STOKES:
      if (compressible) file << "Compressible Laminar Navier-Stokes' equations.\n";
      if (incompressible) file << "Incompressible Laminar Navier-Stokes' equations.\n";
      break;
    case MAIN_SOLVER::RANS:
    case MAIN_SOLVER::INC_RANS:
      if (compressible) file << "Compressible RANS equations.\n";
      if (incompressible) file << "Incompressible RANS equations.\n";
      file << "Turbulence model: ";
      switch (Kind_Turb_Model) {
        case TURB_MODEL::NONE: break;
        case TURB_MODEL::SA:
          /// TODO: add the submodels here
          file << "Spalart Allmaras\n";
          break;
        case TURB_MODEL::SST:
          /// TODO: add the submodels here
          if (config->GetSSTParsedOptions().sust)
            file << "Menter's SST with sustaining terms\n";
          else
            file << "Menter's SST\n";
         break;
      }
      if (transition) {
        file << "Transition model: ";
        switch (Kind_Trans_Model) {
        case TURB_TRANS_MODEL::NONE: break;
        case TURB_TRANS_MODEL::LM:
          file << "Langtry and Menter's transition";
          if (config->GetLMParsedOptions().LM2015) {
            file << " w/ cross-flow corrections (2015)\n";
          } else {
            file << " (2009)\n";
          }
          break;
        }
      }
      break;
    default:
      break;
  }

  /*--- Compressible version of console output ---*/

  if (compressible) {
    file << "Mach number: " << config->GetMach() << ".\n";
    file << "Angle of attack (AoA): " << config->GetAoA() << " deg, and angle of sideslip (AoS): " << config->GetAoS()
         << " deg.\n";
    if (viscous)
      file << "Reynolds number: " << config->GetReynolds() << ".\n";

    if (fixed_cl) {
      file << "Simulation at a cte. CL: " << config->GetTarget_CL() << ".\n";
      file << "Approx. Delta CL / Delta AoA: " << config->GetdCL_dAlpha() << " (1/deg).\n";
      file << "Approx. Delta CD / Delta CL: " << config->GetdCD_dCL() << ".\n";
      if (nDim == 3) {
        file << "Approx. Delta CMx / Delta CL: " << config->GetdCMx_dCL() << ".\n";
        file << "Approx. Delta CMy / Delta CL: " << config->GetdCMy_dCL() << ".\n";
      }
      file << "Approx. Delta CMz / Delta CL: " << config->GetdCMz_dCL() << ".\n";
    }

    if (Ref_NonDim == DIMENSIONAL) {
      file << "Dimensional simulation.\n";
    } else if (Ref_NonDim == FREESTREAM_PRESS_EQ_ONE) {
      file << "Non-Dimensional simulation (P=1.0, Rho=1.0, T=1.0 at the farfield).\n";
    } else if (Ref_NonDim == FREESTREAM_VEL_EQ_MACH) {
      file << "Non-Dimensional simulation (V=Mach, Rho=1.0, T=1.0 at the farfield).\n";
    } else if (Ref_NonDim == FREESTREAM_VEL_EQ_ONE) {
      file << "Non-Dimensional simulation (V=1.0, Rho=1.0, T=1.0 at the farfield).\n";
    }

    if (si_units) {
      file << "The reference area is " << config->GetRefArea() << " m^2.\n";
      file << "The reference length is " << config->GetRefLength() << " m.\n";
    }

    if (us_units) {
      file << "The reference area is " << config->GetRefArea() * 12.0 * 12.0 << " in^2.\n";
      file << "The reference length is " << config->GetRefLength() * 12.0 << " in.\n";
    }
    file << "\n\nProblem definition:\n\n";

    if (viscous) {
      file << "Viscous flow: Computing pressure using the ideal gas law\n";
      file << "based on the free-stream temperature and a density computed\n";
      file << "from the Reynolds number.\n";
    } else {
      file << "Inviscid flow: Computing density based on free-stream\n";
      file << "temperature and pressure using the ideal gas law.\n";
    }

    if (dynamic_grid)
      file << "Force coefficients computed using MACH_MOTION.\n";
    else
      file << "Force coefficients computed using free-stream values.\n";

    file << "-- Input conditions:\n";

    switch (config->GetKind_FluidModel()) {
      case STANDARD_AIR:
        file << "Fluid Model: STANDARD_AIR \n";
        file << "Specific gas constant: " << config->GetGas_Constant();
        if (si_units) file << " N.m/kg.K.\n";
        else file << " lbf.ft/slug.R.\n";
        file << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << "\n";
        file << "Specific Heat Ratio: 1.4000 \n";
        break;

      case IDEAL_GAS:
        file << "Fluid Model: IDEAL_GAS \n";
        file << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K.\n";
        file << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << "\n";
        file << "Specific Heat Ratio: " << config->GetGamma() << "\n";
        break;

      case VW_GAS:
        file << "Fluid Model: Van der Waals \n";
        file << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K.\n";
        file << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << "\n";
        file << "Specific Heat Ratio: " << config->GetGamma() << "\n";
        file << "Critical Pressure:   " << config->GetPressure_Critical() << " Pa.\n";
        file << "Critical Temperature:  " << config->GetTemperature_Critical() << " K.\n";
        file << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() / config->GetPressure_Ref()
             << "\n";
        file << "Critical Temperature (non-dim) :  "
             << config->GetTemperature_Critical() / config->GetTemperature_Ref() << "\n";
        break;

      case PR_GAS:
        file << "Fluid Model: Peng-Robinson \n";
        file << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K.\n";
        file << "Specific gas constant(non-dim): " << config->GetGas_ConstantND() << "\n";
        file << "Specific Heat Ratio: " << config->GetGamma() << "\n";
        file << "Critical Pressure:   " << config->GetPressure_Critical() << " Pa.\n";
        file << "Critical Temperature:  " << config->GetTemperature_Critical() << " K.\n";
        file << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() / config->GetPressure_Ref()
             << "\n";
        file << "Critical Temperature (non-dim) :  "
             << config->GetTemperature_Critical() / config->GetTemperature_Ref() << "\n";
        break;

     case FLUID_FLAMELET:
        file << "Fluid Model: FLAMELET \n";
        break;

      case COOLPROP: {
        CCoolProp auxFluidModel(config->GetFluid_Name());
        file << "Fluid Model: CoolProp library \n";
        file << "Specific gas constant: " << auxFluidModel.GetGas_Constant()<< " N.m/kg.K.\n";
        file << "Specific gas constant(non-dim): " << config->GetGas_ConstantND() << "\n";
        file << "Specific Heat Ratio: "<< auxFluidModel.GetGamma() << "\n";
        file << "Critical Pressure:   " << auxFluidModel.GetPressure_Critical() << " Pa.\n";
        file << "Critical Temperature:  " << auxFluidModel.GetTemperature_Critical()<< " K.\n";
        file << "Critical Pressure (non-dim):   " << auxFluidModel.GetPressure_Critical()/ config->GetPressure_Ref()
            << "\n";
        file << "Critical Temperature (non-dim) :  "
            << auxFluidModel.GetTemperature_Critical() / config->GetTemperature_Ref() << "\n";
        } break;
    }

    if (viscous) {
      switch (config->GetKind_ViscosityModel()) {
        case VISCOSITYMODEL::CONSTANT:
          file << "Viscosity Model: CONSTANT_VISCOSITY  \n";
          file << "Laminar Viscosity: " << config->GetMu_Constant();
          if (si_units) file << " N.s/m^2.\n";
          else file << " lbf.s/ft^2.\n";
          file << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND() << "\n";
          break;

        case VISCOSITYMODEL::COOLPROP:
          file << "Viscosity Model: CoolProp  \n";
          break;

        case VISCOSITYMODEL::SUTHERLAND:
          file << "Viscosity Model: SUTHERLAND \n";
          file << "Ref. Laminar Viscosity: " << config->GetMu_Ref();
          if (si_units) file << " N.s/m^2.\n";
          else file << " lbf.s/ft^2.\n";
          file << "Ref. Temperature: " << config->GetMu_Temperature_Ref();
          if (si_units) file << " K.\n";
          else file << " R.\n";
          file << "Sutherland Constant: " << config->GetMu_S();
          if (si_units) file << " K.\n";
          else file << " R.\n";
          file << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND() << "\n";
          file << "Ref. Temperature (non-dim): " << config->GetMu_Temperature_RefND() << "\n";
          file << "Sutherland constant (non-dim): " << config->GetMu_SND() << "\n";
          break;

        default:
          break;
      }
      switch (config->GetKind_ConductivityModel()) {
        case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
          file << "Conductivity Model: CONSTANT_PRANDTL \n";
          file << "Prandtl: " << config->GetPrandtl_Lam() << "\n";
          break;

        case CONDUCTIVITYMODEL::CONSTANT:
          file << "Conductivity Model: CONSTANT \n";
          file << "Molecular Conductivity: " << config->GetThermal_Conductivity_Constant() << " W/m^2.K.\n";
          file << "Molecular Conductivity (non-dim): " << config->GetThermal_Conductivity_ConstantND() << "\n";
          break;
        case CONDUCTIVITYMODEL::COOLPROP:
          file << "Conductivity Model: COOLPROP \n";
          break;
        default:
          break;
      }

      if (turbulent) {
        switch (config->GetKind_ConductivityModel_Turb()) {
          case CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL:
            file << "Turbulent Conductivity Model: CONSTANT_PRANDTL \n";
            file << "Turbulent Prandtl: " << config->GetPrandtl_Turb() << "\n";
            break;
          case CONDUCTIVITYMODEL_TURB::NONE:
            file << "Turbulent Conductivity Model: NONE \n";
            file << "No turbulent component in effective thermal conductivity.\n";
            break;
        }
      }
    }

    file << "Free-stream static pressure: " << config->GetPressure_FreeStream();
    if (si_units) file << " Pa.\n";
    else file << " psf.\n";

    file << "Free-stream total pressure: "
         << config->GetPressure_FreeStream() *
                pow(1.0 + config->GetMach() * config->GetMach() * 0.5 * (config->GetGamma() - 1.0),
                    config->GetGamma() / (config->GetGamma() - 1.0));
    if (si_units) file << " Pa.\n";
    else file << " psf.\n";

    file << "Free-stream temperature: " << config->GetTemperature_FreeStream();
    if (si_units) file << " K.\n";
    else file << " R.\n";

    file << "Free-stream total temperature: "
         << config->GetTemperature_FreeStream() *
                (1.0 + config->GetMach() * config->GetMach() * 0.5 * (config->GetGamma() - 1.0));
    if (si_units) file << " K.\n";
    else file << " R.\n";

    file << "Free-stream density: " << config->GetDensity_FreeStream();
    if (si_units) file << " kg/m^3.\n";
    else file << " slug/ft^3.\n";

    file << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0];
    file << ", " << config->GetVelocity_FreeStream()[1];
    if (nDim == 3) {
      file << ", " << config->GetVelocity_FreeStream()[2];
    }
    if (si_units) file << ") m/s. ";
    else file << ") ft/s. ";

    file << "Magnitude: " << config->GetModVel_FreeStream();
    if (si_units) file << " m/s.\n";
    else file << " ft/s.\n";

    file << "Free-stream total energy per unit mass: " << config->GetEnergy_FreeStream();
    if (si_units) file << " m^2/s^2.\n";
    else file << " ft^2/s^2.\n";

    if (viscous) {
      file << "Free-stream viscosity: " << config->GetViscosity_FreeStream();
      if (si_units) file << " N.s/m^2.\n";
      else file << " lbf.s/ft^2.\n";
      if (turbulent) {
        file << "Free-stream turb. kinetic energy per unit mass: " << config->GetTke_FreeStream();
        if (si_units) file << " m^2/s^2.\n";
        else file << " ft^2/s^2.\n";
        file << "Free-stream specific dissipation: " << config->GetOmega_FreeStream();
        if (si_units) file << " 1/s.\n";
        else file << " 1/s.\n";
      }
    }

    if (unsteady) {
      file << "Total time: " << config->GetTotal_UnstTime() << " s. Time step: " << config->GetDelta_UnstTime()
           << " s.\n";
    }

    /*--- Print out reference values. ---*/

    file << "-- Reference values:\n";

    file << "Reference specific gas constant: " << config->GetGas_Constant_Ref();
    if (si_units) file << " N.m/kg.K.\n";
    else file << " lbf.ft/slug.R.\n";

    file << "Reference pressure: " << config->GetPressure_Ref();
    if (si_units) file << " Pa.\n";
    else file << " psf.\n";

    file << "Reference temperature: " << config->GetTemperature_Ref();
    if (si_units) file << " K.\n";
    else file << " R.\n";

    file << "Reference density: " << config->GetDensity_Ref();
    if (si_units) file << " kg/m^3.\n";
    else file << " slug/ft^3.\n";

    file << "Reference velocity: " << config->GetVelocity_Ref();
    if (si_units) file << " m/s.\n";
    else file << " ft/s.\n";

    file << "Reference energy per unit mass: " << config->GetEnergy_Ref();
    if (si_units) file << " m^2/s^2.\n";
    else file << " ft^2/s^2.\n";

    if (viscous) {
      file << "Reference viscosity: " << config->GetViscosity_Ref();
      if (si_units) file << " N.s/m^2.\n";
      else file << " lbf.s/ft^2.\n";
      file << "Reference conductivity: " << config->GetThermal_Conductivity_Ref();
      if (si_units) file << " W/m^2.K.\n";
      else file << " lbf/ft.s.R.\n";
    }

    if (unsteady) file << "Reference time: " << config->GetTime_Ref() << " s.\n";

    /*--- Print out resulting non-dim values here. ---*/

    file << "-- Resulting non-dimensional state:\n";
    file << "Mach number (non-dim): " << config->GetMach() << "\n";
    if (viscous) {
      file << "Reynolds number (non-dim): " << config->GetReynolds() << ". Re length: " << config->GetLength_Reynolds();
      if (si_units) file << " m.\n";
      else file << " ft.\n";
    }
    if (gravity) {
      file << "Froude number (non-dim): " << config->GetFroude() << "\n";
      file << "Lenght of the baseline wave (non-dim): " << 2.0 * PI_NUMBER * config->GetFroude() * config->GetFroude()
           << "\n";
    }

    file << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << "\n";
    file << "Free-stream temperature (non-dim): " << config->GetTemperature_FreeStreamND() << "\n";
    file << "Free-stream pressure (non-dim): " << config->GetPressure_FreeStreamND() << "\n";
    file << "Free-stream density (non-dim): " << config->GetDensity_FreeStreamND() << "\n";

    if (nDim == 2) {
      file << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      file << config->GetVelocity_FreeStreamND()[1] << "). ";
    } else {
      file << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      file << config->GetVelocity_FreeStreamND()[1] << ", " << config->GetVelocity_FreeStreamND()[2] << "). ";
    }
    file << "Magnitude: " << config->GetModVel_FreeStreamND() << "\n";
    file << "Free-stream total energy per unit mass (non-dim): " << config->GetEnergy_FreeStreamND() << "\n";

    if (viscous) {
      file << "Free-stream viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << "\n";
      if (turbulent) {
        file << "Free-stream turb. kinetic energy (non-dim): " << config->GetTke_FreeStreamND() << "\n";
        file << "Free-stream specific dissipation (non-dim): " << config->GetOmega_FreeStreamND() << "\n";
      }
    }

    if (unsteady) {
      file << "Total time (non-dim): " << config->GetTotal_UnstTimeND() << "\n";
      file << "Time step (non-dim): " << config->GetDelta_UnstTimeND() << "\n";
    }

  } else {

    /*--- Incompressible version of the console output ---*/

    const bool energy = config->GetEnergy_Equation();
    const bool boussinesq = (config->GetKind_DensityModel() == INC_DENSITYMODEL::BOUSSINESQ);

    if (config->GetRef_Inc_NonDim() == DIMENSIONAL) {
      file << "Viscous and Inviscid flow: rho_ref, vel_ref, temp_ref, p_ref\n";
      file << "are set to 1.0 in order to perform a dimensional calculation.\n";
    } else if (config->GetRef_Inc_NonDim() == INITIAL_VALUES) {
      file << "Viscous and Inviscid flow: rho_ref, vel_ref, and temp_ref\n";
      file << "are based on the initial values, p_ref = rho_ref*vel_ref^2.\n";
    } else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
      file << "Viscous and Inviscid flow: rho_ref, vel_ref, and temp_ref\n";
      file << "are user-provided reference values, p_ref = rho_ref*vel_ref^2.\n";
    }
    if (dynamic_grid)
      file << "Force coefficients computed using MACH_MOTION.\n";
    else
      file << "Force coefficients computed using initial values.\n";

    file << "The reference area for force coeffs. is " << config->GetRefArea() << " m^2.\n";
    file << "The reference length for force coeffs. is " << config->GetRefLength() << " m.\n";

    file << "The pressure is decomposed into thermodynamic and dynamic components.\n";
    file << "The initial value of the dynamic pressure is 0.\n";

    file << "Mach number: " << config->GetMach();
    if (config->GetKind_FluidModel() == CONSTANT_DENSITY) {
      file << ", computed using the Bulk modulus.\n";
    } else {
      file << ", computed using fluid speed of sound.\n";
    }
    file << "For external flows, the initial state is imposed at the far-field.\n";
    file << "Angle of attack (deg): " << config->GetAoA() << ", computed using the initial velocity.\n";
    file << "Side slip angle (deg): " << config->GetAoS() << ", computed using the initial velocity.\n";

    if (viscous) {
      file << "Reynolds number per meter: " << config->GetReynolds() << ", computed using initial values.\n";
      file << "Reynolds number is a byproduct of inputs only (not used internally).\n";
    }
    file << "SI units only. The grid should be dimensional (meters).\n";

    switch (config->GetKind_DensityModel()) {
      case INC_DENSITYMODEL::CONSTANT:
        if (energy)
          file << "Energy equation is active and decoupled.\n";
        else
          file << "No energy equation.\n";
        break;

      case INC_DENSITYMODEL::BOUSSINESQ:
        if (energy) file << "Energy equation is active and coupled through Boussinesq approx.\n";
        break;

      case INC_DENSITYMODEL::VARIABLE:
        if (energy) file << "Energy equation is active and coupled for variable density.\n";
        break;

      case INC_DENSITYMODEL::FLAMELET:
        file << "Density is obtained through flamelet manifold.\n";
        break;
    }

    file << "-- Input conditions:\n";

    switch (config->GetKind_FluidModel()) {
      case CONSTANT_DENSITY:
        file << "Fluid Model: CONSTANT_DENSITY \n";
        if (energy) {
          file << "Specific heat at constant pressure (Cp): " << config->GetSpecific_Heat_Cp() << " N.m/kg.K.\n";
        }
        if (boussinesq) file << "Thermal expansion coefficient: " << config->GetThermal_Expansion_Coeff() << " K^-1.\n";
        file << "Thermodynamic pressure not required.\n";
        break;

      case INC_IDEAL_GAS:
        file << "Fluid Model: INC_IDEAL_GAS \n";
        file << "Variable density incompressible flow using ideal gas law.\n";
        file << "Density is a function of temperature (constant thermodynamic pressure).\n";
        file << "Specific heat at constant pressure (Cp): " << config->GetSpecific_Heat_Cp() << " N.m/kg.K.\n";
        file << "Molecular weight : " << config->GetMolecular_Weight() << " g/mol\n";
        file << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K.\n";
        file << "Thermodynamic pressure: " << config->GetPressure_Thermodynamic();
        if (si_units) file << " Pa.\n";
        else file << " psf.\n";
        break;

      case FLUID_MIXTURE:
        file << "Fluid Model: FLUID_MIXTURE \n";
        file << "Variable density incompressible flow using ideal gas law.\n";
        file << "Density is a function of temperature (constant thermodynamic pressure).\n";
        file << "Specific heat at constant pressure (Cp): " << config->GetSpecific_Heat_Cp() << " N.m/kg.K.\n";
        file << "Molecular weight : " << config->GetMolecular_Weight() << " g/mol\n";
        file << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K.\n";
        file << "Thermodynamic pressure: " << config->GetPressure_Thermodynamic();
        if (si_units) file << " Pa.\n";
        else file << " psf.\n";
        break;

      case FLUID_FLAMELET:
        file << "Fluid model: FLUID_FLAMELET \n";
        if (si_units) file << " Pa.\n";
        else file << " psf.\n";
        break;

      case INC_IDEAL_GAS_POLY:
        file << "Fluid Model: INC_IDEAL_GAS_POLY \n";
        file << "Variable density incompressible flow using ideal gas law.\n";
        file << "Density is a function of temperature (constant thermodynamic pressure).\n";
        file << "Molecular weight: " << config->GetMolecular_Weight() << " g/mol.\n";
        file << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K.\n";
        file << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << "\n";
        file << "Thermodynamic pressure: " << config->GetPressure_Thermodynamic();
        if (si_units) file << " Pa.\n";
        else file << " psf.\n";
        file << "Cp(T) polynomial coefficients: \n  (";
        for (unsigned short iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
          file << config->GetCp_PolyCoeff(iVar);
          if (iVar < config->GetnPolyCoeffs() - 1) file << ", ";
        }
        file << ").\n";
        file << "Cp(T) polynomial coefficients (non-dim.): \n  (";
        for (unsigned short iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
          file << config->GetCp_PolyCoeffND(iVar);
          if (iVar < config->GetnPolyCoeffs() - 1) file << ", ";
        }
        file << ").\n";
        break;
    }
    if (viscous) {
      switch (config->GetKind_ViscosityModel()) {
        case VISCOSITYMODEL::CONSTANT:
          file << "Viscosity Model: CONSTANT_VISCOSITY  \n";
          file << "Constant Laminar Viscosity: " << config->GetMu_Constant();
          if (si_units) file << " N.s/m^2.\n";
          else file << " lbf.s/ft^2.\n";
          file << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND() << "\n";
          break;

        case VISCOSITYMODEL::FLAMELET:
          file << "Viscosity Model: FLAMELET  \n";
          if (si_units) file << " N.s/m^2.\n";
          else file << " lbf.s/ft^2.\n";
          file << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND() << "\n";
          break;

        case VISCOSITYMODEL::COOLPROP:
          file << "Viscosity Model: CoolProp \n";
          break;

        case VISCOSITYMODEL::SUTHERLAND:
          file << "Viscosity Model: SUTHERLAND \n";
          file << "Ref. Laminar Viscosity: " << config->GetMu_Ref();
          if (si_units) file << " N.s/m^2.\n";
          else file << " lbf.s/ft^2.\n";
          file << "Ref. Temperature: " << config->GetMu_Temperature_Ref();
          if (si_units) file << " K.\n";
          else file << " R.\n";
          file << "Sutherland Constant: " << config->GetMu_S();
          if (si_units) file << " K.\n";
          else file << " R.\n";
          file << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND() << "\n";
          file << "Ref. Temperature (non-dim): " << config->GetMu_Temperature_RefND() << "\n";
          file << "Sutherland constant (non-dim): " << config->GetMu_SND() << "\n";
          break;

        case VISCOSITYMODEL::POLYNOMIAL:
          file << "Viscosity Model: POLYNOMIAL_VISCOSITY  \n";
          file << "Mu(T) polynomial coefficients: \n  (";
          for (unsigned short iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
            file << config->GetMu_PolyCoeff(iVar);
            if (iVar < config->GetnPolyCoeffs() - 1) file << ", ";
          }
          file << ").\n";
          file << "Mu(T) polynomial coefficients (non-dim.): \n  (";
          for (unsigned short iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
            file << config->GetMu_PolyCoeffND(iVar);
            if (iVar < config->GetnPolyCoeffs() - 1) file << ", ";
          }
          file << ").\n";
          break;
      }

      if (energy) {
        switch (config->GetKind_ConductivityModel()) {
          case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
            file << "Conductivity Model: CONSTANT_PRANDTL  \n";
            file << "Prandtl (Laminar): " << config->GetPrandtl_Lam() << "\n";
            break;

          case CONDUCTIVITYMODEL::CONSTANT:
            file << "Conductivity Model: CONSTANT \n";
            file << "Molecular Conductivity: " << config->GetThermal_Conductivity_Constant() << " W/m^2.K.\n";
            file << "Molecular Conductivity (non-dim): " << config->GetThermal_Conductivity_ConstantND() << "\n";
            break;
          case CONDUCTIVITYMODEL::COOLPROP:
            file << "Conductivity Model: COOLPROP \n";
            break;

          case CONDUCTIVITYMODEL::FLAMELET:
            file << "Conductivity Model: FLAMELET \n";
            file << "Molecular Conductivity units: "  << " W/m^2.K.\n";
            file << "Molecular Conductivity (non-dim): " << config->GetThermal_Conductivity_ConstantND() << "\n";
            break;

          case CONDUCTIVITYMODEL::POLYNOMIAL:
            file << "Viscosity Model: POLYNOMIAL \n";
            file << "Kt(T) polynomial coefficients: \n  (";
            for (unsigned short iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
              file << config->GetKt_PolyCoeff(iVar);
              if (iVar < config->GetnPolyCoeffs() - 1) file << ", ";
            }
            file << ").\n";
            file << "Kt(T) polynomial coefficients (non-dim.): \n  (";
            for (unsigned short iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
              file << config->GetKt_PolyCoeffND(iVar);
              if (iVar < config->GetnPolyCoeffs() - 1) file << ", ";
            }
            file << ").\n";
            break;
        }

        if (turbulent) {
          switch (config->GetKind_ConductivityModel_Turb()) {
            case CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL:
              file << "Turbulent Conductivity Model: CONSTANT_PRANDTL  \n";
              file << "Turbulent Prandtl: " << config->GetPrandtl_Turb() << "\n";
              break;
            case CONDUCTIVITYMODEL_TURB::NONE:
              file << "Turbulent Conductivity Model: CONDUCTIVITYMODEL_TURB::NONE \n";
              file << "No turbulent component in effective thermal conductivity.\n";
              break;
          }
        }
      }
    }

    if (config->GetKind_FluidModel() == CONSTANT_DENSITY) {
      file << "Bulk modulus: " << config->GetBulk_Modulus();
      if (si_units) file << " Pa.\n";
      else file << " psf.\n";
    }

    file << "Initial dynamic pressure: " << config->GetPressure_FreeStream();
    if (si_units) file << " Pa.\n";
    else file << " psf.\n";

    file << "Initial total pressure: "
         << config->GetPressure_FreeStream() +
                0.5 * config->GetDensity_FreeStream() * config->GetModVel_FreeStream() * config->GetModVel_FreeStream();
    if (si_units) file << " Pa.\n";
    else file << " psf.\n";

    if (energy) {
      file << "Initial temperature: " << config->GetTemperature_FreeStream();
      if (si_units) file << " K.\n";
      else file << " R.\n";
    }

    file << "Initial density: " << config->GetDensity_FreeStream();
    if (si_units) file << " kg/m^3.\n";
    else file << " slug/ft^3.\n";

    file << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0];
    file << ", " << config->GetVelocity_FreeStream()[1];
    if (nDim == 3) {
      file << ", " << config->GetVelocity_FreeStream()[2];
    }
    if (si_units) file << ") m/s. ";
    else file << ") ft/s. ";

    file << "Magnitude: " << config->GetModVel_FreeStream();
    if (si_units) file << " m/s.\n";
    else file << " ft/s.\n";

    if (viscous) {
      file << "Initial laminar viscosity: " << config->GetViscosity_FreeStream();
      if (si_units) file << " N.s/m^2.\n";
      else file << " lbf.s/ft^2.\n";
      if (turbulent) {
        file << "Initial turb. kinetic energy per unit mass: " << config->GetTke_FreeStream();
        if (si_units) file << " m^2/s^2.\n";
        else file << " ft^2/s^2.\n";
        file << "Initial specific dissipation: " << config->GetOmega_FreeStream();
        if (si_units) file << " 1/s.\n";
        else file << " 1/s.\n";
      }
    }

    if (unsteady) {
      file << "Total time: " << config->GetTotal_UnstTime() << " s. Time step: " << config->GetDelta_UnstTime()
           << " s.\n";
    }

    /*--- Print out reference values. ---*/

    file << "-- Reference values:\n";

    if (config->GetKind_FluidModel() != CONSTANT_DENSITY) {
      file << "Reference specific gas constant: " << config->GetGas_Constant_Ref();
      if (si_units) file << " N.m/kg.K.\n";
      else file << " lbf.ft/slug.R.\n";
    } else {
      if (energy) {
        file << "Reference specific heat: " << config->GetGas_Constant_Ref();
        if (si_units) file << " N.m/kg.K.\n";
        else file << " lbf.ft/slug.R.\n";
      }
    }

    file << "Reference pressure: " << config->GetPressure_Ref();
    if (si_units) file << " Pa.\n";
    else file << " psf.\n";

    if (energy) {
      file << "Reference temperature: " << config->GetTemperature_Ref();
      if (si_units) file << " K.\n";
      else file << " R.\n";
    }

    file << "Reference density: " << config->GetDensity_Ref();
    if (si_units) file << " kg/m^3.\n";
    else file << " slug/ft^3.\n";

    file << "Reference velocity: " << config->GetVelocity_Ref();
    if (si_units) file << " m/s.\n";
    else file << " ft/s.\n";

    file << "Reference length: " << config->GetLength_Ref();
    if (si_units) file << " m.\n";
    else file << " in.\n";

    if (viscous) {
      file << "Reference viscosity: " << config->GetViscosity_Ref();
      if (si_units) file << " N.s/m^2.\n";
      else file << " lbf.s/ft^2.\n";
    }

    if (unsteady) file << "Reference time: " << config->GetTime_Ref() << " s.\n";

    /*--- Print out resulting non-dim values here. ---*/

    file << "-- Resulting non-dimensional state:\n";
    file << "Mach number (non-dim): " << config->GetMach() << "\n";
    if (viscous) {
      file << "Reynolds number (per m): " << config->GetReynolds() << "\n";
    }

    if (config->GetKind_FluidModel() != CONSTANT_DENSITY) {
      file << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << "\n";
      file << "Initial thermodynamic pressure (non-dim): " << config->GetPressure_ThermodynamicND() << "\n";
    } else {
      if (energy) {
        file << "Specific heat at constant pressure (non-dim): " << config->GetSpecific_Heat_CpND() << "\n";
        if (boussinesq)
          file << "Thermal expansion coefficient (non-dim.): " << config->GetThermal_Expansion_CoeffND() << " K^-1.\n";
      }
    }

    if (energy) file << "Initial temperature (non-dim): " << config->GetTemperature_FreeStreamND() << "\n";
    file << "Initial pressure (non-dim): " << config->GetPressure_FreeStreamND() << "\n";
    file << "Initial density (non-dim): " << config->GetDensity_FreeStreamND() << "\n";

    file << "Initial velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0];
    file << ", " << config->GetVelocity_FreeStreamND()[1];
    if (nDim == 3) {
      file << ", " << config->GetVelocity_FreeStreamND()[2];
    }
    file << "). Magnitude: " << config->GetModVel_FreeStreamND() << "\n";

    if (viscous) {
      file << "Initial viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << "\n";
      if (turbulent) {
        file << "Initial turb. kinetic energy (non-dim): " << config->GetTke_FreeStreamND() << "\n";
        file << "Initial specific dissipation (non-dim): " << config->GetOmega_FreeStreamND() << "\n";
      }
    }

    if (unsteady) {
      file << "Total time (non-dim): " << config->GetTotal_UnstTimeND() << "\n";
      file << "Time step (non-dim): " << config->GetDelta_UnstTimeND() << "\n";
    }
  }

  /*--- Begin forces breakdown info. ---*/

  file << fixed;
  file << "\n\nForces breakdown:\n\n";

  if (nDim == 3) {
    su2double m = flow_solver->GetTotal_CFz() / flow_solver->GetTotal_CFx();
    su2double term = (Total_CoPz / m) - Total_CoPx;

    if (term > 0)
      file << "Center of Pressure: X=" << 1 / m << "Z-" << term << ".\n\n";
    else
      file << "Center of Pressure: X=" << 1 / m << "Z+" << fabs(term);
    if (si_units) file << " m.\n\n";
    else file << " in.\n\n";
  } else {
    su2double m = flow_solver->GetTotal_CFy() / flow_solver->GetTotal_CFx();
    su2double term = (Total_CoPy / m) - Total_CoPx;
    if (term > 0)
      file << "Center of Pressure: X=" << 1 / m << "Y-" << term << ".\n\n";
    else
      file << "Center of Pressure: X=" << 1 / m << "Y+" << fabs(term);
    if (si_units) file << " m.\n\n";
    else file << " in.\n\n";
  }

  /*--- Reference area and force factors. ---*/

  const su2double Factor = flow_solver->GetAeroCoeffsReferenceForce();
  const su2double Ref = config->GetDensity_Ref() * pow(config->GetVelocity_Ref(), 2);

  file << "NOTE: Multiply forces by the non-dimensional factor: " << Factor << ", and the reference factor: " << Ref
       << "\nto obtain the dimensional force.\n\n";

  file << "Total CL:    ";
  file.width(11);
  file << Total_CL;
  file << " | Pressure (";
  file.width(5);
  file << SU2_TYPE::Int((Inv_CL * 100.0) / (Total_CL + EPS));
  file << "%): ";
  file.width(11);
  file << Inv_CL;
  file << " | Friction (";
  file.width(5);
  file << SU2_TYPE::Int((Visc_CL * 100.0) / (Total_CL + EPS));
  file << "%): ";
  file.width(11);
  file << Visc_CL;
  file << " | Momentum (";
  file.width(5);
  file << SU2_TYPE::Int((Mnt_CL * 100.0) / (Total_CL + EPS));
  file << "%): ";
  file.width(11);
  file << Mnt_CL << "\n";

  file << "Total CD:    ";
  file.width(11);
  file << Total_CD;
  file << " | Pressure (";
  file.width(5);
  file << SU2_TYPE::Int((Inv_CD * 100.0) / (Total_CD + EPS)) << "%): ";
  file.width(11);
  file << Inv_CD;
  file << " | Friction (";
  file.width(5);
  file << SU2_TYPE::Int((Visc_CD * 100.0) / (Total_CD + EPS)) << "%): ";
  file.width(11);
  file << Visc_CD;
  file << " | Momentum (";
  file.width(5);
  file << SU2_TYPE::Int((Mnt_CD * 100.0) / (Total_CD + EPS)) << "%): ";
  file.width(11);
  file << Mnt_CD << "\n";

  if (nDim == 3) {
    file << "Total CSF:   ";
    file.width(11);
    file << Total_CSF;
    file << " | Pressure (";
    file.width(5);
    file << SU2_TYPE::Int((Inv_CSF * 100.0) / (Total_CSF + EPS));
    file << "%): ";
    file.width(11);
    file << Inv_CSF;
    file << " | Friction (";
    file.width(5);
    file << SU2_TYPE::Int((Visc_CSF * 100.0) / (Total_CSF + EPS));
    file << "%): ";
    file.width(11);
    file << Visc_CSF;
    file << " | Momentum (";
    file.width(5);
    file << SU2_TYPE::Int((Mnt_CSF * 100.0) / (Total_CSF + EPS));
    file << "%): ";
    file.width(11);
    file << Mnt_CSF << "\n";
  }

  file << "Total CL/CD: ";
  file.width(11);
  file << Total_CEff;
  file << " | Pressure (";
  file.width(5);
  file << SU2_TYPE::Int((Inv_CEff * 100.0) / (Total_CEff + EPS));
  file << "%): ";
  file.width(11);
  file << Inv_CEff;
  file << " | Friction (";
  file.width(5);
  file << SU2_TYPE::Int((Visc_CEff * 100.0) / (Total_CEff + EPS));
  file << "%): ";
  file.width(11);
  file << Visc_CEff;
  file << " | Momentum (";
  file.width(5);
  file << SU2_TYPE::Int((Mnt_CEff * 100.0) / (Total_CEff + EPS));
  file << "%): ";
  file.width(11);
  file << Mnt_CEff << "\n";

  if (nDim == 3) {
    file << "Total CMx:   ";
    file.width(11);
    file << Total_CMx;
    file << " | Pressure (";
    file.width(5);
    file << SU2_TYPE::Int((Inv_CMx * 100.0) / (Total_CMx + EPS));
    file << "%): ";
    file.width(11);
    file << Inv_CMx;
    file << " | Friction (";
    file.width(5);
    file << SU2_TYPE::Int((Visc_CMx * 100.0) / (Total_CMx + EPS));
    file << "%): ";
    file.width(11);
    file << Visc_CMx;
    file << " | Momentum (";
    file.width(5);
    file << SU2_TYPE::Int((Mnt_CMx * 100.0) / (Total_CMx + EPS));
    file << "%): ";
    file.width(11);
    file << Mnt_CMx << "\n";

    file << "Total CMy:   ";
    file.width(11);
    file << Total_CMy;
    file << " | Pressure (";
    file.width(5);
    file << SU2_TYPE::Int((Inv_CMy * 100.0) / (Total_CMy + EPS));
    file << "%): ";
    file.width(11);
    file << Inv_CMy;
    file << " | Friction (";
    file.width(5);
    file << SU2_TYPE::Int((Visc_CMy * 100.0) / (Total_CMy + EPS));
    file << "%): ";
    file.width(11);
    file << Visc_CMy;
    file << " | Momentum (";
    file.width(5);
    file << SU2_TYPE::Int((Mnt_CMz * 100.0) / (Total_CMz + EPS));
    file << "%): ";
    file.width(11);
    file << Mnt_CMy << "\n";
  }

  file << "Total CMz:   ";
  file.width(11);
  file << Total_CMz;
  file << " | Pressure (";
  file.width(5);
  file << SU2_TYPE::Int((Inv_CMz * 100.0) / (Total_CMz + EPS));
  file << "%): ";
  file.width(11);
  file << Inv_CMz;
  file << " | Friction (";
  file.width(5);
  file << SU2_TYPE::Int((Visc_CMz * 100.0) / (Total_CMz + EPS));
  file << "%): ";
  file.width(11);
  file << Visc_CMz;
  file << " | Momentum (";
  file.width(5);
  file << SU2_TYPE::Int((Mnt_CMz * 100.0) / (Total_CMz + EPS));
  file << "%): ";
  file.width(11);
  file << Mnt_CMz << "\n";

  file << "Total CFx:   ";
  file.width(11);
  file << Total_CFx;
  file << " | Pressure (";
  file.width(5);
  file << SU2_TYPE::Int((Inv_CFx * 100.0) / (Total_CFx + EPS));
  file << "%): ";
  file.width(11);
  file << Inv_CFx;
  file << " | Friction (";
  file.width(5);
  file << SU2_TYPE::Int((Visc_CFx * 100.0) / (Total_CFx + EPS));
  file << "%): ";
  file.width(11);
  file << Visc_CFx;
  file << " | Momentum (";
  file.width(5);
  file << SU2_TYPE::Int((Mnt_CFx * 100.0) / (Total_CFx + EPS));
  file << "%): ";
  file.width(11);
  file << Mnt_CFx << "\n";

  file << "Total CFy:   ";
  file.width(11);
  file << Total_CFy;
  file << " | Pressure (";
  file.width(5);
  file << SU2_TYPE::Int((Inv_CFy * 100.0) / (Total_CFy + EPS));
  file << "%): ";
  file.width(11);
  file << Inv_CFy;
  file << " | Friction (";
  file.width(5);
  file << SU2_TYPE::Int((Visc_CFy * 100.0) / (Total_CFy + EPS));
  file << "%): ";
  file.width(11);
  file << Visc_CFy;
  file << " | Momentum (";
  file.width(5);
  file << SU2_TYPE::Int((Mnt_CFy * 100.0) / (Total_CFy + EPS));
  file << "%): ";
  file.width(11);
  file << Mnt_CFy << "\n";

  if (nDim == 3) {
    file << "Total CFz:   ";
    file.width(11);
    file << Total_CFz;
    file << " | Pressure (";
    file.width(5);
    file << SU2_TYPE::Int((Inv_CFz * 100.0) / (Total_CFz + EPS));
    file << "%): ";
    file.width(11);
    file << Inv_CFz;
    file << " | Friction (";
    file.width(5);
    file << SU2_TYPE::Int((Visc_CFz * 100.0) / (Total_CFz + EPS));
    file << "%): ";
    file.width(11);
    file << Visc_CFz;
    file << " | Momentum (";
    file.width(5);
    file << SU2_TYPE::Int((Mnt_CFz * 100.0) / (Total_CFz + EPS));
    file << "%): ";
    file.width(11);
    file << Mnt_CFz << "\n";
  }

  file << "\n\n";

  for (auto iMarker = 0u; iMarker < nMonitoring; iMarker++) {
    file << "Surface name: " << config->GetMarker_Monitoring_TagBound(iMarker) << "\n\n";

    file << "Total CL    (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CL[iMarker] * 100.0) / (Total_CL + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CL[iMarker];
    file << " | Pressure (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CL_Inv[iMarker] * 100.0) / (Surface_CL[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CL_Inv[iMarker];
    file << " | Friction (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CL_Visc[iMarker] * 100.0) / (Surface_CL[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CL_Visc[iMarker];
    file << " | Momentum (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CL_Mnt[iMarker] * 100.0) / (Surface_CL[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CL_Mnt[iMarker] << "\n";

    file << "Total CD    (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CD[iMarker] * 100.0) / (Total_CD + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CD[iMarker];
    file << " | Pressure (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CD_Inv[iMarker] * 100.0) / (Surface_CD[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CD_Inv[iMarker];
    file << " | Friction (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CD_Visc[iMarker] * 100.0) / (Surface_CD[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CD_Visc[iMarker];
    file << " | Momentum (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CD_Mnt[iMarker] * 100.0) / (Surface_CD[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CD_Mnt[iMarker] << "\n";

    if (nDim == 3) {
      file << "Total CSF   (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CSF[iMarker] * 100.0) / (Total_CSF + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CSF[iMarker];
      file << " | Pressure (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CSF_Inv[iMarker] * 100.0) / (Surface_CSF[iMarker] + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CSF_Inv[iMarker];
      file << " | Friction (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CSF_Visc[iMarker] * 100.0) / (Surface_CSF[iMarker] + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CSF_Visc[iMarker];
      file << " | Momentum (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CSF_Mnt[iMarker] * 100.0) / (Surface_CSF[iMarker] + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CSF_Mnt[iMarker] << "\n";
    }

    file << "Total CL/CD (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CEff[iMarker] * 100.0) / (Total_CEff + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CEff[iMarker];
    file << " | Pressure (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CEff_Inv[iMarker] * 100.0) / (Surface_CEff[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CEff_Inv[iMarker];
    file << " | Friction (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CEff_Visc[iMarker] * 100.0) / (Surface_CEff[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CEff_Visc[iMarker];
    file << " | Momentum (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CEff_Mnt[iMarker] * 100.0) / (Surface_CEff[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CEff_Mnt[iMarker] << "\n";

    if (nDim == 3) {
      file << "Total CMx   (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CMx[iMarker] * 100.0) / (Total_CMx + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CMx[iMarker];
      file << " | Pressure (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CMx_Inv[iMarker] * 100.0) / (Surface_CMx[iMarker] + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CMx_Inv[iMarker];
      file << " | Friction (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CMx_Visc[iMarker] * 100.0) / (Surface_CMx[iMarker] + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CMx_Visc[iMarker];
      file << " | Momentum (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CMx_Mnt[iMarker] * 100.0) / (Surface_CMx[iMarker] + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CMx_Mnt[iMarker] << "\n";

      file << "Total CMy   (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CMy[iMarker] * 100.0) / (Total_CMy + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CMy[iMarker];
      file << " | Pressure (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CMy_Inv[iMarker] * 100.0) / (Surface_CMy[iMarker] + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CMy_Inv[iMarker];
      file << " | Friction (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CMy_Visc[iMarker] * 100.0) / (Surface_CMy[iMarker] + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CMy_Visc[iMarker];
      file << " | Momentum (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CMy_Mnt[iMarker] * 100.0) / (Surface_CMy[iMarker] + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CMy_Mnt[iMarker] << "\n";
    }

    file << "Total CMz   (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CMz[iMarker] * 100.0) / (Total_CMz + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CMz[iMarker];
    file << " | Pressure (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CMz_Inv[iMarker] * 100.0) / (Surface_CMz[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CMz_Inv[iMarker];
    file << " | Friction (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CMz_Visc[iMarker] * 100.0) / (Surface_CMz[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CMz_Visc[iMarker];
    file << " | Momentum (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CMz_Mnt[iMarker] * 100.0) / (Surface_CMz[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CMz_Mnt[iMarker] << "\n";

    file << "Total CFx   (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CFx[iMarker] * 100.0) / (Total_CFx + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CFx[iMarker];
    file << " | Pressure (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CFx_Inv[iMarker] * 100.0) / (Surface_CFx[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CFx_Inv[iMarker];
    file << " | Friction (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CFx_Visc[iMarker] * 100.0) / (Surface_CFx[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CFx_Visc[iMarker];
    file << " | Momentum (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CFx_Mnt[iMarker] * 100.0) / (Surface_CFx[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CFx_Mnt[iMarker] << "\n";

    file << "Total CFy   (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CFy[iMarker] * 100.0) / (Total_CFy + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CFy[iMarker];
    file << " | Pressure (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CFy_Inv[iMarker] * 100.0) / (Surface_CFy[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CFy_Inv[iMarker];
    file << " | Friction (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CFy_Visc[iMarker] * 100.0) / (Surface_CFy[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CFy_Visc[iMarker];
    file << " | Momentum (";
    file.width(5);
    file << SU2_TYPE::Int((Surface_CFy_Mnt[iMarker] * 100.0) / (Surface_CFy[iMarker] + EPS));
    file << "%): ";
    file.width(11);
    file << Surface_CFy_Mnt[iMarker] << "\n";

    if (nDim == 3) {
      file << "Total CFz   (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CFz[iMarker] * 100.0) / (Total_CFz + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CFz[iMarker];
      file << " | Pressure (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CFz_Inv[iMarker] * 100.0) / (Surface_CFz[iMarker] + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CFz_Inv[iMarker];
      file << " | Friction (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CFz_Visc[iMarker] * 100.0) / (Surface_CFz[iMarker] + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CFz_Visc[iMarker];
      file << " | Momentum (";
      file.width(5);
      file << SU2_TYPE::Int((Surface_CFz_Mnt[iMarker] * 100.0) / (Surface_CFz[iMarker] + EPS));
      file << "%): ";
      file.width(11);
      file << Surface_CFz_Mnt[iMarker] << "\n";
    }

    file << "\n";
  }
  // clang-format on
}

bool CFlowOutput::WriteVolumeOutput(CConfig *config, unsigned long Iter, bool force_writing, unsigned short iFile){

  bool writeRestart = false;
  auto FileFormat = config->GetVolumeOutputFiles();

  if (config->GetTime_Domain()){
    if (((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) || (config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING)) &&
        ((Iter == 0) || (Iter % config->GetVolumeOutputFrequency(iFile) == 0))){
      return true;
    }

    /* check if we want to write a restart file*/
    if (FileFormat[iFile] == OUTPUT_TYPE::RESTART_ASCII || FileFormat[iFile] == OUTPUT_TYPE::RESTART_BINARY || FileFormat[iFile] == OUTPUT_TYPE::CSV) {
      writeRestart = true;
    }

    /* only write 'double' files for the restart files */
    if ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND) &&
      ((Iter == 0) || (Iter % config->GetVolumeOutputFrequency(iFile) == 0) ||
      (((Iter+1) % config->GetVolumeOutputFrequency(iFile) == 0) && writeRestart) || // Restarts need 2 old solutions.
      (((Iter+2) == config->GetnTime_Iter()) && writeRestart))){      // The last timestep is written anyway but one needs the step before for restarts.
      return true;
    }
  } else {
    if (config->GetFixed_CL_Mode() && config->GetFinite_Difference_Mode()) return false;
    return ((Iter > 0) && Iter % config->GetVolumeOutputFrequency(iFile) == 0) || force_writing;
  }

  return force_writing;
}

void CFlowOutput::SetTimeAveragedFields() {
  AddVolumeOutput("MEAN_DENSITY", "MeanDensity", "TIME_AVERAGE", "Mean density");
  AddVolumeOutput("MEAN_VELOCITY-X", "MeanVelocity_x", "TIME_AVERAGE", "Mean velocity x-component");
  AddVolumeOutput("MEAN_VELOCITY-Y", "MeanVelocity_y", "TIME_AVERAGE", "Mean velocity y-component");
  if (nDim == 3)
    AddVolumeOutput("MEAN_VELOCITY-Z", "MeanVelocity_z", "TIME_AVERAGE", "Mean velocity z-component");

  AddVolumeOutput("MEAN_PRESSURE", "MeanPressure", "TIME_AVERAGE", "Mean pressure");
  AddVolumeOutput("RMS_U",   "RMS[u]", "TIME_AVERAGE", "RMS u");
  AddVolumeOutput("RMS_V",   "RMS[v]", "TIME_AVERAGE", "RMS v");
  AddVolumeOutput("RMS_UV",  "RMS[uv]", "TIME_AVERAGE", "RMS uv");
  AddVolumeOutput("RMS_P",   "RMS[Pressure]",   "TIME_AVERAGE", "RMS Pressure");
  AddVolumeOutput("UUPRIME", "u'u'", "TIME_AVERAGE", "Mean Reynolds-stress component u'u'");
  AddVolumeOutput("VVPRIME", "v'v'", "TIME_AVERAGE", "Mean Reynolds-stress component v'v'");
  AddVolumeOutput("UVPRIME", "u'v'", "TIME_AVERAGE", "Mean Reynolds-stress component u'v'");
  AddVolumeOutput("PPRIME",  "p'p'",   "TIME_AVERAGE", "Mean pressure fluctuation p'p'");
  if (nDim == 3){
    AddVolumeOutput("RMS_W",   "RMS[w]", "TIME_AVERAGE", "RMS u");
    AddVolumeOutput("RMS_UW", "RMS[uw]", "TIME_AVERAGE", "RMS uw");
    AddVolumeOutput("RMS_VW", "RMS[vw]", "TIME_AVERAGE", "RMS vw");
    AddVolumeOutput("WWPRIME", "w'w'", "TIME_AVERAGE", "Mean Reynolds-stress component w'w'");
    AddVolumeOutput("UWPRIME", "w'u'", "TIME_AVERAGE", "Mean Reynolds-stress component w'u'");
    AddVolumeOutput("VWPRIME", "w'v'", "TIME_AVERAGE", "Mean Reynolds-stress component w'v'");
  }
}

void CFlowOutput::LoadTimeAveragedData(unsigned long iPoint, const CVariable *Node_Flow){
  SetAvgVolumeOutputValue("MEAN_DENSITY", iPoint, Node_Flow->GetDensity(iPoint));
  SetAvgVolumeOutputValue("MEAN_VELOCITY-X", iPoint, Node_Flow->GetVelocity(iPoint,0));
  SetAvgVolumeOutputValue("MEAN_VELOCITY-Y", iPoint, Node_Flow->GetVelocity(iPoint,1));
  if (nDim == 3)
    SetAvgVolumeOutputValue("MEAN_VELOCITY-Z", iPoint, Node_Flow->GetVelocity(iPoint,2));

  SetAvgVolumeOutputValue("MEAN_PRESSURE", iPoint, Node_Flow->GetPressure(iPoint));

  SetAvgVolumeOutputValue("RMS_U", iPoint, pow(Node_Flow->GetVelocity(iPoint,0),2));
  SetAvgVolumeOutputValue("RMS_V", iPoint, pow(Node_Flow->GetVelocity(iPoint,1),2));
  SetAvgVolumeOutputValue("RMS_UV", iPoint, Node_Flow->GetVelocity(iPoint,0) * Node_Flow->GetVelocity(iPoint,1));
  SetAvgVolumeOutputValue("RMS_P", iPoint, pow(Node_Flow->GetPressure(iPoint),2));
  if (nDim == 3){
    SetAvgVolumeOutputValue("RMS_W", iPoint, pow(Node_Flow->GetVelocity(iPoint,2),2));
    SetAvgVolumeOutputValue("RMS_VW", iPoint, Node_Flow->GetVelocity(iPoint,2) * Node_Flow->GetVelocity(iPoint,1));
    SetAvgVolumeOutputValue("RMS_UW", iPoint,  Node_Flow->GetVelocity(iPoint,2) * Node_Flow->GetVelocity(iPoint,0));
  }

  const su2double umean  = GetVolumeOutputValue("MEAN_VELOCITY-X", iPoint);
  const su2double uumean = GetVolumeOutputValue("RMS_U", iPoint);
  const su2double vmean  = GetVolumeOutputValue("MEAN_VELOCITY-Y", iPoint);
  const su2double vvmean = GetVolumeOutputValue("RMS_V", iPoint);
  const su2double uvmean = GetVolumeOutputValue("RMS_UV", iPoint);
  const su2double pmean  = GetVolumeOutputValue("MEAN_PRESSURE", iPoint);
  const su2double ppmean = GetVolumeOutputValue("RMS_P", iPoint);

  SetVolumeOutputValue("UUPRIME", iPoint, -(umean*umean - uumean));
  SetVolumeOutputValue("VVPRIME", iPoint, -(vmean*vmean - vvmean));
  SetVolumeOutputValue("UVPRIME", iPoint, -(umean*vmean - uvmean));
  SetVolumeOutputValue("PPRIME",  iPoint, -(pmean*pmean - ppmean));
  if (nDim == 3){
    const su2double wmean  = GetVolumeOutputValue("MEAN_VELOCITY-Z", iPoint);
    const su2double wwmean = GetVolumeOutputValue("RMS_W", iPoint);
    const su2double uwmean = GetVolumeOutputValue("RMS_UW", iPoint);
    const su2double vwmean = GetVolumeOutputValue("RMS_VW", iPoint);
    SetVolumeOutputValue("WWPRIME", iPoint, -(wmean*wmean - wwmean));
    SetVolumeOutputValue("UWPRIME", iPoint, -(umean*wmean - uwmean));
    SetVolumeOutputValue("VWPRIME",  iPoint, -(vmean*wmean - vwmean));
  }
}

void CFlowOutput::SetFixedCLScreenOutput(const CConfig *config){
  PrintingToolbox::CTablePrinter FixedCLSummary(&cout);

  if (fabs(historyOutput_Map["CL_DRIVER_COMMAND"].value) > 1e-16){
    FixedCLSummary.AddColumn("Fixed CL Mode", 40);
    FixedCLSummary.AddColumn("Value", 30);
    FixedCLSummary.SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    FixedCLSummary.PrintHeader();
    FixedCLSummary << "Current CL" << historyOutput_Map["LIFT"].value;
    FixedCLSummary << "Target CL" << config->GetTarget_CL();
    FixedCLSummary << "Previous AOA" << historyOutput_Map["PREV_AOA"].value;
    if (config->GetFinite_Difference_Mode()){
      FixedCLSummary << "Changed AoA by (Finite Difference step)" << historyOutput_Map["CL_DRIVER_COMMAND"].value;
      lastInnerIter = curInnerIter - 1;
    }
    else
      FixedCLSummary << "Changed AoA by" << historyOutput_Map["CL_DRIVER_COMMAND"].value;
    FixedCLSummary.PrintFooter();
    SetScreenHeader(config);
  }

  else if (config->GetFinite_Difference_Mode() && historyOutput_Map["AOA"].value == historyOutput_Map["PREV_AOA"].value){
    FixedCLSummary.AddColumn("Fixed CL Mode (Finite Difference)", 40);
    FixedCLSummary.AddColumn("Value", 30);
    FixedCLSummary.SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    FixedCLSummary.PrintHeader();
    FixedCLSummary << "Delta CL / Delta AoA" << config->GetdCL_dAlpha();
    FixedCLSummary << "Delta CD / Delta CL" << config->GetdCD_dCL();
    if (nDim == 3){
      FixedCLSummary << "Delta CMx / Delta CL" << config->GetdCMx_dCL();
      FixedCLSummary << "Delta CMy / Delta CL" << config->GetdCMy_dCL();
    }
    FixedCLSummary << "Delta CMz / Delta CL" << config->GetdCMz_dCL();
    FixedCLSummary.PrintFooter();
    curInnerIter = lastInnerIter;
    WriteMetaData(config);
    curInnerIter = config->GetInnerIter();
  }
}
