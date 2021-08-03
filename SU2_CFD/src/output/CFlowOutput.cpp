/*!
 * \file CFlowOutput.cpp
 * \brief Main subroutines for compressible flow output
 * \author R. Sanchez
 * \version 7.1.1 "Blackbird"
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

#include "../../include/output/CFlowOutput.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/solvers/CSolver.hpp"

CFlowOutput::CFlowOutput(CConfig *config, unsigned short nDim, bool fem_output) : CFVMOutput (config, nDim, fem_output){

  lastInnerIter = curInnerIter;
}

void CFlowOutput::AddAnalyzeSurfaceOutput(CConfig *config){


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
  AddHistoryOutput("SURFACE_PRESSURE_DROP",    "Pressure_Drop",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total pressure drop on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average mass fraction of CO    
  AddHistoryOutput("AVG_CO",                   "Avg_CO",                    ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average mass fraction of CO on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average mass fraction of NO    
  AddHistoryOutput("AVG_NOX",                  "Avg_NOx",                   ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average mass fraction of NO on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average temperature
  AddHistoryOutput("AVG_TEMP",                 "Avg_Temp",                  ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average temperature on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
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
  /// DESCRIPTION: Pressure drop
  AddHistoryOutputPerSurface("SURFACE_PRESSURE_DROP",    "Pressure_Drop",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average mass fraction of CO    
  AddHistoryOutputPerSurface("AVG_CO",                   "Avg_CO",                    ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average mass fraction of NO    
  AddHistoryOutputPerSurface("AVG_NOX",                  "Avg_NOx",                   ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average temperature    
  AddHistoryOutputPerSurface("AVG_TEMP",                 "Avg_Temp",                  ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// END_GROUP

}

void CFlowOutput::SetAnalyzeSurface(CSolver **solver, CGeometry *geometry, CConfig *config, bool output){

  unsigned short iDim, iMarker, iMarker_Analyze;
  unsigned long iVertex, iPoint;
  su2double Mach = 0.0, Pressure, Temperature = 0.0, TotalPressure = 0.0, TotalTemperature = 0.0,
  Enthalpy, Velocity[3] = {0.0}, TangVel[3], Vector[3], Velocity2, MassFlow, Density, Area,
  AxiFactor = 1.0, SoundSpeed, Vn, Vn2, Vtang2, Weight = 1.0;

  const su2double Gas_Constant      = config->GetGas_ConstantND();
  const su2double Gamma             = config->GetGamma();
  const unsigned short nMarker      = config->GetnMarker_All();
  const unsigned short nDim         = geometry->GetnDim();
  const unsigned short Kind_Average = config->GetKind_Average();

  const bool compressible   = config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE;
  const bool incompressible = config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE;
  const bool energy         = config->GetEnergy_Equation();
  const bool flamelet_model = config->GetKind_Scalar_Model() == PROGRESS_VARIABLE;
  const bool streamwisePeriodic = (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE);

  const bool axisymmetric               = config->GetAxisymmetric();
  const unsigned short nMarker_Analyze  = config->GetnMarker_Analyze();

  CSolver* flow_solver   = solver[FLOW_SOL];
  CSolver* scalar_solver = solver[SCALAR_SOL];

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
  vector<su2double> Surface_CO                (nMarker,0.0);
  vector<su2double> Surface_NOx               (nMarker,0.0);

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
  su2double  Tot_Surface_PressureDrop      = 0.0;
  su2double  Tot_Surface_CO                = 0.0;
  su2double  Tot_Surface_NOx               = 0.0;
  su2double  Tot_Surface_Temp              = 0.0;
  //su2double  Tot_Surface_Scalar[n_scalars];
  //for (int i_scalar = 0; i_scalar < n_scalars; ++i_scalar)
  //  Tot_Surface_Scalar[i_scalar] = 0.0;

  /*--- Compute the numerical fan face Mach number, and the total area of the inflow ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if (config->GetMarker_All_Analyze(iMarker) == YES) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->nodes->GetDomain(iPoint)) {

          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);

          if (axisymmetric) {
            if (geometry->nodes->GetCoord(iPoint, 1) != 0.0)
              AxiFactor = 2.0*PI_NUMBER*geometry->nodes->GetCoord(iPoint, 1);
            else {
              /*--- Find the point "above" by finding the neighbor of iPoint that is also a vertex of iMarker. ---*/
              AxiFactor = 0.0;
              for (unsigned short iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); ++iNeigh) {
                auto jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
                if (geometry->nodes->GetVertex(jPoint, iMarker) >= 0) {
                  /*--- Not multiplied by two since we need to half the y coordinate. ---*/
                  AxiFactor = PI_NUMBER * geometry->nodes->GetCoord(jPoint, 1);
                  break;
                }
              }
            }
          } else {
            AxiFactor = 1.0;
          }

          Density = flow_solver->GetNodes()->GetDensity(iPoint);
          Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0; Vn = 0.0; Vtang2 = 0.0;

          for (iDim = 0; iDim < nDim; iDim++) {
            Area += (Vector[iDim] * AxiFactor) * (Vector[iDim] * AxiFactor);
            Velocity[iDim] = flow_solver->GetNodes()->GetVelocity(iPoint,iDim);
            Velocity2 += Velocity[iDim] * Velocity[iDim];
            Vn += Velocity[iDim] * Vector[iDim] * AxiFactor;
            MassFlow += Vector[iDim] * AxiFactor * Density * Velocity[iDim];
          }

          Area       = sqrt (Area);
          if (AxiFactor == 0.0) Vn = 0.0; else Vn /= Area;
          Vn2        = Vn * Vn;

          Pressure   = flow_solver->GetNodes()->GetPressure(iPoint);
          /*--- Use recovered pressure here as pressure difference between in and outlet is zero otherwise  ---*/
          if(streamwisePeriodic) Pressure = flow_solver->GetNodes()->GetStreamwise_Periodic_RecoveredPressure(iPoint);
          SoundSpeed = flow_solver->GetNodes()->GetSoundSpeed(iPoint);

          for (iDim = 0; iDim < nDim; iDim++) {
            TangVel[iDim] = Velocity[iDim] - Vn*Vector[iDim]*AxiFactor/Area;
            Vtang2       += TangVel[iDim]*TangVel[iDim];
          }

          if (incompressible){
            if (config->GetKind_DensityModel() == INC_DENSITYMODEL::VARIABLE) {
              Mach = sqrt(flow_solver->GetNodes()->GetVelocity2(iPoint))/
              sqrt(flow_solver->GetNodes()->GetSpecificHeatCp(iPoint)*config->GetPressure_ThermodynamicND()/(flow_solver->GetNodes()->GetSpecificHeatCv(iPoint)*flow_solver->GetNodes()->GetDensity(iPoint)));
            } else {
              Mach = sqrt(flow_solver->GetNodes()->GetVelocity2(iPoint))/
              sqrt(config->GetBulk_Modulus()/(flow_solver->GetNodes()->GetDensity(iPoint)));
            }
            Temperature       = flow_solver->GetNodes()->GetTemperature(iPoint);
            Enthalpy          = flow_solver->GetNodes()->GetSpecificHeatCp(iPoint)*Temperature;
            TotalTemperature  = Temperature + 0.5*Velocity2/flow_solver->GetNodes()->GetSpecificHeatCp(iPoint);
            TotalPressure     = Pressure + 0.5*Density*Velocity2;
          }
          else{
            Mach              = sqrt(Velocity2)/SoundSpeed;
            Temperature       = Pressure / (Gas_Constant * Density);
            Enthalpy          = flow_solver->GetNodes()->GetEnthalpy(iPoint);
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
          if (flamelet_model){
            Surface_CO[iMarker]  += scalar_solver->GetNodes()->GetSolution(iPoint, I_CO) * Weight;
            Surface_NOx[iMarker] += scalar_solver->GetNodes()->GetSolution(iPoint, I_NOX) * Weight;
          }

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
  vector<su2double> Surface_CO_Local               (nMarker_Analyze,0.0);
  vector<su2double> Surface_NOx_Local               (nMarker_Analyze,0.0);
  
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
  vector<su2double> Surface_CO_Total                (nMarker_Analyze,0.0);
  vector<su2double> Surface_NOx_Total               (nMarker_Analyze,0.0);

  vector<su2double> Surface_MomentumDistortion_Total (nMarker_Analyze,0.0);

  /*--- Compute the numerical fan face Mach number, mach number, temperature and the total area ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

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
          Surface_CO_Local[iMarker_Analyze]                += Surface_CO[iMarker];
          Surface_NOx_Local[iMarker_Analyze]               += Surface_NOx[iMarker];
        }

      }

    }

  }

  auto Allreduce = [](const vector<su2double>& src, vector<su2double>& dst) {
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
  Allreduce(Surface_CO_Local, Surface_CO_Total);
  Allreduce(Surface_NOx_Local,Surface_NOx_Total);
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
      Surface_CO_Total[iMarker_Analyze]               /= Weight;
      Surface_NOx_Total[iMarker_Analyze]              /= Weight;
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
      Surface_CO_Total[iMarker_Analyze]               = 0.0;
      Surface_NOx_Total[iMarker_Analyze]              = 0.0;
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
    if (config->GetSystemMeasurements() == US) MassFlow *= 32.174;
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

    su2double y_CO = Surface_CO_Total[iMarker_Analyze];
    SetHistoryOutputPerSurfaceValue("AVG_CO", y_CO, iMarker_Analyze);
    Tot_Surface_CO += y_CO;
    config->SetSurface_CO(iMarker_Analyze, y_CO);
    
    su2double y_NOx = Surface_NOx_Total[iMarker_Analyze];
    SetHistoryOutputPerSurfaceValue("AVG_NOX", y_NOx, iMarker_Analyze);
    Tot_Surface_NOx += y_NOx;
    config->SetSurface_NOx(iMarker_Analyze, y_NOx);

    su2double temp = Surface_TotalTemperature_Total[iMarker_Analyze];
    SetHistoryOutputPerSurfaceValue("AVG_TEMP", temp, iMarker_Analyze);
    Tot_Surface_Temp += temp;
    config->SetSurface_Temperature(iMarker_Analyze, temp);
  }

  /*--- Compute the average static pressure drop between two surfaces. Note
   that this assumes we have two surfaces being analyzed and that the outlet
   is first followed by the inlet. This is because we may also want to choose
   outlet values (temperature, uniformity, etc.) for our design problems,
   which require the outlet to be listed first. This is a simple first version
   that could be generalized to a different orders/lists/etc. ---*/

  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
    su2double Pressure_Drop = 0.0;
    if (nMarker_Analyze == 2) {
      Pressure_Drop = (Surface_Pressure_Total[1]-Surface_Pressure_Total[0]) * config->GetPressure_Ref();
      config->SetSurface_PressureDrop(iMarker_Analyze, Pressure_Drop);
    }
    SetHistoryOutputPerSurfaceValue("SURFACE_PRESSURE_DROP",  Pressure_Drop, iMarker_Analyze);
    Tot_Surface_PressureDrop += Pressure_Drop;
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
  SetHistoryOutputValue("SURFACE_PRESSURE_DROP", Tot_Surface_PressureDrop);
  SetHistoryOutputValue("AVG_CO",   Tot_Surface_CO);
  SetHistoryOutputValue("AVG_NOX",  Tot_Surface_NOx);
  SetHistoryOutputValue("AVG_TEMP", Tot_Surface_Temp);

  if ((rank == MASTER_NODE) && !config->GetDiscrete_Adjoint() && output) {

    cout.precision(6);
    cout.setf(ios::scientific, ios::floatfield);
    cout << endl << "Computing surface mean values." << endl << endl;

    for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
      cout << "Surface "<< config->GetMarker_Analyze_TagBound(iMarker_Analyze) << ":" << endl;

      if (nDim == 3) { if (config->GetSystemMeasurements() == SI) cout << setw(20) << "Area (m^2): "; else cout << setw(20) << "Area (ft^2): "; }
      else { if (config->GetSystemMeasurements() == SI) cout << setw(20) << "Area (m): "; else cout << setw(20) << "Area (ft): "; }

      if (config->GetSystemMeasurements() == SI)      cout << setw(15) << fabs(Surface_Area_Total[iMarker_Analyze]);
      else if (config->GetSystemMeasurements() == US) cout << setw(15) << fabs(Surface_Area_Total[iMarker_Analyze])*12.0*12.0;

      cout << endl;

      su2double MassFlow = config->GetSurface_MassFlow(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "Mf (kg/s): " << setw(15) << MassFlow;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "Mf (lbs/s): " << setw(15) << MassFlow;

      su2double NormalVelocity = config->GetSurface_NormalVelocity(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "Vn (m/s): " << setw(15) << NormalVelocity;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "Vn (ft/s): " << setw(15) << NormalVelocity;

      cout << endl;

      su2double Uniformity = config->GetSurface_Uniformity(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "Uniformity (m/s): " << setw(15) << Uniformity;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "Uniformity (ft/s): " << setw(15) << Uniformity;

      su2double SecondaryStrength = config->GetSurface_SecondaryStrength(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "Secondary (m/s): " << setw(15) << SecondaryStrength;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "Secondary (ft/s): " << setw(15) << SecondaryStrength;

      cout << endl;

      su2double MomentumDistortion = config->GetSurface_MomentumDistortion(iMarker_Analyze);
      cout << setw(20) << "Mom. Distortion: " << setw(15) << MomentumDistortion;

      su2double SecondOverUniform = config->GetSurface_SecondOverUniform(iMarker_Analyze);
      cout << setw(20) << "Second/Uniform: " << setw(15) << SecondOverUniform;

      cout << endl;

      su2double Pressure = config->GetSurface_Pressure(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "P (Pa): " << setw(15) << Pressure;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "P (psf): " << setw(15) << Pressure;

      su2double TotalPressure = config->GetSurface_TotalPressure(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "PT (Pa): " << setw(15) <<TotalPressure;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "PT (psf): " << setw(15) <<TotalPressure;

      cout << endl;

      su2double Mach = config->GetSurface_Mach(iMarker_Analyze);
      cout << setw(20) << "Mach: " << setw(15) << Mach;

      su2double Density = config->GetSurface_Density(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "Rho (kg/m^3): " << setw(15) << Density;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "Rho (lb/ft^3): " << setw(15) << Density*32.174;

      cout << endl;

      if (compressible || energy) {
        su2double Temperature = config->GetSurface_Temperature(iMarker_Analyze);
        if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "T (K): " << setw(15) << Temperature;
        else if (config->GetSystemMeasurements() == US) cout << setw(20) << "T (R): " << setw(15) << Temperature;

        su2double TotalTemperature = config->GetSurface_TotalTemperature(iMarker_Analyze);
        if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "TT (K): " << setw(15) << TotalTemperature;
        else if (config->GetSystemMeasurements() == US) cout << setw(20) << "TT (R): " << setw(15) << TotalTemperature;

        cout << endl;
      }
    }
    cout.unsetf(ios_base::floatfield);

  }

  std::cout << std::resetiosflags(std::cout.flags());
}

void CFlowOutput::AddAerodynamicCoefficients(CConfig *config){

  /// BEGIN_GROUP: AERO_COEFF, DESCRIPTION: Sum of the aerodynamic coefficients and forces on all surfaces (markers) set with MARKER_MONITORING.
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
  /// DESCRIPTION: Custom objective
  AddHistoryOutput("CUSTOM_OBJFUNC", "Custom_ObjFunc", ScreenOutputFormat::FIXED, "AERO_COEFF", "Custom objective function on all surfaces set with MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
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

void CFlowOutput::SetAerodynamicCoefficients(CConfig *config, CSolver *flow_solver){

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
  SetHistoryOutputValue("CUSTOM_OBJFUNC", flow_solver->GetTotal_Custom_ObjFunc());

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

  SetHistoryOutputValue("COMBO", flow_solver->GetTotal_ComboObj());
}

void CFlowOutput::SetRotatingFrameCoefficients(CConfig *config, CSolver *flow_solver) {

  SetHistoryOutputValue("THRUST", flow_solver->GetTotal_CT());
  SetHistoryOutputValue("TORQUE", flow_solver->GetTotal_CQ());
  SetHistoryOutputValue("FIGURE_OF_MERIT", flow_solver->GetTotal_CMerit());
}


void CFlowOutput::Add_CpInverseDesignOutput(){

  AddHistoryOutput("INVERSE_DESIGN_PRESSURE", "Cp_Diff", ScreenOutputFormat::FIXED, "CP_DIFF", "Cp difference for inverse design");

}

void CFlowOutput::Set_CpInverseDesign(CSolver *solver, const CGeometry *geometry, const CConfig *config){

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

void CFlowOutput::WriteAdditionalFiles(CConfig *config, CGeometry *geometry, CSolver **solver_container){

  if (config->GetFixed_CL_Mode() || config->GetFixed_CM_Mode()){
    WriteMetaData(config);
  }

  if (config->GetWrt_ForcesBreakdown()){
    WriteForcesBreakdown(config, geometry, solver_container);
  }

}

void CFlowOutput::WriteMetaData(const CConfig *config){

  ofstream meta_file;

  string filename = "flow";

  filename = config->GetFilename(filename, ".meta", curTimeIter);

  /*--- All processors open the file. ---*/

  if (rank == MASTER_NODE) {
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


    if (( config->GetKind_Solver() == DISC_ADJ_EULER ||
          config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES ||
          config->GetKind_Solver() == DISC_ADJ_RANS )) {
      meta_file << "SENS_AOA=" << GetHistoryFieldValue("SENS_AOA") * PI_NUMBER / 180.0 << endl;
    }
  }

  meta_file.close();
}

void CFlowOutput::WriteForcesBreakdown(CConfig *config, CGeometry *geometry, CSolver **solver_container){

  unsigned short iMarker_Monitoring;

  const bool compressible    = (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);
  const bool incompressible  = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
  const bool unsteady        = config->GetTime_Domain();
  const bool viscous         = config->GetViscous();
  const bool dynamic_grid    = config->GetDynamic_Grid();
  const bool gravity         = config->GetGravityForce();
  const bool turbulent       = config->GetKind_Solver() == RANS;
  const bool fixed_cl        = config->GetFixed_CL_Mode();
  const auto Kind_Solver     = config->GetKind_Solver();
  const auto Kind_Turb_Model = config->GetKind_Turb_Model();
  const auto Ref_NonDim      = config->GetRef_NonDim();

  const auto nDim =  geometry->GetnDim();

  auto fileName = config->GetBreakdown_FileName();
  if (unsteady) {
    const auto lastindex = fileName.find_last_of(".");
    const auto ext = fileName.substr(lastindex, fileName.size());
    fileName = fileName.substr(0, lastindex);
    fileName = config->GetFilename(fileName, ext, curTimeIter);
  }

  /*--- Output the mean flow solution using only the master node ---*/

  if ( rank == MASTER_NODE) {

    cout << endl << "Writing the forces breakdown file ("<< fileName << ")." << endl;

    /*--- Initialize variables to store information from all domains (direct solution) ---*/

    su2double Total_CL = 0.0, Total_CD = 0.0, Total_CSF = 0.0,
    Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0,
    Total_CoPx = 0.0, Total_CoPy = 0.0, Total_CoPz = 0.0,
    Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0, Inv_CL = 0.0,
    Inv_CD = 0.0, Inv_CSF = 0.0, Inv_CMx = 0.0, Inv_CMy = 0.0,
    Inv_CMz = 0.0, Inv_CEff = 0.0, Inv_CFx = 0.0, Inv_CFy = 0.0, Inv_CFz =
    0.0,      Mnt_CL = 0.0,
    Mnt_CD = 0.0, Mnt_CSF = 0.0, Mnt_CMx = 0.0, Mnt_CMy = 0.0,
    Mnt_CMz = 0.0, Mnt_CEff = 0.0, Mnt_CFx = 0.0, Mnt_CFy = 0.0, Mnt_CFz =
    0.0, Visc_CL = 0.0,
    Visc_CD = 0.0, Visc_CSF = 0.0, Visc_CMx = 0.0, Visc_CMy = 0.0,
    Visc_CMz = 0.0, Visc_CEff = 0.0, Visc_CFx = 0.0, Visc_CFy = 0.0, Visc_CFz =
    0.0, *Surface_CL = nullptr, *Surface_CD = nullptr,
    *Surface_CSF = nullptr, *Surface_CEff = nullptr, *Surface_CFx = nullptr,
    *Surface_CFy = nullptr, *Surface_CFz = nullptr,
    *Surface_CMx = nullptr, *Surface_CMy = nullptr, *Surface_CMz = nullptr,
    *Surface_CL_Inv = nullptr,
    *Surface_CD_Inv = nullptr, *Surface_CSF_Inv = nullptr,
    *Surface_CEff_Inv = nullptr, *Surface_CFx_Inv = nullptr, *Surface_CFy_Inv =
    nullptr, *Surface_CFz_Inv = nullptr, *Surface_CMx_Inv = nullptr,
    *Surface_CMy_Inv = nullptr, *Surface_CMz_Inv = nullptr,
    *Surface_CL_Visc = nullptr,
    *Surface_CD_Visc = nullptr, *Surface_CSF_Visc = nullptr,
    *Surface_CEff_Visc = nullptr, *Surface_CFx_Visc = nullptr, *Surface_CFy_Visc =
    nullptr, *Surface_CFz_Visc = nullptr, *Surface_CMx_Visc = nullptr,
    *Surface_CMy_Visc = nullptr, *Surface_CMz_Visc = nullptr,
    *Surface_CL_Mnt = nullptr,
    *Surface_CD_Mnt = nullptr, *Surface_CSF_Mnt = nullptr,
    *Surface_CEff_Mnt = nullptr, *Surface_CFx_Mnt = nullptr, *Surface_CFy_Mnt =
    nullptr, *Surface_CFz_Mnt = nullptr, *Surface_CMx_Mnt = nullptr,
    *Surface_CMy_Mnt = nullptr, *Surface_CMz_Mnt = nullptr;

    /*--- Allocate memory for the coefficients being monitored ---*/

    Surface_CL      = new su2double[config->GetnMarker_Monitoring()];
    Surface_CD      = new su2double[config->GetnMarker_Monitoring()];
    Surface_CSF = new su2double[config->GetnMarker_Monitoring()];
    Surface_CEff       = new su2double[config->GetnMarker_Monitoring()];
    Surface_CFx        = new su2double[config->GetnMarker_Monitoring()];
    Surface_CFy        = new su2double[config->GetnMarker_Monitoring()];
    Surface_CFz        = new su2double[config->GetnMarker_Monitoring()];
    Surface_CMx        = new su2double[config->GetnMarker_Monitoring()];
    Surface_CMy        = new su2double[config->GetnMarker_Monitoring()];
    Surface_CMz        = new su2double[config->GetnMarker_Monitoring()];

    Surface_CL_Inv      = new su2double[config->GetnMarker_Monitoring()];
    Surface_CD_Inv      = new su2double[config->GetnMarker_Monitoring()];
    Surface_CSF_Inv = new su2double[config->GetnMarker_Monitoring()];
    Surface_CEff_Inv       = new su2double[config->GetnMarker_Monitoring()];
    Surface_CFx_Inv        = new su2double[config->GetnMarker_Monitoring()];
    Surface_CFy_Inv        = new su2double[config->GetnMarker_Monitoring()];
    Surface_CFz_Inv        = new su2double[config->GetnMarker_Monitoring()];
    Surface_CMx_Inv        = new su2double[config->GetnMarker_Monitoring()];
    Surface_CMy_Inv        = new su2double[config->GetnMarker_Monitoring()];
    Surface_CMz_Inv        = new su2double[config->GetnMarker_Monitoring()];

    Surface_CL_Visc = new su2double[config->GetnMarker_Monitoring()];
    Surface_CD_Visc = new su2double[config->GetnMarker_Monitoring()];
    Surface_CSF_Visc =
    new su2double[config->GetnMarker_Monitoring()];
    Surface_CEff_Visc = new su2double[config->GetnMarker_Monitoring()];
    Surface_CFx_Visc = new su2double[config->GetnMarker_Monitoring()];
    Surface_CFy_Visc = new su2double[config->GetnMarker_Monitoring()];
    Surface_CFz_Visc = new su2double[config->GetnMarker_Monitoring()];
    Surface_CMx_Visc = new su2double[config->GetnMarker_Monitoring()];
    Surface_CMy_Visc = new su2double[config->GetnMarker_Monitoring()];
    Surface_CMz_Visc = new su2double[config->GetnMarker_Monitoring()];


    Surface_CL_Mnt = new su2double[config->GetnMarker_Monitoring()];
    Surface_CD_Mnt = new su2double[config->GetnMarker_Monitoring()];
    Surface_CSF_Mnt =
    new su2double[config->GetnMarker_Monitoring()];
    Surface_CEff_Mnt = new su2double[config->GetnMarker_Monitoring()];
    Surface_CFx_Mnt = new su2double[config->GetnMarker_Monitoring()];
    Surface_CFy_Mnt = new su2double[config->GetnMarker_Monitoring()];
    Surface_CFz_Mnt = new su2double[config->GetnMarker_Monitoring()];
    Surface_CMx_Mnt = new su2double[config->GetnMarker_Monitoring()];
    Surface_CMy_Mnt = new su2double[config->GetnMarker_Monitoring()];
    Surface_CMz_Mnt = new su2double[config->GetnMarker_Monitoring()];

    /*--- Flow solution coefficients ---*/

    Total_CL       = solver_container[FLOW_SOL]->GetTotal_CL();
    Total_CD       = solver_container[FLOW_SOL]->GetTotal_CD();
    Total_CSF      = solver_container[FLOW_SOL]->GetTotal_CSF();
    Total_CEff        = solver_container[FLOW_SOL]->GetTotal_CEff();
    Total_CMx         = solver_container[FLOW_SOL]->GetTotal_CMx();
    Total_CMy         = solver_container[FLOW_SOL]->GetTotal_CMy();
    Total_CMz         = solver_container[FLOW_SOL]->GetTotal_CMz();
    Total_CFx         = solver_container[FLOW_SOL]->GetTotal_CFx();
    Total_CFy         = solver_container[FLOW_SOL]->GetTotal_CFy();
    Total_CFz         = solver_container[FLOW_SOL]->GetTotal_CFz();

    if (nDim == 2) {
      Total_CoPx = solver_container[FLOW_SOL]->GetTotal_CoPx() / solver_container[FLOW_SOL]->GetTotal_CFy();
      Total_CoPy = solver_container[FLOW_SOL]->GetTotal_CoPy() / solver_container[FLOW_SOL]->GetTotal_CFx();
      Total_CoPz = 0.0;
    }
    if (nDim == 3) {
      Total_CoPx = solver_container[FLOW_SOL]->GetTotal_CoPx() / solver_container[FLOW_SOL]->GetTotal_CFz();
      Total_CoPy = 0.0;
      Total_CoPz = solver_container[FLOW_SOL]->GetTotal_CoPz() / solver_container[FLOW_SOL]->GetTotal_CFx();
    }

    if (config->GetSystemMeasurements() == US) { Total_CoPx *= 12.0; Total_CoPy *= 12.0; Total_CoPz *= 12.0; }

    /*--- Flow inviscid solution coefficients ---*/

    Inv_CL =
    solver_container[FLOW_SOL]->GetAllBound_CL_Inv();
    Inv_CD =
    solver_container[FLOW_SOL]->GetAllBound_CD_Inv();
    Inv_CSF =
    solver_container[FLOW_SOL]->GetAllBound_CSF_Inv();
    Inv_CEff =
    solver_container[FLOW_SOL]->GetAllBound_CEff_Inv();
    Inv_CMx =
    solver_container[FLOW_SOL]->GetAllBound_CMx_Inv();
    Inv_CMy =
    solver_container[FLOW_SOL]->GetAllBound_CMy_Inv();
    Inv_CMz =
    solver_container[FLOW_SOL]->GetAllBound_CMz_Inv();
    Inv_CFx =
    solver_container[FLOW_SOL]->GetAllBound_CFx_Inv();
    Inv_CFy =
    solver_container[FLOW_SOL]->GetAllBound_CFy_Inv();
    Inv_CFz =
    solver_container[FLOW_SOL]->GetAllBound_CFz_Inv();

    /*--- Flow viscous solution coefficients ---*/

    Visc_CL =
    solver_container[FLOW_SOL]->GetAllBound_CL_Visc();
    Visc_CD =
    solver_container[FLOW_SOL]->GetAllBound_CD_Visc();
    Visc_CSF =
    solver_container[FLOW_SOL]->GetAllBound_CSF_Visc();
    Visc_CEff =
    solver_container[FLOW_SOL]->GetAllBound_CEff_Visc();
    Visc_CMx =
    solver_container[FLOW_SOL]->GetAllBound_CMx_Visc();
    Visc_CMy =
    solver_container[FLOW_SOL]->GetAllBound_CMy_Visc();
    Visc_CMz =
    solver_container[FLOW_SOL]->GetAllBound_CMz_Visc();
    Visc_CFx =
    solver_container[FLOW_SOL]->GetAllBound_CFx_Visc();
    Visc_CFy =
    solver_container[FLOW_SOL]->GetAllBound_CFy_Visc();
    Visc_CFz =
    solver_container[FLOW_SOL]->GetAllBound_CFz_Visc();

    /*--- Flow momentum solution coefficients ---*/

    Mnt_CL =
    solver_container[FLOW_SOL]->GetAllBound_CL_Mnt();
    Mnt_CD =
    solver_container[FLOW_SOL]->GetAllBound_CD_Mnt();
    Mnt_CSF =
    solver_container[FLOW_SOL]->GetAllBound_CSF_Mnt();
    Mnt_CEff =
    solver_container[FLOW_SOL]->GetAllBound_CEff_Mnt();
    Mnt_CMx =
    solver_container[FLOW_SOL]->GetAllBound_CMx_Mnt();
    Mnt_CMy =
    solver_container[FLOW_SOL]->GetAllBound_CMy_Mnt();
    Mnt_CMz =
    solver_container[FLOW_SOL]->GetAllBound_CMz_Mnt();
    Mnt_CFx =
    solver_container[FLOW_SOL]->GetAllBound_CFx_Mnt();
    Mnt_CFy =
    solver_container[FLOW_SOL]->GetAllBound_CFy_Mnt();
    Mnt_CFz =
    solver_container[FLOW_SOL]->GetAllBound_CFz_Mnt();


    /*--- Look over the markers being monitored and get the desired values ---*/

    for (iMarker_Monitoring = 0;
         iMarker_Monitoring < config->GetnMarker_Monitoring();
         iMarker_Monitoring++) {
      Surface_CL[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CL(
                                                             iMarker_Monitoring);
      Surface_CD[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CD(
                                                             iMarker_Monitoring);
      Surface_CSF[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CSF(
                                                              iMarker_Monitoring);
      Surface_CEff[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CEff(
                                                               iMarker_Monitoring);
      Surface_CMx[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CMx(
                                                              iMarker_Monitoring);
      Surface_CMy[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CMy(
                                                              iMarker_Monitoring);
      Surface_CMz[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CMz(
                                                              iMarker_Monitoring);
      Surface_CFx[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CFx(
                                                              iMarker_Monitoring);
      Surface_CFy[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CFy(
                                                              iMarker_Monitoring);
      Surface_CFz[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CFz(
                                                              iMarker_Monitoring);

      Surface_CL_Inv[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CL_Inv(
                                                                 iMarker_Monitoring);
      Surface_CD_Inv[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CD_Inv(
                                                                 iMarker_Monitoring);
      Surface_CSF_Inv[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CSF_Inv(
                                                                  iMarker_Monitoring);
      Surface_CEff_Inv[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CEff_Inv(
                                                                   iMarker_Monitoring);
      Surface_CMx_Inv[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CMx_Inv(
                                                                  iMarker_Monitoring);
      Surface_CMy_Inv[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CMy_Inv(
                                                                  iMarker_Monitoring);
      Surface_CMz_Inv[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CMz_Inv(
                                                                  iMarker_Monitoring);
      Surface_CFx_Inv[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CFx_Inv(
                                                                  iMarker_Monitoring);
      Surface_CFy_Inv[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CFy_Inv(
                                                                  iMarker_Monitoring);
      Surface_CFz_Inv[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CFz_Inv(
                                                                  iMarker_Monitoring);
      Surface_CL_Visc[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CL_Visc(
                                                                  iMarker_Monitoring);
      Surface_CD_Visc[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CD_Visc(
                                                                  iMarker_Monitoring);
      Surface_CSF_Visc[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CSF_Visc(
                                                                   iMarker_Monitoring);
      Surface_CEff_Visc[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CEff_Visc(
                                                                    iMarker_Monitoring);
      Surface_CMx_Visc[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CMx_Visc(
                                                                   iMarker_Monitoring);
      Surface_CMy_Visc[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CMy_Visc(
                                                                   iMarker_Monitoring);
      Surface_CMz_Visc[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CMz_Visc(
                                                                   iMarker_Monitoring);
      Surface_CFx_Visc[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CFx_Visc(
                                                                   iMarker_Monitoring);
      Surface_CFy_Visc[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CFy_Visc(
                                                                   iMarker_Monitoring);
      Surface_CFz_Visc[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CFz_Visc(
                                                                   iMarker_Monitoring);

      Surface_CL_Mnt[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CL_Mnt(
                                                                 iMarker_Monitoring);
      Surface_CD_Mnt[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CD_Mnt(
                                                                 iMarker_Monitoring);
      Surface_CSF_Mnt[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CSF_Mnt(
                                                                  iMarker_Monitoring);
      Surface_CEff_Mnt[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CEff_Mnt(
                                                                   iMarker_Monitoring);
      Surface_CMx_Mnt[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CMx_Mnt(
                                                                  iMarker_Monitoring);
      Surface_CMy_Mnt[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CMy_Mnt(
                                                                  iMarker_Monitoring);
      Surface_CMz_Mnt[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CMz_Mnt(
                                                                  iMarker_Monitoring);
      Surface_CFx_Mnt[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CFx_Mnt(
                                                                  iMarker_Monitoring);
      Surface_CFy_Mnt[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CFy_Mnt(
                                                                  iMarker_Monitoring);
      Surface_CFz_Mnt[iMarker_Monitoring] =
      solver_container[FLOW_SOL]->GetSurface_CFz_Mnt(
                                                                  iMarker_Monitoring);

    }

    /*--- Write file name with extension ---*/

    ofstream Breakdown_file;
    Breakdown_file.open(fileName);

    Breakdown_file << "\n" <<"-------------------------------------------------------------------------" << "\n";
    Breakdown_file << "|    ___ _   _ ___                                                      |" << "\n";
    Breakdown_file << "|   / __| | | |_  )   Release 7.1.1 \"Blackbird\"                       |" << "\n";
    Breakdown_file << "|   \\__ \\ |_| |/ /                                                    |" << "\n";
    Breakdown_file << "|   |___/\\___//___|   Suite (Computational Fluid Dynamics Code)        |" << "\n";
    Breakdown_file << "|                                                                       |" << "\n";
    //Breakdown_file << "|   Local date and time: " << dt << "                      |" << "\n";
    Breakdown_file << "-------------------------------------------------------------------------" << "\n";
    Breakdown_file << "| SU2 Project Website: https://su2code.github.io                        |" << "\n";
    Breakdown_file << "|                                                                       |" << "\n";
    Breakdown_file << "| The SU2 Project is maintained by the SU2 Foundation                   |" << "\n";
    Breakdown_file << "| (http://su2foundation.org)                                            |" << "\n";
    Breakdown_file << "-------------------------------------------------------------------------" << "\n";
    Breakdown_file << "| Copyright 2012-2021, SU2 Contributors                                 |" << "\n";
    Breakdown_file << "|                                                                       |" << "\n";
    Breakdown_file << "| SU2 is free software; you can redistribute it and/or                  |" << "\n";
    Breakdown_file << "| modify it under the terms of the GNU Lesser General Public            |" << "\n";
    Breakdown_file << "| License as published by the Free Software Foundation; either          |" << "\n";
    Breakdown_file << "| version 2.1 of the License, or (at your option) any later version.    |" << "\n";
    Breakdown_file << "|                                                                       |" << "\n";
    Breakdown_file << "| SU2 is distributed in the hope that it will be useful,                |" << "\n";
    Breakdown_file << "| but WITHOUT ANY WARRANTY; without even the implied warranty of        |" << "\n";
    Breakdown_file << "| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |" << "\n";
    Breakdown_file << "| Lesser General Public License for more details.                       |" << "\n";
    Breakdown_file << "|                                                                       |" << "\n";
    Breakdown_file << "| You should have received a copy of the GNU Lesser General Public      |" << "\n";
    Breakdown_file << "| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |" << "\n";
    Breakdown_file << "-------------------------------------------------------------------------" << "\n";

    Breakdown_file.precision(6);

    Breakdown_file << "\n" << "\n" << "Problem definition:" << "\n" << "\n";

    switch (Kind_Solver) {
      case EULER: case INC_EULER:
        if (compressible) Breakdown_file << "Compressible Euler equations." << "\n";
        if (incompressible) Breakdown_file << "Incompressible Euler equations." << "\n";
        break;
      case NAVIER_STOKES: case INC_NAVIER_STOKES:
        if (compressible) Breakdown_file << "Compressible Laminar Navier-Stokes' equations." << "\n";
        if (incompressible) Breakdown_file << "Incompressible Laminar Navier-Stokes' equations." << "\n";
        break;
      case RANS: case INC_RANS:
        if (compressible) Breakdown_file << "Compressible RANS equations." << "\n";
        if (incompressible) Breakdown_file << "Incompressible RANS equations." << "\n";
        Breakdown_file << "Turbulence model: ";
        switch (Kind_Turb_Model) {
          case SA:        Breakdown_file << "Spalart Allmaras" << "\n"; break;
          case SA_NEG:    Breakdown_file << "Negative Spalart Allmaras" << "\n"; break;
          case SA_E:      Breakdown_file << "Edwards Spalart Allmaras" << "\n"; break;
          case SA_COMP:   Breakdown_file << "Compressibility Correction Spalart Allmaras" << "\n"; break;
          case SA_E_COMP: Breakdown_file << "Compressibility Correction Edwards Spalart Allmaras" << "\n"; break;
          case SST:       Breakdown_file << "Menter's SST"     << "\n"; break;
          case SST_SUST:  Breakdown_file << "Menter's SST with sustaining terms" << "\n"; break;
        }
        break;
    }


    /*--- Compressible version of console output ---*/

    if (compressible) {


    if (compressible) {
      Breakdown_file << "Mach number: " << config->GetMach() <<"."<< "\n";
      Breakdown_file << "Angle of attack (AoA): " << config->GetAoA() <<" deg, and angle of sideslip (AoS): " << config->GetAoS() <<" deg."<< "\n";
      if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == INC_NAVIER_STOKES) ||
          (Kind_Solver == RANS) || (Kind_Solver == INC_RANS))
        Breakdown_file << "Reynolds number: " << config->GetReynolds() <<"."<< "\n";
    }

    if (fixed_cl) {
      Breakdown_file << "Simulation at a cte. CL: " << config->GetTarget_CL() << ".\n";
      Breakdown_file << "Approx. Delta CL / Delta AoA: " << config->GetdCL_dAlpha() << " (1/deg).\n";
      Breakdown_file << "Approx. Delta CD / Delta CL: " << config->GetdCD_dCL() << ".\n";
      if (nDim == 3 ) {
        Breakdown_file << "Approx. Delta CMx / Delta CL: " << config->GetdCMx_dCL() << ".\n";
        Breakdown_file << "Approx. Delta CMy / Delta CL: " << config->GetdCMy_dCL() << ".\n";
      }
      Breakdown_file << "Approx. Delta CMz / Delta CL: " << config->GetdCMz_dCL() << ".\n";
    }

    if (Ref_NonDim == DIMENSIONAL) { Breakdown_file << "Dimensional simulation." << "\n"; }
    else if (Ref_NonDim == FREESTREAM_PRESS_EQ_ONE) { Breakdown_file << "Non-Dimensional simulation (P=1.0, Rho=1.0, T=1.0 at the farfield)." << "\n"; }
    else if (Ref_NonDim == FREESTREAM_VEL_EQ_MACH) { Breakdown_file << "Non-Dimensional simulation (V=Mach, Rho=1.0, T=1.0 at the farfield)." << "\n"; }
    else if (Ref_NonDim == FREESTREAM_VEL_EQ_ONE) { Breakdown_file << "Non-Dimensional simulation (V=1.0, Rho=1.0, T=1.0 at the farfield)." << "\n"; }

    if (config->GetSystemMeasurements() == SI) {
      Breakdown_file << "The reference area is " << config->GetRefArea() << " m^2." << "\n";
      Breakdown_file << "The reference length is " << config->GetRefLength() << " m." << "\n";
    }

    if (config->GetSystemMeasurements() == US) {
      Breakdown_file << "The reference area is " << config->GetRefArea()*12.0*12.0 << " in^2." << "\n";
      Breakdown_file << "The reference length is " << config->GetRefLength()*12.0 << " in." << "\n";
    }
    Breakdown_file << "\n" << "\n" <<"Problem definition:" << "\n" << "\n";
    if (compressible) {
      if (viscous) {
        Breakdown_file << "Viscous flow: Computing pressure using the ideal gas law" << "\n";
        Breakdown_file << "based on the free-stream temperature and a density computed" << "\n";
        Breakdown_file << "from the Reynolds number." << "\n";
      } else {
        Breakdown_file << "Inviscid flow: Computing density based on free-stream" << "\n";
        Breakdown_file << "temperature and pressure using the ideal gas law." << "\n";
      }
    }

    if (dynamic_grid) Breakdown_file << "Force coefficients computed using MACH_MOTION." << "\n";
    else Breakdown_file << "Force coefficients computed using free-stream values." << "\n";

    if (incompressible) {
      Breakdown_file << "Viscous and Inviscid flow: rho_ref, and vel_ref" << "\n";
      Breakdown_file << "are based on the free-stream values, p_ref = rho_ref*vel_ref^2." << "\n";
      Breakdown_file << "The free-stream value of the pressure is 0." << "\n";
      Breakdown_file << "Mach number: "<< config->GetMach() << ", computed using the Bulk modulus." << "\n";
      Breakdown_file << "Angle of attack (deg): "<< config->GetAoA() << ", computed using the the free-stream velocity." << "\n";
      Breakdown_file << "Side slip angle (deg): "<< config->GetAoS() << ", computed using the the free-stream velocity." << "\n";
      if (viscous) Breakdown_file << "Reynolds number: " << config->GetReynolds() << ", computed using free-stream values."<< "\n";
      Breakdown_file << "Only dimensional computation, the grid should be dimensional." << "\n";
    }

    Breakdown_file <<"-- Input conditions:"<< "\n";

    if (compressible) {
      switch (config->GetKind_FluidModel()) {

        case STANDARD_AIR:
          Breakdown_file << "Fluid Model: STANDARD_AIR "<< "\n";
          Breakdown_file << "Specific gas constant: " << config->GetGas_Constant();
          if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.m/kg.K." << "\n";
          else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.ft/slug.R." << "\n";
          Breakdown_file << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< "\n";
          Breakdown_file << "Specific Heat Ratio: 1.4000 "<< "\n";
          break;

        case IDEAL_GAS:
          Breakdown_file << "Fluid Model: IDEAL_GAS "<< "\n";
          Breakdown_file << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << "\n";
          Breakdown_file << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< "\n";
          Breakdown_file << "Specific Heat Ratio: "<< config->GetGamma() << "\n";
          break;

        case VW_GAS:
          Breakdown_file << "Fluid Model: Van der Waals "<< "\n";
          Breakdown_file << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << "\n";
          Breakdown_file << "Specific gas constant (non-dim): " << config->GetGas_ConstantND()<< "\n";
          Breakdown_file << "Specific Heat Ratio: "<< config->GetGamma() << "\n";
          Breakdown_file << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << "\n";
          Breakdown_file << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << "\n";
          Breakdown_file << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << "\n";
          Breakdown_file << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << "\n";
          break;

        case PR_GAS:
          Breakdown_file << "Fluid Model: Peng-Robinson "<< "\n";
          Breakdown_file << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << "\n";
          Breakdown_file << "Specific gas constant(non-dim): " << config->GetGas_ConstantND()<< "\n";
          Breakdown_file << "Specific Heat Ratio: "<< config->GetGamma() << "\n";
          Breakdown_file << "Critical Pressure:   " << config->GetPressure_Critical()  << " Pa." << "\n";
          Breakdown_file << "Critical Temperature:  " << config->GetTemperature_Critical() << " K." << "\n";
          Breakdown_file << "Critical Pressure (non-dim):   " << config->GetPressure_Critical() /config->GetPressure_Ref() << "\n";
          Breakdown_file << "Critical Temperature (non-dim) :  " << config->GetTemperature_Critical() /config->GetTemperature_Ref() << "\n";
          break;
      }

      if (viscous) {

        switch (config->GetKind_ViscosityModel()) {

          case VISCOSITYMODEL::CONSTANT:
            Breakdown_file << "Viscosity Model: CONSTANT_VISCOSITY  "<< "\n";
            Breakdown_file << "Laminar Viscosity: " << config->GetMu_Constant();
            if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
            else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
            Breakdown_file << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< "\n";
            break;

          case VISCOSITYMODEL::SUTHERLAND:
            Breakdown_file << "Viscosity Model: SUTHERLAND "<< "\n";
            Breakdown_file << "Ref. Laminar Viscosity: " << config->GetMu_Ref();
            if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
            else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
            Breakdown_file << "Ref. Temperature: " << config->GetMu_Temperature_Ref();
            if (config->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
            else if (config->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
            Breakdown_file << "Sutherland Constant: "<< config->GetMu_S();
            if (config->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
            else if (config->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
            Breakdown_file << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< "\n";
            Breakdown_file << "Ref. Temperature (non-dim): " << config->GetMu_Temperature_RefND()<< "\n";
            Breakdown_file << "Sutherland constant (non-dim): "<< config->GetMu_SND()<< "\n";
            break;

          default:
            break;

        }
        switch (config->GetKind_ConductivityModel()) {

          case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
            Breakdown_file << "Conductivity Model: CONSTANT_PRANDTL "<< "\n";
            Breakdown_file << "Prandtl: " << config->GetPrandtl_Lam()<< "\n";
            break;

          case CONDUCTIVITYMODEL::CONSTANT:
            Breakdown_file << "Conductivity Model: CONSTANT "<< "\n";
            Breakdown_file << "Molecular Conductivity: " << config->GetThermal_Conductivity_Constant()<< " W/m^2.K." << "\n";
            Breakdown_file << "Molecular Conductivity (non-dim): " << config->GetThermal_Conductivity_ConstantND()<< "\n";
            break;

          default:
            break;

        }

        if ((Kind_Solver == RANS) || (Kind_Solver == INC_RANS)) {
          switch (config->GetKind_ConductivityModel_Turb()) {
            case CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL:
              Breakdown_file << "Turbulent Conductivity Model: CONSTANT_PRANDTL "<< "\n";
              Breakdown_file << "Turbulent Prandtl: " << config->GetPrandtl_Turb()<< "\n";
              break;
            case CONDUCTIVITYMODEL_TURB::NONE:
              Breakdown_file << "Turbulent Conductivity Model: NONE "<< "\n";
              Breakdown_file << "No turbulent component in effective thermal conductivity." << "\n";
              break;
          }
        }

      }
    }

    if (incompressible) {
      Breakdown_file << "Bulk modulus: " << config->GetBulk_Modulus();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";
      Breakdown_file << "Epsilon^2 multiplier of Beta for incompressible preconditioner: " << config->GetBeta_Factor();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";
    }

    Breakdown_file << "Free-stream static pressure: " << config->GetPressure_FreeStream();
    if (config->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
    else if (config->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";

    Breakdown_file << "Free-stream total pressure: " << config->GetPressure_FreeStream() * pow( 1.0+config->GetMach()*config->GetMach()*0.5*(config->GetGamma()-1.0), config->GetGamma()/(config->GetGamma()-1.0) );
    if (config->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
    else if (config->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";

    if (compressible) {
      Breakdown_file << "Free-stream temperature: " << config->GetTemperature_FreeStream();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";

      Breakdown_file << "Free-stream total temperature: " << config->GetTemperature_FreeStream() * (1.0 + config->GetMach() * config->GetMach() * 0.5 * (config->GetGamma() - 1.0));
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
    }

    Breakdown_file << "Free-stream density: " << config->GetDensity_FreeStream();
    if (config->GetSystemMeasurements() == SI) Breakdown_file << " kg/m^3." << "\n";
    else if (config->GetSystemMeasurements() == US) Breakdown_file << " slug/ft^3." << "\n";

    if (nDim == 2) {
      Breakdown_file << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      Breakdown_file << config->GetVelocity_FreeStream()[1] << ")";
    }
    if (nDim == 3) {
      Breakdown_file << "Free-stream velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
      Breakdown_file << config->GetVelocity_FreeStream()[1] << ", " << config->GetVelocity_FreeStream()[2] << ")";
    }
    if (config->GetSystemMeasurements() == SI) Breakdown_file << " m/s. ";
    else if (config->GetSystemMeasurements() == US) Breakdown_file << " ft/s. ";

    Breakdown_file << "Magnitude: "  << config->GetModVel_FreeStream();
    if (config->GetSystemMeasurements() == SI) Breakdown_file << " m/s." << "\n";
    else if (config->GetSystemMeasurements() == US) Breakdown_file << " ft/s." << "\n";

    if (compressible) {
      Breakdown_file << "Free-stream total energy per unit mass: " << config->GetEnergy_FreeStream();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " m^2/s^2." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " ft^2/s^2." << "\n";
    }

    if (viscous) {
      Breakdown_file << "Free-stream viscosity: " << config->GetViscosity_FreeStream();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
      if (turbulent) {
        Breakdown_file << "Free-stream turb. kinetic energy per unit mass: " << config->GetTke_FreeStream();
        if (config->GetSystemMeasurements() == SI) Breakdown_file << " m^2/s^2." << "\n";
        else if (config->GetSystemMeasurements() == US) Breakdown_file << " ft^2/s^2." << "\n";
        Breakdown_file << "Free-stream specific dissipation: " << config->GetOmega_FreeStream();
        if (config->GetSystemMeasurements() == SI) Breakdown_file << " 1/s." << "\n";
        else if (config->GetSystemMeasurements() == US) Breakdown_file << " 1/s." << "\n";
      }
    }

    if (unsteady) { Breakdown_file << "Total time: " << config->GetTotal_UnstTime() << " s. Time step: " << config->GetDelta_UnstTime() << " s." << "\n"; }

    /*--- Print out reference values. ---*/

    Breakdown_file <<"-- Reference values:"<< "\n";

    if (compressible) {
      Breakdown_file << "Reference specific gas constant: " << config->GetGas_Constant_Ref();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.m/kg.K." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.ft/slug.R." << "\n";
    }

    Breakdown_file << "Reference pressure: " << config->GetPressure_Ref();
    if (config->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
    else if (config->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";

    if (compressible) {
      Breakdown_file << "Reference temperature: " << config->GetTemperature_Ref();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
    }

    Breakdown_file << "Reference density: " << config->GetDensity_Ref();
    if (config->GetSystemMeasurements() == SI) Breakdown_file << " kg/m^3." << "\n";
    else if (config->GetSystemMeasurements() == US) Breakdown_file << " slug/ft^3." << "\n";

    Breakdown_file << "Reference velocity: " << config->GetVelocity_Ref();
    if (config->GetSystemMeasurements() == SI) Breakdown_file << " m/s." << "\n";
    else if (config->GetSystemMeasurements() == US) Breakdown_file << " ft/s." << "\n";

    if (compressible) {
      Breakdown_file << "Reference energy per unit mass: " << config->GetEnergy_Ref();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " m^2/s^2." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " ft^2/s^2." << "\n";
    }

    if (incompressible) {
      Breakdown_file << "Reference length: " << config->GetLength_Ref();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " m." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " in." << "\n";
    }

    if (viscous) {
      Breakdown_file << "Reference viscosity: " << config->GetViscosity_Ref();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
      if (compressible){
        Breakdown_file << "Reference conductivity: " << config->GetConductivity_Ref();
        if (config->GetSystemMeasurements() == SI) Breakdown_file << " W/m^2.K." << "\n";
        else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf/ft.s.R." << "\n";
      }
    }


    if (unsteady) Breakdown_file << "Reference time: " << config->GetTime_Ref() <<" s." << "\n";

    /*--- Print out resulting non-dim values here. ---*/

    Breakdown_file << "-- Resulting non-dimensional state:" << "\n";
    Breakdown_file << "Mach number (non-dim): " << config->GetMach() << "\n";
    if (viscous) {
      Breakdown_file << "Reynolds number (non-dim): " << config->GetReynolds() <<". Re length: " << config->GetLength_Reynolds();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " m." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " ft." << "\n";
    }
    if (gravity) {
      Breakdown_file << "Froude number (non-dim): " << config->GetFroude() << "\n";
      Breakdown_file << "Lenght of the baseline wave (non-dim): " << 2.0*PI_NUMBER*config->GetFroude()*config->GetFroude() << "\n";
    }

    if (compressible) {
      Breakdown_file << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << "\n";
      Breakdown_file << "Free-stream temperature (non-dim): " << config->GetTemperature_FreeStreamND() << "\n";
    }

    Breakdown_file << "Free-stream pressure (non-dim): " << config->GetPressure_FreeStreamND() << "\n";

    Breakdown_file << "Free-stream density (non-dim): " << config->GetDensity_FreeStreamND() << "\n";

    if (nDim == 2) {
      Breakdown_file << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      Breakdown_file << config->GetVelocity_FreeStreamND()[1] << "). ";
    } else {
      Breakdown_file << "Free-stream velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
      Breakdown_file << config->GetVelocity_FreeStreamND()[1] << ", " << config->GetVelocity_FreeStreamND()[2] << "). ";
    }
    Breakdown_file << "Magnitude: "   << config->GetModVel_FreeStreamND() << "\n";

    if (compressible)
      Breakdown_file << "Free-stream total energy per unit mass (non-dim): " << config->GetEnergy_FreeStreamND() << "\n";

    if (viscous) {
      Breakdown_file << "Free-stream viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << "\n";
      if (turbulent) {
        Breakdown_file << "Free-stream turb. kinetic energy (non-dim): " << config->GetTke_FreeStreamND() << "\n";
        Breakdown_file << "Free-stream specific dissipation (non-dim): " << config->GetOmega_FreeStreamND() << "\n";
      }
    }

    if (unsteady) {
      Breakdown_file << "Total time (non-dim): " << config->GetTotal_UnstTimeND() << "\n";
      Breakdown_file << "Time step (non-dim): " << config->GetDelta_UnstTimeND() << "\n";
    }

    } else {

    /*--- Incompressible version of the console output ---*/

      bool energy     = config->GetEnergy_Equation();
      bool boussinesq = (config->GetKind_DensityModel() == INC_DENSITYMODEL::BOUSSINESQ);

      if (config->GetRef_Inc_NonDim() == DIMENSIONAL) {
        Breakdown_file << "Viscous and Inviscid flow: rho_ref, vel_ref, temp_ref, p_ref" << "\n";
        Breakdown_file << "are set to 1.0 in order to perform a dimensional calculation." << "\n";
        if (dynamic_grid) Breakdown_file << "Force coefficients computed using MACH_MOTION." << "\n";
        else Breakdown_file << "Force coefficients computed using initial values." << "\n";
      }
      else if (config->GetRef_Inc_NonDim() == INITIAL_VALUES) {
        Breakdown_file << "Viscous and Inviscid flow: rho_ref, vel_ref, and temp_ref" << "\n";
        Breakdown_file << "are based on the initial values, p_ref = rho_ref*vel_ref^2." << "\n";
        if (dynamic_grid) Breakdown_file << "Force coefficients computed using MACH_MOTION." << "\n";
        else Breakdown_file << "Force coefficients computed using initial values." << "\n";
      }
      else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
        Breakdown_file << "Viscous and Inviscid flow: rho_ref, vel_ref, and temp_ref" << "\n";
        Breakdown_file << "are user-provided reference values, p_ref = rho_ref*vel_ref^2." << "\n";
        if (dynamic_grid) Breakdown_file << "Force coefficients computed using MACH_MOTION." << "\n";
        else Breakdown_file << "Force coefficients computed using reference values." << "\n";
      }
      Breakdown_file << "The reference area for force coeffs. is " << config->GetRefArea() << " m^2." << "\n";
      Breakdown_file << "The reference length for force coeffs. is " << config->GetRefLength() << " m." << "\n";

      Breakdown_file << "The pressure is decomposed into thermodynamic and dynamic components." << "\n";
      Breakdown_file << "The initial value of the dynamic pressure is 0." << "\n";

      Breakdown_file << "Mach number: "<< config->GetMach();
      if (config->GetKind_FluidModel() == CONSTANT_DENSITY) {
        Breakdown_file << ", computed using the Bulk modulus." << "\n";
      } else {
        Breakdown_file << ", computed using fluid speed of sound." << "\n";
      }

      Breakdown_file << "For external flows, the initial state is imposed at the far-field." << "\n";
      Breakdown_file << "Angle of attack (deg): "<< config->GetAoA() << ", computed using the initial velocity." << "\n";
      Breakdown_file << "Side slip angle (deg): "<< config->GetAoS() << ", computed using the initial velocity." << "\n";

      if (viscous) {
        Breakdown_file << "Reynolds number per meter: " << config->GetReynolds() << ", computed using initial values."<< "\n";
        Breakdown_file << "Reynolds number is a byproduct of inputs only (not used internally)." << "\n";
      }
      Breakdown_file << "SI units only. The grid should be dimensional (meters)." << "\n";

      switch (config->GetKind_DensityModel()) {

        case INC_DENSITYMODEL::CONSTANT:
          if (energy) Breakdown_file << "Energy equation is active and decoupled." << "\n";
          else Breakdown_file << "No energy equation." << "\n";
          break;

        case INC_DENSITYMODEL::BOUSSINESQ:
          if (energy) Breakdown_file << "Energy equation is active and coupled through Boussinesq approx." << "\n";
          break;

        case INC_DENSITYMODEL::VARIABLE:
          if (energy) Breakdown_file << "Energy equation is active and coupled for variable density." << "\n";
          break;

      }

      Breakdown_file <<"-- Input conditions:"<< "\n";

      switch (config->GetKind_FluidModel()) {

        case CONSTANT_DENSITY:
          Breakdown_file << "Fluid Model: CONSTANT_DENSITY "<< "\n";
          if (energy) {
            Breakdown_file << "Specific heat at constant pressure (Cp): " << config->GetSpecific_Heat_Cp() << " N.m/kg.K." << "\n";
          }
          if (boussinesq) Breakdown_file << "Thermal expansion coefficient: " << config->GetThermal_Expansion_Coeff() << " K^-1." << "\n";
          Breakdown_file << "Thermodynamic pressure not required." << "\n";
          break;

        case INC_IDEAL_GAS:
          Breakdown_file << "Fluid Model: INC_IDEAL_GAS "<< endl;
          Breakdown_file << "Variable density incompressible flow using ideal gas law." << endl;
          Breakdown_file << "Density is a function of temperature (constant thermodynamic pressure)." << endl;
          Breakdown_file << "Specific heat at constant pressure (Cp): " << config->GetSpecific_Heat_Cp() << " N.m/kg.K." << endl;
          Breakdown_file << "Molecular weight : "<< config->GetMolecular_Weight() << " g/mol" << endl;
          Breakdown_file << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
          Breakdown_file << "Thermodynamic pressure: " << config->GetPressure_Thermodynamic();
          if (config->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << endl;
          else if (config->GetSystemMeasurements() == US) Breakdown_file << " psf." << endl;
          break;

        case INC_IDEAL_GAS_POLY:
          Breakdown_file << "Fluid Model: INC_IDEAL_GAS_POLY "<< endl;
          Breakdown_file << "Variable density incompressible flow using ideal gas law." << endl;
          Breakdown_file << "Density is a function of temperature (constant thermodynamic pressure)." << endl;
          Breakdown_file << "Molecular weight: " << config->GetMolecular_Weight() << " g/mol." << endl;
          Breakdown_file << "Specific gas constant: " << config->GetGas_Constant() << " N.m/kg.K." << endl;
          Breakdown_file << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << endl;
          Breakdown_file << "Thermodynamic pressure: " << config->GetPressure_Thermodynamic();
          if (config->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << endl;
          else if (config->GetSystemMeasurements() == US) Breakdown_file << " psf." << endl;
          Breakdown_file << "Cp(T) polynomial coefficients: \n  (";
          for (unsigned short iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
            Breakdown_file << config->GetCp_PolyCoeff(iVar);
            if (iVar < config->GetnPolyCoeffs()-1) Breakdown_file << ", ";
          }
          Breakdown_file << ")." << endl;
          Breakdown_file << "Cp(T) polynomial coefficients (non-dim.): \n  (";
          for (unsigned short iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
            Breakdown_file << config->GetCp_PolyCoeffND(iVar);
            if (iVar < config->GetnPolyCoeffs()-1) Breakdown_file << ", ";
          }
          Breakdown_file << ")." << endl;
          break;

      }
      if (viscous) {
        switch (config->GetKind_ViscosityModel()) {

          case VISCOSITYMODEL::CONSTANT:
            Breakdown_file << "Viscosity Model: CONSTANT_VISCOSITY  "<< "\n";
            Breakdown_file << "Constant Laminar Viscosity: " << config->GetMu_Constant();
            if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
            else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
            Breakdown_file << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< "\n";
            break;

          case VISCOSITYMODEL::SUTHERLAND:
            Breakdown_file << "Viscosity Model: SUTHERLAND "<< "\n";
            Breakdown_file << "Ref. Laminar Viscosity: " << config->GetMu_Ref();
            if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
            else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
            Breakdown_file << "Ref. Temperature: " << config->GetMu_Temperature_Ref();
            if (config->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
            else if (config->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
            Breakdown_file << "Sutherland Constant: "<< config->GetMu_S();
            if (config->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
            else if (config->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
            Breakdown_file << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< "\n";
            Breakdown_file << "Ref. Temperature (non-dim): " << config->GetMu_Temperature_RefND()<< "\n";
            Breakdown_file << "Sutherland constant (non-dim): "<< config->GetMu_SND()<< "\n";
            break;

          case VISCOSITYMODEL::POLYNOMIAL:
            Breakdown_file << "Viscosity Model: POLYNOMIAL_VISCOSITY  "<< endl;
            Breakdown_file << "Mu(T) polynomial coefficients: \n  (";
            for (unsigned short iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
              Breakdown_file << config->GetMu_PolyCoeff(iVar);
              if (iVar < config->GetnPolyCoeffs()-1) Breakdown_file << ", ";
            }
            Breakdown_file << ")." << endl;
            Breakdown_file << "Mu(T) polynomial coefficients (non-dim.): \n  (";
            for (unsigned short iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
              Breakdown_file << config->GetMu_PolyCoeffND(iVar);
              if (iVar < config->GetnPolyCoeffs()-1) Breakdown_file << ", ";
            }
            Breakdown_file << ")." << endl;
            break;

          case VISCOSITYMODEL::FLAMELET:
            break;
        }

        if (energy) {
          switch (config->GetKind_ConductivityModel()) {

            case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
              Breakdown_file << "Conductivity Model: CONSTANT_PRANDTL  "<< "\n";
              Breakdown_file << "Prandtl (Laminar): " << config->GetPrandtl_Lam()<< "\n";
              break;

            case CONDUCTIVITYMODEL::CONSTANT:
              Breakdown_file << "Conductivity Model: CONSTANT "<< "\n";
              Breakdown_file << "Molecular Conductivity: " << config->GetThermal_Conductivity_Constant()<< " W/m^2.K." << "\n";
              Breakdown_file << "Molecular Conductivity (non-dim): " << config->GetThermal_Conductivity_ConstantND()<< "\n";
              break;

            case CONDUCTIVITYMODEL::POLYNOMIAL:
              Breakdown_file << "Viscosity Model: POLYNOMIAL "<< endl;
              Breakdown_file << "Kt(T) polynomial coefficients: \n  (";
              for (unsigned short iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
                Breakdown_file << config->GetKt_PolyCoeff(iVar);
                if (iVar < config->GetnPolyCoeffs()-1) Breakdown_file << ", ";
              }
              Breakdown_file << ")." << endl;
              Breakdown_file << "Kt(T) polynomial coefficients (non-dim.): \n  (";
              for (unsigned short iVar = 0; iVar < config->GetnPolyCoeffs(); iVar++) {
                Breakdown_file << config->GetKt_PolyCoeffND(iVar);
                if (iVar < config->GetnPolyCoeffs()-1) Breakdown_file << ", ";
              }
              Breakdown_file << ")." << endl;
              break;

            case CONDUCTIVITYMODEL::FLAMELET:
              break;

          }

          if ((Kind_Solver == RANS) || (Kind_Solver == ADJ_RANS) || (Kind_Solver == DISC_ADJ_RANS)) {
            switch (config->GetKind_ConductivityModel_Turb()) {
              case CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL:
                Breakdown_file << "Turbulent Conductivity Model: CONSTANT_PRANDTL  "<< "\n";
                Breakdown_file << "Turbulent Prandtl: " << config->GetPrandtl_Turb()<< "\n";
                break;
              case CONDUCTIVITYMODEL_TURB::NONE:
                Breakdown_file << "Turbulent Conductivity Model: CONDUCTIVITYMODEL_TURB::NONE "<< "\n";
                Breakdown_file << "No turbulent component in effective thermal conductivity." << "\n";
                break;
            }
          }

        }

      }

      if (config->GetKind_FluidModel() == CONSTANT_DENSITY) {
        Breakdown_file << "Bulk modulus: " << config->GetBulk_Modulus();
        if (config->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
        else if (config->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";
      }

      Breakdown_file << "Initial dynamic pressure: " << config->GetPressure_FreeStream();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";

      Breakdown_file << "Initial total pressure: " << config->GetPressure_FreeStream() + 0.5*config->GetDensity_FreeStream()*config->GetModVel_FreeStream()*config->GetModVel_FreeStream();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";

      if (energy) {
        Breakdown_file << "Initial temperature: " << config->GetTemperature_FreeStream();
        if (config->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
        else if (config->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
      }

      Breakdown_file << "Initial density: " << config->GetDensity_FreeStream();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " kg/m^3." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " slug/ft^3." << "\n";

      if (nDim == 2) {
        Breakdown_file << "Initial velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
        Breakdown_file << config->GetVelocity_FreeStream()[1] << ")";
      }
      if (nDim == 3) {
        Breakdown_file << "Initial velocity: (" << config->GetVelocity_FreeStream()[0] << ", ";
        Breakdown_file << config->GetVelocity_FreeStream()[1] << ", " << config->GetVelocity_FreeStream()[2] << ")";
      }
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " m/s. ";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " ft/s. ";

      Breakdown_file << "Magnitude: "  << config->GetModVel_FreeStream();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " m/s." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " ft/s." << "\n";

      if (viscous) {
        Breakdown_file << "Initial laminar viscosity: " << config->GetViscosity_FreeStream();
        if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
        else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
        if (turbulent) {
          Breakdown_file << "Initial turb. kinetic energy per unit mass: " << config->GetTke_FreeStream();
          if (config->GetSystemMeasurements() == SI) Breakdown_file << " m^2/s^2." << "\n";
          else if (config->GetSystemMeasurements() == US) Breakdown_file << " ft^2/s^2." << "\n";
          Breakdown_file << "Initial specific dissipation: " << config->GetOmega_FreeStream();
          if (config->GetSystemMeasurements() == SI) Breakdown_file << " 1/s." << "\n";
          else if (config->GetSystemMeasurements() == US) Breakdown_file << " 1/s." << "\n";
        }
      }

      if (unsteady) { Breakdown_file << "Total time: " << config->GetTotal_UnstTime() << " s. Time step: " << config->GetDelta_UnstTime() << " s." << "\n"; }

      /*--- Print out reference values. ---*/

      Breakdown_file <<"-- Reference values:"<< "\n";

      if (config->GetKind_FluidModel() != CONSTANT_DENSITY) {
        Breakdown_file << "Reference specific gas constant: " << config->GetGas_Constant_Ref();
        if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.m/kg.K." << "\n";
        else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.ft/slug.R." << "\n";
      } else {
        if (energy) {
          Breakdown_file << "Reference specific heat: " << config->GetGas_Constant_Ref();
          if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.m/kg.K." << "\n";
          else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.ft/slug.R." << "\n";
        }
      }

      Breakdown_file << "Reference pressure: " << config->GetPressure_Ref();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";

      if (energy) {
        Breakdown_file << "Reference temperature: " << config->GetTemperature_Ref();
        if (config->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
        else if (config->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
      }

      Breakdown_file << "Reference density: " << config->GetDensity_Ref();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " kg/m^3." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " slug/ft^3." << "\n";

      Breakdown_file << "Reference velocity: " << config->GetVelocity_Ref();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " m/s." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " ft/s." << "\n";

      Breakdown_file << "Reference length: " << config->GetLength_Ref();
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " m." << "\n";
      else if (config->GetSystemMeasurements() == US) Breakdown_file << " in." << "\n";

      if (viscous) {
        Breakdown_file << "Reference viscosity: " << config->GetViscosity_Ref();
        if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
        else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
      }

      if (unsteady) Breakdown_file << "Reference time: " << config->GetTime_Ref() <<" s." << "\n";

      /*--- Print out resulting non-dim values here. ---*/

      Breakdown_file << "-- Resulting non-dimensional state:" << "\n";
      Breakdown_file << "Mach number (non-dim): " << config->GetMach() << "\n";
      if (viscous) {
        Breakdown_file << "Reynolds number (per m): " << config->GetReynolds() << "\n";
      }

      if (config->GetKind_FluidModel() != CONSTANT_DENSITY) {
        Breakdown_file << "Specific gas constant (non-dim): " << config->GetGas_ConstantND() << "\n";
        Breakdown_file << "Initial thermodynamic pressure (non-dim): " << config->GetPressure_ThermodynamicND() << "\n";
      } else {
        if (energy) {
          Breakdown_file << "Specific heat at constant pressure (non-dim): " << config->GetSpecific_Heat_CpND() << "\n";
          if (boussinesq) Breakdown_file << "Thermal expansion coefficient (non-dim.): " << config->GetThermal_Expansion_CoeffND() << " K^-1." << "\n";
        }
      }

      if (energy) Breakdown_file << "Initial temperature (non-dim): " << config->GetTemperature_FreeStreamND() << "\n";
      Breakdown_file << "Initial pressure (non-dim): " << config->GetPressure_FreeStreamND() << "\n";
      Breakdown_file << "Initial density (non-dim): " << config->GetDensity_FreeStreamND() << "\n";

      if (nDim == 2) {
        Breakdown_file << "Initial velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
        Breakdown_file << config->GetVelocity_FreeStreamND()[1] << "). ";
      } else {
        Breakdown_file << "Initial velocity (non-dim): (" << config->GetVelocity_FreeStreamND()[0] << ", ";
        Breakdown_file << config->GetVelocity_FreeStreamND()[1] << ", " << config->GetVelocity_FreeStreamND()[2] << "). ";
      }
      Breakdown_file << "Magnitude: "   << config->GetModVel_FreeStreamND() << "\n";

      if (viscous) {
        Breakdown_file << "Initial viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << "\n";
        if (turbulent) {
          Breakdown_file << "Initial turb. kinetic energy (non-dim): " << config->GetTke_FreeStreamND() << "\n";
          Breakdown_file << "Initial specific dissipation (non-dim): " << config->GetOmega_FreeStreamND() << "\n";
        }
      }

      if (unsteady) {
        Breakdown_file << "Total time (non-dim): " << config->GetTotal_UnstTimeND() << "\n";
        Breakdown_file << "Time step (non-dim): " << config->GetDelta_UnstTimeND() << "\n";
      }

    }

    /*--- Begin forces breakdown info. ---*/

    Breakdown_file << fixed;
    Breakdown_file << "\n" << "\n" <<"Forces breakdown:" << "\n" << "\n";

    if (nDim == 3) {
      su2double m = solver_container[FLOW_SOL]->GetTotal_CFz()/solver_container[FLOW_SOL]->GetTotal_CFx();
      su2double term = (Total_CoPz/m)-Total_CoPx;

      if (term > 0) Breakdown_file << "Center of Pressure: X="  << 1/m <<"Z-"<< term << "." << "\n\n";
      else Breakdown_file << "Center of Pressure: X="  << 1/m <<"Z+"<< fabs(term);
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " m." << "\n\n";
      else Breakdown_file << " in." << "\n\n";
    }
    else {
      su2double m = solver_container[FLOW_SOL]->GetTotal_CFy()/solver_container[FLOW_SOL]->GetTotal_CFx();
      su2double term = (Total_CoPy/m)-Total_CoPx;
      if (term > 0) Breakdown_file << "Center of Pressure: X="  << 1/m <<"Y-"<< term << "." << "\n\n";
      else Breakdown_file << "Center of Pressure: X="  << 1/m <<"Y+"<< fabs(term);
      if (config->GetSystemMeasurements() == SI) Breakdown_file << " m." << "\n\n";
      else Breakdown_file << " in." << "\n\n";
    }

    /*--- Reference area and force factors. ---*/

    const su2double Factor = solver_container[FLOW_SOL]->GetAeroCoeffsReferenceForce();
    const su2double Ref = config->GetDensity_Ref() * pow(config->GetVelocity_Ref(),2);

    Breakdown_file << "NOTE: Multiply forces by the non-dimensional factor: " << Factor << ", and the reference factor: " << Ref  << "\n";
    Breakdown_file << "to obtain the dimensional force."  << "\n" << "\n";

    Breakdown_file << "Total CL:    ";
    Breakdown_file.width(11);
    Breakdown_file << Total_CL;
    Breakdown_file << " | Pressure (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Inv_CL * 100.0) / (Total_CL + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Inv_CL;
    Breakdown_file << " | Friction (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Visc_CL * 100.0) / (Total_CL + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Visc_CL;
    Breakdown_file << " | Momentum (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Mnt_CL * 100.0) / (Total_CL + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Mnt_CL << "\n";

    Breakdown_file << "Total CD:    ";
    Breakdown_file.width(11);
    Breakdown_file << Total_CD;
    Breakdown_file << " | Pressure (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Inv_CD * 100.0) / (Total_CD + EPS)) << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Inv_CD;
    Breakdown_file << " | Friction (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Visc_CD * 100.0) / (Total_CD + EPS)) << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Visc_CD;
    Breakdown_file << " | Momentum (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Mnt_CD * 100.0) / (Total_CD + EPS)) << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Mnt_CD << "\n";

    if (nDim == 3) {
      Breakdown_file << "Total CSF:   ";
      Breakdown_file.width(11);
      Breakdown_file << Total_CSF;
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Inv_CSF * 100.0) / (Total_CSF + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Inv_CSF;
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file <<  SU2_TYPE::Int((Visc_CSF * 100.0) / (Total_CSF + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Visc_CSF;
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Mnt_CSF * 100.0) / (Total_CSF + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Mnt_CSF << "\n";
    }

    Breakdown_file << "Total CL/CD: ";
    Breakdown_file.width(11);
    Breakdown_file << Total_CEff;
    Breakdown_file << " | Pressure (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Inv_CEff * 100.0) / (Total_CEff + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Inv_CEff;
    Breakdown_file << " | Friction (";
    Breakdown_file.width(5);
    Breakdown_file <<  SU2_TYPE::Int((Visc_CEff * 100.0) / (Total_CEff + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Visc_CEff;
    Breakdown_file << " | Momentum (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Mnt_CEff * 100.0) / (Total_CEff + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Mnt_CEff << "\n";

    if (nDim == 3) {
      Breakdown_file << "Total CMx:   ";
      Breakdown_file.width(11);
      Breakdown_file << Total_CMx;
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Inv_CMx * 100.0) / (Total_CMx + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Inv_CMx;
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Visc_CMx * 100.0) / (Total_CMx + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Visc_CMx;
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Mnt_CMx * 100.0) / (Total_CMx + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Mnt_CMx << "\n";

      Breakdown_file << "Total CMy:   ";
      Breakdown_file.width(11);
      Breakdown_file << Total_CMy;
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Inv_CMy * 100.0) / (Total_CMy + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Inv_CMy;
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Visc_CMy * 100.0) / (Total_CMy + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Visc_CMy;
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Mnt_CMz * 100.0) / (Total_CMz + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Mnt_CMy << "\n";
    }

    Breakdown_file << "Total CMz:   ";
    Breakdown_file.width(11);
    Breakdown_file << Total_CMz;
    Breakdown_file << " | Pressure (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Inv_CMz * 100.0) / (Total_CMz + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Inv_CMz;
    Breakdown_file << " | Friction (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Visc_CMz * 100.0) / (Total_CMz + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Visc_CMz;
    Breakdown_file << " | Momentum (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Mnt_CMz * 100.0) / (Total_CMz + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Mnt_CMz << "\n";

    Breakdown_file << "Total CFx:   ";
    Breakdown_file.width(11);
    Breakdown_file << Total_CFx;
    Breakdown_file << " | Pressure (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Inv_CFx * 100.0) / (Total_CFx + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Inv_CFx;
    Breakdown_file << " | Friction (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Visc_CFx * 100.0) / (Total_CFx + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Visc_CFx;
    Breakdown_file << " | Momentum (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Mnt_CFx * 100.0) / (Total_CFx + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Mnt_CFx << "\n";

    Breakdown_file << "Total CFy:   ";
    Breakdown_file.width(11);
    Breakdown_file << Total_CFy;
    Breakdown_file << " | Pressure (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Inv_CFy * 100.0) / (Total_CFy + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Inv_CFy;
    Breakdown_file << " | Friction (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Visc_CFy * 100.0) / (Total_CFy + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Visc_CFy;
    Breakdown_file << " | Momentum (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Mnt_CFy * 100.0) / (Total_CFy + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Mnt_CFy << "\n";

    if (nDim == 3) {
      Breakdown_file << "Total CFz:   ";
      Breakdown_file.width(11);
      Breakdown_file << Total_CFz;
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Inv_CFz * 100.0) / (Total_CFz + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Inv_CFz;
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Visc_CFz * 100.0) / (Total_CFz + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Visc_CFz;
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Mnt_CFz * 100.0) / (Total_CFz + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Mnt_CFz << "\n";
    }

    Breakdown_file << "\n" << "\n";

    for (iMarker_Monitoring = 0;
         iMarker_Monitoring < config->GetnMarker_Monitoring();
         iMarker_Monitoring++) {

      Breakdown_file << "Surface name: "
      << config->GetMarker_Monitoring_TagBound(
                                                          iMarker_Monitoring) << "\n" << "\n";

      Breakdown_file << "Total CL    (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CL[iMarker_Monitoring] * 100.0)
                       / (Total_CL + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CL[iMarker_Monitoring];
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CL_Inv[iMarker_Monitoring] * 100.0)
                       / (Surface_CL[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CL_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CL_Visc[iMarker_Monitoring] * 100.0)
                       / (Surface_CL[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CL_Visc[iMarker_Monitoring];
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CL_Mnt[iMarker_Monitoring] * 100.0)
                       / (Surface_CL[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CL_Mnt[iMarker_Monitoring] << "\n";

      Breakdown_file << "Total CD    (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CD[iMarker_Monitoring] * 100.0)
                       / (Total_CD + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CD[iMarker_Monitoring];
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CD_Inv[iMarker_Monitoring] * 100.0)
                       / (Surface_CD[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CD_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CD_Visc[iMarker_Monitoring] * 100.0)
                       / (Surface_CD[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CD_Visc[iMarker_Monitoring];
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CD_Mnt[iMarker_Monitoring] * 100.0)
                       / (Surface_CD[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CD_Mnt[iMarker_Monitoring] << "\n";

      if (nDim == 3) {
        Breakdown_file << "Total CSF   (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CSF[iMarker_Monitoring] * 100.0)
                         / (Total_CSF + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CSF[iMarker_Monitoring];
        Breakdown_file << " | Pressure (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CSF_Inv[iMarker_Monitoring] * 100.0)
                         / (Surface_CSF[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CSF_Inv[iMarker_Monitoring];
        Breakdown_file << " | Friction (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CSF_Visc[iMarker_Monitoring] * 100.0)
                         / (Surface_CSF[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CSF_Visc[iMarker_Monitoring];
        Breakdown_file << " | Momentum (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CSF_Mnt[iMarker_Monitoring] * 100.0)
                         / (Surface_CSF[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CSF_Mnt[iMarker_Monitoring] << "\n";
      }

      Breakdown_file << "Total CL/CD (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CEff[iMarker_Monitoring] * 100.0) / (Total_CEff + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CEff[iMarker_Monitoring];
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CEff_Inv[iMarker_Monitoring] * 100.0)
                       / (Surface_CEff[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CEff_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CEff_Visc[iMarker_Monitoring] * 100.0)
                       / (Surface_CEff[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CEff_Visc[iMarker_Monitoring];
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CEff_Mnt[iMarker_Monitoring] * 100.0)
                       / (Surface_CEff[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CEff_Mnt[iMarker_Monitoring] << "\n";

      if (nDim == 3) {

        Breakdown_file << "Total CMx   (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMx[iMarker_Monitoring] * 100.0) / (Total_CMx + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CMx[iMarker_Monitoring];
        Breakdown_file << " | Pressure (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMx_Inv[iMarker_Monitoring] * 100.0)
                         / (Surface_CMx[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CMx_Inv[iMarker_Monitoring];
        Breakdown_file << " | Friction (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMx_Visc[iMarker_Monitoring] * 100.0)
                         / (Surface_CMx[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CMx_Visc[iMarker_Monitoring];
        Breakdown_file << " | Momentum (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMx_Mnt[iMarker_Monitoring] * 100.0)
                         / (Surface_CMx[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CMx_Mnt[iMarker_Monitoring] << "\n";

        Breakdown_file << "Total CMy   (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMy[iMarker_Monitoring] * 100.0) / (Total_CMy + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CMy[iMarker_Monitoring];
        Breakdown_file << " | Pressure (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMy_Inv[iMarker_Monitoring] * 100.0)
                         / (Surface_CMy[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CMy_Inv[iMarker_Monitoring];
        Breakdown_file << " | Friction (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMy_Visc[iMarker_Monitoring] * 100.0)
                         / (Surface_CMy[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CMy_Visc[iMarker_Monitoring];
        Breakdown_file << " | Momentum (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMy_Mnt[iMarker_Monitoring] * 100.0)
                         / (Surface_CMy[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CMy_Mnt[iMarker_Monitoring] << "\n";
      }

      Breakdown_file << "Total CMz   (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int((Surface_CMz[iMarker_Monitoring] * 100.0) / (Total_CMz + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CMz[iMarker_Monitoring];
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CMz_Inv[iMarker_Monitoring] * 100.0)
                       / (Surface_CMz[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CMz_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CMz_Visc[iMarker_Monitoring] * 100.0)
                       / (Surface_CMz[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CMz_Visc[iMarker_Monitoring];
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CMz_Mnt[iMarker_Monitoring] * 100.0)
                       / (Surface_CMz[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CMz_Mnt[iMarker_Monitoring] << "\n";

      Breakdown_file << "Total CFx   (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int((Surface_CFx[iMarker_Monitoring] * 100.0) / (Total_CFx + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CFx[iMarker_Monitoring];
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CFx_Inv[iMarker_Monitoring] * 100.0)
                       / (Surface_CFx[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CFx_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CFx_Visc[iMarker_Monitoring] * 100.0)
                       / (Surface_CFx[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CFx_Visc[iMarker_Monitoring];
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CFx_Mnt[iMarker_Monitoring] * 100.0)
                       / (Surface_CFx[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CFx_Mnt[iMarker_Monitoring] << "\n";

      Breakdown_file << "Total CFy   (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int((Surface_CFy[iMarker_Monitoring] * 100.0) / (Total_CFy + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CFy[iMarker_Monitoring];
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CFy_Inv[iMarker_Monitoring] * 100.0)
                       / (Surface_CFy[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CFy_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CFy_Visc[iMarker_Monitoring] * 100.0)
                       / (Surface_CFy[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CFy_Visc[iMarker_Monitoring];
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CFy_Mnt[iMarker_Monitoring] * 100.0)
                       / (Surface_CFy[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CFy_Mnt[iMarker_Monitoring] << "\n";

      if (nDim == 3) {
        Breakdown_file << "Total CFz   (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CFz[iMarker_Monitoring] * 100.0) / (Total_CFz + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CFz[iMarker_Monitoring];
        Breakdown_file << " | Pressure (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CFz_Inv[iMarker_Monitoring] * 100.0)
                         / (Surface_CFz[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CFz_Inv[iMarker_Monitoring];
        Breakdown_file << " | Friction (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CFz_Visc[iMarker_Monitoring] * 100.0)
                         / (Surface_CFz[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CFz_Visc[iMarker_Monitoring];
        Breakdown_file << " | Momentum (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CFz_Mnt[iMarker_Monitoring] * 100.0)
                         / (Surface_CFz[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CFz_Mnt[iMarker_Monitoring] << "\n";

      }

      Breakdown_file << "\n";


    }

    delete [] Surface_CL;
    delete [] Surface_CD;
    delete [] Surface_CSF;
    delete [] Surface_CEff;
    delete [] Surface_CFx;
    delete [] Surface_CFy;
    delete [] Surface_CFz;
    delete [] Surface_CMx;
    delete [] Surface_CMy;
    delete [] Surface_CMz;

    delete [] Surface_CL_Inv;
    delete [] Surface_CD_Inv;
    delete [] Surface_CSF_Inv;
    delete [] Surface_CEff_Inv;
    delete [] Surface_CFx_Inv;
    delete [] Surface_CFy_Inv;
    delete [] Surface_CFz_Inv;
    delete [] Surface_CMx_Inv;
    delete [] Surface_CMy_Inv;
    delete [] Surface_CMz_Inv;

    delete [] Surface_CL_Visc;
    delete [] Surface_CD_Visc;
    delete [] Surface_CSF_Visc;
    delete [] Surface_CEff_Visc;
    delete [] Surface_CFx_Visc;
    delete [] Surface_CFy_Visc;
    delete [] Surface_CFz_Visc;
    delete [] Surface_CMx_Visc;
    delete [] Surface_CMy_Visc;
    delete [] Surface_CMz_Visc;

    delete [] Surface_CL_Mnt;
    delete [] Surface_CD_Mnt;
    delete [] Surface_CSF_Mnt;
    delete [] Surface_CEff_Mnt;
    delete [] Surface_CFx_Mnt;
    delete [] Surface_CFy_Mnt;
    delete [] Surface_CFz_Mnt;
    delete [] Surface_CMx_Mnt;
    delete [] Surface_CMy_Mnt;
    delete [] Surface_CMz_Mnt;

  }

}

bool CFlowOutput::WriteVolume_Output(CConfig *config, unsigned long Iter, bool force_writing){

  if (config->GetTime_Domain()){
    if (((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) || (config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING)) &&
        ((Iter == 0) || (Iter % config->GetVolume_Wrt_Freq() == 0))){
      return true;
    }

    if ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND) &&
        ((Iter == 0) ||
         (Iter % config->GetVolume_Wrt_Freq() == 0) ||
         ((Iter+1) % config->GetVolume_Wrt_Freq() == 0) || // Restarts need 2 old solution.
         ((Iter+2) == config->GetnTime_Iter()))){ // The last timestep is written anyways but again one needs the step before for restarts.
      return true;
    }
  } else {
    if (config->GetFixed_CL_Mode() && config->GetFinite_Difference_Mode()) return false;
    return ((Iter > 0) && Iter % config->GetVolume_Wrt_Freq() == 0) || force_writing;
  }

  return false || force_writing;
}

void CFlowOutput::SetTimeAveragedFields(){
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

void CFlowOutput::LoadTimeAveragedData(unsigned long iPoint, CVariable *Node_Flow){
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
    SetScreen_Header(config);
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
