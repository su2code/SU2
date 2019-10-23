/*!
 * \file output_flow.cpp
 * \brief Main subroutines for compressible flow output
 * \author R. Sanchez
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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
#include "../../include/output/CFlowOutput.hpp"
#include "../../../Common/include/geometry_structure.hpp"
#include "../../include/solver_structure.hpp"

CFlowOutput::CFlowOutput(CConfig *config, unsigned short nDim, bool fem_output) : COutput (config, nDim, fem_output){
  
}


CFlowOutput::~CFlowOutput(void){}

void CFlowOutput::AddAnalyzeSurfaceOutput(CConfig *config){
  
  
  /// DESCRIPTION: Average mass flow    
  AddHistoryOutput("AVG_MASSFLOW",             "Avg_Massflow",              ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average mass flow on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Mach number      
  AddHistoryOutput("AVG_MACH",                 "Avg_Mach",                  ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average mach number on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Temperature        
  AddHistoryOutput("AVG_TEMP",                 "Avg_Temp",                  ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average temperature on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Pressure  
  AddHistoryOutput("AVG_PRESS",                "Avg_Press",                 ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average pressure on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Density  
  AddHistoryOutput("AVG_DENSITY",              "Avg_Density",               ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average density on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Enthalpy  
  AddHistoryOutput("AVG_ENTHALPY",             "Avg_Enthalpy",              ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average enthalpy on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average velocity in normal direction of the surface
  AddHistoryOutput("AVG_NORMALVEL",            "Avg_NormalVel",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average normal velocity on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Flow uniformity 
  AddHistoryOutput("UNIFORMITY",               "Uniformity",                ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total flow uniformity on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Secondary strength
  AddHistoryOutput("SECONDARY_STRENGTH",       "Secondary_Strength",        ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total secondary strength on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Momentum distortion  
  AddHistoryOutput("MOMENTUM_DISTORTION",      "Momentum_Distortion",       ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total momentum distortion on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Secondary over uniformity 
  AddHistoryOutput("SECONDARY_OVER_UNIFORMITY", "Secondary_Over_Uniformity", ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total secondary over uniformity on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average total temperature  
  AddHistoryOutput("AVG_TOTALTEMP",            "Avg_TotalTemp",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average total temperature all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average total pressure   
  AddHistoryOutput("AVG_TOTALPRESS",           "Avg_TotalPress",            ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total average total pressure on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Pressure drop    
  AddHistoryOutput("PRESSURE_DROP",            "Pressure_Drop",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF", "Total pressure drop on all markers set in MARKER_ANALYZE", HistoryFieldType::COEFFICIENT);
  /// END_GROUP
  
  
  /// BEGIN_GROUP: AERO_COEFF_SURF, DESCRIPTION: Surface values on non-solid markers.
  vector<string> Marker_Analyze;
  for (unsigned short iMarker_Analyze = 0; iMarker_Analyze < config->GetnMarker_Analyze(); iMarker_Analyze++){
    Marker_Analyze.push_back(config->GetMarker_Analyze_TagBound(iMarker_Analyze));
  }  
  
  /// DESCRIPTION: Average mass flow    
  AddHistoryOutputPerSurface("AVG_MASSFLOW",             "Avg_Massflow",              ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Mach number      
  AddHistoryOutputPerSurface("AVG_MACH",                 "Avg_Mach",                  ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Temperature        
  AddHistoryOutputPerSurface("AVG_TEMP",                 "Avg_Temp",                  ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Pressure  
  AddHistoryOutputPerSurface("AVG_PRESS",                "Avg_Press",                 ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Density  
  AddHistoryOutputPerSurface("AVG_DENSITY",              "Avg_Density",               ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average Enthalpy  
  AddHistoryOutputPerSurface("AVG_ENTHALPY",             "Avg_Enthalpy",              ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average velocity in normal direction of the surface
  AddHistoryOutputPerSurface("AVG_NORMALVEL",            "Avg_NormalVel",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Flow uniformity 
  AddHistoryOutputPerSurface("UNIFORMITY",               "Uniformity",                ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Secondary strength
  AddHistoryOutputPerSurface("SECONDARY_STRENGTH",       "Secondary_Strength",        ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Momentum distortion  
  AddHistoryOutputPerSurface("MOMENTUM_DISTORTION",      "Momentum_Distortion",       ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Secondary over uniformity 
  AddHistoryOutputPerSurface("SECONDARY_OVER_UNIFORMITY", "Secondary_Over_Uniformity", ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average total temperature  
  AddHistoryOutputPerSurface("AVG_TOTALTEMP",            "Avg_TotalTemp",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Average total pressure   
  AddHistoryOutputPerSurface("AVG_TOTALPRESS",           "Avg_TotalPress",            ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// DESCRIPTION: Pressure drop    
  AddHistoryOutputPerSurface("PRESSURE_DROP",            "Pressure_Drop",             ScreenOutputFormat::SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, HistoryFieldType::COEFFICIENT);
  /// END_GROUP
  
}

void CFlowOutput::SetAnalyzeSurface(CSolver *solver, CGeometry *geometry, CConfig *config, bool output){
  
  unsigned short iDim, iMarker, iMarker_Analyze;
  unsigned long iVertex, iPoint;
  su2double Mach = 0.0, Pressure, Temperature = 0.0, TotalPressure = 0.0, TotalTemperature = 0.0,
  Enthalpy, Velocity[3] = {}, TangVel[3], Velocity2, MassFlow, Density, Area,
  AxiFactor = 1.0, SoundSpeed, Vn, Vn2, Vtang2, Weight = 1.0;

  su2double Gas_Constant      = config->GetGas_ConstantND();
  su2double Gamma             = config->GetGamma();
  unsigned short nMarker      = config->GetnMarker_All();
  unsigned short nDim         = geometry->GetnDim();
  unsigned short Kind_Average = config->GetKind_Average();

  bool compressible   = config->GetKind_Regime() == COMPRESSIBLE;
  bool incompressible = config->GetKind_Regime() == INCOMPRESSIBLE;
  bool energy         = config->GetEnergy_Equation();


  bool axisymmetric               = config->GetAxisymmetric();
  unsigned short nMarker_Analyze  = config->GetnMarker_Analyze();
  
  su2double  *Vector                    = new su2double[nDim];
  su2double  *Surface_MassFlow          = new su2double[nMarker];
  su2double  *Surface_Mach              = new su2double[nMarker];
  su2double  *Surface_Temperature       = new su2double[nMarker];
  su2double  *Surface_Density           = new su2double[nMarker];
  su2double  *Surface_Enthalpy          = new su2double[nMarker];
  su2double  *Surface_NormalVelocity    = new su2double[nMarker];
  su2double  *Surface_StreamVelocity2   = new su2double[nMarker];
  su2double  *Surface_TransvVelocity2   = new su2double[nMarker];
  su2double  *Surface_Pressure          = new su2double[nMarker];
  su2double  *Surface_TotalTemperature  = new su2double[nMarker];
  su2double  *Surface_TotalPressure     = new su2double[nMarker];
  su2double  *Surface_VelocityIdeal     = new su2double[nMarker];
  su2double  *Surface_Area              = new su2double[nMarker];
  su2double  *Surface_MassFlow_Abs      = new su2double[nMarker];
  
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
  
  /*--- Compute the numerical fan face Mach number, and the total area of the inflow ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    Surface_MassFlow[iMarker]          = 0.0;
    Surface_Mach[iMarker]              = 0.0;
    Surface_Temperature[iMarker]       = 0.0;
    Surface_Density[iMarker]           = 0.0;
    Surface_Enthalpy[iMarker]          = 0.0;
    Surface_NormalVelocity[iMarker]    = 0.0;
    Surface_StreamVelocity2[iMarker]   = 0.0;
    Surface_TransvVelocity2[iMarker]   = 0.0;
    Surface_Pressure[iMarker]          = 0.0;
    Surface_TotalTemperature[iMarker]  = 0.0;
    Surface_TotalPressure[iMarker]     = 0.0;
    Surface_VelocityIdeal[iMarker]     = 0.0;
    Surface_Area[iMarker]              = 0.0;
    Surface_MassFlow_Abs[iMarker]      = 0.0;

    if (config->GetMarker_All_Analyze(iMarker) == YES) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);
          
          if (axisymmetric) {
            if (geometry->node[iPoint]->GetCoord(1) != 0.0)
              AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
            else
              AxiFactor = 1.0;
          } else {
            AxiFactor = 1.0;
          }

          Density = solver->GetNodes()->GetDensity(iPoint);
          Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0; Vn = 0.0; Vtang2 = 0.0;

          for (iDim = 0; iDim < nDim; iDim++) {
            Area += (Vector[iDim] * AxiFactor) * (Vector[iDim] * AxiFactor);
            Velocity[iDim] = solver->GetNodes()->GetVelocity(iPoint,iDim);
            Velocity2 += Velocity[iDim] * Velocity[iDim];
            Vn += Velocity[iDim] * Vector[iDim] * AxiFactor;
            MassFlow += Vector[iDim] * AxiFactor * Density * Velocity[iDim];
          }
          
          Area       = sqrt (Area);
          if (AxiFactor == 0.0) Vn = 0.0; else Vn /= Area;
          Vn2        = Vn * Vn;
          Pressure   = solver->GetNodes()->GetPressure(iPoint);
          SoundSpeed = solver->GetNodes()->GetSoundSpeed(iPoint);

          for (iDim = 0; iDim < nDim; iDim++) {
            TangVel[iDim] = Velocity[iDim] - Vn*Vector[iDim]*AxiFactor/Area;
            Vtang2       += TangVel[iDim]*TangVel[iDim];
          }

          if (incompressible){
            if (config->GetKind_DensityModel() == VARIABLE) {
              Mach = sqrt(solver->GetNodes()->GetVelocity2(iPoint))/
              sqrt(solver->GetNodes()->GetSpecificHeatCp(iPoint)*config->GetPressure_ThermodynamicND()/(solver->GetNodes()->GetSpecificHeatCv(iPoint)*solver->GetNodes()->GetDensity(iPoint)));
            } else {
              Mach = sqrt(solver->GetNodes()->GetVelocity2(iPoint))/
              sqrt(config->GetBulk_Modulus()/(solver->GetNodes()->GetDensity(iPoint)));
            }
            Temperature       = solver->GetNodes()->GetTemperature(iPoint);
            Enthalpy          = solver->GetNodes()->GetSpecificHeatCp(iPoint)*Temperature;
            TotalTemperature  = Temperature + 0.5*Velocity2/solver->GetNodes()->GetSpecificHeatCp(iPoint);
            TotalPressure     = Pressure + 0.5*Density*Velocity2;
          }
          else{
            Mach              = sqrt(Velocity2)/SoundSpeed;
            Temperature       = Pressure / (Gas_Constant * Density);
            Enthalpy          = solver->GetNodes()->GetEnthalpy(iPoint);
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

          /*--- For now, always used the area to weight the uniformities. ---*/

          Weight = abs(Area);

          Surface_StreamVelocity2[iMarker]   += Vn2*Weight;
          Surface_TransvVelocity2[iMarker]   += Vtang2*Weight;

        }
      }
      
    }
    
  }
  
  /*--- Copy to the appropriate structure ---*/
  
  su2double *Surface_MassFlow_Local          = new su2double [nMarker_Analyze];
  su2double *Surface_Mach_Local              = new su2double [nMarker_Analyze];
  su2double *Surface_Temperature_Local       = new su2double [nMarker_Analyze];
  su2double *Surface_Density_Local           = new su2double [nMarker_Analyze];
  su2double *Surface_Enthalpy_Local          = new su2double [nMarker_Analyze];
  su2double *Surface_NormalVelocity_Local    = new su2double [nMarker_Analyze];
  su2double *Surface_StreamVelocity2_Local   = new su2double [nMarker_Analyze];
  su2double *Surface_TransvVelocity2_Local   = new su2double [nMarker_Analyze];
  su2double *Surface_Pressure_Local          = new su2double [nMarker_Analyze];
  su2double *Surface_TotalTemperature_Local  = new su2double [nMarker_Analyze];
  su2double *Surface_TotalPressure_Local     = new su2double [nMarker_Analyze];
  su2double *Surface_Area_Local              = new su2double [nMarker_Analyze];
  su2double *Surface_MassFlow_Abs_Local      = new su2double [nMarker_Analyze];
  
  su2double *Surface_MassFlow_Total          = new su2double [nMarker_Analyze];
  su2double *Surface_Mach_Total              = new su2double [nMarker_Analyze];
  su2double *Surface_Temperature_Total       = new su2double [nMarker_Analyze];
  su2double *Surface_Density_Total           = new su2double [nMarker_Analyze];
  su2double *Surface_Enthalpy_Total          = new su2double [nMarker_Analyze];
  su2double *Surface_NormalVelocity_Total    = new su2double [nMarker_Analyze];
  su2double *Surface_StreamVelocity2_Total   = new su2double [nMarker_Analyze];
  su2double *Surface_TransvVelocity2_Total   = new su2double [nMarker_Analyze];
  su2double *Surface_Pressure_Total          = new su2double [nMarker_Analyze];
  su2double *Surface_TotalTemperature_Total  = new su2double [nMarker_Analyze];
  su2double *Surface_TotalPressure_Total     = new su2double [nMarker_Analyze];
  su2double *Surface_Area_Total              = new su2double [nMarker_Analyze];
  su2double *Surface_MassFlow_Abs_Total      = new su2double [nMarker_Analyze];

  su2double *Surface_MomentumDistortion_Total = new su2double [nMarker_Analyze];

  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
    Surface_MassFlow_Local[iMarker_Analyze]          = 0.0;
    Surface_Mach_Local[iMarker_Analyze]              = 0.0;
    Surface_Temperature_Local[iMarker_Analyze]       = 0.0;
    Surface_Density_Local[iMarker_Analyze]           = 0.0;
    Surface_Enthalpy_Local[iMarker_Analyze]          = 0.0;
    Surface_NormalVelocity_Local[iMarker_Analyze]    = 0.0;
    Surface_StreamVelocity2_Local[iMarker_Analyze]   = 0.0;
    Surface_TransvVelocity2_Local[iMarker_Analyze]   = 0.0;
    Surface_Pressure_Local[iMarker_Analyze]          = 0.0;
    Surface_TotalTemperature_Local[iMarker_Analyze]  = 0.0;
    Surface_TotalPressure_Local[iMarker_Analyze]     = 0.0;
    Surface_Area_Local[iMarker_Analyze]              = 0.0;
    Surface_MassFlow_Abs_Local[iMarker_Analyze]      = 0.0;
    
    Surface_MassFlow_Total[iMarker_Analyze]          = 0.0;
    Surface_Mach_Total[iMarker_Analyze]              = 0.0;
    Surface_Temperature_Total[iMarker_Analyze]       = 0.0;
    Surface_Density_Total[iMarker_Analyze]           = 0.0;
    Surface_Enthalpy_Total[iMarker_Analyze]          = 0.0;
    Surface_NormalVelocity_Total[iMarker_Analyze]    = 0.0;
    Surface_StreamVelocity2_Total[iMarker_Analyze]   = 0.0;
    Surface_TransvVelocity2_Total[iMarker_Analyze]   = 0.0;
    Surface_Pressure_Total[iMarker_Analyze]          = 0.0;
    Surface_TotalTemperature_Total[iMarker_Analyze]  = 0.0;
    Surface_TotalPressure_Total[iMarker_Analyze]     = 0.0;
    Surface_Area_Total[iMarker_Analyze]              = 0.0;
    Surface_MassFlow_Abs_Total[iMarker_Analyze]      = 0.0;

    Surface_MomentumDistortion_Total[iMarker_Analyze] = 0.0;

  }
  
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
        }
        
      }
      
    }
    
  }
  
#ifdef HAVE_MPI
  
  SU2_MPI::Allreduce(Surface_MassFlow_Local, Surface_MassFlow_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_Mach_Local, Surface_Mach_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_Temperature_Local, Surface_Temperature_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_Density_Local, Surface_Density_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_Enthalpy_Local, Surface_Enthalpy_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_NormalVelocity_Local, Surface_NormalVelocity_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_StreamVelocity2_Local, Surface_StreamVelocity2_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_TransvVelocity2_Local, Surface_TransvVelocity2_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_Pressure_Local, Surface_Pressure_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_TotalTemperature_Local, Surface_TotalTemperature_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_TotalPressure_Local, Surface_TotalPressure_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_Area_Local, Surface_Area_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(Surface_MassFlow_Abs_Local, Surface_MassFlow_Abs_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
  
  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
    Surface_MassFlow_Total[iMarker_Analyze]          = Surface_MassFlow_Local[iMarker_Analyze];
    Surface_Mach_Total[iMarker_Analyze]              = Surface_Mach_Local[iMarker_Analyze];
    Surface_Temperature_Total[iMarker_Analyze]       = Surface_Temperature_Local[iMarker_Analyze];
    Surface_Density_Total[iMarker_Analyze]           = Surface_Density_Local[iMarker_Analyze];
    Surface_Enthalpy_Total[iMarker_Analyze]          = Surface_Enthalpy_Local[iMarker_Analyze];
    Surface_NormalVelocity_Total[iMarker_Analyze]    = Surface_NormalVelocity_Local[iMarker_Analyze];
    Surface_StreamVelocity2_Total[iMarker_Analyze]   = Surface_StreamVelocity2_Local[iMarker_Analyze];
    Surface_TransvVelocity2_Total[iMarker_Analyze]   = Surface_TransvVelocity2_Local[iMarker_Analyze];
    Surface_Pressure_Total[iMarker_Analyze]          = Surface_Pressure_Local[iMarker_Analyze];
    Surface_TotalTemperature_Total[iMarker_Analyze]  = Surface_TotalTemperature_Local[iMarker_Analyze];
    Surface_TotalPressure_Total[iMarker_Analyze]     = Surface_TotalPressure_Local[iMarker_Analyze];
    Surface_Area_Total[iMarker_Analyze]              = Surface_Area_Local[iMarker_Analyze];
    Surface_MassFlow_Abs_Total[iMarker_Analyze]      = Surface_MassFlow_Abs_Local[iMarker_Analyze];
  }
  
#endif
  
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
    }

    /*--- Compute flow uniformity parameters separately (always area for now). ---*/

    Area = fabs(Surface_Area_Total[iMarker_Analyze]);

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
    SetHistoryOutputPerSurfaceValue("AVG_MASSFLOW", MassFlow, iMarker_Analyze);
    Tot_Surface_MassFlow += MassFlow;
    
    su2double Mach = Surface_Mach_Total[iMarker_Analyze];
    SetHistoryOutputPerSurfaceValue("AVG_MACH", Mach, iMarker_Analyze);
    Tot_Surface_Mach += Mach;
    
    su2double Temperature = Surface_Temperature_Total[iMarker_Analyze] * config->GetTemperature_Ref();
    SetHistoryOutputPerSurfaceValue("AVG_TEMP", Temperature, iMarker_Analyze);
    Tot_Surface_Temperature += Temperature;
    
    su2double Pressure = Surface_Pressure_Total[iMarker_Analyze] * config->GetPressure_Ref();
    SetHistoryOutputPerSurfaceValue("AVG_PRESS", Pressure, iMarker_Analyze);
    Tot_Surface_Pressure += Pressure;
    
    su2double Density = Surface_Density_Total[iMarker_Analyze] * config->GetDensity_Ref();
    SetHistoryOutputPerSurfaceValue("AVG_DENSITY", Density, iMarker_Analyze);
    Tot_Surface_Density += Density;
    
    su2double Enthalpy = Surface_Enthalpy_Total[iMarker_Analyze];
    SetHistoryOutputPerSurfaceValue("AVG_ENTHALPY", Enthalpy, iMarker_Analyze);
    Tot_Surface_Enthalpy += Enthalpy;
    
    su2double NormalVelocity = Surface_NormalVelocity_Total[iMarker_Analyze] * config->GetVelocity_Ref();
    SetHistoryOutputPerSurfaceValue("AVG_NORMALVEL", NormalVelocity, iMarker_Analyze);
    Tot_Surface_NormalVelocity += NormalVelocity;
    
    su2double Uniformity = sqrt(Surface_StreamVelocity2_Total[iMarker_Analyze]) * config->GetVelocity_Ref();
    SetHistoryOutputPerSurfaceValue("UNIFORMITY", Uniformity, iMarker_Analyze);
    Tot_Surface_StreamVelocity2 += Uniformity;
    
    su2double SecondaryStrength = sqrt(Surface_TransvVelocity2_Total[iMarker_Analyze]) * config->GetVelocity_Ref();
    SetHistoryOutputPerSurfaceValue("SECONDARY_STRENGTH", SecondaryStrength, iMarker_Analyze);
    Tot_Surface_TransvVelocity2 += SecondaryStrength;
    
    su2double MomentumDistortion = Surface_MomentumDistortion_Total[iMarker_Analyze];
    SetHistoryOutputPerSurfaceValue("MOMENTUM_DISTORTION", MomentumDistortion, iMarker_Analyze);
    Tot_Momentum_Distortion += MomentumDistortion;
    
    su2double SecondOverUniform = SecondaryStrength/Uniformity;
    SetHistoryOutputPerSurfaceValue("SECONDARY_OVER_UNIFORMITY", SecondOverUniform, iMarker_Analyze);
    Tot_SecondOverUniformity += SecondOverUniform;
    
    su2double TotalTemperature = Surface_TotalTemperature_Total[iMarker_Analyze] * config->GetTemperature_Ref();
    SetHistoryOutputPerSurfaceValue("AVG_TOTALTEMP", TotalTemperature, iMarker_Analyze);
    Tot_Surface_TotalTemperature += TotalTemperature;
    
    su2double TotalPressure = Surface_TotalPressure_Total[iMarker_Analyze] * config->GetPressure_Ref();
    SetHistoryOutputPerSurfaceValue("AVG_TOTALPRESS", TotalPressure, iMarker_Analyze);
    Tot_Surface_TotalPressure += TotalPressure;
    
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
    SetHistoryOutputPerSurfaceValue("PRESSURE_DROP",  Pressure_Drop, iMarker_Analyze);
    Tot_Surface_PressureDrop += Pressure_Drop;
  }
  
  SetHistoryOutputValue("AVG_MASSFLOW", Tot_Surface_MassFlow);
  SetHistoryOutputValue("AVG_MACH", Tot_Surface_Mach);
  SetHistoryOutputValue("AVG_TEMP", Tot_Surface_Temperature);  
  SetHistoryOutputValue("AVG_PRESS", Tot_Surface_Pressure);
  SetHistoryOutputValue("AVG_DENSITY", Tot_Surface_Density);
  SetHistoryOutputValue("AVG_ENTHALPY", Tot_Surface_Enthalpy);
  SetHistoryOutputValue("AVG_NORMALVEL", Tot_Surface_Enthalpy);
  SetHistoryOutputValue("UNIFORMITY", Tot_Surface_StreamVelocity2);
  SetHistoryOutputValue("SECONDARY_STRENGTH", Tot_Surface_TransvVelocity2);
  SetHistoryOutputValue("MOMENTUM_DISTORTION", Tot_Momentum_Distortion);
  SetHistoryOutputValue("SECONDARY_OVER_UNIFORMITY", Tot_SecondOverUniformity);
  SetHistoryOutputValue("AVG_TOTALTEMP", Tot_Surface_TotalTemperature);
  SetHistoryOutputValue("AVG_TOTALPRESS", Tot_Surface_TotalPressure);
  SetHistoryOutputValue("PRESSURE_DROP",  Tot_Surface_PressureDrop);

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
  
  delete [] Surface_MassFlow_Local;
  delete [] Surface_Mach_Local;
  delete [] Surface_Temperature_Local;
  delete [] Surface_Density_Local;
  delete [] Surface_Enthalpy_Local;
  delete [] Surface_NormalVelocity_Local;
  delete [] Surface_StreamVelocity2_Local;
  delete [] Surface_TransvVelocity2_Local;
  delete [] Surface_Pressure_Local;
  delete [] Surface_TotalTemperature_Local;
  delete [] Surface_TotalPressure_Local;
  delete [] Surface_Area_Local;
  delete [] Surface_MassFlow_Abs_Local;
  
  delete [] Surface_MassFlow_Total;
  delete [] Surface_Mach_Total;
  delete [] Surface_Temperature_Total;
  delete [] Surface_Density_Total;
  delete [] Surface_Enthalpy_Total;
  delete [] Surface_NormalVelocity_Total;
  delete [] Surface_StreamVelocity2_Total;
  delete [] Surface_TransvVelocity2_Total;
  delete [] Surface_Pressure_Total;
  delete [] Surface_TotalTemperature_Total;
  delete [] Surface_TotalPressure_Total;
  delete [] Surface_Area_Total;
  delete [] Surface_MassFlow_Abs_Total;
  delete [] Surface_MomentumDistortion_Total;

  delete [] Surface_MassFlow;
  delete [] Surface_Mach;
  delete [] Surface_Temperature;
  delete [] Surface_Density;
  delete [] Surface_Enthalpy;
  delete [] Surface_NormalVelocity;
  delete [] Surface_StreamVelocity2;
  delete [] Surface_TransvVelocity2;
  delete [] Surface_Pressure;
  delete [] Surface_TotalTemperature;
  delete [] Surface_TotalPressure;
  delete [] Surface_Area;
  delete [] Vector;
  delete [] Surface_VelocityIdeal;
  delete [] Surface_MassFlow_Abs;
  
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


void CFlowOutput::Add_CpInverseDesignOutput(CConfig *config){
  
  AddHistoryOutput("CP_DIFF", "Cp_Diff", ScreenOutputFormat::FIXED, "CP_DIFF", "Cp difference for inverse design");
  
}

void CFlowOutput::Set_CpInverseDesign(CSolver *solver, CGeometry *geometry, CConfig *config){
  
  unsigned short iMarker, icommas, Boundary, iDim;
  unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
  su2double XCoord, YCoord, ZCoord, Pressure, PressureCoeff = 0, Cp, CpTarget, *Normal = NULL, Area, PressDiff = 0.0;
  bool *PointInDomain;
  string text_line, surfCp_filename;
  ifstream Surface_file;
  char cstr[200];
  
  /*--- Prepare to read the surface pressure files (CSV) ---*/
  
  surfCp_filename = "TargetCp";
  
  surfCp_filename = config->GetUnsteady_FileName(surfCp_filename, (int)curTimeIter, ".dat");
  
  strcpy (cstr, surfCp_filename.c_str());
    
  /*--- Read the surface pressure file ---*/
  
  string::size_type position;
  
  Surface_file.open(cstr, ios::in);
  
  if (!(Surface_file.fail())) {
    
    nPointLocal = geometry->GetnPoint();
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
    nPointGlobal = nPointLocal;
#endif
    
    Point2Vertex = new unsigned long[nPointGlobal][2];
    PointInDomain = new bool[nPointGlobal];
    
    for (iPoint = 0; iPoint < nPointGlobal; iPoint ++)
      PointInDomain[iPoint] = false;
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      Boundary   = config->GetMarker_All_KindBC(iMarker);
      
      if ((Boundary == EULER_WALL             ) ||
          (Boundary == HEAT_FLUX              ) ||
          (Boundary == ISOTHERMAL             ) ||
          (Boundary == NEARFIELD_BOUNDARY)) {
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          
          /*--- The Pressure file uses the global numbering ---*/
          
          iPoint = geometry->node[geometry->vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();
          
          if (geometry->vertex[iMarker][iVertex]->GetNode() < geometry->GetnPointDomain()) {
            Point2Vertex[iPoint][0] = iMarker;
            Point2Vertex[iPoint][1] = iVertex;
            PointInDomain[iPoint] = true;
            solver->SetCPressureTarget(iMarker, iVertex, 0.0);
          }
          
        }
      }
    }
    
    getline(Surface_file, text_line);
    
    while (getline(Surface_file, text_line)) {
      for (icommas = 0; icommas < 50; icommas++) {
        position = text_line.find( ",", 0 );
        if (position!=string::npos) text_line.erase (position,1);
      }
      stringstream  point_line(text_line);
      
      if (geometry->GetnDim() == 2) point_line >> iPoint >> XCoord >> YCoord >> Pressure >> PressureCoeff;
      if (geometry->GetnDim() == 3) point_line >> iPoint >> XCoord >> YCoord >> ZCoord >> Pressure >> PressureCoeff;
      
      if (PointInDomain[iPoint]) {
        
        /*--- Find the vertex for the Point and Marker ---*/
        
        iMarker = Point2Vertex[iPoint][0];
        iVertex = Point2Vertex[iPoint][1];
        
        solver->SetCPressureTarget(iMarker, iVertex, PressureCoeff);
        
      }
      
    }
    
    Surface_file.close();
    
    delete [] Point2Vertex;
    delete [] PointInDomain;
    
    /*--- Compute the pressure difference ---*/
    
    PressDiff = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      Boundary   = config->GetMarker_All_KindBC(iMarker);
      
      if ((Boundary == EULER_WALL             ) ||
          (Boundary == HEAT_FLUX              ) ||
          (Boundary == ISOTHERMAL             ) ||
          (Boundary == NEARFIELD_BOUNDARY)) {
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          
          Cp = solver->GetCPressure(iMarker, iVertex);
          CpTarget = solver->GetCPressureTarget(iMarker, iVertex);
          
          Area = 0.0;
          for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
          
          PressDiff += Area * (CpTarget - Cp) * (CpTarget - Cp);
        }
        
      }
    }
    
#ifdef HAVE_MPI
    su2double MyPressDiff = PressDiff;    
    SU2_MPI::Allreduce(&MyPressDiff, &PressDiff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    
  }
  
  /*--- Update the total Cp difference coeffient ---*/
  
  solver->SetTotal_CpDiff(PressDiff);
  
  SetHistoryOutputValue("CP_DIFF", PressDiff);

}

su2double CFlowOutput::GetQ_Criterion(su2double** VelocityGradient) const {

  /*--- Make a 3D copy of the gradient so we do not have worry about nDim ---*/

  su2double Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};

  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
      Grad_Vel[iDim][jDim] = VelocityGradient[iDim][jDim];

  /*--- Q Criterion Eq 1.2 of HALLER, G. (2005). An objective definition of a vortex.
   Journal of Fluid Mechanics, 525, 1-26. doi:10.1017/S0022112004002526 ---*/

  /*--- Components of the strain rate tensor (symmetric) ---*/
  su2double s11 = Grad_Vel[0][0];
  su2double s12 = 0.5 * (Grad_Vel[0][1] + Grad_Vel[1][0]);
  su2double s13 = 0.5 * (Grad_Vel[0][2] + Grad_Vel[2][0]);
  su2double s22 = Grad_Vel[1][1];
  su2double s23 = 0.5 * (Grad_Vel[1][2] + Grad_Vel[2][1]);
  su2double s33 = Grad_Vel[2][2];

  /*--- Components of the spin tensor (skew-symmetric) ---*/
  su2double omega12 = 0.5 * (Grad_Vel[0][1] - Grad_Vel[1][0]);
  su2double omega13 = 0.5 * (Grad_Vel[0][2] - Grad_Vel[2][0]);
  su2double omega23 = 0.5 * (Grad_Vel[1][2] - Grad_Vel[2][1]);

  /*--- Q = ||Omega|| - ||Strain|| ---*/
  su2double Q = 2*(pow(omega12,2) + pow(omega13,2) + pow(omega23,2)) - 
    (pow(s11,2) + pow(s22,2) + pow(s33,2) + 2*(pow(s12,2) + pow(s13,2) + pow(s23,2)));

  return Q;
}

void CFlowOutput::WriteAdditionalFiles(CConfig *config, CGeometry *geometry, CSolver **solver_container){
  
  if (config->GetFixed_CL_Mode() || config->GetFixed_CM_Mode()){
    WriteMetaData(config);
  }
  
  if (config->GetWrt_ForcesBreakdown()){
    WriteForcesBreakdown(config, geometry, solver_container);
  }
  
}

void CFlowOutput::WriteMetaData(CConfig *config){
    
  ofstream meta_file;
  
  string filename = "flow";
  
  filename = config->GetFilename(filename, ".meta", curTimeIter);
  
  /*--- All processors open the file. ---*/

  if (rank == MASTER_NODE) {
    meta_file.open(filename.c_str(), ios::out);
    meta_file.precision(15);
    
    if (config->GetTime_Marching() == DT_STEPPING_1ST || config->GetTime_Marching() == DT_STEPPING_2ND)
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
  
  char cstr[200];
  unsigned short iDim, iMarker_Monitoring;
  ofstream Breakdown_file;
  
  bool compressible       = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool unsteady           = (config->GetTime_Marching() != NO);
  bool viscous            = config->GetViscous();
  bool dynamic_grid       = config->GetDynamic_Grid();
  bool gravity            = config->GetGravityForce();
  bool turbulent          = config->GetKind_Solver() == RANS;
  bool fixed_cl           = config->GetFixed_CL_Mode();
  unsigned short Kind_Solver = config->GetKind_Solver();
  unsigned short Kind_Turb_Model = config->GetKind_Turb_Model();
  unsigned short Ref_NonDim = config->GetRef_NonDim();

  unsigned short nDim =  geometry->GetnDim();
  
  /*--- Output the mean flow solution using only the master node ---*/
  
  if ( rank == MASTER_NODE) {
    
    cout << endl << "Writing the forces breakdown file ("<< config->GetBreakdown_FileName() << ")." << endl;
    
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
    0.0, *Surface_CL = NULL, *Surface_CD = NULL,
    *Surface_CSF = NULL, *Surface_CEff = NULL, *Surface_CFx = NULL,
    *Surface_CFy = NULL, *Surface_CFz = NULL,
    *Surface_CMx = NULL, *Surface_CMy = NULL, *Surface_CMz = NULL,
    *Surface_CL_Inv = NULL,
    *Surface_CD_Inv = NULL, *Surface_CSF_Inv = NULL,
    *Surface_CEff_Inv = NULL, *Surface_CFx_Inv = NULL, *Surface_CFy_Inv =
    NULL, *Surface_CFz_Inv = NULL, *Surface_CMx_Inv = NULL,
    *Surface_CMy_Inv = NULL, *Surface_CMz_Inv = NULL,
    *Surface_CL_Visc = NULL,
    *Surface_CD_Visc = NULL, *Surface_CSF_Visc = NULL,
    *Surface_CEff_Visc = NULL, *Surface_CFx_Visc = NULL, *Surface_CFy_Visc =
    NULL, *Surface_CFz_Visc = NULL, *Surface_CMx_Visc = NULL,
    *Surface_CMy_Visc = NULL, *Surface_CMz_Visc = NULL,
    *Surface_CL_Mnt = NULL,
    *Surface_CD_Mnt = NULL, *Surface_CSF_Mnt = NULL,
    *Surface_CEff_Mnt = NULL, *Surface_CFx_Mnt = NULL, *Surface_CFy_Mnt =
    NULL, *Surface_CFz_Mnt = NULL, *Surface_CMx_Mnt = NULL,
    *Surface_CMy_Mnt = NULL, *Surface_CMz_Mnt = NULL;
    
    /*--- WARNING: when compiling on Windows, ctime() is not available. Comment out
     the two lines below that use the dt variable. ---*/
    //time_t now = time(0);
    //string dt = ctime(&now); dt[24] = '.';
    
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
    
    string filename = config->GetBreakdown_FileName();
    strcpy (cstr, filename.data());
    
    Breakdown_file.open(cstr, ios::out);
    
    Breakdown_file << "\n" <<"-------------------------------------------------------------------------" << "\n";
    Breakdown_file <<"|    ___ _   _ ___                                                      |" << "\n";
    Breakdown_file <<"|   / __| | | |_  )   Release 6.1.0  \"Falcon\"                           |" << "\n";
    Breakdown_file <<"|   \\__ \\ |_| |/ /                                                      |" << "\n";
    Breakdown_file <<"|   |___/\\___//___|   Suite (Computational Fluid Dynamics Code)         |" << "\n";
    Breakdown_file << "|                                                                       |" << "\n";
    //Breakdown_file << "|   Local date and time: " << dt << "                      |" << "\n";
    Breakdown_file <<"-------------------------------------------------------------------------" << "\n";
    Breakdown_file << "| The current SU2 release has been coordinated by the                   |" << "\n";
    Breakdown_file << "| SU2 International Developers Society <www.su2devsociety.org>          |" << "\n";
    Breakdown_file << "| with selected contributions from the open-source community            |" << "\n";
    Breakdown_file <<"-------------------------------------------------------------------------" << "\n";
    Breakdown_file << "| The main research teams contributing to the current release are:      |" << "\n";
    Breakdown_file << "| - Prof. Juan J. Alonso's group at Stanford University.                |" << "\n";
    Breakdown_file << "| - Prof. Piero Colonna's group at Delft University of Technology.      |" << "\n";
    Breakdown_file << "| - Prof. Nicolas R. Gauger's group at Kaiserslautern U. of Technology. |" << "\n";
    Breakdown_file << "| - Prof. Alberto Guardone's group at Polytechnic University of Milan.  |" << "\n";
    Breakdown_file << "| - Prof. Rafael Palacios' group at Imperial College London.            |" << "\n";
    Breakdown_file << "| - Prof. Vincent Terrapon's group at the University of Liege.          |" << "\n";
    Breakdown_file << "| - Prof. Edwin van der Weide's group at the University of Twente.      |" << "\n";
    Breakdown_file << "| - Lab. of New Concepts in Aeronautics at Tech. Inst. of Aeronautics.  |" << "\n";
    Breakdown_file <<"-------------------------------------------------------------------------" << "\n";
    Breakdown_file << "| Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,       |" << "\n";
    Breakdown_file << "|                      Tim Albring, and the SU2 contributors.           |" << "\n";
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
    Breakdown_file <<"-------------------------------------------------------------------------" << "\n";
    
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
            
          case CONSTANT_VISCOSITY:
            Breakdown_file << "Viscosity Model: CONSTANT_VISCOSITY  "<< "\n";
            Breakdown_file << "Laminar Viscosity: " << config->GetMu_Constant();
            if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
            else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
            Breakdown_file << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< "\n";
            break;
            
          case SUTHERLAND:
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
            
        }
        switch (config->GetKind_ConductivityModel()) {
            
          case CONSTANT_PRANDTL:
            Breakdown_file << "Conductivity Model: CONSTANT_PRANDTL  "<< "\n";
            Breakdown_file << "Prandtl: " << config->GetPrandtl_Lam()<< "\n";
            break;
            
          case CONSTANT_CONDUCTIVITY:
            Breakdown_file << "Conductivity Model: CONSTANT_CONDUCTIVITY "<< "\n";
            Breakdown_file << "Molecular Conductivity: " << config->GetKt_Constant()<< " W/m^2.K." << "\n";
            Breakdown_file << "Molecular Conductivity (non-dim): " << config->GetKt_ConstantND()<< "\n";
            break;
            
        }
        
        if ((Kind_Solver == RANS) || (Kind_Solver == INC_RANS)) {
          switch (config->GetKind_ConductivityModel_Turb()) {
            case CONSTANT_PRANDTL_TURB:
              Breakdown_file << "Turbulent Conductivity Model: CONSTANT_PRANDTL_TURB  "<< "\n";
              Breakdown_file << "Turbulent Prandtl: " << config->GetPrandtl_Turb()<< "\n";
              break;
            case NO_CONDUCTIVITY_TURB:
              Breakdown_file << "Turbulent Conductivity Model: NO_CONDUCTIVITY_TURB "<< "\n";
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
      bool boussinesq = (config->GetKind_DensityModel() == BOUSSINESQ);

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

        case CONSTANT:
          if (energy) Breakdown_file << "Energy equation is active and decoupled." << "\n";
          else Breakdown_file << "No energy equation." << "\n";
          break;

        case BOUSSINESQ:
          if (energy) Breakdown_file << "Energy equation is active and coupled through Boussinesq approx." << "\n";
          break;

        case VARIABLE:
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

          case CONSTANT_VISCOSITY:
            Breakdown_file << "Viscosity Model: CONSTANT_VISCOSITY  "<< "\n";
            Breakdown_file << "Constant Laminar Viscosity: " << config->GetMu_Constant();
            if (config->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
            else if (config->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
            Breakdown_file << "Laminar Viscosity (non-dim): " << config->GetMu_ConstantND()<< "\n";
            break;

          case SUTHERLAND:
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
            
          case POLYNOMIAL_VISCOSITY:
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

        }

        if (energy) {
          switch (config->GetKind_ConductivityModel()) {

            case CONSTANT_PRANDTL:
              Breakdown_file << "Conductivity Model: CONSTANT_PRANDTL  "<< "\n";
              Breakdown_file << "Prandtl (Laminar): " << config->GetPrandtl_Lam()<< "\n";
              break;

            case CONSTANT_CONDUCTIVITY:
              Breakdown_file << "Conductivity Model: CONSTANT_CONDUCTIVITY "<< "\n";
              Breakdown_file << "Molecular Conductivity: " << config->GetKt_Constant()<< " W/m^2.K." << "\n";
              Breakdown_file << "Molecular Conductivity (non-dim): " << config->GetKt_ConstantND()<< "\n";
              break;

            case POLYNOMIAL_CONDUCTIVITY:
              Breakdown_file << "Viscosity Model: POLYNOMIAL_CONDUCTIVITY "<< endl;
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
              
          }
          
          if ((Kind_Solver == RANS) || (Kind_Solver == ADJ_RANS) || (Kind_Solver == DISC_ADJ_RANS)) {
            switch (config->GetKind_ConductivityModel_Turb()) {
              case CONSTANT_PRANDTL_TURB:
                Breakdown_file << "Turbulent Conductivity Model: CONSTANT_PRANDTL_TURB  "<< "\n";
                Breakdown_file << "Turbulent Prandtl: " << config->GetPrandtl_Turb()<< "\n";
                break;
              case NO_CONDUCTIVITY_TURB:
                Breakdown_file << "Turbulent Conductivity Model: NO_CONDUCTIVITY_TURB "<< "\n";
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

    su2double RefDensity, RefArea, RefVel, Factor, Ref;
    RefArea     = config->GetRefArea();
    if (compressible) {
      RefDensity  = solver_container[FLOW_SOL]->GetDensity_Inf();
      RefVel = solver_container[FLOW_SOL]->GetModVelocity_Inf();
    } else {
      if ((config->GetRef_Inc_NonDim() == DIMENSIONAL) ||
          (config->GetRef_Inc_NonDim() == INITIAL_VALUES)) {
        RefDensity  = solver_container[FLOW_SOL]->GetDensity_Inf();
        RefVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          RefVel  += solver_container[FLOW_SOL]->GetVelocity_Inf(iDim)*solver_container[FLOW_SOL]->GetVelocity_Inf(iDim);
        RefVel = sqrt(RefVel);
      } else {
        RefDensity = config->GetInc_Density_Ref();
        RefVel    = config->GetInc_Velocity_Ref();
      }
    }
    Factor = (0.5*RefDensity*RefArea*RefVel*RefVel);
    Ref = config->GetDensity_Ref() * config->GetVelocity_Ref() * config->GetVelocity_Ref() * 1.0 * 1.0;

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

    Breakdown_file.close();
    
  }
  
}


bool CFlowOutput::WriteVolume_Output(CConfig *config, unsigned long Iter, bool force_writing){
  
  if (config->GetTime_Domain()){
    if (((config->GetTime_Marching() == DT_STEPPING_1ST) ||
         (config->GetTime_Marching() == TIME_STEPPING)) &&
        ((Iter == 0) || (Iter % config->GetVolume_Wrt_Freq() == 0))){
      return true;
    }
    
    if ((config->GetTime_Marching() == DT_STEPPING_2ND) &&
        ((Iter == 0) || (Iter    % config->GetVolume_Wrt_Freq() == 0) ||
         ((Iter+1) % config->GetVolume_Wrt_Freq() == 0) || 
         ((Iter+2 == config->GetnTime_Iter())))){
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
