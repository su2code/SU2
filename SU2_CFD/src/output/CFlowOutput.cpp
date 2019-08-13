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

CFlowOutput::CFlowOutput(CConfig *config, unsigned short nDim) : COutput (config, nDim){
  
}


CFlowOutput::~CFlowOutput(void){}

void CFlowOutput::AddAnalyzeSurfaceOutput(CConfig *config){
  
  
  /// DESCRIPTION: Average mass flow    
  AddHistoryOutput("AVG_MASSFLOW",             "Avg_Massflow",              FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total average mass flow on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Mach number      
  AddHistoryOutput("AVG_MACH",                 "Avg_Mach",                  FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total average mach number on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Temperature        
  AddHistoryOutput("AVG_TEMP",                 "Avg_Temp",                  FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total average temperature on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Pressure  
  AddHistoryOutput("AVG_PRESS",                "Avg_Press",                 FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total average pressure on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Density  
  AddHistoryOutput("AVG_DENSITY",              "Avg_Density",               FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total average density on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Enthalpy  
  AddHistoryOutput("AVG_ENTHALPY",             "Avg_Enthalpy",              FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total average enthalpy on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Average velocity in normal direction of the surface
  AddHistoryOutput("AVG_NORMALVEL",            "Avg_NormalVel",             FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total average normal velocity on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Flow uniformity 
  AddHistoryOutput("UNIFORMITY",               "Uniformity",                FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total flow uniformity on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Secondary strength
  AddHistoryOutput("SECONDARY_STRENGTH",       "Secondary_Strength",        FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total secondary strength on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Momentum distortion  
  AddHistoryOutput("MOMENTUM_DISTORTION",      "Momentum_Distortion",       FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total momentum distortion on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Secondary over uniformity 
  AddHistoryOutput("SECONDARY_OVER_UNIFORMITY", "Secondary_Over_Uniformity", FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total secondary over uniformity on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Average total temperature  
  AddHistoryOutput("AVG_TOTALTEMP",            "Avg_TotalTemp",             FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total average total temperature all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Average total pressure   
  AddHistoryOutput("AVG_TOTALPRESS",           "Avg_TotalPress",            FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total average total pressure on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// DESCRIPTION: Pressure drop    
  AddHistoryOutput("PRESSURE_DROP",            "Pressure_Drop",             FORMAT_SCIENTIFIC, "FLOW_COEFF", "Total pressure drop on all markers set in MARKER_ANALYZE", TYPE_COEFFICIENT);
  /// END_GROUP
  
  
  /// BEGIN_GROUP: AERO_COEFF_SURF, DESCRIPTION: Surface values on non-solid markers.
  vector<string> Marker_Analyze;
  for (unsigned short iMarker_Analyze = 0; iMarker_Analyze < config->GetnMarker_Analyze(); iMarker_Analyze++){
    Marker_Analyze.push_back(config->GetMarker_Analyze_TagBound(iMarker_Analyze));
  }  
  
  /// DESCRIPTION: Average mass flow    
  AddHistoryOutputPerSurface("AVG_MASSFLOW",             "Avg_Massflow",              FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Mach number      
  AddHistoryOutputPerSurface("AVG_MACH",                 "Avg_Mach",                  FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Temperature        
  AddHistoryOutputPerSurface("AVG_TEMP",                 "Avg_Temp",                  FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Pressure  
  AddHistoryOutputPerSurface("AVG_PRESS",                "Avg_Press",                 FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Density  
  AddHistoryOutputPerSurface("AVG_DENSITY",              "Avg_Density",               FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average Enthalpy  
  AddHistoryOutputPerSurface("AVG_ENTHALPY",             "Avg_Enthalpy",              FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average velocity in normal direction of the surface
  AddHistoryOutputPerSurface("AVG_NORMALVEL",            "Avg_NormalVel",             FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Flow uniformity 
  AddHistoryOutputPerSurface("UNIFORMITY",               "Uniformity",                FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Secondary strength
  AddHistoryOutputPerSurface("SECONDARY_STRENGTH",       "Secondary_Strength",        FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Momentum distortion  
  AddHistoryOutputPerSurface("MOMENTUM_DISTORTION",      "Momentum_Distortion",       FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Secondary over uniformity 
  AddHistoryOutputPerSurface("SECONDARY_OVER_UNIFORMITY", "Secondary_Over_Uniformity", FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average total temperature  
  AddHistoryOutputPerSurface("AVG_TOTALTEMP",            "Avg_TotalTemp",             FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Average total pressure   
  AddHistoryOutputPerSurface("AVG_TOTALPRESS",           "Avg_TotalPress",            FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// DESCRIPTION: Pressure drop    
  AddHistoryOutputPerSurface("PRESSURE_DROP",            "Pressure_Drop",             FORMAT_SCIENTIFIC, "FLOW_COEFF_SURF", Marker_Analyze, TYPE_COEFFICIENT);
  /// END_GROUP
  
}

void CFlowOutput::SetAnalyzeSurface(CSolver *solver, CGeometry *geometry, CConfig *config, bool output){
  
  unsigned short iDim, iMarker, iMarker_Analyze;
  unsigned long iVertex, iPoint;
  su2double Mach = 0.0, Pressure, Temperature = 0.0, TotalPressure = 0.0, TotalTemperature = 0.0,
  Enthalpy, Velocity[3], TangVel[3], Velocity2, MassFlow, Density, Area,
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

          Density = solver->node[iPoint]->GetDensity();
          Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0; Vn = 0.0; Vtang2 = 0.0;

          for (iDim = 0; iDim < nDim; iDim++) {
            Area += (Vector[iDim] * AxiFactor) * (Vector[iDim] * AxiFactor);
            Velocity[iDim] = solver->node[iPoint]->GetVelocity(iDim);
            Velocity2 += Velocity[iDim] * Velocity[iDim];
            Vn += Velocity[iDim] * Vector[iDim] * AxiFactor;
            MassFlow += Vector[iDim] * AxiFactor * Density * Velocity[iDim];
          }
          
          Area       = sqrt (Area);
          if (AxiFactor == 0.0) Vn = 0.0; else Vn /= Area;
          Vn2        = Vn * Vn;
          Pressure   = solver->node[iPoint]->GetPressure();
          SoundSpeed = solver->node[iPoint]->GetSoundSpeed();

          for (iDim = 0; iDim < nDim; iDim++) {
            TangVel[iDim] = Velocity[iDim] - Vn*Vector[iDim]*AxiFactor/Area;
            Vtang2       += TangVel[iDim]*TangVel[iDim];
          }

          if (incompressible){
            if (config->GetKind_DensityModel() == VARIABLE) {
              Mach = sqrt(solver->node[iPoint]->GetVelocity2())/
              sqrt(solver->node[iPoint]->GetSpecificHeatCp()*config->GetPressure_ThermodynamicND()/(solver->node[iPoint]->GetSpecificHeatCv()*solver->node[iPoint]->GetDensity()));
            } else {
              Mach = sqrt(solver->node[iPoint]->GetVelocity2())/
              sqrt(config->GetBulk_Modulus()/(solver->node[iPoint]->GetDensity()));
            }
            Temperature       = solver->node[iPoint]->GetTemperature();
            Enthalpy          = solver->node[iPoint]->GetSpecificHeatCp()*Temperature;
            TotalTemperature  = Temperature + 0.5*Velocity2/solver->node[iPoint]->GetSpecificHeatCp();
            TotalPressure     = Pressure + 0.5*Density*Velocity2;
          }
          else{
            Mach              = sqrt(Velocity2)/SoundSpeed;
            Temperature       = Pressure / (Gas_Constant * Density);
            Enthalpy          = solver->node[iPoint]->GetEnthalpy();
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
  AddHistoryOutput("DRAG",       "CD",   FORMAT_FIXED, "AERO_COEFF", "Total drag coefficient on all surfaces set with MARKER_MONITORING", TYPE_COEFFICIENT);
  /// DESCRIPTION: Lift coefficient 
  AddHistoryOutput("LIFT",       "CL",   FORMAT_FIXED, "AERO_COEFF", "Total lift coefficient on all surfaces set with MARKER_MONITORING", TYPE_COEFFICIENT);
  /// DESCRIPTION: Sideforce coefficient   
  AddHistoryOutput("SIDEFORCE",  "CSF",  FORMAT_FIXED, "AERO_COEFF", "Total sideforce coefficient on all surfaces set with MARKER_MONITORING", TYPE_COEFFICIENT);
  /// DESCRIPTION: Moment around the x-axis    
  AddHistoryOutput("MOMENT-X",   "CMx",  FORMAT_FIXED, "AERO_COEFF", "Total momentum x-component on all surfaces set with MARKER_MONITORING", TYPE_COEFFICIENT);
  /// DESCRIPTION: Moment around the y-axis    
  AddHistoryOutput("MOMENT-Y",   "CMy",  FORMAT_FIXED, "AERO_COEFF", "Total momentum y-component on all surfaces set with MARKER_MONITORING", TYPE_COEFFICIENT);
  /// DESCRIPTION: Moment around the z-axis      
  AddHistoryOutput("MOMENT-Z",   "CMz",  FORMAT_FIXED, "AERO_COEFF", "Total momentum z-component on all surfaces set with MARKER_MONITORING", TYPE_COEFFICIENT);
  /// DESCRIPTION: Force in x direction    
  AddHistoryOutput("FORCE-X",    "CFx",  FORMAT_FIXED, "AERO_COEFF", "Total force x-component on all surfaces set with MARKER_MONITORING", TYPE_COEFFICIENT);
  /// DESCRIPTION: Force in y direction    
  AddHistoryOutput("FORCE-Y",    "CFy",  FORMAT_FIXED, "AERO_COEFF", "Total force y-component on all surfaces set with MARKER_MONITORING", TYPE_COEFFICIENT);
  /// DESCRIPTION: Force in z direction      
  AddHistoryOutput("FORCE-Z",    "CFz",  FORMAT_FIXED, "AERO_COEFF", "Total force z-component on all surfaces set with MARKER_MONITORING", TYPE_COEFFICIENT);
  /// DESCRIPTION: Lift-to-drag ratio
  AddHistoryOutput("EFFICIENCY", "CEff", FORMAT_FIXED, "AERO_COEFF", "Total lift-to-drag ratio on all surfaces set with MARKER_MONITORING", TYPE_COEFFICIENT);
  /// END_GROUP  
  
  /// BEGIN_GROUP: AERO_COEFF_SURF, DESCRIPTION: Aerodynamic coefficients and forces per surface.
  vector<string> Marker_Monitoring;
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++){
    Marker_Monitoring.push_back(config->GetMarker_Monitoring_TagBound(iMarker_Monitoring));
  }  
  /// DESCRIPTION: Drag coefficient   
  AddHistoryOutputPerSurface("DRAG_ON_SURFACE",       "CD",   FORMAT_FIXED, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Lift coefficient   
  AddHistoryOutputPerSurface("LIFT_ON_SURFACE",       "CL",   FORMAT_FIXED, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Sideforce coefficient     
  AddHistoryOutputPerSurface("SIDEFORCE_ON_SURFACE",  "CSF",  FORMAT_FIXED, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Moment around the x-axis      
  AddHistoryOutputPerSurface("MOMENT-X_ON_SURFACE",   "CMx",  FORMAT_FIXED, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Moment around the y-axis      
  AddHistoryOutputPerSurface("MOMENT-Y_ON_SURFACE",   "CMy",  FORMAT_FIXED, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Moment around the z-axis        
  AddHistoryOutputPerSurface("MOMENT-Z_ON_SURFACE",   "CMz",  FORMAT_FIXED, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Force in x direction      
  AddHistoryOutputPerSurface("FORCE-X_ON_SURFACE",    "CFx",  FORMAT_FIXED, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Force in y direction      
  AddHistoryOutputPerSurface("FORCE-Y_ON_SURFACE",    "CFy",  FORMAT_FIXED, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Force in z direction        
  AddHistoryOutputPerSurface("FORCE-Z_ON_SURFACE",    "CFz",  FORMAT_FIXED, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// DESCRIPTION: Lift-to-drag ratio  
  AddHistoryOutputPerSurface("EFFICIENCY_ON_SURFACE", "CEff", FORMAT_FIXED, "AERO_COEFF_SURF", Marker_Monitoring, TYPE_COEFFICIENT);
  /// END_GROUP 

  /// DESCRIPTION: Angle of attack  
  AddHistoryOutput("AOA", "AoA", FORMAT_SCIENTIFIC, "AOA", "Angle of attack");
}

void CFlowOutput::SetAerodynamicCoefficients(CConfig *config, CSolver *flow_solver){
  
  SetHistoryOutputValue("DRAG", flow_solver->GetTotal_CD());
  SetHistoryOutputValue("LIFT", flow_solver->GetTotal_CL());
  if (nDim == 3)
    SetHistoryOutputValue("SIDEFORCE", flow_solver->GetTotal_CSF());
  if (nDim == 3){
    SetHistoryOutputValue("MOMENT-X", flow_solver->GetTotal_CMx());
    SetHistoryOutputValue("MOMENT-Y", flow_solver->GetTotal_CMy());
  }
  SetHistoryOutputValue("MOMENT-Z", flow_solver->GetTotal_CMz());
  SetHistoryOutputValue("FORCE-X", flow_solver->GetTotal_CFx());
  SetHistoryOutputValue("FORCE-Y", flow_solver->GetTotal_CFy());
  if (nDim == 3)
    SetHistoryOutputValue("FORCE-Z", flow_solver->GetTotal_CFz());
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
  
  AddHistoryOutput("CP_DIFF", "Cp_Diff", FORMAT_FIXED, "CP_DIFF", "Cp difference for inverse design");
  
}

void CFlowOutput::Set_CpInverseDesign(CSolver *solver_container, CGeometry *geometry, CConfig *config){
  
  unsigned short iMarker, icommas, Boundary, iDim;
  unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
  su2double XCoord, YCoord, ZCoord, Pressure, PressureCoeff = 0, Cp, CpTarget, *Normal = NULL, Area, PressDiff;
  bool *PointInDomain;
  string text_line, surfCp_filename;
  ifstream Surface_file;
  char cstr[200];
  
  
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
        
#ifndef HAVE_MPI
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
#else
        iPoint = geometry->node[geometry->vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();
#endif
        
        if (geometry->vertex[iMarker][iVertex]->GetNode() < geometry->GetnPointDomain()) {
          Point2Vertex[iPoint][0] = iMarker;
          Point2Vertex[iPoint][1] = iVertex;
          PointInDomain[iPoint] = true;
          solver_container->SetCPressureTarget(iMarker, iVertex, 0.0);
        }
        
      }
    }
  }
  
  /*--- Prepare to read the surface pressure files (CSV) ---*/
  
  surfCp_filename = "TargetCp";
  
  surfCp_filename = config->GetUnsteady_FileName(surfCp_filename, (int)curr_TimeIter, ".dat");
  
  strcpy (cstr, surfCp_filename.c_str());
    
  /*--- Read the surface pressure file ---*/
  
  string::size_type position;
  
  Surface_file.open(cstr, ios::in);
  
  if (!(Surface_file.fail())) {
    
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
        
        solver_container->SetCPressureTarget(iMarker, iVertex, PressureCoeff);
        
      }
      
    }
    
    Surface_file.close();
    
  }
  
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
        
        Cp = solver_container->GetCPressure(iMarker, iVertex);
        CpTarget = solver_container->GetCPressureTarget(iMarker, iVertex);
        
        Area = 0.0;
        for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
          Area += Normal[iDim]*Normal[iDim];
        Area = sqrt(Area);
        
        PressDiff += Area * (CpTarget - Cp) * (CpTarget - Cp);
      }
      
    }
  }
  
#ifdef HAVE_MPI
  su2double MyPressDiff = PressDiff;   PressDiff = 0.0;
  SU2_MPI::Allreduce(&MyPressDiff, &PressDiff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Update the total Cp difference coeffient ---*/
  
  solver_container->SetTotal_CpDiff(PressDiff);
  
  SetHistoryOutputValue("CP_DIFF", PressDiff);
  
  delete [] Point2Vertex;
  delete [] PointInDomain;
}
