/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information
 * \author F. Palacios, T. Economon
 * \version 6.1.0 "Falcon"
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

#include "../include/output_structure.hpp"

COutput::COutput(CConfig *config) {

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  field_width = 12;
  
  ConvergenceTable = new PrintingToolbox::CTablePrinter(&std::cout);
  MultiZoneHeaderTable = new PrintingToolbox::CTablePrinter(&std::cout);

  unsigned short iDim, iZone, iSpan, iMarker;
  
  /*--- Initialize point and connectivity counters to zero. ---*/
  
  nGlobal_Poin      = 0;
  nSurf_Poin        = 0;
  nGlobal_Elem      = 0;
  nSurf_Elem        = 0;
  nGlobal_Tria      = 0;
  nGlobal_Quad      = 0;
  nGlobal_Tetr      = 0;
  nGlobal_Hexa      = 0;
  nGlobal_Pris      = 0;
  nGlobal_Pyra      = 0;
  nGlobal_Line      = 0;
  nGlobal_BoundTria = 0;
  nGlobal_BoundQuad = 0;

  /*--- Initialize pointers to NULL ---*/
  
  Coords = NULL;
  Conn_Line = NULL;     Conn_BoundTria = NULL;  Conn_BoundQuad = NULL;
  Conn_Tria = NULL;     Conn_Quad = NULL;       Conn_Tetr = NULL;
  Conn_Hexa = NULL;     Conn_Pris = NULL;       Conn_Pyra = NULL;
  Data = NULL;
  
  /*--- Initialize parallel pointers to NULL ---*/
  
  nGlobal_Poin_Par    = 0;
  nGlobal_Elem_Par    = 0;
  nGlobal_Surf_Poin   = 0;
  nParallel_Poin      = 0;
  nSurf_Poin_Par      = 0;
  nSurf_Elem_Par      = 0;
  nParallel_Tria      = 0;
  nParallel_Quad      = 0;
  nParallel_Tetr      = 0;
  nParallel_Hexa      = 0;
  nParallel_Pris      = 0;
  nParallel_Pyra      = 0;
  nParallel_Line      = 0;
  nParallel_BoundTria = 0;
  nParallel_BoundQuad = 0;
  
  /*--- Initialize pointers to NULL ---*/
  
  Conn_BoundLine_Par = NULL;  Conn_BoundTria_Par = NULL;  Conn_BoundQuad_Par = NULL;
  Conn_Tria_Par = NULL;  Conn_Quad_Par = NULL;       Conn_Tetr_Par = NULL;
  Conn_Hexa_Par = NULL;  Conn_Pris_Par = NULL;       Conn_Pyra_Par = NULL;
  
  Local_Data         = NULL;
  Local_Data_Copy    = NULL;
  Parallel_Data      = NULL;
  Parallel_Surf_Data = NULL;

  /*--- Initialize structures for storing linear partitioning offsets ---*/

  nGlobalPoint_Sort = 0;
  nLocalPoint_Sort  = 0;
  nPoint_Restart    = 0;

  Local_Halo_Sort = NULL;

  beg_node = NULL;
  end_node = NULL;

  nPoint_Lin = NULL;
  nPoint_Cum = NULL;
  
  /*--- Inlet profile data structures. ---*/

  nRow_InletFile    = NULL;
  nRowCum_InletFile = NULL;
  InletCoords       = NULL;

  Marker_Tags_InletFile.clear();

  /*--- Initialize CGNS write flag ---*/
  
  wrote_base_file = false;
  
  /*--- Initialize CGNS write flag ---*/
  
  wrote_CGNS_base = false;
  
  /*--- Initialize Tecplot surface flag ---*/
  
  wrote_surf_file = false;
  
  /*--- Initialize Paraview write flag ---*/
  
  wrote_Paraview_base = false;
  
  /*--- Initialize residual ---*/
  
  RhoRes_New = EPS;
  RhoRes_Old = new su2double[config->GetnZone()];
  for (iZone = 0; iZone < config->GetnZone(); iZone++) RhoRes_Old[iZone] = EPS;
  
  wrote_Paraview_base = false;

  /*--- Initialize turbo flag ---*/
  turbo = config->GetBoolTurbomachinery();

  if(turbo){
    /*--- Initializate quantities for turboperformace ---*/
    nSpanWiseSections = config->GetnSpanMaxAllZones();
    nMarkerTurboPerf  = config->GetnMarker_TurboPerformance();


    TotalStaticEfficiency         = new su2double*[nMarkerTurboPerf];
    TotalTotalEfficiency          = new su2double*[nMarkerTurboPerf];
    KineticEnergyLoss             = new su2double*[nMarkerTurboPerf];
    TRadius                       = new su2double*[nMarkerTurboPerf];
    TotalPressureLoss             = new su2double*[nMarkerTurboPerf];
    MassFlowIn                    = new su2double*[nMarkerTurboPerf];
    MassFlowOut                   = new su2double*[nMarkerTurboPerf];
    FlowAngleIn                   = new su2double*[nMarkerTurboPerf];
    FlowAngleIn_BC                = new su2double*[nMarkerTurboPerf];
    FlowAngleOut                  = new su2double*[nMarkerTurboPerf];
    EulerianWork                  = new su2double*[nMarkerTurboPerf];
    TotalEnthalpyIn               = new su2double*[nMarkerTurboPerf];
    TotalEnthalpyIn_BC            = new su2double*[nMarkerTurboPerf];
    EntropyIn                     = new su2double*[nMarkerTurboPerf];
    EntropyOut                    = new su2double*[nMarkerTurboPerf];
    EntropyIn_BC                  = new su2double*[nMarkerTurboPerf];
    PressureRatio                 = new su2double*[nMarkerTurboPerf];
    TotalTemperatureIn            = new su2double*[nMarkerTurboPerf];
    EnthalpyOut                   = new su2double*[nMarkerTurboPerf];
    MachIn                        = new su2double**[nMarkerTurboPerf];
    MachOut                       = new su2double**[nMarkerTurboPerf];
    VelocityOutIs                 = new su2double*[nMarkerTurboPerf];
    DensityIn                     = new su2double*[nMarkerTurboPerf];
    PressureIn                    = new su2double*[nMarkerTurboPerf];
    TurboVelocityIn               = new su2double**[nMarkerTurboPerf];
    DensityOut                    = new su2double*[nMarkerTurboPerf];
    PressureOut                   = new su2double*[nMarkerTurboPerf];
    TurboVelocityOut              = new su2double**[nMarkerTurboPerf];
    EnthalpyOutIs                 = new su2double*[nMarkerTurboPerf];
    EntropyGen                    = new su2double*[nMarkerTurboPerf];
    AbsFlowAngleIn                = new su2double*[nMarkerTurboPerf];
    TotalEnthalpyOut              = new su2double*[nMarkerTurboPerf];
    TotalEnthalpyOutIs            = new su2double*[nMarkerTurboPerf];
    RothalpyIn                    = new su2double*[nMarkerTurboPerf];
    RothalpyOut                   = new su2double*[nMarkerTurboPerf];
    AbsFlowAngleOut               = new su2double*[nMarkerTurboPerf];
    PressureOut_BC                = new su2double*[nMarkerTurboPerf];
    TemperatureIn                 = new su2double*[nMarkerTurboPerf];
    TemperatureOut                = new su2double*[nMarkerTurboPerf];
    TotalPressureIn               = new su2double*[nMarkerTurboPerf];
    TotalPressureOut              = new su2double*[nMarkerTurboPerf];
    TotalTemperatureOut           = new su2double*[nMarkerTurboPerf];
    EnthalpyIn                    = new su2double*[nMarkerTurboPerf];
    TurbIntensityIn               = new su2double*[nMarkerTurboPerf];
    Turb2LamViscRatioIn           = new su2double*[nMarkerTurboPerf];
    TurbIntensityOut              = new su2double*[nMarkerTurboPerf];
    Turb2LamViscRatioOut          = new su2double*[nMarkerTurboPerf];
    NuFactorIn                    = new su2double*[nMarkerTurboPerf];
    NuFactorOut                   = new su2double*[nMarkerTurboPerf];

    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
      TotalStaticEfficiency   [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalTotalEfficiency    [iMarker] = new su2double [nSpanWiseSections + 1];
      KineticEnergyLoss       [iMarker] = new su2double [nSpanWiseSections + 1];
      TRadius                 [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalPressureLoss       [iMarker] = new su2double [nSpanWiseSections + 1];
      MassFlowIn              [iMarker] = new su2double [nSpanWiseSections + 1];
      MassFlowOut             [iMarker] = new su2double [nSpanWiseSections + 1];
      FlowAngleIn             [iMarker] = new su2double [nSpanWiseSections + 1];
      FlowAngleIn_BC          [iMarker] = new su2double [nSpanWiseSections + 1];
      FlowAngleOut            [iMarker] = new su2double [nSpanWiseSections + 1];
      EulerianWork            [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalEnthalpyIn         [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalEnthalpyIn_BC      [iMarker] = new su2double [nSpanWiseSections + 1];
      EntropyIn               [iMarker] = new su2double [nSpanWiseSections + 1];
      EntropyOut              [iMarker] = new su2double [nSpanWiseSections + 1];
      EntropyIn_BC            [iMarker] = new su2double [nSpanWiseSections + 1];
      PressureRatio           [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalTemperatureIn      [iMarker] = new su2double [nSpanWiseSections + 1];
      EnthalpyOut             [iMarker] = new su2double [nSpanWiseSections + 1];
      MachIn                  [iMarker] = new su2double*[nSpanWiseSections + 1];
      MachOut                 [iMarker] = new su2double*[nSpanWiseSections + 1];
      VelocityOutIs           [iMarker] = new su2double [nSpanWiseSections + 1];
      DensityIn               [iMarker] = new su2double [nSpanWiseSections + 1];
      PressureIn              [iMarker] = new su2double [nSpanWiseSections + 1];
      TurboVelocityIn         [iMarker] = new su2double*[nSpanWiseSections + 1];
      DensityOut              [iMarker] = new su2double [nSpanWiseSections + 1];
      PressureOut             [iMarker] = new su2double [nSpanWiseSections + 1];
      TurboVelocityOut        [iMarker] = new su2double*[nSpanWiseSections + 1];
      EnthalpyOutIs           [iMarker] = new su2double [nSpanWiseSections + 1];
      EntropyGen              [iMarker] = new su2double [nSpanWiseSections + 1];
      AbsFlowAngleIn          [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalEnthalpyOut        [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalEnthalpyOutIs      [iMarker] = new su2double [nSpanWiseSections + 1];
      RothalpyIn              [iMarker] = new su2double [nSpanWiseSections + 1];
      RothalpyOut             [iMarker] = new su2double [nSpanWiseSections + 1];
      AbsFlowAngleOut         [iMarker] = new su2double [nSpanWiseSections + 1];
      PressureOut_BC          [iMarker] = new su2double [nSpanWiseSections + 1];
      TemperatureIn           [iMarker] = new su2double [nSpanWiseSections + 1];
      TemperatureOut          [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalPressureIn         [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalPressureOut        [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalTemperatureOut     [iMarker] = new su2double [nSpanWiseSections + 1];
      EnthalpyIn              [iMarker] = new su2double [nSpanWiseSections + 1];
      TurbIntensityIn         [iMarker] = new su2double [nSpanWiseSections + 1];
      Turb2LamViscRatioIn     [iMarker] = new su2double [nSpanWiseSections + 1];
      TurbIntensityOut        [iMarker] = new su2double [nSpanWiseSections + 1];
      Turb2LamViscRatioOut    [iMarker] = new su2double [nSpanWiseSections + 1];
      NuFactorIn              [iMarker] = new su2double [nSpanWiseSections + 1];
      NuFactorOut             [iMarker] = new su2double [nSpanWiseSections + 1];


      for (iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++){
        TotalStaticEfficiency   [iMarker][iSpan] = 0.0;
        TotalTotalEfficiency    [iMarker][iSpan] = 0.0;
        KineticEnergyLoss       [iMarker][iSpan] = 0.0;
        TRadius                 [iMarker][iSpan] = 0.0;
        TotalPressureLoss       [iMarker][iSpan] = 0.0;
        MassFlowIn              [iMarker][iSpan] = 0.0;
        MassFlowOut             [iMarker][iSpan] = 0.0;
        FlowAngleIn             [iMarker][iSpan] = 0.0;
        FlowAngleIn_BC          [iMarker][iSpan] = config->GetFlowAngleIn_BC();
        FlowAngleOut            [iMarker][iSpan] = 0.0;
        EulerianWork            [iMarker][iSpan] = 0.0;
        TotalEnthalpyIn         [iMarker][iSpan] = 0.0;
        TotalEnthalpyIn_BC      [iMarker][iSpan] = 0.0;
        EntropyIn               [iMarker][iSpan] = 0.0;
        EntropyOut              [iMarker][iSpan] = 0.0;
        EntropyIn_BC            [iMarker][iSpan] = 0.0;
        PressureRatio           [iMarker][iSpan] = 0.0;
        TotalTemperatureIn      [iMarker][iSpan] = 0.0;
        EnthalpyOut             [iMarker][iSpan] = 0.0;


        VelocityOutIs           [iMarker][iSpan] = 0.0;
        DensityIn               [iMarker][iSpan] = 0.0;
        PressureIn              [iMarker][iSpan] = 0.0;

        DensityOut              [iMarker][iSpan] = 0.0;
        PressureOut             [iMarker][iSpan] = 0.0;

        EnthalpyOutIs           [iMarker][iSpan] = 0.0;
        EntropyGen              [iMarker][iSpan] = 0.0;
        AbsFlowAngleIn          [iMarker][iSpan] = 0.0;
        TotalEnthalpyOut        [iMarker][iSpan] = 0.0;
        TotalEnthalpyOutIs      [iMarker][iSpan] = 0.0;
        RothalpyIn              [iMarker][iSpan] = 0.0;
        RothalpyOut             [iMarker][iSpan] = 0.0;
        AbsFlowAngleOut         [iMarker][iSpan] = 0.0;
        PressureOut_BC          [iMarker][iSpan] = config->GetPressureOut_BC();

        TemperatureIn           [iMarker][iSpan] = 0.0;
        TemperatureOut          [iMarker][iSpan] = 0.0;
        TotalPressureIn         [iMarker][iSpan] = 0.0;
        TotalPressureOut        [iMarker][iSpan] = 0.0;
        TotalTemperatureOut     [iMarker][iSpan] = 0.0;
        EnthalpyIn              [iMarker][iSpan] = 0.0;
        TurbIntensityIn         [iMarker][iSpan] = 0.0;
        Turb2LamViscRatioIn     [iMarker][iSpan] = 0.0;
        TurbIntensityOut        [iMarker][iSpan] = 0.0;
        Turb2LamViscRatioOut    [iMarker][iSpan] = 0.0;
        NuFactorIn              [iMarker][iSpan] = 0.0;
        NuFactorOut             [iMarker][iSpan] = 0.0;
        MachIn                  [iMarker][iSpan] = new su2double[4];
        MachOut                 [iMarker][iSpan] = new su2double[4];
        TurboVelocityIn         [iMarker][iSpan] = new su2double[4];
        TurboVelocityOut        [iMarker][iSpan] = new su2double[4];

        for (iDim = 0; iDim < 4; iDim++){
          MachIn           [iMarker][iSpan][iDim]   = 0.0;
          MachOut          [iMarker][iSpan][iDim]   = 0.0;
          TurboVelocityIn  [iMarker][iSpan][iDim]   = 0.0;
          TurboVelocityOut [iMarker][iSpan][iDim]   = 0.0;
        }
      }
    }
  }
  
  nRequestedHistoryFields = config->GetnHistoryOutput();
  for (unsigned short iField = 0; iField < nRequestedHistoryFields; iField++){
    RequestedHistoryFields.push_back(config->GetHistoryOutput_Field(iField));
  }
  
  nRequestedScreenFields = config->GetnScreenOutput();
  for (unsigned short iField = 0; iField < nRequestedScreenFields; iField++){
    RequestedScreenFields.push_back(config->GetScreenOutput_Field(iField));
  }
  
  nRequestedVolumeFields = config->GetnVolumeOutput();
  for (unsigned short iField = 0; iField < nRequestedVolumeFields; iField++){
    RequestedVolumeFields.push_back(config->GetVolumeOutput_Field(iField));
  }
  
  grid_movement = config->GetGrid_Movement(); 
  
  multizone     = config->GetMultizone_Problem();
  
}

COutput::~COutput(void) {
  /* delete pointers initialized at construction*/
  /* Coords and Conn_*(Connectivity) have their own dealloc functions */
  /* Data is taken care of in DeallocateSolution function */

  if (RhoRes_Old != NULL) delete [] RhoRes_Old;

  /*--- Delete turboperformance pointers initiliazed at constrction  ---*/
  unsigned short iMarker, iSpan;
  if(turbo){
    for(iMarker = 0; iMarker< nMarkerTurboPerf; iMarker++){
      for(iSpan=0; iSpan<nSpanWiseSections+1; iSpan++){
        delete [] MachIn          [iMarker][iSpan];
        delete [] MachOut         [iMarker][iSpan];
        delete [] TurboVelocityIn [iMarker][iSpan];
        delete [] TurboVelocityOut[iMarker][iSpan];
      }
    }
    for(iMarker = 0; iMarker< nMarkerTurboPerf; iMarker++){
      delete [] TotalStaticEfficiency[iMarker];
      delete [] TotalTotalEfficiency [iMarker];
      delete [] KineticEnergyLoss    [iMarker];
      delete [] TRadius              [iMarker];
      delete [] TotalPressureLoss    [iMarker];
      delete [] MassFlowIn           [iMarker];
      delete [] MassFlowOut          [iMarker];
      delete [] FlowAngleIn          [iMarker];
      delete [] FlowAngleOut         [iMarker];
      delete [] EulerianWork         [iMarker];
      delete [] TotalEnthalpyIn      [iMarker];
      delete [] TotalEnthalpyOut     [iMarker];
      delete [] TotalEnthalpyOutIs   [iMarker];
      delete [] PressureRatio        [iMarker];
      delete [] EnthalpyOut          [iMarker];
      delete [] VelocityOutIs        [iMarker];
      delete [] TotalTemperatureIn   [iMarker];
      delete [] FlowAngleIn_BC       [iMarker];
      delete [] EntropyIn            [iMarker];
      delete [] EntropyIn_BC         [iMarker];
      delete [] EntropyOut           [iMarker];
      delete [] TotalEnthalpyIn_BC   [iMarker];
      delete [] DensityIn            [iMarker];
      delete [] PressureIn           [iMarker];
      delete [] DensityOut           [iMarker];
      delete [] PressureOut          [iMarker];
      delete [] EnthalpyOutIs        [iMarker];
      delete [] EntropyGen           [iMarker];
      delete [] AbsFlowAngleIn       [iMarker];
      delete [] RothalpyIn           [iMarker];
      delete [] RothalpyOut          [iMarker];
      delete [] AbsFlowAngleOut      [iMarker];
      delete [] PressureOut_BC       [iMarker];
      delete [] MachIn               [iMarker];
      delete [] MachOut              [iMarker];
      delete [] TurboVelocityIn      [iMarker];
      delete [] TurboVelocityOut     [iMarker];
      delete [] TemperatureIn        [iMarker];
      delete [] TemperatureOut       [iMarker];
      delete [] TotalPressureIn      [iMarker];
      delete [] TotalPressureOut     [iMarker];
      delete [] TotalTemperatureOut  [iMarker];
      delete [] EnthalpyIn           [iMarker];
      delete [] TurbIntensityIn      [iMarker];
      delete [] Turb2LamViscRatioIn  [iMarker];
      delete [] TurbIntensityOut     [iMarker];
      delete [] Turb2LamViscRatioOut [iMarker];
      delete [] NuFactorIn           [iMarker];
      delete [] NuFactorOut          [iMarker];


    }
    delete [] TotalStaticEfficiency;
    delete [] TotalTotalEfficiency;
    delete [] KineticEnergyLoss;
    delete [] TRadius;
    delete [] TotalPressureLoss;
    delete [] MassFlowIn;
    delete [] MassFlowOut;
    delete [] FlowAngleIn;
    delete [] FlowAngleOut;
    delete [] EulerianWork;
    delete [] TotalEnthalpyIn;
    delete [] TotalEnthalpyOut;
    delete [] TotalEnthalpyOutIs;
    delete [] PressureRatio;
    delete [] EnthalpyOut;
    delete [] VelocityOutIs;
    delete [] TotalTemperatureIn;
    delete [] FlowAngleIn_BC;
    delete [] EntropyIn;
    delete [] EntropyIn_BC;
    delete [] EntropyOut;
    delete [] TotalEnthalpyIn_BC;
    delete [] DensityIn;
    delete [] PressureIn;
    delete [] DensityOut;
    delete [] PressureOut;
    delete [] EnthalpyOutIs;
    delete [] EntropyGen;
    delete [] AbsFlowAngleIn;
    delete [] RothalpyIn;
    delete [] RothalpyOut;
    delete [] AbsFlowAngleOut;
    delete [] PressureOut_BC;
    delete [] MachIn;
    delete [] MachOut;
    delete [] TurboVelocityIn;
    delete [] TurboVelocityOut;
    delete [] TemperatureIn;
    delete [] TemperatureOut;
    delete [] TotalPressureIn;
    delete [] TotalPressureOut;
    delete [] TotalTemperatureOut;
    delete [] EnthalpyIn;
    delete [] TurbIntensityIn;
    delete [] Turb2LamViscRatioIn;
    delete [] TurbIntensityOut;
    delete [] Turb2LamViscRatioOut;
    delete [] NuFactorIn;
    delete [] NuFactorOut;
  }
}


void COutput::MergeConnectivity(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
  /*--- Flags identifying the types of files to be written. ---*/
  
  bool Wrt_Vol = config->GetWrt_Vol_Sol();
  bool Wrt_Srf = config->GetWrt_Srf_Sol();
  
  /*--- Merge connectivity for each type of element (excluding halos). Note
   that we only need to merge the connectivity once, as it does not change
   during computation. Check whether the base file has been written. ---*/
  
  /*--- Merge volumetric grid. ---*/
  
  if (Wrt_Vol) {
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tria != 0))
      cout <<"Merging volumetric triangle grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, TRIANGLE    );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Quad != 0))
      cout <<"Merging volumetric quadrilateral grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, QUADRILATERAL   );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tetr != 0))
      cout <<"Merging volumetric tetrahedron grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, TETRAHEDRON );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Hexa != 0))
      cout <<"Merging volumetric hexahedron grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, HEXAHEDRON  );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Pris != 0))
      cout <<"Merging volumetric prism grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, PRISM       );
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Pyra != 0))
      cout <<"Merging volumetric pyramid grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, PYRAMID     );
    
  }
  
  /*--- Merge surface grid. ---*/
  
  if (Wrt_Srf) {
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Line != 0))
      cout <<"Merging surface line grid connectivity." << endl;
    MergeSurfaceConnectivity(config, geometry, LINE);
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundTria != 0))
      cout <<"Merging surface triangle grid connectivity." << endl;
    MergeSurfaceConnectivity(config, geometry, TRIANGLE);
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundQuad != 0))
      cout <<"Merging surface quadrilateral grid connectivity." << endl;
    MergeSurfaceConnectivity(config, geometry, QUADRILATERAL);
    
  }
  
  /*--- Update total number of volume elements after merge. ---*/
  
  nGlobal_Elem = nGlobal_Tria + nGlobal_Quad + nGlobal_Tetr +
  nGlobal_Hexa + nGlobal_Pyra + nGlobal_Pris;
  
  /*--- Update total number of surface elements after merge. ---*/
  
  nSurf_Elem = nGlobal_Line + nGlobal_BoundTria + nGlobal_BoundQuad;
  
}

void COutput::MergeCoordinates(CConfig *config, CGeometry *geometry) {
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint;
  
  unsigned short kind_SU2 = config->GetKind_SU2();
  
#ifndef HAVE_MPI
  
  /*--- In serial, the single process has access to all geometry, so simply
   load the coordinates into the data structure. ---*/
  
  unsigned short iMarker;
  unsigned long iVertex, nTotalPoints = 0;
  int SendRecv;
  
  bool isPeriodic;
  
  /*--- First, create a structure to locate any periodic halo nodes ---*/
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- For SU2_CFD and SU2_SOL we want to remove the periodic halo nodes,
         * but for SU2_DEF we want them to be included, therefore the definition of a periodic point
         * is different in each case ---*/

        if (kind_SU2 == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }

        if (isPeriodic && (SendRecv < 0)) {
          Local_Halo[iPoint] = false;
        }
      }
      
    }
  }
  
  /*--- Total number of points in the mesh (this might include periodic points). ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (!Local_Halo[iPoint]) nTotalPoints++;
  
  nGlobal_Poin = nTotalPoints;
  nGlobal_Doma = geometry->GetnPointDomain();
  
  /*--- Allocate the coordinates data structure. ---*/
  
  Coords = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Coords[iDim] = new su2double[nGlobal_Poin];
  }
  
  /*--- Loop over the mesh to collect the coords of the local points ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node).
     Sort by the global index, even in serial there is a renumbering (e.g. RCM). ---*/
    
    if (!Local_Halo[iPoint]) {
      
      /*--- Retrieve the current coordinates at this node. ---*/
      
      unsigned long iGlobal_Index = geometry->node[iPoint]->GetGlobalIndex();
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Coords[iDim][iGlobal_Index] = geometry->node[iPoint]->GetCoord(iDim);
        
        /*--- If US system, the output should be in inches ---*/
        
        if ((config->GetSystemMeasurements() == US) && (config->GetKind_SU2() != SU2_DEF)) {
          Coords[iDim][iGlobal_Index] *= 12.0;
        }
        
      }
      
    }
  }
  
  
  delete [] Local_Halo;
  
#else
  
  /*--- MPI preprocessing ---*/
  int iProcessor, nProcessor = size;
  unsigned long jPoint;

  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;
  
  /*--- Local variables needed for merging the geometry with MPI. ---*/
  
  unsigned long iVertex, iMarker;
  unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/
  
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        
        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          /*--- For SU2_CFD and SU2_SOL we want to remove the periodic halo nodes,
           * but for SU2_DEF we want them to be included, therefore the definition of a periodic point
           * is different in each case ---*/

          if (kind_SU2 == SU2_DEF) {
            isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
          }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          }
          if (isPeriodic) {
            Local_Halo[iPoint] = false;
          }
        }
      }
    }
    
    /*--- Sum total number of nodes that belong to the domain ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;
  }
  Buffer_Send_nPoin[0] = nLocalPoint;
  
  /*--- Communicate the total number of nodes on this domain. ---*/
  
  SU2_MPI::Gather(&Buffer_Send_nPoin, 1, MPI_UNSIGNED_LONG,
                  Buffer_Recv_nPoin, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  
  if (rank == MASTER_NODE) {
    nGlobal_Doma = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Doma += Buffer_Recv_nPoin[iProcessor];
    }
  }
  nBuffer_Scalar = MaxLocalPoint;
  
  /*--- Send and Recv buffers. ---*/
  
  su2double *Buffer_Send_X = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_X = NULL;
  
  su2double *Buffer_Send_Y = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Y = NULL;
  
  su2double *Buffer_Send_Z = NULL, *Buffer_Recv_Z = NULL;
  if (nDim == 3) Buffer_Send_Z = new su2double[MaxLocalPoint];
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_X = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Y = new su2double[nProcessor*MaxLocalPoint];
    if (nDim == 3) Buffer_Recv_Z = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoin[iProcessor];
    }
    Coords = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Coords[iDim] = new su2double[nGlobal_Poin];
    }
  }
  
  /*--- Main communication routine. Loop over each coordinate and perform
   the MPI comm. Temporary 1-D buffers are used to send the coordinates at
   all nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/
  
  /*--- Loop over this partition to collect the coords of the local points. ---*/
  su2double *Coords_Local; jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check for halos and write only if requested ---*/
    if (!Local_Halo[iPoint] || Wrt_Halo) {
      
      /*--- Retrieve local coordinates at this node. ---*/
      Coords_Local = geometry->node[iPoint]->GetCoord();
      
      /*--- Load local coords into the temporary send buffer. ---*/
      Buffer_Send_X[jPoint] = Coords_Local[0];
      Buffer_Send_Y[jPoint] = Coords_Local[1];
      if (nDim == 3) Buffer_Send_Z[jPoint] = Coords_Local[2];
      
      /*--- If US system, the output should be in inches ---*/
      
      if ((config->GetSystemMeasurements() == US) && (config->GetKind_SU2() != SU2_DEF)) {
        Buffer_Send_X[jPoint] *= 12.0;
        Buffer_Send_Y[jPoint] *= 12.0;
        if (nDim == 3) Buffer_Send_Z[jPoint] *= 12.0;
      }
      
      /*--- Store the global index for this local node. ---*/
      Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
      
      /*--- Increment jPoint as the counter. We need this because iPoint
       may include halo nodes that we skip over during this loop. ---*/
      jPoint++;
    }
  }
  
  /*--- Gather the coordinate data on the master node using MPI. ---*/
  
  SU2_MPI::Gather(Buffer_Send_X, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_X, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) {
    SU2_MPI::Gather(Buffer_Send_Z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  }
  SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  /*--- The master node unpacks and sorts this variable by global index ---*/
  
  if (rank == MASTER_NODE) {
    jPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        /*--- Get global index, then loop over each variable and store ---*/
        iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
        if (iGlobal_Index >= nGlobal_Poin) {
          cout << iGlobal_Index << " " << nGlobal_Poin << endl;
        }
        Coords[0][iGlobal_Index] = Buffer_Recv_X[jPoint];
        Coords[1][iGlobal_Index] = Buffer_Recv_Y[jPoint];
        if (nDim == 3) Coords[2][iGlobal_Index] = Buffer_Recv_Z[jPoint];
        jPoint++;
      }
      /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
      jPoint = (iProcessor+1)*nBuffer_Scalar;
    }
  }
  
  /*--- Immediately release the temporary data buffers. ---*/
  
  delete [] Local_Halo;
  delete [] Buffer_Send_X;
  delete [] Buffer_Send_Y;
  if (Buffer_Send_Z != NULL) delete [] Buffer_Send_Z;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_X;
    delete [] Buffer_Recv_Y;
    if (Buffer_Recv_Z != NULL)  delete [] Buffer_Recv_Z;
    delete [] Buffer_Recv_GlobalIndex;
    delete [] Buffer_Recv_nPoin;
  }
  
#endif
  
}

void COutput::MergeVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  int iProcessor;
  unsigned short NODES_PER_ELEMENT;
  unsigned long iPoint, iNode, jNode;
  unsigned long iElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;
  
  unsigned long iVertex, iMarker;
  unsigned long jElem;
  int SendRecv, RecvFrom;
  
  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0;
  unsigned long MaxLocalElem = 0, iGlobal_Index, jPoint, kPoint;
  
  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem = NULL, notPeriodic, notHalo, addedPeriodic, isPeriodic;
  
  unsigned short kind_SU2 = config->GetKind_SU2();
  
  int *Conn_Elem = NULL;
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  switch (Elem_Type) {
    case TRIANGLE:
      nLocalElem = geometry->GetnElemTria();
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      nLocalElem = geometry->GetnElemQuad();
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      nLocalElem = geometry->GetnElemTetr();
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      nLocalElem = geometry->GetnElemHexa();
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case PRISM:
      nLocalElem = geometry->GetnElemPris();
      NODES_PER_ELEMENT = N_POINTS_PRISM;
      break;
    case PYRAMID:
      nLocalElem = geometry->GetnElemPyra();
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      nLocalElem = 0;
      NODES_PER_ELEMENT = 0;
      break;
  }
  
  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/
  
  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[size];
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalElem, &MaxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalElem = nLocalElem;
  Buffer_Recv_nElem[0] = Buffer_Send_nElem[0];
#endif
  
  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;
  
  /*--- Send and Recv buffers ---*/
  
  unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
  unsigned long *Buffer_Recv_Elem = NULL;
  
  unsigned short *Buffer_Send_Halo = new unsigned short[MaxLocalElem];
  unsigned short *Buffer_Recv_Halo = NULL;
  
  /*--- Prepare the receive buffers on the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new unsigned long[size*nBuffer_Scalar];
    Buffer_Recv_Halo = new unsigned short[size*MaxLocalElem];
    if (MaxLocalElem > 0) Conn_Elem = new int[size*MaxLocalElem*NODES_PER_ELEMENT];
  }
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that we part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  
  if (kind_SU2 != SU2_DEF) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
              (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
              (SendRecv < 0)) {
            Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
          }
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found ---*/
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];
  
  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++) Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iGlobal_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain.
         For SU2_DEF we want to keep all periodic nodes. ---*/
        
        if (kind_SU2 == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }
        
        notPeriodic = (isPeriodic && (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (iGlobal_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          kPoint = (iProcessor+1)*maxAddedPeriodic;
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
      }
    }
  }
  
  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/
  
  jNode = 0; jElem = 0;
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if (geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {
      
      /*--- Loop over all nodes in this element and load the
       connectivity into the send buffer. ---*/
      
      Buffer_Send_Halo[jElem] = false;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        
        /*--- Store the global index values directly. ---*/
        
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        Buffer_Send_Elem[jNode] = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- Check if this is a halo node. If so, flag this element
         as a halo cell. We will use this later to sort and remove
         any duplicates from the connectivity list. ---*/
        
        if (Local_Halo[iPoint]) {
          Buffer_Send_Halo[jElem] = true;
        }
        
        /*--- Increment jNode as the counter. We need this because iElem
         may include other elements that we skip over during this loop. ---*/
        
        jNode++;
      }
      jElem++;
    }
  }
  
  /*--- Gather the element connectivity information. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, Buffer_Recv_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Elem[iPoint] = Buffer_Send_Elem[iPoint];
  for (iPoint = 0; iPoint < MaxLocalElem; iPoint++) Buffer_Recv_Halo[iPoint] = Buffer_Send_Halo[iPoint];
#endif
  
  /*--- The master node unpacks and sorts the connectivity. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/
    
    Write_Elem = new bool[size*MaxLocalElem];
    for (iElem = 0; iElem < size*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }
    
    /*--- Remove the rind layer from the solution only if requested ---*/
    
    if (!Wrt_Halo) {
      
      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/
      
      kElem = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
          
          /*--- Check if this element was marked as a halo. ---*/
          if (Buffer_Recv_Halo[kElem+iElem])
            Write_Elem[kElem+iElem] = false;
          
        }
        kElem = (iProcessor+1)*MaxLocalElem;
      }
    }
    
    /*--- Store the unique connectivity list for this element type. ---*/
    
    jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {
          
          /*--- Increment total count for this element type ---*/
          nElem_Total++;
          
          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/
          
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = (int)Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
            kNode++;
          }
        }
      }
      /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
      jElem = (iProcessor+1)*MaxLocalElem;
      jNode = (iProcessor+1)*nBuffer_Scalar;
    }
  }
  
  /*--- Immediately release the temporary buffers. ---*/
  delete [] Buffer_Send_Elem;
  delete [] Buffer_Send_Halo;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] Local_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }
  
  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case TRIANGLE:
        nGlobal_Tria = nElem_Total;
        if (nGlobal_Tria > 0) Conn_Tria = Conn_Elem;
        break;
      case QUADRILATERAL:
        nGlobal_Quad = nElem_Total;
        if (nGlobal_Quad > 0) Conn_Quad = Conn_Elem;
        break;
      case TETRAHEDRON:
        nGlobal_Tetr = nElem_Total;
        if (nGlobal_Tetr > 0) Conn_Tetr = Conn_Elem;
        break;
      case HEXAHEDRON:
        nGlobal_Hexa = nElem_Total;
        if (nGlobal_Hexa > 0) Conn_Hexa = Conn_Elem;
        break;
      case PRISM:
        nGlobal_Pris = nElem_Total;
        if (nGlobal_Pris > 0) Conn_Pris = Conn_Elem;
        break;
      case PYRAMID:
        nGlobal_Pyra = nElem_Total;
        if (nGlobal_Pyra > 0) Conn_Pyra = Conn_Elem;
        break;
      default:
        SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
        break;
    }
  }
  
}

void COutput::MergeSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  unsigned short NODES_PER_ELEMENT;
  
  unsigned short iMarker;
  unsigned long iPoint, iNode, jNode;
  unsigned long iElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;
  
  int iProcessor;
  unsigned long jElem;
  
  unsigned long iVertex;
  
  int SendRecv, RecvFrom;
  
  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0;
  unsigned long MaxLocalElem = 0, iGlobal_Index, jPoint, kPoint;
  
  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem = NULL, notPeriodic, notHalo, addedPeriodic;
  
  
  int *Conn_Elem = NULL;
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  nLocalElem = 0;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          nLocalElem++;
        }
      }
    }
  }
  
  switch (Elem_Type) {
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      NODES_PER_ELEMENT = 0;
      break;
  }
  
  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/
  
  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[size];
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalElem, &MaxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalElem = nLocalElem;
  Buffer_Recv_nElem[0] = Buffer_Send_nElem[0];
#endif
  
  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;
  
  /*--- Send and Recv buffers ---*/
  
  unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
  unsigned long *Buffer_Recv_Elem = NULL;
  
  unsigned short *Buffer_Send_Halo = new unsigned short[MaxLocalElem];
  unsigned short *Buffer_Recv_Halo = NULL;
  
  /*--- Prepare the receive buffers on the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new unsigned long[size*nBuffer_Scalar];
    Buffer_Recv_Halo = new unsigned short[size*MaxLocalElem];
    if (MaxLocalElem > 0) Conn_Elem = new int[size*MaxLocalElem*NODES_PER_ELEMENT];
  }
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that we part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
            (SendRecv < 0)) {
          Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found ---*/
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];
  
  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++) Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iGlobal_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain ---*/
        notPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                       (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
                       (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (iGlobal_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          kPoint = (iProcessor+1)*maxAddedPeriodic;
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
      }
    }
  }
  
  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/
  jNode = 0; jElem = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        
        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          
          /*--- Loop over all nodes in this element and load the
           connectivity into the send buffer. ---*/
          
          Buffer_Send_Halo[jElem] = false;
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            
            /*--- Store the global index values directly. ---*/
            
            iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
            Buffer_Send_Elem[jNode] = geometry->node[iPoint]->GetGlobalIndex();
            
            /*--- Check if this is a halo node. If so, flag this element
             as a halo cell. We will use this later to sort and remove
             any duplicates from the connectivity list. ---*/
            
            if (Local_Halo[iPoint])
              Buffer_Send_Halo[jElem] = true;
            
            /*--- Increment jNode as the counter. We need this because iElem
             may include other elements that we skip over during this loop. ---*/
            
            jNode++;
          }
          jElem++;
        }
      }
  
  /*--- Gather the element connectivity information. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, Buffer_Recv_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Elem[iPoint] = Buffer_Send_Elem[iPoint];
  for (iPoint = 0; iPoint < MaxLocalElem; iPoint++) Buffer_Recv_Halo[iPoint] = Buffer_Send_Halo[iPoint];
#endif
  
  /*--- The master node unpacks and sorts the connectivity. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/
    
    Write_Elem = new bool[size*MaxLocalElem];
    for (iElem = 0; iElem < size*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }
    
    /*--- Remove the rind layer from the solution only if requested ---*/
    
    if (!Wrt_Halo) {
      
      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/
      
      kElem = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
          
          /*--- Check if this element was marked as a halo. ---*/
          if (Buffer_Recv_Halo[kElem+iElem])
            Write_Elem[kElem+iElem] = false;
          
        }
        kElem = (iProcessor+1)*MaxLocalElem;
      }
    }
    
    /*--- Store the unique connectivity list for this element type. ---*/
    
    jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {
        
        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {
          
          /*--- Increment total count for this element type ---*/
          nElem_Total++;
          
          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/
          
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = (int)Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
            kNode++;
          }
        }
      }
      /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
      jElem = (iProcessor+1)*MaxLocalElem;
      jNode = (iProcessor+1)*nBuffer_Scalar;
    }
  }
  
  /*--- Immediately release the temporary buffers. ---*/
  delete [] Buffer_Send_Elem;
  delete [] Buffer_Send_Halo;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] Local_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }
  
  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case LINE:
        nGlobal_Line = nElem_Total;
        if (nGlobal_Line > 0) Conn_Line = Conn_Elem;
        break;
      case TRIANGLE:
        nGlobal_BoundTria = nElem_Total;
        if (nGlobal_BoundTria > 0) Conn_BoundTria = Conn_Elem;
        break;
      case QUADRILATERAL:
        nGlobal_BoundQuad = nElem_Total;
        if (nGlobal_BoundQuad > 0) Conn_BoundQuad = Conn_Elem;
        break;
      default:
        SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
        break;
    }
  }
  
}

void COutput::MergeSolution(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar = 0, jVar = 0, FirstIndex = NONE, SecondIndex = NONE, ThirdIndex = NONE;
  unsigned short nVar_First = 0, nVar_Second = 0, nVar_Third = 0;
  unsigned short iVar_GridVel = 0, iVar_PressCp = 0, iVar_Lam = 0, iVar_MachMean = 0,
  iVar_ViscCoeffs = 0, iVar_HeatCoeffs = 0, iVar_Sens = 0, iVar_Extra = 0, iVar_Eddy = 0, iVar_Sharp = 0,
  iVar_FEA_Vel = 0, iVar_FEA_Accel = 0, iVar_FEA_Stress = 0, iVar_FEA_Stress_3D = 0,
  iVar_FEA_Extra = 0, iVar_SensDim = 0;
  unsigned long iPoint = 0, jPoint = 0, iVertex = 0, iMarker = 0;
  su2double Gas_Constant, Mach2Vel, Mach_Motion, RefDensity, RefPressure = 0.0, factor = 0.0;
  
  su2double *Aux_Frict_x = NULL, *Aux_Frict_y = NULL, *Aux_Frict_z = NULL, *Aux_Heat = NULL, *Aux_yPlus = NULL, *Aux_Sens = NULL;
  
  unsigned short CurrentIndex;
  int *Local_Halo;
  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;
  
  int iProcessor;
  
  bool grid_movement  = (config->GetGrid_Movement());
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool transition     = (config->GetKind_Trans_Model() == LM);
  bool flow           = (( config->GetKind_Solver() == EULER             ) ||
                         ( config->GetKind_Solver() == NAVIER_STOKES     ) ||
                         ( config->GetKind_Solver() == RANS              ) ||
                         ( config->GetKind_Solver() == FEM_EULER         ) ||
                         ( config->GetKind_Solver() == FEM_NAVIER_STOKES ) ||
                         ( config->GetKind_Solver() == FEM_RANS          ) ||
                         ( config->GetKind_Solver() == FEM_LES           ) ||
                         ( config->GetKind_Solver() == ADJ_EULER         ) ||
                         ( config->GetKind_Solver() == ADJ_NAVIER_STOKES ) ||
                         ( config->GetKind_Solver() == ADJ_RANS          )   );
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);
  
  unsigned short iDim;
  unsigned short nDim = geometry->GetnDim();
  su2double RefArea = config->GetRefArea();
  su2double Gamma = config->GetGamma();
  su2double RefVel2, *Normal, Area;
  
  /*--- Set the non-dimensionalization ---*/
  if (flow) {
    if (grid_movement) {
      Gas_Constant = config->GetGas_ConstantND();
      Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
      Mach_Motion = config->GetMach_Motion();
      RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
    }
    else {
      RefVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        RefVel2  += solver[FLOW_SOL]->GetVelocity_Inf(iDim)*solver[FLOW_SOL]->GetVelocity_Inf(iDim);
    }
    RefDensity  = solver[FLOW_SOL]->GetDensity_Inf();
    RefPressure = solver[FLOW_SOL]->GetPressure_Inf();
    factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);
  }
  
  /*--- Prepare send buffers for the conservative variables. Need to
   find the total number of conservative variables and also the
   index for their particular solution container. ---*/
  
  switch (Kind_Solver) {
    case EULER : case NAVIER_STOKES: FirstIndex = FLOW_SOL; if(config->GetWeakly_Coupled_Heat()) SecondIndex = HEAT_SOL; else SecondIndex = NONE; ThirdIndex = NONE; break;
    case RANS : FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; if (transition) ThirdIndex=TRANS_SOL; else ThirdIndex = NONE;  if(config->GetWeakly_Coupled_Heat()) ThirdIndex = HEAT_SOL; else ThirdIndex = NONE; break;
    case FEM_EULER : case FEM_NAVIER_STOKES: case FEM_LES: FirstIndex = FLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case FEM_RANS : FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; break;
    case FEM_ELASTICITY: FirstIndex = FEA_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case ADJ_EULER : case ADJ_NAVIER_STOKES : FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case ADJ_RANS : FirstIndex = ADJFLOW_SOL; if (config->GetFrozen_Visc_Cont()) SecondIndex = NONE; else SecondIndex = ADJTURB_SOL; ThirdIndex = NONE; break;
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case DISC_ADJ_RANS: FirstIndex = ADJFLOW_SOL; if (config->GetFrozen_Visc_Disc()) SecondIndex = NONE; else SecondIndex = ADJTURB_SOL; ThirdIndex = NONE; break;
    case DISC_ADJ_FEM: FirstIndex = ADJFEA_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    default: SecondIndex = NONE; ThirdIndex = NONE; break;
  }
  
  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  if (ThirdIndex != NONE) nVar_Third = solver[ThirdIndex]->GetnVar();
  nVar_Consv = nVar_First + nVar_Second + nVar_Third;
  nVar_Total = nVar_Consv;
  
  if (!config->GetLow_MemoryOutput()) {
    
    /*--- Add the limiters ---*/
    
    if (config->GetWrt_Limiters()) nVar_Total += nVar_Consv;
    
    /*--- Add the residuals ---*/
    
    if (config->GetWrt_Residuals()) nVar_Total += nVar_Consv;
    
    /*--- Add the grid velocity to the restart file for the unsteady adjoint ---*/
    
    if (grid_movement && !fem) {
      iVar_GridVel = nVar_Total;
      if (geometry->GetnDim() == 2) nVar_Total += 2;
      else if (geometry->GetnDim() == 3) nVar_Total += 3;
    }
    
    /*--- Add Pressure, Temperature, Cp, Mach to the restart file ---*/
    
    if (Kind_Solver == EULER     || Kind_Solver == NAVIER_STOKES     || Kind_Solver == RANS ||
        Kind_Solver == FEM_EULER || Kind_Solver == FEM_NAVIER_STOKES || Kind_Solver == FEM_RANS ||
        Kind_Solver == FEM_LES) {
      iVar_PressCp = nVar_Total; nVar_Total += 3;
      iVar_MachMean = nVar_Total; nVar_Total += 1;
    }
    
    /*--- Add Laminar Viscosity, Skin Friction, Heat Flux, & yPlus to the restart file ---*/
    
    if (Kind_Solver == NAVIER_STOKES     || Kind_Solver == RANS     ||
        Kind_Solver == FEM_NAVIER_STOKES || Kind_Solver == FEM_RANS || Kind_Solver == FEM_LES) {
      iVar_Lam = nVar_Total;
      nVar_Total += 1;
      iVar_ViscCoeffs = nVar_Total;
      if (geometry->GetnDim() == 2) nVar_Total += 2;
      else if (geometry->GetnDim() == 3) nVar_Total += 3;
      iVar_HeatCoeffs = nVar_Total;
      nVar_Total += 2;
    }
    
    /*--- Add Eddy Viscosity to the restart file ---*/
    
    if (Kind_Solver == RANS || Kind_Solver == FEM_RANS || Kind_Solver == FEM_LES) {
      iVar_Eddy = nVar_Total; nVar_Total += 1;
    }
    
    /*--- Add Sharp edges to the restart file ---*/
    
    if (config->GetWrt_SharpEdges()) {
      if (((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
        ((Kind_Solver == FEM_EULER) || (Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {
        iVar_Sharp = nVar_Total; nVar_Total += 1;
      }
    }
    
    
    if (( Kind_Solver == ADJ_EULER              ) ||
        ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
        ( Kind_Solver == ADJ_RANS               )) {
      iVar_Sens   = nVar_Total; nVar_Total += 2;
    }
    
    if (Kind_Solver == FEM_ELASTICITY)  {
      /*--- If the analysis is dynamic... ---*/
      if (config->GetDynamic_Analysis() == DYNAMIC) {
        /*--- Velocities ---*/
        iVar_FEA_Vel = nVar_Total;
        if (geometry->GetnDim() == 2) nVar_Total += 2;
        else if (geometry->GetnDim() == 3) nVar_Total += 3;
        /*--- Accelerations ---*/
        iVar_FEA_Accel = nVar_Total;
        if (geometry->GetnDim() == 2) nVar_Total += 2;
        else if (geometry->GetnDim() == 3) nVar_Total += 3;
      }
      iVar_FEA_Stress  = nVar_Total; nVar_Total += 3;
      if (geometry->GetnDim() == 3) {iVar_FEA_Stress_3D = nVar_Total; nVar_Total += 3;}
      iVar_FEA_Extra = nVar_Total; nVar_Total += 1;
    }
    
    if ((Kind_Solver == DISC_ADJ_EULER)         ||
        (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
        (Kind_Solver == DISC_ADJ_RANS)) {
      iVar_Sens    = nVar_Total; nVar_Total += 1;
      iVar_SensDim = nVar_Total; nVar_Total += nDim;
    }
    
    if ((Kind_Solver == DISC_ADJ_FEM) && (config->GetFSI_Simulation())){
      iVar_FEA_Extra = nVar_Total; nVar_Total += 2;
    }

    if (config->GetExtraOutput()) {
      if (Kind_Solver == RANS) {
        iVar_Extra  = nVar_Total; nVar_Extra  = solver[TURB_SOL]->GetnOutputVariables(); nVar_Total += nVar_Extra;
      }
    }
    
  }
  
  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/
  
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        
        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
    
    /*--- Sum total number of nodes that belong to the domain ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;
    
  }
  Buffer_Send_nPoint[0] = nLocalPoint;
  
  /*--- Each processor sends its local number of nodes to the master. ---*/
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[size];
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalPoint = nLocalPoint;
  Buffer_Recv_nPoint[0] = Buffer_Send_nPoint[0];
#endif
  
  nBuffer_Scalar = MaxLocalPoint;
  
  /*--- Send and Recv buffers. ---*/
  
  su2double *Buffer_Send_Var = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Var = NULL;
  
  su2double *Buffer_Send_Res = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Res = NULL;
  
  su2double *Buffer_Send_Vol = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Vol = NULL;
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Auxiliary vectors for surface coefficients ---*/
  
  if (((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))  ||
      ((Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {
    Aux_Frict_x = new su2double[geometry->GetnPoint()];
    Aux_Frict_y = new su2double[geometry->GetnPoint()];
    Aux_Frict_z = new su2double[geometry->GetnPoint()];
    Aux_Heat  = new su2double[geometry->GetnPoint()];
    Aux_yPlus = new su2double[geometry->GetnPoint()];
  }
  
  if ((Kind_Solver == ADJ_EULER) ||
      (Kind_Solver == ADJ_NAVIER_STOKES) ||
      (Kind_Solver == ADJ_RANS)  ||
      (Kind_Solver == DISC_ADJ_EULER) ||
      (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
      (Kind_Solver == DISC_ADJ_RANS)) {
    Aux_Sens = new su2double[geometry->GetnPoint()];
  }
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_Var = new su2double[size*MaxLocalPoint];
    Buffer_Recv_Res = new su2double[size*MaxLocalPoint];
    Buffer_Recv_Vol = new su2double[size*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[size*MaxLocalPoint];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new su2double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new su2double[nGlobal_Poin];
    }
  }
  
  /*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/
  
  for (iVar = 0; iVar < nVar_Consv; iVar++) {
    
    /*--- Logic for which solution class to draw from. ---*/
    
    jVar = iVar;
    CurrentIndex = FirstIndex;
    if ((SecondIndex != NONE) && (iVar > nVar_First-1)) {
      jVar = iVar - nVar_First;
      CurrentIndex = SecondIndex;
    }
    if ((SecondIndex != NONE) && (ThirdIndex != NONE) && (iVar > (nVar_First + nVar_Second-1))) {
      jVar = iVar - nVar_First - nVar_Second;
      CurrentIndex = ThirdIndex;
    }
    
    /*--- Loop over this partition to collect the current variable ---*/
    
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos & write only if requested ---*/
      
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Get this variable into the temporary send buffer. ---*/
        
        Buffer_Send_Var[jPoint] = solver[CurrentIndex]->node[iPoint]->GetSolution(jVar);
        
        if (!config->GetLow_MemoryOutput()) {
          
          if (config->GetWrt_Limiters()) {
            Buffer_Send_Vol[jPoint] = solver[CurrentIndex]->node[iPoint]->GetLimiter_Primitive(jVar);
          }
          
          if (config->GetWrt_Residuals()) {
            if (!config->GetDiscrete_Adjoint()) {
              Buffer_Send_Res[jPoint] = solver[CurrentIndex]->LinSysRes.GetBlock(iPoint, jVar);
            } else {
              Buffer_Send_Res[jPoint] = solver[CurrentIndex]->node[iPoint]->GetSolution(jVar) -
              solver[CurrentIndex]->node[iPoint]->GetSolution_Old(jVar);
            }
          }
          
        }
        
        /*--- Only send/recv the volumes & global indices during the first loop ---*/
        
        if (iVar == 0) {
          Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
        }
        
        jPoint++;
        
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    
#ifdef HAVE_MPI
    SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
    if (!config->GetLow_MemoryOutput()) {
      
      if (config->GetWrt_Limiters()) {
#ifdef HAVE_MPI
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif
      }
      
      if (config->GetWrt_Residuals()) {
#ifdef HAVE_MPI
        SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
#endif
      }
      
    }
    
    if (iVar == 0) {
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_GlobalIndex[iPoint] = Buffer_Send_GlobalIndex[iPoint];
#endif
    }
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    
    if (rank == MASTER_NODE) {
      jPoint = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          
          if (!config->GetLow_MemoryOutput()) {
            
            if (config->GetWrt_Limiters()) {
              Data[iVar+nVar_Consv][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            }
            
            if (config->GetWrt_Residuals()) {
              unsigned short ExtraIndex;
              ExtraIndex = nVar_Consv;
              if (config->GetWrt_Limiters()) ExtraIndex = 2*nVar_Consv;
              Data[iVar+ExtraIndex][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            }
            
          }
          
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
    
  }
  
  if (!config->GetLow_MemoryOutput()) {
    
    /*--- Additional communication routine for the grid velocity. Note that
     we are reusing the same temporary buffers from above for efficiency.
     Also, in the future more routines like this could be used to write
     an arbitrary number of additional variables to the file. ---*/
    
    if (grid_movement && !fem) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0; su2double *Grid_Vel;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the three grid velocity components. ---*/
          
          Grid_Vel = geometry->node[iPoint]->GetGridVel();
          Buffer_Send_Var[jPoint] = Grid_Vel[0];
          Buffer_Send_Res[jPoint] = Grid_Vel[1];
          if (geometry->GetnDim() == 3) Buffer_Send_Vol[jPoint] = Grid_Vel[2];
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (geometry->GetnDim() == 3) {
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      }
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (geometry->GetnDim() == 3) {
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
      }
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_GridVel;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            if (geometry->GetnDim() == 3)
              Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }
    
    /*--- Communicate Pressure, Cp, and Mach ---*/
    
    if (((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
      ((Kind_Solver == FEM_EULER) || (Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {
      
      /*--- First, loop through the mesh in order to find and store the
       value of the coefficient of pressure at any surface nodes. They
       will be placed in an auxiliary vector and then communicated like
       all other volumetric variables. ---*/
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the pressure, Cp, and mach variables. ---*/

          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure();
          if (compressible){
            Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature();
          } else{
            Buffer_Send_Res[jPoint] =  0.0;
          }
          Buffer_Send_Vol[jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure() - RefPressure)*factor*RefArea;

          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_PressCp;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }
    
    /*--- Communicate Mach---*/
    
    if (((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
      ((Kind_Solver == FEM_EULER) || (Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
          
          if (compressible) {
            Buffer_Send_Var[jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
            solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed();
          }
          if (incompressible) {
            Buffer_Send_Var[jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())*config->GetVelocity_Ref()/
            sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->node[iPoint]->GetDensity()*config->GetDensity_Ref()));
          }
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_MachMean;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }
    
    /*--- Laminar Viscosity ---*/
    
    if (((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
      ((Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
          
          Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();

          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Lam;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
      
      /*--- Communicate skin friction ---*/
      
      /*--- First, loop through the mesh in order to find and store the
       value of the viscous coefficients at any surface nodes. They
       will be placed in an auxiliary vector and then communicated like
       all other volumetric variables. ---*/
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        Aux_Frict_x[iPoint] = 0.0;
        Aux_Frict_y[iPoint] = 0.0;
        Aux_Frict_z[iPoint] = 0.0;
      }
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
        if (config->GetMarker_All_Plotting(iMarker) == YES) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            Aux_Frict_x[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 0);
            Aux_Frict_y[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 1);
            if (geometry->GetnDim() == 3) Aux_Frict_z[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 2);
          }
        }
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the three grid velocity components. ---*/
          
          Buffer_Send_Var[jPoint] = Aux_Frict_x[iPoint];
          Buffer_Send_Res[jPoint] = Aux_Frict_y[iPoint];
          if (geometry->GetnDim() == 3)
            Buffer_Send_Vol[jPoint] = Aux_Frict_z[iPoint];
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (geometry->GetnDim() == 3) {
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      }
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
        Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (geometry->GetnDim() == 3) {
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
          Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
      }
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0;
        iVar = iVar_ViscCoeffs;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            Data[iVar + 1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            if (geometry->GetnDim() == 3)
              Data[iVar + 2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor + 1) * nBuffer_Scalar;
        }
      }
      
      /*--- Communicate heat transfer, y+ ---*/
      
      /*--- First, loop through the mesh in order to find and store the
       value of the viscous coefficients at any surface nodes. They
       will be placed in an auxiliary vector and then communicated like
       all other volumetric variables. ---*/
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        Aux_Heat[iPoint] = 0.0;
        Aux_yPlus[iPoint] = 0.0;
      }
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
        if (config->GetMarker_All_Plotting(iMarker) == YES) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            Aux_Heat[iPoint] = solver[FLOW_SOL]->GetHeatFlux(iMarker, iVertex);
            Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker, iVertex);
          }
        }
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/
          
          if (compressible) {
            Buffer_Send_Res[jPoint] = Aux_Heat[iPoint];
            Buffer_Send_Vol[jPoint] = Aux_yPlus[iPoint];
          }
          if (incompressible) {
            Buffer_Send_Res[jPoint] = Aux_Heat[iPoint];
            Buffer_Send_Vol[jPoint] = Aux_yPlus[iPoint];
          }
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
        Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
        Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0;
        iVar = iVar_HeatCoeffs;

        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar + 0][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            Data[iVar + 1][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor + 1) * nBuffer_Scalar;
        }
      }
    }
    
    
    /*--- Communicate the Eddy Viscosity ---*/
    
    if (Kind_Solver == RANS || Kind_Solver == FEM_RANS || Kind_Solver == FEM_LES) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the pressure and mach variables. ---*/
          
          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosity();

          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Eddy;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
      
    }
    
    /*--- Communicate the Sharp Edges ---*/
    
    if (config->GetWrt_SharpEdges()) {
      
      if (((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
        ((Kind_Solver == FEM_EULER) || (Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {
        
        /*--- Loop over this partition to collect the current variable ---*/
        jPoint = 0;
        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
          
          /*--- Check for halos & write only if requested ---*/
          
          if (!Local_Halo[iPoint] || Wrt_Halo) {
            
            /*--- Load buffers with the pressure and mach variables. ---*/
            
            Buffer_Send_Var[jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance();
            jPoint++;
          }
        }
        
        /*--- Gather the data on the master node. ---*/
        
#ifdef HAVE_MPI
        SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
        
        /*--- The master node unpacks and sorts this variable by global index ---*/
        
        if (rank == MASTER_NODE) {
          jPoint = 0; iVar = iVar_Sharp;
          for (iProcessor = 0; iProcessor < size; iProcessor++) {
            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
              
              /*--- Get global index, then loop over each variable and store ---*/
              
              iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
              Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
              jPoint++;
            }
            
            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
            
            jPoint = (iProcessor+1)*nBuffer_Scalar;
          }
        }
      }
    }
    
    /*--- Communicate the surface sensitivity ---*/
    
    if ((Kind_Solver == ADJ_EULER)         ||
        (Kind_Solver == ADJ_NAVIER_STOKES) ||
        (Kind_Solver == ADJ_RANS)          ||
        (Kind_Solver == DISC_ADJ_EULER)    ||
        (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
        (Kind_Solver == DISC_ADJ_RANS)) {
      
      /*--- First, loop through the mesh in order to find and store the
       value of the surface sensitivity at any surface nodes. They
       will be placed in an auxiliary vector and then communicated like
       all other volumetric variables. ---*/
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Aux_Sens[iPoint] = 0.0;
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
        if (config->GetMarker_All_Plotting(iMarker) == YES) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
            Area = sqrt (Area);
            Aux_Sens[iPoint] = solver[ADJFLOW_SOL]->GetCSensitivity(iMarker, iVertex)/Area;
          }
        }
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/
          
          Buffer_Send_Var[jPoint] = Aux_Sens[iPoint];
          if ((config->GetKind_ConvNumScheme() == SPACE_CENTERED) && (!config->GetDiscrete_Adjoint()))
            Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensor(iPoint);
          if ((config->GetKind_ConvNumScheme() == SPACE_UPWIND) && (!config->GetDiscrete_Adjoint()))
            Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetLimiter(0);
          
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (!config->GetDiscrete_Adjoint())
        SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      if (!config->GetDiscrete_Adjoint())
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Sens;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar+0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            if (!config->GetDiscrete_Adjoint())
              Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }
    
    if ((Kind_Solver == DISC_ADJ_EULER)    ||
        (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
        (Kind_Solver == DISC_ADJ_RANS)) {
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/
          
          Buffer_Send_Var[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensitivity(0);
          Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensitivity(1);
          if (nDim == 3)
            Buffer_Send_Vol[jPoint] = solver[ADJFLOW_SOL]->node[iPoint]->GetSensitivity(2);
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (nDim == 3)
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (nDim == 3)
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_SensDim;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar+0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            if (nDim == 3)
              Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }
    
    
    /*--- Communicate the Velocities for dynamic FEM problem ---*/
    
    if ((Kind_Solver == FEM_ELASTICITY) && (config->GetDynamic_Analysis() == DYNAMIC)) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0; su2double *Node_Vel;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the three grid velocity components. ---*/
          
          Node_Vel = solver[FEA_SOL]->node[iPoint]->GetSolution_Vel();
          Buffer_Send_Var[jPoint] = Node_Vel[0];
          Buffer_Send_Res[jPoint] = Node_Vel[1];
          if (geometry->GetnDim() == 3) Buffer_Send_Vol[jPoint] = Node_Vel[2];
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (geometry->GetnDim() == 3) {
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      }
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (geometry->GetnDim() == 3) {
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
      }
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Vel;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            if (geometry->GetnDim() == 3)
              Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }
    
    /*--- Communicate the Accelerations for dynamic FEM problem ---*/
    
    if ((Kind_Solver == FEM_ELASTICITY) && (config->GetDynamic_Analysis() == DYNAMIC)) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0; su2double *Node_Accel;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the three grid velocity components. ---*/
          
          Node_Accel = solver[FEA_SOL]->node[iPoint]->GetSolution_Accel();
          Buffer_Send_Var[jPoint] = Node_Accel[0];
          Buffer_Send_Res[jPoint] = Node_Accel[1];
          if (geometry->GetnDim() == 3) Buffer_Send_Vol[jPoint] = Node_Accel[2];
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (geometry->GetnDim() == 3) {
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      }
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (geometry->GetnDim() == 3) {
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
      }
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Accel;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            if (geometry->GetnDim() == 3)
              Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }

    /*--- Communicate the FEM elasticity stresses (2D) - New elasticity solver---*/
    
    if (Kind_Solver == FEM_ELASTICITY) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0; su2double *Stress;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the three grid velocity components. ---*/
          
          Stress = solver[FEA_SOL]->node[iPoint]->GetStress_FEM();
          /*--- Sigma xx ---*/
          Buffer_Send_Var[jPoint] = Stress[0];
          /*--- Sigma yy ---*/
          Buffer_Send_Res[jPoint] = Stress[1];
          /*--- Sigma xy ---*/
          Buffer_Send_Vol[jPoint] = Stress[2];
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Stress;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }
    
    /*--- Communicate the FEM elasticity stresses (3D) - New elasticity solver---*/
    
    if ((Kind_Solver == FEM_ELASTICITY) && (geometry->GetnDim() == 3)) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      
      jPoint = 0; su2double *Stress;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the three grid velocity components. ---*/
          
          Stress = solver[FEA_SOL]->node[iPoint]->GetStress_FEM();
          /*--- Sigma zz ---*/
          Buffer_Send_Var[jPoint] = Stress[3];
          /*--- Sigma xz ---*/
          Buffer_Send_Res[jPoint] = Stress[4];
          /*--- Sigma yz ---*/
          Buffer_Send_Vol[jPoint] = Stress[5];
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
      
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Stress_3D;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }
    
    
    /*--- Communicate the Linear elasticity ---*/
    
    if ( Kind_Solver == FEM_ELASTICITY ) {
      
      /*--- Loop over this partition to collect the current variable ---*/
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Check for halos & write only if requested ---*/
        
        if (!Local_Halo[iPoint] || Wrt_Halo) {
          
          /*--- Load buffers with the temperature and laminar viscosity variables. ---*/
          
          Buffer_Send_Var[jPoint] = solver[FEA_SOL]->node[iPoint]->GetVonMises_Stress();
          jPoint++;
        }
      }
      
      /*--- Gather the data on the master node. ---*/
      
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
      
      /*--- The master node unpacks and sorts this variable by global index ---*/
      
      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Extra;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
            
            /*--- Get global index, then loop over each variable and store ---*/
            
            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            jPoint++;
          }
          
          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
          
          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }
    
    if ((Kind_Solver == DISC_ADJ_FEM) && (config->GetFSI_Simulation())) {
      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/

          Buffer_Send_Var[jPoint] = solver[ADJFEA_SOL]->node[iPoint]->GetGeometry_CrossTerm_Derivative(0);
          Buffer_Send_Res[jPoint] = solver[ADJFEA_SOL]->node[iPoint]->GetGeometry_CrossTerm_Derivative(1);
          if (geometry->GetnDim() == 3)
            Buffer_Send_Vol[jPoint] = solver[ADJFEA_SOL]->node[iPoint]->GetGeometry_CrossTerm_Derivative(2);
          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (nDim == 3)
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (nDim == 3)
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Extra;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar+0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            if (nDim == 3)
              Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }

    if (config->GetExtraOutput()) {
      
      for (jVar = 0; jVar < nVar_Extra; jVar++) {
        
        /*--- Loop over this partition to collect the current variable ---*/
        
        jPoint = 0;
        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
          
          /*--- Check for halos & write only if requested ---*/
          
          if (!Local_Halo[iPoint] || Wrt_Halo) {
            
            /*--- Get this variable into the temporary send buffer. ---*/
            
            if (Kind_Solver == RANS) {
              Buffer_Send_Var[jPoint] = solver[TURB_SOL]->OutputVariables[iPoint*nVar_Extra+jVar];
            }
            jPoint++;
            
          }
        }
        
        /*--- Gather the data on the master node. ---*/
        
#ifdef HAVE_MPI
        SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif
        
        /*--- The master node unpacks and sorts this variable by global index ---*/
        
        if (rank == MASTER_NODE) {
          jPoint = 0; iVar = iVar_Extra;
          for (iProcessor = 0; iProcessor < size; iProcessor++) {
            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
              
              /*--- Get global index, then loop over each variable and store ---*/
              
              iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
              Data[iVar+jVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
              jPoint++;
            }
            
            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
            
            jPoint = (iProcessor+1)*nBuffer_Scalar;
          }
        }
      }
    }
    
  }
  
  /*--- Immediately release the temporary buffers. ---*/
  
  delete [] Buffer_Send_Var;
  delete [] Buffer_Send_Res;
  delete [] Buffer_Send_Vol;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nPoint;
    delete [] Buffer_Recv_Var;
    delete [] Buffer_Recv_Res;
    delete [] Buffer_Recv_Vol;
    delete [] Buffer_Recv_GlobalIndex;
  }
  
  /*--- Release memory needed for surface coefficients ---*/
  
  delete [] Local_Halo;
  
  if (((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
    ((Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {
    delete[] Aux_Frict_x; delete[] Aux_Frict_y; delete[] Aux_Frict_z;
    delete [] Aux_Heat; delete [] Aux_yPlus;
  }
  if (( Kind_Solver == ADJ_EULER              ) ||
      ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
      ( Kind_Solver == ADJ_RANS               ) ||
      ( Kind_Solver == DISC_ADJ_EULER         ) ||
      ( Kind_Solver == DISC_ADJ_NAVIER_STOKES ) ||
      ( Kind_Solver == DISC_ADJ_RANS          )) {
    delete [] Aux_Sens;
  }
  
}

void COutput::MergeBaselineSolution(CConfig *config, CGeometry *geometry, CSolver *solver, unsigned short val_iZone) {
  
  /*--- Local variables needed on all processors ---*/
  unsigned short iVar;
  unsigned long iPoint = 0, jPoint = 0;
  
  nVar_Total = config->fields.size() - 1;
  
  /*--- Merge the solution either in serial or parallel. ---*/
  
#ifndef HAVE_MPI
  
  /*--- In serial, the single process has access to all solution data,
   so it is simple to retrieve and store inside Solution_Data. ---*/
  
  unsigned short iMarker;
  unsigned long iVertex, nTotalPoints = 0;
  int SendRecv;
  
  /*--- First, create a structure to locate any periodic halo nodes ---*/
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
            (SendRecv < 0)) {
          Local_Halo[iPoint] = false;
        }
      }
      
    }
  }
  
  /*--- Total number of points in the mesh (this might include periodic points). ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (!Local_Halo[iPoint]) nTotalPoints++;
  
  nGlobal_Poin = nTotalPoints;
  Data = new su2double*[nVar_Total];
  for (iVar = 0; iVar < nVar_Total; iVar++) {
    Data[iVar] = new su2double[nGlobal_Poin];
  }
  
  /*--- Loop over all points in the mesh, but only write data
   for nodes in the domain (ignore periodic halo nodes). ---*/
  
  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    if (!Local_Halo[iPoint]) {
      
      /*--- Solution (first, and second system of equations) ---*/
      
      unsigned short jVar = 0;
      for (iVar = 0; iVar < nVar_Total; iVar++) {
        Data[jVar][jPoint] = solver->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
    }
    
    /*--- Increment jPoint as the counter. We need this because iPoint
     may include halo nodes that we skip over during this loop. ---*/
    
    jPoint++;
    
  }
  
#else
  
  /*--- MPI preprocessing ---*/
  
  int nProcessor = size, iProcessor;
  
  /*--- Local variables needed for merging with MPI ---*/
  
  unsigned long iVertex, iMarker;
  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/
  
  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        
        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
    
    /*--- Sum total number of nodes that belong to the domain ---*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;
    
  }
  Buffer_Send_nPoint[0] = nLocalPoint;
  
  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];
  
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  
  nBuffer_Scalar = MaxLocalPoint;
  
  /*--- Send and Recv buffers. ---*/
  
  su2double *Buffer_Send_Var = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Var = NULL;
  
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;
  
  /*--- Prepare the receive buffers in the master node only. ---*/
  if (rank == MASTER_NODE) {
    
    Buffer_Recv_Var = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];
    
    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new su2double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new su2double[nGlobal_Poin];
    }
    
  }
  
  /*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/
  
  for (iVar = 0; iVar < nVar_Total; iVar++) {
    
    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Check for halos and write only if requested ---*/
      if (!Local_Halo[iPoint] || Wrt_Halo) {
        
        /*--- Get this variable into the temporary send buffer. ---*/
        Buffer_Send_Var[jPoint] = solver->node[iPoint]->GetSolution(iVar);
        
        /*--- Only send/recv the volumes & global indices during the first loop ---*/
        if (iVar == 0) {
          Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
        }
        jPoint++;
      }
    }
    
    /*--- Gather the data on the master node. ---*/
    
    SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    if (iVar == 0) {
      SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    }
    
    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {
          
          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }
  
  /*--- Immediately release the temporary buffers. ---*/
  
  delete [] Buffer_Send_Var;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_Var;
    delete [] Buffer_Recv_GlobalIndex;
  }
  
#endif
  
  delete [] Local_Halo;
  
}

void COutput::SetRestart(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables ---*/
  
  unsigned short nZone = geometry->GetnZone();
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iExtIter = config->GetExtIter();
  bool grid_movement = config->GetGrid_Movement();
  bool dynamic_fem = (config->GetDynamic_Analysis() == DYNAMIC);
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);
  bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);
  ofstream restart_file;
  ofstream meta_file;
  string filename, meta_filename;
  bool adjoint = config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Retrieve filename from config ---*/
  
  if (((config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint())) && ((config->GetKind_Solver() != DISC_ADJ_FEM)))  {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else if (fem) {
    filename = config->GetRestart_FEMFileName();
  } else if (disc_adj_fem){
    filename = config->GetRestart_AdjFEMFileName();
  } else {
    filename = config->GetRestart_FlowFileName();
  }
  
  /*--- Append the zone number if multizone problems ---*/
  if (nZone > 1)
    filename= config->GetMultizone_FileName(filename, val_iZone);
  
  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(config->GetiInst()));
  } else if (config->GetWrt_Unsteady()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  } else if ((fem || disc_adj_fem) && (config->GetWrt_Dynamic())) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  }
  
  /*--- Open the restart file and write the solution. ---*/
  
  restart_file.open(filename.c_str(), ios::out);
  restart_file.precision(15);
  
  /*--- Write the header line based on the particular solver ----*/
  
  restart_file << "\"PointID\"";
  
  /*--- Mesh coordinates are always written to the restart first ---*/
  
  if (nDim == 2) {
    restart_file << "\t\"x\"\t\"y\"";
  } else {
    restart_file << "\t\"x\"\t\"y\"\t\"z\"";
  }
  
  for (iVar = 0; iVar < nVar_Consv; iVar++) {
  if (( Kind_Solver == FEM_ELASTICITY ) || ( Kind_Solver == DISC_ADJ_FEM))
      restart_file << "\t\"Displacement_" << iVar+1<<"\"";
    else
      restart_file << "\t\"Conservative_" << iVar+1<<"\"";
  }
  
  if (!config->GetLow_MemoryOutput()) {
    
    if (config->GetWrt_Limiters()) {
      for (iVar = 0; iVar < nVar_Consv; iVar++) {
        restart_file << "\t\"Limiter_" << iVar+1<<"\"";
      }
    }
    if (config->GetWrt_Residuals()) {
      for (iVar = 0; iVar < nVar_Consv; iVar++) {
        restart_file << "\t\"Residual_" << iVar+1<<"\"";
      }
    }
    
    /*--- Mesh velocities for dynamic mesh cases ---*/
    
    if (grid_movement && !fem) {
      if (nDim == 2) {
        restart_file << "\t\"Grid_Velx\"\t\"Grid_Vely\"";
      } else {
        restart_file << "\t\"Grid_Velx\"\t\"Grid_Vely\"\t\"Grid_Velz\"";
      }
    }
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      if ((config->GetOutput_FileFormat() == PARAVIEW) || (config->GetOutput_FileFormat() == PARAVIEW_BINARY)) {
        restart_file << "\t\"Pressure\"\t\"Temperature\"\t\"Pressure_Coefficient\"\t\"Mach\"";
      } else
        restart_file << "\t\"Pressure\"\t\"Temperature\"\t\"C<sub>p</sub>\"\t\"Mach\"";
    }
    
    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      if ((config->GetOutput_FileFormat() == PARAVIEW) || (config->GetOutput_FileFormat() == PARAVIEW_BINARY)) {
        if (nDim == 2) restart_file << "\t\"Laminar_Viscosity\"\t\"Skin_Friction_Coefficient_X\"\t\"Skin_Friction_Coefficient_Y\"\t\"Heat_Flux\"\t\"Y_Plus\"";
        if (nDim == 3) restart_file << "\t\"Laminar_Viscosity\"\t\"Skin_Friction_Coefficient_X\"\t\"Skin_Friction_Coefficient_Y\"\t\"Skin_Friction_Coefficient_Z\"\t\"Heat_Flux\"\t\"Y_Plus\"";
      } else {
        if (nDim == 2) restart_file << "\t\"<greek>m</greek>\"\t\"C<sub>f</sub>_x\"\t\"C<sub>f</sub>_y\"\t\"h\"\t\"y<sup>+</sup>\"";
        if (nDim == 3) restart_file << "\t\"<greek>m</greek>\"\t\"C<sub>f</sub>_x\"\t\"C<sub>f</sub>_y\"\t\"C<sub>f</sub>_z\"\t\"h\"\t\"y<sup>+</sup>\"";
      }
    }
    
    if (Kind_Solver == RANS) {
      if ((config->GetOutput_FileFormat() == PARAVIEW) || (config->GetOutput_FileFormat() == PARAVIEW_BINARY)) {
        restart_file << "\t\"Eddy_Viscosity\"";
      } else
        restart_file << "\t\"<greek>m</greek><sub>t</sub>\"";
    }
    
    if (config->GetWrt_SharpEdges()) {
      if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
        restart_file << "\t\"Sharp_Edge_Dist\"";
      }
    }
    
    if ((Kind_Solver == ADJ_EULER              ) ||
        (Kind_Solver == ADJ_NAVIER_STOKES      ) ||
        (Kind_Solver == ADJ_RANS               )   ) {
      restart_file << "\t\"Surface_Sensitivity\"\t\"Solution_Sensor\"";
    }
    if (( Kind_Solver == DISC_ADJ_EULER              ) ||
        ( Kind_Solver == DISC_ADJ_NAVIER_STOKES      ) ||
        ( Kind_Solver == DISC_ADJ_RANS               )) {
      restart_file << "\t\"Surface_Sensitivity\"\t\"Sensitivity_x\"\t\"Sensitivity_y\"";
      if (geometry->GetnDim() == 3) {
        restart_file << "\t\"Sensitivity_z\"";
      }
    }
    
    if (Kind_Solver == FEM_ELASTICITY) {
      if (!dynamic_fem) {
        if (geometry->GetnDim() == 2)
          restart_file << "\t\"Sxx\"\t\"Syy\"\t\"Sxy\"\t\"Von_Mises_Stress\"";
        if (geometry->GetnDim() == 3)
          restart_file << "\t\"Sxx\"\t\"Syy\"\t\"Sxy\"\t\"Szz\"\t\"Sxz\"\t\"Syz\"\t\"Von_Mises_Stress\"";
      }
      else if (dynamic_fem) {
        if (geometry->GetnDim() == 2) {
          restart_file << "\t\"Velocity_1\"\t\"Velocity_2\"\t\"Acceleration_1\"\t\"Acceleration_2\"";
          restart_file << "\t\"Sxx\"\t\"Syy\"\t\"Sxy\"\t\"Von_Mises_Stress\"";
        }
        if (geometry->GetnDim() == 3) {
          restart_file << "\t\"Velocity_1\"\t\"Velocity_2\"\t\"Velocity_3\"\t\"Acceleration_1\"\t\"Acceleration_2\"\t\"Acceleration_3\"";
          restart_file << "\t\"Sxx\"\t\"Syy\"\t\"Sxy\"\t\"Szz\"\t\"Sxz\"\t\"Syz\"\t\"Von_Mises_Stress\"";
        }
      }
    }
    
    if ((Kind_Solver == DISC_ADJ_FEM) && (config->GetFSI_Simulation())){
      if (geometry->GetnDim() == 2)
        restart_file << "\t\"CrossTerm_1\"\t\"CrossTerm_2\"";
      if (geometry->GetnDim() == 3)
        restart_file << "\t\"CrossTerm_1\"\t\"CrossTerm_2\"\t\"CrossTerm_3\"";
    }

    
    if (config->GetExtraOutput()) {
      string *headings = NULL;
      //if (Kind_Solver == RANS) {
      headings = solver[TURB_SOL]->OutputHeadingNames;
      //}
      
      for (iVar = 0; iVar < nVar_Extra; iVar++) {
        if (headings == NULL) {
          restart_file << "\t\"ExtraOutput_" << iVar+1<<"\"";
        } else {
          restart_file << "\t\""<< headings[iVar] <<"\"";
        }
      }
    }
  }
  
  restart_file << "\n";
  
  /*--- Write the restart file ---*/
  
  for (iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
    
    /*--- Index of the point ---*/
    restart_file << iPoint << "\t";
    
    /*--- Write the grid coordinates first ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      restart_file << scientific << Coords[iDim][iPoint] << "\t";
    }
    
    /*--- Loop over the variables and write the values to file ---*/
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      restart_file << scientific << Data[iVar][iPoint] << "\t";
    }
    restart_file << "\n";
  }

  /*--- Write the general header and flow conditions ----*/

  if (dual_time)
    restart_file <<"EXT_ITER= " << config->GetExtIter() + 1 << endl;
  else
    restart_file <<"EXT_ITER= " << config->GetExtIter() + config->GetExtIter_OffSet() + 1 << endl;
  restart_file <<"AOA= " << config->GetAoA() - config->GetAoA_Offset() << endl;
  restart_file <<"SIDESLIP_ANGLE= " << config->GetAoS() - config->GetAoS_Offset() << endl;
  restart_file <<"INITIAL_BCTHRUST= " << config->GetInitial_BCThrust() << endl;
  restart_file <<"DCD_DCL_VALUE= " << config->GetdCD_dCL() << endl;
  restart_file <<"DCMX_DCL_VALUE= " << config->GetdCMx_dCL() << endl;
  restart_file <<"DCMY_DCL_VALUE= " << config->GetdCMy_dCL() << endl;
  restart_file <<"DCMZ_DCL_VALUE= " << config->GetdCMz_dCL() << endl;
  if (adjoint) restart_file << "SENS_AOA=" << solver[ADJFLOW_SOL]->GetTotal_Sens_AoA() * PI_NUMBER / 180.0 << endl;

  /*--- Close the data portion of the restart file. ---*/

  restart_file.close();

}

void COutput::DeallocateCoordinates(CConfig *config, CGeometry *geometry) {
  
  unsigned short iDim, nDim = geometry->GetnDim();
  
  /*--- The master node alone owns all data found in this routine. ---*/
  
  if (rank == MASTER_NODE) {
    
    /*--- Deallocate memory for coordinate data ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      delete [] Coords[iDim];
    }
    delete [] Coords;
  }
}

void COutput::DeallocateConnectivity(CConfig *config, CGeometry *geometry, bool surf_sol) {

  /*--- The master node alone owns all data found in this routine. ---*/
  if (rank == MASTER_NODE) {
    
    /*--- Deallocate memory for connectivity data ---*/
    if (surf_sol) {
      if (nGlobal_Line > 0      && Conn_Line      != NULL) delete [] Conn_Line;
      if (nGlobal_BoundTria > 0 && Conn_BoundTria != NULL) delete [] Conn_BoundTria;
      if (nGlobal_BoundQuad > 0 && Conn_BoundQuad != NULL) delete [] Conn_BoundQuad;
    }
    else {
      if (nGlobal_Tria > 0 && Conn_Tria != NULL) delete [] Conn_Tria;
      if (nGlobal_Quad > 0 && Conn_Quad != NULL) delete [] Conn_Quad;
      if (nGlobal_Tetr > 0 && Conn_Tetr != NULL) delete [] Conn_Tetr;
      if (nGlobal_Hexa > 0 && Conn_Hexa != NULL) delete [] Conn_Hexa;
      if (nGlobal_Pris > 0 && Conn_Pris != NULL) delete [] Conn_Pris;
      if (nGlobal_Pyra > 0 && Conn_Pyra != NULL) delete [] Conn_Pyra;
      
    }
    
  }
}

void COutput::DeallocateSolution(CConfig *config, CGeometry *geometry) {

  /*--- The master node alone owns all data found in this routine. ---*/
  if (rank == MASTER_NODE) {
    
    /*--- Deallocate memory for solution data ---*/
    for (unsigned short iVar = 0; iVar < nVar_Total; iVar++) {
      delete [] Data[iVar];
    }
    delete [] Data;
    
  }
}

void COutput::SetConvHistory_Body(CGeometry ****geometry,
                                     CSolver *****solver_container,
                                     CConfig **config,
                                     CIntegration ****integration,
                                     bool DualTime_Iteration,
                                     su2double timeused,
                                     unsigned short val_iZone,
                                     unsigned short val_iInst) {


  /*--- Output using only the master node ---*/

  if (rank == MASTER_NODE) {

    bool write_header, write_history, write_screen;

    /*--- Retrieve residual and extra data -----------------------------------------------------------------*/

    LoadHistoryData(geometry, solver_container, config, integration, DualTime_Iteration,
                    timeused, val_iZone, val_iInst);
    
    Postprocess_HistoryData(config[val_iZone], DualTime_Iteration);
        
    /*--- Write the history file ---------------------------------------------------------------------------*/
    write_history = WriteHistoryFile_Output(config[val_iZone], DualTime_Iteration);
    if (write_history) SetHistoryFile_Output(config[val_iZone]);

    /*--- Write the screen header---------------------------------------------------------------------------*/
    write_header = WriteScreen_Header(config[val_iZone]);
    if (write_header) SetScreen_Header(config[val_iZone]);

    /*--- Write the screen output---------------------------------------------------------------------------*/
    write_screen = WriteScreen_Output(config[val_iZone], DualTime_Iteration);
    if (write_screen) SetScreen_Output(config[val_iZone]);

  }

}

void COutput::SetCFL_Number(CSolver *****solver_container, CConfig **config, unsigned short val_iZone) {
  
  su2double CFLFactor = 1.0, power = 1.0, CFL = 0.0, CFLMin = 0.0, CFLMax = 0.0, Div = 1.0, Diff = 0.0, MGFactor[100];
  unsigned short iMesh;
  
  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  unsigned long ExtIter = config[val_iZone]->GetExtIter();
  unsigned short nVar = 1;

  bool energy = config[val_iZone]->GetEnergy_Equation();
  bool weakly_coupled_heat = config[val_iZone]->GetWeakly_Coupled_Heat();

  switch( config[val_iZone]->GetKind_Solver()) {
    case EULER : case NAVIER_STOKES : case RANS:
      if (energy) {
        nVar = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetnVar();
        RhoRes_New = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetRes_RMS(nVar-1);
      }
      else if (weakly_coupled_heat) {
        RhoRes_New = solver_container[val_iZone][INST_0][FinestMesh][HEAT_SOL]->GetRes_RMS(0);
      }
      else {
        RhoRes_New = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetRes_RMS(0);
      }
      break;
    case ADJ_EULER : case ADJ_NAVIER_STOKES: case ADJ_RANS:
      RhoRes_New = solver_container[val_iZone][INST_0][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(0);
      break;
    case HEAT_EQUATION_FVM:
      RhoRes_New = solver_container[val_iZone][INST_0][FinestMesh][HEAT_SOL]->GetRes_RMS(0);
      break;
  }
  
  if (RhoRes_New < EPS) RhoRes_New = EPS;
  if (RhoRes_Old[val_iZone] < EPS) RhoRes_Old[val_iZone] = RhoRes_New;

  Div = RhoRes_Old[val_iZone]/RhoRes_New;
  Diff = RhoRes_New-RhoRes_Old[val_iZone];

  /*--- Compute MG factor ---*/

  for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
    if (iMesh == MESH_0) MGFactor[iMesh] = 1.0;
    else MGFactor[iMesh] = MGFactor[iMesh-1] * config[val_iZone]->GetCFL(iMesh)/config[val_iZone]->GetCFL(iMesh-1);
  }

  if (Div < 1.0) power = config[val_iZone]->GetCFL_AdaptParam(0);
  else power = config[val_iZone]->GetCFL_AdaptParam(1);

  /*--- Detect a stall in the residual ---*/

  if ((fabs(Diff) <= RhoRes_New*1E-8) && (ExtIter != 0)) { Div = 0.1; power = config[val_iZone]->GetCFL_AdaptParam(1); }

  CFLMin = config[val_iZone]->GetCFL_AdaptParam(2);
  CFLMax = config[val_iZone]->GetCFL_AdaptParam(3);

  CFLFactor = pow(Div, power);

  for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
    CFL = config[val_iZone]->GetCFL(iMesh);
    CFL *= CFLFactor;

    if ((iMesh == MESH_0) && (CFL <= CFLMin)) {
      for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
        config[val_iZone]->SetCFL(iMesh, 1.001*CFLMin*MGFactor[iMesh]);
      }
      break;
    }
    if ((iMesh == MESH_0) && (CFL >= CFLMax)) {
      for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++)
        config[val_iZone]->SetCFL(iMesh, 0.999*CFLMax*MGFactor[iMesh]);
      break;
    }

    config[val_iZone]->SetCFL(iMesh, CFL);
  }

  switch( config[val_iZone]->GetKind_Solver()) {
  case EULER : case NAVIER_STOKES : case RANS:
    nVar = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetnVar();
    if (energy) RhoRes_Old[val_iZone] = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetRes_RMS(nVar-1);
    else if (weakly_coupled_heat) RhoRes_Old[val_iZone] = solver_container[val_iZone][INST_0][FinestMesh][HEAT_SOL]->GetRes_RMS(0);
    else RhoRes_Old[val_iZone] = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetRes_RMS(0);
    break;
  case ADJ_EULER : case ADJ_NAVIER_STOKES: case ADJ_RANS:
    RhoRes_Old[val_iZone] = solver_container[val_iZone][INST_0][FinestMesh][ADJFLOW_SOL]->GetRes_RMS(0);
    break;
  case HEAT_EQUATION_FVM:
    RhoRes_Old[val_iZone] = solver_container[val_iZone][INST_0][FinestMesh][HEAT_SOL]->GetRes_RMS(0);
    break;
  }
  
}

void COutput::SetBaselineResult_Files(CSolver ***solver, CGeometry ***geometry, CConfig **config,
                                      unsigned long iExtIter, unsigned short val_nZone) {
  
  unsigned short iZone, iInst, nInst;
  
  for (iZone = 0; iZone < val_nZone; iZone++) {
    
    nInst = config[iZone]->GetnTimeInstances();

    for (iInst = 0; iInst < nInst; iInst++) {

      config[iZone]->SetiInst(iInst);

      /*--- Flags identifying the types of files to be written. ---*/

      bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
      if (config[iZone]->GetKind_SU2() == SU2_DOT) { Wrt_Vol = false; }
      bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();

      /*--- Get the file output format ---*/

      unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();

      /*--- Merge the node coordinates and connectivity if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/

      if ((Wrt_Vol || Wrt_Srf)) {
        if (rank == MASTER_NODE) cout << "Merging connectivities in the Master node." << endl;
        MergeConnectivity(config[iZone], geometry[iZone][iInst], iZone);
      }

      /*--- Merge the solution data needed for volume solutions and restarts ---*/

      if ((Wrt_Vol || Wrt_Srf)) {
        if (rank == MASTER_NODE) cout << "Merging solution in the Master node." << endl;
        MergeBaselineSolution(config[iZone], geometry[iZone][iInst], solver[iZone][iInst], iZone);
      }

      /*--- Write restart, Tecplot or Paraview files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/


      if (rank == MASTER_NODE) {

        if (Wrt_Vol) {

          switch (FileFormat) {

          case TECPLOT:

            /*--- Write a Tecplot ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (volume grid)." << endl;
            SetTecplotASCII(config[iZone], geometry[iZone][iInst], &solver[iZone][iInst], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
            break;

          case FIELDVIEW:

            /*--- Write a FieldView ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing FieldView ASCII file (volume grid)." << endl;
            SetFieldViewASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone);
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
            break;

          case TECPLOT_BINARY:

            /*--- Write a Tecplot binary solution file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot Binary file (volume grid)." << endl;
            SetTecplotBinary_DomainMesh(config[iZone], geometry[iZone][iInst], iZone);
            SetTecplotBinary_DomainSolution(config[iZone], geometry[iZone][iInst], iZone);
            break;

          case FIELDVIEW_BINARY:

            /*--- Write a binary binary file ---*/

            if (rank == MASTER_NODE) cout << "Writing FieldView ASCII file (volume grid)." << endl;
            SetFieldViewBinary(config[iZone], geometry[iZone][iInst], iZone, val_nZone);
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
            break;

          case PARAVIEW:

            /*--- Write a Paraview ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Paraview ASCII file (volume grid)." << endl;
            SetParaview_ASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
            break;

          default:
            break;
          }

        }

        if (Wrt_Srf) {

          switch (FileFormat) {

          case TECPLOT:

            /*--- Write a Tecplot ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (surface grid)." << endl;
            SetTecplotASCII(config[iZone], geometry[iZone][iInst], &solver[iZone][iInst], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], true);
            break;

          case TECPLOT_BINARY:

            /*--- Write a Tecplot binary solution file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot Binary file (surface grid)." << endl;
            SetTecplotBinary_SurfaceMesh(config[iZone], geometry[iZone][iInst], iZone);
            SetTecplotBinary_SurfaceSolution(config[iZone], geometry[iZone][iInst], iZone);
            break;

          case PARAVIEW:

            /*--- Write a Paraview ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Paraview ASCII file (surface grid)." << endl;
            SetParaview_ASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], true);
            break;

          default:
            break;
          }
        }

        if (FileFormat == TECPLOT_BINARY) {
          if (!wrote_base_file)
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
          if (!wrote_surf_file)
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], wrote_surf_file);
        }

        if (Wrt_Vol || Wrt_Srf)
          DeallocateSolution(config[iZone], geometry[iZone][iInst]);
      }



      /*--- Final broadcast (informing other procs that the base output
     file was written). ---*/

#ifdef HAVE_MPI
      SU2_MPI::Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif

    }

  }
}

void COutput::SetMesh_Files(CGeometry **geometry, CConfig **config, unsigned short val_nZone, bool new_file, bool su2_file) {

  char cstr[MAX_STRING_SIZE], out_file[MAX_STRING_SIZE];
  unsigned short iZone;
  ofstream output_file;
  string str;

  /*--- Read the name of the output and input file ---*/

  if (su2_file) {
  if (rank == MASTER_NODE) {
      str = config[ZONE_0]->GetMesh_Out_FileName();
      strcpy (out_file, str.c_str());
      strcpy (cstr, out_file);
      output_file.precision(15);
      output_file.open(cstr, ios::out);
      if (val_nZone > 1){
        output_file << "NZONE= " << val_nZone << endl;
      }
    }
  }

  for (iZone = 0; iZone < val_nZone; iZone++) {
    
    /*--- Flags identifying the types of files to be written. ---*/
    
    bool Wrt_Vol = config[iZone]->GetVisualize_Volume_Def();
    bool Wrt_Srf = config[iZone]->GetVisualize_Surface_Def();
    bool Wrt_Crd = config[iZone]->GetWrt_Crd_Sol();
    
    /*--- Merge the node coordinates and connectivity if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/
    
    if (rank == MASTER_NODE) cout <<"Merging grid connectivity." << endl;
    MergeConnectivity(config[iZone], geometry[iZone], iZone);
    
    /*--- Merge coordinates of all grid nodes (excluding ghost points).
     The grid coordinates are always merged and included first in the
     restart files. ---*/
    
    if (rank == MASTER_NODE) cout <<"Merging grid coordinates." << endl;
    MergeCoordinates(config[iZone], geometry[iZone]);
    
    /*--- Write restart, Tecplot or Paraview files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/
    
    if (rank == MASTER_NODE) {
      
      if (Wrt_Vol) {
        
        if (rank == MASTER_NODE) cout <<"Writing volume mesh file." << endl;
        
        /*--- Write a Tecplot ASCII file ---*/

        if (config[iZone]->GetOutput_FileFormat() == PARAVIEW) SetParaview_MeshASCII(config[iZone], geometry[iZone], iZone,  val_nZone, false,new_file);
        else if (config[iZone]->GetOutput_FileFormat() == PARAVIEW_BINARY) {
          if (rank == MASTER_NODE) cout <<"Writing ASCII paraview volume mesh file by default." << endl;
          SetParaview_MeshASCII(config[iZone], geometry[iZone], iZone,  val_nZone, false, new_file);
        }
        else if (config[iZone]->GetOutput_FileFormat() == TECPLOT) SetTecplotASCII_Mesh(config[iZone], geometry[iZone], iZone, false, new_file);
        else if (config[iZone]->GetOutput_FileFormat() == TECPLOT_BINARY) {
          if (rank == MASTER_NODE) cout <<"Writing ASCII tecplot volume mesh file by default." << endl;
          SetTecplotASCII_Mesh(config[iZone], geometry[iZone], iZone, false, new_file);
        }
        
      }
      
      if (Wrt_Srf) {
        
        if (rank == MASTER_NODE) cout <<"Writing surface mesh file." << endl;
        
        /*--- Write a Tecplot ASCII file ---*/
        
        if (config[iZone]->GetOutput_FileFormat()==PARAVIEW) SetParaview_MeshASCII(config[iZone], geometry[iZone], iZone,  val_nZone, true,new_file);
        else if (config[iZone]->GetOutput_FileFormat() == PARAVIEW_BINARY) {
          if (rank == MASTER_NODE) cout <<"Writing ASCII paraview surface mesh file by default." << endl;
          SetParaview_MeshASCII(config[iZone], geometry[iZone], iZone,  val_nZone, true, new_file);
        }
        else if (config[iZone]->GetOutput_FileFormat() == TECPLOT) SetTecplotASCII_Mesh(config[iZone], geometry[iZone], iZone, true, new_file);
        else if (config[iZone]->GetOutput_FileFormat() == TECPLOT_BINARY) {
          if (rank == MASTER_NODE) cout <<"Writing ASCII tecplot surface mesh file by default." << endl;
          SetTecplotASCII_Mesh(config[iZone], geometry[iZone], iZone, true, new_file);
        }
        
      }
      
      /*--- Write a .su2 ASCII file ---*/

      if (su2_file) {

      if (rank == MASTER_NODE) cout <<"Writing .su2 file." << endl;

        SetSU2_MeshASCII(config[iZone], geometry[iZone], iZone, output_file);
      
      /*--- Write an stl surface file ---*/

      if (rank == MASTER_NODE) cout <<"Writing .stl surface file." << endl;
      
        SetSTL_MeshASCII(config[iZone], geometry[iZone]);
        
      }
      
      /*--- Write a binary file with the grid coordinates alone. ---*/
      
      if (Wrt_Crd) {
        if (rank == MASTER_NODE) cout <<"Writing .dat binary coordinates file." << endl;
        WriteCoordinates_Binary(config[iZone], geometry[iZone], iZone);
      }

      
      /*--- Deallocate connectivity ---*/
      
      DeallocateConnectivity(config[iZone], geometry[iZone], true);
      DeallocateConnectivity(config[iZone], geometry[iZone], false);
      DeallocateCoordinates(config[iZone], geometry[iZone]);
      
    }

    /*--- Final broadcast (informing other procs that the base output
     file was written). ---*/
    
#ifdef HAVE_MPI
    SU2_MPI::Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif
    
    /*--- Write an csv surface file, done in parallel ---*/

    if (rank == MASTER_NODE) cout <<"Writing .csv surface file." << endl;

    if (su2_file) SetCSV_MeshASCII(config[iZone], geometry[iZone]);

  }

  if (rank == MASTER_NODE) {
    if (su2_file){
      output_file.close();
    }
  }
}

void COutput::SetSensitivity_Files(CGeometry ***geometry, CConfig **config, unsigned short val_nZone) {

  unsigned short iMarker,iDim, nDim, iVar, nMarker, nVar;
  unsigned long iVertex, iPoint, nPoint, nVertex;
  su2double *Normal, Prod, Sens = 0.0, SensDim, Area;

  unsigned short iZone;

  CSolver ***solver = new CSolver**[val_nZone];
  for (iZone = 0; iZone < val_nZone; iZone++) {
    solver[iZone] = new CSolver*[1];
  }

  for (iZone = 0; iZone < val_nZone; iZone++) {

    nPoint = geometry[iZone][INST_0]->GetnPoint();
    nDim   = geometry[iZone][INST_0]->GetnDim();
    nMarker = config[iZone]->GetnMarker_All();
    nVar = nDim + 1;

    /*--- We create a baseline solver to easily merge the sensitivity information ---*/

    vector<string> fieldnames;
    fieldnames.push_back("\"Point\"");
    fieldnames.push_back("\"x\"");
    fieldnames.push_back("\"y\"");
    if (nDim == 3) {
      fieldnames.push_back("\"z\"");
    }
    fieldnames.push_back("\"Sensitivity_x\"");
    fieldnames.push_back("\"Sensitivity_y\"");
    if (nDim == 3) {
      fieldnames.push_back("\"Sensitivity_z\"");
    }
    fieldnames.push_back("\"Surface_Sensitivity\"");

    solver[iZone][INST_0] = new CBaselineSolver(geometry[iZone][INST_0], config[iZone], nVar+nDim, fieldnames);

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        solver[iZone][INST_0]->node[iPoint]->SetSolution(iDim, geometry[iZone][INST_0]->node[iPoint]->GetCoord(iDim));
      }
      for (iVar = 0; iVar < nDim; iVar++) {
        solver[iZone][INST_0]->node[iPoint]->SetSolution(iVar+nDim, geometry[iZone][INST_0]->GetSensitivity(iPoint, iVar));
      }
    }

    /*--- Compute the sensitivity in normal direction ---*/

    for (iMarker = 0; iMarker < nMarker; iMarker++) {

      if((config[iZone]->GetMarker_All_KindBC(iMarker) == HEAT_FLUX )   ||
         (config[iZone]->GetMarker_All_KindBC(iMarker) == EULER_WALL )  ||
         (config[iZone]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL )  ||
         (config[iZone]->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE )) {
        
        nVertex = geometry[iZone][INST_0]->GetnVertex(iMarker);

        for (iVertex = 0; iVertex < nVertex; iVertex++) {
          iPoint = geometry[iZone][INST_0]->vertex[iMarker][iVertex]->GetNode();
          Normal = geometry[iZone][INST_0]->vertex[iMarker][iVertex]->GetNormal();
          Prod = 0.0;
          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {

            /*--- Retrieve the gradient calculated with discrete adjoint method ---*/

            SensDim = geometry[iZone][INST_0]->GetSensitivity(iPoint, iDim);

            /*--- Calculate scalar product for projection onto the normal vector ---*/

            Prod += Normal[iDim]*SensDim;

            Area += Normal[iDim]*Normal[iDim];
          }

          Area = sqrt(Area);

          /*--- Projection of the gradient onto the normal vector of the surface ---*/

          Sens = Prod/Area;

          solver[iZone][INST_0]->node[iPoint]->SetSolution(2*nDim, Sens);

        }
      }
    }
  }

  /*--- Merge the information and write the output files ---*/

  SetBaselineResult_Files(solver, geometry, config, 0, val_nZone);

  for (iZone = 0; iZone < val_nZone; iZone++) {
    delete solver[iZone][0];
    delete solver[iZone];
  }
  delete [] solver;
}


void COutput::SetResult_Files_Parallel(CSolver *****solver_container,
                                       CGeometry ****geometry,
                                       CConfig **config,
                                       unsigned long iExtIter,
                                       unsigned short iZone,
                                       unsigned short val_nZone) {
  
  unsigned short iVar, iInst;
  unsigned long iPoint;
  bool compressible = true;
  
  unsigned short nInst = config[iZone]->GetnTimeInstances();

    /*--- Get the file output format ---*/
    
    unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();
    
    for (iInst = 0; iInst < nInst; iInst++){

      config[iZone]->SetiInst(iInst);
      
      bool cont_adj = config[iZone]->GetContinuous_Adjoint();
      bool disc_adj = config[iZone]->GetDiscrete_Adjoint();

      /*--- Flags identifying the types of files to be written. ---*/
      /*--- For now, we are disabling the parallel writers for Tecplot
          ASCII until we have parallel versions of all file formats
          available. SU2_SOL will remain intact for writing files
          until this capability is completed. ---*/

      bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
      bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();
      bool Wrt_Csv = config[iZone]->GetWrt_Csv_Sol();

#ifdef HAVE_MPI
      /*--- Do not merge the connectivity or write the visualization files
     if we are running in parallel, unless we are using ParaView binary.
       Force the use of SU2_SOL to merge and write the viz. files in this
       case to save overhead. ---*/

      if ((size > SINGLE_NODE) && (FileFormat != PARAVIEW_BINARY)) {
        Wrt_Vol = false;
        Wrt_Srf = false;
      }
#endif

    /*--- Check for compressible/incompressible flow problems. ---*/

    compressible = (config[iZone]->GetKind_Regime() == COMPRESSIBLE);
    
    /*--- First, prepare the offsets needed throughout below. ---*/    
    
    PrepareOffsets(config[iZone], geometry[iZone][iInst][MESH_0]);

//    /*--- Write out CSV files in parallel for flow and adjoint. ---*/
    
//    if (rank == MASTER_NODE) cout << endl << "Writing comma-separated values (CSV) surface files." << endl;
    
//    switch (config[iZone]->GetKind_Solver()) {
//      case EULER : case NAVIER_STOKES : case RANS :
//        if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][iInst][MESH_0],
//            solver_container[iZone][iInst][MESH_0][FLOW_SOL], iExtIter, iZone, iInst);
//        break;
//      case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS :
//      case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
//        if (Wrt_Csv) SetSurfaceCSV_Adjoint(config[iZone], geometry[iZone][iInst][MESH_0],
//            solver_container[iZone][iInst][MESH_0][ADJFLOW_SOL],
//            solver_container[iZone][iInst][MESH_0][FLOW_SOL], iExtIter, iZone, iInst);
//        break;
//      default: break;
//    }

    /*--- Write a template inlet profile file if requested. ---*/

    if (config[iZone]->GetWrt_InletFile()) {
      MergeInletCoordinates(config[iZone], geometry[iZone][iInst][MESH_0]);

      if (rank == MASTER_NODE) {
        Write_InletFile_Flow(config[iZone], geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0]);
        DeallocateInletCoordinates(config[iZone], geometry[iZone][iInst][MESH_0]);
      }
      config[iZone]->SetWrt_InletFile(false);
    }

    /*--- This switch statement will become a call to a virtual function
     defined within each of the "physics" output child classes that loads
     the local data for that particular problem alone. ---*/

      if (rank == MASTER_NODE)
        cout << endl << "Loading solution output data locally on each rank." << endl;

      CollectVolumeData(config[iZone], geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0]);
         
    /*--- Store the solution to be used on the final iteration with cte. lift mode. ---*/

    if ((!cont_adj) && (!disc_adj) && (config[iZone]->GetFixed_CL_Mode()) &&
        (config[iZone]->GetnExtIter()-config[iZone]->GetIter_dCL_dAlpha() -1 == iExtIter)) {

      if (rank == MASTER_NODE)
        cout << "Storing solution output data locally on each rank (cte. CL mode)." << endl;
      
      Local_Data_Copy = new su2double*[geometry[iZone][iInst][MESH_0]->GetnPoint()];
      for (iPoint = 0; iPoint < geometry[iZone][iInst][MESH_0]->GetnPoint(); iPoint++) {
        Local_Data_Copy[iPoint] = new su2double[GlobalField_Counter];
      }
      
      for (iPoint = 0; iPoint < geometry[iZone][iInst][MESH_0]->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
          Local_Data_Copy[iPoint][iVar] = Local_Data[iPoint][iVar];
        }
      }
      
    }
    
    /*--- Recover the solution to be used on the final iteration with cte. lift mode. ---*/

    if ((!cont_adj) && (!disc_adj) && (config[iZone]->GetFixed_CL_Mode()) &&
        (config[iZone]->GetnExtIter() - 1 == iExtIter) && (Local_Data_Copy != NULL)) {

      if (rank == MASTER_NODE)
        cout << "Recovering solution output data locally on each rank (cte. CL mode)." << endl;
      
      for (iPoint = 0; iPoint < geometry[iZone][iInst][MESH_0]->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
          Local_Data[iPoint][iVar] = Local_Data_Copy[iPoint][iVar];
        }
      }
      
      for (iPoint = 0; iPoint < geometry[iZone][iInst][MESH_0]->GetnPoint(); iPoint++)
        delete [] Local_Data_Copy[iPoint];
      delete [] Local_Data_Copy;
      
    }
    
    /*--- After loading the data local to a processor, we perform a sorting,
     i.e., a linear partitioning of the data across all ranks in the communicator. ---*/
    
    if (rank == MASTER_NODE) cout << "Sorting output data across all ranks." << endl;
    SortOutputData(config[iZone], geometry[iZone][iInst][MESH_0]);
    
    /*--- Write either a binary or ASCII restart file in parallel. ---*/

    if (config[iZone]->GetWrt_Binary_Restart()) {
      if (rank == MASTER_NODE) cout << "Writing binary SU2 native restart file." << endl;
      WriteRestart_Parallel_Binary(config[iZone], geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], iZone, iInst);
    } else {
      if (rank == MASTER_NODE) cout << "Writing ASCII SU2 native restart file." << endl;
      WriteRestart_Parallel_ASCII(config[iZone], geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], iZone, iInst);
    }

    /*--- Write a slice on a structured mesh if requested. ---*/

    if (config[iZone]->GetWrt_Slice()) {
      WriteCSV_Slice(config[iZone], geometry[iZone][iInst][MESH_0],
                     solver_container[iZone][iInst][MESH_0][FLOW_SOL], iExtIter, iZone, 0);
      WriteCSV_Slice(config[iZone], geometry[iZone][iInst][MESH_0],
                     solver_container[iZone][iInst][MESH_0][FLOW_SOL], iExtIter, iZone, 1);
    }

    /*--- Write the solution files if they are requested and we are executing
     with a single rank (all data on one proc and no comm. overhead). Once we
     have parallel binary versions of Tecplot / ParaView / CGNS / etc., we
     can allow the write of the viz. files as well. ---*/

    if ((Wrt_Vol || Wrt_Srf)) {
      
      /*--- First, sort all connectivity into linearly partitioned chunks of elements. ---*/

      if (rank == MASTER_NODE)
        cout << "Preparing element connectivity across all ranks." << endl;
      SortConnectivity(config[iZone], geometry[iZone][iInst][MESH_0], iZone);

      /*--- Sort the surface data and renumber if for writing. ---*/

      SortOutputData_Surface(config[iZone], geometry[iZone][iInst][MESH_0]);

      /*--- Write Tecplot/ParaView ASCII files for the volume and/or surface solutions. ---*/

      if (Wrt_Vol) {

        switch (FileFormat) {

          case TECPLOT:

            /*--- Write a Tecplot ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file volume solution file." << endl;
            WriteTecplotASCII_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
                solver_container[iZone][iInst][MESH_0], iZone, val_nZone, iInst, nInst, false);
            break;

          case FIELDVIEW:

            /*--- We do not yet have a version of FieldView ASCII for new parallel output. ---*/

            if (rank == MASTER_NODE) cout << "FieldView ASCII volume files not available in serial with SU2_CFD." << endl;
            if (rank == MASTER_NODE) cout << "  Run SU2_SOL to generate FieldView ASCII." << endl;

            break;

          case TECPLOT_BINARY:

            /*--- Write a Tecplot ASCII file instead for now in serial. ---*/

            if (rank == MASTER_NODE) cout << "Tecplot binary volume files not available in serial with SU2_CFD." << endl;
            if (rank == MASTER_NODE) cout << "  Run SU2_SOL to generate Tecplot binary." << endl;
            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file volume solution file instead." << endl;
            WriteTecplotASCII_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
                solver_container[iZone][iInst][MESH_0], iZone, val_nZone, iInst, nInst, false);
            break;

          case FIELDVIEW_BINARY:

            /*--- FieldView binary files not yet available for parallel output. ---*/

            if (rank == MASTER_NODE) cout << "FieldView ASCII volume files not available in serial with SU2_CFD." << endl;
            if (rank == MASTER_NODE) cout << "  Run SU2_SOL to generate FieldView ASCII." << endl;
            break;

          case PARAVIEW:

            /*--- Write a Paraview ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Paraview ASCII volume solution file." << endl;
            WriteParaViewASCII_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
                solver_container[iZone][iInst][MESH_0], iZone, val_nZone, iInst, nInst, false);
            break;
            
          case PARAVIEW_BINARY:
            
            /*--- Write a Paraview binary file ---*/
            
            if (rank == MASTER_NODE) cout << "Writing Paraview binary volume solution file." << endl;
            WriteParaViewBinary_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
                                        solver_container[iZone][iInst][MESH_0], iZone, val_nZone, nInst, false);
            break;

          default:
            break;
          }

        }

        if (Wrt_Srf) {

          switch (FileFormat) {

          case TECPLOT:

            /*--- Write a Tecplot ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file surface solution file." << endl;
            WriteTecplotASCII_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
                solver_container[iZone][iInst][MESH_0], iZone, val_nZone, iInst, nInst, true);
            break;

          case TECPLOT_BINARY:

            /*--- Write a Tecplot ASCII file instead for now in serial. ---*/

            if (rank == MASTER_NODE) cout << "Tecplot binary surface files not available in serial with SU2_CFD." << endl;
            if (rank == MASTER_NODE) cout << "  Run SU2_SOL to generate Tecplot binary." << endl;
            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file surface solution file instead." << endl;
            WriteTecplotASCII_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
                solver_container[iZone][iInst][MESH_0], iZone, val_nZone, iInst, nInst, true);
            break;

          case PARAVIEW:

            /*--- Write a Paraview ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Paraview ASCII surface solution file." << endl;
            WriteParaViewASCII_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
                solver_container[iZone][iInst][MESH_0], iZone, val_nZone, iInst, nInst, true);
            break;

          case PARAVIEW_BINARY:
            
            /*--- Write a Paraview binary file ---*/
            
            if (rank == MASTER_NODE) cout << "Writing Paraview binary surface solution file." << endl;
            WriteParaViewBinary_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
                                         solver_container[iZone][iInst][MESH_0], iZone, val_nZone, nInst, true);
            break;
            

          default:
            break;
          }

        }

        /*--- Clean up the connectivity data that was allocated for output. ---*/

        DeallocateConnectivity_Parallel(config[iZone], geometry[iZone][iInst][MESH_0], false);
        DeallocateConnectivity_Parallel(config[iZone], geometry[iZone][iInst][MESH_0], true);

        /*--- Clean up the surface data that was only needed for output. ---*/

        DeallocateSurfaceData_Parallel(config[iZone], geometry[iZone][iInst][MESH_0]);

      }

      /*--- Deallocate the nodal data needed for writing restarts. ---*/

      DeallocateData_Parallel(config[iZone], geometry[iZone][iInst][MESH_0]);

    }

}

void COutput::SortConnectivity(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {

  /*--- Flags identifying the types of files to be written. ---*/
  
  bool Wrt_Vol = config->GetWrt_Vol_Sol();
  bool Wrt_Srf = config->GetWrt_Srf_Sol();
  
  /*--- Sort connectivity for each type of element (excluding halos). Note
   In these routines, we sort the connectivity into a linear partitioning
   across all processors based on the global index of the grid nodes. ---*/
  
  /*--- Sort volumetric grid connectivity. ---*/
  
  if (Wrt_Vol) {
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
      cout <<"Sorting volumetric grid connectivity." << endl;
    
    SortVolumetricConnectivity(config, geometry, TRIANGLE     );
    SortVolumetricConnectivity(config, geometry, QUADRILATERAL);
    SortVolumetricConnectivity(config, geometry, TETRAHEDRON  );
    SortVolumetricConnectivity(config, geometry, HEXAHEDRON   );
    SortVolumetricConnectivity(config, geometry, PRISM        );
    SortVolumetricConnectivity(config, geometry, PYRAMID      );
    
  }
  
  /*--- Sort surface grid connectivity. ---*/
  
  if (Wrt_Srf) {
    
    if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
      cout <<"Sorting surface grid connectivity." << endl;
    
    SortSurfaceConnectivity(config, geometry, LINE         );
    SortSurfaceConnectivity(config, geometry, TRIANGLE     );
    SortSurfaceConnectivity(config, geometry, QUADRILATERAL);
    
  }
  
  /*--- Reduce the total number of cells we will be writing in the output files. ---*/
  
  unsigned long nTotal_Elem = nParallel_Tria + nParallel_Quad + nParallel_Tetr + nParallel_Hexa + nParallel_Pris + nParallel_Pyra;
  unsigned long nTotal_Surf_Elem = nParallel_Line + nParallel_BoundTria + nParallel_BoundQuad;
#ifndef HAVE_MPI
  nGlobal_Elem_Par = nTotal_Elem;
  nSurf_Elem_Par   = nTotal_Surf_Elem;
#else
  SU2_MPI::Allreduce(&nTotal_Elem, &nGlobal_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nTotal_Surf_Elem, &nSurf_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
}

void COutput::SortVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  unsigned long iProcessor;
  unsigned short NODES_PER_ELEMENT = 0;
  unsigned long iPoint, jPoint, kPoint, nLocalPoint, nTotalPoint;
  unsigned long nElem_Total = 0, Global_Index;
  
  unsigned long iVertex, iMarker;
  int SendRecv, RecvFrom;
  
  bool notPeriodic, notHalo, addedPeriodic, isPeriodic;
  
  int *Local_Halo = NULL;
  int *Conn_Elem  = NULL;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  switch (Elem_Type) {
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case PRISM:
      NODES_PER_ELEMENT = N_POINTS_PRISM;
      break;
    case PYRAMID:
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
  }
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that were part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  
  if (config->GetKind_SU2() != SU2_DEF) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
              (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
              (SendRecv < 0)) {
            Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
          }
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found. ---*/
  
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];
  
  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++)
    Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain.
         For SU2_DEF we want to keep all periodic nodes. ---*/
        
        if (config->GetKind_SU2() == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }
        
        notPeriodic = (isPeriodic && (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < (unsigned long)size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (Global_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          
          kPoint = (iProcessor+1)*maxAddedPeriodic;
          
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
        
      }
    }
  }
  
  /*--- Now that we've done the gymnastics to find any periodic points,
   compute the total number of local and global points for the output. ---*/
  
  nLocalPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif
  
  /*--- Compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  unsigned long *nPoint_Linear = new unsigned long[size+1];
  
  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }
  
  /*--- Get the number of remainder points after the even division. ---*/
  
  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }
  
  /*--- Store the local number of nodes and the beginning/end index ---*/
  
  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[size] = nTotalPoint;
  
  /*--- We start with the connectivity distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many elements we must send to each other rank in order to
   have all elements sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first nPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/
  
  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];
  
  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;
  
  for (int ii = 0; ii < (int)geometry->GetnElem(); ii++ ) {
    if (geometry->elem[ii]->GetVTK_Type() == Elem_Type) {
      for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {
        
        /*--- Get the index of the current point. ---*/
        
        iPoint = geometry->elem[ii]->GetNode(jj);
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- Search for the lowest global index in this element. We
         send the element to the processor owning the range that includes
         the lowest global index value. ---*/
        
        for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
          jPoint = geometry->elem[ii]->GetNode(kk);
          unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
          if (newID < Global_Index) Global_Index = newID;
        }
        
        /*--- Search for the processor that owns this point ---*/
        
        iProcessor = Global_Index/npoint_procs[0];
        if (iProcessor >= (unsigned long)size)
          iProcessor = (unsigned long)size-1;
        if (Global_Index >= nPoint_Linear[iProcessor])
          while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
        else
          while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
        
        /*--- If we have not visited this element yet, increment our
         number of elements that must be sent to a particular proc. ---*/
        
        if ((nElem_Flag[iProcessor] != ii)) {
          nElem_Flag[iProcessor] = ii;
          nElem_Send[iProcessor+1]++;
        }
        
      }
    }
  }
  
  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif
  
  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
  
  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;
    
    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }
  
  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/
  
  unsigned long *connSend = NULL;
  connSend = new unsigned long[NODES_PER_ELEMENT*nElem_Send[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Send[size]; ii++)
    connSend[ii] = 0;
  
  /*--- Allocate arrays for storing halo flags. ---*/
  
  unsigned short *haloSend = new unsigned short[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    haloSend[ii] = false;
  
  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/
  
  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = NODES_PER_ELEMENT*nElem_Send[ii];
  
  unsigned long *haloIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) haloIndex[ii] = nElem_Send[ii];
  
  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/
  
  for (int ii = 0; ii < (int)geometry->GetnElem(); ii++) {
    if (geometry->elem[ii]->GetVTK_Type() == Elem_Type) {
      for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {
        
        /*--- Get the index of the current point. ---*/
        
        iPoint = geometry->elem[ii]->GetNode(jj);
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- Search for the lowest global index in this element. We
         send the element to the processor owning the range that includes
         the lowest global index value. ---*/
        
        for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
          jPoint = geometry->elem[ii]->GetNode(kk);
          unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
          if (newID < Global_Index) Global_Index = newID;
        }
        
        /*--- Search for the processor that owns this point ---*/
        
        iProcessor = Global_Index/npoint_procs[0];
        if (iProcessor >= (unsigned long)size)
          iProcessor = (unsigned long)size-1;
        if (Global_Index >= nPoint_Linear[iProcessor])
          while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
        else
          while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
        
        /*--- Load connectivity into the buffer for sending ---*/
        
        if (nElem_Flag[iProcessor] != ii) {
          
          nElem_Flag[iProcessor] = ii;
          unsigned long nn = index[iProcessor];
          unsigned long mm = haloIndex[iProcessor];
          
          /*--- Load the connectivity values. ---*/
          
          for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
            iPoint = geometry->elem[ii]->GetNode(kk);
            connSend[nn] = geometry->node[iPoint]->GetGlobalIndex(); nn++;
            
            /*--- Check if this is a halo node. If so, flag this element
             as a halo cell. We will use this later to sort and remove
             any duplicates from the connectivity list. ---*/
            
            if (Local_Halo[iPoint]) haloSend[mm] = true;
            
          }
          
          /*--- Increment the index by the message length ---*/
          
          index[iProcessor]    += NODES_PER_ELEMENT;
          haloIndex[iProcessor]++;
          
        }
      }
    }
  }
  
  /*--- Free memory after loading up the send buffer. ---*/
  
  delete [] index;
  delete [] haloIndex;
  
  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  unsigned long *connRecv = NULL;
  connRecv = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Recv[size]; ii++)
    connRecv[ii] = 0;
  
  unsigned short *haloRecv = new unsigned short[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    haloRecv[ii] = false;
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/
  
  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];
  
  /*--- Launch the non-blocking recv's for the connectivity. ---*/
  
  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = NODES_PER_ELEMENT*nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = NODES_PER_ELEMENT*nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Repeat the process to communicate the halo flags. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(haloRecv[ll]), count, MPI_UNSIGNED_SHORT, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the halo flags. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(haloSend[ll]), count, MPI_UNSIGNED_SHORT, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif
  
  /*--- Copy my own rank's data into the recv buffer directly. ---*/
  
  int mm = NODES_PER_ELEMENT*nElem_Recv[rank];
  int ll = NODES_PER_ELEMENT*nElem_Send[rank];
  int kk = NODES_PER_ELEMENT*nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];
  
  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) haloRecv[mm] = haloSend[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. Note that we add 1 here
   to the connectivity for vizualization packages. First, allocate
   appropriate amount of memory for this section. ---*/
  
  if (nElem_Recv[size] > 0) Conn_Elem = new int[NODES_PER_ELEMENT*nElem_Recv[size]];
  int count = 0; nElem_Total = 0;
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    if (!haloRecv[ii]) {
      nElem_Total++;
      for (int jj = 0; jj < NODES_PER_ELEMENT; jj++) {
        Conn_Elem[count] = (int)connRecv[ii*NODES_PER_ELEMENT+jj] + 1;
        count++;
      }
    }
  }
  
  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
  switch (Elem_Type) {
    case TRIANGLE:
      nParallel_Tria = nElem_Total;
      if (nParallel_Tria > 0) Conn_Tria_Par = Conn_Elem;
      break;
    case QUADRILATERAL:
      nParallel_Quad = nElem_Total;
      if (nParallel_Quad > 0) Conn_Quad_Par = Conn_Elem;
      break;
    case TETRAHEDRON:
      nParallel_Tetr = nElem_Total;
      if (nParallel_Tetr > 0) Conn_Tetr_Par = Conn_Elem;
      break;
    case HEXAHEDRON:
      nParallel_Hexa = nElem_Total;
      if (nParallel_Hexa > 0) Conn_Hexa_Par = Conn_Elem;
      break;
    case PRISM:
      nParallel_Pris = nElem_Total;
      if (nParallel_Pris > 0) Conn_Pris_Par = Conn_Elem;
      break;
    case PYRAMID:
      nParallel_Pyra = nElem_Total;
      if (nParallel_Pyra > 0) Conn_Pyra_Par = Conn_Elem;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      break;
  }
  
  /*--- Free temporary memory from communications ---*/
  
  delete [] connSend;
  delete [] connRecv;
  delete [] haloSend;
  delete [] haloRecv;
  delete [] Local_Halo;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic; 
  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear;

}

void COutput::SortSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  unsigned long iProcessor;
  unsigned short NODES_PER_ELEMENT;
  unsigned long iPoint, jPoint, kPoint, nLocalPoint, nTotalPoint;
  unsigned long nElem_Total = 0, Global_Index;
  
  unsigned long iVertex, iMarker;
  int SendRecv, RecvFrom;
  
  bool notPeriodic, notHalo, addedPeriodic, isPeriodic;
  
  int *Local_Halo = NULL;
  int *Conn_Elem  = NULL;
  
#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  switch (Elem_Type) {
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      NODES_PER_ELEMENT = 0;
      break;
  }
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that were part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  
  if (config->GetKind_SU2() != SU2_DEF) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
              (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
              (SendRecv < 0)) {
            Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
          }
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found. ---*/
  
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];
  
  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++)
    Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain.
         For SU2_DEF we want to keep all periodic nodes. ---*/
        
        if (config->GetKind_SU2() == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }
        
        notPeriodic = (isPeriodic && (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < (unsigned long)size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (Global_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          
          kPoint = (iProcessor+1)*maxAddedPeriodic;
          
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
        
      }
    }
  }
  
  /*--- Now that we've done the gymnastics to find any periodic points,
   compute the total number of local and global points for the output. ---*/
  
  nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif
  
  /*--- Compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  unsigned long *nPoint_Linear = new unsigned long[size+1];
  
  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }
  
  /*--- Get the number of remainder points after the even division. ---*/
  
  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }
  
  /*--- Store the local number of nodes and the beginning/end index ---*/
  
  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[size] = nTotalPoint;
  
  /*--- We start with the connectivity distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many elements we must send to each other rank in order to
   have all elements sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first nPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/
  
  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];
  
  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      
      for (int ii = 0; ii < (int)geometry->GetnElem_Bound(iMarker); ii++) {
        
        if (geometry->bound[iMarker][ii]->GetVTK_Type() == Elem_Type) {
          for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {
            
            /*--- Get the index of the current point. ---*/
            
            iPoint = geometry->bound[iMarker][ii]->GetNode(jj);
            Global_Index = geometry->node[iPoint]->GetGlobalIndex();
            
            /*--- Search for the lowest global index in this element. We
             send the element to the processor owning the range that includes
             the lowest global index value. ---*/
            
            for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
              jPoint = geometry->bound[iMarker][ii]->GetNode(kk);
              unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
              if (newID < Global_Index) Global_Index = newID;
            }
            
            /*--- Search for the processor that owns this point ---*/
            
            iProcessor = Global_Index/npoint_procs[0];
            if (iProcessor >= (unsigned long)size)
              iProcessor = (unsigned long)size-1;
            if (Global_Index >= nPoint_Linear[iProcessor])
              while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
            else
              while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
            
            /*--- If we have not visited this element yet, increment our
             number of elements that must be sent to a particular proc. ---*/
            
            if ((nElem_Flag[iProcessor] != ii)) {
              nElem_Flag[iProcessor] = ii;
              nElem_Send[iProcessor+1]++;
            }
            
          }
        }
      }
    }
  }
  
  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif
  
  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
  
  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;
    
    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }
  
  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/
  
  unsigned long *connSend = NULL;
  connSend = new unsigned long[NODES_PER_ELEMENT*nElem_Send[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Send[size]; ii++)
    connSend[ii] = 0;
  
  /*--- Allocate arrays for storing halo flags. ---*/
  
  unsigned short *haloSend = new unsigned short[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    haloSend[ii] = false;
  
  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/
  
  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = NODES_PER_ELEMENT*nElem_Send[ii];
  
  unsigned long *haloIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) haloIndex[ii] = nElem_Send[ii];
  
  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      
      for (int ii = 0; ii < (int)geometry->GetnElem_Bound(iMarker); ii++) {
        
        if (geometry->bound[iMarker][ii]->GetVTK_Type() == Elem_Type) {
          for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {
            
            /*--- Get the index of the current point. ---*/
            
            iPoint = geometry->bound[iMarker][ii]->GetNode(jj);
            Global_Index = geometry->node[iPoint]->GetGlobalIndex();
            
            /*--- Search for the lowest global index in this element. We
             send the element to the processor owning the range that includes
             the lowest global index value. ---*/
            
            for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
              jPoint = geometry->bound[iMarker][ii]->GetNode(kk);
              unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
              if (newID < Global_Index) Global_Index = newID;
            }
            
            /*--- Search for the processor that owns this point ---*/
            
            iProcessor = Global_Index/npoint_procs[0];
            if (iProcessor >= (unsigned long)size)
              iProcessor = (unsigned long)size-1;
            if (Global_Index >= nPoint_Linear[iProcessor])
              while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
            else
              while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
            
            /*--- Load connectivity into the buffer for sending ---*/
            
            if (nElem_Flag[iProcessor] != ii) {
              
              nElem_Flag[iProcessor] = ii;
              unsigned long nn = index[iProcessor];
              unsigned long mm = haloIndex[iProcessor];
              
              /*--- Load the connectivity values. ---*/
              
              for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
                iPoint = geometry->bound[iMarker][ii]->GetNode(kk);
                connSend[nn] = geometry->node[iPoint]->GetGlobalIndex(); nn++;
                
                /*--- Check if this is a halo node. If so, flag this element
                 as a halo cell. We will use this later to sort and remove
                 any duplicates from the connectivity list. ---*/
                
                if (Local_Halo[iPoint]) haloSend[mm] = true;
                
              }
              
              /*--- Increment the index by the message length ---*/
              
              index[iProcessor]    += NODES_PER_ELEMENT;
              haloIndex[iProcessor]++;
              
            }
          }
        }
      }
    }
  }
  
  /*--- Free memory after loading up the send buffer. ---*/
  
  delete [] index;
  delete [] haloIndex;
  
  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  unsigned long *connRecv = NULL;
  connRecv = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Recv[size]; ii++)
    connRecv[ii] = 0;
  
  unsigned short *haloRecv = new unsigned short[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    haloRecv[ii] = false;
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/
  
  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];
  
  /*--- Launch the non-blocking recv's for the connectivity. ---*/
  
  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = NODES_PER_ELEMENT*nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = NODES_PER_ELEMENT*nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Repeat the process to communicate the halo flags. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(haloRecv[ll]), count, MPI_UNSIGNED_SHORT, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the halo flags. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(haloSend[ll]), count, MPI_UNSIGNED_SHORT, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif
  
  /*--- Copy my own rank's data into the recv buffer directly. ---*/
  
  int mm = NODES_PER_ELEMENT*nElem_Recv[rank];
  int ll = NODES_PER_ELEMENT*nElem_Send[rank];
  int kk = NODES_PER_ELEMENT*nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];
  
  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) haloRecv[mm] = haloSend[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. Note that we add 1 here
   to the connectivity for vizualization packages. First, allocate
   appropriate amount of memory for this section. ---*/
  
  if (nElem_Recv[size] > 0) Conn_Elem = new int[NODES_PER_ELEMENT*nElem_Recv[size]];
  int count = 0; nElem_Total = 0;
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    if (!haloRecv[ii]) {
      nElem_Total++;
      for (int jj = 0; jj < NODES_PER_ELEMENT; jj++) {
        Conn_Elem[count] = (int)connRecv[ii*NODES_PER_ELEMENT+jj] + 1;
        count++;
      }
    }
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
  switch (Elem_Type) {
    case LINE:
      nParallel_Line = nElem_Total;
      if (nParallel_Line > 0) Conn_BoundLine_Par = Conn_Elem;
      break;
    case TRIANGLE:
      nParallel_BoundTria = nElem_Total;
      if (nParallel_BoundTria > 0) Conn_BoundTria_Par = Conn_Elem;
      break;
    case QUADRILATERAL:
      nParallel_BoundQuad = nElem_Total;
      if (nParallel_BoundQuad > 0) Conn_BoundQuad_Par = Conn_Elem;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      break;
  }
  
  /*--- Free temporary memory from communications ---*/
  
  delete [] connSend;
  delete [] connRecv;
  delete [] haloSend;
  delete [] haloRecv;
  delete [] Local_Halo;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear;
  
}

void COutput::SortOutputData(CConfig *config, CGeometry *geometry) {
  
  unsigned short iMarker;
  unsigned long iProcessor;
  unsigned long iPoint, Global_Index, nLocalPoint, nTotalPoint, iVertex;
  
  int VARS_PER_POINT = GlobalField_Counter;
  int *Local_Halo = NULL;

  bool isPeriodic;
  
#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/
  
  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      
      /*--- Checking for less than or equal to the rank, because there may
       be some periodic halo nodes that send info to the same rank. ---*/
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                      (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        if (isPeriodic) Local_Halo[iPoint] = false;
      }
    }
  }
  
  /*--- Sum total number of nodes that belong to the domain ---*/
  
  nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif
  
  /*--- Now that we know the actual number of points we need to output,
   compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  unsigned long *nPoint_Linear = new unsigned long[size+1];
  
  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }
  
  /*--- Get the number of remainder points after the even division. ---*/
  
  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }
  
  /*--- Store the local number of nodes and the beginning/end index ---*/
  
  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[size] = nTotalPoint;
  
  /*--- We start with the grid nodes distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many nodes we must send to each other rank in order to
   have all nodes sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first ~ nGlobalPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/
  
  int *nPoint_Send = new int[size+1]; nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1]; nPoint_Recv[0] = 0;
  int *nPoint_Flag = new int[size];
  
  for (int ii=0; ii < size; ii++) {
    nPoint_Send[ii] = 0;
    nPoint_Recv[ii] = 0;
    nPoint_Flag[ii]= -1;
  }
  nPoint_Send[size] = 0; nPoint_Recv[size] = 0;
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++ ) {
    
    /*--- We only write interior points and recovered periodic points. ---*/
    
    if (!Local_Halo[iPoint]) {
      
      /*--- Get the global index of the current point. ---*/
      
      Global_Index = geometry->node[iPoint]->GetGlobalIndex();
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear[iProcessor])
        while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
      
      /*--- If we have not visited this node yet, increment our
       number of elements that must be sent to a particular proc. ---*/
      
      if (nPoint_Flag[iProcessor] != (int)iPoint) {
        nPoint_Flag[iProcessor] = (int)iPoint;
        nPoint_Send[iProcessor+1]++;
      }
      
    }
  }
  
  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nPoint_Recv[1] = nPoint_Send[1];
#endif
  
  /*--- Prepare to send coordinates. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
  
  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nPoint_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nPoint_Recv[ii+1] > 0)) nRecvs++;
    
    nPoint_Send[ii+1] += nPoint_Send[ii];
    nPoint_Recv[ii+1] += nPoint_Recv[ii];
  }
  
  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/
  
  su2double *connSend = NULL;
  connSend = new su2double[VARS_PER_POINT*nPoint_Send[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Send[size]; ii++)
    connSend[ii] = 0;
  
  /*--- Allocate arrays for sending the global ID. ---*/
  
  unsigned long *idSend = new unsigned long[nPoint_Send[size]];
  for (int ii = 0; ii < nPoint_Send[size]; ii++)
    idSend[ii] = 0;
  
  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/
  
  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = VARS_PER_POINT*nPoint_Send[ii];
  
  unsigned long *idIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) idIndex[ii] = nPoint_Send[ii];
  
  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- We only write interior points and recovered periodic points. ---*/
    
    if (!Local_Halo[iPoint]) {
      
      /*--- Get the index of the current point. ---*/
      
      Global_Index = geometry->node[iPoint]->GetGlobalIndex();
      
      /*--- Search for the processor that owns this point. ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear[iProcessor])
        while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;
      
      /*--- Load node coordinates into the buffer for sending. ---*/
      
      if (nPoint_Flag[iProcessor] != (int)iPoint) {
        
        nPoint_Flag[iProcessor] = (int)iPoint;
        unsigned long nn = index[iProcessor];
        
        /*--- Load the data values. ---*/
        
        for (unsigned short kk = 0; kk < VARS_PER_POINT; kk++) {
          connSend[nn] = Local_Data[iPoint][kk]; nn++;
        }
        
        /*--- Load the global ID (minus offset) for sorting the
         points once they all reach the correct processor. ---*/
        
        nn = idIndex[iProcessor];
        idSend[nn] = Global_Index - starting_node[iProcessor];
        
        /*--- Increment the index by the message length ---*/
        
        index[iProcessor]  += VARS_PER_POINT;
        idIndex[iProcessor]++;
        
      }
    }
  }
  
  /*--- Free memory after loading up the send buffer. ---*/
  
  delete [] index;
  delete [] idIndex;
  
  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  su2double *connRecv = NULL;
  connRecv = new su2double[VARS_PER_POINT*nPoint_Recv[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Recv[size]; ii++)
    connRecv[ii] = 0;
  
  unsigned long *idRecv = new unsigned long[nPoint_Recv[size]];
  for (int ii = 0; ii < nPoint_Recv[size]; ii++)
    idRecv[ii] = 0;
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the global IDs. ---*/
  
  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];
  
  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = VARS_PER_POINT*nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = VARS_PER_POINT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_DOUBLE, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = VARS_PER_POINT*nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = VARS_PER_POINT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_DOUBLE, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Repeat the process to communicate the global IDs. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the global IDs. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif
  
  /*--- Copy my own rank's data into the recv buffer directly. ---*/
  
  int mm = VARS_PER_POINT*nPoint_Recv[rank];
  int ll = VARS_PER_POINT*nPoint_Send[rank];
  int kk = VARS_PER_POINT*nPoint_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];
  
  mm = nPoint_Recv[rank];
  ll = nPoint_Send[rank];
  kk = nPoint_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. First, allocate the
   appropriate amount of memory for this section. ---*/
  
  Parallel_Data = new su2double*[VARS_PER_POINT];
  for (int jj = 0; jj < VARS_PER_POINT; jj++) {
    Parallel_Data[jj] = new su2double[nPoint_Recv[size]];
    for (int ii = 0; ii < nPoint_Recv[size]; ii++) {
      Parallel_Data[jj][idRecv[ii]] = connRecv[ii*VARS_PER_POINT+jj];
    }
  }
  
  /*--- Store the total number of local points my rank has for
   the current section after completing the communications. ---*/
  
  nParallel_Poin = nPoint_Recv[size];
  
  /*--- Reduce the total number of points we will write in the output files. ---*/

#ifndef HAVE_MPI
  nGlobal_Poin_Par = nParallel_Poin;
#else
  SU2_MPI::Allreduce(&nParallel_Poin, &nGlobal_Poin_Par, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Free temporary memory from communications ---*/
  
  delete [] connSend;
  delete [] connRecv;
  delete [] idSend;
  delete [] idRecv;
  delete [] nPoint_Recv;
  delete [] nPoint_Send;
  delete [] nPoint_Flag;

  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear;
  
}

void COutput::SortOutputData_Surface(CConfig *config, CGeometry *geometry) {
  
  unsigned short iMarker;
  unsigned long iProcessor;
  unsigned long iPoint, jPoint, kPoint, iElem;
  unsigned long Global_Index, nLocalPoint, nTotalPoint, iVertex;
  
  int VARS_PER_POINT = GlobalField_Counter;
  int *Local_Halo = NULL;
  int iNode, count;
  int SendRecv, RecvFrom;
  
  bool notPeriodic, notHalo, addedPeriodic, isPeriodic;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 1: We already have the surface connectivity spread out in     ---*/
  /*---         linear partitions across all procs and the output data     ---*/
  /*---         for the entire field is similarly linearly partitioned.    ---*/
  /*---         We need to identify which nodes in the volume data are     ---*/
  /*---         also surface points. Our first step is to loop over all    ---*/
  /*---         of the sorted surface connectivity and create a data       ---*/
  /*---         structure on each proc that can identify the local surf    ---*/
  /*---         points. Note that the linear partitioning is slightly      ---*/
  /*---         different between the nodes and elements, so we will       ---*/
  /*---         have to move between the two systems in this routine.      ---*/
  /*--------------------------------------------------------------------------*/
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. This is the linear partitioning for nodes. ---*/
  
  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      
      /*--- Checking for less than or equal to the rank, because there may
       be some periodic halo nodes that send info to the same rank. ---*/
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                      (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        if (isPeriodic) Local_Halo[iPoint] = false;
      }
    }
  }
  
  /*--- Sum total number of nodes that belong to the domain ---*/
  
  nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif
  
  /*--- Now that we know the actual number of points we need to output,
   compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  
  unsigned long *nPoint_Linear_Nodes = new unsigned long[size+1];
  unsigned long *nPoint_Linear_Elems = new unsigned long[size+1];
  
  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }
  
  /*--- Get the number of remainder points after the even division. ---*/
  
  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }
  
  /*--- Store the local number of nodes and the beginning/end index ---*/
  
  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear_Nodes[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear_Nodes[ii] = nPoint_Linear_Nodes[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear_Nodes[size] = nTotalPoint;
  
  /*--- Prepare to check and communicate the nodes that each proc has
   locally from the surface connectivity. ---*/
  
  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];
  
  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;
  
  /*--- Loop through our local line elements and check where each
   of the grid nodes resides based on global index. ---*/
  
  for (int ii = 0; ii < (int)nParallel_Line; ii++) {
    for ( int jj = 0; jj < N_POINTS_LINE; jj++ ) {
      
      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/
      
      iNode = ii*N_POINTS_LINE+jj;
      Global_Index = Conn_BoundLine_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
      
      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/
      
      if ((nElem_Flag[iProcessor] != iNode)) {
        nElem_Flag[iProcessor] = iNode;
        nElem_Send[iProcessor+1]++;
      }
      
    }
  }
  
  /*--- Reset out flags and then loop through our local triangle surface
   elements performing the same check for where each grid node resides. ---*/
  
  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;
  
  for (int ii = 0; ii < (int)nParallel_BoundTria; ii++) {
    for ( int jj = 0; jj < N_POINTS_TRIANGLE; jj++ ) {
      
      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/
      
      iNode = ii*N_POINTS_TRIANGLE + jj;
      Global_Index = Conn_BoundTria_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
      
      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/
      
      if ((nElem_Flag[iProcessor] != iNode)) {
        nElem_Flag[iProcessor] = iNode;
        nElem_Send[iProcessor+1]++;
      }
      
    }
  }
  
  /*--- Reset out flags and then loop through our local quad surface
   elements performing the same check for where each grid node resides. ---*/
  
  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;
  
  for (int ii = 0; ii < (int)nParallel_BoundQuad; ii++) {
    for ( int jj = 0; jj < N_POINTS_QUADRILATERAL; jj++ ) {
      
      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/
      
      iNode = ii*N_POINTS_QUADRILATERAL+jj;
      Global_Index = Conn_BoundQuad_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
      
      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/
      
      if ((nElem_Flag[iProcessor] != iNode)) {
        nElem_Flag[iProcessor] = iNode;
        nElem_Send[iProcessor+1]++;
      }
      
    }
  }
  
  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many nodes it will receive from each other processor. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif
  
  /*--- Prepare to send. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
  
  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;
    
    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }
  
  /*--- Allocate arrays for sending the global ID. ---*/
  
  unsigned long *idSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++) idSend[ii] = 0;
  
  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/
  
  unsigned long *idIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) idIndex[ii] = nElem_Send[ii];
  
  /*--- Now loop back through the local connectivities for the surface
   elements and load up the global IDs for sending to their home proc. ---*/
  
  for (int ii = 0; ii < (int)nParallel_Line; ii++) {
    for ( int jj = 0; jj < N_POINTS_LINE; jj++ ) {
      
      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/
      
      iNode = ii*N_POINTS_LINE+jj;
      Global_Index = Conn_BoundLine_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
      
      /*--- Load global ID into the buffer for sending ---*/
      
      if (nElem_Flag[iProcessor] != iNode) {
        
        nElem_Flag[iProcessor] = iNode;
        unsigned long nn = idIndex[iProcessor];
        
        /*--- Load the connectivity values. ---*/
        
        idSend[nn] = Global_Index; nn++;
        
        /*--- Increment the index by the message length ---*/
        
        idIndex[iProcessor]++;
        
      }
      
    }
  }
  
  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;
  
  for (int ii = 0; ii < (int)nParallel_BoundTria; ii++) {
    for ( int jj = 0; jj < N_POINTS_TRIANGLE; jj++ ) {
      
      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/
      
      iNode = ii*N_POINTS_TRIANGLE + jj;
      Global_Index = Conn_BoundTria_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
      
      /*--- Load global ID into the buffer for sending ---*/
      
      if (nElem_Flag[iProcessor] != iNode) {
        
        nElem_Flag[iProcessor] = iNode;
        unsigned long nn = idIndex[iProcessor];
        
        /*--- Load the connectivity values. ---*/
        
        idSend[nn] = Global_Index; nn++;
        
        /*--- Increment the index by the message length ---*/
        
        idIndex[iProcessor]++;
        
      }
      
    }
  }
  
  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;
  
  for (int ii = 0; ii < (int)nParallel_BoundQuad; ii++) {
    for ( int jj = 0; jj < N_POINTS_QUADRILATERAL; jj++ ) {
      
      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/
      
      iNode = ii*N_POINTS_QUADRILATERAL+jj;
      Global_Index = Conn_BoundQuad_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
      
      /*--- Load global ID into the buffer for sending ---*/
      
      if (nElem_Flag[iProcessor] != iNode) {
        
        nElem_Flag[iProcessor] = iNode;
        unsigned long nn = idIndex[iProcessor];
        
        /*--- Load the connectivity values. ---*/
        
        idSend[nn] = Global_Index; nn++;
        
        /*--- Increment the index by the message length ---*/
        
        idIndex[iProcessor]++;
        
      }
      
    }
  }
  
  /*--- Allocate the memory that we need for receiving the global IDs
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  unsigned long *idRecv = NULL;
  idRecv = new unsigned long[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    idRecv[ii] = 0;
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/
  
  send_req = new SU2_MPI::Request[nSends];
  recv_req = new SU2_MPI::Request[nRecvs];
  
  /*--- Launch the non-blocking recv's for the global IDs. ---*/
  
  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the global IDs. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
#endif
  
  /*--- Copy my own rank's data into the recv buffer directly. ---*/
  
  int mm = nElem_Recv[rank];
  int ll = nElem_Send[rank];
  int kk = nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  int number = nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Each proc now knows which is its local grid nodes from     ---*/
  /*---         the entire volume solution are part of the surface. We     ---*/
  /*---         now apply a mask to extract just those points on the       ---*/
  /*---         surface. We also need to perform a renumbering so that     ---*/
  /*---         the surface data (nodes and connectivity) have their       ---*/
  /*---         own global numbering. This is important for writing        ---*/
  /*---         output files in a later routine.                           ---*/
  /*--------------------------------------------------------------------------*/
  
  /*--- Create a local data structure that acts as a mask to extract the
   set of points within the local set that are on the surface. ---*/
  
  int *surfPoint = new int[nParallel_Poin];
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) surfPoint[iPoint] = -1;
  
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    surfPoint[(int)idRecv[ii]- starting_node[rank]] = (int)idRecv[ii];
  }
  
  /*--- First, add up the number of surface points I have on my rank. ---*/
  
  nSurf_Poin_Par = 0;
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
    if (surfPoint[iPoint] != -1) {
      nSurf_Poin_Par++;
    }
  }
  
  /*--- Communicate this number of local surface points to all other
   processors so that it can be used to create offsets for the new
   global numbering for the surface points. ---*/
  
  int *nPoint_Send = new int[size+1]; nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1]; nPoint_Recv[0] = 0;
  
  for (int ii=1; ii < size+1; ii++) nPoint_Send[ii]= (int)nSurf_Poin_Par;
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nPoint_Recv[1] = nPoint_Send[1];
#endif
  
  /*--- Go to cumulative storage format to compute the offsets. ---*/
  
  for (int ii = 0; ii < size; ii++) {
    nPoint_Send[ii+1] += nPoint_Send[ii];
    nPoint_Recv[ii+1] += nPoint_Recv[ii];
  }
  
  /*--- Now that we know the number of local surface points that we have,
   we can allocate the new data structure to hold these points alone. Here,
   we also copy the data for those points from our volume data structure. ---*/
  
  Parallel_Surf_Data = new su2double*[VARS_PER_POINT];
  for (int jj = 0; jj < VARS_PER_POINT; jj++) {
    Parallel_Surf_Data[jj] = new su2double[nSurf_Poin_Par];
    count = 0;
    for (int ii = 0; ii < (int)nParallel_Poin; ii++) {
      if (surfPoint[ii] !=-1) {
        Parallel_Surf_Data[jj][count] = Parallel_Data[jj][ii];
        count++;
      }
    }
  }
  
  /*--- Reduce the total number of surf points we have. This will be
   needed for writing the surface solution files later. ---*/
  
#ifndef HAVE_MPI
  nGlobal_Surf_Poin = nSurf_Poin_Par;
#else
  SU2_MPI::Allreduce(&nSurf_Poin_Par, &nGlobal_Surf_Poin, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Now that we know every proc's global offset for the number of
   surface points, we can create the new global numbering. Here, we
   create a new mapping using two arrays, which will need to be
   communicated. We use our mask again here.  ---*/
  
  unsigned long *globalP = new unsigned long[nSurf_Poin_Par];
  unsigned long *renumbP = new unsigned long[nSurf_Poin_Par];
  
  count = 0;
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
    if (surfPoint[iPoint] != -1) {
      globalP[count] = surfPoint[iPoint];
      renumbP[count] = count + nPoint_Recv[rank];
      count++;
    }
  }
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Communicate the arrays with the new global surface point   ---*/
  /*---         numbering to the procs that hold the connectivity for      ---*/
  /*---         each element. This will be done in two phases. First,      ---*/
  /*---         we send the arrays around to the other procs based on      ---*/
  /*---         the linear partitioning for the elems. This gets us        ---*/
  /*---         most of the way there, however, due to the type of         ---*/
  /*---         linear partitioning for the elements, there may exist      ---*/
  /*---         elements that have nodes outside of the linear part.       ---*/
  /*---         bounds. This is because the elems are distributed based    ---*/
  /*---         on the node with the smallest global ID.                   ---*/
  /*--------------------------------------------------------------------------*/
  
  /*--- First, we perform the linear partitioning again as it is done
   for elements, which is slightly different than for nodes (above). ---*/
  
  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that were part of the original domain (we want these to be in the
   output files). ---*/
  
  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  
  if (config->GetKind_SU2() != SU2_DEF) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
              (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
              (SendRecv < 0)) {
            Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
          }
        }
      }
    }
  }
  
  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/
  
  /*--- First communicate the number of points that each rank has found. ---*/
  
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];
  
  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif
  
  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];
  
  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }
  
  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++)
    Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif
  
  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        
        /*--- We need to keep one copy of overlapping halo cells. ---*/
        
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));
        
        /*--- We want to keep the periodic nodes that were part of the original domain.
         For SU2_DEF we want to keep all periodic nodes. ---*/
        
        if (config->GetKind_SU2() == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }
        
        notPeriodic = (isPeriodic && (SendRecv < 0));
        
        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < (unsigned long)size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (Global_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          
          kPoint = (iProcessor+1)*maxAddedPeriodic;
          
        }
        
        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
        
      }
    }
  }
  
  /*--- Now that we've done the gymnastics to find any periodic points,
   compute the total number of local and global points for the output. ---*/
  
  nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;
  
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif
  
  /*--- Compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }
  
  /*--- Get the number of remainder points after the even division. ---*/
  
  rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }
  
  /*--- Store the local number of nodes and the beginning/end index ---*/
  
  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear_Elems[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear_Elems[ii] = nPoint_Linear_Elems[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear_Elems[size] = nTotalPoint;
  
  /*--- Reset our flags and counters ---*/
  
  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;
  
  /*--- Loop through my local surface nodes, find which proc the global
   value lives on, then communicate the global ID and remumbered value. ---*/
  
  for (int ii = 0; ii < (int)nSurf_Poin_Par; ii++) {
    
    Global_Index = globalP[ii];
    
    /*--- Search for the processor that owns this point ---*/
    
    iProcessor = Global_Index/npoint_procs[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Linear_Elems[iProcessor])
      while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;
    
    /*--- If we have not visited this element yet, increment our
     number of elements that must be sent to a particular proc. ---*/
    
    if ((nElem_Flag[iProcessor] != ii)) {
      nElem_Flag[iProcessor] = ii;
      nElem_Send[iProcessor+1]++;
    }
    
  }
  
  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif
  
  /*--- Prepare to send. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
  
  nSends = 0; nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;
    
    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }
  
  /*--- Allocate memory to hold the globals that we are
   sending. ---*/
  
  unsigned long *globalSend = NULL;
  globalSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    globalSend[ii] = 0;
  
  /*--- Allocate memory to hold the renumbering that we are
   sending. ---*/
  
  unsigned long *renumbSend = NULL;
  renumbSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    renumbSend[ii] = 0;
  
  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/
  
  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = nElem_Send[ii];
  
  /*--- Loop back through and load up the buffers for the global IDs
   and their new renumbering values. ---*/
  
  for (int ii = 0; ii < (int)nSurf_Poin_Par; ii++) {
    
    Global_Index = globalP[ii];
    
    /*--- Search for the processor that owns this point ---*/
    
    iProcessor = Global_Index/npoint_procs[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Linear_Elems[iProcessor])
      while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;
    
    
    if (nElem_Flag[iProcessor] != ii) {
      
      nElem_Flag[iProcessor] = ii;
      unsigned long nn = index[iProcessor];
      
      globalSend[nn] = Global_Index;
      renumbSend[nn] = renumbP[ii];
      
      /*--- Increment the index by the message length ---*/
      
      index[iProcessor]++;
      
    }
  }
  
  /*--- Free memory after loading up the send buffer. ---*/
  
  delete [] index;
  
  /*--- Allocate the memory that we need for receiving the
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  unsigned long *globalRecv = NULL;
  globalRecv = new unsigned long[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    globalRecv[ii] = 0;
  
  unsigned long *renumbRecv = NULL;
  renumbRecv = new unsigned long[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    renumbRecv[ii] = 0;
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/
  
  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];
  
  /*--- Launch the non-blocking recv's for the global ID. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(globalRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the global ID. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(globalSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking recv's for the renumbered ID. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(renumbRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the renumbered ID. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(renumbSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
  
#endif
  
  /*--- Load our own procs data into the buffers directly. ---*/
  
  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) globalRecv[mm] = globalSend[nn];
  
  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) renumbRecv[mm] = renumbSend[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*-- Now update my local connectivitiy for the surface with the new
   numbering. Create a new mapping for global -> renumber for nodes. Note
   the adding of 1 back in here for the eventual viz. purposes. ---*/
  
  map<unsigned long,unsigned long> Global2Renumber;
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    Global2Renumber[globalRecv[ii]] = renumbRecv[ii] + 1;
  }
  
  
  /*--- The final step is one last pass over all elements to check
   for points outside of the linear partitions of the elements. Again,
   note that elems were distributed based on their smallest global ID,
   so some nodes of the elem may have global IDs lying outside of the
   linear partitioning. We need to recover the mapping for these
   outliers. We loop over all local surface elements to find these. ---*/
  
  vector<unsigned long>::iterator it;
  vector<unsigned long> outliers;
  
  for (int ii = 0; ii < (int)nParallel_Line; ii++) {
    for ( int jj = 0; jj < N_POINTS_LINE; jj++ ) {
      
      iNode = ii*N_POINTS_LINE+jj;
      Global_Index = Conn_BoundLine_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Elems[iProcessor])
        while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;
      
      /*--- Store the global ID if it is outside our own linear partition. ---*/
      
      if ((iProcessor != (unsigned long)rank)) {
        outliers.push_back(Global_Index);
      }
      
    }
  }
  
  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;
  
  for (int ii = 0; ii < (int)nParallel_BoundTria; ii++) {
    for ( int jj = 0; jj < N_POINTS_TRIANGLE; jj++ ) {
      
      iNode = ii*N_POINTS_TRIANGLE + jj;
      Global_Index = Conn_BoundTria_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Elems[iProcessor])
        while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;
      
      /*--- Store the global ID if it is outside our own linear partition. ---*/
      
      if ((iProcessor != (unsigned long)rank)) {
        outliers.push_back(Global_Index);
      }
      
    }
  }
  
  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;
  
  for (int ii = 0; ii < (int)nParallel_BoundQuad; ii++) {
    for ( int jj = 0; jj < N_POINTS_QUADRILATERAL; jj++ ) {
      
      iNode = ii*N_POINTS_QUADRILATERAL+jj;
      Global_Index = Conn_BoundQuad_Par[iNode]-1;
      
      /*--- Search for the processor that owns this point ---*/
      
      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Elems[iProcessor])
        while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;
      
      /*--- Store the global ID if it is outside our own linear partition. ---*/
      
      if ((iProcessor != (unsigned long)rank)) {
        outliers.push_back(Global_Index);
      }
      
    }
  }
  
  /*--- Create a unique list of global IDs that fall outside of our procs
   linear partition. ---*/
  
  sort(outliers.begin(), outliers.end());
  it = unique(outliers.begin(), outliers.end());
  outliers.resize(it - outliers.begin());
  
  /*--- Now loop over the outliers and communicate to those procs that
   hold the new numbering for our outlier points. We need to ask for the
   new numbering from these procs. ---*/
  
  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;
  
  for (int ii = 0; ii < (int)outliers.size(); ii++) {
    
    Global_Index = outliers[ii];
    
    /*--- Search for the processor that owns this point ---*/
    
    iProcessor = Global_Index/npoint_procs[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
      while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
    
    /*--- If we have not visited this element yet, increment our
     number of elements that must be sent to a particular proc. ---*/
    
    if ((nElem_Flag[iProcessor] != ii)) {
      nElem_Flag[iProcessor] = ii;
      nElem_Send[iProcessor+1]++;
    }
    
  }
  
  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/
  
#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif
  
  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/
  
  nSends = 0; nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;
  
  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;
    
    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }
  
  delete [] idSend;
  idSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    idSend[ii] = 0;
  
  /*--- Reset our index variable for reuse. ---*/
  
  for (int ii=0; ii < size; ii++) idIndex[ii] = nElem_Send[ii];
  
  /*--- Loop over the outliers again and load up the global IDs. ---*/
  
  for (int ii = 0; ii < (int)outliers.size(); ii++) {
    
    Global_Index = outliers[ii];
    
    /*--- Search for the processor that owns this point ---*/
    
    iProcessor = Global_Index/npoint_procs[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
      while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;
    
    /*--- If we have not visited this element yet, increment our
     number of elements that must be sent to a particular proc. ---*/
    
    if ((nElem_Flag[iProcessor] != ii)) {
      
      nElem_Flag[iProcessor] = ii;
      unsigned long nn = idIndex[iProcessor];
      
      /*--- Load the global ID values. ---*/
      
      idSend[nn] = Global_Index; nn++;
      
      /*--- Increment the index by the message length ---*/
      
      idIndex[iProcessor]++;
      
    }
  }
  
  /*--- Allocate the memory that we need for receiving the
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/
  
  delete [] idRecv;
  idRecv = new unsigned long[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    idRecv[ii] = 0;
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/
  
  send_req = new SU2_MPI::Request[nSends];
  recv_req = new SU2_MPI::Request[nRecvs];
  
  /*--- Launch the non-blocking recv's for the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking sends of the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
#endif
  
  /*--- Copy my own rank's data into the recv buffer directly. ---*/
  
  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  number = nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*--- The procs holding the outlier grid nodes now have the global IDs
   that they need to have their renumbering shared. ---*/
  
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    for (iPoint = 0; iPoint < nSurf_Poin_Par; iPoint++) {
      if (idRecv[ii] == globalP[iPoint]) {
        idRecv[ii] = renumbP[iPoint];
      }
    }
  }
  
  /*--- Now simply reverse the last communication to give the renumbered IDs
   back to the owner of the outlier points. Note everything is flipped. ---*/
  
#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/
  
  send_req = new SU2_MPI::Request[nRecvs];
  recv_req = new SU2_MPI::Request[nSends];
  
  /*--- Launch the non-blocking sends of the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }
  
  /*--- Launch the non-blocking recv's for the connectivity. ---*/
  
  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
#endif
  
  /*--- Copy my own rank's data into the recv buffer directly. ---*/
  
  mm = nElem_Send[rank];
  ll = nElem_Recv[rank];
  kk = nElem_Recv[rank+1];
  
  for (int nn=ll; nn<kk; nn++, mm++) idSend[mm] = idRecv[nn];
  
  /*--- Wait for the non-blocking sends and recvs to complete. ---*/
  
#ifdef HAVE_MPI
  number = nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);
  
  number = nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);
  
  delete [] send_req;
  delete [] recv_req;
#endif
  
  /*--- Add the renumbering for the outliers to the map from before carrying
   the global -> renumber transformation. Note that by construction, 
   nElem_Send[ii] == outliers.size(). We also add in the 1 for viz. here. ---*/
  
  for (int ii = 0; ii < nElem_Send[size]; ii++) {
    Global2Renumber[outliers[ii]] = idSend[ii] + 1;
  }
  
  /*--- We can now overwrite the local connectivity for our surface elems
   using our completed map with the new global renumbering. Whew!! Note
   the -1 when accessing the conn from the map. ---*/
  
  for (iElem = 0; iElem < nParallel_Line; iElem++) {
    iNode = (int)iElem*N_POINTS_LINE;
    Conn_BoundLine_Par[iNode+0] = (int)Global2Renumber[Conn_BoundLine_Par[iNode+0]-1];
    Conn_BoundLine_Par[iNode+1] = (int)Global2Renumber[Conn_BoundLine_Par[iNode+1]-1];
  }
  
  for (iElem = 0; iElem < nParallel_BoundTria; iElem++) {
    iNode = (int)iElem*N_POINTS_TRIANGLE;
    Conn_BoundTria_Par[iNode+0] = (int)Global2Renumber[Conn_BoundTria_Par[iNode+0]-1];
    Conn_BoundTria_Par[iNode+1] = (int)Global2Renumber[Conn_BoundTria_Par[iNode+1]-1];
    Conn_BoundTria_Par[iNode+2] = (int)Global2Renumber[Conn_BoundTria_Par[iNode+2]-1];
  }
  
  for (iElem = 0; iElem < nParallel_BoundQuad; iElem++) {
    iNode = (int)iElem*N_POINTS_QUADRILATERAL;
    Conn_BoundQuad_Par[iNode+0] = (int)Global2Renumber[Conn_BoundQuad_Par[iNode+0]-1];
    Conn_BoundQuad_Par[iNode+1] = (int)Global2Renumber[Conn_BoundQuad_Par[iNode+1]-1];
    Conn_BoundQuad_Par[iNode+2] = (int)Global2Renumber[Conn_BoundQuad_Par[iNode+2]-1];
    Conn_BoundQuad_Par[iNode+3] = (int)Global2Renumber[Conn_BoundQuad_Par[iNode+3]-1];
  }
  
  /*--- Free temporary memory ---*/
  
  delete [] idIndex;
  delete [] surfPoint;
  delete [] globalP;
  delete [] renumbP;
  
  delete [] idSend;
  delete [] idRecv;
  delete [] globalSend;
  delete [] globalRecv;
  delete [] renumbSend;
  delete [] renumbRecv;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;
  delete [] Local_Halo;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear_Elems;
  delete [] nPoint_Linear_Nodes;
  delete [] nPoint_Send;
  delete [] nPoint_Recv;
  
}

void COutput::WriteRestart_Parallel_ASCII(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_iInst) {
  
  /*--- Local variables ---*/
  
  unsigned short nZone = geometry->GetnZone(), nInst = config->GetnTimeInstances();
  unsigned short iVar;
  unsigned long iPoint, iExtIter = config->GetExtIter();
  bool fem       = (config->GetKind_Solver() == FEM_ELASTICITY);
  bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);
  bool adjoint   = (config->GetContinuous_Adjoint() ||
                    config->GetDiscrete_Adjoint());
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  ofstream restart_file;
  string filename;
  
  int iProcessor;

  /*--- Retrieve filename from config ---*/
  
  if ((config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint())) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else if (fem) {
    filename = config->GetRestart_FEMFileName();
  } else if (disc_adj_fem){
    filename = config->GetRestart_AdjFEMFileName();
  } else {
    filename = config->GetRestart_FlowFileName();
  }
  
  /*--- Append the zone number if multizone problems ---*/
  if (nZone > 1)
    filename= config->GetMultizone_FileName(filename, val_iZone);
  
  /*--- Append the zone number if multiple instance problems ---*/
  if (nInst > 1)
    filename= config->GetMultiInstance_FileName(filename, val_iInst);

  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(val_iInst));
  } else if (config->GetWrt_Unsteady()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  } else if ((fem || disc_adj_fem) && (config->GetWrt_Dynamic())) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  }
  
  /*--- Only the master node writes the header. ---*/
  
  if (rank == MASTER_NODE) {
    restart_file.open(filename.c_str(), ios::out);
    restart_file.precision(15);
    restart_file << "\"PointID\"";
    for (iVar = 0; iVar < Variable_Names.size()-1; iVar++)
      restart_file << "\t\"" << Variable_Names[iVar] << "\"";
    restart_file << "\t\"" << Variable_Names[Variable_Names.size()-1] << "\"" << endl;
    restart_file.close();
  }
  
#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
  
  /*--- All processors open the file. ---*/
  
  restart_file.open(filename.c_str(), ios::out | ios::app);
  restart_file.precision(15);
  
  /*--- Write the restart file in parallel, processor by processor. ---*/
  
  unsigned long myPoint = 0, offset = 0, Global_Index;
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
        
        /*--- Global Index of the current point. (note outer loop over procs) ---*/
        
        Global_Index = iPoint + offset;
        
        /*--- Only write original domain points, i.e., exclude any periodic
         or halo nodes, even if they are output in the viz. files. ---*/
        
        if (Global_Index < geometry->GetGlobal_nPointDomain()) {
          
          /*--- Write global index. (note outer loop over procs) ---*/
          
          restart_file << Global_Index << "\t";
          myPoint++;
          
          /*--- Loop over the variables and write the values to file ---*/
          
          for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
            restart_file << scientific << Parallel_Data[iVar][iPoint] << "\t";
          }
          restart_file << "\n";
        }
      }
    }
    /*--- Flush the file and wait for all processors to arrive. ---*/
    restart_file.flush();
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&myPoint, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif
    
  }

  /*--- Write the metadata (master rank alone) ----*/

  if (rank == MASTER_NODE) {
    if (dual_time)
      restart_file <<"EXT_ITER= " << config->GetExtIter() + 1 << endl;
    else
      restart_file <<"EXT_ITER= " << config->GetExtIter() + config->GetExtIter_OffSet() + 1 << endl;
    restart_file <<"AOA= " << config->GetAoA() - config->GetAoA_Offset() << endl;
    restart_file <<"SIDESLIP_ANGLE= " << config->GetAoS() - config->GetAoS_Offset() << endl;
    restart_file <<"INITIAL_BCTHRUST= " << config->GetInitial_BCThrust() << endl;
    restart_file <<"DCD_DCL_VALUE= " << config->GetdCD_dCL() << endl;
    restart_file <<"DCMX_DCL_VALUE= " << config->GetdCMx_dCL() << endl;
    restart_file <<"DCMY_DCL_VALUE= " << config->GetdCMy_dCL() << endl;
    restart_file <<"DCMZ_DCL_VALUE= " << config->GetdCMz_dCL() << endl;
    if (adjoint) restart_file << "SENS_AOA=" << solver[ADJFLOW_SOL]->GetTotal_Sens_AoA() * PI_NUMBER / 180.0 << endl;
  }

  /*--- All processors close the file. ---*/

  restart_file.close();
  
}

void COutput::WriteRestart_Parallel_Binary(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_iInst) {

  /*--- Local variables ---*/

  unsigned short iVar, nZone = geometry->GetnZone(), nInst = config->GetnTimeInstances();
  unsigned long iPoint, iExtIter = config->GetExtIter();
  bool fem       = (config->GetKind_Solver() == FEM_ELASTICITY);
  bool adjoint   = (config->GetContinuous_Adjoint() ||
                    config->GetDiscrete_Adjoint());
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool wrt_perf  = config->GetWrt_Performance();
  ofstream restart_file;
  string filename;
  char str_buf[CGNS_STRING_SIZE], fname[100];
  su2double file_size = 0.0, StartTime, StopTime, UsedTime, Bandwidth;

  /*--- Retrieve filename from config ---*/

  if ((config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint())) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else if (fem) {
    filename = config->GetRestart_FEMFileName();
  } else {
    filename = config->GetRestart_FlowFileName();
  }

  /*--- Append the zone number if multizone problems ---*/
  if (nZone > 1)
    filename= config->GetMultizone_FileName(filename, val_iZone);

  /*--- Append the zone number if multiple instance problems ---*/
  if (nInst > 1)
    filename= config->GetMultiInstance_FileName(filename, val_iInst);

  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(val_iInst));
  } else if (config->GetWrt_Unsteady()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  } else if ((fem) && (config->GetWrt_Dynamic())) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  }

  strcpy(fname, filename.c_str());

  /*--- Prepare the first ints containing the counts. The first is a
   magic number that we can use to check for binary files (it is the hex
   representation for "SU2"). The second two values are number of variables
   and number of points (DoFs). The last two values are for metadata: 
   one int for ExtIter and 8 su2doubles. ---*/

  int var_buf_size = 5;
  int var_buf[5] = {535532, GlobalField_Counter, (int)nGlobalPoint_Sort, 1, 8};

  /*--- Prepare the 1D data buffer on this rank. ---*/

  passivedouble *buf = new passivedouble[nParallel_Poin*GlobalField_Counter];

  /*--- For now, create a temp 1D buffer to load up the data for writing.
   This will be replaced with a derived data type most likely. ---*/

  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++)
    for (iVar = 0; iVar < GlobalField_Counter; iVar++)
      buf[iPoint*GlobalField_Counter+iVar] = SU2_TYPE::GetValue(Parallel_Data[iVar][iPoint]);

  /*--- Prepare metadata. ---*/

  int Restart_ExtIter;
  if (dual_time)
    Restart_ExtIter= (int)config->GetExtIter() + 1;
  else
    Restart_ExtIter = (int)config->GetExtIter() + (int)config->GetExtIter_OffSet() + 1;

  passivedouble Restart_Metadata[8] = {
    SU2_TYPE::GetValue(config->GetAoA() - config->GetAoA_Offset()),
    SU2_TYPE::GetValue(config->GetAoS() - config->GetAoS_Offset()),
    SU2_TYPE::GetValue(config->GetInitial_BCThrust()),
    SU2_TYPE::GetValue(config->GetdCD_dCL()),
    SU2_TYPE::GetValue(config->GetdCMx_dCL()),
    SU2_TYPE::GetValue(config->GetdCMy_dCL()),
    SU2_TYPE::GetValue(config->GetdCMz_dCL()),
    0.0
  };

  if (adjoint) Restart_Metadata[4] = SU2_TYPE::GetValue(solver[ADJFLOW_SOL]->GetTotal_Sens_AoA() * PI_NUMBER / 180.0);

  /*--- Set a timer for the binary file writing. ---*/
  
#ifndef HAVE_MPI
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StartTime = MPI_Wtime();
#endif
  
#ifndef HAVE_MPI

  FILE* fhw;
  fhw = fopen(fname, "wb");

  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
  }

  /*--- First, write the number of variables and points. ---*/

  fwrite(var_buf, var_buf_size, sizeof(int), fhw);
  file_size += (su2double)var_buf_size*sizeof(int);
  
  /*--- Write the variable names to the file. Note that we are adopting a
   fixed length of 33 for the string length to match with CGNS. This is 
   needed for when we read the strings later. ---*/

  for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
    strncpy(str_buf, Variable_Names[iVar].c_str(), CGNS_STRING_SIZE);
    fwrite(str_buf, CGNS_STRING_SIZE, sizeof(char), fhw);
    file_size += (su2double)CGNS_STRING_SIZE*sizeof(char);
  }

  /*--- Call to write the entire restart file data in binary in one shot. ---*/

  fwrite(buf, GlobalField_Counter*nParallel_Poin, sizeof(passivedouble), fhw);
  file_size += (su2double)GlobalField_Counter*nParallel_Poin*sizeof(passivedouble);

  /*--- Write the external iteration. ---*/

  fwrite(&Restart_ExtIter, 1, sizeof(int), fhw);
  file_size += (su2double)sizeof(int);

  /*--- Write the metadata. ---*/

  fwrite(Restart_Metadata, 8, sizeof(passivedouble), fhw);
  file_size += (su2double)8*sizeof(passivedouble);

  /*--- Close the file. ---*/

  fclose(fhw);

#else

  /*--- Parallel binary output using MPI I/O. ---*/

  MPI_File fhw;
  SU2_MPI::Status status;
  MPI_Datatype etype, filetype;
  MPI_Offset disp;
  int ierr;

  /*--- We're writing only su2doubles in the data portion of the file. ---*/

  etype = MPI_DOUBLE;

  /*--- Define a derived datatype for this ranks contiguous chunk of data
   that will be placed in the restart (1D array size = num points * num vars). ---*/

  MPI_Type_contiguous(GlobalField_Counter*nParallel_Poin, MPI_DOUBLE, &filetype);
  MPI_Type_commit(&filetype);

  /*--- All ranks open the file using MPI. Here, we try to open the file with
   exclusive so that an error is generated if the file exists. We always want
   to write a fresh restart file, so we delete any existing files and create
   a new one. ---*/

  ierr = MPI_File_open(MPI_COMM_WORLD, fname,
                       MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                       MPI_INFO_NULL, &fhw);
  if (ierr != MPI_SUCCESS)  {
    MPI_File_close(&fhw);
    if (rank == 0)
      MPI_File_delete(fname, MPI_INFO_NULL);
    ierr = MPI_File_open(MPI_COMM_WORLD, fname,
                         MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                         MPI_INFO_NULL, &fhw);
  }

  /*--- Error check opening the file. ---*/

  if (ierr) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
  }

  /*--- First, write the number of variables and points (i.e., cols and rows),
   which we will need in order to read the file later. Also, write the 
   variable string names here. Only the master rank writes the header. ---*/

  if (rank == MASTER_NODE) {
    MPI_File_write(fhw, var_buf, var_buf_size, MPI_INT, MPI_STATUS_IGNORE);
    file_size += (su2double)var_buf_size*sizeof(int);

    /*--- Write the variable names to the file. Note that we are adopting a
     fixed length of 33 for the string length to match with CGNS. This is
     needed for when we read the strings later. ---*/

    for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
      disp = var_buf_size*sizeof(int) + iVar*CGNS_STRING_SIZE*sizeof(char);
      strcpy(str_buf, Variable_Names[iVar].c_str());
      MPI_File_write_at(fhw, disp, str_buf, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);
      file_size += (su2double)CGNS_STRING_SIZE*sizeof(char);
    }
  }

  /*--- Compute the offset for this rank's linear partition of the data in bytes.
   After the calculations above, we have the partition sizes store in nPoint_Linear
   in cumulative storage format. ---*/

  disp = (var_buf_size*sizeof(int) + GlobalField_Counter*CGNS_STRING_SIZE*sizeof(char) +
          GlobalField_Counter*nPoint_Cum[rank]*sizeof(passivedouble));

  /*--- Set the view for the MPI file write, i.e., describe the location in
   the file that this rank "sees" for writing its piece of the restart file. ---*/

  MPI_File_set_view(fhw, disp, etype, filetype, (char*)"native", MPI_INFO_NULL);

  /*--- Collective call for all ranks to write to their view simultaneously. ---*/

  MPI_File_write_all(fhw, buf, GlobalField_Counter*nParallel_Poin, MPI_DOUBLE, &status);
  file_size += (su2double)GlobalField_Counter*nParallel_Poin*sizeof(passivedouble);

  /*--- Free the derived datatype. ---*/

  MPI_Type_free(&filetype);

  /*--- Reset the file view before writing the metadata. ---*/

  MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE, (char*)"native", MPI_INFO_NULL);

  /*--- Finally, the master rank writes the metadata. ---*/

  if (rank == MASTER_NODE) {

    /*--- External iteration. ---*/

    disp = (var_buf_size*sizeof(int) + GlobalField_Counter*CGNS_STRING_SIZE*sizeof(char) +
            GlobalField_Counter*nGlobalPoint_Sort*sizeof(passivedouble));
    MPI_File_write_at(fhw, disp, &Restart_ExtIter, 1, MPI_INT, MPI_STATUS_IGNORE);
    file_size += (su2double)sizeof(int);

    /*--- Additional doubles for AoA, AoS, etc. ---*/

    disp = (var_buf_size*sizeof(int) + GlobalField_Counter*CGNS_STRING_SIZE*sizeof(char) +
            GlobalField_Counter*nGlobalPoint_Sort*sizeof(passivedouble) + 1*sizeof(int));
    MPI_File_write_at(fhw, disp, Restart_Metadata, 8, MPI_DOUBLE, MPI_STATUS_IGNORE);
    file_size += (su2double)8*sizeof(passivedouble);

  }

  /*--- All ranks close the file after writing. ---*/

  MPI_File_close(&fhw);

#endif

  /*--- Compute and store the write time. ---*/
  
#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  UsedTime = StopTime-StartTime;
  
  /*--- Communicate the total file size for the restart ---*/
  
#ifdef HAVE_MPI
  su2double my_file_size = file_size;
  SU2_MPI::Allreduce(&my_file_size, &file_size, 1,
                     MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Compute and store the bandwidth ---*/
  
  Bandwidth = file_size/(1.0e6)/UsedTime;
  config->SetRestart_Bandwidth_Agg(config->GetRestart_Bandwidth_Agg()+Bandwidth);
  
  if ((rank == MASTER_NODE) && (wrt_perf)) {
    cout << "Wrote " << file_size/1.0e6 << " MB to disk in ";
    cout << UsedTime << " s. (" << Bandwidth << " MB/s)." << endl;
  }
  
  /*--- Free temporary data buffer for writing the binary file. ---*/

  delete [] buf;

}

void COutput::WriteCSV_Slice(CConfig *config, CGeometry *geometry,
                             CSolver *FlowSolver, unsigned long iExtIter,
                             unsigned short val_iZone, unsigned short val_direction) {

  /*--- This routine is for exporting slices of field data for 2D cartesian
   grids. It assumes that the grid points lie on lines of constant x-
   or y-coordinates. It is a simple way to export slices or profiles on
   these meshes for use in verification and validation work. It will
   eventually be replaced by general routines for probing/slicing.  ---*/

  int DIRECTION = (int)val_direction;

  su2double coordMin, coordMax;
  coordMin = config->GetStations_Bounds(0);
  coordMax = config->GetStations_Bounds(1);

  unsigned short iVar;
  unsigned long iPoint, iVertex, Global_Index;
  char cstr[200];

  int rank = MASTER_NODE, iProcessor, nProcessor = SINGLE_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

  bool isPeriodic;

  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

      /*--- Checking for less than or equal to the rank, because there may
       be some periodic halo nodes that send info to the same rank. ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                      (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        if (isPeriodic) Local_Halo[iPoint] = false;
      }
    }
  }

  /*--- Sum total number of nodes that belong to the domain ---*/

  unsigned long nTotalPoint;
  unsigned long nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif

  unsigned long *npoint_procs  = new unsigned long[nProcessor];
  unsigned long *nPoint_Linear = new unsigned long[nProcessor+1];

  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < nProcessor; ii++) {
    npoint_procs[ii] = nTotalPoint/nProcessor;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }

  /*--- Get the number of remainder points after the even division. ---*/

  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }

  /*--- Store the point offsets for each rank. ---*/

  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < nProcessor; ii++) {
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[nProcessor] = nTotalPoint;

  unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = NULL;
  unsigned long nLocalVertex_Surface = 0;
  unsigned long MaxLocalVertex_Surface = 0;

  /*--- Find the max number of vertices we will send from among all
   partitions and set up buffers. The master node will handle the
   writing of the CSV file after gathering all of the data. ---*/

  nLocalVertex_Surface = 0;
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {

    /*--- Global Index of the current point. (note outer loop over procs) ---*/

    Global_Index = iPoint + nPoint_Linear[rank];

    /*--- Only write original domain points, i.e., exclude any periodic
     or halo nodes, even if they are output in the viz. files. ---*/

    if (Global_Index < geometry->GetGlobal_nPointDomain()) {
      if ((Parallel_Data[DIRECTION][iPoint] > coordMin) &&
          (Parallel_Data[DIRECTION][iPoint] < coordMax)) {
        nLocalVertex_Surface++;
      }
    }
  }

  /*--- Communicate the number of local vertices on each partition
   to the master node ---*/

  Buffer_Send_nVertex[0] = nLocalVertex_Surface;
  if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long [nProcessor];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalVertex_Surface = nLocalVertex_Surface;
  Buffer_Recv_nVertex[0] = Buffer_Send_nVertex[0];
#endif

  /*--- Send and Recv buffers ---*/

  su2double *Buffer_Send_Data = new su2double [MaxLocalVertex_Surface*GlobalField_Counter];
  su2double *Buffer_Recv_Data = NULL;

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long [MaxLocalVertex_Surface];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Prepare the receive buffers on the master node only. ---*/

  if (rank == MASTER_NODE) {
    Buffer_Recv_Data = new su2double [nProcessor*MaxLocalVertex_Surface*GlobalField_Counter];
    Buffer_Recv_GlobalIndex  = new unsigned long [nProcessor*MaxLocalVertex_Surface];
  }

  /*--- Loop over all vertices in this partition and load the
   data of the specified type into the buffer to be sent to
   the master node. ---*/

  nLocalVertex_Surface = 0;
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {

    /*--- Global Index of the current point. (note outer loop over procs) ---*/

    Global_Index = iPoint + nPoint_Linear[rank];

    /*--- Only write original domain points, i.e., exclude any periodic
     or halo nodes, even if they are output in the viz. files. ---*/

    if (Global_Index < geometry->GetGlobal_nPointDomain()) {
      if ((Parallel_Data[DIRECTION][iPoint] > coordMin) &&
          (Parallel_Data[DIRECTION][iPoint] < coordMax)) {
        Buffer_Send_GlobalIndex[nLocalVertex_Surface] = Global_Index;
        for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
          Buffer_Send_Data[nLocalVertex_Surface*GlobalField_Counter+iVar] = Parallel_Data[iVar][iPoint];
        }
        nLocalVertex_Surface++;
      }
    }
  }

  /*--- Send the information to the master node ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Data, MaxLocalVertex_Surface*GlobalField_Counter, MPI_DOUBLE, Buffer_Recv_Data, MaxLocalVertex_Surface*GlobalField_Counter, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (iVertex = 0; iVertex < Buffer_Recv_nVertex[0]; iVertex++) {
    Buffer_Recv_GlobalIndex[iVertex] = Buffer_Send_GlobalIndex[iVertex];
    for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
      Buffer_Recv_Data[iVertex*GlobalField_Counter+iVar] = Buffer_Send_Data[iVertex*GlobalField_Counter+iVar];
    }
  }
#endif

  /*--- The master node unpacks the data and writes the surface CSV file ---*/

  if (rank == MASTER_NODE) {

    /*--- Write file name with extension if unsteady ---*/
    char buffer[50];
    string filename = "slice";
    if (DIRECTION == 0) {
      SPRINTF (buffer, "_vert.csv");
    } else if (DIRECTION == 1) {
      SPRINTF (buffer, "_hori.csv");
    }
    ofstream SurfFlow_file;

    /*--- Write file name with extension if unsteady ---*/
    strcpy (cstr, filename.c_str());
    strcat (cstr, buffer);
    SurfFlow_file.precision(15);
    SurfFlow_file.open(cstr, ios::out);

    /*--- Global index is first, then the rest of the data. We already have the names. ---*/

    SurfFlow_file << "\"Global_Index\",";
    for (iVar = 0; iVar < Variable_Names.size()-1; iVar++) {
      SurfFlow_file << "\"" << Variable_Names[iVar] << "\",";
    }
    SurfFlow_file << "\"" << Variable_Names[Variable_Names.size()-1] << "\"" << endl;

    /*--- Loop through all of the collected data and write each node's values ---*/

    unsigned long Total_Index;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {

        /*--- Current index position and global index ---*/

        Total_Index  = iProcessor*MaxLocalVertex_Surface+iVertex;
        Global_Index = Buffer_Recv_GlobalIndex[Total_Index];

        /*--- Write the the data ---*/

        SurfFlow_file << scientific << Global_Index;
        Total_Index  = iProcessor*MaxLocalVertex_Surface*GlobalField_Counter+iVertex*GlobalField_Counter;
        for (iVar = 0; iVar < GlobalField_Counter; iVar++) {
          SurfFlow_file << scientific << ", " << Buffer_Recv_Data[Total_Index+iVar];
        }
        SurfFlow_file << endl;

      }
    }

    /*--- Close the CSV file ---*/

    SurfFlow_file.close();

    /*--- Release the recv buffers on the master node ---*/

    delete [] Buffer_Recv_Data;
    delete [] Buffer_Recv_GlobalIndex;
    delete [] Buffer_Recv_nVertex;

  }

  /*--- Release the memory for the remaining buffers and exit ---*/

  delete [] Buffer_Send_Data;
  delete [] Buffer_Send_GlobalIndex;
  
}

void COutput::DeallocateConnectivity_Parallel(CConfig *config, CGeometry *geometry, bool surf_sol) {
  
  /*--- Deallocate memory for connectivity data on each processor. ---*/
  
  if (surf_sol) {
    if (nParallel_Line > 0      && Conn_BoundLine_Par      != NULL)
      delete [] Conn_BoundLine_Par;
    if (nParallel_BoundTria > 0 && Conn_BoundTria_Par != NULL)
      delete [] Conn_BoundTria_Par;
    if (nParallel_BoundQuad > 0 && Conn_BoundQuad_Par != NULL)
      delete [] Conn_BoundQuad_Par;
  }
  else {
    if (nParallel_Tria > 0 && Conn_Tria_Par != NULL) delete [] Conn_Tria_Par;
    if (nParallel_Quad > 0 && Conn_Quad_Par != NULL) delete [] Conn_Quad_Par;
    if (nParallel_Tetr > 0 && Conn_Tetr_Par != NULL) delete [] Conn_Tetr_Par;
    if (nParallel_Hexa > 0 && Conn_Hexa_Par != NULL) delete [] Conn_Hexa_Par;
    if (nParallel_Pris > 0 && Conn_Pris_Par != NULL) delete [] Conn_Pris_Par;
    if (nParallel_Pyra > 0 && Conn_Pyra_Par != NULL) delete [] Conn_Pyra_Par;
  }
  
}

void COutput::DeallocateData_Parallel(CConfig *config, CGeometry *geometry) {
  
  /*--- Deallocate memory for solution data ---*/
  
  for (unsigned short iVar = 0; iVar < GlobalField_Counter; iVar++) {
    if (Parallel_Data[iVar] != NULL) delete [] Parallel_Data[iVar];
  }
  if (Parallel_Data != NULL) delete [] Parallel_Data;

  /*--- Deallocate the structures holding the linear partitioning ---*/

  if (Local_Halo_Sort != NULL) delete [] Local_Halo_Sort;

  if (beg_node != NULL) delete [] beg_node;
  if (end_node != NULL) delete [] end_node;

  if (nPoint_Lin != NULL) delete [] nPoint_Lin;
  if (nPoint_Cum != NULL) delete [] nPoint_Cum;

}

void COutput::DeallocateSurfaceData_Parallel(CConfig *config, CGeometry *geometry) {
  
  /*--- Deallocate memory for surface solution data ---*/

  for (unsigned short iVar = 0; iVar < GlobalField_Counter; iVar++) {
    if (Parallel_Surf_Data[iVar] != NULL) delete [] Parallel_Surf_Data[iVar];
  }
  if (Parallel_Surf_Data != NULL) delete [] Parallel_Surf_Data;
  
}

void COutput::MergeInletCoordinates(CConfig *config, CGeometry *geometry) {

  /*--- Local variables needed on all processors ---*/

  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, jPoint, kPoint;

  int iProcessor, nProcessor = size;

  unsigned long iVertex, iMarker;
  unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;

  unsigned long index, iChar;

  char str_buf[MAX_STRING_SIZE];
  vector<string> Marker_Tags;
  vector<string>::iterator it;

  unsigned long *nRowCum_Counter = NULL;

  if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];

  /*--- Search all boundaries on the present rank to count the number
   of nodes found on inlet markers. ---*/

  nLocalPoint = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Only communicate owned nodes to avoid duplicates. ---*/

        if (geometry->node[iPoint]->GetDomain())
          nLocalPoint++;

      }
    }
  }
  Buffer_Send_nPoin[0] = nLocalPoint;

  /*--- Communicate the total number of nodes on this domain. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(&Buffer_Send_nPoin, 1, MPI_UNSIGNED_LONG,
                  Buffer_Recv_nPoin, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
#else
  Buffer_Recv_nPoin[0] = Buffer_Send_nPoin[0];
  MaxLocalPoint        = nLocalPoint;
#endif

  /*--- Send and Recv buffers. ---*/

  su2double *Buffer_Send_X = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_X = NULL;

  su2double *Buffer_Send_Y = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Y = NULL;

  su2double *Buffer_Send_Z = NULL, *Buffer_Recv_Z = NULL;
  if (nDim == 3) Buffer_Send_Z = new su2double[MaxLocalPoint];

  char *Buffer_Send_Str = new char[MaxLocalPoint*MAX_STRING_SIZE];
  char *Buffer_Recv_Str = NULL;

  /*--- Prepare the receive buffers in the master node only. ---*/

  if (rank == MASTER_NODE) {

    Buffer_Recv_X = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Y = new su2double[nProcessor*MaxLocalPoint];
    if (nDim == 3) Buffer_Recv_Z = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Str = new char[nProcessor*MaxLocalPoint*MAX_STRING_SIZE];

    /*--- Sum total number of nodes to be written and allocate arrays ---*/

    unsigned long nGlobal_InletPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_InletPoint += Buffer_Recv_nPoin[iProcessor];
    }
    InletCoords = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      InletCoords[iDim] = new su2double[nGlobal_InletPoint];
    }
  }

  /*--- Main communication routine. Loop over each coordinate and perform
   the MPI comm. Temporary 1-D buffers are used to send the coordinates at
   all nodes on each partition to the master node. These are then unpacked
   by the master and sorted by marker tag in one large n-dim. array. ---*/

  /*--- Loop over this partition to collect the coords of the local points. ---*/

  su2double *Coords_Local; jPoint = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Only communicate owned nodes to avoid duplicates. ---*/

        if (geometry->node[iPoint]->GetDomain()) {

          /*--- Retrieve local coordinates at this node. ---*/

          Coords_Local = geometry->node[iPoint]->GetCoord();

          /*--- Load local coords into the temporary send buffer. ---*/

          Buffer_Send_X[jPoint] = Coords_Local[0];
          Buffer_Send_Y[jPoint] = Coords_Local[1];
          if (nDim == 3) Buffer_Send_Z[jPoint] = Coords_Local[2];

          /*--- If US system, the output should be in inches ---*/

          if (config->GetSystemMeasurements() == US) {
            Buffer_Send_X[jPoint] *= 12.0;
            Buffer_Send_Y[jPoint] *= 12.0;
            if (nDim == 3) Buffer_Send_Z[jPoint] *= 12.0;
          }

          /*--- Store the marker tag for this particular node. ---*/

          SPRINTF(&Buffer_Send_Str[jPoint*MAX_STRING_SIZE], "%s",
                  config->GetMarker_All_TagBound(iMarker).c_str());

          /*--- Increment jPoint as the counter. We need this because iPoint
           may include halo nodes that we skip over during this loop. ---*/
          
          jPoint++;
          
        }
      }
    }
  }

  /*--- Gather the coordinate data on the master node using MPI. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_X, MaxLocalPoint, MPI_DOUBLE, Buffer_Recv_X, MaxLocalPoint, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Y, MaxLocalPoint, MPI_DOUBLE, Buffer_Recv_Y, MaxLocalPoint, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) {
    SU2_MPI::Gather(Buffer_Send_Z, MaxLocalPoint, MPI_DOUBLE, Buffer_Recv_Z, MaxLocalPoint, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  }
  SU2_MPI::Gather(Buffer_Send_Str, MaxLocalPoint*MAX_STRING_SIZE, MPI_CHAR, Buffer_Recv_Str, MaxLocalPoint*MAX_STRING_SIZE, MPI_CHAR, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < MaxLocalPoint; iPoint++) {
    Buffer_Recv_X[iPoint] = Buffer_Send_X[iPoint];
    Buffer_Recv_Y[iPoint] = Buffer_Send_Y[iPoint];
    if (nDim == 3) Buffer_Recv_Z[iPoint] = Buffer_Send_Z[iPoint];
    index = iPoint*MAX_STRING_SIZE;
    for (iChar = 0; iChar < MAX_STRING_SIZE; iChar++) {
      Buffer_Recv_Str[index + iChar] = Buffer_Send_Str[index + iChar];
    }
  }
#endif

  /*--- The master node unpacks and sorts this variable by marker tag. ---*/

  if (rank == MASTER_NODE) {

    Marker_Tags_InletFile.clear();

    /*--- First, parse the marker tags to count how many total inlet markers
     we have now on the master. ---*/

    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        index = (iProcessor*MaxLocalPoint + iPoint)*MAX_STRING_SIZE;
        for (iChar = 0; iChar < MAX_STRING_SIZE; iChar++) {
          str_buf[iChar] = Buffer_Recv_Str[index + iChar];
        }
        Marker_Tags.push_back(str_buf);
        Marker_Tags_InletFile.push_back(str_buf);
      }
    }

    /*--- Sort and remove the duplicate inlet marker strings. ---*/

    sort(Marker_Tags_InletFile.begin(), Marker_Tags_InletFile.end());
    Marker_Tags_InletFile.erase(unique(Marker_Tags_InletFile.begin(),
                                       Marker_Tags_InletFile.end()),
                                Marker_Tags_InletFile.end());

    /*--- Store the unique number of markers for writing later. ---*/

    nMarker_InletFile = Marker_Tags_InletFile.size();

    /*--- Count the number of rows (nodes) per marker. ---*/

    nRow_InletFile = new unsigned long[nMarker_InletFile];
    for (iMarker = 0; iMarker < nMarker_InletFile; iMarker++) {
      nRow_InletFile[iMarker] = 0;
    }

    /*--- Now count the number of points per marker. ---*/

    jPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        for (iMarker = 0; iMarker < nMarker_InletFile; iMarker++) {
          if (Marker_Tags_InletFile[iMarker] == Marker_Tags[jPoint]) {
            nRow_InletFile[iMarker]++;
          }
        }
        jPoint++;
      }
    }

    /*--- Now put the number of points per marker into cumulative storage.
     We will also create an extra counter to make sorting easier. ---*/

    nRowCum_InletFile = new unsigned long[nMarker_InletFile+1];
    nRowCum_Counter   = new unsigned long[nMarker_InletFile+1];

    nRowCum_InletFile[0] = 0; nRowCum_Counter[0] = 0;
    for (iMarker = 0; iMarker < nMarker_InletFile; iMarker++) {
      nRowCum_InletFile[iMarker+1] = nRowCum_InletFile[iMarker] + nRow_InletFile[iMarker];
      nRowCum_Counter[iMarker+1]   = nRowCum_Counter[iMarker]   + nRow_InletFile[iMarker];
    }

    /*--- Load up the coordinates, sorted into chunks per marker. ---*/

    jPoint = 0; kPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        for (iMarker = 0; iMarker < nMarker_InletFile; iMarker++) {
          if (Marker_Tags_InletFile[iMarker] == Marker_Tags[kPoint]) {

            /*--- Find our current index for this marker and store coords. ---*/

            index = nRowCum_Counter[iMarker];
            InletCoords[0][index] = Buffer_Recv_X[jPoint];
            InletCoords[1][index] = Buffer_Recv_Y[jPoint];
            if (nDim == 3) InletCoords[2][index] = Buffer_Recv_Z[jPoint];

            /*--- Increment the counter for this marker. ---*/

            nRowCum_Counter[iMarker]++;

          }
        }

        /*--- Increment point counter for marker tags and data. ---*/

        kPoint++;
        jPoint++;
        
      }

      /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

      jPoint = (iProcessor+1)*MaxLocalPoint;

    }
  }

  /*--- Immediately release the temporary data buffers. ---*/

  delete [] Buffer_Send_X;
  delete [] Buffer_Send_Y;
  if (Buffer_Send_Z != NULL) delete [] Buffer_Send_Z;
  delete [] Buffer_Send_Str;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_X;
    delete [] Buffer_Recv_Y;
    if (Buffer_Recv_Z != NULL)  delete [] Buffer_Recv_Z;
    delete [] Buffer_Recv_nPoin;
    delete [] Buffer_Recv_Str;
    delete [] nRowCum_Counter;
  }

}

void COutput::Write_InletFile_Flow(CConfig *config, CGeometry *geometry, CSolver **solver) {

  unsigned short iMarker, iDim, iVar;
  unsigned long iPoint;
  su2double turb_val[2] = {0.0,0.0};

  const unsigned short nDim = geometry->GetnDim();

  bool turbulent = (config->GetKind_Solver() == RANS ||
                    config->GetKind_Solver() == ADJ_RANS ||
                    config->GetKind_Solver() == DISC_ADJ_RANS);

  unsigned short nVar_Turb = 0;
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
        nVar_Turb = 1;
        turb_val[0] = solver[TURB_SOL]->GetNuTilde_Inf();
        break;
      case SST:
        nVar_Turb = 2;
        turb_val[0] = solver[TURB_SOL]->GetTke_Inf();
        turb_val[1] = solver[TURB_SOL]->GetOmega_Inf();
        break;
      default:
        SU2_MPI::Error("Specified turbulence model unavailable or none selected", CURRENT_FUNCTION);
        break;
    }

  /*--- Count the number of columns that we have for this flow case.
   Here, we have nDim entries for node coordinates, 2 entries for the total
   conditions or mass flow, another nDim for the direction vector, and
   finally entries for the number of turbulence variables. ---*/

  unsigned short nCol_InletFile = nDim + 2 + nDim + nVar_Turb;

  /*--- Write the inlet profile file. Note that we have already merged
   all of the information for the markers and coordinates previously
   in the MergeInletCoordinates() routine. ---*/

  ofstream node_file("inlet_example.dat");

  node_file << "NMARK= " << nMarker_InletFile << endl;

  for (iMarker = 0; iMarker < nMarker_InletFile; iMarker++) {

    /*--- Access the default data for this marker. ---*/

    string Marker_Tag   = Marker_Tags_InletFile[iMarker];
    su2double p_total   = config->GetInlet_Ptotal(Marker_Tag);
    su2double t_total   = config->GetInlet_Ttotal(Marker_Tag);
    su2double* flow_dir = config->GetInlet_FlowDir(Marker_Tag);

    /*--- Header information for this marker. ---*/

    node_file << "MARKER_TAG= " << Marker_Tag              << endl;
    node_file << "NROW="        << nRow_InletFile[iMarker] << endl;
    node_file << "NCOL="        << nCol_InletFile          << endl;

    node_file << setprecision(15);
    node_file << std::scientific;

    /*--- Loop over the data structure and write the coords and vars to file. ---*/

    for (iPoint = nRowCum_InletFile[iMarker]; iPoint < nRowCum_InletFile[iMarker+1]; iPoint++) {

      for (iDim = 0; iDim < nDim; iDim++) {
        node_file << InletCoords[iDim][iPoint] << "\t";
      }
      node_file << t_total << "\t" << p_total;
      for (iDim = 0; iDim < nDim; iDim++) {
        node_file << "\t" << flow_dir[iDim];
      }
      for (iVar = 0; iVar < nVar_Turb; iVar++) {
        node_file << "\t" << turb_val[iVar];
      }
      node_file << endl;
    }

  }
  node_file.close();

  /*--- Print a message to inform the user about the template file. ---*/
  
  stringstream err;
  err << endl;
  err << "  Created a template inlet profile file with node coordinates" << endl;
  err << "  and solver variables at `inlet_example.dat`." << endl;
  err << "  You can use this file as a guide for making your own inlet" << endl;
  err << "  specification." << endl << endl;
  SU2_MPI::Error(err.str(), CURRENT_FUNCTION);

}

void COutput::DeallocateInletCoordinates(CConfig *config, CGeometry *geometry) {

  unsigned short iDim, nDim = geometry->GetnDim();

  /*--- The master node alone owns all data found in this routine. ---*/

  if (rank == MASTER_NODE) {

    /*--- Deallocate memory for inlet coordinate data ---*/

    if (nRow_InletFile != NULL)    delete [] nRow_InletFile;
    if (nRowCum_InletFile != NULL) delete [] nRowCum_InletFile;

    Marker_Tags_InletFile.clear();

    for (iDim = 0; iDim < nDim; iDim++) {
      if (InletCoords[iDim] != NULL) delete [] InletCoords[iDim];
    }
    if (InletCoords != NULL) delete [] InletCoords;
  }

}



void COutput::MergeConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {

  /*--- Flags identifying the types of files to be written. ---*/

  bool Wrt_Vol = config->GetWrt_Vol_Sol();
  bool Wrt_Srf = config->GetWrt_Srf_Sol();

  /*--- Merge connectivity for each type of element (excluding halos). Note
   that we only need to merge the connectivity once, as it does not change
   during computation. Check whether the base file has been written. ---*/

  /*--- Merge volumetric grid. ---*/

  if (Wrt_Vol) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tria != 0))
      cout <<"Merging volumetric triangle grid connectivity." << endl;
    MergeVolumetricConnectivity_FEM(config, geometry, TRIANGLE    );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Quad != 0))
      cout <<"Merging volumetric quadrilateral grid connectivity." << endl;
    MergeVolumetricConnectivity_FEM(config, geometry, QUADRILATERAL   );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tetr != 0))
      cout <<"Merging volumetric tetrahedron grid connectivity." << endl;
    MergeVolumetricConnectivity_FEM(config, geometry, TETRAHEDRON );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Hexa != 0))
      cout <<"Merging volumetric hexahedron grid connectivity." << endl;
    MergeVolumetricConnectivity_FEM(config, geometry, HEXAHEDRON  );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Pris != 0))
      cout <<"Merging volumetric prism grid connectivity." << endl;
    MergeVolumetricConnectivity_FEM(config, geometry, PRISM       );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Pyra != 0))
      cout <<"Merging volumetric pyramid grid connectivity." << endl;
    MergeVolumetricConnectivity_FEM(config, geometry, PYRAMID     );

  }

  /*--- Merge surface grid. ---*/

  if (Wrt_Srf) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Line != 0))
      cout <<"Merging surface line grid connectivity." << endl;
    MergeSurfaceConnectivity_FEM(config, geometry, LINE);

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundTria != 0))
      cout <<"Merging surface triangle grid connectivity." << endl;
    MergeSurfaceConnectivity_FEM(config, geometry, TRIANGLE);

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundQuad != 0))
      cout <<"Merging surface quadrilateral grid connectivity." << endl;
    MergeSurfaceConnectivity_FEM(config, geometry, QUADRILATERAL);

  }

  /*--- Update total number of volume elements after merge. ---*/

  nGlobal_Elem = nGlobal_Tria + nGlobal_Quad + nGlobal_Tetr +
  nGlobal_Hexa + nGlobal_Pyra + nGlobal_Pris;

  /*--- Update total number of surface elements after merge. ---*/

  nSurf_Elem = nGlobal_Line + nGlobal_BoundTria + nGlobal_BoundQuad;

}

void COutput::MergeCoordinates_FEM(CConfig *config, CGeometry *geometry) {

  /*--- Local variables needed on all processors ---*/

  unsigned short iDim;
  unsigned long iPoint;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned short nDim = DGGeometry->GetnDim();

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem = DGGeometry->GetVolElem();

  /*--- Create the map from the global DOF ID to the local index. ---*/
  vector<su2double> DOFsCoords;
  vector<unsigned long> globalID;

  /*--- Set the global ID and the coordinates of the local DOFs
        by looping over the owned volume elements. ---*/
  unsigned long nLocalPoint = 0;
  for(unsigned long l=0; l<nVolElemOwned; ++l) {
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {

      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      globalID.push_back(globalIndex);

      const su2double *coor = volElem[l].coorSolDOFs.data() + j*nDim;
      for(unsigned short k=0; k<nDim; ++k) DOFsCoords.push_back(coor[k]);

      nLocalPoint++;
    }
  }

#ifndef HAVE_MPI

  /*--- In serial, the single process has access to all geometry, so simply
   load the coordinates into the data structure. ---*/

  // need to double check for halos in parallel

  nGlobal_Poin = nLocalPoint;
  nGlobal_Doma = nLocalPoint;

  /*--- Allocate the coordinates data structure. ---*/

  Coords = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Coords[iDim] = new su2double[nGlobal_Poin];
  }

  /*--- Loop over the mesh to collect the coords of the local points ---*/

  for (iPoint = 0; iPoint < nLocalPoint; iPoint++) {

    /*--- Check if the node belongs to the domain (i.e, not a halo node).
     Sort by the global index, even in serial there is a renumbering (e.g. RCM). ---*/

    /*--- Retrieve the current coordinates at this node. ---*/

    unsigned long iGlobal_Index = globalID[iPoint];

    for (iDim = 0; iDim < nDim; iDim++) {
      Coords[iDim][iGlobal_Index] = DOFsCoords[iPoint*nDim+iDim];

      /*--- If US system, the output should be in inches ---*/

      if ((config->GetSystemMeasurements() == US) && (config->GetKind_SU2() != SU2_DEF)) {
        Coords[iDim][iGlobal_Index] *= 12.0;
      }
    }
  }

#else

  /*--- MPI preprocessing ---*/
  int iProcessor, nProcessor = size;
  unsigned long jPoint;

  /*--- Local variables needed for merging the geometry with MPI. ---*/

  unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
  unsigned long MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;

  if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];

  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/

  // Need to double check halo layers

  Buffer_Send_nPoin[0] = nLocalPoint;

  /*--- Communicate the total number of nodes on this domain. ---*/

  SU2_MPI::Gather(&Buffer_Send_nPoin, 1, MPI_UNSIGNED_LONG,
                  Buffer_Recv_nPoin, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

  if (rank == MASTER_NODE) {
    nGlobal_Doma = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Doma += Buffer_Recv_nPoin[iProcessor];
    }
  }
  nBuffer_Scalar = MaxLocalPoint;

  /*--- Send and Recv buffers. ---*/

  su2double *Buffer_Send_X = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_X = NULL;

  su2double *Buffer_Send_Y = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Y = NULL;

  su2double *Buffer_Send_Z = NULL, *Buffer_Recv_Z = NULL;
  if (nDim == 3) Buffer_Send_Z = new su2double[MaxLocalPoint];

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Prepare the receive buffers in the master node only. ---*/

  if (rank == MASTER_NODE) {

    Buffer_Recv_X = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Y = new su2double[nProcessor*MaxLocalPoint];
    if (nDim == 3) Buffer_Recv_Z = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];

    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoin[iProcessor];
    }
    Coords = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Coords[iDim] = new su2double[nGlobal_Poin];
    }
  }

  /*--- Main communication routine. Loop over each coordinate and perform
   the MPI comm. Temporary 1-D buffers are used to send the coordinates at
   all nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/

  /*--- Loop over this partition to collect the coords of the local points. ---*/
  jPoint = 0;
  for (iPoint = 0; iPoint < nLocalPoint; iPoint++) {

    /*--- Check for halos and write only if requested ---*/
    //    if (!Local_Halo[iPoint] || Wrt_Halo) {

    /*--- Load local coords into the temporary send buffer. These were stored above ---*/
    Buffer_Send_X[jPoint] = DOFsCoords[iPoint*nDim+0];
    Buffer_Send_Y[jPoint] = DOFsCoords[iPoint*nDim+1];
    if (nDim == 3) Buffer_Send_Z[jPoint] = DOFsCoords[iPoint*nDim+2];

    /*--- If US system, the output should be in inches ---*/

    if ((config->GetSystemMeasurements() == US) && (config->GetKind_SU2() != SU2_DEF)) {
      Buffer_Send_X[jPoint] *= 12.0;
      Buffer_Send_Y[jPoint] *= 12.0;
      if (nDim == 3) Buffer_Send_Z[jPoint] *= 12.0;
    }

    /*--- Store the global index for this local node. ---*/
    Buffer_Send_GlobalIndex[jPoint] = globalID[iPoint];

    /*--- Increment jPoint as the counter. We need this because iPoint
     may include halo nodes that we skip over during this loop. ---*/
    jPoint++;
    //    }
  }

  /*--- Gather the coordinate data on the master node using MPI. ---*/

  SU2_MPI::Gather(Buffer_Send_X, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_X, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) {
    SU2_MPI::Gather(Buffer_Send_Z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  }
  SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

  /*--- The master node unpacks and sorts this variable by global index ---*/

  if (rank == MASTER_NODE) {
    jPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        /*--- Get global index, then loop over each variable and store ---*/
        iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
        if (iGlobal_Index >= nGlobal_Poin){
          cout << iGlobal_Index << " " << nGlobal_Poin << endl;
        }
        Coords[0][iGlobal_Index] = Buffer_Recv_X[jPoint];
        Coords[1][iGlobal_Index] = Buffer_Recv_Y[jPoint];
        if (nDim == 3) Coords[2][iGlobal_Index] = Buffer_Recv_Z[jPoint];
        jPoint++;
      }
      /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
      jPoint = (iProcessor+1)*nBuffer_Scalar;
    }
  }

  /*--- Immediately release the temporary data buffers. ---*/

  //  delete [] Local_Halo;
  delete [] Buffer_Send_X;
  delete [] Buffer_Send_Y;
  if (Buffer_Send_Z != NULL) delete [] Buffer_Send_Z;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_X;
    delete [] Buffer_Recv_Y;
    if (Buffer_Recv_Z != NULL)  delete [] Buffer_Recv_Z;
    delete [] Buffer_Recv_GlobalIndex;
    delete [] Buffer_Recv_nPoin;
  }

#endif

}

void COutput::MergeVolumetricConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

  int iProcessor;
  unsigned short NODES_PER_ELEMENT = 0;
  unsigned long iNode, jNode;
  unsigned long iElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;

  unsigned long jElem;

  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0;
  unsigned long MaxLocalElem = 0;

  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem = NULL;

  int *Conn_Elem = NULL;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();

  CVolumeElementFEM *volElem = DGGeometry->GetVolElem();

  const CFEMStandardElement *standardElementsSol = DGGeometry->GetStandardElementsSol();

  /* Define the vectors for the connectivity of the local linear subelements. */
  vector<unsigned long> volumeConn;

  /*--- Create the map from the global DOF ID to the local index. ---*/
  map<unsigned long, unsigned long> mapLocal2Global;
  unsigned long ii = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j, ++ii) {
      mapLocal2Global[ii] = volElem[i].offsetDOFsSolGlobal+j;
    }
  }

  /*--- Counter for keeping track of the number of this element locally. ---*/
  nLocalElem = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /* Determine the necessary data from the corresponding standard elem,
     such as the number of linear sub elements, the number of DOFs per
     linear sub element and the corresponding local connectivity. */
    const unsigned short ind       = volElem[i].indStandardElement;
    const unsigned short VTK_Type1 = standardElementsSol[ind].GetVTK_Type1();
    const unsigned short VTK_Type2 = standardElementsSol[ind].GetVTK_Type2();

    /*--- Only store the linear sub elements if they are of
     the current type that we are merging. ---*/
    if (Elem_Type == VTK_Type1) {

      const unsigned short nSubElems       = standardElementsSol[ind].GetNSubElemsType1();
      const unsigned short nDOFsPerSubElem = standardElementsSol[ind].GetNDOFsPerSubElem(VTK_Type1);
      const unsigned short *connSubElems   = standardElementsSol[ind].GetSubConnType1();

      /* Loop over the number of subfaces and store the required data. */
      unsigned short kk = 0;
      for(unsigned short j=0; j<nSubElems; ++j) {

        /*--- Store the global index for the surface conn ---*/
        for(unsigned short k=0; k<nDOFsPerSubElem; ++k, ++kk)
          volumeConn.push_back(mapLocal2Global[volElem[i].offsetDOFsSolLocal+connSubElems[kk]]);

        nLocalElem++;
      }

    } else if (Elem_Type == VTK_Type2) {

      const unsigned short nSubElems       = standardElementsSol[ind].GetNSubElemsType2();
      const unsigned short nDOFsPerSubElem = standardElementsSol[ind].GetNDOFsPerSubElem(VTK_Type2);
      const unsigned short *connSubElems   = standardElementsSol[ind].GetSubConnType2();

      /* Loop over the number of subfaces and store the required data. */
      unsigned short kk = 0;
      for(unsigned short j=0; j<nSubElems; ++j) {

        /*--- Store the global index for the surface conn ---*/
        for(unsigned short k=0; k<nDOFsPerSubElem; ++k, ++kk)
          volumeConn.push_back(mapLocal2Global[volElem[i].offsetDOFsSolLocal+connSubElems[kk]]);

        nLocalElem++;
      }
    }
  }

  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/

  switch (Elem_Type) {
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case PRISM:
      NODES_PER_ELEMENT = N_POINTS_PRISM;
      break;
    case PYRAMID:
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
  }

  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/

  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[size];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalElem, &MaxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalElem = nLocalElem;
  Buffer_Recv_nElem[0] = Buffer_Send_nElem[0];
#endif

  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;

  /*--- Send and Recv buffers ---*/

  unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
  unsigned long *Buffer_Recv_Elem = NULL;

  unsigned short *Buffer_Send_Halo = new unsigned short[MaxLocalElem];
  unsigned short *Buffer_Recv_Halo = NULL;

  /*--- Prepare the receive buffers on the master node only. ---*/

  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new unsigned long[size*nBuffer_Scalar];
    Buffer_Recv_Halo = new unsigned short[size*MaxLocalElem];
    if (MaxLocalElem > 0) Conn_Elem = new int[size*MaxLocalElem*NODES_PER_ELEMENT];
  }

  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/

  jNode = 0;
  for (iElem = 0; iElem < nLocalElem; iElem++) {

    /*--- Loop over all nodes in this element and load the
     connectivity into the send buffer. ---*/
    Buffer_Send_Halo[iElem] = false;
    for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {

      /*--- Store the global index values directly. ---*/

      Buffer_Send_Elem[jNode] = volumeConn[iElem*NODES_PER_ELEMENT+iNode];

      /*--- Increment jNode as the counter. We need this because iElem
       may include other elements that we skip over during this loop. ---*/

      jNode++;
    }
  }

  /*--- Gather the element connectivity information. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, Buffer_Recv_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (unsigned long iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Elem[iPoint] = Buffer_Send_Elem[iPoint];
  for (unsigned long iPoint = 0; iPoint < MaxLocalElem; iPoint++) Buffer_Recv_Halo[iPoint] = Buffer_Send_Halo[iPoint];
#endif

  /*--- The master node unpacks and sorts the connectivity. ---*/

  if (rank == MASTER_NODE) {

    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/

    Write_Elem = new bool[size*MaxLocalElem];
    for (iElem = 0; iElem < size*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }

    /*--- Remove the rind layer from the solution only if requested ---*/

    if (!Wrt_Halo) {

      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/

      kElem = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

          /*--- Check if this element was marked as a halo. ---*/
          if (Buffer_Recv_Halo[kElem+iElem])
            Write_Elem[kElem+iElem] = false;

        }
        kElem = (iProcessor+1)*MaxLocalElem;
      }
    }

    /*--- Store the unique connectivity list for this element type. ---*/

    jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {

          /*--- Increment total count for this element type ---*/
          nElem_Total++;

          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/

          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = (int)Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
            kNode++;
          }
        }
      }
      /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
      jElem = (iProcessor+1)*MaxLocalElem;
      jNode = (iProcessor+1)*nBuffer_Scalar;
    }
  }

  /*--- Immediately release the temporary buffers. ---*/
  delete [] Buffer_Send_Elem;
  delete [] Buffer_Send_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case TRIANGLE:
        nGlobal_Tria = nElem_Total;
        if (nGlobal_Tria > 0) Conn_Tria = Conn_Elem;
        break;
      case QUADRILATERAL:
        nGlobal_Quad = nElem_Total;
        if (nGlobal_Quad > 0) Conn_Quad = Conn_Elem;
        break;
      case TETRAHEDRON:
        nGlobal_Tetr = nElem_Total;
        if (nGlobal_Tetr > 0) Conn_Tetr = Conn_Elem;
        break;
      case HEXAHEDRON:
        nGlobal_Hexa = nElem_Total;
        if (nGlobal_Hexa > 0) Conn_Hexa = Conn_Elem;
        break;
      case PRISM:
        nGlobal_Pris = nElem_Total;
        if (nGlobal_Pris > 0) Conn_Pris = Conn_Elem;
        break;
      case PYRAMID:
        nGlobal_Pyra = nElem_Total;
        if (nGlobal_Pyra > 0) Conn_Pyra = Conn_Elem;
        break;
      default:
        SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
    }
  }
  
}

void COutput::MergeSurfaceConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

  unsigned short NODES_PER_ELEMENT = 0;

  unsigned long iNode, jNode;
  unsigned long iElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;

  int iProcessor;
  unsigned long jElem;

  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0;
  unsigned long MaxLocalElem = 0;

  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem = NULL;

  int *Conn_Elem = NULL;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();

  CVolumeElementFEM *volElem = DGGeometry->GetVolElem();

  const CBoundaryFEM *boundaries = DGGeometry->GetBoundaries();

  const CFEMStandardBoundaryFace *standardBoundaryFacesSol = DGGeometry->GetStandardBoundaryFacesSol();

  /*--- Create the map from the global DOF ID to the local index. ---*/
  map<unsigned long, unsigned long> mapLocal2Global;
  unsigned long ii = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j, ++ii) {
      mapLocal2Global[ii] = volElem[i].offsetDOFsSolGlobal+j;
    }
  }

  /* Define the vectors for the connectivity of the local linear subelements. */
  vector<unsigned long>  surfaceConn;

  /*--- Counter for keeping track of the number of this element locally. ---*/
  nLocalElem = 0;

  /* Loop over the boundary markers. */
  for(unsigned short iMarker=0; iMarker < config->GetnMarker_All(); ++iMarker) {
    if( !boundaries[iMarker].periodicBoundary ) {

      /* Check for markers to be plotted. */
      if (config->GetMarker_All_Plotting(iMarker) == YES) {

        /* Loop over the surface elements of this marker. */
        const vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
        for(unsigned long i=0; i<surfElem.size(); ++i) {

          /* Determine the necessary data from the corresponding standard face,
           such as the number of linear subfaces, the number of DOFs per
           linear subface and the corresponding local connectivity. */
          const unsigned short ind           = surfElem[i].indStandardElement;
          const unsigned short VTK_Type      = standardBoundaryFacesSol[ind].GetVTK_Type();

          /*--- Only store the linear sub-elements if they are of
           the current type that we are merging. ---*/
          if (VTK_Type == Elem_Type) {

            const unsigned short nSubFaces     = standardBoundaryFacesSol[ind].GetNSubFaces();
            const unsigned short nDOFsPerFace  = standardBoundaryFacesSol[ind].GetNDOFsPerSubFace();
            const unsigned short *connSubFaces = standardBoundaryFacesSol[ind].GetSubFaceConn();

            /* Loop over the number of subfaces and store the required data. */
            unsigned short ii = 0;
            for(unsigned short j=0; j<nSubFaces; ++j) {

              /*--- Store the global index for the surface conn ---*/
              for(unsigned short k=0; k<nDOFsPerFace; ++k, ++ii)
                surfaceConn.push_back(mapLocal2Global[surfElem[i].DOFsSolFace[connSubFaces[ii]]]);

              nLocalElem++;
            }
          }
        }
      }
    }
  }

  switch (Elem_Type) {
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
  }

  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/

  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[size];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalElem, &MaxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalElem = nLocalElem;
  Buffer_Recv_nElem[0] = Buffer_Send_nElem[0];
#endif

  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;

  /*--- Send and Recv buffers ---*/

  unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
  unsigned long *Buffer_Recv_Elem = NULL;

  unsigned short *Buffer_Send_Halo = new unsigned short[MaxLocalElem];
  unsigned short *Buffer_Recv_Halo = NULL;

  /*--- Prepare the receive buffers on the master node only. ---*/

  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new unsigned long[size*nBuffer_Scalar];
    Buffer_Recv_Halo = new unsigned short[size*MaxLocalElem];
    if (MaxLocalElem > 0) Conn_Elem = new int[size*MaxLocalElem*NODES_PER_ELEMENT];
  }

  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/
  jNode = 0;
  for (iElem = 0; iElem < nLocalElem; iElem++) {

    /*--- Loop over all nodes in this element and load the
     connectivity into the send buffer. ---*/
    Buffer_Send_Halo[iElem] = false;
    for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {

      /*--- Store the global index values directly. ---*/
      Buffer_Send_Elem[jNode] = surfaceConn[iElem*NODES_PER_ELEMENT+iNode];

      /*--- Increment jNode as the counter for the buffer array. ---*/
      jNode++;
    }
  }

  /*--- Gather the element connectivity information. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, Buffer_Recv_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (unsigned long iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Elem[iPoint] = Buffer_Send_Elem[iPoint];
  for (unsigned long iPoint = 0; iPoint < MaxLocalElem; iPoint++) Buffer_Recv_Halo[iPoint] = Buffer_Send_Halo[iPoint];
#endif

  /*--- The master node unpacks and sorts the connectivity. ---*/

  if (rank == MASTER_NODE) {

    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/

    Write_Elem = new bool[size*MaxLocalElem];
    for (iElem = 0; iElem < size*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }

    /*--- Remove the rind layer from the solution only if requested ---*/

    if (!Wrt_Halo) {

      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/

      kElem = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

          /*--- Check if this element was marked as a halo. ---*/
          if (Buffer_Recv_Halo[kElem+iElem])
            Write_Elem[kElem+iElem] = false;

        }
        kElem = (iProcessor+1)*MaxLocalElem;
      }
    }

    /*--- Store the unique connectivity list for this element type. ---*/

    jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {

          /*--- Increment total count for this element type ---*/
          nElem_Total++;

          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/

          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = (int)Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
            kNode++;
          }
        }
      }
      /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
      jElem = (iProcessor+1)*MaxLocalElem;
      jNode = (iProcessor+1)*nBuffer_Scalar;
    }
  }

  /*--- Immediately release the temporary buffers. ---*/
  delete [] Buffer_Send_Elem;
  delete [] Buffer_Send_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case LINE:
        nGlobal_Line = nElem_Total;
        if (nGlobal_Line > 0) Conn_Line = Conn_Elem;
        break;
      case TRIANGLE:
        nGlobal_BoundTria = nElem_Total;
        if (nGlobal_BoundTria > 0) Conn_BoundTria = Conn_Elem;
        break;
      case QUADRILATERAL:
        nGlobal_BoundQuad = nElem_Total;
        if (nGlobal_BoundQuad > 0) Conn_BoundQuad = Conn_Elem;
        break;
      default:
        SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
    }
  }
}

void COutput::MergeSolution_FEM(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {

  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar = 0, FirstIndex = NONE, SecondIndex = NONE, ThirdIndex = NONE;
  unsigned short nVar_First = 0, nVar_Second = 0, nVar_Third = 0;

  unsigned long iPoint = 0, jPoint = 0;

  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;

  int iProcessor;

  /*--- Prepare send buffers for the conservative variables. Need to
   find the total number of conservative variables and also the
   index for their particular solution container. ---*/

  switch (Kind_Solver) {
    case FEM_EULER : case FEM_NAVIER_STOKES: case FEM_LES: FirstIndex = FLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    default: SecondIndex = NONE; ThirdIndex = NONE; break;
  }

  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  if (ThirdIndex != NONE) nVar_Third = solver[ThirdIndex]->GetnVar();
  nVar_Consv = nVar_First + nVar_Second + nVar_Third;
  nVar_Total = nVar_Consv;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem = DGGeometry->GetVolElem();

  /*--- Create the map from the global DOF ID to the local index. ---*/
  map<unsigned long, unsigned long> mapLocal2Global;
  vector<su2double> DOFsSol;
  vector<unsigned long> globalID;

  /*--- Update the solution by looping over the owned volume elements. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Set the pointers for the residual and solution for this element. */
    const unsigned long offset  = nVar_Consv*volElem[l].offsetDOFsSolLocal;
    su2double *solDOFs          = solver[FirstIndex]->GetVecSolDOFs() + offset;

    /* Loop over the DOFs for this element and store the solution. */

    unsigned int i = 0;
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {

      nLocalPoint++;
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      globalID.push_back(globalIndex);

      for(unsigned short iVar=0; iVar<nVar_Consv; ++iVar, ++i) {
        DOFsSol.push_back(solDOFs[i]);
      }
    }
  }
  Buffer_Send_nPoint[0] = nLocalPoint;

  /*--- Each processor sends its local number of nodes to the master. ---*/

  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[size];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalPoint = nLocalPoint;
  Buffer_Recv_nPoint[0] = Buffer_Send_nPoint[0];
#endif

  nBuffer_Scalar = MaxLocalPoint;

  /*--- Send and Recv buffers. ---*/

  su2double *Buffer_Send_Var = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Var = NULL;

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Prepare the receive buffers in the master node only. ---*/

  if (rank == MASTER_NODE) {

    Buffer_Recv_Var = new su2double[size*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[size*MaxLocalPoint];

    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new su2double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new su2double[nGlobal_Poin];
    }
  }

  /*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/

  for (iVar = 0; iVar < nVar_Consv; iVar++) {

    /*--- Loop over this partition to collect the current variable ---*/

    jPoint = 0;
    for (iPoint = 0; iPoint < nLocalPoint; iPoint++) {

      /*--- Get this variable into the temporary send buffer. ---*/

      Buffer_Send_Var[jPoint] = DOFsSol[iPoint*nVar_Consv+iVar];

      /*--- Only send/recv the volumes & global indices during the first loop ---*/

      if (iVar == 0) {
        Buffer_Send_GlobalIndex[jPoint] = globalID[iPoint];
      }

      jPoint++;

    }

    /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
    SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

    if (iVar == 0) {
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_GlobalIndex[iPoint] = Buffer_Send_GlobalIndex[iPoint];
#endif
    }

    /*--- The master node unpacks and sorts this variable by global index ---*/

    if (rank == MASTER_NODE) {
      jPoint = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

          /*--- Get global index, then loop over each variable and store ---*/

          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];

          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];

          jPoint++;

        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }

  }

  /*--- Immediately release the temporary buffers. ---*/

  delete [] Buffer_Send_Var;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nPoint;
    delete [] Buffer_Recv_Var;
    delete [] Buffer_Recv_GlobalIndex;
  }

}

void COutput::MergeBaselineSolution_FEM(CConfig *config, CGeometry *geometry, CSolver *solver, unsigned short val_iZone) {

  /*--- Local variables needed on all processors ---*/
  unsigned short iVar;
  unsigned long iPoint = 0, jPoint = 0;

  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;

  int iProcessor;

  /*--- We know the number of fields from reading the restart. ---*/

  nVar_Total = config->fields.size() - 1;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem = DGGeometry->GetVolElem();

  /*--- Create the map from the global DOF ID to the local index. ---*/
  vector<su2double> DOFsSol;
  vector<unsigned long> globalID;

  /*--- Update the solution by looping over the owned volume elements. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Set the pointers for the residual and solution for this element. */
    const unsigned long offset  = nVar_Total*volElem[l].offsetDOFsSolLocal;
    su2double *solDOFs          = solver->GetVecSolDOFs() + offset;

    /* Loop over the DOFs for this element and store the solution. */

    unsigned int i = 0;
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {

      nLocalPoint++;
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      globalID.push_back(globalIndex);

      for(unsigned short iVar=0; iVar<nVar_Total; ++iVar, ++i) {
        DOFsSol.push_back(solDOFs[i]);
      }
    }
  }
  Buffer_Send_nPoint[0] = nLocalPoint;

  /*--- Each processor sends its local number of nodes to the master. ---*/

  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[size];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalPoint = nLocalPoint;
  Buffer_Recv_nPoint[0] = Buffer_Send_nPoint[0];
#endif

  nBuffer_Scalar = MaxLocalPoint;

  /*--- Send and Recv buffers. ---*/

  su2double *Buffer_Send_Var = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Var = NULL;

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Prepare the receive buffers in the master node only. ---*/

  if (rank == MASTER_NODE) {

    Buffer_Recv_Var = new su2double[size*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[size*MaxLocalPoint];

    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new su2double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new su2double[nGlobal_Poin];
    }
  }

  /*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/

  for (iVar = 0; iVar < nVar_Total; iVar++) {

    /*--- Loop over this partition to collect the current variable ---*/

    jPoint = 0;
    for (iPoint = 0; iPoint < nLocalPoint; iPoint++) {

      /*--- Get this variable into the temporary send buffer. ---*/

      Buffer_Send_Var[jPoint] = DOFsSol[iPoint*nVar_Total+iVar];

      /*--- Only send/recv the volumes & global indices during the first loop ---*/

      if (iVar == 0) {
        Buffer_Send_GlobalIndex[jPoint] = globalID[iPoint];
      }

      jPoint++;

    }

    /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
    SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

    if (iVar == 0) {
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_GlobalIndex[iPoint] = Buffer_Send_GlobalIndex[iPoint];
#endif
    }

    /*--- The master node unpacks and sorts this variable by global index ---*/

    if (rank == MASTER_NODE) {
      jPoint = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

          /*--- Get global index, then loop over each variable and store ---*/

          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];

          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];

          jPoint++;

        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }

  }

  /*--- Immediately release the temporary buffers. ---*/

  delete [] Buffer_Send_Var;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nPoint;
    delete [] Buffer_Recv_Var;
    delete [] Buffer_Recv_GlobalIndex;
  }
  
}

void COutput::SetResult_Files_FEM(CSolver ****solver_container, CGeometry ***geometry, CConfig **config,
                                  unsigned long iExtIter, unsigned short val_nZone) {

  unsigned short iZone;

  for (iZone = 0; iZone < val_nZone; iZone++) {

    /*--- Check whether we are writing any output at all. We can
     disable all output to avoid serial bottleneck at scale. ---*/
    if (config[iZone]->GetWrt_Output()) {

      /*--- Flags identifying the types of files to be written. ---*/

      bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
      bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();

      /*--- Get the file output format ---*/

      unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();

      /*--- Merge the node coordinates and connectivity, if necessary. This
       is only performed if a volume solution file is requested, and it
       is active by default. ---*/

      if (Wrt_Vol || Wrt_Srf) {
        if (rank == MASTER_NODE) cout << endl << "Merging connectivities in the Master node." << endl;
        MergeConnectivity_FEM(config[iZone], geometry[iZone][MESH_0], iZone);
      }

      /*--- Merge coordinates of all grid nodes (excluding ghost points).
       The grid coordinates are always merged and included first in the
       restart files. ---*/

      if (rank == MASTER_NODE) cout << "Merging coordinates in the Master node." << endl;
      MergeCoordinates_FEM(config[iZone], geometry[iZone][MESH_0]);

      if ((rank == MASTER_NODE) && (Wrt_Vol || Wrt_Srf)) {
        if (FileFormat == TECPLOT_BINARY) {
          if (rank == MASTER_NODE) cout << "Writing Tecplot binary volume and surface mesh files." << endl;
          SetTecplotBinary_DomainMesh(config[iZone], geometry[iZone][MESH_0], iZone);
          SetTecplotBinary_SurfaceMesh(config[iZone], geometry[iZone][MESH_0], iZone);
        }
      }

      /*--- Merge the solution data needed for volume solutions and restarts ---*/

      if (rank == MASTER_NODE) cout << "Merging solution in the Master node." << endl;
      MergeSolution_FEM(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0], iZone);

      /*--- Write restart, or Tecplot files using the merged data.
       This data lives only on the master, and these routines are currently
       executed by the master proc alone (as if in serial). ---*/

      if (rank == MASTER_NODE) {

        /*--- Write a native restart file ---*/

        if (rank == MASTER_NODE) cout << "Writing SU2 native restart file." << endl;
        SetRestart(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0] , iZone);

        if (Wrt_Vol) {

          switch (FileFormat) {

            case TECPLOT:

              /*--- Write a Tecplot ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file volume solution file." << endl;
              SetTecplotASCII(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0], iZone, val_nZone, false);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
              break;

            case FIELDVIEW:

              /*--- Write a FieldView ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing FieldView ASCII file volume solution file." << endl;
              SetFieldViewASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
              break;

            case TECPLOT_BINARY:

              /*--- Write a Tecplot binary solution file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot binary volume solution file." << endl;
              SetTecplotBinary_DomainSolution(config[iZone], geometry[iZone][MESH_0], iZone);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
              break;

            case FIELDVIEW_BINARY:

              /*--- Write a FieldView binary file ---*/

              if (rank == MASTER_NODE) cout << "Writing FieldView binary file volume solution file." << endl;
              SetFieldViewBinary(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
              break;

            case PARAVIEW:

              /*--- Write a Paraview ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Paraview ASCII volume solution file." << endl;
              SetParaview_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, false);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
              break;

            default:
              break;
          }

        }

        if (Wrt_Srf) {

          switch (FileFormat) {

            case TECPLOT:

              /*--- Write a Tecplot ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII surface solution file." << endl;
              SetTecplotASCII(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0] , iZone, val_nZone, true);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
              break;

            case TECPLOT_BINARY:

              /*--- Write a Tecplot binary solution file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot binary surface solution file." << endl;
              SetTecplotBinary_SurfaceSolution(config[iZone], geometry[iZone][MESH_0], iZone);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
              break;

            case PARAVIEW:

              /*--- Write a Paraview ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Paraview ASCII surface solution file." << endl;
              SetParaview_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, true);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
              break;

            default:
              break;
          }

        }

        /*--- Release memory needed for merging the solution data. ---*/

        DeallocateCoordinates(config[iZone], geometry[iZone][MESH_0]);
        DeallocateSolution(config[iZone], geometry[iZone][MESH_0]);

      }

      /*--- Final broadcast (informing other procs that the base output
       file was written). ---*/

#ifdef HAVE_MPI
      SU2_MPI::Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Bcast(&wrote_surf_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif

    } else {
      if (rank == MASTER_NODE) cout << endl << "Restart and solution output disabled." << endl;
    }
  }
}

void COutput::SetBaselineResult_Files_FEM(CSolver ***solver, CGeometry ***geometry, CConfig **config,
                                          unsigned long iExtIter, unsigned short val_nZone) {

  unsigned short iZone, iInst, nInst;

  for (iZone = 0; iZone < val_nZone; iZone++) {

    nInst = config[iZone]->GetnTimeInstances();

    for (iInst = 0; iInst < nInst; iInst++) {

      config[iZone]->SetiInst(iInst);

      /*--- Flags identifying the types of files to be written. ---*/

      bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
      bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();

      /*--- Get the file output format ---*/

      unsigned short FileFormat = config[iZone]->GetOutput_FileFormat();

      /*--- Merge the node coordinates and connectivity if necessary. This
       is only performed if a volume solution file is requested, and it
       is active by default. ---*/

      if ((Wrt_Vol || Wrt_Srf)) {
        if (rank == MASTER_NODE) cout << "Merging connectivities in the Master node." << endl;
        MergeConnectivity_FEM(config[iZone], geometry[iZone][iInst], iZone);
      }

      /*--- Merge the solution data needed for volume solutions and restarts ---*/

      if ((Wrt_Vol || Wrt_Srf)) {
        if (rank == MASTER_NODE) cout << "Merging solution in the Master node." << endl;
        MergeBaselineSolution_FEM(config[iZone], geometry[iZone][iInst], solver[iZone][iInst], iZone);
      }

      /*--- Write restart, Tecplot or Paraview files using the merged data.
       This data lives only on the master, and these routines are currently
       executed by the master proc alone (as if in serial). ---*/

      if (rank == MASTER_NODE) {

        if (Wrt_Vol) {

          switch (FileFormat) {

            case TECPLOT:

              /*--- Write a Tecplot ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (volume grid)." << endl;
              SetTecplotASCII(config[iZone], geometry[iZone][iInst], &solver[iZone][iInst], iZone, val_nZone, false);
              DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
              break;

            case FIELDVIEW:

              /*--- Write a FieldView ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing FieldView ASCII file (volume grid)." << endl;
              SetFieldViewASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone);
              DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
              break;

            case TECPLOT_BINARY:

              /*--- Write a Tecplot binary solution file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot Binary file (volume grid)." << endl;
              SetTecplotBinary_DomainMesh(config[iZone], geometry[iZone][iInst], iZone);
              SetTecplotBinary_DomainSolution(config[iZone], geometry[iZone][iInst], iZone);
              break;

            case FIELDVIEW_BINARY:

              /*--- Write a binary binary file ---*/

              if (rank == MASTER_NODE) cout << "Writing FieldView ASCII file (volume grid)." << endl;
              SetFieldViewBinary(config[iZone], geometry[iZone][iInst], iZone, val_nZone);
              DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
              break;

            case PARAVIEW:

              /*--- Write a Paraview ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Paraview ASCII file (volume grid)." << endl;
              SetParaview_ASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone, false);
              DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
              break;

            default:
              break;
          }

        }

        if (Wrt_Srf) {

          switch (FileFormat) {

            case TECPLOT:

              /*--- Write a Tecplot ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (surface grid)." << endl;
              SetTecplotASCII(config[iZone], geometry[iZone][iInst], &solver[iZone][iInst], iZone, val_nZone, true);
              DeallocateConnectivity(config[iZone], geometry[iZone][iInst], true);
              break;

            case TECPLOT_BINARY:

              /*--- Write a Tecplot binary solution file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot Binary file (surface grid)." << endl;
              SetTecplotBinary_SurfaceMesh(config[iZone], geometry[iZone][iInst], iZone);
              SetTecplotBinary_SurfaceSolution(config[iZone], geometry[iZone][iInst], iZone);
              break;

            case PARAVIEW:

              /*--- Write a Paraview ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Paraview ASCII file (surface grid)." << endl;
              SetParaview_ASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone, true);
              DeallocateConnectivity(config[iZone], geometry[iZone][iInst], true);
              break;

            default:
              break;
          }
        }

        if (FileFormat == TECPLOT_BINARY) {
          if (!wrote_base_file)
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
          if (!wrote_surf_file)
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], wrote_surf_file);
        }

        if (Wrt_Vol || Wrt_Srf)
          DeallocateSolution(config[iZone], geometry[iZone][iInst]);
      }
    
      /*--- Final broadcast (informing other procs that the base output
       file was written). ---*/
    
#ifdef HAVE_MPI
      SU2_MPI::Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif
    }
  }
}

void COutput::LoadLocalData_FEM(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {

  unsigned short iDim;
  unsigned short Kind_Solver = config->GetKind_Solver();
  unsigned short nDim = geometry->GetnDim();

  unsigned long iVar, jVar;
  unsigned long iPoint, jPoint, FirstIndex = NONE, SecondIndex = NONE;
  unsigned long nVar_First = 0, nVar_Second = 0, nVar_Consv_Par = 0;

  stringstream varname;

  /*--- Use a switch statement to decide how many solver containers we have
   in this zone for output. ---*/

  switch (Kind_Solver) {
    case FEM_EULER : case FEM_NAVIER_STOKES: case FEM_LES: FirstIndex = FLOW_SOL; SecondIndex = NONE; break;
    case FEM_RANS : FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; break;
    default: SecondIndex = NONE; break;
  }

  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  nVar_Consv_Par = nVar_First + nVar_Second;

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Register the variables that will be output. To register a  ---*/
  /*---         variable, two things are required. First, increment the    ---*/
  /*---         counter for the number of variables (nVar_Par), which      ---*/
  /*---         controls the size of the data structure allocation, i.e.,  ---*/
  /*---         the number of columns in an nPoint x nVar structure.       ---*/
  /*---         Second, add a name for the variable to the vector that     ---*/
  /*---         holds the string names.                                    ---*/
  /*--------------------------------------------------------------------------*/

  /*--- All output files first need the grid coordinates. ---*/

  GlobalField_Counter  = 1; Variable_Names.push_back("x");
  GlobalField_Counter += 1; Variable_Names.push_back("y");
  if (geometry->GetnDim() == 3) {
    GlobalField_Counter += 1; Variable_Names.push_back("z");
  }

  /*--- At a mininum, the restarts and visualization files need the
   conservative variables, so these follow next. ---*/

  GlobalField_Counter += nVar_Consv_Par;

  Variable_Names.push_back("Density");
  Variable_Names.push_back("X-Momentum");
  Variable_Names.push_back("Y-Momentum");
  if (geometry->GetnDim() == 3) Variable_Names.push_back("Z-Momentum");
  Variable_Names.push_back("Energy");

  /*--- Eventually, turbulence model from second container goes here. ---*/

  if (!config->GetLow_MemoryOutput()) {

    /*--- Add Pressure, Temperature, Cp, Mach. ---*/

    GlobalField_Counter += 1;
    Variable_Names.push_back("Pressure");

    GlobalField_Counter += 3;
    Variable_Names.push_back("Temperature");
    if (config->GetOutput_FileFormat() == PARAVIEW){
      Variable_Names.push_back("Pressure_Coefficient");
    } else {
      Variable_Names.push_back("C<sub>p</sub>");
    }
    Variable_Names.push_back("Mach");

    /*--- New variables get registered here before the end of the loop. ---*/
    
    if (Kind_Solver == FEM_NAVIER_STOKES){
      GlobalField_Counter += 1;
      Variable_Names.push_back("Laminar_Viscosity");
    }
    if ((Kind_Solver == FEM_LES) && (config->GetKind_SGS_Model() != IMPLICIT_LES)){
      GlobalField_Counter += 1;
      Variable_Names.push_back("Eddy_Viscosity");
    }
  }

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/

  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

  /*--- Get a pointer to the fluid model class from the DG-FEM solver
   so that we can access the states below. ---*/

  CFluidModel *DGFluidModel = solver[FLOW_SOL]->GetFluidModel();

  /*--- Allocate the local data structure now that we know how many
   variables are in the output. ---*/

  Local_Data = new su2double*[nLocalPoint_Sort];
  for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++) {
    Local_Data[iPoint] = new su2double[GlobalField_Counter];
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Loop over all grid nodes and load up the desired data for  ---*/
  /*---         the restart and vizualization files. Note that we need to  ---*/
  /*---         increment the iVar variable after each variable load.      ---*/
  /*---         The idea is that we're filling up the columns of field     ---*/
  /*---         data for each iPoint (row) of the data structure.          ---*/
  /*---         This data will then be sorted, communicated, and written   ---*/
  /*---         to files automatically after this routine. Note that the   ---*/
  /*---         ordering of the data loading MUST match the order of the   ---*/
  /*---         variable registration above for the files to be correct.   ---*/
  /*--------------------------------------------------------------------------*/

  jPoint = 0;

  /*--- Access the solution by looping over the owned volume elements. ---*/

  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Set the pointers for the solution for this element. */

    const unsigned long offset = nVar_First*volElem[l].offsetDOFsSolLocal;
    su2double *solDOFs         = solver[FirstIndex]->GetVecSolDOFs() + offset;

    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {

      /*--- Restart the column index with each new point. ---*/

      iVar = 0;

      /*--- Get the conservative variables for this particular DOF. ---*/

      const su2double *U = solDOFs + j*nVar_First;

      /*--- Load the coordinate values of the solution DOFs. ---*/

      const su2double *coor = volElem[l].coorSolDOFs.data() + j*nDim;
      for(unsigned short k=0; k<nDim; ++k) {
        Local_Data[jPoint][iVar] = coor[k];
        iVar++;
      }

      /* Load the conservative variables. */

      for(jVar=0; jVar < nVar_First; ++jVar) {
        Local_Data[jPoint][iVar] = U[jVar];
        iVar++;
      }

      /*--- Prepare the primitive states. ---*/

      const su2double DensityInv = 1.0/U[0];
      su2double vel[3], Velocity2 = 0.0;
      for(iDim=0; iDim<nDim; ++iDim) {
        vel[iDim] = U[iDim+1]*DensityInv;
        Velocity2 += vel[iDim]*vel[iDim];
      }
      su2double StaticEnergy = U[nDim+1]*DensityInv - 0.5*Velocity2;
      DGFluidModel->SetTDState_rhoe(U[0], StaticEnergy);

      /*--- Load data for the pressure, temperature, Cp, and Mach variables. ---*/

      Local_Data[jPoint][iVar] = DGFluidModel->GetPressure(); iVar++;
      Local_Data[jPoint][iVar] = DGFluidModel->GetTemperature(); iVar++;
      Local_Data[jPoint][iVar] = DGFluidModel->GetCp(); iVar++;
      Local_Data[jPoint][iVar] = sqrt(Velocity2)/DGFluidModel->GetSoundSpeed(); iVar++;

      /*--- New variables can be loaded to the Local_Data structure here,
       assuming they were registered above correctly. ---*/

      if (Kind_Solver == FEM_NAVIER_STOKES){
        Local_Data[jPoint][iVar] = DGFluidModel->GetLaminarViscosity(); iVar++;
      }
      if ((Kind_Solver == FEM_LES) && (config->GetKind_SGS_Model() != IMPLICIT_LES)){
        // todo: Export Eddy instead of Laminar viscosity
        Local_Data[jPoint][iVar] = DGFluidModel->GetLaminarViscosity(); iVar++;
      }

      /*--- Increment the point counter. ---*/

      jPoint++;
    }
  }

}

void COutput::PrepareOffsets(CConfig *config, CGeometry *geometry) {

  unsigned long iPoint;

  /*--- Bool to distinguish between the FVM and FEM solvers. ---*/

  unsigned short KindSolver = config->GetKind_Solver();
  bool fem_solver = ((KindSolver == FEM_EULER) ||
                     (KindSolver == FEM_NAVIER_STOKES) ||
                     (KindSolver == FEM_RANS) ||
                     (KindSolver == FEM_LES));

  /*--- Reset point sorting counters ---*/

  nGlobalPoint_Sort = 0;
  nLocalPoint_Sort  = 0;

  /*--- Prepare the offsets for the FV solver ---*/

  if (!fem_solver) {

    /*--- Search all send/recv boundaries on this partition for any periodic
     nodes that were part of the original domain. We want to recover these
     for visualization purposes. ---*/

    unsigned long iVertex;
    bool isPeriodic;

    Local_Halo_Sort = new int[geometry->GetnPoint()];
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      Local_Halo_Sort[iPoint] = !geometry->node[iPoint]->GetDomain();

    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo_Sort[iPoint] = false;
        }
      }
    }

    /*--- Sum total number of nodes that belong to the domain ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo_Sort[iPoint] == false)
        nLocalPoint_Sort++;

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalPoint_Sort, &nGlobalPoint_Sort, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
    nGlobalPoint_Sort = nLocalPoint_Sort;
#endif

    /*--- Set a final variable for the number of points in the restart
     file. We do not write the periodic points for the FV solver, even if
     they show up in the viz. files. ---*/

    nPoint_Restart = geometry->GetGlobal_nPointDomain();

  } else {

    /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
     geometrical information for the FEM DG solver. ---*/

    CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

    unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
    CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

    /*--- Update the solution by looping over the owned volume elements. ---*/

    for(unsigned long l=0; l<nVolElemOwned; ++l) {

      /* Count up the number of local points we have for allocating storage. */

      for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
        nLocalPoint_Sort++;
      }
    }

    Local_Halo_Sort = new int[nLocalPoint_Sort];
    for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++)
      Local_Halo_Sort[iPoint] = false;

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalPoint_Sort, &nGlobalPoint_Sort, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
    nGlobalPoint_Sort = nLocalPoint_Sort;
#endif

    /*--- Set a final variable for the number of points in the restart
     file. We do not write the periodic points for the FV solver, even if
     they show up in the viz. files. ---*/

    nPoint_Restart = nGlobalPoint_Sort;

  }

  /*--- Now that we know the actual number of points we need to output,
   compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/

  beg_node = new unsigned long[size];
  end_node = new unsigned long[size];

  nPoint_Lin = new unsigned long[size];
  nPoint_Cum = new unsigned long[size+1];

  unsigned long total_points = 0;
  for (int ii = 0; ii < size; ii++) {
    nPoint_Lin[ii] = nGlobalPoint_Sort/size;
    total_points  += nPoint_Lin[ii];
  }

  /*--- Get the number of remainder points after the even division. ---*/

  unsigned long remainder = nGlobalPoint_Sort - total_points;
  for (unsigned long ii = 0; ii < remainder; ii++) {
    nPoint_Lin[ii]++;
  }

  /*--- Store the local number of nodes on each proc in the linear
   partitioning, the beginning/end index, and the linear partitioning
   within an array in cumulative storage format. ---*/

  beg_node[0] = 0;
  end_node[0] = beg_node[0] + nPoint_Lin[0];
  nPoint_Cum[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    beg_node[ii]   = end_node[ii-1];
    end_node[ii]   = beg_node[ii] + nPoint_Lin[ii];
    nPoint_Cum[ii] = nPoint_Cum[ii-1] + nPoint_Lin[ii-1];
  }
  nPoint_Cum[size] = nGlobalPoint_Sort;
  
}

void COutput::SortConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {

  /*--- Flags identifying the types of files to be written. ---*/

  bool Wrt_Vol = config->GetWrt_Vol_Sol();
  bool Wrt_Srf = config->GetWrt_Srf_Sol();

  /*--- Sort connectivity for each type of element (excluding halos). Note
   In these routines, we sort the connectivity into a linear partitioning
   across all processors based on the global index of the grid nodes. ---*/

  /*--- Sort volumetric grid connectivity. ---*/

  if (Wrt_Vol) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
      cout <<"Sorting volumetric grid connectivity." << endl;

    SortVolumetricConnectivity_FEM(config, geometry, TRIANGLE     );
    SortVolumetricConnectivity_FEM(config, geometry, QUADRILATERAL);
    SortVolumetricConnectivity_FEM(config, geometry, TETRAHEDRON  );
    SortVolumetricConnectivity_FEM(config, geometry, HEXAHEDRON   );
    SortVolumetricConnectivity_FEM(config, geometry, PRISM        );
    SortVolumetricConnectivity_FEM(config, geometry, PYRAMID      );

  }

  /*--- Sort surface grid connectivity. ---*/

  if (Wrt_Srf) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
      cout <<"Sorting surface grid connectivity." << endl;

    SortSurfaceConnectivity_FEM(config, geometry, LINE         );
    SortSurfaceConnectivity_FEM(config, geometry, TRIANGLE     );
    SortSurfaceConnectivity_FEM(config, geometry, QUADRILATERAL);

  }

  /*--- Reduce the total number of cells we will be writing in the output files. ---*/

  unsigned long nTotal_Elem = nParallel_Tria + nParallel_Quad + nParallel_Tetr + nParallel_Hexa + nParallel_Pris + nParallel_Pyra;
  unsigned long nTotal_Surf_Elem = nParallel_Line + nParallel_BoundTria + nParallel_BoundQuad;
#ifndef HAVE_MPI
  nGlobal_Elem_Par = nTotal_Elem;
  nSurf_Elem_Par   = nTotal_Surf_Elem;
#else
  SU2_MPI::Allreduce(&nTotal_Elem, &nGlobal_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nTotal_Surf_Elem, &nSurf_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
}

void COutput::SortVolumetricConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

}

void COutput::SortSurfaceConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

}

void COutput::SortOutputData_FEM(CConfig *config, CGeometry *geometry) {

  unsigned long iProcessor;
  unsigned long iPoint, Global_Index;

  /* For convenience, set the total number of variables stored at each DOF. */

  int VARS_PER_POINT = GlobalField_Counter;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/

  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

  /*--- Create the map from the global DOF ID to the local index. ---*/

  //map<unsigned long, unsigned long> mapLocal2Global;
  vector<unsigned long> globalID;

  /*--- Update the solution by looping over the owned volume elements. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Loop over the DOFs for this element and store the solution. */
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      globalID.push_back(globalIndex);
    }
  }

  /*--- We start with the grid nodes distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many nodes we must send to each other rank in order to
   have all nodes sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first ~ nGlobalPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/

  int *nPoint_Send = new int[size+1]; nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1]; nPoint_Recv[0] = 0;
  int *nPoint_Flag = new int[size];

  for (int ii=0; ii < size; ii++) {
    nPoint_Send[ii] = 0;
    nPoint_Recv[ii] = 0;
    nPoint_Flag[ii]= -1;
  }
  nPoint_Send[size] = 0; nPoint_Recv[size] = 0;

  for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++ ) {

    /*--- Get the global index of the current point. ---*/

    Global_Index = globalID[iPoint];

    /*--- Search for the processor that owns this point ---*/

    iProcessor = Global_Index/nPoint_Lin[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Cum[iProcessor])
      while(Global_Index >= nPoint_Cum[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Cum[iProcessor])   iProcessor--;

    /*--- If we have not visted this node yet, increment our
     number of elements that must be sent to a particular proc. ---*/

    if (nPoint_Flag[iProcessor] != (int)iPoint) {
      nPoint_Flag[iProcessor] = (int)iPoint;
      nPoint_Send[iProcessor+1]++;
    }

  }

  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nPoint_Recv[1] = nPoint_Send[1];
#endif

  /*--- Prepare to send coordinates. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nPoint_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nPoint_Recv[ii+1] > 0)) nRecvs++;

    nPoint_Send[ii+1] += nPoint_Send[ii];
    nPoint_Recv[ii+1] += nPoint_Recv[ii];
  }

  /*--- Allocate memory to hold the connectivity that we are sending. ---*/

  su2double *connSend = NULL;
  connSend = new su2double[VARS_PER_POINT*nPoint_Send[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Send[size]; ii++)
    connSend[ii] = 0;

  /*--- Allocate arrays for sending the global ID. ---*/

  unsigned long *idSend = new unsigned long[nPoint_Send[size]];
  for (int ii = 0; ii < nPoint_Send[size]; ii++)
    idSend[ii] = 0;

  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = VARS_PER_POINT*nPoint_Send[ii];

  unsigned long *idIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) idIndex[ii] = nPoint_Send[ii];

  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/

  for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++) {

    /*--- Get the index of the current point. ---*/

    Global_Index = globalID[iPoint];

    /*--- Search for the processor that owns this point. ---*/

    iProcessor = Global_Index/nPoint_Lin[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Cum[iProcessor])
      while(Global_Index >= nPoint_Cum[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Cum[iProcessor])   iProcessor--;

    /*--- Load data into the buffer for sending. ---*/

    if (nPoint_Flag[iProcessor] != (int)iPoint) {

      nPoint_Flag[iProcessor] = (int)iPoint;
      unsigned long nn = index[iProcessor];

      /*--- Load the data values. ---*/

      for (unsigned short kk = 0; kk < VARS_PER_POINT; kk++) {
        connSend[nn] = Local_Data[iPoint][kk]; nn++;
      }

      /*--- Load the global ID (minus offset) for sorting the
       points once they all reach the correct processor. ---*/

      nn = idIndex[iProcessor];
      idSend[nn] = Global_Index - beg_node[iProcessor];

      /*--- Increment the index by the message length ---*/

      index[iProcessor]  += VARS_PER_POINT;
      idIndex[iProcessor]++;

    }

  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;
  delete [] idIndex;

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  su2double *connRecv = NULL;
  connRecv = new su2double[VARS_PER_POINT*nPoint_Recv[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Recv[size]; ii++)
    connRecv[ii] = 0;

  unsigned long *idRecv = new unsigned long[nPoint_Recv[size]];
  for (int ii = 0; ii < nPoint_Recv[size]; ii++)
    idRecv[ii] = 0;

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the global IDs. ---*/

  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];

  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = VARS_PER_POINT*nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = VARS_PER_POINT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_DOUBLE, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = VARS_PER_POINT*nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = VARS_PER_POINT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_DOUBLE, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Repeat the process to communicate the global IDs. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the global IDs. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  int mm = VARS_PER_POINT*nPoint_Recv[rank];
  int ll = VARS_PER_POINT*nPoint_Send[rank];
  int kk = VARS_PER_POINT*nPoint_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];

  mm = nPoint_Recv[rank];
  ll = nPoint_Send[rank];
  kk = nPoint_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. First, allocate the
   appropriate amount of memory for this section. ---*/

  Parallel_Data = new su2double*[VARS_PER_POINT];
  for (int jj = 0; jj < VARS_PER_POINT; jj++) {
    Parallel_Data[jj] = new su2double[nPoint_Recv[size]];
    for (int ii = 0; ii < nPoint_Recv[size]; ii++) {
      Parallel_Data[jj][idRecv[ii]] = connRecv[ii*VARS_PER_POINT+jj];
    }
  }

  /*--- Store the total number of local points my rank has for
   the current section after completing the communications. ---*/

  nParallel_Poin = nPoint_Recv[size];

  /*--- Reduce the total number of points we will write in the output files. ---*/

#ifndef HAVE_MPI
  nGlobal_Poin_Par = nParallel_Poin;
#else
  SU2_MPI::Allreduce(&nParallel_Poin, &nGlobal_Poin_Par, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Free temporary memory from communications ---*/

  delete [] connSend;
  delete [] connRecv;
  delete [] idSend;
  delete [] idRecv;
  delete [] nPoint_Recv;
  delete [] nPoint_Send;
  delete [] nPoint_Flag;
  
  for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++)
    delete [] Local_Data[iPoint];
  delete [] Local_Data;
  
}

void COutput::SortOutputData_Surface_FEM(CConfig *config, CGeometry *geometry) {

}

void COutput::SetHistoryFile_Header(CConfig *config) { 

  
  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
      (config->GetOutput_FileFormat() == FIELDVIEW) ||
      (config->GetOutput_FileFormat() == FIELDVIEW_BINARY)) {
    HistFile << "TITLE = \"SU2 Simulation\"" << endl;
    HistFile << "VARIABLES = ";
  }
  
  stringstream out;
  string RequestedField;
  std::vector<bool> found_field(nRequestedHistoryFields, false);
    
  for (unsigned short iField_Output = 0; iField_Output < HistoryOutput_List.size(); iField_Output++){
    HistoryOutputField &Field = HistoryOutput_Map[HistoryOutput_List[iField_Output]];
    for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      RequestedField = RequestedHistoryFields[iReqField];   
      if (RequestedField == Field.OutputGroup || (RequestedField == HistoryOutput_List[iField_Output])){
        found_field[iReqField] = true;        
        out << "\"" << Field.FieldName << "\"" << HistorySep;
      }
    }  
  }
  
  for (unsigned short iField_Output = 0; iField_Output < HistoryOutputPerSurface_List.size(); iField_Output++){
    for (unsigned short iMarker = 0; iMarker < HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField_Output]].size(); iMarker++){  
      HistoryOutputField &Field = HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField_Output]][iMarker];
      for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
        RequestedField = RequestedHistoryFields[iReqField];   
        if (RequestedField == Field.OutputGroup || (RequestedField == HistoryOutputPerSurface_List[iField_Output])){
          found_field[iReqField] = true;          
          out << "\"" << Field.FieldName << "\"" << HistorySep;        
        }
      }
    }
  }
  
  
  for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
    if (!found_field[iReqField]){
      SU2_MPI::Error("Requested history field " + RequestedHistoryFields[iReqField] + " not defined in the current solver.", CURRENT_FUNCTION);
    }
  }
  
  /*--- Print the string to file and remove the last character (a separator) ---*/
  HistFile << out.str().substr(0, out.str().size() - 1);
  HistFile << endl;
  if (config->GetOutput_FileFormat() == TECPLOT ||
      config->GetOutput_FileFormat() == TECPLOT_BINARY ||
      config->GetOutput_FileFormat() == FIELDVIEW ||
      config->GetOutput_FileFormat() == FIELDVIEW_BINARY) {
    HistFile << "ZONE T= \"Convergence history\"" << endl;
  }
  HistFile.flush();
  
}


void COutput::SetHistoryFile_Output(CConfig *config) { 
  
  stringstream out;
  string RequestedField;
  
  
  for (unsigned short iField_Output = 0; iField_Output < HistoryOutput_List.size(); iField_Output++){
    HistoryOutputField &Field = HistoryOutput_Map[HistoryOutput_List[iField_Output]];
    for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
      RequestedField = RequestedHistoryFields[iReqField];   
      if (RequestedField == Field.OutputGroup){
        out << std::setprecision(10) << Field.Value << HistorySep << " ";
      }
    }
  }
  
  for (unsigned short iField_Output = 0; iField_Output < HistoryOutputPerSurface_List.size(); iField_Output++){
    for (unsigned short iMarker = 0; iMarker < HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField_Output]].size(); iMarker++){
      HistoryOutputField &Field = HistoryOutputPerSurface_Map[HistoryOutputPerSurface_List[iField_Output]][iMarker];
      for (unsigned short iReqField = 0; iReqField < nRequestedHistoryFields; iReqField++){
        RequestedField = RequestedHistoryFields[iReqField];   
        if (RequestedField == Field.OutputGroup){
          out << std::setprecision(10) << Field.Value << HistorySep << " ";
        }
      }
    }
  }
  
  /*--- Print the string to file and remove the last two characters (a separator and a space) ---*/
  
  HistFile << out.str().substr(0, out.str().size()-2);
  HistFile << endl;
  HistFile.flush();
}

void COutput::SetScreen_Header(CConfig *config) {
  if (config->GetMultizone_Problem()) 
    MultiZoneHeaderTable->PrintHeader();
  ConvergenceTable->PrintHeader();
}


void COutput::SetScreen_Output(CConfig *config) {
  
  string RequestedField;
  
  for (unsigned short iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
    stringstream out;
    RequestedField = RequestedScreenFields[iReqField]; 
    if (HistoryOutput_Map.count(RequestedField) > 0){  
      switch (HistoryOutput_Map[RequestedField].ScreenFormat) {
        case FORMAT_INTEGER:
          PrintScreenInteger(out, SU2_TYPE::Int(HistoryOutput_Map[RequestedField].Value));
          break;
        case FORMAT_FIXED:
          PrintScreenFixed(out, HistoryOutput_Map[RequestedField].Value);
          break;
        case FORMAT_SCIENTIFIC:
          PrintScreenScientific(out, HistoryOutput_Map[RequestedField].Value);
          break;      
      }
    }
    if (HistoryOutputPerSurface_Map.count(RequestedField) > 0){
      switch (HistoryOutputPerSurface_Map[RequestedField][0].ScreenFormat) {
        case FORMAT_INTEGER:
          PrintScreenInteger(out, SU2_TYPE::Int(HistoryOutputPerSurface_Map[RequestedField][0].Value));
          break;
        case FORMAT_FIXED:
          PrintScreenFixed(out, HistoryOutputPerSurface_Map[RequestedField][0].Value);
          break;
        case FORMAT_SCIENTIFIC:
          PrintScreenScientific(out, HistoryOutputPerSurface_Map[RequestedField][0].Value);
          break;   
      }
    }      
    (*ConvergenceTable) << out.str();
  }
}


void COutput::PreprocessHistoryOutput(CConfig *config){
  
  if (rank == MASTER_NODE){
   
    HistorySep = ",";
      
    char buffer[50];
    
     /*--- Retrieve the history filename ---*/
    
    string history_filename = config->GetConv_FileName();
    
     /*--- Append the zone ID ---*/
    
    if(config->GetnZone() > 1){
      history_filename = config->GetMultizone_HistoryFileName(history_filename, config->GetiZone());
    }
    strcpy (char_histfile, history_filename.data());
    
     /*--- Append the restart iteration: if dynamic problem and restart ---*/
    
    if (config->GetTime_Domain() && config->GetRestart()) {
      long iExtIter = config->GetRestart_Iter();
      if (SU2_TYPE::Int(iExtIter) < 10) SPRINTF (buffer, "_0000%d", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 10) && (SU2_TYPE::Int(iExtIter) < 100)) SPRINTF (buffer, "_000%d", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 100) && (SU2_TYPE::Int(iExtIter) < 1000)) SPRINTF (buffer, "_00%d", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d", SU2_TYPE::Int(iExtIter));
      if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d", SU2_TYPE::Int(iExtIter));
      strcat(char_histfile, buffer);
    }
    
    /*--- Add the correct file extension depending on the file format ---*/
    
    if ((config->GetOutput_FileFormat() == TECPLOT) ||
        (config->GetOutput_FileFormat() == FIELDVIEW)) SPRINTF (buffer, ".dat");
    else if ((config->GetOutput_FileFormat() == TECPLOT_BINARY) ||
             (config->GetOutput_FileFormat() == FIELDVIEW_BINARY))  SPRINTF (buffer, ".plt");
    else if (config->GetOutput_FileFormat() == PARAVIEW)  SPRINTF (buffer, ".csv");
    strcat(char_histfile, buffer);
    
//    if(!(std::find(RequestedHistoryFields.begin(), RequestedHistoryFields.end(), "EXT_ITER") != RequestedHistoryFields.end())) {
//      RequestedHistoryFields.push_back("EXT_ITER");
//      nRequestedHistoryFields++;
//    }
    
    /*--- Set the History output fields using a virtual function call to the child implementation ---*/
    
    SetHistoryOutputFields(config);
   
    Postprocess_HistoryFields(config);
    
    /*--- Open the history file ---*/
    
    cout << "History filename: " << char_histfile << endl;
    HistFile.open(char_histfile, ios::out);
    HistFile.precision(15);
    
    /*--- Add the header to the history file. ---*/
    
    SetHistoryFile_Header(config);    
    
    string RequestedField;
  
    /*--- Set screen convergence output header ---*/
    
    // Evaluate the requested output
    for (unsigned short iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
      RequestedField = RequestedScreenFields[iReqField];  
      if (HistoryOutput_Map.count(RequestedField) > 0){ 
        ConvergenceTable->AddColumn(HistoryOutput_Map[RequestedField].FieldName, field_width);
      } else {
  //      SU2_MPI::Error(string("Requested screen output field not found: ") + currentField, CURRENT_FUNCTION);
      }
      if (HistoryOutputPerSurface_Map.count(RequestedField) > 0){
        ConvergenceTable->AddColumn(HistoryOutputPerSurface_Map[RequestedField][0].FieldName, field_width);
      }
    }
    
    if (config->GetMultizone_Problem()){
      MultiZoneHeaderTable->AddColumn(MultiZoneHeaderString, nRequestedScreenFields*field_width + (nRequestedScreenFields-1));      
      MultiZoneHeaderTable->SetAlign(PrintingToolbox::CTablePrinter::CENTER);
      MultiZoneHeaderTable->SetPrintHeaderBottomLine(false);
    }
    
  }
  
}

void COutput::PreprocessVolumeOutput(CConfig *config, CGeometry *geometry){
  
  /*--- Make sure that coordinates are always in the volume output --- */
  
  if(!(std::find(RequestedVolumeFields.begin(), RequestedVolumeFields.end(), "COORDINATES") != RequestedVolumeFields.end())) {
    RequestedVolumeFields.push_back("COORDINATES");
    nRequestedVolumeFields++;
  }
  
  /*--- Set the volume output fields using a virtual function call to the child implementation ---*/  
  
  SetVolumeOutputFields(config);
 
  // TODO: Add check for coordinates in solver output 
  
  GlobalField_Counter = 0;
  
  string RequestedField;
  std::vector<bool> found_field(nRequestedVolumeFields, false);
  
  
  /*--- Loop through all fields defined in the corresponding SetVolumeOutputFields(). 
 * If it is also defined in the config (either as part of a group or a single field), the field 
 * object gets an offset so that we know where to find the data in the Local_Data() array.
 *  Note that the default offset is -1. An index !=-1 defines this field as part of the output. ---*/

  for (unsigned short iField_Output = 0; iField_Output < VolumeOutput_List.size(); iField_Output++){
    
    VolumeOutputField &Field = VolumeOutput_Map[VolumeOutput_List[iField_Output]];
    
    /*--- Loop through all fields specified in the config ---*/
    
    for (unsigned short iReqField = 0; iReqField < nRequestedVolumeFields; iReqField++){
      
      RequestedField = RequestedVolumeFields[iReqField];  
            
      if (((RequestedField == Field.OutputGroup) || (RequestedField == VolumeOutput_List[iField_Output])) && (Field.Offset == -1)){
        Field.Offset = GlobalField_Counter;
        Variable_Names.push_back(Field.FieldName);
        GlobalField_Counter++;
        
        found_field[iReqField] = true;
      }
    }
    
    
//    if ((std::find(found_field.begin(), found_field.end(), false) != found_field.end())){
////      SU2_MPI::Error(string("There is no volume output field/group with name ") + std::find(found_field.begin(), found_field.end(), false) + string(" defined in the current solver."), CURRENT_FUNCTION);
//    }
  }
  
  
  /*--- Now that we know the number of fields, create the local data array to temporarily store the volume output 
   * before writing it to file ---*/
    
  unsigned long iPoint, iVertex;
  bool Wrt_Halo, isPeriodic;
  
  unsigned short iMarker, iField;
  
  Wrt_Halo = config->GetWrt_Halo();
  
  Local_Data = new su2double*[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    Local_Data[iPoint] = new su2double[GlobalField_Counter];
    for (iField = 0; iField < GlobalField_Counter; iField++){
      Local_Data[iPoint][iField] = 0.0;
    }
  }
  
  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();
  
  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/
  
  if (!Wrt_Halo) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        
        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/
        
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
  }
}

void COutput::CollectVolumeData(CConfig* config, CGeometry* geometry, CSolver** solver){
  
  bool Wrt_Halo = config->GetWrt_Halo();
  unsigned short iMarker = 0;
  unsigned long iPoint = 0, jPoint = 0;
  long iVertex = 0;
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check for halos & write only if requested ---*/
    
    if (!Local_Halo[iPoint] || Wrt_Halo) {
      
      /*--- Load the volume data into the Local_Data() array. --- */
      
      LoadVolumeData(config, geometry, solver, jPoint);
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_Plotting(iMarker) == YES) {
          iVertex = geometry->node[iPoint]->GetVertex(iMarker);
          if (iVertex != -1){ 
            
            /*--- Load the surface data into the Local_Data() array. --- */
            
            LoadSurfaceData(config, geometry, solver, jPoint, iMarker, iVertex);
          }
        }
      }
      jPoint++;
    }
  }
}

void COutput::Postprocess_HistoryData(CConfig *config, bool dualtime){
  
  for (unsigned short iField = 0; iField < HistoryOutput_List.size(); iField++){
    HistoryOutputField &currentField = HistoryOutput_Map[HistoryOutput_List[iField]];
    if (currentField.FieldType == TYPE_RESIDUAL){
      if ( SetInit_Residuals(config) ) {
        Init_Residuals[HistoryOutput_List[iField]] = currentField.Value;
      }
      SetHistoryOutputValue("REL_" + HistoryOutput_List[iField], currentField.Value - Init_Residuals[HistoryOutput_List[iField]]);
    }
    if (currentField.FieldType == TYPE_COEFFICIENT){
      if(SetUpdate_Averages(config, dualtime)){
        SetHistoryOutputValue("TAVG_" + HistoryOutput_List[iField], RunningAverages[HistoryOutput_List[iField]].Update(currentField.Value));
      }
      if (config->GetDirectDiff() != NO_DERIVATIVE){
        SetHistoryOutputValue("D_" + HistoryOutput_List[iField], SU2_TYPE::GetDerivative(currentField.Value));      
        SetHistoryOutputValue("D_TAVG_" + HistoryOutput_List[iField], SU2_TYPE::GetDerivative(RunningAverages[HistoryOutput_List[iField]].Get()));
        
      }
    }
  }
}

void COutput::Postprocess_HistoryFields(CConfig *config){
  for (unsigned short iField = 0; iField < HistoryOutput_List.size(); iField++){
    HistoryOutputField &currentField = HistoryOutput_Map[HistoryOutput_List[iField]];
    if (currentField.FieldType == TYPE_RESIDUAL){
      AddHistoryOutput("REL_" + HistoryOutput_List[iField], "rel" + currentField.FieldName, currentField.ScreenFormat, "REL_" + currentField.OutputGroup);
    }
    if (currentField.FieldType == TYPE_COEFFICIENT){
      AddHistoryOutput("TAVG_"   + HistoryOutput_List[iField], "tavg["  + currentField.FieldName + "]", currentField.ScreenFormat, "TAVG_"   + currentField.OutputGroup);
      AddHistoryOutput("D_"      + HistoryOutput_List[iField], "d["     + currentField.FieldName + "]", currentField.ScreenFormat, "D_"      + currentField.OutputGroup);  
      AddHistoryOutput("D_TAVG_" + HistoryOutput_List[iField], "dtavg[" + currentField.FieldName + "]", currentField.ScreenFormat, "D_TAVG_" + currentField.OutputGroup);  
    }
  }
}
