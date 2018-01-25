/*!
 * \file transfer_structure.cpp
 * \brief Main subroutines for physics of the information transfer between zones
 * \author R. Sanchez
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/transfer_structure.hpp"

CTransfer_FlowTraction::CTransfer_FlowTraction(void) : CTransfer() {

}

CTransfer_FlowTraction::CTransfer_FlowTraction(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_FlowTraction::~CTransfer_FlowTraction(void) {

}

void CTransfer_FlowTraction::GetPhysical_Constants(CSolver *flow_solution, CSolver *struct_solution,
                                                         CGeometry *flow_geometry, CGeometry *struct_geometry,
                           CConfig *flow_config, CConfig *struct_config) {

  unsigned short iVar;

  /*--- We have to clear the traction before applying it, because we are "adding" to node and not "setting" ---*/

  for (unsigned long iPoint = 0; iPoint < struct_geometry->GetnPoint(); iPoint++) 
    struct_solution->node[iPoint]->Clear_FlowTraction();

  /*--- Redimensionalize the pressure ---*/

  su2double *Velocity_ND, *Velocity_Real;
  su2double Density_ND,  Density_Real, Velocity2_Real, Velocity2_ND;

  Velocity_Real = flow_config->GetVelocity_FreeStream();
  Density_Real  = flow_config->GetDensity_FreeStream();

  Velocity_ND = flow_config->GetVelocity_FreeStreamND();
  Density_ND  = flow_config->GetDensity_FreeStreamND();

  Velocity2_Real = 0.0;
  Velocity2_ND   = 0.0;
  for (iVar = 0; iVar < nVar; iVar++) {
    Velocity2_Real += Velocity_Real[iVar]*Velocity_Real[iVar];
    Velocity2_ND   += Velocity_ND[iVar]*Velocity_ND[iVar];
  }

  Physical_Constants[0] = Density_Real * Velocity2_Real / ( Density_ND * Velocity2_ND );

  /*--- Apply a ramp to the transfer of the fluid loads ---*/

  su2double ModAmpl = 0.0;
  su2double CurrentTime = struct_config->GetCurrent_DynTime();
  su2double Static_Time = struct_config->GetStatic_Time();
  su2double Transfer_Time = 0.0;

  bool Ramp_Load = struct_config->GetRamp_Load();
  su2double Ramp_Time = struct_config->GetRamp_Time();

  /*--- Polynomial functions from https://en.wikipedia.org/wiki/Smoothstep ---*/
  if (CurrentTime < Static_Time){
    ModAmpl = 0.0;
  }
  else if (Ramp_Load){
    Transfer_Time = (CurrentTime - Static_Time) / Ramp_Time;
    switch (struct_config->GetDynamic_LoadTransfer()){
    case INSTANTANEOUS:
      ModAmpl = 1.0;
      break;
    case POL_ORDER_1:
      ModAmpl = Transfer_Time;
      break;
    case POL_ORDER_3:
      ModAmpl = -2.0 * pow(Transfer_Time,3.0) + 3.0 * pow(Transfer_Time,2.0);
      break;
    case POL_ORDER_5:
      ModAmpl = 6.0 * pow(Transfer_Time, 5.0) - 15.0 * pow(Transfer_Time, 4.0) + 10 * pow(Transfer_Time, 3.0);
      break;
    case SIGMOID_10:
      ModAmpl = (1 / (1+exp(-1.0 * 10.0 * (Transfer_Time - 0.5)) ) );
      break;
    case SIGMOID_20:
      ModAmpl = (1 / (1+exp(-1.0 * 20.0 * (Transfer_Time - 0.5)) ) );
      break;
    }
    ModAmpl = max(ModAmpl,0.0);
    ModAmpl = min(ModAmpl,1.0);
  }
  else{
    ModAmpl = 1.0;
  }
  
  Physical_Constants[1] = ModAmpl;
  
  /*--- For static FSI, we cannot apply the ramp like this ---*/
  if ((flow_config->GetUnsteady_Simulation() == STEADY) && (struct_config->GetDynamic_Analysis() == STATIC)){
    Physical_Constants[1] = 1.0;
    if (Ramp_Load){
      CurrentTime = static_cast<su2double>(struct_config->GetFSIIter());
      Ramp_Time = static_cast<su2double>(struct_config->GetnIterFSI_Ramp() - 1);
      if (Ramp_Time != 0.0) Transfer_Time = CurrentTime / Ramp_Time;
      switch (struct_config->GetDynamic_LoadTransfer()){
      case INSTANTANEOUS:
        ModAmpl = 1.0;
        break;
      case POL_ORDER_1:
        ModAmpl = Transfer_Time;
        break;
      case POL_ORDER_3:
        ModAmpl = -2.0 * pow(Transfer_Time,3.0) + 3.0 * pow(Transfer_Time,2.0);
        break;
      case POL_ORDER_5:
        ModAmpl = 6.0 * pow(Transfer_Time, 5.0) - 15.0 * pow(Transfer_Time, 4.0) + 10 * pow(Transfer_Time, 3.0);
        break;
      case SIGMOID_10:
        ModAmpl = (1 / (1+exp(-1.0 * 10.0 * (Transfer_Time - 0.5)) ) );
        break;
      case SIGMOID_20:
        ModAmpl = (1 / (1+exp(-1.0 * 20.0 * (Transfer_Time - 0.5)) ) );
        break;
      }
      ModAmpl = max(ModAmpl,0.0);
      ModAmpl = min(ModAmpl,1.0);
      if (CurrentTime > Ramp_Time) ModAmpl = 1.0;
      Physical_Constants[1] = ModAmpl;
    }
  }

}

void CTransfer_FlowTraction::GetDonor_Variable(CSolver *flow_solution, CGeometry *flow_geometry, CConfig *flow_config,
                                           unsigned long Marker_Flow, unsigned long Vertex_Flow, unsigned long Point_Struct) {


  unsigned short iVar, jVar;
  unsigned long Point_Flow;
  su2double *Normal_Flow;

  // Check the kind of fluid problem
  bool compressible       = (flow_config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (flow_config->GetKind_Regime() == INCOMPRESSIBLE);
  bool viscous_flow       = ((flow_config->GetKind_Solver() == NAVIER_STOKES) ||
                              (flow_config->GetKind_Solver() == RANS) ||
                              (flow_config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES) ||
                              (flow_config->GetKind_Solver() == DISC_ADJ_RANS));

  // Parameters for the calculations
  // Pn: Pressure
  // Pinf: Pressure_infinite
  // div_vel: Velocity divergence
  // Dij: Dirac delta
  su2double Pn = 0.0, div_vel = 0.0;
  su2double Viscosity = 0.0;
  su2double Tau[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  su2double Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  su2double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

  su2double Pinf = flow_solution->GetPressure_Inf();

  Point_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNode();
  // Get the normal at the vertex: this normal goes inside the fluid domain.
  Normal_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNormal();

  // Retrieve the values of pressure

  Pn = flow_solution->node[Point_Flow]->GetPressure();

  // Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
  for (iVar = 0; iVar < nVar; iVar++) 
    Donor_Variable[iVar] = -(Pn-Pinf)*Normal_Flow[iVar];

  // Calculate tn in the fluid nodes for the viscous term

  if ((incompressible || compressible) && viscous_flow) {

    Viscosity = flow_solution->node[Point_Flow]->GetLaminarViscosity();

    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0 ; jVar < nVar; jVar++) {
        Grad_Vel[iVar][jVar] = flow_solution->node[Point_Flow]->GetGradient_Primitive(iVar+1, jVar);
      }
    }

    // Divergence of the velocity
    div_vel = 0.0; for (iVar = 0; iVar < nVar; iVar++) div_vel += Grad_Vel[iVar][iVar];
  
    for (iVar = 0; iVar < nVar; iVar++) {
      
      for (jVar = 0 ; jVar < nVar; jVar++) {
      
        // Viscous stress
        Tau[iVar][jVar] = Viscosity*(Grad_Vel[jVar][iVar] + Grad_Vel[iVar][jVar]) - TWO3*Viscosity*div_vel*delta[iVar][jVar];

        // Viscous component in the tn vector --> Units of force (non-dimensional).
        Donor_Variable[iVar] += Tau[iVar][jVar]*Normal_Flow[jVar];
      }
    }
  }

  // Redimensionalize and take into account ramp transfer of the loads
  for (iVar = 0; iVar < nVar; iVar++) {
    Donor_Variable[iVar] = Donor_Variable[iVar] * Physical_Constants[0] * Physical_Constants[1];
  }

}

void CTransfer_FlowTraction::SetTarget_Variable(CSolver *fea_solution, CGeometry *fea_geometry,
                        CConfig *fea_config, unsigned long Marker_Struct,
                        unsigned long Vertex_Struct, unsigned long Point_Struct) {

  /*--- Add to the Flow traction ---*/
  fea_solution->node[Point_Struct]->Add_FlowTraction(Target_Variable);

}



CTransfer_StructuralDisplacements::CTransfer_StructuralDisplacements(void) : CTransfer() {

}

CTransfer_StructuralDisplacements::CTransfer_StructuralDisplacements(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_StructuralDisplacements::~CTransfer_StructuralDisplacements(void) {

}


void CTransfer_StructuralDisplacements::GetPhysical_Constants(CSolver *struct_solution, CSolver *flow_solution,
                                                         CGeometry *struct_geometry, CGeometry *flow_geometry,
                           CConfig *struct_config, CConfig *flow_config) {

}

void CTransfer_StructuralDisplacements::GetDonor_Variable(CSolver *struct_solution, CGeometry *struct_geometry, CConfig *struct_config,
                                                       unsigned long Marker_Struct, unsigned long Vertex_Struct, unsigned long Point_Struct) {

  su2double *DisplacementDonor, *DisplacementDonor_Prev;
  unsigned short iVar;

  /*--- The displacements come from the predicted solution ---*/
  DisplacementDonor = struct_solution->node[Point_Struct]->GetSolution_Pred();

  DisplacementDonor_Prev = struct_solution->node[Point_Struct]->GetSolution_Pred_Old();

  for (iVar = 0; iVar < nVar; iVar++) 
  Donor_Variable[iVar] = DisplacementDonor[iVar] - DisplacementDonor_Prev[iVar];
}

void CTransfer_StructuralDisplacements::SetTarget_Variable(CSolver *flow_solution, CGeometry *flow_geometry,
                               CConfig *flow_config, unsigned long Marker_Flow,
                               unsigned long Vertex_Flow, unsigned long Point_Flow) {

  su2double VarCoord[3] = {0.0, 0.0, 0.0};
  unsigned short iVar;

  for (iVar = 0; iVar < nVar; iVar++)
    VarCoord[iVar] = Target_Variable[iVar];

  flow_geometry->vertex[Marker_Flow][Vertex_Flow]->SetVarCoord(VarCoord);

}

CTransfer_StructuralDisplacements_DiscAdj::CTransfer_StructuralDisplacements_DiscAdj(void) : CTransfer() {

}

CTransfer_StructuralDisplacements_DiscAdj::CTransfer_StructuralDisplacements_DiscAdj(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_StructuralDisplacements_DiscAdj::~CTransfer_StructuralDisplacements_DiscAdj(void) {

}


void CTransfer_StructuralDisplacements_DiscAdj::GetPhysical_Constants(CSolver *struct_solution, CSolver *flow_solution,
                                                         CGeometry *struct_geometry, CGeometry *flow_geometry,
                           CConfig *struct_config, CConfig *flow_config) {

}

void CTransfer_StructuralDisplacements_DiscAdj::GetDonor_Variable(CSolver *struct_solution, CGeometry *struct_geometry, CConfig *struct_config,
                                                       unsigned long Marker_Struct, unsigned long Vertex_Struct, unsigned long Point_Struct) {


  su2double *Coord_Struct, *Displacement_Struct;
  unsigned short iVar;

  Coord_Struct = struct_geometry->node[Point_Struct]->GetCoord();

  /*--- The displacements come from the predicted solution ---*/
  Displacement_Struct = struct_solution->node[Point_Struct]->GetSolution();

  for (iVar = 0; iVar < nVar; iVar++) 
    Donor_Variable[iVar] = Coord_Struct[iVar] + Displacement_Struct[iVar];
}

void CTransfer_StructuralDisplacements_DiscAdj::SetTarget_Variable(CSolver *flow_solution, CGeometry *flow_geometry,
                               CConfig *flow_config, unsigned long Marker_Flow,
                               unsigned long Vertex_Flow, unsigned long Point_Flow) {

  su2double *Coord, VarCoord[3] = {0.0, 0.0, 0.0};
  unsigned short iVar;

  Coord = flow_geometry->node[Point_Flow]->GetCoord();

  for (iVar = 0; iVar < nVar; iVar++)
    VarCoord[iVar] = Target_Variable[iVar]-Coord[iVar];

  flow_geometry->vertex[Marker_Flow][Vertex_Flow]->SetVarCoord(VarCoord);
}

CTransfer_FlowTraction_DiscAdj::CTransfer_FlowTraction_DiscAdj(void) : CTransfer() {

}

CTransfer_FlowTraction_DiscAdj::CTransfer_FlowTraction_DiscAdj(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_FlowTraction_DiscAdj::~CTransfer_FlowTraction_DiscAdj(void) {

}

void CTransfer_FlowTraction_DiscAdj::GetPhysical_Constants(CSolver *flow_solution, CSolver *struct_solution,
                                               CGeometry *flow_geometry, CGeometry *struct_geometry,
                           CConfig *flow_config, CConfig *struct_config){

  unsigned short iVar;

  /*--- We have to clear the traction before applying it, because we are "adding" to node and not "setting" ---*/

  for (unsigned long iPoint = 0; iPoint < struct_geometry->GetnPoint(); iPoint++){
    struct_solution->node[iPoint]->Clear_FlowTraction();
  }

  /*--- Redimensionalize the pressure ---*/

  su2double *Velocity_ND, *Velocity_Real;
  su2double Density_ND,  Density_Real, Velocity2_Real, Velocity2_ND;

  Velocity_Real = flow_config->GetVelocity_FreeStream();
  Density_Real = flow_config->GetDensity_FreeStream();

  Velocity_ND = flow_config->GetVelocity_FreeStreamND();
  Density_ND = flow_config->GetDensity_FreeStreamND();

  Velocity2_Real = 0.0;
  Velocity2_ND = 0.0;
  for (iVar = 0; iVar < nVar; iVar++){
    Velocity2_Real += Velocity_Real[iVar]*Velocity_Real[iVar];
    Velocity2_ND += Velocity_ND[iVar]*Velocity_ND[iVar];
  }

  Physical_Constants[0] = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);

}

void CTransfer_FlowTraction_DiscAdj::GetDonor_Variable(CSolver *flow_solution, CGeometry *flow_geometry, CConfig *flow_config,
                                     unsigned long Marker_Flow, unsigned long Vertex_Flow, unsigned long Point_Struct){


  unsigned short iVar, jVar;
  unsigned long Point_Flow;
  su2double *Normal_Flow;

  // Check the kind of fluid problem
  bool compressible       = (flow_config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (flow_config->GetKind_Regime() == INCOMPRESSIBLE);
  bool viscous_flow       = ((flow_config->GetKind_Solver() == NAVIER_STOKES) ||
                              (flow_config->GetKind_Solver() == RANS) ||
                              (flow_config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES) ||
                              (flow_config->GetKind_Solver() == DISC_ADJ_RANS));

  // Parameters for the calculations
  // Pn: Pressure
  // Pinf: Pressure_infinite
  // div_vel: Velocity divergence
  // Dij: Dirac delta
  su2double Pn = 0.0, div_vel = 0.0;
  su2double Viscosity = 0.0;
  su2double Tau[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  su2double Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  su2double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

  su2double Pinf = flow_solution->GetPressure_Inf();

  Point_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNode();
    // Get the normal at the vertex: this normal goes inside the fluid domain.
  Normal_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNormal();

  // Retrieve the values of pressure
  Pn = flow_solution->node[Point_Flow]->GetPressure();

  // Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
  for (iVar = 0; iVar < nVar; iVar++) {
    Donor_Variable[iVar] = -(Pn-Pinf)*Normal_Flow[iVar];
  }

  // Calculate tn in the fluid nodes for the viscous term

  if ((incompressible || compressible) && viscous_flow){

    Viscosity = flow_solution->node[Point_Flow]->GetLaminarViscosity();

    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0 ; jVar < nVar; jVar++) {
        Grad_Vel[iVar][jVar] = flow_solution->node[Point_Flow]->GetGradient_Primitive(iVar+1, jVar);
      }
    }

    // Divergence of the velocity
    div_vel = 0.0; for (iVar = 0; iVar < nVar; iVar++) div_vel += Grad_Vel[iVar][iVar];

    for (iVar = 0; iVar < nVar; iVar++) {

      for (jVar = 0 ; jVar < nVar; jVar++) {

        // Viscous stress
        Tau[iVar][jVar] = Viscosity*(Grad_Vel[jVar][iVar] + Grad_Vel[iVar][jVar]) - TWO3*Viscosity*div_vel*delta[iVar][jVar];

        // Viscous component in the tn vector --> Units of force (non-dimensional).
        Donor_Variable[iVar] += Tau[iVar][jVar]*Normal_Flow[jVar];
      }
    }
  }

  // Redimensionalize
  for (iVar = 0; iVar < nVar; iVar++){
    Donor_Variable[iVar] = Donor_Variable[iVar] * Physical_Constants[0];
  }

}

void CTransfer_FlowTraction_DiscAdj::SetTarget_Variable(CSolver *fea_solution, CGeometry *fea_geometry,
                        CConfig *fea_config, unsigned long Marker_Struct,
                        unsigned long Vertex_Struct, unsigned long Point_Struct){

  /*--- Add to the Flow traction ---*/
  fea_solution->node[Point_Struct]->Add_FlowTraction(Target_Variable);

}



CTransfer_ConservativeVars::CTransfer_ConservativeVars(void) : CTransfer() {

}

CTransfer_ConservativeVars::CTransfer_ConservativeVars(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_ConservativeVars::~CTransfer_ConservativeVars(void) {

}


void CTransfer_ConservativeVars::GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                                                             CGeometry *donor_geometry, CGeometry *target_geometry,
                             CConfig *donor_config, CConfig *target_config) {

}

void CTransfer_ConservativeVars::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry, CConfig *donor_config,
                                                unsigned long Marker_Donor, unsigned long Vertex_Donor, unsigned long Point_Donor) {

  su2double *Solution;
  unsigned short iVar;

  /*--- Retrieve solution and set it as the donor variable ---*/
  Solution = donor_solution->node[Point_Donor]->GetSolution();

  for (iVar = 0; iVar < nVar; iVar++)
    Donor_Variable[iVar] = Solution[iVar];
}

void CTransfer_ConservativeVars::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                          CConfig *target_config, unsigned long Marker_Target,
                          unsigned long Vertex_Target, unsigned long Point_Target) {

  /*--- Set the target solution with the value of the Target Variable ---*/
  target_solution->node[Point_Target]->SetSolution(Target_Variable);

}

CTransfer_MixingPlaneInterface::CTransfer_MixingPlaneInterface(void) : CTransfer() {

}

CTransfer_MixingPlaneInterface::CTransfer_MixingPlaneInterface(unsigned short val_nVar, unsigned short val_nConst, CConfig *donor_config, CConfig *target_config){
  unsigned short iVar;
  nVar = val_nVar;


  Donor_Variable     = new su2double[nVar + 5];
  Target_Variable    = new su2double[nVar + 5];



  for (iVar = 0; iVar < nVar + 5; iVar++){
    Donor_Variable[iVar]  = 0.0;
    Target_Variable[iVar] = 0.0;
  }
}

CTransfer_MixingPlaneInterface::~CTransfer_MixingPlaneInterface(void) {
}



void CTransfer_MixingPlaneInterface::SetSpanWiseLevels(CConfig *donor_config, CConfig *target_config){

  unsigned short iSpan;
  nSpanMaxAllZones = donor_config->GetnSpanMaxAllZones();


  SpanValueCoeffTarget = new su2double[target_config->GetnSpanWiseSections() + 1];
  SpanLevelDonor       = new unsigned short[target_config->GetnSpanWiseSections() + 1];


  for (iSpan = 0; iSpan < target_config->GetnSpanWiseSections() + 1;iSpan++){
    SpanValueCoeffTarget[iSpan] = 0.0;
    SpanLevelDonor[iSpan]       = 1;
  }

}

void CTransfer_MixingPlaneInterface::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry,
    CConfig *donor_config, unsigned long Marker_Donor,
    unsigned long iSpan, unsigned long rank) {

  unsigned short nDim = nVar - 2;
  bool turbulent = ((donor_config->GetKind_Solver() == RANS) || (donor_config->GetKind_Solver() == DISC_ADJ_RANS));



  Donor_Variable[0] = donor_solution->GetAverageDensity(Marker_Donor, iSpan);
  Donor_Variable[1]	= donor_solution->GetAveragePressure(Marker_Donor, iSpan);
  Donor_Variable[2] = donor_solution->GetAverageTurboVelocity(Marker_Donor, iSpan)[0];
  Donor_Variable[3] = donor_solution->GetAverageTurboVelocity(Marker_Donor, iSpan)[1];

  if(nDim == 3){
    Donor_Variable[4] = donor_solution->GetAverageTurboVelocity(Marker_Donor, iSpan)[2];
  }
  else{
    Donor_Variable[4] = -1.0;
  }

  if(turbulent){
    Donor_Variable[5] = donor_solution->GetAverageNu(Marker_Donor, iSpan);
    Donor_Variable[6] = donor_solution->GetAverageKine(Marker_Donor, iSpan);
    Donor_Variable[7] = donor_solution->GetAverageOmega(Marker_Donor, iSpan);
  }
  else{
    Donor_Variable[5] = -1.0;
    Donor_Variable[6] = -1.0;
    Donor_Variable[7] = -1.0;
  }

}


void CTransfer_MixingPlaneInterface::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
    CConfig *target_config, unsigned long Marker_Target,
    unsigned long iSpan, unsigned long rank) {

  unsigned short nDim = nVar - 2;
  bool turbulent = ((target_config->GetKind_Solver() == RANS) || (target_config->GetKind_Solver() == DISC_ADJ_RANS));


  target_solution->SetExtAverageDensity(Marker_Target, iSpan, Target_Variable[0]);
  target_solution->SetExtAveragePressure(Marker_Target, iSpan, Target_Variable[1]);
  target_solution->SetExtAverageTurboVelocity(Marker_Target, iSpan, 0, Target_Variable[2]);
  target_solution->SetExtAverageTurboVelocity(Marker_Target, iSpan, 1, Target_Variable[3]);

  if(nDim == 3){
    target_solution->SetExtAverageTurboVelocity(Marker_Target, iSpan, 2, Target_Variable[4]);
  }

  if(turbulent){
    target_solution->SetExtAverageNu(Marker_Target, iSpan, Target_Variable[5]);
    target_solution->SetExtAverageKine(Marker_Target, iSpan, Target_Variable[6]);
    target_solution->SetExtAverageOmega(Marker_Target, iSpan,  Target_Variable[7]);
  }

}

void CTransfer_MixingPlaneInterface::SetAverageValues(CSolver *donor_solution, CSolver *target_solution, unsigned short donorZone){
  unsigned short iSpan;

  for(iSpan = 0; iSpan<nSpanMaxAllZones +1; iSpan++){
    /*--- trnasfer inviscid quantities ---*/
    target_solution->SetDensityIn(donor_solution->GetDensityIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetPressureIn(donor_solution->GetPressureIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetTurboVelocityIn(donor_solution->GetTurboVelocityIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetDensityOut(donor_solution->GetDensityOut(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetPressureOut(donor_solution->GetPressureOut(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetTurboVelocityOut(donor_solution->GetTurboVelocityOut(donorZone, iSpan), donorZone, iSpan);

    /*--- transfer turbulent quantities ---*/
    target_solution->SetKineIn(donor_solution->GetKineIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetOmegaIn(donor_solution->GetOmegaIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetNuIn(donor_solution->GetNuIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetKineOut(donor_solution->GetKineOut(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetOmegaOut(donor_solution->GetOmegaOut(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetNuOut(donor_solution->GetNuOut(donorZone, iSpan), donorZone, iSpan);

  }
}

void CTransfer_MixingPlaneInterface::SetAverageTurboGeoValues(CGeometry *donor_geometry, CGeometry *target_geometry, unsigned short donorZone){
  unsigned short iSpan;

  for(iSpan = 0; iSpan<nSpanMaxAllZones+1; iSpan++){
    target_geometry->SetTurboRadiusIn(donor_geometry->GetTurboRadiusIn(donorZone, iSpan), donorZone, iSpan);
    target_geometry->SetSpanAreaIn(donor_geometry->GetSpanAreaIn(donorZone, iSpan), donorZone, iSpan);
    target_geometry->SetTangGridVelIn(donor_geometry->GetTangGridVelIn(donorZone, iSpan), donorZone, iSpan);
    target_geometry->SetTurboRadiusOut(donor_geometry->GetTurboRadiusOut(donorZone, iSpan), donorZone, iSpan);
    target_geometry->SetSpanAreaOut(donor_geometry->GetSpanAreaOut(donorZone, iSpan), donorZone, iSpan);
    target_geometry->SetTangGridVelOut(donor_geometry->GetTangGridVelOut(donorZone, iSpan), donorZone, iSpan);
  }

}



CTransfer_SlidingInterface::CTransfer_SlidingInterface(void) : CTransfer() {

}

CTransfer_SlidingInterface::CTransfer_SlidingInterface(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_SlidingInterface::~CTransfer_SlidingInterface(void) {

}


void CTransfer_SlidingInterface::GetPhysical_Constants(CSolver *donor_solution, CSolver *target_solution,
                                                             CGeometry *donor_geometry, CGeometry *target_geometry,
                             CConfig *donor_config, CConfig *target_config) {

}

void CTransfer_SlidingInterface::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry, CConfig *donor_config,
                                                unsigned long Marker_Donor, unsigned long Vertex_Donor, unsigned long Point_Donor) {

  unsigned short iVar, nDonorVar;
  nDonorVar = donor_solution->GetnPrimVar();

  /*---  the number of primitive variables is set to two by default for the turbulent solver ---*/
  bool turbulent = (nDonorVar == 2) ;

  if (turbulent){

    /*---  for turbulent solver retrieve solution and set it as the donor variable ---*/
    Donor_Variable[0] = donor_solution->node[Point_Donor]->GetSolution(0);
    Donor_Variable[1] = donor_solution->node[Point_Donor]->GetSolution(1);

  } else{

    /*---  Retrieve primitive variables and set them as the donor variables ---*/
    for (iVar = 0; iVar < nDonorVar; iVar++)
    Donor_Variable[iVar] = donor_solution->node[Point_Donor]->GetPrimitive(iVar);

  }
}

void CTransfer_SlidingInterface::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                          CConfig *target_config, unsigned long Marker_Target,
                          unsigned long Vertex_Target, unsigned long Point_Target) {

  unsigned short iVar, iDonorVertex, nTargetVar;
  nTargetVar = target_solution->GetnPrimVar();
  /*--- Set the Sliding solution with the value of the Target Variable ---*/

  iDonorVertex = target_solution->GetnSlidingStates(Marker_Target, Vertex_Target);

  for (iVar = 0; iVar < nTargetVar+1; iVar++)
    target_solution->SetSlidingState(Marker_Target, Vertex_Target, iVar, iDonorVertex, Target_Variable[iVar]);

  target_solution->SetnSlidingStates( Marker_Target, Vertex_Target, iDonorVertex + 1 );
}
