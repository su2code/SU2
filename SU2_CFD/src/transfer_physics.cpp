/*!
 * \file transfer_structure.cpp
 * \brief Main subroutines for physics of the information transfer between zones
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

#include "../include/transfer_structure.hpp"

CTransfer_FlowTraction::CTransfer_FlowTraction(void) : CTransfer() {

}

CTransfer_FlowTraction::CTransfer_FlowTraction(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_FlowTraction::~CTransfer_FlowTraction(void) {

}

void CTransfer_FlowTraction::Preprocess(CConfig *flow_config) {

  /*--- Store if consistent interpolation is in use, in which case we need to transfer stresses
        and integrate on the structural side rather than directly transferring forces. ---*/
  consistent_interpolation = (!flow_config->GetConservativeInterpolation() ||
                             (flow_config->GetKindInterpolation() == WEIGHTED_AVERAGE));

  /*--- Compute the constant factor to dimensionalize pressure and shear stress. ---*/
  su2double *Velocity_ND, *Velocity_Real;
  su2double Density_ND,  Density_Real, Velocity2_Real, Velocity2_ND;

  Velocity_Real = flow_config->GetVelocity_FreeStream();
  Density_Real  = flow_config->GetDensity_FreeStream();

  Velocity_ND = flow_config->GetVelocity_FreeStreamND();
  Density_ND  = flow_config->GetDensity_FreeStreamND();

  Velocity2_Real = 0.0;
  Velocity2_ND   = 0.0;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Velocity2_Real += Velocity_Real[iVar]*Velocity_Real[iVar];
    Velocity2_ND   += Velocity_ND[iVar]*Velocity_ND[iVar];
  }

  Physical_Constants[0] = Density_Real * Velocity2_Real / ( Density_ND * Velocity2_ND );
}

void CTransfer_FlowTraction::GetPhysical_Constants(CSolver *flow_solution, CSolver *struct_solution,
                                                   CGeometry *flow_geometry, CGeometry *struct_geometry,
                                                   CConfig *flow_config, CConfig *struct_config) {

  /*--- We have to clear the traction before applying it, because we are "adding" to node and not "setting" ---*/

  for (unsigned long iPoint = 0; iPoint < struct_geometry->GetnPoint(); iPoint++) 
    struct_solution->node[iPoint]->Clear_FlowTraction();

  Preprocess(flow_config);

  /*--- Apply a ramp to the transfer of the fluid loads ---*/

  su2double ModAmpl = 0.0;
  su2double CurrentTime = struct_config->GetCurrent_DynTime();

  bool Ramp_Load = struct_config->GetRamp_Load();
  su2double Ramp_Time = struct_config->GetRamp_Time();

  ModAmpl = struct_solution->Compute_LoadCoefficient(CurrentTime, Ramp_Time, struct_config);

  Physical_Constants[1] = ModAmpl;
  
  /*--- For static FSI, we cannot apply the ramp like this ---*/
  if ((flow_config->GetUnsteady_Simulation() == STEADY) && (struct_config->GetDynamic_Analysis() == STATIC)){
    Physical_Constants[1] = 1.0;
    if (Ramp_Load){
      CurrentTime = static_cast<su2double>(struct_config->GetOuterIter());
      Ramp_Time = static_cast<su2double>(struct_config->GetnIterFSI_Ramp() - 1);

      ModAmpl = struct_solution->Compute_LoadCoefficient(CurrentTime, Ramp_Time, struct_config);
      Physical_Constants[1] = ModAmpl;
    }
  }

  /*--- Store the force coefficient ---*/
  struct_solution->SetForceCoeff(Physical_Constants[1]);

}

void CTransfer_FlowTraction::GetDonor_Variable(CSolver *flow_solution, CGeometry *flow_geometry, CConfig *flow_config,
                                           unsigned long Marker_Flow, unsigned long Vertex_Flow, unsigned long Point_Struct) {


  unsigned short iVar, jVar;
  unsigned long Point_Flow;
  su2double const *Normal_Flow;

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
  // area: area of the face, needed if we transfer stress instead of force
  su2double Pn = 0.0, div_vel = 0.0;
  su2double Viscosity = 0.0;
  su2double Tau[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  su2double Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  su2double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  su2double area = 0.0;

  su2double Pinf = flow_solution->GetPressure_Inf();

  Point_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNode();
  // Get the normal at the vertex: this normal goes inside the fluid domain.
  Normal_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNormal();
  
  if (consistent_interpolation)
    for (iVar = 0; iVar < nVar; ++iVar) area += Normal_Flow[iVar]*Normal_Flow[iVar];
  else
    area = 1.0;
  area = sqrt(area);

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
    Donor_Variable[iVar] *= Physical_Constants[0] * Physical_Constants[1] / area;
  }

}

void CTransfer_FlowTraction::SetTarget_Variable(CSolver *fea_solution, CGeometry *fea_geometry,
                        CConfig *fea_config, unsigned long Marker_Struct,
                        unsigned long Vertex_Struct, unsigned long Point_Struct) {

  /*--- Add to the Flow traction. If nonconservative interpolation is in use,
        this is a stress and is integrated by the structural solver later on. ---*/
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

  flow_geometry->vertex[Marker_Flow][Vertex_Flow]->SetVarCoord(Target_Variable);
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

CTransfer_FlowTraction_DiscAdj::CTransfer_FlowTraction_DiscAdj(void) : CTransfer_FlowTraction() {

}

CTransfer_FlowTraction_DiscAdj::CTransfer_FlowTraction_DiscAdj(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer_FlowTraction(val_nVar, val_nConst, config) {

}

CTransfer_FlowTraction_DiscAdj::~CTransfer_FlowTraction_DiscAdj(void) {

}

void CTransfer_FlowTraction_DiscAdj::GetPhysical_Constants(CSolver *flow_solution, CSolver *struct_solution,
                                                           CGeometry *flow_geometry, CGeometry *struct_geometry,
                                                           CConfig *flow_config, CConfig *struct_config){

  /*--- We have to clear the traction before applying it, because we are "adding" to node and not "setting" ---*/

  for (unsigned long iPoint = 0; iPoint < struct_geometry->GetnPoint(); iPoint++)
    struct_solution->node[iPoint]->Clear_FlowTraction();

  Preprocess(flow_config);

  /*--- No ramp applied ---*/
  Physical_Constants[1] = 1.0;
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

CTransfer_SlidingInterface::CTransfer_SlidingInterface(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer() {

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  Physical_Constants = NULL;
  Donor_Variable     = NULL;
  Target_Variable    = NULL;

  unsigned short iVar;

  Physical_Constants = new su2double[val_nConst];
  Donor_Variable     = new su2double[val_nVar];

  Target_Variable    = new su2double[val_nVar+1];

  valAggregated      = false;

  nVar = val_nVar;

  for (iVar = 0; iVar < nVar; iVar++) {
    Donor_Variable[iVar]  = 0.0;
    Target_Variable[iVar] = 0.0;
  }

  for (iVar = 0; iVar < val_nConst; iVar++) {
    Physical_Constants[iVar] = 0.0;
  }

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

void CTransfer_SlidingInterface::InitializeTarget_Variable(CSolver *target_solution, unsigned long Marker_Target,
                          unsigned long Vertex_Target, unsigned short nDonorPoints) {

  target_solution->SetnSlidingStates(Marker_Target, Vertex_Target, nDonorPoints); // This is to allocate
  target_solution->SetSlidingStateStructure(Marker_Target, Vertex_Target);
  target_solution->SetnSlidingStates(Marker_Target, Vertex_Target, 0); // Reset counter to 0

}

void CTransfer_SlidingInterface::RecoverTarget_Variable(long indexPoint_iVertex, su2double *Buffer_Bcast_Variables,
                                                        su2double donorCoeff){
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Target_Variable[iVar] = Buffer_Bcast_Variables[ indexPoint_iVertex*nVar + iVar ];

  Target_Variable[nVar] = donorCoeff;
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

CTransfer_ConjugateHeatVars::CTransfer_ConjugateHeatVars(void) : CTransfer() {

}

CTransfer_ConjugateHeatVars::CTransfer_ConjugateHeatVars(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) : CTransfer(val_nVar, val_nConst, config) {

}

CTransfer_ConjugateHeatVars::~CTransfer_ConjugateHeatVars(void) {

}

void CTransfer_ConjugateHeatVars::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry, CConfig *donor_config,
                                                unsigned long Marker_Donor, unsigned long Vertex_Donor, unsigned long Point_Donor) {

  unsigned long iPoint;
  unsigned long PointNormal;
  unsigned short nDim, iDim;

  nDim = donor_geometry->GetnDim();

  su2double *Coord, *Coord_Normal, *Normal, *Edge_Vector, dist, dist2, Area, Twall = 0.0, Tnormal = 0.0,
      dTdn, cp_fluid, rho_cp_solid, Prandtl_Lam, laminar_viscosity,
      thermal_diffusivity, thermal_conductivity, thermal_conductivityND, heat_flux_density, conductivity_over_dist, Temperature_Ref;
  su2double Gamma = donor_config->GetGamma();
  su2double Gas_Constant = donor_config->GetGas_ConstantND();
  su2double Cp = (Gamma / (Gamma - 1.0)) * Gas_Constant;
  Edge_Vector = new su2double[nDim];

  /*--- Check whether the current zone is a solid zone or a fluid zone ---*/
  bool flow = ((donor_config->GetKind_Solver() == NAVIER_STOKES)
               || (donor_config->GetKind_Solver() == RANS)
               || (donor_config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (donor_config->GetKind_Solver() == DISC_ADJ_RANS));
  bool compressible_flow  = (donor_config->GetKind_Regime() == COMPRESSIBLE) && flow;
  bool incompressible_flow = (donor_config->GetEnergy_Equation()) && flow;
  bool heat_equation      = donor_config->GetKind_Solver() == HEAT_EQUATION_FVM;

  Temperature_Ref   = donor_config->GetTemperature_Ref();
  Prandtl_Lam       = donor_config->GetPrandtl_Lam();
  laminar_viscosity = donor_config->GetMu_ConstantND(); // TDE check for consistency
  cp_fluid          = donor_config->GetSpecific_Heat_Cp();
  rho_cp_solid      = donor_config->GetSpecific_Heat_Cp_Solid()*donor_config->GetDensity_Solid();

  PointNormal   = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNormal_Neighbor();
  Coord         = donor_geometry->node[Point_Donor]->GetCoord();
  Coord_Normal  = donor_geometry->node[PointNormal]->GetCoord();
  Normal        = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNormal();

  dist2 = 0.0;
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_Normal[iDim] - Coord[iDim];
    dist2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  dist = sqrt(dist2);
  Area = sqrt(Area);

  /*--- Retrieve temperature solution (later set is as the first donor variable) and its gradient ---*/

  dTdn = 0.0;

  if (compressible_flow) {

    Twall   = donor_solution->node[Point_Donor]->GetPrimitive(0)*Temperature_Ref;
    Tnormal = donor_solution->node[PointNormal]->GetPrimitive(0)*Temperature_Ref;

    dTdn = (Twall - Tnormal)/dist;
  }
  else if (incompressible_flow) {

    Twall   = donor_solution->node[Point_Donor]->GetTemperature()*Temperature_Ref;
    Tnormal = donor_solution->node[PointNormal]->GetTemperature()*Temperature_Ref;

    dTdn = (Twall - Tnormal)/dist;
  }
  else if (flow || heat_equation) {
    Twall   = donor_solution->node[Point_Donor]->GetSolution(0)*Temperature_Ref;
    Tnormal = donor_solution->node[PointNormal]->GetSolution(0)*Temperature_Ref;

//    for (iDim = 0; iDim < nDim; iDim++) {
//      dTdn += (Twall - Tnormal)/dist * (Edge_Vector[iDim]/dist) * (Normal[iDim]/Area);
//    }

    dTdn = (Twall - Tnormal)/dist;
  }
  else {
    cout << "WARNING: Transfer of conjugate heat variables is called with non-supported donor solver!" << endl;
  }

  /*--- Calculate the heat flux density (temperature gradient times thermal conductivity) and set it as second donor variable ---*/
  if (compressible_flow) {

    iPoint = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNode();

    thermal_conductivityND  = Cp*(laminar_viscosity/Prandtl_Lam);
    thermal_conductivity    = thermal_conductivityND*donor_config->GetViscosity_Ref();

    heat_flux_density       = thermal_conductivity*dTdn;
    conductivity_over_dist  = thermal_conductivity/dist;
  }
  else if (incompressible_flow) {

    iPoint = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNode();

    thermal_conductivityND  = donor_solution->node[iPoint]->GetThermalConductivity();
    thermal_conductivity = thermal_conductivityND*donor_config->GetConductivity_Ref();

    switch (donor_config->GetKind_ConductivityModel()) {

      case CONSTANT_CONDUCTIVITY:
        thermal_conductivity = thermal_conductivityND*donor_config->GetConductivity_Ref();
        break;

      case CONSTANT_PRANDTL:
        thermal_conductivity = thermal_conductivityND*donor_config->GetGas_Constant_Ref()*donor_config->GetViscosity_Ref();
        break;
    }

    heat_flux_density       = thermal_conductivity*dTdn;
    conductivity_over_dist  = thermal_conductivity/dist;
  }
  else if (flow) {

    iPoint = donor_geometry->vertex[Marker_Donor][Vertex_Donor]->GetNode();

    thermal_conductivityND  = laminar_viscosity/Prandtl_Lam;
    thermal_conductivity    = thermal_conductivityND*donor_config->GetViscosity_Ref()*cp_fluid;

    heat_flux_density       = thermal_conductivity*dTdn;
    conductivity_over_dist  = thermal_conductivity/dist;
  }
  else {

    thermal_diffusivity     = donor_config->GetThermalDiffusivity_Solid();
    heat_flux_density       = (thermal_diffusivity*dTdn)*rho_cp_solid;
    conductivity_over_dist  = thermal_diffusivity*rho_cp_solid/dist;
  }

  Donor_Variable[0] = Twall;
  Donor_Variable[1] = heat_flux_density;
  Donor_Variable[2] = conductivity_over_dist;
  Donor_Variable[3] = Tnormal;

  delete [] Edge_Vector;
}

void CTransfer_ConjugateHeatVars::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                          CConfig *target_config, unsigned long Marker_Target,
                          unsigned long Vertex_Target, unsigned long Point_Target) {

  target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 0, target_config->GetRelaxation_Factor_CHT(), Target_Variable[0]);
  target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 1, target_config->GetRelaxation_Factor_CHT(), Target_Variable[1]);
  target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 2, target_config->GetRelaxation_Factor_CHT(), Target_Variable[2]);
  target_solution->SetConjugateHeatVariable(Marker_Target, Vertex_Target, 3, target_config->GetRelaxation_Factor_CHT(), Target_Variable[3]);
}
