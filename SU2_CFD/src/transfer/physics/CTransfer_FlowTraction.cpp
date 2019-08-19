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

#include "../include/transfer/CTransfer.hpp"

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


