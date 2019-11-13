/*!
 * \file CFlowTractionInterface.cpp
 * \brief Declaration and inlines of the class to transfer flow tractions
 *        from a fluid zone into a structural zone.
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

#include "../../../include/interfaces/fsi/CFlowTractionInterface.hpp"

CFlowTractionInterface::CFlowTractionInterface(void) : CInterface() {

}

CFlowTractionInterface::CFlowTractionInterface(unsigned short val_nVar, unsigned short val_nConst, CConfig *config) :
  CInterface(val_nVar, val_nConst, config) {

}

CFlowTractionInterface::~CFlowTractionInterface(void) {

}

void CFlowTractionInterface::Preprocess(CConfig *flow_config) {

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

void CFlowTractionInterface::GetPhysical_Constants(CSolver *flow_solution, CSolver *struct_solution,
                                                   CGeometry *flow_geometry, CGeometry *struct_geometry,
                                                   CConfig *flow_config, CConfig *struct_config) {

  /*--- We have to clear the traction before applying it, because we are "adding" to node and not "setting" ---*/

  struct_solution->GetNodes()->Clear_FlowTraction();

  Preprocess(flow_config);

  /*--- Apply a ramp to the transfer of the fluid loads ---*/

  su2double ModAmpl = 0.0;
  su2double CurrentTime = struct_config->GetCurrent_DynTime();

  bool Ramp_Load = struct_config->GetRamp_Load();
  su2double Ramp_Time = struct_config->GetRamp_Time();

  ModAmpl = struct_solution->Compute_LoadCoefficient(CurrentTime, Ramp_Time, struct_config);

  Physical_Constants[1] = ModAmpl;
  
  /*--- For static FSI, we cannot apply the ramp like this ---*/
  if ((!flow_config->GetTime_Domain())){
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

void CFlowTractionInterface::GetDonor_Variable(CSolver *flow_solution, CGeometry *flow_geometry,
                                               CConfig *flow_config, unsigned long Marker_Flow,
                                               unsigned long Vertex_Flow, unsigned long Point_Struct) {


  unsigned short iVar, jVar;
  unsigned long Point_Flow;
  su2double const *Normal_Flow;

  // Check the kind of fluid problem
  bool compressible       = (flow_config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (flow_config->GetKind_Regime() == INCOMPRESSIBLE);
  bool viscous_flow       = ((flow_config->GetKind_Solver() == NAVIER_STOKES) ||
                             (flow_config->GetKind_Solver() == RANS) ||
                             (flow_config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES) ||
                             (flow_config->GetKind_Solver() == DISC_ADJ_RANS) ||
                             (flow_config->GetKind_Solver() == DISC_ADJ_INC_NAVIER_STOKES) ||
                             (flow_config->GetKind_Solver() == DISC_ADJ_INC_RANS));

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

  Pn = flow_solution->GetNodes()->GetPressure(Point_Flow);

  // Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
  for (iVar = 0; iVar < nVar; iVar++)
    Donor_Variable[iVar] = -(Pn-Pinf)*Normal_Flow[iVar];

  // Calculate tn in the fluid nodes for the viscous term

  if ((incompressible || compressible) && viscous_flow) {

    Viscosity = flow_solution->GetNodes()->GetLaminarViscosity(Point_Flow);

    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0 ; jVar < nVar; jVar++) {
        Grad_Vel[iVar][jVar] = flow_solution->GetNodes()->GetGradient_Primitive(Point_Flow, iVar+1, jVar);
      }
    }

    // Divergence of the velocity
    div_vel = 0.0; for (iVar = 0; iVar < nVar; iVar++) div_vel += Grad_Vel[iVar][iVar];

    for (iVar = 0; iVar < nVar; iVar++) {

      for (jVar = 0 ; jVar < nVar; jVar++) {

        // Viscous stress
        Tau[iVar][jVar] = Viscosity*(Grad_Vel[jVar][iVar] + Grad_Vel[iVar][jVar])
                          - TWO3*Viscosity*div_vel*delta[iVar][jVar];

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

void CFlowTractionInterface::SetTarget_Variable(CSolver *fea_solution, CGeometry *fea_geometry,
                                                CConfig *fea_config, unsigned long Marker_Struct,
                                                unsigned long Vertex_Struct, unsigned long Point_Struct) {

  /*--- Add to the Flow traction. If nonconservative interpolation is in use,
        this is a stress and is integrated by the structural solver later on. ---*/
  fea_solution->GetNodes()->Add_FlowTraction(Point_Struct,Target_Variable);

}
