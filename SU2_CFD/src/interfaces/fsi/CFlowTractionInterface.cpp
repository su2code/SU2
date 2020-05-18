/*!
 * \file CFlowTractionInterface.cpp
 * \brief Declaration and inlines of the class to transfer flow tractions
 *        from a fluid zone into a structural zone.
 * \author R. Sanchez
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


CFlowTractionInterface::CFlowTractionInterface(unsigned short val_nVar, unsigned short val_nConst,
                                               CConfig *config, bool integrate_tractions_) :
  CInterface(val_nVar, val_nConst, config),
  integrate_tractions(integrate_tractions_) {

}

void CFlowTractionInterface::Preprocess(CConfig *flow_config) {

  /*--- Compute the constant factor to dimensionalize pressure and shear stress. ---*/
  const su2double *Velocity_ND, *Velocity_Real;
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

  // Check the kind of fluid problem
  bool viscous_flow;
  switch (flow_config->GetKind_Solver()) {
    case RANS: case INC_RANS:
    case NAVIER_STOKES: case INC_NAVIER_STOKES:
    case DISC_ADJ_RANS: case DISC_ADJ_INC_RANS:
    case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_INC_NAVIER_STOKES:
      viscous_flow = true; break;
    default:
      viscous_flow = false; break;
  }

  const auto Point_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNode();

  // Get the normal at the vertex: this normal goes inside the fluid domain.
  const su2double* Normal_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNormal();

  // If we do not want integrated tractions, i.e. forces, we will need to divide by area.
  su2double oneOnArea = 1.0;
  if (!integrate_tractions) {
    oneOnArea = 0.0;
    for (iVar = 0; iVar < nVar; ++iVar) oneOnArea += pow(Normal_Flow[iVar], 2);
    oneOnArea = 1.0 / sqrt(oneOnArea);
  }

  // Retrieve the values of pressure

  CVariable* flow_nodes = flow_solution->GetNodes();
  su2double Pinf = flow_solution->GetPressure_Inf();
  su2double Pn = flow_nodes->GetPressure(Point_Flow);

  // Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).

  for (iVar = 0; iVar < nVar; iVar++)
    Donor_Variable[iVar] = -(Pn-Pinf)*Normal_Flow[iVar];

  // Calculate tn in the fluid nodes for the viscous term

  if (viscous_flow) {

    su2double Viscosity = flow_nodes->GetLaminarViscosity(Point_Flow);

    const su2double* const* GradVel = &flow_nodes->GetGradient_Primitive(Point_Flow)[1];

    // Divergence of the velocity
    su2double DivVel = 0.0;
    for (iVar = 0; iVar < nVar; iVar++) DivVel += GradVel[iVar][iVar];

    for (iVar = 0; iVar < nVar; iVar++) {

      for (jVar = 0 ; jVar < nVar; jVar++) {

        // Viscous stress
        su2double delta_ij = (iVar == jVar);
        su2double tau_ij = Viscosity*(GradVel[jVar][iVar] + GradVel[iVar][jVar] - TWO3*DivVel*delta_ij);

        // Viscous component in the tn vector --> Units of force (non-dimensional).
        Donor_Variable[iVar] += tau_ij * Normal_Flow[jVar];
      }
    }
  }

  // Redimensionalize and take into account ramp transfer of the loads
  for (iVar = 0; iVar < nVar; iVar++) {
    Donor_Variable[iVar] *= Physical_Constants[0] * Physical_Constants[1] * oneOnArea;
  }

}

void CFlowTractionInterface::SetTarget_Variable(CSolver *fea_solution, CGeometry *fea_geometry,
                                                CConfig *fea_config, unsigned long Marker_Struct,
                                                unsigned long Vertex_Struct, unsigned long Point_Struct) {

  /*--- Add to the Flow traction. If nonconservative interpolation is in use,
        this is a stress and is integrated by the structural solver later on. ---*/
  fea_solution->GetNodes()->Add_FlowTraction(Point_Struct,Target_Variable);

}
