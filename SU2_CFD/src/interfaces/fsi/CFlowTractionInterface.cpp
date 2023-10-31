/*!
 * \file CFlowTractionInterface.cpp
 * \brief Declaration and inlines of the class to transfer flow tractions
 *        from a fluid zone into a structural zone.
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

#include "../../../include/interfaces/fsi/CFlowTractionInterface.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../include/solvers/CSolver.hpp"
#include "../../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../../Common/include/interface_interpolation/CInterpolator.hpp"
#include <unordered_set>

CFlowTractionInterface::CFlowTractionInterface(unsigned short val_nVar, unsigned short val_nConst,
                                               const CConfig *config, bool conservative_) :
  CInterface(val_nVar, val_nConst),
  conservative(conservative_) {
}

void CFlowTractionInterface::Preprocess(const CConfig *flow_config, const CConfig *struct_config,
                                        CGeometry *struct_geometry, CSolver *struct_solution) {

  /*--- Clear the tractions only on the markers involved in interface, fluid tractions
   * on other markers can be specified via e.g. the python wrapper. ---*/

  for (auto iMarkerInt = 0u; iMarkerInt < struct_config->GetMarker_n_ZoneInterface()/2; iMarkerInt++) {
    const auto markDonor = flow_config->FindInterfaceMarker(iMarkerInt);
    const auto markTarget = struct_config->FindInterfaceMarker(iMarkerInt);

    if(!CInterpolator::CheckInterfaceBoundary(markDonor, markTarget)) continue;
    if (markTarget < 0) continue;

    for (auto iVertex = 0ul; iVertex < struct_geometry->GetnVertex(markTarget); iVertex++) {
      const auto iPoint = struct_geometry->vertex[markTarget][iVertex]->GetNode();
      if (!struct_geometry->nodes->GetDomain(iPoint)) continue;
      su2double zeros[3] = {};
      struct_solution->GetNodes()->Set_FlowTraction(iPoint, zeros);
    }
  }

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

void CFlowTractionInterface::ComputeVertexAreas(const CConfig *config,
                                                CGeometry *geometry,
                                                CSolver *solution) {
  constexpr size_t MAXNDIM = 3;
  const auto nodes = solution->GetNodes();
  const auto nDim = nVar;

  /*--- Compute current area associated with each vertex. ---*/

  vertexArea.clear();
  unordered_set<short> processedMarkers({-1});

  for (auto iMarkerInt = 0; iMarkerInt < config->GetMarker_n_ZoneInterface()/2; ++iMarkerInt) {
    /*--- Find the marker index associated with the pair. ---*/
    const auto iMarker = config->FindInterfaceMarker(iMarkerInt);
    /*--- The current mpi rank may not have this marker, or it may have been processed already. ---*/
    if (processedMarkers.count(iMarker) > 0) continue;
    processedMarkers.insert(iMarker);

    for (auto iElem = 0u; iElem < geometry->GetnElem_Bound(iMarker); ++iElem) {
      /*--- Define the boundary element. ---*/
      unsigned long nodeList[N_POINTS_MAXIMUM] = {0};
      su2double coords[N_POINTS_MAXIMUM][MAXNDIM] = {{0.0}};
      bool quad = geometry->bound[iMarker][iElem]->GetVTK_Type() == QUADRILATERAL;
      auto nNode = quad? 4u : nDim;

      for (auto iNode = 0u; iNode < nNode; ++iNode) {
        nodeList[iNode] = geometry->bound[iMarker][iElem]->GetNode(iNode);
        for (auto iDim = 0u; iDim < nDim; ++iDim)
          coords[iNode][iDim] = geometry->nodes->GetCoord(nodeList[iNode], iDim)+
                                nodes->GetSolution(nodeList[iNode],iDim);
      }

      /*--- Compute the area contribution to each node. ---*/
      su2double normal[MAXNDIM] = {0.0};

      switch (nNode) {
        case 2: GeometryToolbox::LineNormal(coords, normal); break;
        case 3: GeometryToolbox::TriangleNormal(coords, normal); break;
        case 4: GeometryToolbox::QuadrilateralNormal(coords, normal); break;
      }

      su2double area = GeometryToolbox::Norm(MAXNDIM, normal) / nNode;

      /*--- Update area of nodes. ---*/
      for (auto iNode = 0u; iNode < nNode; ++iNode) {
        auto iPoint = nodeList[iNode];
        if (vertexArea.count(iPoint) == 0)
          vertexArea[iPoint] = area;
        else
          vertexArea[iPoint] += area;
      }
    }
  }
}

void CFlowTractionInterface::GetPhysical_Constants(CSolver *flow_solution, CSolver *struct_solution,
                                                   CGeometry *flow_geometry, CGeometry *struct_geometry,
                                                   const CConfig *flow_config, const CConfig *struct_config) {

  Preprocess(flow_config, struct_config, struct_geometry, struct_solution);

  if (!conservative) ComputeVertexAreas(struct_config, struct_geometry, struct_solution);

  /*--- Apply a ramp to the transfer of the fluid loads ---*/

  su2double ModAmpl = 0.0;
  su2double CurrentTime = struct_config->GetCurrent_UnstTime();

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
                                               const CConfig *flow_config, unsigned long Marker_Flow,
                                               unsigned long Vertex_Flow, unsigned long Point_Struct) {

  const auto Point_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNode();

  // Get the normal at the vertex: this normal goes inside the fluid domain.
  const su2double* Normal_Flow = flow_geometry->vertex[Marker_Flow][Vertex_Flow]->GetNormal();

  // If we do not want integrated tractions, i.e. forces, we will need to divide by area.
  su2double oneOnArea = 1.0;
  if (!conservative) oneOnArea = 1.0 / GeometryToolbox::Norm(nVar, Normal_Flow);

  // Retrieve the values of pressure

  CVariable* flow_nodes = flow_solution->GetNodes();
  su2double Pinf = flow_solution->GetPressure_Inf();
  su2double Pn = flow_nodes->GetPressure(Point_Flow);

  // Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).

  for (auto iVar = 0u; iVar < nVar; iVar++)
    Donor_Variable[iVar] = -(Pn-Pinf)*Normal_Flow[iVar];

  // Calculate tn in the fluid nodes for the viscous term

  if (flow_config->GetViscous()) {

    su2double Viscosity = flow_nodes->GetLaminarViscosity(Point_Flow);

    su2double tau[3][3];
    CNumerics::ComputeStressTensor(nVar, tau, flow_nodes->GetVelocityGradient(Point_Flow), Viscosity);
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      for (auto jVar = 0u; jVar < nVar; jVar++) {
        // Viscous component in the tn vector --> Units of force (non-dimensional).
        Donor_Variable[iVar] += tau[iVar][jVar] * Normal_Flow[jVar];
      }
    }

  }

  // Redimensionalize and take into account ramp transfer of the loads
  for (auto iVar = 0u; iVar < nVar; iVar++) {
    Donor_Variable[iVar] *= Physical_Constants[0] * Physical_Constants[1] * oneOnArea;
  }

}

void CFlowTractionInterface::SetTarget_Variable(CSolver *fea_solution, CGeometry *fea_geometry,
                                                const CConfig *fea_config, unsigned long Marker_Struct,
                                                unsigned long Vertex_Struct, unsigned long Point_Struct) {

  /*--- Add the Flow traction, if consistent interpolation is in use the traction needs to be integrated. ---*/
  if (!conservative) {
    assert(vertexArea.count(Point_Struct));
    const auto area = vertexArea.at(Point_Struct);
    for (auto iVar = 0u; iVar < nVar; ++iVar) Target_Variable[iVar] *= area;
  }

  fea_solution->GetNodes()->Add_FlowTraction(Point_Struct, Target_Variable);
}
