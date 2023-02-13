/*!
 * \file python_wrapper_structure.cpp
 * \brief Driver subroutines that are used by the Python wrapper. Those routines are usually called from an external Python environment.
 * \author D. Thomas
 * \version 7.5.1 "Blackbird"
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


#include "../include/drivers/CDriver.hpp"
#include "../include/drivers/CSinglezoneDriver.hpp"
#include "../../Common/include/toolboxes/geometry_toolbox.hpp"

void CDriver::PythonInterface_Preprocessing(CConfig **config, CGeometry ****geometry, CSolver *****solver){

  int rank = MASTER_NODE;
  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);

  /* --- Initialize boundary conditions customization, this is achieve through the Python wrapper --- */
  for(iZone=0; iZone < nZone; iZone++){

    if (config[iZone]->GetnMarker_PyCustom() > 0){

      if (rank == MASTER_NODE) cout << endl << "----------------- Python Interface Preprocessing ( Zone "<< iZone <<" ) -----------------" << endl;

      if (rank == MASTER_NODE) cout << "Setting customized boundary conditions for zone " << iZone << endl;
      for (iMesh = 0; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
        geometry[iZone][INST_0][iMesh]->SetCustomBoundary(config[iZone]);
      }
      geometry[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry[iZone][INST_0], config[iZone]);

      if ((config[iZone]->GetKind_Solver() == MAIN_SOLVER::EULER) ||
          (config[iZone]->GetKind_Solver() == MAIN_SOLVER::NAVIER_STOKES) ||
          (config[iZone]->GetKind_Solver() == MAIN_SOLVER::RANS)) {

        solver[iZone][INST_0][MESH_0][FLOW_SOL]->UpdateCustomBoundaryConditions(geometry[iZone][INST_0], config[iZone]);
      }
    }
  }

}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to the global performance indices (Lift, Drag, ecc..) */
/////////////////////////////////////////////////////////////////////////////

passivedouble CDriver::Get_Drag() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CDrag, factor, val_Drag;

  /*--- Calculate drag force based on drag coefficient ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CDrag = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CD();

  val_Drag = CDrag*factor;

  return SU2_TYPE::GetValue(val_Drag);
}

passivedouble CDriver::Get_Lift() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CLift, factor, val_Lift;

  /*--- Calculate drag force based on drag coefficient ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CLift = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CL();

  val_Lift = CLift*factor;

  return SU2_TYPE::GetValue(val_Lift);
}

passivedouble CDriver::Get_Mx() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMx, RefLengthCoeff, factor, val_Mx;

  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate moment around x-axis based on coefficients ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CMx = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMx();

  val_Mx = CMx*factor*RefLengthCoeff;

  return SU2_TYPE::GetValue(val_Mx);
}

passivedouble CDriver::Get_My() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMy, RefLengthCoeff, factor, val_My;

  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate moment around x-axis based on coefficients ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CMy = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMy();

  val_My = CMy*factor*RefLengthCoeff;

  return SU2_TYPE::GetValue(val_My);
}

passivedouble CDriver::Get_Mz() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMz, RefLengthCoeff, factor, val_Mz;

  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate moment around z-axis based on coefficients ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CMz = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMz();

  val_Mz = CMz*factor*RefLengthCoeff;

  return SU2_TYPE::GetValue(val_Mz);
}

passivedouble CDriver::Get_DragCoeff() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CDrag;

  CDrag = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CD();

  return SU2_TYPE::GetValue(CDrag);
}

passivedouble CDriver::Get_LiftCoeff() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CLift;

  CLift = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CL();

  return SU2_TYPE::GetValue(CLift);
}

/////////////////////////////////////////////////////////////////////////////
/* Functions to obtain information from the geometry/mesh                  */
/////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetNumberVertices(unsigned short iMarker) const {

  return geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker];

}

unsigned long CDriver::GetNumberHaloVertices(unsigned short iMarker) const {

  unsigned long nHaloVertices, iVertex, iPoint;

  nHaloVertices = 0;
  for(iVertex = 0; iVertex < geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker]; iVertex++){
    iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
    if(!(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint))) nHaloVertices += 1;
  }

  return nHaloVertices;

}

unsigned long CDriver::GetVertexGlobalIndex(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint, GlobalIndex;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  GlobalIndex = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetGlobalIndex(iPoint);

  return GlobalIndex;

}

bool CDriver::IsAHaloNode(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  if(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint)) return false;
  else return true;

}

vector<passivedouble> CDriver::GetInitialMeshCoord(unsigned short iMarker, unsigned long iVertex) const {

  vector<su2double> coord(3,0.0);
  vector<passivedouble> coord_passive(3, 0.0);

  auto iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  for (auto iDim = 0 ; iDim < nDim ; iDim++){
   coord[iDim] = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->GetMesh_Coord(iPoint,iDim);
  }

  coord_passive[0] = SU2_TYPE::GetValue(coord[0]);
  coord_passive[1] = SU2_TYPE::GetValue(coord[1]);
  coord_passive[2] = SU2_TYPE::GetValue(coord[2]);

  return coord_passive;
}

vector<passivedouble> CDriver::GetVertexNormal(unsigned short iMarker, unsigned long iVertex, bool unitNormal) const {

  su2double *Normal;
  su2double Area;
  vector<su2double> ret_Normal(3, 0.0);
  vector<passivedouble> ret_Normal_passive(3, 0.0);

  Normal = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();

  if (!unitNormal) {

    ret_Normal_passive[0] = SU2_TYPE::GetValue(Normal[0]);
    ret_Normal_passive[1] = SU2_TYPE::GetValue(Normal[1]);
    if(nDim>2) ret_Normal_passive[2] = SU2_TYPE::GetValue(Normal[2]);

    return ret_Normal_passive;
  }

  Area = GeometryToolbox::Norm(nDim, Normal);

  ret_Normal[0] = Normal[0]/Area;
  ret_Normal[1] = Normal[1]/Area;
  if(nDim>2) ret_Normal[2] = Normal[2]/Area;

  ret_Normal_passive[0] = SU2_TYPE::GetValue(ret_Normal[0]);
  ret_Normal_passive[1] = SU2_TYPE::GetValue(ret_Normal[1]);
  ret_Normal_passive[2] = SU2_TYPE::GetValue(ret_Normal[2]);

  return ret_Normal_passive;
}

//////////////////////////////////////////////////////////////////////////////////
/* Functions to obtain global parameters from SU2 (time steps, delta t, ecc...) */
//////////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetnTimeIter() const {

  return config_container[ZONE_0]->GetnTime_Iter();
}

unsigned long CDriver::GetTime_Iter() const{

  return TimeIter;
}

passivedouble CDriver::GetUnsteady_TimeStep() const {

  return SU2_TYPE::GetValue(config_container[ZONE_0]->GetTime_Step());
}

string CDriver::GetSurfaceFileName() const {

  return config_container[ZONE_0]->GetSurfCoeff_FileName();
}

///////////////////////////////////////////////////////////////////////////////
/* Functions related to CHT solver                                           */
///////////////////////////////////////////////////////////////////////////////

passivedouble CDriver::GetVertexTemperature(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  su2double vertexWallTemp(0.0);

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint) && compressible){
    vertexWallTemp = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
  }

  return SU2_TYPE::GetValue(vertexWallTemp);

}

void CDriver::SetVertexTemperature(unsigned short iMarker, unsigned long iVertex, passivedouble val_WallTemp){

  geometry_container[ZONE_0][INST_0][MESH_0]->SetCustomBoundaryTemperature(iMarker, iVertex, val_WallTemp);
}

vector<passivedouble> CDriver::GetVertexHeatFluxes(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  unsigned short iDim;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double laminar_viscosity, thermal_conductivity;
  vector<su2double> GradT (3,0.0);
  vector<su2double> HeatFlux (3,0.0);
  vector<passivedouble> HeatFluxPassive (3,0.0);

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(compressible){
    laminar_viscosity    = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);
    for(iDim=0; iDim < nDim; iDim++){
      GradT[iDim] = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint, 0, iDim);
      HeatFlux[iDim] = -thermal_conductivity*GradT[iDim];
    }
  }

  HeatFluxPassive[0] = SU2_TYPE::GetValue(HeatFlux[0]);
  HeatFluxPassive[1] = SU2_TYPE::GetValue(HeatFlux[1]);
  HeatFluxPassive[2] = SU2_TYPE::GetValue(HeatFlux[2]);

  return HeatFluxPassive;
}

passivedouble CDriver::GetVertexNormalHeatFlux(unsigned short iMarker, unsigned long iVertex) const{

  unsigned long iPoint;
  unsigned short iDim;
  su2double vertexWallHeatFlux;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double Area;
  su2double laminar_viscosity, thermal_conductivity, dTdn;
  su2double *Normal, GradT[3] = {0.0,0.0,0.0}, UnitNormal[3] = {0.0,0.0,0.0};

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  vertexWallHeatFlux = 0.0;
  dTdn = 0.0;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint) && compressible){
    Normal = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();

    Area = GeometryToolbox::Norm(nDim, Normal);

    for (iDim = 0; iDim < nDim; iDim++)
      UnitNormal[iDim] = Normal[iDim]/Area;

    laminar_viscosity    = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);
    /*Compute wall heat flux (normal to the wall) based on computed temperature gradient*/
    for(iDim=0; iDim < nDim; iDim++){
      GradT[iDim] = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint, 0, iDim);
      dTdn += GradT[iDim]*UnitNormal[iDim];
    }

    vertexWallHeatFlux = -thermal_conductivity*dTdn;
  }

  return SU2_TYPE::GetValue(vertexWallHeatFlux);
}

void CDriver::SetVertexNormalHeatFlux(unsigned short iMarker, unsigned long iVertex, passivedouble val_WallHeatFlux){

  geometry_container[ZONE_0][INST_0][MESH_0]->SetCustomBoundaryHeatFlux(iMarker, iVertex, val_WallHeatFlux);
}

passivedouble CDriver::GetThermalConductivity(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double laminar_viscosity, thermal_conductivity;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  laminar_viscosity    = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
  thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);

  return SU2_TYPE::GetValue(thermal_conductivity);

}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to the management of markers                             */
////////////////////////////////////////////////////////////////////////////////

vector<string> CDriver::GetAllBoundaryMarkersTag() const {

  vector<string> boundariesTagList;
  unsigned short iMarker,nBoundariesMarkers;
  string Marker_Tag;

  nBoundariesMarkers = config_container[ZONE_0]->GetnMarker_All();
  boundariesTagList.resize(nBoundariesMarkers);

  for(iMarker=0; iMarker < nBoundariesMarkers; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
    boundariesTagList[iMarker] = Marker_Tag;
  }

  return boundariesTagList;
}

vector<string> CDriver::GetAllDeformMeshMarkersTag() const {

  vector<string> interfaceBoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_Deform_Mesh();
  interfaceBoundariesTagList.resize(nBoundariesMarker);

  for(iMarker=0; iMarker < nBoundariesMarker; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_Deform_Mesh_TagBound(iMarker);
    interfaceBoundariesTagList[iMarker] = Marker_Tag;
  }

  return interfaceBoundariesTagList;
}

vector<string> CDriver::GetAllCHTMarkersTag() const {

  vector<string> CHTBoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_All();

  //The CHT markers can be identified as the markers that are customizable with a BC type HEAT_FLUX or ISOTHERMAL.
  for(iMarker=0; iMarker<nBoundariesMarker; iMarker++){
    if((config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == HEAT_FLUX || config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) && config_container[ZONE_0]->GetMarker_All_PyCustom(iMarker)){
      Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
      CHTBoundariesTagList.push_back(Marker_Tag);
    }
  }

  return CHTBoundariesTagList;
}

vector<string> CDriver::GetAllInletMarkersTag() const {

  vector<string> BoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_All();

  for(iMarker=0; iMarker<nBoundariesMarker; iMarker++){
    bool isCustomizable = config_container[ZONE_0]->GetMarker_All_PyCustom(iMarker);
    bool isInlet = (config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == INLET_FLOW);
    if(isCustomizable && isInlet) {
      Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
      BoundariesTagList.push_back(Marker_Tag);
    }
  }

  return BoundariesTagList;
}

map<string, int> CDriver::GetAllBoundaryMarkers() const {

  map<string, int>  allBoundariesMap;
  unsigned short iMarker, nBoundaryMarkers;
  string Marker_Tag;

  nBoundaryMarkers = config_container[ZONE_0]->GetnMarker_All();

  for(iMarker=0; iMarker < nBoundaryMarkers; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
    allBoundariesMap[Marker_Tag] = iMarker;
  }

  return allBoundariesMap;
}

map<string, string> CDriver::GetAllBoundaryMarkersType() const {

  map<string, string> allBoundariesTypeMap;
  unsigned short iMarker, KindBC;
  string Marker_Tag, Marker_Type;

  for(iMarker=0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
    KindBC = config_container[ZONE_0]->GetMarker_All_KindBC(iMarker);
    switch(KindBC){
      case EULER_WALL:
        Marker_Type = "EULER_WALL";
        break;
      case FAR_FIELD:
        Marker_Type = "FARFIELD";
        break;
      case ISOTHERMAL:
        Marker_Type = "ISOTHERMAL";
        break;
      case HEAT_FLUX:
        Marker_Type = "HEATFLUX";
        break;
      case INLET_FLOW:
        Marker_Type = "INLET_FLOW";
        break;
      case OUTLET_FLOW:
        Marker_Type = "OUTLET_FLOW";
        break;
      case SYMMETRY_PLANE:
        Marker_Type = "SYMMETRY";
        break;
      case SEND_RECEIVE:
        Marker_Type = "SEND_RECEIVE";
        break;
      default:
        Marker_Type = "UNKNOWN_TYPE";
    }
    allBoundariesTypeMap[Marker_Tag] = Marker_Type;
  }

  return allBoundariesTypeMap;
}

void CDriver::SetHeatSource_Position(passivedouble alpha, passivedouble pos_x, passivedouble pos_y, passivedouble pos_z){

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][RAD_SOL];

  config_container[ZONE_0]->SetHeatSource_Rot_Z(alpha);
  config_container[ZONE_0]->SetHeatSource_Center(pos_x, pos_y, pos_z);

  solver->SetVolumetricHeatSource(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0]);

}

void CDriver::SetInlet_Angle(unsigned short iMarker, passivedouble alpha){

  su2double alpha_rad = alpha * PI_NUMBER/180.0;

  unsigned long iVertex;

  for (iVertex = 0; iVertex < geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker]; iVertex++){
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInlet_FlowDir(iMarker, iVertex, 0, cos(alpha_rad));
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInlet_FlowDir(iMarker, iVertex, 1, sin(alpha_rad));
  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Functions related to simulation control, high level functions (reset convergence, set initial mesh, ecc...) */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CDriver::ResetConvergence() {

  for(iZone = 0; iZone < nZone; iZone++) {
    switch (config_container[iZone]->GetKind_Solver()) {

    case MAIN_SOLVER::EULER: case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS:
    case MAIN_SOLVER::INC_EULER: case MAIN_SOLVER::INC_NAVIER_STOKES: case MAIN_SOLVER::INC_RANS:
      integration_container[iZone][INST_0][FLOW_SOL]->SetConvergence(false);
      if (config_container[iZone]->GetKind_Solver() == MAIN_SOLVER::RANS) integration_container[iZone][INST_0][TURB_SOL]->SetConvergence(false);
      if(config_container[iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) integration_container[iZone][INST_0][TRANS_SOL]->SetConvergence(false);
      break;

    case MAIN_SOLVER::FEM_ELASTICITY:
      integration_container[iZone][INST_0][FEA_SOL]->SetConvergence(false);
      break;

    case MAIN_SOLVER::ADJ_EULER: case MAIN_SOLVER::ADJ_NAVIER_STOKES: case MAIN_SOLVER::ADJ_RANS: case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
    case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
      integration_container[iZone][INST_0][ADJFLOW_SOL]->SetConvergence(false);
      if( (config_container[iZone]->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) || (config_container[iZone]->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_RANS) )
        integration_container[iZone][INST_0][ADJTURB_SOL]->SetConvergence(false);
      break;

    default:
      break;
    }
  }

}

void CSinglezoneDriver::SetInitialMesh() {

  DynamicMeshUpdate(0);

  SU2_OMP_PARALLEL {
    // Overwrite fictious velocities
    for (iMesh = 0u; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++) {
      SU2_OMP_FOR_STAT(roundUpDiv(geometry_container[ZONE_0][INST_0][iMesh]->GetnPoint(),omp_get_max_threads()))
      for (unsigned long iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][iMesh]->GetnPoint(); iPoint++) {

        /*--- Overwrite fictitious velocities ---*/
        su2double Grid_Vel[3] = {0.0, 0.0, 0.0};

        /*--- Set the grid velocity for this coarse node. ---*/
        geometry_container[ZONE_0][INST_0][iMesh]->nodes->SetGridVel(iPoint, Grid_Vel);
      }
      END_SU2_OMP_FOR
      /*--- Push back the volume. ---*/
      geometry_container[ZONE_0][INST_0][iMesh]->nodes->SetVolume_n();
      geometry_container[ZONE_0][INST_0][iMesh]->nodes->SetVolume_nM1();
    }
    /*--- Push back the solution so that there is no fictious velocity at the next step. ---*/
    solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->Set_Solution_time_n();
    solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->Set_Solution_time_n1();
  }
  END_SU2_OMP_PARALLEL
}

void CDriver::BoundaryConditionsUpdate(){

  int rank = MASTER_NODE;
  unsigned short iZone;

  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);

  if(rank == MASTER_NODE) cout << "Updating boundary conditions." << endl;
  for(iZone = 0; iZone < nZone; iZone++){
    geometry_container[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry_container[iZone][INST_0], config_container[iZone]);
  }
}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to finite elements                                       */
////////////////////////////////////////////////////////////////////////////////

void CDriver::SetFEA_Loads(unsigned short iMarker, unsigned long iVertex, passivedouble LoadX,
                       passivedouble LoadY, passivedouble LoadZ) {

  unsigned long iPoint;
  su2double NodalForce[3] = {0.0,0.0,0.0};
  NodalForce[0] = LoadX;
  NodalForce[1] = LoadY;
  NodalForce[2] = LoadZ;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->GetNodes()->Set_FlowTraction(iPoint,NodalForce);

}

vector<passivedouble> CDriver::GetFEA_Displacements(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  vector<su2double> Displacements(3, 0.0);
  vector<passivedouble> Displacements_passive(3, 0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  Displacements[0] = solver->GetNodes()->GetSolution(iPoint, 0);
  Displacements[1] = solver->GetNodes()->GetSolution(iPoint, 1);
  if (geometry->GetnDim() == 3)
    Displacements[2] = solver->GetNodes()->GetSolution(iPoint, 2);
  else
    Displacements[2] = 0.0;

  Displacements_passive[0] = SU2_TYPE::GetValue(Displacements[0]);
  Displacements_passive[1] = SU2_TYPE::GetValue(Displacements[1]);
  Displacements_passive[2] = SU2_TYPE::GetValue(Displacements[2]);

  return Displacements_passive;
}


vector<passivedouble> CDriver::GetFEA_Velocity(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  vector<su2double> Velocity(3, 0.0);
  vector<passivedouble> Velocity_passive(3,0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  if (config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC){
    Velocity[0] = solver->GetNodes()->GetSolution_Vel(iPoint, 0);
    Velocity[1] = solver->GetNodes()->GetSolution_Vel(iPoint, 1);
    if (geometry->GetnDim() == 3)
      Velocity[2] = solver->GetNodes()->GetSolution_Vel(iPoint, 2);
    else
      Velocity[2] = 0.0;
  }

  Velocity_passive[0] = SU2_TYPE::GetValue(Velocity[0]);
  Velocity_passive[1] = SU2_TYPE::GetValue(Velocity[1]);
  Velocity_passive[2] = SU2_TYPE::GetValue(Velocity[2]);

  return Velocity_passive;
}

vector<passivedouble> CDriver::GetFEA_Velocity_n(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  vector<su2double> Velocity_n(3, 0.0);
  vector<passivedouble> Velocity_n_passive(3, 0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  if (config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC){
    Velocity_n[0] = solver->GetNodes()->GetSolution_Vel_time_n(iPoint, 0);
    Velocity_n[1] = solver->GetNodes()->GetSolution_Vel_time_n(iPoint, 1);
    if (geometry->GetnDim() == 3)
      Velocity_n[2] = solver->GetNodes()->GetSolution_Vel_time_n(iPoint, 2);
    else
      Velocity_n[2] = 0.0;
  }

  Velocity_n_passive[0] = SU2_TYPE::GetValue(Velocity_n[0]);
  Velocity_n_passive[1] = SU2_TYPE::GetValue(Velocity_n[1]);
  Velocity_n_passive[2] = SU2_TYPE::GetValue(Velocity_n[2]);

  return Velocity_n_passive;

}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to adjoint simulations                                   */
////////////////////////////////////////////////////////////////////////////////

vector<passivedouble> CDriver::GetMeshDisp_Sensitivity(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  vector<su2double> Disp_Sens(3, 0.0);
  vector<passivedouble> Disp_Sens_passive(3, 0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  CSolver *solver =  solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  Disp_Sens[0] = solver->GetNodes()->GetBoundDisp_Sens(iPoint, 0);
  Disp_Sens[1] = solver->GetNodes()->GetBoundDisp_Sens(iPoint, 1);
  if (geometry->GetnDim() == 3)
    Disp_Sens[2] = solver->GetNodes()->GetBoundDisp_Sens(iPoint, 2);
  else
    Disp_Sens[2] = 0.0;

  Disp_Sens_passive[0] = SU2_TYPE::GetValue(Disp_Sens[0]);
  Disp_Sens_passive[1] = SU2_TYPE::GetValue(Disp_Sens[1]);
  Disp_Sens_passive[2] = SU2_TYPE::GetValue(Disp_Sens[2]);

  return Disp_Sens_passive;

}

vector<passivedouble> CDriver::GetFlowLoad_Sensitivity(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  vector<su2double> FlowLoad_Sens(3, 0.0);
  vector<passivedouble> FlowLoad_Sens_passive(3, 0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  FlowLoad_Sens[0] = solver->GetNodes()->GetFlowTractionSensitivity(iPoint, 0);
  FlowLoad_Sens[1] = solver->GetNodes()->GetFlowTractionSensitivity(iPoint, 1);
  if (geometry->GetnDim() == 3)
    FlowLoad_Sens[2] = solver->GetNodes()->GetFlowTractionSensitivity(iPoint, 2);
  else
    FlowLoad_Sens[2] = 0.0;

  FlowLoad_Sens_passive[0] = SU2_TYPE::GetValue(FlowLoad_Sens[0]);
  FlowLoad_Sens_passive[1] = SU2_TYPE::GetValue(FlowLoad_Sens[1]);
  FlowLoad_Sens_passive[2] = SU2_TYPE::GetValue(FlowLoad_Sens[2]);

  return FlowLoad_Sens_passive;

}

void CDriver::SetFlowLoad_Adjoint(unsigned short iMarker, unsigned long iVertex, passivedouble val_AdjointX,
                                  passivedouble val_AdjointY, passivedouble val_AdjointZ) {

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  solver->StoreVertexTractionsAdjoint(iMarker, iVertex, 0, val_AdjointX);
  solver->StoreVertexTractionsAdjoint(iMarker, iVertex, 1, val_AdjointY);
  if (geometry->GetnDim() == 3)
    solver->StoreVertexTractionsAdjoint(iMarker, iVertex, 2, val_AdjointZ);

}

void CDriver::SetSourceTerm_DispAdjoint(unsigned short iMarker, unsigned long iVertex, passivedouble val_AdjointX,
                                        passivedouble val_AdjointY, passivedouble val_AdjointZ) {

  unsigned long iPoint;

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];
  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  solver->GetNodes()->SetSourceTerm_DispAdjoint(iPoint, 0, val_AdjointX);
  solver->GetNodes()->SetSourceTerm_DispAdjoint(iPoint, 1, val_AdjointY);
  if (geometry->GetnDim() == 3)
    solver->GetNodes()->SetSourceTerm_DispAdjoint(iPoint, 2, val_AdjointZ);

}

void CDriver::SetSourceTerm_VelAdjoint(unsigned short iMarker, unsigned long iVertex, passivedouble val_AdjointX,
                                        passivedouble val_AdjointY, passivedouble val_AdjointZ) {

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];
  const auto iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  solver->GetNodes()->SetSourceTerm_VelAdjoint(iPoint, 0, val_AdjointX);
  solver->GetNodes()->SetSourceTerm_VelAdjoint(iPoint, 1, val_AdjointY);
  if (geometry->GetnDim() == 3)
    solver->GetNodes()->SetSourceTerm_VelAdjoint(iPoint, 2, val_AdjointZ);

}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to mesh deformation */
////////////////////////////////////////////////////////////////////////////////

void CDriver::SetMeshDisplacement(unsigned short iMarker, unsigned long iVertex, passivedouble DispX, passivedouble DispY, passivedouble DispZ) {

  unsigned long iPoint;
  su2double MeshDispl[3] =  {0.0,0.0,0.0};

  MeshDispl[0] = DispX;
  MeshDispl[1] = DispY;
  MeshDispl[2] = DispZ;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->GetNodes()->SetBound_Disp(iPoint,MeshDispl);

}

void CDriver::CommunicateMeshDisplacement(void) {

  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->InitiateComms(geometry_container[ZONE_0][INST_0][MESH_0],
                                                                    config_container[ZONE_0], MESH_DISPLACEMENTS);
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->CompleteComms(geometry_container[ZONE_0][INST_0][MESH_0],
                                                                    config_container[ZONE_0], MESH_DISPLACEMENTS);

}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to flow loads */
////////////////////////////////////////////////////////////////////////////////

vector<passivedouble> CDriver::GetFlowLoad(unsigned short iMarker, unsigned long iVertex) const {

  vector<su2double> FlowLoad(3, 0.0);
  vector<passivedouble> FlowLoad_passive(3, 0.0);

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  if (config_container[ZONE_0]->GetSolid_Wall(iMarker)) {
    FlowLoad[0] = solver->GetVertexTractions(iMarker, iVertex, 0);
    FlowLoad[1] = solver->GetVertexTractions(iMarker, iVertex, 1);
    if (geometry->GetnDim() == 3)
      FlowLoad[2] = solver->GetVertexTractions(iMarker, iVertex, 2);
    else
      FlowLoad[2] = 0.0;
  }

  FlowLoad_passive[0] = SU2_TYPE::GetValue(FlowLoad[0]);
  FlowLoad_passive[1] = SU2_TYPE::GetValue(FlowLoad[1]);
  FlowLoad_passive[2] = SU2_TYPE::GetValue(FlowLoad[2]);

  return FlowLoad_passive;

}
