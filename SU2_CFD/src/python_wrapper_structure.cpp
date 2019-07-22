/*!
 * \file python_wrapper_structure.cpp
 * \brief Driver subroutines that are used by the Python wrapper. Those routines are usually called from an external Python environment.
 * \author D. Thomas
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

 #include "../include/driver_structure.hpp"

void CDriver::PythonInterface_Preprocessing(){

  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /* --- Initialize boundary conditions customization, this is achieve through the Python wrapper --- */
  for(iZone=0; iZone < nZone; iZone++){
    if (rank == MASTER_NODE) cout << "Setting customized boundary conditions for zone " << iZone << endl;
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
      geometry_container[iZone][INST_0][iMesh]->SetCustomBoundary(config_container[iZone]);
    }
    geometry_container[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry_container[iZone][INST_0], config_container[iZone]);

    if ((config_container[iZone]->GetKind_Solver() == EULER) ||
        (config_container[iZone]->GetKind_Solver() == NAVIER_STOKES) ||
        (config_container[iZone]->GetKind_Solver() == RANS)) {

          solver_container[iZone][INST_0][MESH_0][FLOW_SOL]->UpdateCustomBoundaryConditions(geometry_container[iZone][INST_0], config_container[iZone]);
    }
  }
  /*--- Initialize some variables used for external communications trough the Py wrapper. ---*/
  PyWrapVarCoord[0] = 0.0;
  PyWrapVarCoord[1] = 0.0;
  PyWrapVarCoord[2] = 0.0;
  PyWrapNodalForce[0] = 0.0;
  PyWrapNodalForce[1] = 0.0;
  PyWrapNodalForce[2] = 0.0;
  PyWrapNodalForceDensity[0] = 0.0;
  PyWrapNodalForceDensity[1] = 0.0;
  PyWrapNodalForceDensity[2] = 0.0;
  PyWrapNodalHeatFlux[0] = 0.0;
  PyWrapNodalHeatFlux[1] = 0.0;
  PyWrapNodalHeatFlux[2] = 0.0;

}

su2double CDriver::Get_Drag() {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CDrag, RefDensity, RefArea, RefVel2, factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefArea = config_container[val_iZone]->GetRefArea();

  /*--- Calculate free-stream velocity (squared) ---*/
  RefVel2 = 0.0;
  for(unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate drag force based on drag coefficient ---*/
  factor = 0.5*RefDensity*RefArea*RefVel2;
  CDrag = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CD();

  return CDrag*factor;
}

su2double CDriver::Get_Lift() {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CLift, RefDensity, RefArea, RefVel2, factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefArea = config_container[val_iZone]->GetRefArea();

  /*--- Calculate free-stream velocity (squared) ---*/
  RefVel2 = 0.0;
  for(unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate drag force based on drag coefficient ---*/
  factor = 0.5*RefDensity*RefArea*RefVel2;
  CLift = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CL();

  return CLift*factor;
}

su2double CDriver::Get_Mx(){

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMx, RefDensity, RefArea, RefLengthCoeff, RefVel2, factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefArea = config_container[val_iZone]->GetRefArea();
  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate free-stream velocity (squared) ---*/
  RefVel2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate moment around x-axis based on coefficients ---*/
  factor = 0.5*RefDensity*RefArea*RefVel2;
  CMx = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMx();

  return CMx*factor*RefLengthCoeff;

}

su2double CDriver::Get_My(){

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMy, RefDensity, RefArea, RefLengthCoeff, RefVel2, factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefArea = config_container[val_iZone]->GetRefArea();
  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate free-stream velocity (squared) ---*/
  RefVel2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate moment around x-axis based on coefficients ---*/
  factor = 0.5*RefDensity*RefArea*RefVel2;
  CMy = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMy();

  return CMy*factor*RefLengthCoeff;

}

su2double CDriver::Get_Mz() {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMz, RefDensity, RefArea, RefLengthCoeff, RefVel2, factor;

  /*--- Export free-stream density and reference area ---*/
  RefDensity = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetDensity_Inf();
  RefArea = config_container[val_iZone]->GetRefArea();
  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate free-stream velocity (squared) ---*/
  RefVel2 = 0.0;
  for(unsigned short iDim = 0; iDim < nDim; iDim++)
    RefVel2 += pow(solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim),2);

  /*--- Calculate moment around z-axis based on coefficients ---*/
  factor = 0.5*RefDensity*RefArea*RefVel2;
  CMz = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMz();

  return CMz*factor*RefLengthCoeff;

}

su2double CDriver::Get_DragCoeff() {

    unsigned short val_iZone = ZONE_0;
    unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
    su2double CDrag;

    CDrag = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CD();

    return CDrag;
}

su2double CDriver::Get_LiftCoeff() {

    unsigned short val_iZone = ZONE_0;
    unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
    su2double CLift;

    CLift = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CL();

    return CLift;
}

unsigned short CDriver::GetMovingMarker() {

  unsigned short IDtoSend,iMarker, jMarker, Moving;
  string Marker_Tag, Moving_Tag;

  IDtoSend = 0;
  for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
    Moving = config_container[ZONE_0]->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config_container[ZONE_0]->GetnMarker_Moving(); jMarker++) {
        Moving_Tag = config_container[ZONE_0]->GetMarker_Moving_TagBound(jMarker);
        Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Moving_Tag) {
          IDtoSend = iMarker;
          break;
        }
      }
    }
  }

  return IDtoSend;

}

unsigned long CDriver::GetNumberVertices(unsigned short iMarker){

  return geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker];

}

unsigned long CDriver::GetNumberHaloVertices(unsigned short iMarker){

  unsigned long nHaloVertices, iVertex, iPoint;

  nHaloVertices = 0;
  for(iVertex = 0; iVertex < geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker]; iVertex++){
    iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
    if(!(geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetDomain())) nHaloVertices += 1;
  }

  return nHaloVertices;

}

unsigned long CDriver::GetVertexGlobalIndex(unsigned short iMarker, unsigned short iVertex) {

  unsigned long iPoint, GlobalIndex;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  GlobalIndex = geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetGlobalIndex();

  return GlobalIndex;

}

bool CDriver::IsAHaloNode(unsigned short iMarker, unsigned short iVertex) {

  unsigned long iPoint;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  if(geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetDomain()) return false;
  else return true;

}

unsigned long CDriver::GetnExtIter() {

    return config_container[ZONE_0]->GetnExtIter();
}

unsigned long CDriver::GetExtIter(){

  return ExtIter;
}

su2double CDriver::GetUnsteady_TimeStep(){

  return config_container[ZONE_0]->GetDelta_UnstTime();
}

su2double CDriver::GetVertexCoordX(unsigned short iMarker, unsigned short iVertex) {

  su2double* Coord;
  unsigned long iPoint;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetCoord();
  return Coord[0];

}

su2double CDriver::GetVertexCoordY(unsigned short iMarker, unsigned short iVertex) {

  su2double* Coord;
  unsigned long iPoint;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetCoord();
  return Coord[1];
}

su2double CDriver::GetVertexCoordZ(unsigned short iMarker, unsigned short iVertex) {

  su2double* Coord;
  unsigned long iPoint;

  if(nDim == 3) {
    iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
    Coord = geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetCoord();
    return Coord[2];
  }
  else {
    return 0.0;
  }

}

bool CDriver::ComputeVertexForces(unsigned short iMarker, unsigned short iVertex) {

  unsigned long iPoint;
  unsigned short iDim, jDim;
  su2double *Normal, AreaSquare, Area;
  bool halo;

  unsigned short FinestMesh = config_container[ZONE_0]->GetFinestMesh();

  /*--- Check the kind of fluid problem ---*/
  bool compressible       = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (config_container[ZONE_0]->GetKind_Regime() == INCOMPRESSIBLE);
  bool viscous_flow       = ((config_container[ZONE_0]->GetKind_Solver() == NAVIER_STOKES) ||
                 (config_container[ZONE_0]->GetKind_Solver() == RANS) );

  /*--- Parameters for the calculations ---*/
  // Pn: Pressure
  // Pinf: Pressure_infinite
  // div_vel: Velocity divergence
  // Dij: Dirac delta
  su2double Pn = 0.0, div_vel = 0.0, Dij = 0.0;
  su2double Viscosity = 0.0;
  su2double Grad_Vel[3][3] = { {0.0, 0.0, 0.0} ,
              {0.0, 0.0, 0.0} ,
              {0.0, 0.0, 0.0} } ;
  su2double Tau[3][3] = { {0.0, 0.0, 0.0} ,
              {0.0, 0.0, 0.0} ,
              {0.0, 0.0, 0.0} } ;

  su2double Pinf = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->GetPressure_Inf();

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  /*--- It is necessary to distinguish the halo nodes from the others, since they introduce non physical forces. ---*/
  if(geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetDomain()) {
    /*--- Get the normal at the vertex: this normal goes inside the fluid domain. ---*/
    Normal = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();
    AreaSquare = 0.0;
    for(iDim = 0; iDim < nDim; iDim++) {
      AreaSquare += Normal[iDim]*Normal[iDim];
    }
    Area = sqrt(AreaSquare);

    /*--- Get the values of pressure and viscosity ---*/
    Pn = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetPressure();
    if (viscous_flow) {
      for(iDim=0; iDim<nDim; iDim++) {
        for(jDim=0; jDim<nDim; jDim++) {
          Grad_Vel[iDim][jDim] = solver_container[ZONE_0][INST_0][FinestMesh][FLOW_SOL]->node[iPoint]->GetGradient_Primitive(iDim+1, jDim);
        }
      }
      Viscosity = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    }

    /*--- Calculate the inviscid (pressure) part of tn in the fluid nodes (force units) ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
     PyWrapNodalForce[iDim] = -(Pn-Pinf)*Normal[iDim];     //NB : norm(Normal) = Area
    }

    /*--- Calculate the viscous (shear stress) part of tn in the fluid nodes (force units) ---*/
    if ((incompressible || compressible) && viscous_flow) {
      div_vel = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        div_vel += Grad_Vel[iDim][iDim];
     if (incompressible) div_vel = 0.0;

      for (iDim = 0; iDim < nDim; iDim++) {
       for (jDim = 0 ; jDim < nDim; jDim++) {
         Dij = 0.0; if (iDim == jDim) Dij = 1.0;
         Tau[iDim][jDim] = Viscosity*(Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim]) - TWO3*Viscosity*div_vel*Dij;
         PyWrapNodalForce[iDim] += Tau[iDim][jDim]*Normal[jDim];
        }
      }
    }

    //Divide by local are in case of force density communication.
   for(iDim = 0; iDim < nDim; iDim++) {
     PyWrapNodalForceDensity[iDim] = PyWrapNodalForce[iDim]/Area;
    }

    halo = false;
  }
  else {
    halo = true;
  }

  return halo;

}

su2double CDriver::GetVertexForceX(unsigned short iMarker, unsigned short iVertex) {

  return PyWrapNodalForce[0];

}

su2double CDriver::GetVertexForceY(unsigned short iMarker, unsigned short iVertex) {

  return PyWrapNodalForce[1];

}

su2double CDriver::GetVertexForceZ(unsigned short iMarker, unsigned short iVertex) {

  return PyWrapNodalForce[2];

}

su2double CDriver::GetVertexForceDensityX(unsigned short iMarker, unsigned short iVertex) {
  return PyWrapNodalForceDensity[0];
}

su2double CDriver::GetVertexForceDensityY(unsigned short iMarker, unsigned short iVertex) {
  return PyWrapNodalForceDensity[1];
}

su2double CDriver::GetVertexForceDensityZ(unsigned short iMarker, unsigned short iVertex) {
  return PyWrapNodalForceDensity[2];
}

void CDriver::SetVertexCoordX(unsigned short iMarker, unsigned short iVertex, su2double newPosX) {

  unsigned long iPoint;
  su2double *Coord;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetCoord();

  PyWrapVarCoord[0] = newPosX - Coord[0];

}

void CDriver::SetVertexCoordY(unsigned short iMarker, unsigned short iVertex, su2double newPosY) {

  unsigned long iPoint;
  su2double *Coord;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetCoord();

  PyWrapVarCoord[1] = newPosY - Coord[1];
}

void CDriver::SetVertexCoordZ(unsigned short iMarker, unsigned short iVertex, su2double newPosZ) {

  unsigned long iPoint;
  su2double *Coord;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetCoord();

  if(nDim > 2) {
    PyWrapVarCoord[2] = newPosZ - Coord[2];
  }
  else {
    PyWrapVarCoord[2] = 0.0;
  }
}

su2double CDriver::SetVertexVarCoord(unsigned short iMarker, unsigned short iVertex) {

  su2double nodalVarCoordNorm;

    geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(PyWrapVarCoord);
    nodalVarCoordNorm = sqrt((PyWrapVarCoord[0])*(PyWrapVarCoord[0]) + (PyWrapVarCoord[1])*(PyWrapVarCoord[1]) + (PyWrapVarCoord[2])*(PyWrapVarCoord[2]));

  return nodalVarCoordNorm;

}

su2double CDriver::GetVertexTemperature(unsigned short iMarker, unsigned short iVertex){

  unsigned long iPoint;
  su2double vertexWallTemp(0.0);

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetDomain() && compressible){
    vertexWallTemp = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetTemperature();
  }

  return vertexWallTemp;

}

void CDriver::SetVertexTemperature(unsigned short iMarker, unsigned short iVertex, su2double val_WallTemp){

  geometry_container[ZONE_0][INST_0][MESH_0]->SetCustomBoundaryTemperature(iMarker, iVertex, val_WallTemp);
}

bool CDriver::ComputeVertexHeatFluxes(unsigned short iMarker, unsigned short iVertex){

  unsigned long iPoint;
  unsigned short iDim;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double laminar_viscosity, thermal_conductivity;
  su2double GradT[3] = {0.0,0.0,0.0};

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);
  bool halo;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetDomain()){
    halo = false;
  }
  else{
    halo = true;
  }

  if(!halo && compressible){
    laminar_viscosity    = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);
    for(iDim=0; iDim < nDim; iDim++){
      GradT[iDim] = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetGradient_Primitive(0, iDim);
      PyWrapNodalHeatFlux[iDim] = -thermal_conductivity*GradT[iDim];
    }
  }

  return halo;
}

su2double CDriver::GetVertexHeatFluxX(unsigned short iMarker, unsigned short iVertex){

  return PyWrapNodalHeatFlux[0];
}

su2double CDriver::GetVertexHeatFluxY(unsigned short iMarker, unsigned short iVertex){

  return PyWrapNodalHeatFlux[1];
}

su2double CDriver::GetVertexHeatFluxZ(unsigned short iMarker, unsigned short iVertex){

  return PyWrapNodalHeatFlux[2];
}

su2double CDriver::GetVertexNormalHeatFlux(unsigned short iMarker, unsigned short iVertex){

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

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);

  vertexWallHeatFlux = 0.0;
  dTdn = 0.0;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetDomain() && compressible){
    Normal = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();
    Area = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Area += Normal[iDim]*Normal[iDim];
    Area = sqrt(Area);

    for (iDim = 0; iDim < nDim; iDim++)
      UnitNormal[iDim] = Normal[iDim]/Area;

    laminar_viscosity    = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);
    /*Compute wall heat flux (normal to the wall) based on computed temperature gradient*/
    for(iDim=0; iDim < nDim; iDim++){
      GradT[iDim] = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetGradient_Primitive(0, iDim);
      dTdn += GradT[iDim]*UnitNormal[iDim];
    }

    vertexWallHeatFlux = -thermal_conductivity*dTdn;
  }

  return vertexWallHeatFlux;
}

void CDriver::SetVertexNormalHeatFlux(unsigned short iMarker, unsigned short iVertex, su2double val_WallHeatFlux){

  geometry_container[ZONE_0][INST_0][MESH_0]->SetCustomBoundaryHeatFlux(iMarker, iVertex, val_WallHeatFlux);
}

su2double CDriver::GetThermalConductivity(unsigned short iMarker, unsigned short iVertex){

  unsigned long iPoint;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double laminar_viscosity, thermal_conductivity;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  laminar_viscosity    = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
  thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);

  return thermal_conductivity;

}

vector<su2double> CDriver::GetVertexUnitNormal(unsigned short iMarker, unsigned short iVertex){

  unsigned short iDim;
  su2double *Normal;
  su2double Area;
  vector<su2double> ret_Normal(3, 0.0);

  Normal = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
      Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  ret_Normal[0] = Normal[0]/Area;
  ret_Normal[1] = Normal[1]/Area;
  if(nDim>2) ret_Normal[2] = Normal[2]/Area;

  return ret_Normal;


}

vector<string> CDriver::GetAllBoundaryMarkersTag(){

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

vector<string> CDriver::GetAllMovingMarkersTag(){

  vector<string> movingBoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_Moving();
  movingBoundariesTagList.resize(nBoundariesMarker);

  for(iMarker=0; iMarker < nBoundariesMarker; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_Moving_TagBound(iMarker);
    movingBoundariesTagList[iMarker] = Marker_Tag;
  }

  return movingBoundariesTagList;
}

vector<string> CDriver::GetAllCHTMarkersTag(){

  vector<string> CHTBoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_All();
  //CHTBoundariesTagList.resize(nBoundariesMarker);

  //The CHT markers can be identified as the markers that are customizable with a BC type HEAT_FLUX or ISOTHERMAL.
  for(iMarker=0; iMarker<nBoundariesMarker; iMarker++){
    if((config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == HEAT_FLUX || config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) && config_container[ZONE_0]->GetMarker_All_PyCustom(iMarker)){
      Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
      CHTBoundariesTagList.push_back(Marker_Tag);
    }
  }

  return CHTBoundariesTagList;
}

vector<string> CDriver::GetAllInletMarkersTag(){

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

map<string, int> CDriver::GetAllBoundaryMarkers(){

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

map<string, string> CDriver::GetAllBoundaryMarkersType(){

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

void CFluidDriver::ResetConvergence() {

  for(iZone = 0; iZone < nZone; iZone++) {
    switch (config_container[iZone]->GetKind_Solver()) {

    case EULER: case NAVIER_STOKES: case RANS:
      integration_container[iZone][INST_0][FLOW_SOL]->SetConvergence(false);
      if (config_container[iZone]->GetKind_Solver() == RANS) integration_container[iZone][INST_0][TURB_SOL]->SetConvergence(false);
      if(config_container[iZone]->GetKind_Trans_Model() == LM) integration_container[iZone][INST_0][TRANS_SOL]->SetConvergence(false);
      break;

    case FEM_ELASTICITY:
      integration_container[iZone][INST_0][FEA_SOL]->SetConvergence(false);
      break;

    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS: case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
      integration_container[iZone][INST_0][ADJFLOW_SOL]->SetConvergence(false);
      if( (config_container[iZone]->GetKind_Solver() == ADJ_RANS) || (config_container[iZone]->GetKind_Solver() == DISC_ADJ_RANS) )
        integration_container[iZone][INST_0][ADJTURB_SOL]->SetConvergence(false);
      break;
    }
  }

}

void CFluidDriver::StaticMeshUpdate() {

  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  for(iZone = 0; iZone < nZone; iZone++) {
    if(rank == MASTER_NODE) cout << " Deforming the volume grid." << endl;
    grid_movement[iZone][INST_0]->SetVolume_Deformation(geometry_container[iZone][INST_0][MESH_0], config_container[iZone], true);

    if(rank == MASTER_NODE) cout << "No grid velocity to be computde : static grid deformation." << endl;

    if(rank == MASTER_NODE) cout << " Updating multigrid structure." << endl;
    grid_movement[iZone][INST_0]->UpdateMultiGrid(geometry_container[iZone][INST_0], config_container[iZone]);
  }
}

void CFluidDriver::SetInitialMesh() {

  unsigned long iPoint;

  StaticMeshUpdate();

  /*--- Propagate the initial deformation to the past ---*/
  //if (!restart) {
    for(iZone = 0; iZone < nZone; iZone++) {
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++) {
        for(iPoint = 0; iPoint < geometry_container[iZone][INST_0][iMesh]->GetnPoint(); iPoint++) {
        //solver_container[iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
        //solver_container[iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
        geometry_container[iZone][INST_0][iMesh]->node[iPoint]->SetVolume_n();
        geometry_container[iZone][INST_0][iMesh]->node[iPoint]->SetVolume_nM1();
        geometry_container[iZone][INST_0][iMesh]->node[iPoint]->SetCoord_n();
        geometry_container[iZone][INST_0][iMesh]->node[iPoint]->SetCoord_n1();
      }
    }
  }
  //}
}

void CFluidDriver::SetVertexTtotal(unsigned short iMarker, unsigned short iVertex, su2double val_Ttotal){

  solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInlet_Ttotal(iMarker, iVertex, val_Ttotal);

}

void CFluidDriver::SetVertexPtotal(unsigned short iMarker, unsigned short iVertex, su2double val_Ptotal){

  solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInlet_Ptotal(iMarker, iVertex, val_Ptotal);

}

void CFluidDriver::SetVertexFlowDir(unsigned short iMarker, unsigned short iVertex, unsigned short iDim, su2double val_FlowDir){

  solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInlet_FlowDir(iMarker, iVertex, iDim, val_FlowDir);

}

void CFluidDriver::SetVertexTurbVar(unsigned short iMarker, unsigned short iVertex, unsigned short iDim, su2double val_turb_var){

  if (solver_container[ZONE_0][INST_0] == NULL ||
      solver_container[ZONE_0][INST_0][MESH_0] ==  NULL) {
    SU2_MPI::Error("Could not find an appropriate solver.", CURRENT_FUNCTION);
  } else if (solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL] == NULL) {
    SU2_MPI::Error("Tried to set turbulence variables without a turbulence solver.", CURRENT_FUNCTION);
  }
  solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->SetInlet_TurbVar(iMarker, iVertex, iDim, val_turb_var);

}

void CFluidDriver::BoundaryConditionsUpdate(){

  int rank = MASTER_NODE;
  unsigned short iZone;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if(rank == MASTER_NODE) cout << "Updating boundary conditions." << endl;
  for(iZone = 0; iZone < nZone; iZone++){
    geometry_container[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry_container[iZone][INST_0], config_container[iZone]);
  }
}
