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

 #include "../include/drivers/CDiscAdjSinglezoneDriver.hpp"

void CDriver::PythonInterface_Preprocessing(CConfig **config, CGeometry ****geometry, CSolver *****solver){

  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /* --- Initialize boundary conditions customization, this is achieve through the Python wrapper --- */
  for(iZone=0; iZone < nZone; iZone++){
    
    if (config[iZone]->GetnMarker_PyCustom() > 0){
      
      if (rank == MASTER_NODE) cout << endl << "----------------- Python Interface Preprocessing ( Zone "<< iZone <<" ) -----------------" << endl;
      
      if (rank == MASTER_NODE) cout << "Setting customized boundary conditions for zone " << iZone << endl;
      for (iMesh = 0; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
        geometry[iZone][INST_0][iMesh]->SetCustomBoundary(config[iZone]);
      }
      geometry[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry[iZone][INST_0], config[iZone]);
      
      if ((config[iZone]->GetKind_Solver() == EULER) ||
          (config[iZone]->GetKind_Solver() == NAVIER_STOKES) ||
          (config[iZone]->GetKind_Solver() == RANS)) {
        
        solver[iZone][INST_0][MESH_0][FLOW_SOL]->UpdateCustomBoundaryConditions(geometry[iZone][INST_0], config[iZone]);
      }
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

passivedouble CDriver::Get_Drag() {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CDrag, RefDensity, RefArea, RefVel2, factor, val_Drag;

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

  val_Drag = CDrag*factor;

  return SU2_TYPE::GetValue(val_Drag);
}

passivedouble CDriver::Get_Lift() {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CLift, RefDensity, RefArea, RefVel2, factor, val_Lift;

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

  val_Lift = CLift*factor;

  return SU2_TYPE::GetValue(val_Lift);
}

passivedouble CDriver::Get_Mx(){

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMx, RefDensity, RefArea, RefLengthCoeff, RefVel2, factor, val_Mx;

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

  val_Mx = CMx*factor*RefLengthCoeff;

  return SU2_TYPE::GetValue(val_Mx);

}

passivedouble CDriver::Get_My(){

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMy, RefDensity, RefArea, RefLengthCoeff, RefVel2, factor, val_My;

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

  val_My = CMy*factor*RefLengthCoeff;

  return SU2_TYPE::GetValue(val_My);

}

passivedouble CDriver::Get_Mz() {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMz, RefDensity, RefArea, RefLengthCoeff, RefVel2, factor, val_Mz;

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

  val_Mz = CMz*factor*RefLengthCoeff;

  return SU2_TYPE::GetValue(val_Mz);

}

passivedouble CDriver::Get_DragCoeff() {

    unsigned short val_iZone = ZONE_0;
    unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
    su2double CDrag;

    CDrag = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CD();

    return SU2_TYPE::GetValue(CDrag);
}

passivedouble CDriver::Get_LiftCoeff() {

    unsigned short val_iZone = ZONE_0;
    unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
    su2double CLift;

    CLift = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CL();

    return SU2_TYPE::GetValue(CLift);
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

unsigned long CDriver::GetnIter() {

    return config_container[ZONE_0]->GetnIter();
}

unsigned long CDriver::GetIter(){

  return ExtIter;
}

passivedouble CDriver::GetUnsteady_TimeStep(){

  return SU2_TYPE::GetValue(config_container[ZONE_0]->GetDelta_UnstTime());
}

passivedouble CDriver::GetVertexCoordX(unsigned short iMarker, unsigned short iVertex) {

  su2double* Coord;
  unsigned long iPoint;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetCoord();
  return SU2_TYPE::GetValue(Coord[0]);

}

passivedouble CDriver::GetVertexCoordY(unsigned short iMarker, unsigned short iVertex) {

  su2double* Coord;
  unsigned long iPoint;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetCoord();
  return SU2_TYPE::GetValue(Coord[1]);
}

passivedouble CDriver::GetVertexCoordZ(unsigned short iMarker, unsigned short iVertex) {

  su2double* Coord;
  unsigned long iPoint;

  if(nDim == 3) {
    iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
    Coord = geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetCoord();
    return SU2_TYPE::GetValue(Coord[2]);
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

  su2double Pinf = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetPressure_Inf();

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

passivedouble CDriver::GetVertexForceX(unsigned short iMarker, unsigned short iVertex) {

  return SU2_TYPE::GetValue(PyWrapNodalForce[0]);

}

passivedouble CDriver::GetVertexForceY(unsigned short iMarker, unsigned short iVertex) {

  return SU2_TYPE::GetValue(PyWrapNodalForce[1]);

}

passivedouble CDriver::GetVertexForceZ(unsigned short iMarker, unsigned short iVertex) {

  return SU2_TYPE::GetValue(PyWrapNodalForce[2]);

}

passivedouble CDriver::GetVertexForceDensityX(unsigned short iMarker, unsigned short iVertex) {
  return SU2_TYPE::GetValue(PyWrapNodalForceDensity[0]);
}

passivedouble CDriver::GetVertexForceDensityY(unsigned short iMarker, unsigned short iVertex) {
  return SU2_TYPE::GetValue(PyWrapNodalForceDensity[1]);
}

passivedouble CDriver::GetVertexForceDensityZ(unsigned short iMarker, unsigned short iVertex) {
  return SU2_TYPE::GetValue(PyWrapNodalForceDensity[2]);
}

void CDriver::SetVertexCoordX(unsigned short iMarker, unsigned short iVertex, passivedouble newPosX) {

  unsigned long iPoint;
  su2double *Coord;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetCoord();

  PyWrapVarCoord[0] = newPosX - Coord[0];

}

void CDriver::SetVertexCoordY(unsigned short iMarker, unsigned short iVertex, passivedouble newPosY) {

  unsigned long iPoint;
  su2double *Coord;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  Coord = geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetCoord();

  PyWrapVarCoord[1] = newPosY - Coord[1];
}

void CDriver::SetVertexCoordZ(unsigned short iMarker, unsigned short iVertex, passivedouble newPosZ) {

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

passivedouble CDriver::SetVertexVarCoord(unsigned short iMarker, unsigned short iVertex) {

  su2double nodalVarCoordNorm;

    geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->SetVarCoord(PyWrapVarCoord);
    nodalVarCoordNorm = sqrt((PyWrapVarCoord[0])*(PyWrapVarCoord[0]) + (PyWrapVarCoord[1])*(PyWrapVarCoord[1]) + (PyWrapVarCoord[2])*(PyWrapVarCoord[2]));

  return SU2_TYPE::GetValue(nodalVarCoordNorm);

}

passivedouble CDriver::GetVertexTemperature(unsigned short iMarker, unsigned short iVertex){

  unsigned long iPoint;
  su2double vertexWallTemp(0.0);

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);
  
  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][INST_0][MESH_0]->node[iPoint]->GetDomain() && compressible){
    vertexWallTemp = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetTemperature();
  }

  return SU2_TYPE::GetValue(vertexWallTemp);

}

void CDriver::SetVertexTemperature(unsigned short iMarker, unsigned short iVertex, passivedouble val_WallTemp){

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

passivedouble CDriver::GetVertexHeatFluxX(unsigned short iMarker, unsigned short iVertex){

  return SU2_TYPE::GetValue(PyWrapNodalHeatFlux[0]);
}

passivedouble CDriver::GetVertexHeatFluxY(unsigned short iMarker, unsigned short iVertex){

  return SU2_TYPE::GetValue(PyWrapNodalHeatFlux[1]);
}

passivedouble CDriver::GetVertexHeatFluxZ(unsigned short iMarker, unsigned short iVertex){

  return SU2_TYPE::GetValue(PyWrapNodalHeatFlux[2]);
}

passivedouble CDriver::GetVertexNormalHeatFlux(unsigned short iMarker, unsigned short iVertex){

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

  return SU2_TYPE::GetValue(vertexWallHeatFlux);
}

void CDriver::SetVertexNormalHeatFlux(unsigned short iMarker, unsigned short iVertex, passivedouble val_WallHeatFlux){

  geometry_container[ZONE_0][INST_0][MESH_0]->SetCustomBoundaryHeatFlux(iMarker, iVertex, val_WallHeatFlux);
}

passivedouble CDriver::GetThermalConductivity(unsigned short iMarker, unsigned short iVertex){

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

  return SU2_TYPE::GetValue(thermal_conductivity);

}

vector<passivedouble> CDriver::GetVertexUnitNormal(unsigned short iMarker, unsigned short iVertex){

  unsigned short iDim;
  su2double *Normal;
  su2double Area;
  vector<su2double> ret_Normal(3, 0.0);
  vector<passivedouble> ret_Normal_passive(3, 0.0);

  Normal = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
      Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  ret_Normal[0] = Normal[0]/Area;
  ret_Normal[1] = Normal[1]/Area;
  if(nDim>2) ret_Normal[2] = Normal[2]/Area;

  ret_Normal_passive[0] = SU2_TYPE::GetValue(ret_Normal[0]);
  ret_Normal_passive[1] = SU2_TYPE::GetValue(ret_Normal[1]);
  ret_Normal_passive[2] = SU2_TYPE::GetValue(ret_Normal[2]);

  return ret_Normal_passive;


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

vector<string> CDriver::GetAllDeformMeshMarkersTag(){

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

vector<string> CDriver::GetAllFluidLoadMarkersTag(){

  vector<string> interfaceBoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_Fluid_Load();
  interfaceBoundariesTagList.resize(nBoundariesMarker);

  for(iMarker=0; iMarker < nBoundariesMarker; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_Fluid_Load_TagBound(iMarker);
    interfaceBoundariesTagList[iMarker] = Marker_Tag;
  }

  return interfaceBoundariesTagList;
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

void CDriver::ResetConvergence() {

  for(iZone = 0; iZone < nZone; iZone++) {
    switch (config_container[iZone]->GetKind_Solver()) {

    case EULER: case NAVIER_STOKES: case RANS:
    case INC_EULER: case INC_NAVIER_STOKES: case INC_RANS:
      integration_container[iZone][INST_0][FLOW_SOL]->SetConvergence(false);
      if (config_container[iZone]->GetKind_Solver() == RANS) integration_container[iZone][INST_0][TURB_SOL]->SetConvergence(false);
      if(config_container[iZone]->GetKind_Trans_Model() == LM) integration_container[iZone][INST_0][TRANS_SOL]->SetConvergence(false);
      break;

    case FEM_ELASTICITY:
      integration_container[iZone][INST_0][FEA_SOL]->SetConvergence(false);
      break;

    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS: case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
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

void CFluidDriver::SetVertexTtotal(unsigned short iMarker, unsigned short iVertex, passivedouble val_Ttotal_passive){

  su2double val_Ttotal = val_Ttotal_passive;

  solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInlet_Ttotal(iMarker, iVertex, val_Ttotal);

}

void CFluidDriver::SetVertexPtotal(unsigned short iMarker, unsigned short iVertex, passivedouble val_Ptotal_passive){

  su2double val_Ptotal = val_Ptotal_passive;

  solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInlet_Ptotal(iMarker, iVertex, val_Ptotal);

}

void CFluidDriver::SetVertexFlowDir(unsigned short iMarker, unsigned short iVertex, unsigned short iDim, passivedouble val_FlowDir_passive){

  su2double val_FlowDir = val_FlowDir_passive;

  solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInlet_FlowDir(iMarker, iVertex, iDim, val_FlowDir);

}

void CFluidDriver::SetVertexTurbVar(unsigned short iMarker, unsigned short iVertex, unsigned short iDim, passivedouble val_turb_var_passive){

  su2double val_turb_var = val_turb_var_passive;

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

void CDriver::SetMeshDisplacement(unsigned short iMarker, unsigned long iVertex, passivedouble DispX, passivedouble DispY, passivedouble DispZ) {

  unsigned long iPoint;
  PyWrapVarCoord[0] = DispX;
  PyWrapVarCoord[1] = DispY;
  PyWrapVarCoord[2] = DispZ;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->node[iPoint]->SetBound_Disp(PyWrapVarCoord);

}

void CDriver::CommunicateMeshDisplacement(void) {

  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->InitiateComms(geometry_container[ZONE_0][INST_0][MESH_0],
                                                                    config_container[ZONE_0], MESH_DISPLACEMENTS);
  solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL]->CompleteComms(geometry_container[ZONE_0][INST_0][MESH_0],
                                                                    config_container[ZONE_0], MESH_DISPLACEMENTS);

}

vector<passivedouble> CDriver::GetMeshDisp_Sensitivity(unsigned short iMarker, unsigned short iVertex) {

  unsigned long iPoint;
  vector<su2double> Disp_Sens(3, 0.0);
  vector<passivedouble> Disp_Sens_passive(3, 0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  Disp_Sens[0] = solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL]->node[iPoint]->GetBoundDisp_Sens(0);
  Disp_Sens[1] = solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL]->node[iPoint]->GetBoundDisp_Sens(1);
  if (solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL]->GetnVar() == 3)
    Disp_Sens[2] = solver_container[ZONE_0][INST_0][MESH_0][ADJMESH_SOL]->node[iPoint]->GetBoundDisp_Sens(2);
  else
    Disp_Sens[2] = 0.0;

  Disp_Sens_passive[0] = SU2_TYPE::GetValue(Disp_Sens[0]);
  Disp_Sens_passive[1] = SU2_TYPE::GetValue(Disp_Sens[1]);
  Disp_Sens_passive[2] = SU2_TYPE::GetValue(Disp_Sens[2]);

  return Disp_Sens_passive;

}

void CDriver::SetFEA_Loads(unsigned short iMarker, unsigned short iVertex, passivedouble LoadX,
                       passivedouble LoadY, passivedouble LoadZ) {

  unsigned long iPoint;
  PyWrapNodalForce[0] = LoadX;
  PyWrapNodalForce[1] = LoadY;
  PyWrapNodalForce[2] = LoadZ;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->node[iPoint]->Set_FlowTraction(PyWrapNodalForce);

}

vector<passivedouble> CDriver::GetFEA_Displacements(unsigned short iMarker, unsigned short iVertex) {

  unsigned long iPoint;
  vector<su2double> Displacements(3, 0.0);
  vector<passivedouble> Displacements_passive(3, 0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  Displacements[0] = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->node[iPoint]->GetSolution(0);
  Displacements[1] = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->node[iPoint]->GetSolution(1);
  if (solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->GetnVar() == 3)
    Displacements[2] = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->node[iPoint]->GetSolution(2);
  else
    Displacements[2] = 0.0;

  Displacements_passive[0] = SU2_TYPE::GetValue(Displacements[0]);
  Displacements_passive[1] = SU2_TYPE::GetValue(Displacements[1]);
  Displacements_passive[2] = SU2_TYPE::GetValue(Displacements[2]);

  return Displacements_passive;
}


vector<passivedouble> CDriver::GetFEA_Velocity(unsigned short iMarker, unsigned short iVertex) {

  unsigned long iPoint;
  vector<su2double> Velocity(3, 0.0);
  vector<passivedouble> Velocity_passive(3,0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if (config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC){
    Velocity[0] = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel(0);
    Velocity[1] = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel(1);
    if (solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->GetnVar() == 3)
      Velocity[2] = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel(2);
    else
      Velocity[2] = 0.0;
  }
  else{
    Velocity[0] = 0.0;
    Velocity[1] = 0.0;
    Velocity[2] = 0.0;
  }

  Velocity_passive[0] = SU2_TYPE::GetValue(Velocity[0]);
  Velocity_passive[1] = SU2_TYPE::GetValue(Velocity[1]);
  Velocity_passive[2] = SU2_TYPE::GetValue(Velocity[2]);

  return Velocity_passive;
}

vector<passivedouble> CDriver::GetFEA_Velocity_n(unsigned short iMarker, unsigned short iVertex) {

  unsigned long iPoint;
  vector<su2double> Velocity_n(3, 0.0);
  vector<passivedouble> Velocity_n_passive(3, 0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if (config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC){
    Velocity_n[0] = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n(0);
    Velocity_n[1] = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n(1);
    if (solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->GetnVar() == 3)
      Velocity_n[2] = solver_container[ZONE_0][INST_0][MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n(2);
    else
      Velocity_n[2] = 0.0;
  }
  else{
    Velocity_n[0] = 0.0;
    Velocity_n[1] = 0.0;
    Velocity_n[2] = 0.0;
  }

  Velocity_n_passive[0] = SU2_TYPE::GetValue(Velocity_n[0]);
  Velocity_n_passive[1] = SU2_TYPE::GetValue(Velocity_n[1]);
  Velocity_n_passive[2] = SU2_TYPE::GetValue(Velocity_n[2]);

  return Velocity_n_passive;

}

vector<passivedouble> CDriver::GetFlowLoad_Sensitivity(unsigned short iMarker, unsigned short iVertex) {

  unsigned long iPoint;
  vector<su2double> FlowLoad_Sens(3, 0.0);
  vector<passivedouble> FlowLoad_Sens_passive(3, 0.0);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  FlowLoad_Sens[0] = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL]->node[iPoint]->GetFlowTractionSensitivity(0);
  FlowLoad_Sens[1] = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL]->node[iPoint]->GetFlowTractionSensitivity(1);
  if (solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL]->GetnVar() == 3)
    FlowLoad_Sens[2] = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL]->node[iPoint]->GetFlowTractionSensitivity(2);
  else
    FlowLoad_Sens[2] = 0.0;

  FlowLoad_Sens_passive[0] = SU2_TYPE::GetValue(FlowLoad_Sens[0]);
  FlowLoad_Sens_passive[1] = SU2_TYPE::GetValue(FlowLoad_Sens[1]);
  FlowLoad_Sens_passive[2] = SU2_TYPE::GetValue(FlowLoad_Sens[2]);

  return FlowLoad_Sens_passive;

}

vector<passivedouble> CDriver::GetFlowLoad(unsigned short iMarker, unsigned short iVertex) {

  vector<su2double> FlowLoad(3, 0.0);
  vector<passivedouble> FlowLoad_passive(3, 0.0);

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];

  if (config_container[ZONE_0]->GetMarker_All_Fluid_Load(iMarker) == YES) {
    FlowLoad[0] = solver->GetVertexTractions(iMarker, iVertex, 0);
    FlowLoad[1] = solver->GetVertexTractions(iMarker, iVertex, 1);
    if (solver->GetnVar() == 3)
      FlowLoad[2] = solver->GetVertexTractions(iMarker, iVertex, 2);
    else
      FlowLoad[2] = 0.0;
  }
  else{
    FlowLoad[0] = 0.0;
    FlowLoad[1] = 0.0;
    FlowLoad[2] = 0.0;
  }

  FlowLoad_passive[0] = SU2_TYPE::GetValue(FlowLoad[0]);
  FlowLoad_passive[1] = SU2_TYPE::GetValue(FlowLoad[1]);
  FlowLoad_passive[2] = SU2_TYPE::GetValue(FlowLoad[2]);

  return FlowLoad_passive;

}

void CDriver::SetFlowLoad_Adjoint(unsigned short iMarker, unsigned short iVertex, passivedouble val_AdjointX,
                                  passivedouble val_AdjointY, passivedouble val_AdjointZ) {

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];

  solver->StoreVertexTractionsAdjoint(iMarker, iVertex, 0, val_AdjointX);
  solver->StoreVertexTractionsAdjoint(iMarker, iVertex, 1, val_AdjointY);
  if (solver->GetnVar() == 3)
    solver->StoreVertexTractionsAdjoint(iMarker, iVertex, 2, val_AdjointZ);

}

void CDriver::SetSourceTerm_DispAdjoint(unsigned short iMarker, unsigned short iVertex, passivedouble val_AdjointX,
                                        passivedouble val_AdjointY, passivedouble val_AdjointZ) {

  unsigned long iPoint;

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][ADJFEA_SOL];
  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  solver->node[iPoint]->SetSourceTerm_DispAdjoint(0, val_AdjointX);
  solver->node[iPoint]->SetSourceTerm_DispAdjoint(1, val_AdjointY);
  if (solver->GetnVar() == 3)
    solver->node[iPoint]->SetSourceTerm_DispAdjoint(2, val_AdjointZ);

}

vector<passivedouble> CDriver::GetVertex_UndeformedCoord(unsigned short iMarker, unsigned short iVertex) {

  unsigned long iPoint;
  vector<su2double> MeshCoord(3, 0.0);
  vector<passivedouble> MeshCoord_passive(3, 0.0);

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][MESH_SOL];

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if (solver != NULL) {
    MeshCoord[0] = solver->node[iPoint]->GetMesh_Coord(0);
    MeshCoord[1] = solver->node[iPoint]->GetMesh_Coord(1);
    if (solver->GetnVar() == 3)
      MeshCoord[2] = solver->node[iPoint]->GetMesh_Coord(2);
    else
      MeshCoord[2] = 0.0;
  }
  else{
    MeshCoord[0] = 0.0;
    MeshCoord[1] = 0.0;
    MeshCoord[2] = 0.0;
  }

  MeshCoord_passive[0] = SU2_TYPE::GetValue(MeshCoord[0]);
  MeshCoord_passive[1] = SU2_TYPE::GetValue(MeshCoord[1]);
  MeshCoord_passive[2] = SU2_TYPE::GetValue(MeshCoord[2]);

  return MeshCoord_passive;

}

bool CDiscAdjSinglezoneDriver::DirectIteration(unsigned long TimeIter) {

  bool steady = (config->GetUnsteady_Simulation() == STEADY);

  config->SetDiscrete_Adjoint(false);

  switch (config->GetKind_Solver()) {
    case DISC_ADJ_EULER: config->SetKind_Solver(EULER); break;
    case DISC_ADJ_NAVIER_STOKES: config->SetKind_Solver(NAVIER_STOKES); break;
    case DISC_ADJ_RANS: config->SetKind_Solver(RANS); break;
  }

  if(TimeIter == 0) {
    /*--- Create an output structure for the direct solver ---*/
    direct_output = new COutput(config);

    /*--- Open the convergence history file ---*/
    ConvHist_file_direct = NULL;
    if (rank == MASTER_NODE){
      ConvHist_file_direct = new ofstream[nInst[ZONE_0]];
      for (iInst = 0; iInst < nInst[ZONE_0]; iInst++) {
        direct_output->SetConvHistory_Header(&ConvHist_file_direct[iInst], config, ZONE_0, iInst);
        config->SetHistFile(&ConvHist_file_direct[INST_0]);
      }
    }

    /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
    if (config->GetWrt_Unsteady() && config->GetRestart())
      ExtIter = config->GetUnst_RestartIter();

    /*--- Check for a dynamic restart (structural analysis). Update ExtIter if necessary. ---*/
    if (config->GetKind_Solver() == FEM_ELASTICITY
        && config->GetWrt_Dynamic() && config->GetRestart())
      ExtIter = config->GetDyn_RestartIter();

    /*--- Zone preprocessing ---*/
    direct_iteration->Preprocess(direct_output, integration_container, geometry_container, solver_container,
                                 numerics_container, config_container, surface_movement, grid_movement, 
                                 FFDBox, ZONE_0, INST_0);
  }

  /*--- For steady-state flow simulations, we need to loop over ExtIter for the number of time steps ---*/
  /*--- However, ExtIter is the number of FSI iterations, so nIntIter is used in this case ---*/

  /*--- For steady-state flow simulations, we need to loop over ExtIter for the number of time steps ---*/
  if (steady) config->SetExtIter(TimeIter);
  /*--- For unsteady flow simulations, we need to loop over IntIter for the number of time steps ---*/
  else config->SetIntIter(TimeIter);
  /*--- If only one internal iteration is required, the ExtIter/IntIter is the OuterIter of the block structure ---*/
  if (config->GetnIter() == 1) {
    if (steady) config->SetExtIter(config->GetOuterIter());
    else config->SetIntIter(config->GetOuterIter());
  }

  /*--- Iterate the direct solver ---*/

  direct_iteration->Iterate(direct_output, integration_container, geometry_container, solver_container, 
                            numerics_container, config_container, surface_movement, grid_movement, 
                            FFDBox, ZONE_0, INST_0);

  /*--- Postprocess the direct solver ---*/

  direct_iteration->Postprocess(direct_output, integration_container, geometry_container, solver_container,
                                numerics_container, config_container, surface_movement, grid_movement,
                                FFDBox, ZONE_0, INST_0);

  /*--- A corrector step can help preventing numerical instabilities ---*/

  if (config->GetRelaxation())
    direct_iteration->Relaxation(direct_output, integration_container, geometry_container, solver_container,
                                 numerics_container, config_container, surface_movement, grid_movement, 
                                 FFDBox, ZONE_0, INST_0);

   /*--- Update the direct solver ---*/

  direct_iteration->Update(direct_output, integration_container, geometry_container, solver_container,
                           numerics_container, config_container, surface_movement, grid_movement,
                           FFDBox, ZONE_0, INST_0);

  /*--- Monitor the direct solver ---*/
  bool StopCalc = false;

#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  UsedTime = StopTime - StartTime;

  /*--- If convergence was reached --*/
  StopCalc = integration[FLOW_SOL]->GetConvergence();

  /*--- Write the convergence history for the fluid (only screen output) ---*/

  /*--- The logic is right now case dependent ----*/
  /*--- This needs to be generalized when the new output structure comes ---*/

  if (steady) direct_output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, false, UsedTime, ZONE_0, INST_0);

  /*--- File output ---*/

  bool output_files = false;

  /*--- Determine whether a solution needs to be written
   after the current iteration ---*/

  if (

      /*--- Fixed CL problem ---*/

      ((config->GetFixed_CL_Mode()) &&
       (config->GetnExtIter()-config->GetIter_dCL_dAlpha() - 1 == TimeIter)) ||

      /*--- Steady problems ---*/

      ((((TimeIter % config->GetWrt_Sol_Freq() == 0) && (TimeIter != 0)) || (config->GetnExtIter() - 1 == TimeIter)) &&
       ((config->GetUnsteady_Simulation() == STEADY) ||
        (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) ||
        (config->GetUnsteady_Simulation() == ROTATIONAL_FRAME))) ||

      /*--- No inlet profile file found. Print template. ---*/

      (config->GetWrt_InletFile())

      ) {

    output_files = true;

  }

  /*--- Determine whether a solution doesn't need to be written
   after the current iteration ---*/

  if (config->GetFixed_CL_Mode()) {
    if (config->GetnExtIter()-config->GetIter_dCL_dAlpha() - 1 < TimeIter) output_files = false;
    if (config->GetnExtIter() - 1 == TimeIter) output_files = true;
  }

  /*--- write the solution ---*/

  if (output_files) {

    if (rank == MASTER_NODE) cout << endl << "-------------------------- File Output Summary --------------------------";

    /*--- Execute the routine for writing restart, volume solution,
     surface solution, and surface comma-separated value files. ---*/

    direct_output->SetResult_Files_Parallel(solver_container, geometry_container, config_container, TimeIter, nZone);

    /*--- Execute the routine for writing special output. ---*/
    direct_output->SetSpecial_Output(solver_container, geometry_container, config_container, TimeIter, nZone);


    if (rank == MASTER_NODE) cout << "-------------------------------------------------------------------------" << endl << endl;

  }

  config->SetDiscrete_Adjoint(true);

  switch (config->GetKind_Solver()) {
    case EULER: config->SetKind_Solver(DISC_ADJ_EULER); break;
    case NAVIER_STOKES: config->SetKind_Solver(DISC_ADJ_NAVIER_STOKES); break;
    case RANS: config->SetKind_Solver(DISC_ADJ_RANS); break;
  }

  return StopCalc;
}

CConfig* CDriver::GetConfig(unsigned short val_iZone) {
  return config_container[val_iZone];
}

CGeometry* CDriver::GetGeometry(unsigned short val_iZone, 
                                unsigned short val_iInst, 
                                unsigned short val_iMesh) {
  return geometry_container[val_iZone][val_iInst][val_iMesh];
}

CSolver* CDriver::GetSolver(unsigned short val_iZone, 
                            unsigned short val_iInst, 
                            unsigned short val_iMesh, 
                            unsigned short val_iSol) {
  return solver_container[val_iZone][val_iInst][val_iMesh][val_iSol];
}

void CDriver::Adapted_Input_Preprocessing(SU2_Comm MPICommunicator, char* confFile, vector<vector<passivedouble> > const &SolAdap,
                                          vector<vector<passivedouble> > const &PoiAdap, vector<vector<unsigned long> > const &EdgAdap, 
                                          vector<vector<unsigned long> > const &TriAdap, vector<vector<unsigned long> > const &TetAdap,
                                          unsigned short val_iZone, unsigned short val_nZone) {

  unsigned short iMesh, requestedMGlevels = config_container[val_iZone]->GetnMGLevels();
  bool fea = false;

  /*--- Determine whether or not the FEM solver is used, which decides the
   type of geometry classes that are instantiated. ---*/
  fem_solver = ((config_container[val_iZone]->GetKind_Solver() == FEM_EULER)         ||
                (config_container[val_iZone]->GetKind_Solver() == FEM_NAVIER_STOKES) ||
                (config_container[val_iZone]->GetKind_Solver() == FEM_RANS)          ||
                (config_container[val_iZone]->GetKind_Solver() == FEM_LES)           ||
                (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM_EULER) ||
                (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM_NS)    ||
                (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM_RANS));

  /*--- Read the number of instances for each zone ---*/

  nInst[val_iZone] = config_container[val_iZone]->GetnTimeInstances();

  for (iInst = 0; iInst < nInst[val_iZone]; iInst++){

    config_container[val_iZone]->SetiInst(iInst);

    /*--- Set number of markers to number in config before partitioning. ---*/

    config_container[val_iZone]->SetnMarker_All(config_container[val_iZone]->GetnMarker_CfgFile());

    /*--- De-allocate the memory of the current domain and solver, and divide the grid
     between the ranks. ---*/

    delete [] geometry_container[val_iZone][iInst];
    delete [] solver_container[val_iZone][iInst];

    Adapted_Geometrical_Preprocessing(config_container[val_iZone], geometry_container[val_iZone][iInst], 
                                      PoiAdap, EdgAdap, TriAdap, TetAdap, val_nZone);

    Adapted_Solver_Preprocessing(config_container[val_iZone], geometry_container[val_iZone][iInst], solver_container[val_iZone][iInst], SolAdap);
    
  }

}

void CDriver::Adapted_Geometrical_Preprocessing(CConfig* config, CGeometry **&geometry, vector<vector<passivedouble> > const &PoiAdap, 
                                                vector<vector<unsigned long> > const &EdgAdap, vector<vector<unsigned long> > const &TriAdap, 
                                                vector<vector<unsigned long> > const &TetAdap, unsigned short val_nZone) {

  unsigned short iMesh, requestedMGlevels = config->GetnMGLevels();
  bool fea = false;

  /*--- Definition of the geometry class to store the primal grid in the
   partitioning process. ---*/

  CGeometry *geometry_aux = NULL;

  /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

  geometry_aux = new CPhysicalGeometry(PoiAdap, EdgAdap, TriAdap, TetAdap, config, nDim, val_nZone);

  /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

  geometry_aux->SetColorGrid_Parallel(config);


  geometry = new CGeometry*[config->GetnMGLevels()+1];

  /*--- Build the grid data structures using the ParMETIS coloring. ---*/
    
  geometry[MESH_0] = new CPhysicalGeometry(geometry_aux, config);

  /*--- Deallocate the memory of geometry_aux and solver_aux ---*/

  delete geometry_aux;

  /*--- Add the Send/Receive boundaries ---*/
  geometry[MESH_0]->SetSendReceive(config);

  /*--- Add the Send/Receive boundaries ---*/
  geometry[MESH_0]->SetBoundaries(config);

  fea = ((config->GetKind_Solver() == FEM_ELASTICITY) ||
        (config->GetKind_Solver() == DISC_ADJ_FEM));

  /*--- Compute elements surrounding points, points surrounding points ---*/
  
  if (rank == MASTER_NODE) cout << "Setting point connectivity." << endl;
  geometry[MESH_0]->SetPoint_Connectivity();
  
  /*--- Renumbering points using Reverse Cuthill McKee ordering ---*/
  
  if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
  geometry[MESH_0]->SetRCM_Ordering(config);
  
  /*--- recompute elements surrounding points, points surrounding points ---*/
  
  if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
  geometry[MESH_0]->SetPoint_Connectivity();
  
  /*--- Compute elements surrounding elements ---*/
  
  if (rank == MASTER_NODE) cout << "Setting element connectivity." << endl;
  geometry[MESH_0]->SetElement_Connectivity();
  
  /*--- Check the orientation before computing geometrical quantities ---*/
  
  geometry[MESH_0]->SetBoundVolume();
  if (config->GetReorientElements()) {
    if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." << endl;
    geometry[MESH_0]->Check_IntElem_Orientation(config);
    geometry[MESH_0]->Check_BoundElem_Orientation(config);
  }
  
  /*--- Create the edge structure ---*/
  
  if (rank == MASTER_NODE) cout << "Identifying edges and vertices." << endl;
  geometry[MESH_0]->SetEdges();
  geometry[MESH_0]->SetVertex(config);
  
  /*--- Compute cell center of gravity ---*/
  
  if ((rank == MASTER_NODE) && (!fea)) cout << "Computing centers of gravity." << endl;
  geometry[MESH_0]->SetCoord_CG();
  
  /*--- Create the control volume structures ---*/
  
  if ((rank == MASTER_NODE) && (!fea)) cout << "Setting the control volume structure." << endl;
  geometry[MESH_0]->SetControlVolume(config, ALLOCATE);
  geometry[MESH_0]->SetBoundControlVolume(config, ALLOCATE);
  
  /*--- Visualize a dual control volume if requested ---*/
  
  if ((config->GetVisualize_CV() >= 0) &&
      (config->GetVisualize_CV() < (long)geometry[MESH_0]->GetnPointDomain()))
    geometry[MESH_0]->VisualizeControlVolume(config, UPDATE);
  
  /*--- Identify closest normal neighbor ---*/
  
  if (rank == MASTER_NODE) cout << "Searching for the closest normal neighbors to the surfaces." << endl;
  geometry[MESH_0]->FindNormal_Neighbor(config);
  
  /*--- Store the global to local mapping. ---*/
  
  if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
  geometry[MESH_0]->SetGlobal_to_Local_Point();
  
  /*--- Compute the surface curvature ---*/
  
  if ((rank == MASTER_NODE) && (!fea)) cout << "Compute the surface curvature." << endl;
  geometry[MESH_0]->ComputeSurf_Curvature(config);
  
  /*--- Check for periodicity and disable MG if necessary. ---*/
  
  if (rank == MASTER_NODE) cout << "Checking for periodicity." << endl;
  geometry[MESH_0]->Check_Periodicity(config);
  
  geometry[MESH_0]->SetMGLevel(MESH_0);
  if ((config->GetnMGLevels() != 0) && (rank == MASTER_NODE))
    cout << "Setting the multigrid structure." << endl;
  
  /*--- Loop over all the new grid ---*/
  
  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    
    /*--- Create main agglomeration structure ---*/
    
    geometry[iMesh] = new CMultiGridGeometry(geometry, config, iMesh);
    
    /*--- Compute points surrounding points. ---*/
    
    geometry[iMesh]->SetPoint_Connectivity(geometry[iMesh-1]);
    
    /*--- Create the edge structure ---*/
    
    geometry[iMesh]->SetEdges();
    geometry[iMesh]->SetVertex(geometry[iMesh-1], config);
    
    /*--- Create the control volume structures ---*/
    
    geometry[iMesh]->SetControlVolume(config, geometry[iMesh-1], ALLOCATE);
    geometry[iMesh]->SetBoundControlVolume(config, geometry[iMesh-1], ALLOCATE);
    geometry[iMesh]->SetCoord(geometry[iMesh-1]);
    
    /*--- Find closest neighbor to a surface point ---*/
    
    geometry[iMesh]->FindNormal_Neighbor(config);
    
    /*--- Store our multigrid index. ---*/
    
    geometry[iMesh]->SetMGLevel(iMesh);
    
    /*--- Protect against the situation that we were not able to complete
       the agglomeration for this level, i.e., there weren't enough points.
       We need to check if we changed the total number of levels and delete
       the incomplete CMultiGridGeometry object. ---*/
    
    if (config->GetnMGLevels() != requestedMGlevels) {
      delete geometry[iMesh];
      break;
    }
  }

  /*--- Create the data structure for MPI point-to-point communications. ---*/
  
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    geometry[iMesh]->PreprocessP2PComms(geometry[iMesh], config);

  /*--- Perform a few preprocessing routines and communications. ---*/
  
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    
    /*--- Compute the max length. ---*/
    
    if ((rank == MASTER_NODE) && (!fea) && (iMesh == MESH_0)) cout << "Finding max control volume width." << endl;
    geometry[iMesh]->SetMaxLength(config);
    
    /*--- Communicate the number of neighbors. This is needed for
         some centered schemes and for multigrid in parallel. ---*/
    
    if ((rank == MASTER_NODE) && (size > SINGLE_NODE) && (!fea) && (iMesh == MESH_0)) cout << "Communicating number of neighbors." << endl;
    geometry[iMesh]->InitiateComms(geometry[iMesh], config, NEIGHBORS);
    geometry[iMesh]->CompleteComms(geometry[iMesh], config, NEIGHBORS);
  }

}

void CDriver::Adapted_Solver_Preprocessing(CConfig* config, CGeometry **geometry, CSolver ***&solver, vector<vector<passivedouble> > const &SolAdap) {

  unsigned short iSol;
  int val_iter = 0;
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------ Solver Preprocessing ( Zone " << config->GetiZone() <<" ) ------------------" << endl;

  solver = new CSolver**[config->GetnMGLevels()+1];
  
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    solver[iMesh] = NULL;
  
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    solver[iMesh] = new CSolver* [MAX_SOLS];
    for (iSol = 0; iSol < MAX_SOLS; iSol++)
      solver[iMesh][iSol] = NULL;
  }
  
  unsigned short iMGlevel;
  bool euler, ns, turbulent,
      fem_euler, fem_ns, fem_turbulent, fem_transition,
      adj_euler, adj_ns, adj_turb,
      heat_fvm,
      fem, disc_adj_fem,
      spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
      template_solver, disc_adj, disc_adj_turb, disc_adj_heat,
      fem_dg_flow, fem_dg_shock_persson,
      e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras;
  
  /*--- Count the number of DOFs per solution point. ---*/
  
  DOFsPerPoint = 0;
  
  /*--- Initialize some useful booleans ---*/
  
  euler            = false;  ns              = false;  turbulent     = false;
  fem_euler        = false;  fem_ns          = false;  fem_turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb      = false;
  spalart_allmaras = false;  menter_sst      = false;  disc_adj_turb = false;
  neg_spalart_allmaras = false;
  disc_adj         = false;
  fem              = false;  disc_adj_fem     = false;
  heat_fvm         = false;  disc_adj_heat    = false;
  transition       = false;  fem_transition   = false;
  template_solver  = false;
  fem_dg_flow      = false;  fem_dg_shock_persson = false;
  e_spalart_allmaras = false; comp_spalart_allmaras = false; e_comp_spalart_allmaras = false;
  
  bool compressible   = false;
  bool incompressible = false;

  /*--- Adjust iteration number for unsteady restarts. ---*/

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  bool adjoint = (config->GetDiscrete_Adjoint() || config->GetContinuous_Adjoint());
  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC); // Dynamic simulation (FSI).

  if (dual_time) {
    if (adjoint) val_iter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
    else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
      val_iter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
    else val_iter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
  }

  if (time_stepping) {
    if (adjoint) val_iter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
    else val_iter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
  }
  
  /*--- Assign booleans ---*/
  
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; compressible = true; break;
    case NAVIER_STOKES: ns = true; compressible = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case RANS : ns = true; turbulent = true; compressible = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case INC_EULER : euler = true; incompressible = true; break;
    case INC_NAVIER_STOKES: ns = true; incompressible = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case INC_RANS : ns = true; turbulent = true; incompressible = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : fem_euler = true; compressible = true; break;
    case FEM_NAVIER_STOKES: fem_ns = true; compressible = true; break;
    case FEM_RANS : fem_ns = true; fem_turbulent = true; compressible = true; if(config->GetKind_Trans_Model() == LM) fem_transition = true; break;
    case FEM_LES : fem_ns = true; compressible = true; break;
    case HEAT_EQUATION_FVM: heat_fvm = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; compressible = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); compressible = true; adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; compressible = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; compressible = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; compressible = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; compressible = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_INC_EULER: euler = true; disc_adj = true; incompressible = true; break;
    case DISC_ADJ_INC_NAVIER_STOKES: ns = true; disc_adj = true; incompressible = true; heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_INC_RANS: ns = true; turbulent = true; disc_adj = true; incompressible = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat_fvm = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: fem_euler = true; disc_adj = true; compressible = true; break;
    case DISC_ADJ_FEM_NS: fem_ns = true; disc_adj = true; compressible = true; break;
    case DISC_ADJ_FEM_RANS: fem_ns = true; fem_turbulent = true; disc_adj = true; compressible = true; if(config->GetKind_Trans_Model() == LM) fem_transition = true; break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; compressible = true; break;
    case DISC_ADJ_HEAT: heat_fvm = true; disc_adj_heat = true; break;
  }
  
  /*--- Determine the kind of FEM solver used for the flow. ---*/
  
  switch( config->GetKind_FEM_Flow() ) {
    case DG: fem_dg_flow = true; break;
  }
  
  /*--- Determine the kind of shock capturing method for FEM DG solver. ---*/
  
  switch( config->GetKind_FEM_DG_Shock() ) {
    case PERSSON: fem_dg_shock_persson = true; break;
  }
  
  /*--- Assign turbulence model booleans ---*/
  
  if (turbulent || fem_turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:        spalart_allmaras = true;        break;
      case SA_NEG:    neg_spalart_allmaras = true;    break;
      case SA_E:      e_spalart_allmaras = true;      break;
      case SA_COMP:   comp_spalart_allmaras = true;   break;
      case SA_E_COMP: e_comp_spalart_allmaras = true; break;
      case SST:       menter_sst = true;              break;
      case SST_SUST:  menter_sst = true;              break;
      default: SU2_MPI::Error("Specified turbulence model unavailable or none selected", CURRENT_FUNCTION); break;
    }
  
  /*--- Definition of the Class for the solution: solver[DOMAIN][INSTANCE][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/
  
  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    
    /*--- Allocate solution for a template problem ---*/
    
    if (template_solver) {
      solver[iMGlevel][TEMPLATE_SOL] = new CTemplateSolver(geometry[iMGlevel], config);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][TEMPLATE_SOL]->GetnVar();
    }
    
    /*--- Allocate solution for direct problem, and run the preprocessing and postprocessing ---*/
    
    if (euler) {
      if (compressible) {
        solver[iMGlevel][FLOW_SOL] = new CEulerSolver(geometry[iMGlevel], config, iMGlevel);
        solver[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      }
      if (incompressible) {
        solver[iMGlevel][FLOW_SOL] = new CIncEulerSolver(geometry[iMGlevel], config, iMGlevel);
        solver[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      }
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][FLOW_SOL]->GetnVar();
    }
    if (ns) {
      if (compressible) {
        solver[iMGlevel][FLOW_SOL] = new CNSSolver(geometry[iMGlevel], config, iMGlevel);
      }
      if (incompressible) {
        solver[iMGlevel][FLOW_SOL] = new CIncNSSolver(geometry[iMGlevel], config, iMGlevel);
      }
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][FLOW_SOL]->GetnVar();
    }
    if (turbulent) {
      if (spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras || neg_spalart_allmaras) {
        solver[iMGlevel][TURB_SOL] = new CTurbSASolver(geometry[iMGlevel], config, iMGlevel, solver[iMGlevel][FLOW_SOL]->GetFluidModel() );
        solver[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel);
      }
      else if (menter_sst) {
        solver[iMGlevel][TURB_SOL] = new CTurbSSTSolver(geometry[iMGlevel], config, iMGlevel);
        solver[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel);
        solver[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      }
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][TURB_SOL]->GetnVar();
      if (transition) {
        solver[iMGlevel][TRANS_SOL] = new CTransLMSolver(geometry[iMGlevel], config, iMGlevel);
        if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][TRANS_SOL]->GetnVar();
      }
    }
    if (fem_euler) {
      if( fem_dg_flow ) {
        if( fem_dg_shock_persson ) {
          solver[iMGlevel][FLOW_SOL] = new CFEM_DG_NSSolver(geometry[iMGlevel], config, iMGlevel);
        }
        else {
          solver[iMGlevel][FLOW_SOL] = new CFEM_DG_EulerSolver(geometry[iMGlevel], config, iMGlevel);
        }
      }
    }
    if (fem_ns) {
      if( fem_dg_flow )
        solver[iMGlevel][FLOW_SOL] = new CFEM_DG_NSSolver(geometry[iMGlevel], config, iMGlevel);
    }
    if (fem_turbulent) {
      SU2_MPI::Error("Finite element turbulence model not yet implemented.", CURRENT_FUNCTION);
      
      if(fem_transition)
        SU2_MPI::Error("Finite element transition model not yet implemented.", CURRENT_FUNCTION);
    }
    if (heat_fvm) {
      solver[iMGlevel][HEAT_SOL] = new CHeatSolverFVM(geometry[iMGlevel], config, iMGlevel);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][HEAT_SOL]->GetnVar();
    }
    if (fem) {
      solver[iMGlevel][FEA_SOL] = new CFEASolver(geometry[iMGlevel], config);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][FEA_SOL]->GetnVar();
    }
    
    /*--- Allocate solution for adjoint problem ---*/
    
    if (adj_euler) {
      if (compressible) {
        solver[iMGlevel][ADJFLOW_SOL] = new CAdjEulerSolver(geometry[iMGlevel], config, iMGlevel);
      }
      if (incompressible) {
        SU2_MPI::Error("Continuous adjoint for the incompressible solver is not currently available.", CURRENT_FUNCTION);
      }
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJFLOW_SOL]->GetnVar();
    }
    if (adj_ns) {
      if (compressible) {
        solver[iMGlevel][ADJFLOW_SOL] = new CAdjNSSolver(geometry[iMGlevel], config, iMGlevel);
      }
      if (incompressible) {
        SU2_MPI::Error("Continuous adjoint for the incompressible solver is not currently available.", CURRENT_FUNCTION);
      }
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJFLOW_SOL]->GetnVar();
    }
    if (adj_turb) {
      solver[iMGlevel][ADJTURB_SOL] = new CAdjTurbSolver(geometry[iMGlevel], config, iMGlevel);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJTURB_SOL]->GetnVar();
    }
    
    if (disc_adj) {
      solver[iMGlevel][ADJFLOW_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver[iMGlevel][FLOW_SOL], RUNTIME_FLOW_SYS, iMGlevel);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJFLOW_SOL]->GetnVar();
      if (disc_adj_turb) {
        solver[iMGlevel][ADJTURB_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver[iMGlevel][TURB_SOL], RUNTIME_TURB_SYS, iMGlevel);
        if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJTURB_SOL]->GetnVar();
      }
      if (heat_fvm) {
        solver[iMGlevel][ADJHEAT_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver[iMGlevel][HEAT_SOL], RUNTIME_HEAT_SYS, iMGlevel);
        if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJHEAT_SOL]->GetnVar();
      }
    }
    
    if (disc_adj_fem) {
      solver[iMGlevel][ADJFEA_SOL] = new CDiscAdjFEASolver(geometry[iMGlevel], config, solver[iMGlevel][FEA_SOL], RUNTIME_FEA_SYS, iMGlevel);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJFEA_SOL]->GetnVar();
    }
    
    if (disc_adj_heat) {
      solver[iMGlevel][ADJHEAT_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver[iMGlevel][HEAT_SOL], RUNTIME_HEAT_SYS, iMGlevel);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJHEAT_SOL]->GetnVar();
    }
  }

  /*--- Use the SortAdaptedSolution routines. ---*/
  
  bool update_geo = true;
  if (config->GetFSI_Simulation()) update_geo = false;
  
  /*--- Load adapted solutions for any of the active solver containers. Note that
   these restart routines fill the fine grid and interpolate to all MG levels. ---*/

  if (euler || ns) {
    solver[MESH_0][FLOW_SOL]->SortAdaptedSolution(geometry, solver, config, SolAdap, val_iter, update_geo);
  }
  if (turbulent) {
    solver[MESH_0][TURB_SOL]->SortAdaptedSolution(geometry, solver, config, SolAdap, val_iter, update_geo);
  }
  if (fem) {
    if (dynamic) val_iter = SU2_TYPE::Int(config->GetDyn_RestartIter())-1;
    solver[MESH_0][FEA_SOL]->SortAdaptedSolution(geometry, solver, config, SolAdap, val_iter, update_geo);
  }
  if (fem_euler || fem_ns) {
    if (fem_dg_flow)
      solver[MESH_0][FLOW_SOL]->SortAdaptedSolution(geometry, solver, config, SolAdap, val_iter, update_geo);
  }
  if (heat_fvm) {
    solver[MESH_0][HEAT_SOL]->SortAdaptedSolution(geometry, solver, config, SolAdap, val_iter, update_geo);
  }
  if (adj_euler || adj_ns) {
    solver[MESH_0][ADJFLOW_SOL]->SortAdaptedSolution(geometry, solver, config, SolAdap, val_iter, update_geo);
  }
  if (disc_adj) {
    solver[MESH_0][ADJFLOW_SOL]->SortAdaptedSolution(geometry, solver, config, SolAdap, val_iter, update_geo);
    if (disc_adj_turb)
      solver[MESH_0][ADJTURB_SOL]->SortAdaptedSolution(geometry, solver, config, SolAdap, val_iter, update_geo);
    if (disc_adj_heat)
      solver[MESH_0][ADJHEAT_SOL]->SortAdaptedSolution(geometry, solver, config, SolAdap, val_iter, update_geo);
  }
  if (disc_adj_fem) {
      if (dynamic) val_iter = SU2_TYPE::Int(config->GetDyn_RestartIter())-1;
      solver[MESH_0][ADJFEA_SOL]->SortAdaptedSolution(geometry, solver, config, SolAdap, val_iter, update_geo);
  }
  if (disc_adj_heat) {
    solver[MESH_0][ADJHEAT_SOL]->SortAdaptedSolution(geometry, solver, config, SolAdap, val_iter, update_geo);
  }
  
  /*--- Set up any necessary inlet profiles ---*/
  
  Inlet_Preprocessing(solver, geometry, config);

}
